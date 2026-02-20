"""Contains utility functions for the GSPulseApp() GUI"""
import numpy as np
import json
from pathlib import Path
import copy
from contextlib import suppress

def get_state(state_tk):
    """
    All state definitions are a dict of dicts, with upper level dict keys corresponding to different
    input groupings (e.g. "settings", "profiles", "targets", etc.). The lower level dicts contain individual 
    input variables. 

    state_tk        := variables are all tkinter variables (e.g. tk.StringVar, tk.DoubleVar, etc.) corresponding to GUI entries 
    state           := variables are the result of tkvar.get(). This is pulse configuration that is saved to json file. 
    state_parsed    := variables are all parsed into python objects (strings, arrays, floats, etc.)
    state_processed := parsed + additional processing such as reading shape data from files, replacing inf times with edge values, etc. This is what is consumed by GSPulse (after saving to json). 

    Functional relationships: 
    state_tk        = GSPulseApp().state_tk, see GSPulseApp.get_default_state_tk()
    state           = get_state(state_tk)
    state_parsed    = parse_state(state)
    state_processed = process_state_parsed(state_parsed)
    """
    state = {}
    for k, v in state_tk.items():
        state[k] = {subk: subv.get() for subk, subv in v.items()}    
    return state

def parse_state(state):    
    state_parsed = {}
    for k, v in state.items():
        state_parsed[k] = {}
        for subk, subv in v.items():
            state_parsed[k][subk] = parse_variable(subv, subk)
    return state_parsed

def check_replace_infs_all_time_signals(state):

    # For any variables that represent time, replace infs with edge values
    for k, v in state.items():
        for subk, subv in v.items():
            if subk.endswith("time") and subv is not None:
                state[k][subk] = replace_inf_times_with_edge_values(state['settings'], subv)

    return state

def process_state_parsed(state_parsed):
    
    state_processed = copy.deepcopy(state_parsed)
    state_processed = check_replace_infs_all_time_signals(state_processed)

    # pulse_id should be a string reprsenting an int (or just some naming string), but not a string that represents a float or some other expression
    with suppress(ValueError):
        state_processed["settings"]["pulse_id"] = str(int(state_processed["settings"]["pulse_id"]))


    # read shape array data and add it to state
    times, shape_array_data = read_multishape(state_processed)
    state_processed["shape_array_data"] = shape_array_data
    state_processed["shape_array_time"] = times

    # Convert shape-related inputs from "channel-mode" to "manual_mode".
    # That is, if (r,z) data is specified via reference to the shape file quantities
    # like 'cp_r', 'cp_z', 'R_x_pt1', 'Z_x_pt1', etc., then read the shape files
    # and write the actual (r,z) data into the inputs dict
    def _convert_channel_mode_to_manual_mode(  
        state_processed: dict, groupname: str, prefix: str, times: np.ndarray, shape_array_data: np.ndarray
    ) -> dict:
        if not state_processed[groupname][f"{prefix}_is_manual_mode"]:
            channel_eval_str = state_processed[groupname][f"{prefix}_channels"]
            if channel_eval_str is not None:
                channels = eval(channel_eval_str)  # noqa: S307
                data = np.vstack([shape_array_data[ch] for ch in channels])
                state_processed[groupname][f"{prefix}_time"] = np.array(times).tolist()
                state_processed[groupname][f"{prefix}_data"] = np.array(data).tolist()
                state_processed[groupname][f"{prefix}_is_manual_mode"] = True
        return state_processed

    for groupname in [f"rel_flux{i}" for i in range(1, 4)]:
        for prefix in ["r", "z", "r_multiref", "z_multiref", "wt_multiref"]:
            state_processed = _convert_channel_mode_to_manual_mode(state_processed, groupname, prefix, times, shape_array_data)

    for groupname in ["flux_abs_avg", "vac_flux_abs_avg", "Br", "Bz", "Br_vac", "Bz_vac"]:
        for prefix in ["r", "z"]:
            state_processed = _convert_channel_mode_to_manual_mode(state_processed, groupname, prefix, times, shape_array_data)
    
    # Convert profile editor inputs from coeffs to full profiles
    list_p_coeffs = dict_of_arrays_to_list_of_dicts(state_processed["pprime_profiles"])
    list_ff_coeffs = dict_of_arrays_to_list_of_dicts(state_processed["ffprime1_profiles"])

    # pprime profile: read Time vs Data
    pprime = [None] * len(list_p_coeffs)
    for i, p_coeffs in enumerate(list_p_coeffs):
        psin, pprime[i] = pprime_from_coeffs(p_coeffs)
    state_processed["pprime"] = {"Time": state_processed["pprime_profiles"]["time"], "Data": np.vstack(pprime), "psin": psin}

    # ffprime1 profile: read Time vs Data
    ffprime = [None] * len(list_ff_coeffs)
    for i, ff_coeffs in enumerate(list_ff_coeffs):
        psin, ffprime[i] = ffprime_from_coeffs(ff_coeffs)
    state_processed["ffprime1"] = {"Time": state_processed["ffprime1_profiles"]["time"], "Data": np.vstack(ffprime), "psin": psin}

    # ffprime2 profile: convert single expression into Time vs Data format
    x = state_processed["ffprime2_editor"]["psin"]  # TODO(jwai): safe eval, don't rely on psin variable being called 'x'
    ffprime2 = eval(state_processed["ffprime2_editor"]["ffprime2_expression"])  # noqa: S307
    state_processed["ffprime2"] = {
        "Time": replace_inf_times_with_edge_values(state_processed["settings"], [-np.inf, np.inf]),
        "Data": np.vstack([ffprime2, ffprime2]),
        "psin": x,
    }

    return state_processed

def parse_variable(val, nam):    
    if isinstance(val, bool) or isinstance(val, int):
        pass  # do nothing, val is already formatted as boolean or int
    elif isinstance(val, str):
        if len(val) == 0:
            val = None   # val was not inputted (empty string)
        else:
            try:
                val = float(val)  # val represents a numeric scalar, e.g. val='10.1'
            except ValueError:
                try: # noqa: SIM105
                    # val represents a list/array, e.g. val='np.array([0,1,2])'
                    val = np.array(eval(val)).astype(float)  # noqa: S307
                except (NameError, ValueError, SyntaxError):
                     pass  # keep as-is, val actually does represent a string, e.g. val='python_cvxopt'
    else:
        raise TypeError(f"Variable '{nam}' with type '{type(val)}' could not be parsed")
    
    return val

# NOTE(mveldhoen): It is assumed that the user will never change the root,
# if that becomes a problem, we need to add a setter too.
def get_gspulse_root() -> Path:
    """Get the root directory of the GSPulse package."""
    # GSPulse root directory (3 levels up)
    return Path(__file__).parents[3].resolve()

def ffprime_from_coeffs(coeffs: dict) -> tuple[np.ndarray, np.ndarray]:
    """Generate ffprime profile from coefficient dict."""
    psin = np.linspace(0, 1, 101)
    ffprime = (1.0 - psin) ** coeffs["exp"]
    flat_idx = psin < coeffs["flat_val"]
    if np.any(flat_idx):
        flat_val = ffprime[flat_idx][-1]
        ffprime[flat_idx] = flat_val + np.linspace(coeffs["core_val"], 0, np.sum(flat_idx))
    with np.errstate(over="ignore"):
        coeffs["ped_width"] = max(coeffs["ped_width"], 1e-3)
        pedestal_ff = coeffs["ped_height"] * (1 + np.tanh(-1 * (psin - coeffs["ped_center"]) / coeffs["ped_width"]))
    ffprime = ffprime + pedestal_ff
    ffprime[-1] = 0

    return psin, ffprime


def pprime_from_coeffs(coeffs: dict) -> tuple[np.ndarray, np.ndarray]:
    """Generate pprime profile from coefficient dict."""
    psin = np.linspace(0, 1, 101)
    pprime = (1.0 - psin) ** coeffs["exp"]
    pprime = pprime * psin ** coeffs["exp2"]
    with np.errstate(over="ignore"):
        coeffs["ped_width"] = max(coeffs["ped_width"], 1e-3)
        pedestal = coeffs["ped_height"] * (1 / np.cosh((psin - coeffs["ped_center"]) / coeffs["ped_width"])) ** 2
    pprime = pprime + pedestal
    pprime[-1] = 0

    return psin, pprime


def dict_of_arrays_to_list_of_dicts(din: dict) -> list[dict]:
    """Convert a dict of arrays to a list of dicts."""
    l_ = []
    n = len(next(iter(din.values())))
    for i in range(n):
        d = {}
        for k, v in din.items():
            d[k] = v[i]
        l_.append(d)
    return l_


def read_multishape(state_processed: dict) -> tuple[np.ndarray, dict]:
    """Read shape data from multiple shape files and marshal into arrays."""
    times, shape_fps = get_multishape_fps(state_processed)
    shape_array_data = {}
    for i, shape_fp in enumerate(shape_fps):
        shape = json.load(shape_fp.open())
        tmp = {**shape["manual_points"], **shape["derived_shape_data"]}
        tmp.pop("segs", None)

        for k, v in tmp.items():
            if i == 0:
                val = np.asarray(v).T
                shape_array_data[k] = np.full((val.size, len(times)), np.nan)
            shape_array_data[k][:, i] = v

    n = shape_array_data["cp_r"].shape[0]
    for i in range(n):
        shape_array_data[f"cp_r{i + 1}"] = shape_array_data["cp_r"][i, np.newaxis]
        shape_array_data[f"cp_z{i + 1}"] = shape_array_data["cp_z"][i, np.newaxis]

    # Marshal strike_pt, x_pt, control_pt_ref
    nref = 3
    nx = 4
    nstrike = 8
    shape_array_data["R_strike_pt"] = np.full((nstrike, len(times)), np.nan)
    shape_array_data["Z_strike_pt"] = np.full((nstrike, len(times)), np.nan)
    shape_array_data["R_x_pt"] = np.full((nx, len(times)), np.nan)
    shape_array_data["Z_x_pt"] = np.full((nx, len(times)), np.nan)
    shape_array_data["R_control_pt_ref"] = np.full((nref, len(times)), np.nan)
    shape_array_data["Z_control_pt_ref"] = np.full((nref, len(times)), np.nan)
    for i in range(nref):
        shape_array_data["R_control_pt_ref"][i, :] = shape_array_data[f"R_control_pt_ref{i + 1}"]
        shape_array_data["Z_control_pt_ref"][i, :] = shape_array_data[f"Z_control_pt_ref{i + 1}"]
    for i in range(nx):
        shape_array_data["R_x_pt"][i, :] = shape_array_data[f"R_x_pt{i + 1}"]
        shape_array_data["Z_x_pt"][i, :] = shape_array_data[f"Z_x_pt{i + 1}"]
    for i in range(nstrike):
        shape_array_data["R_strike_pt"][i, :] = shape_array_data[f"R_strike_pt{i + 1}"]
        shape_array_data["Z_strike_pt"][i, :] = shape_array_data[f"Z_strike_pt{i + 1}"]

    return times, shape_array_data


def get_multishape_fps(state_processed: dict) -> tuple[np.ndarray, list[Path]]:
    """Read the shape filepaths and times from processed state and return as arrays."""
   
    times = state_processed["multishape"]["shape_time"]
    shape_fns = state_processed["multishape"]["shape_filenames"].split(",\n")    
    multishape_dir = get_gspulse_root() / Path(state_processed["multishape"]["multishape_dir"])
    shape_fps = [multishape_dir / fn for fn in shape_fns]

    if len(times) != len(shape_fps):
        raise ValueError(f"Number of times ({len(times)}) does not match number shape files ({len(shape_fps)}).")

    return times, shape_fps


def replace_inf_times_with_edge_values(settings: dict, times: np.ndarray) -> np.ndarray:
    """Replace -inf and +inf in time arrays with global min/max time values from the settings."""
    times = np.array(times).astype(float)

    if not np.all(np.diff(times) > 0):
        raise ValueError("Times are not sorted in ascending order")

    # read global times
    global_t_max = -np.inf
    global_t_min = np.inf    
    for i in range(4):
        k = f"interval_t{i + 1}"
        v = settings[k]
        if v is not None:
            global_t_min = np.min([global_t_min, np.min(v)])
            global_t_max = np.max([global_t_max, np.max(v)])

    replace_max = np.isposinf(times[-1])
    replace_min = np.isneginf(times[0])

    replace_max_val = np.max([global_t_max, times[-2]])
    replace_min_val = np.min([global_t_min, times[1]])

    if replace_max:
        times[-1] = replace_max_val
    if replace_min:
        times[0] = replace_min_val

    return times

def recursive_ndarray_to_list(d):
    for key in d:
        if isinstance(d[key], np.ndarray):
            d[key] = d[key].tolist()
        elif isinstance(d[key], dict):
            d[key] = recursive_ndarray_to_list(d[key])
    return d
