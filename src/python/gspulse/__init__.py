"""GSPulse Python API for configuration and inputs."""

import json
import logging
from pathlib import Path

from gspulse.lib import get_gspulse_root, parse_state, process_state_parsed, recursive_ndarray_to_list

# Set up package logger
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def load_and_process(gui_state_fp: Path) -> dict:
    """Load saved GUI state from a json and process it into the final format consumed by GSPulse.

    Some background info on GUI "state" definitions. All state definitions are a dict of dicts, with
    upper level dict keys corresponding to different input groupings (e.g. "settings", "profiles", "targets", etc.).
    The lower level dicts contain individual input variables.

    Summary:
    -------
    state_tk        := variables are all tkinter variables (e.g. tk.StringVar, tk.DoubleVar, etc.) corresponding to GUI entries
    state           := variables are the result of tkvar.get(). This is GUI state configuration that is saved to json file.
    state_parsed    := variables are all parsed into python objects (strings, arrays, floats, etc.)
    state_processed := parsed + additional processing such as reading shape data from files, replacing inf times with edge values, etc.
                       This is what is consumed by GSPulse (after saving to json).

    Functional relationships:
    -------------------------
    state_tk        = GSPulseApp().state_tk, see GSPulseApp.get_default_state_tk()
    state           = get_state(state_tk)
    state_parsed    = parse_state(state)
    state_processed = process_state_parsed(state_parsed)
    """
    with gui_state_fp.open() as fid:
        state = json.load(fid)
    return process_state_parsed(parse_state(state))


def startup_octave() -> object:
    """Start an octave/Oct2Py instance."""
    import oct2py  # noqa: PLC0415, import is verrrry slow don't import unless necessary

    return oct2py.Oct2Py()


def run_gspulse_from_config_file(gui_state_fp: Path, output_dir: Path, oc_instance: object = None) -> None:
    """Load GSPulse configuration from a file and run GSPulse."""
    state_processed = load_and_process(gui_state_fp)
    return run_gspulse_from_config_dict(state_processed, output_dir, oc_instance)


def run_gspulse_from_config_dict(gspulse_config: dict, output_dir: Path, oc_instance: object = None, save: bool = True) -> None:  # noqa: FBT002, FBT001
    """Run GSPulse from a configuration dictionary."""
    if oc_instance is None:
        oc_instance = startup_octave()

    tokamak = gspulse_config["settings"]["tokamak"]
    pulse_id = gspulse_config["settings"]["pulse_id"]

    # save processed state to file
    output_dir.mkdir(parents=True, exist_ok=True)
    inputs_fp = output_dir / f"{pulse_id}_input_config.json"
    fid = inputs_fp.open("w")
    data = recursive_ndarray_to_list(gspulse_config)  # for json compatibility
    fid.write(json.dumps(data, indent=4))
    logger.info(f"GSPulse inputs saved to {fid.name}")
    fid.close()

    # construct Octave commands
    cmds = []
    cmds = [
        f"addpath('{get_gspulse_root() / 'src' / 'matlab'}');",
        f"startup_gspulse('{tokamak}');",
        f"[soln, gsp_inputs] = run_pulse('{tokamak}', '{pulse_id}', 'gui_data_fp', '{inputs_fp}');",
    ]
    if save:
        soln_fp = output_dir / f"{pulse_id}_soln.mat"
        cmds.append("save('-v7', '{soln_fp}', 'soln');")
    else:
        soln_fp = None

    # run GSPulse through Octave
    oc_instance.eval("\n".join(cmds))
    logger.info("GSPulse run complete.")
    if save:
        logger.info("Solution saved to:", soln_fp)

    return soln_fp, inputs_fp


__all__ = [
    "get_gspulse_root",
    "load_and_process",
    "parse_state_raw",
    "process_state_parsed",
    "run_gspulse_from_config_dict",
    "run_gspulse_from_config_file",
    "startup_octave",
]
