import copy
import json
import tkinter as tk
import numpy as np
import argparse
import textwrap
from functools import partial
from pathlib import Path
from tkinter import Scale, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from gspulse.shape_funs import get_segs, seg_intersections, shape_create_deadstart
from gspulse.lib import (
    get_gspulse_root,
    get_state,
    parse_state,
    process_state_parsed,
    replace_inf_times_with_edge_values,
    read_multishape,
    pprime_from_coeffs,
    ffprime_from_coeffs,
    get_multishape_fps,
    recursive_ndarray_to_list,
    check_replace_infs_all_time_signals
)

gspulse_root = get_gspulse_root()

class GspulseApp:
    def __init__(self, tokamak):
        # define root window
        self.root = tk.Tk()
        self.define_root_window()
        self.initialize_scrollable_frame()

        # initialize global data
        self.tokamak = tokamak.lower()
        fn = gspulse_root / "tokamaks" / self.tokamak / "shapes" / f"{self.tokamak}_gui_settings.json"
        with open(fn) as fid:
            self.tokamak_settings = json.load(fid)
        self.coils = self.tokamak_settings["coils"]
        self.ncoils = len(self.coils)
        self.ncombos = self.tokamak_settings["ncombos"]
        self.shape_dir = gspulse_root / "tokamaks" / self.tokamak / "shapes"
        self.runs_dir = gspulse_root / "tokamaks" / self.tokamak / "runs"
        self.shape_dir.mkdir(parents=True, exist_ok=True)
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self._oc = None
        self.axs = {}
        self.figs = {}
        self.canvases = {}
        self.plot_callbacks_list = []

        # load app inputs
        self.state_tk = self.get_default_state_tk()

        # create notebook with tabs
        tabs, notebook = self.make_tabs(
            self.root_widget, ["settings", "plasma_scalars", "profiles", "shape_editor", "optimization_signals"], expand=1, fill="both"
        )
        
        # Settings tab
        # -----------------
        self.add_settings_panel(tabs["settings"])
        self.add_coil_limits_panel(tabs["settings"])
        self.add_fbt_params_panel(tabs["settings"])

        frame = tk.Frame(tabs["settings"], highlightbackground="gray", highlightthickness=0)
        frame.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
        buttons = [
            ("Load", self.load_all_inputs),
            ("Save", self.save_all_inputs),            
            (
                "Validate inputs",
                partial(self.run_gspulse, inputs_only=True, save_soln=False, success_msg="All inputs passed validation check."),
            ),
            ("Run GSPulse", partial(self.run_gspulse, inputs_only=False, save_soln=True)),
        ]
        for i, (text, callback) in enumerate(buttons):
            B = tk.Button(frame, text=text, command=callback, width=12, height=1)
            B.grid(row=0, column=i, padx=2, pady=2)

        # Plasma scalars tab
        # -------------------------------
        signal_names = ["Ip", "Wk", "qA", "ag"]
        xvars = ["Ip_time", "Wk_time", "qA_time", "ag_time"]
        yvars = ["Ip_data", "Wk_data", "qA_data", "ag_data"]
        ylabels = ["Ip [MA]", "Wk [MJ]", "qA", "ag"]
        for i, (signal_name, xvar, yvar, ylabel) in enumerate(zip(signal_names, xvars, yvars, ylabels)):
            nrows = 2
            frame = tk.LabelFrame(tabs["plasma_scalars"], text=signal_name, highlightbackground="gray", highlightthickness=2)
            frame.grid(row=i % nrows, column=i // nrows, padx=0, pady=0, sticky="nsew")
            self.add_signal_panel(
                frame,
                fig_name=signal_name,
                xvar=self.state_tk["plasma_scalars"][xvar],
                yvar=self.state_tk["plasma_scalars"][yvar],
                xlabel="Time [s]",
                ylabel=ylabel,
            )

        # Plasma profiles tab
        # --------------------

        # pprime editor
        fig_name = "pprime"
        app_inputs_key = "pprime_editor"
        frame = tk.LabelFrame(
            tabs["profiles"], text="pprime editor: P'=(1-psin)^exp * psin^exp2 + pedestal", highlightbackground="gray", highlightthickness=2
        )
        frame.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
        self.add_profile_editor(frame, fig_name, app_inputs_key, pprime_from_coeffs)

        # pprime time-dependent profiles
        frame = tk.LabelFrame(tabs["profiles"], text="pprime time-dependent coeffs", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=1, column=0, padx=0, pady=0, sticky="nsew")
        self.add_inputs_to_panel(
            frame, self.state_tk["pprime_profiles"], callback=partial(self.update_profile_slider_plots, 0), width=20
        )
        frame = tk.Frame(tabs["profiles"], highlightbackground="gray", highlightthickness=2)
        frame.grid(row=2, column=0, padx=0, pady=0, sticky="nsew")
        self.add_plot_to_panel(frame, name="pprime_slider", figsize=(5, 3), row=0, column=0)

        # ffprime1 editor
        fig_name = "ffprime"
        app_inputs_key = "ffprime1_editor"
        frame = tk.LabelFrame(tabs["profiles"], text="ffprime editor: FF'= ...", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=0, column=1, padx=0, pady=0, sticky="nsew")
        self.add_profile_editor(frame, fig_name, app_inputs_key, ffprime_from_coeffs)

        # ffprime time-dependent profiles
        frame = tk.LabelFrame(tabs["profiles"], text="ffprime time-dependent coeffs", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=1, column=1, padx=0, pady=0, sticky="nsew")
        self.add_inputs_to_panel(
            frame, self.state_tk["ffprime1_profiles"], callback=partial(self.update_profile_slider_plots, 0), width=20
        )
        frame = tk.Frame(tabs["profiles"], highlightbackground="gray", highlightthickness=2)
        frame.grid(row=2, column=1, padx=0, pady=0, sticky="nsew")
        self.add_plot_to_panel(frame, name="ffprime1_slider", figsize=(5, 3), row=0, column=0)
        self.add_time_slider(tabs["profiles"], slider_callback=self.update_profile_slider_plots, length=400, row=3, column=0, columnspan=2)

        # ffprime2 editor
        frame = tk.LabelFrame(tabs["profiles"], text="ffprime2 editor", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=0, column=2, padx=0, pady=0, sticky="nsew")
        self.add_signal_panel(
            frame,
            fig_name="ffprime2",
            xvar=self.state_tk["ffprime2_editor"]["psin"],
            yvar=self.state_tk["ffprime2_editor"]["ffprime2_expression"],
            xlabel="psin (:=x)",
            ylabel="ffprime2(x)",
        )

        # Shape editor tab
        # -----------------

        # shape params
        frame = tk.LabelFrame(tabs["shape_editor"], text="Shape params", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
        self.add_inputs_to_panel(frame, self.state_tk["shape_params"], callback=self.update_shape_plot, width=6)

        # manual pts
        frame = tk.LabelFrame(tabs["shape_editor"], text="Manual points", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=1, column=0, padx=0, pady=0, sticky="nsew")
        label = tk.Label(frame, text="R").grid(row=0, column=1)
        label = tk.Label(frame, text="Z").grid(row=0, column=2)

        for i, name in enumerate(self.manual_pt_names):
            label = tk.Label(frame, text=name)
            label.grid(row=i + 1, column=0)
            for j, prefix in enumerate(["R", "Z"]):
                var_name = f"{prefix}_{name}"
                entry = tk.Entry(frame, bd=0, width=4, textvariable=self.state_tk["manual_points"][var_name])
                entry.bind("<Return>", self.update_shape_plot)
                entry.grid(row=i + 1, column=j + 1)

        # control segments
        frame = tk.LabelFrame(tabs["shape_editor"], text="Control segments", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=0, column=1, padx=0, pady=0, sticky="nsew")
        self.add_inputs_to_panel(frame, self.state_tk["control_segments"], callback=self.update_shape_plot, width=6)

        # plot options
        frame = tk.LabelFrame(tabs["shape_editor"], text="Plot options", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=1, column=1, padx=0, pady=0, sticky="nsew")
        self.add_inputs_to_panel(frame, self.state_tk["shape_plot_options"], callback=self.update_shape_plot, width=6)

        # load/save buttons
        frame = tk.Frame(tabs["shape_editor"], highlightbackground="gray", highlightthickness=0)
        frame.grid(row=2, column=0, columnspan=3, padx=0, pady=0, sticky="nsew")

        buttons = [("Load default shape", self.load_default_shape), ("Load shape", self.load_shape), ("Save shape", self.save_shape)]

        for i, (text, callback) in enumerate(buttons):
            B = tk.Button(frame, text=text, command=callback, width=11, height=1)
            B.grid(row=0, column=i, padx=0, pady=0)

        # shape evolution panel
        frame = tk.LabelFrame(tabs["shape_editor"], text="Shape evolution", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=0, column=2, padx=0, pady=0, sticky="nsew")

        label = tk.Label(frame, text="Time [s]")
        label.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
        entry = tk.Entry(frame, bd=0, width=20, textvariable=self.state_tk["multishape"]["shape_time"])
        entry.grid(row=0, column=1, padx=0, pady=0, sticky="nsew")
        entry.bind("<Return>", partial(self.update_shape_plot_from_slider, 0))
        callback = partial(
            self.ask_get_filename,
            self.state_tk["multishape"]["shape_filenames"],
            self.state_tk["multishape"]["multishape_dir"],
        )
        button = tk.Button(frame, text="Select shape files", command=callback)
        button.grid(row=1, column=1, padx=0, pady=0, sticky="nsew")

        # text widget to display selected files (read-only)
        dir_text = tk.Text(frame, height=2, width=30, wrap="word")
        dir_text.insert("1.0", self.state_tk["multishape"]["multishape_dir"].get())
        dir_text.configure(state="disabled")
        dir_text.grid(row=2, column=1, columnspan=1, padx=0, pady=0, sticky="nsew")
        var = self.state_tk["multishape"]["multishape_dir"]
        var.trace_add("write", partial(self.update_text_widget, dir_text, var))

        fn_text = tk.Text(frame, height=16, width=30, wrap="word")
        fn_text.insert("1.0", self.state_tk["multishape"]["shape_filenames"].get())
        fn_text.configure(state="disabled")
        fn_text.grid(row=3, column=1, columnspan=1, padx=0, pady=0, sticky="nsew")
        var = self.state_tk["multishape"]["shape_filenames"]
        var.trace_add("write", partial(self.update_text_widget, fn_text, var))

        # add shape plot figure
        frame = tk.Frame(tabs["shape_editor"], highlightbackground="gray", highlightthickness=0)
        frame.grid(row=0, column=3, rowspan=2, padx=0, pady=0, sticky="nsew")
        self.add_plot_to_panel(frame, "shape", figsize=(6, 8), row=0, column=0)
        self.axs["shape"].set_aspect("equal")
        self.axs["shape"].set_xlim(self.tokamak_settings["xlim"])
        self.axs["shape"].set_ylim(self.tokamak_settings["ylim"])
        self.add_toolbar(frame, "shape", row=1, column=0)
        self.add_time_slider(frame, slider_callback=self.update_shape_plot_from_slider, length=400, row=2, column=0, sticky="nsew")


        # Optimization signals
        # -------------------------------
        sig_tabs, sig_notebook = self.make_tabs(
            tabs["optimization_signals"],
            [
                "voltage",
                "current",
                "current_combos",
                "rel_flux1",
                "rel_flux2",
                "rel_flux3",
                "flux_abs_avg",
                "vac_flux_abs_avg",
                "Br",
                "Bz",
                "Br_vac",
                "Bz_vac",
            ],
            expand=1,
            fill="both",
            side="left",
        )

        for name in ["rel_flux1", "rel_flux2", "rel_flux3", "flux_abs_avg", "vac_flux_abs_avg", "Br", "Bz", "Br_vac", "Bz_vac"]:
            
            # add use_in_optimization button
            tab = sig_tabs[name]
            frame = tk.Frame(tab, highlightbackground="gray", highlightthickness=0)
            frame.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
            self.add_inputs_to_panel(frame, {"use_in_optimization": self.state_tk[name]["use_in_optimization"]}, width=5)
            
            # add info button
            match name:
                case "rel_flux1" | "rel_flux2" | "rel_flux3":
                    info_text = (
                        "Relative flux signal error. The error is calculated as \n\ne = target - (flux @ (r1,z1) - flux @ (r2,z2))\n\nThe (r,z) spatial locations "
                        "can either be specified by selecting a corresponding channel from the shape data files, or by manually entering the coordinates. \n\nThe "
                        "available shape channels are: \n\ncp_r, cp_z, cp_r<#>, cp_z<#>, R_x_pt<#>, Z_x_pt<#>, R_strike_pt<#>, Z_strike_pt<#>, R_control_pt_ref<#>, Z_control_pt_ref<#>")
                    rz_prefixes = ["r", "z", "r_multiref", "z_multiref", "wt_multiref"]
                case "flux_abs_avg" | "vac_flux_abs_avg" | "Br" | "Bz" | "Br_vac" | "Bz_vac":
                    info_text = "To be added ... "
                    rz_prefixes = ["r", "z"]
    
            self.add_info_button(frame, text=info_text, row=1, column=0)
            frame = tk.Frame(tab, highlightbackground="gray", highlightthickness=0)
            frame.grid(row=1, column=0, padx=0, pady=0, sticky="nsew")

            # add rz panels and targ/wt panels
            self.add_rz_panels(frame, name, rz_prefixes)
            self.add_targ_wt_input_plot_panels(frame, name)
            
        # voltage and currents tabs                
        self.add_coil_signals(sig_tabs["voltage"], "voltage")        
        self.add_coil_signals(sig_tabs["current"], "current")

        # current combos tab
        frame = tk.Frame(sig_tabs["current_combos"], highlightbackground="gray", highlightthickness=0)
        frame.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
        self.add_inputs_to_panel(frame, {"use_in_optimization": self.state_tk["current_combos"]["use_in_optimization"]}, width=5)
        self.add_info_button(frame, text="To be added ...", row=1, column=0)
        
        frame = tk.LabelFrame(sig_tabs["current_combos"], text="Combo Matrix", highlightbackground="gray", highlightthickness=2)
        frame.grid(row=1, column=0, padx=0, pady=0, sticky="nsew")
        for i in range(self.ncombos):
            for j in range(self.ncoils):
                entry = tk.Entry(frame, bd=0, width=5, textvariable=self.state_tk["current_combo_matrix"][f"combo_{i}_{j}"])
                entry.grid(row=i + 1, column=j + 1, padx=0, pady=0, sticky="nsew")

        for i in range(self.ncombos):
            label = tk.Label(frame, text=f"combo{i}")
            label.grid(row=i + 1, column=0)

        for j, coil in enumerate(self.coils):
            label = tk.Label(frame, text=coil)
            label.grid(row=0, column=j + 1)

        frame = tk.LabelFrame(sig_tabs["current_combos"], text="Combo weights", highlightbackground="gray", highlightthickness=2)        
        frame.grid(row=2, column=0, padx=0, pady=0, sticky="nsew")
        self.add_targ_wt_input_plot_panels(frame, "current_combos")


        # Open to settings tab
        notebook.select(tabs["settings"])

    def startup_octave_gspulse(self):
        if self._oc is None:
            print("Launching Octave...")
            import oct2py  # very slow, dont import at top of file

            self._oc = oct2py.Oct2Py()
            cmds = [f"addpath('{gspulse_root / 'src' / 'matlab'}');", f"startup_gspulse('{self.tokamak}');"]
            self._oc.eval("\n".join(cmds))
    
    def run_gspulse(self, event=None, inputs_only=False, save_soln=True, success_msg=None):

        # Process state and save to state_processed to file                
        state_processed = process_state_parsed(parse_state(get_state(self.state_tk)))
        tokamak = state_processed["settings"]["tokamak"]
        id = state_processed["settings"]["pulse_id"]
        if id is None:
            id = "empty_pulse_id"

        initialfilepath = self.runs_dir / "processed_inputs" / f"{id}_gui_state_processed.json"
        initialfilepath.parent.mkdir(parents=True, exist_ok=True)
        gui_data_fp = self.save_data_via_filedialog(state_processed, askfilepath=False, initialfilepath=initialfilepath, initialdir=self.runs_dir)
    
        # Construct Octave commands
        cmds = []
        cmds.append(f"[soln, gsp_inputs] = run_pulse('{tokamak}', '{id}', 'gui_data_fp', '{gui_data_fp}', 'inputs_only', {int(inputs_only)});")
        if save_soln:
            soln_fn = self.runs_dir / "soln" / Path(id + "_soln.mat")
            soln_fn.parent.mkdir(parents=True, exist_ok=True)
            cmds.append(f"save('-v7', '{soln_fn}', 'soln');")

        # Run GSPulse in Octave
        self.startup_octave_gspulse()
        self._oc.eval("\n".join(cmds))

        # Print messages
        if not inputs_only and save_soln:
            print("GSPulse run complete. Solution saved to:", soln_fn)
        if success_msg is not None:
            print(success_msg)

    def update_all_plots(self):
        [callback() for callback in self.plot_callbacks_list]
        self.update_profile_slider_plots(0)
        self.update_shape_plot()

    def add_coil_signals(self, parent, groupname):
        
        # make frames to organize inputs
        topframe = tk.Frame(parent, highlightbackground="gray", highlightthickness=0)
        topframe.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
        botframe = tk.Frame(parent, highlightbackground="gray", highlightthickness=0)
        botframe.grid(row=1, column=0, padx=0, pady=0, sticky="nsew")
                
        # coil tabs
        coil_tabs, coil_tabs_nb = self.make_tabs(botframe, ["All (read-only)"] + self.coils)
        for coil in self.coils:
            frame = coil_tabs[coil]
            for i, var_type in enumerate(["targ", "constrain", "wt", "dwt", "d2wt"]):
                fig_name = f"{coil}_{groupname}_{var_type}"
                keys = [f"{coil}_{var_type}_time", f"{coil}_{var_type}_data"]
                vars_dict = self.subdict(self.state_tk[groupname], keys)
                vars_dict_values = list(vars_dict.values())
                input_frame = tk.Frame(frame, highlightbackground="gray", highlightthickness=0)
                input_frame.grid(row=2 * (i // 2), column=2 * (i % 2), padx=0, pady=0, sticky="nsew")
                callback = partial(self.read_vars_update_plot, fig_name, [vars_dict_values[0]], [vars_dict_values[1]], xlabel="Time[s]")
                self.plot_callbacks_list.append(callback)
                self.add_inputs_to_panel(input_frame, vars_dict, callback=callback, width=24)
                plot_frame = tk.Frame(frame, highlightbackground="gray", highlightthickness=0)
                plot_frame.grid(row=2 * (i // 2), column=2 * (i % 2) + 1, rowspan=2, padx=0, pady=0, sticky="nsew")
                self.add_plot_to_panel(plot_frame, fig_name, figsize=(4.5, 3), row=0, column=0)
                toolbar_frame = tk.Frame(frame, highlightbackground="gray", highlightthickness=0)
                toolbar_frame.grid(row=2 * (i // 2) + 1, column=2 * (i % 2), sticky="nsew")
                toolbar = NavigationToolbar2Tk(self.canvases[fig_name], toolbar_frame)
                toolbar.update()

        # header buttons and info
        info_text = "To be added ..."
        self.add_inputs_to_panel(topframe, {"use_in_optimization": self.state_tk[groupname]["use_in_optimization"]}, width=5)
        self.add_info_button(topframe, text=info_text, row=1, column=0)
        B = tk.Button(topframe, 
                      text = "Apply to all coils",
                      command=partial(self.apply_to_all_coils, coil_tabs, coil_tabs_nb, groupname),
                      width=15,
                      height=1,
                      )
        B.grid(row=2, column=0, padx=0, pady=0, sticky="nsew")

        # "all-coils" plot
        frame = coil_tabs["All (read-only)"]
        for i, var_type in enumerate(["targ", "constrain", "wt", "dwt", "d2wt"]):
            fig_name = f"All_coils_{groupname}_{var_type}"

            # add empty plots to tab
            plot_frame = tk.Frame(frame, highlightbackground="gray", highlightthickness=0)
            plot_frame.grid(row=i // 3, column=i % 3, padx=2, pady=2, sticky="nsew")
            self.add_plot_to_panel(plot_frame, fig_name, figsize=(6, 4), row=0, column=0)
            toolbar_frame = tk.Frame(plot_frame, highlightbackground="gray", highlightthickness=0)
            toolbar_frame.grid(row=1, column=0, sticky="nsew")
            toolbar = NavigationToolbar2Tk(self.canvases[fig_name], toolbar_frame)
            toolbar.update()

        # plotting callback
        def update_all_coils_plot():
            for var_type in ["targ", "constrain", "wt", "dwt", "d2wt"]:
                fig_name = f"All_coils_{groupname}_{var_type}"
                xvars = [self.state_tk[groupname][f"{coil}_{var_type}_time"] for coil in self.coils]
                yvars = [self.state_tk[groupname][f"{coil}_{var_type}_data"] for coil in self.coils]
                self.read_vars_update_plot(fig_name, xvars, yvars, xlabel="Time[s]", ylabel="", legend_labels=self.coils)

        self.plot_callbacks_list.append(update_all_coils_plot)

        # bind tab change to update plot
        def on_tab_changed(event):
            selected_tab = event.widget.tab(event.widget.select(), "text")
            if selected_tab == "All (read-only)":
                update_all_coils_plot()

        coil_tabs_nb.bind("<<NotebookTabChanged>>", on_tab_changed)

    def apply_to_all_coils(self, coil_tabs, coil_tabs_nb, groupname, event=None):
        tab_id = coil_tabs_nb.select()
        selected_coil = coil_tabs_nb.tab(tab_id, "text")

        for coil in self.coils:
            if coil != selected_coil:
                for var_type in ["targ", "constrain", "wt", "dwt", "d2wt"]:
                    for suffix in ["time", "data"]:
                        # set values
                        k2copy = f"{selected_coil}_{var_type}_{suffix}"
                        v2copy = self.state_tk[groupname][k2copy].get()
                        k = f"{coil}_{var_type}_{suffix}"
                        self.state_tk[groupname][k].set(v2copy)

                    # now update plot
                    keys = [f"{coil}_{var_type}_time", f"{coil}_{var_type}_data"]
                    vars_dict = self.subdict(self.state_tk[groupname], keys)
                    vars_dict_values = list(vars_dict.values())
                    fig_name = f"{coil}_{groupname}_{var_type}"
                    self.read_vars_update_plot(fig_name, [vars_dict_values[0]], [vars_dict_values[1]], xlabel="Time[s]")

    def add_targ_wt_input_plot_panels(self, tab, input_group_name, width=22):
        # targets, wt, dwt, d2wt
        for i, prefix in enumerate(["targ", "wt", "dwt", "d2wt"]):
            fig_name = f"{prefix}_{input_group_name}"
            vars_dict = self.subdict(
                self.state_tk[input_group_name], [f"{prefix}_time", f"{prefix}_time_val", f"{prefix}_channel_val"]
            )
            frame = tk.LabelFrame(tab, text=prefix, highlightbackground="gray", highlightthickness=2)
            frame.grid(row=i, column=1, padx=0, pady=0, sticky="nsew")
            self.add_targ_wt_input_plot_panel(frame, fig_name, vars_dict, width)

    def add_rz_panels(self, tab, input_group_name, rz_prefixes):
        for i, prefix in enumerate(rz_prefixes):
            frame = tk.LabelFrame(tab, text=prefix, highlightbackground="gray", highlightthickness=2)
            frame.grid(row=i, column=0, padx=0, pady=0, sticky="nsew")
            d = self.state_tk[input_group_name]
            channel_var_dict = self.subdict(d, [f"{prefix}_channels"])
            manual_vars_dict = self.subdict(d, [f"{prefix}_time", f"{prefix}_data"])
            channel_or_manual_var = d[f"{prefix}_is_manual_mode"]
            fig_name = f"{prefix}_{input_group_name}"
            self.add_channel_mode_input_plot_panel(frame, fig_name, channel_var_dict, manual_vars_dict, channel_or_manual_var)

    @staticmethod
    def subdict(d, keys):
        return {k: d[k] for k in keys}

    def update_targ_wts_plot(self, fig_name, vars_dict, event=None):
        vars = list(vars_dict.values())
        try:
            settings = get_vars_format_type(self.state_tk['settings'])
            t = np.array(eval(vars[0].get()))
            t = replace_inf_times_with_edge_values(settings, t)
            t_wt = np.array(eval(vars[1].get()))
            y_wt = np.array(eval(vars[2].get()))
        except Exception:
            ax = self.axs[fig_name]
            ax.clear()
            ax.text(
                0.5,
                0.5,
                "Invalid input data",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                fontsize=12,
            )
            self.figs[fig_name].tight_layout()
            self.canvases[fig_name].draw()
            return

        y = np.outer(t_wt, y_wt)
        ylist = [y[:, i] for i in range(y.shape[1])]
        xlist = [t] * len(ylist)
        legend_labels = [f"{i}" for i in range(len(ylist))]
        self.update_xy_plot(fig_name, xlist=xlist, ylist=ylist, xlabel="Time [s]", ylabel="value", legend_labels=legend_labels)

    def update_xy_plot_channel_mode(self, fig_name, channel_var_dict, event=None):
        try:
            state_processed = process_state_parsed(parse_state(get_state(self.state_tk)))
            times, shape_array_data = read_multishape(state_processed)
            channel_var = list(channel_var_dict.values())[0]
            channels = eval(channel_var.get())
            xlist = [times] * len(channels)
            ylist = [shape_array_data[k] for k in channels]
            self.update_xy_plot(fig_name, xlist=xlist, ylist=ylist, xlabel="Time [s]", ylabel="[m]", legend_labels=channels)
        except:
            self.display_invalid_plot(fig_name)

    def add_targ_wt_input_plot_panel(self, parent, fig_name, vars_dict, width):
        input_frame = tk.Frame(parent, highlightbackground="gray", highlightthickness=0)
        input_frame.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
        callback = partial(self.update_targ_wts_plot, fig_name, vars_dict)
        self.plot_callbacks_list.append(callback)
        self.add_inputs_to_panel(input_frame, vars_dict, callback=callback, width=width)
        plot_frame = tk.Frame(parent, highlightbackground="gray", highlightthickness=0)
        plot_frame.grid(row=0, column=1, rowspan=2, padx=0, pady=0, sticky="nsew")
        self.add_plot_to_panel(plot_frame, fig_name, figsize=(5, 3), row=0, column=0)
        toolbar_frame = tk.Frame(parent, highlightbackground="gray", highlightthickness=0)
        toolbar_frame.grid(row=1, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(self.canvases[fig_name], toolbar_frame)
        toolbar.update()

    def add_channel_mode_input_plot_panel(self, parent, fig_name, channel_var_dict, manual_vars_dict, channel_or_manual_var):
        input_frame = tk.Frame(parent, highlightbackground="gray", highlightthickness=0)
        input_frame.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
        tabs, _ = self.make_tabs(input_frame, ["channel_mode", "manual_mode"], expand=1, fill="both")

        # --------------------------------------------------------------
        # Bind tab selection to channel_or_manual_var
        def on_tab_changed(event):
            selected_tab = event.widget.tab(event.widget.select(), "text")
            if selected_tab == "channel_mode":
                channel_or_manual_var.set(False)
            else:
                channel_or_manual_var.set(True)

        notebook = tabs["channel_mode"].master  # The ttk.Notebook instance
        notebook.bind("<<NotebookTabChanged>>", on_tab_changed)

        # Also, update tab selection if the variable changes
        def on_var_changed(*args):
            if channel_or_manual_var.get():
                notebook.select(tabs["manual_mode"])
            else:
                notebook.select(tabs["channel_mode"])

        channel_or_manual_var.trace_add("write", on_var_changed)
        # --------------------------------------------------------------

        callback = partial(
            self.decide_channel_or_manual_and_update_plot, fig_name, channel_or_manual_var, channel_var_dict, manual_vars_dict
        )
        self.plot_callbacks_list.append(callback)
        self.add_inputs_to_panel(tabs["channel_mode"], channel_var_dict, callback=callback, width=19)
        self.add_inputs_to_panel(tabs["manual_mode"], manual_vars_dict, callback=callback, width=20)

        plot_frame = tk.Frame(parent, highlightbackground="gray", highlightthickness=0)
        plot_frame.grid(row=0, column=1, rowspan=2, padx=0, pady=0, sticky="nsew")

        self.add_plot_to_panel(plot_frame, fig_name, figsize=(5, 3), row=0, column=0)

        toolbar_frame = tk.Frame(parent, highlightbackground="gray", highlightthickness=0)
        toolbar_frame.grid(row=1, column=0, sticky="nsew")
        toolbar = NavigationToolbar2Tk(self.canvases[fig_name], toolbar_frame)
        toolbar.update()

    def decide_channel_or_manual_and_update_plot(self, fig_name, channel_or_manual_var, channel_var_dict, manual_vars_dict, event=None):
        if channel_or_manual_var.get():
            # manual mode
            manual_vars_list = list(manual_vars_dict.values())
            self.read_vars_update_plot(fig_name, [manual_vars_list[0]], [manual_vars_list[1]], xlabel="Time[s]", ylabel="[m]")
        else:
            # channel mode
            self.update_xy_plot_channel_mode(fig_name, channel_var_dict)

    def add_info_button(self, parent, text, row, column):
        # show_info function written by AI
        def show_info():
            # Create a custom top-level window for a wider info box
            info_win = tk.Toplevel()
            info_win.title("Info")
            info_win.geometry("600x300")  # Set width and height as needed

            text_widget = tk.Text(info_win, wrap="word", width=70, height=15)
            text_widget.insert("1.0", text)
            text_widget.configure(state="disabled")
            text_widget.pack(expand=True, fill="both", padx=10, pady=10)

            ok_button = tk.Button(info_win, text="OK", command=info_win.destroy)
            # Bind Enter key to close the info window
            info_win.bind("<Return>", lambda event: info_win.destroy())
            ok_button.focus_set()
            ok_button.pack(pady=(0, 10))

        info_button = tk.Button(parent, text="Info", command=show_info, width=0, height=0)
        info_button.grid(row=row, column=column, padx=0, pady=0, sticky="nsew")

    def add_toolbar(self, parent, fig_name, row, column):
        toolbar_frame = tk.Frame(parent)
        toolbar_frame.grid(row=row, column=column, sticky="nsew")
        toolbar = NavigationToolbar2Tk(self.canvases[fig_name], toolbar_frame)
        toolbar.update()        

    def update_xy_plot(self, fig_name, xlist, ylist, xlabel=None, ylabel=None, title=None, legend_labels=None):
        do_legend = True
        if legend_labels is None:
            do_legend = False
            legend_labels = [None] * len(ylist)

        ax = self.axs[fig_name]
        ax.clear()
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_title(title, fontsize=12)
        for i, (x, y, label) in enumerate(zip(xlist, ylist, legend_labels)):
            args = {"linestyle": "-", "marker": "o", "markersize": 3, "color": f"C{i}"}
            if len(y.shape) < 2:
                ax.plot(x, y, label=label, **args)
            else:
                ax.plot(x, y[0], label=label, **args)  # use legend label only on the first line
                ax.plot(x, y.T, label="_", **args)

        if do_legend:
            leg = ax.legend(fontsize=8)
            leg.set_draggable(True)

        self.figs[fig_name].tight_layout()
        self.canvases[fig_name].draw()

    def make_tabs(self, parent, tabnames, **pack_kwargs):
        notebook = ttk.Notebook(parent)
        tabs = {}
        for name in tabnames:
            tabs[name] = ttk.Frame(notebook)
            notebook.add(tabs[name], text=name)
        notebook.pack(**pack_kwargs)
        return tabs, notebook

    def update_shape_plot_from_slider(self, slider_value, event=None):
        # read times and shape filenames
        try:
            state_parsed = check_replace_infs_all_time_signals(parse_state(get_state(self.state_tk)))
            # state_processed = process_state_parsed(parse_state(get_state(self.state_tk)))
            times, shape_fps = get_multishape_fps(state_parsed)
        except Exception as e:
            self.display_invalid_plot("shape", msg=e)
            return            
            
        slider_value = float(slider_value)
        slider_time = slider_value * (max(times) - min(times)) + min(times)
        idx0 = np.searchsorted(times, slider_time + np.finfo(float).eps, side="left") - 1
        if idx0 < len(times)-1:
            idx1 = idx0 + 1
            alpha = (slider_time - times[idx0]) / (times[idx1] - times[idx0])
        else:
            idx1 = idx0
            alpha = 0

        shape0 = json.load(open(shape_fps[idx0]))
        shape1 = json.load(open(shape_fps[idx1]))

        # Interpolate between shape0 and shape1
        shape = copy.deepcopy(shape0)
        shape["shape_plot_options"] = state_parsed["shape_plot_options"] # use current GUI plot view options, not whatever was saved to file 
        for groupname in ["shape_params", "manual_points", "derived_shape_data"]:
            for varname in shape[groupname].keys():
                val0 = shape0[groupname][varname]
                val1 = shape1[groupname][varname]
                if val0 is None or val1 is None:
                    shape[groupname][varname] = None
                else:
                    shape[groupname][varname] = (1 - alpha) * np.array(val0) + alpha * np.array(val1)

        self.update_shape_plot(shape=shape)
        self.axs["shape"].set_title(f"Time = {slider_time:.3f} s")
        self.canvases["shape"].draw()

    def update_shape_plot(self, event=None, shape=None):
        if shape is None:
            shape = self.read_shape_from_app_inputs()

        ax = self.axs["shape"]
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.clear()
        ax.plot(self.tokamak_settings["limiter_r"], self.tokamak_settings["limiter_z"], linewidth=1.5, color="black")

        for k in self.manual_pt_names:
            r = shape["manual_points"][f"R_{k}"]
            z = shape["manual_points"][f"Z_{k}"]
            if r is not None and z is not None:
                if "control_pt_ref" in k:
                    ax.scatter(r, z, s=70, facecolors="none", edgecolors="g", linewidths=1, marker="s")
                    if shape["shape_plot_options"]["label_cp_ref_pts"]:
                        ax.annotate(k, (r, z))
                if "x_pt" in k:
                    ax.scatter(r, z, s=30, c="b", alpha=1, linewidths=1, marker="x")
                    if shape["shape_plot_options"]["label_xpts"]:
                        ax.annotate(k, (r, z))
                if "strike_pt" in k:
                    ax.scatter(r, z, s=15, c="blue", alpha=1, marker="o")
                    if shape["shape_plot_options"]["label_strike_points"]:
                        ax.annotate(k, (r, z))

        s = shape["derived_shape_data"]
        if s["rbbbs"] is not None and s["zbbbs"] is not None:
            ax.plot(s["rbbbs"], s["zbbbs"], linewidth=1, color="blue", linestyle="-")

            if shape["shape_plot_options"]["plot_segments"] and s["segs"] is not None:
                ax.plot(s["segs"][:, [0, 2]].T, s["segs"][:, [1, 3]].T, c="blue", alpha=0.3, linewidth=0.5)

            if s["cp_r"] is not None and s["cp_z"] is not None:
                ax.scatter(s["cp_r"], s["cp_z"], s=15, c="b", alpha=1, marker=".")
                if shape["shape_plot_options"]["label_control_points"]:
                    for i in range(len(s["cp_r"])):
                        txt = str(i + 1)
                        ax.annotate(txt, (s["cp_r"][i], s["cp_z"][i]))

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_title("Time-independent")
        self.figs["shape"].tight_layout()
        self.canvases["shape"].draw()

    def load_default_shape(self, event=None):
        fn = gspulse_root / "tokamaks" / self.tokamak / "shapes" / f"{self.tokamak}_default_shape.json"
        self.load_shape(filename=fn)

    def load_shape(self, filename=None, event=None):
        self.load_data(
            groupnames=["shape_params", "control_segments", "manual_points", "shape_plot_options"],
            initialdir=self.shape_dir,
            filename=filename,
        )
        self.update_shape_plot()

    def load_data(self, groupnames, initialdir=None, filename=None, event=None):
        if filename is None:
            filetypes = [("JSON files", "*.json"), ("Text Documents", "*.txt"), ("All Files", "*.*")]
            fid = tk.filedialog.askopenfile(filetypes=filetypes, initialdir=str(initialdir))
        else:
            fid = open(filename)
        loaded_data = json.load(fid)
        for groupname in groupnames:
            for varname in self.state_tk[groupname].keys():
                if varname not in loaded_data[groupname]:
                    print(f"Warning: {groupname}['{varname}'] not found in loaded data. Skipping.")
                    continue
                val = loaded_data[groupname][varname]
                self.set_tkvar_value(self.state_tk[groupname][varname], val)
        fid.close()

    def set_tkvar_value(self, tkvar, val):
        if isinstance(tkvar, tk.BooleanVar):
            val = bool(val)
        elif isinstance(tkvar, tk.IntVar):
            val = int(val)
        elif isinstance(tkvar, tk.StringVar):
            if val is None:
                val = ""
            else:
                val = str(val)
        tkvar.set(val)

    def save_data_via_filedialog(self, data, askfilepath=True, initialfilepath=None, initialdir=None):
        data = recursive_ndarray_to_list(data)  # for json compatibility
        if askfilepath:
            if initialfilepath is None:
                initialfile = ".json"
                initialdir = initialdir
            else:
                initialdir = Path(initialfilepath).parent
                initialfile = Path(initialfilepath).name

            fid = tk.filedialog.asksaveasfile(
                initialfile=str(initialfile),
                initialdir=str(initialdir),
                defaultextension=".json",
                filetypes=[("All Files", "*.*"), ("Text Documents", "*.txt")],
            )

        else:
            fid = open(initialfilepath, "w")

        fid.write(json.dumps(data, indent=4))
        print(f"Data saved to {fid.name}")
        fid.close()

        return fid.name        
    
    def save_all_inputs(self, event=None):
        # read pulse id
        id = self.state_tk["settings"]["pulse_id"].get()
        if len(id) == 0:
            id = "empty_pulse_id"
        fp = self.runs_dir / f"{id}_gui_state.json"

        # read GUI inputs
        inputs = {}
        for k, v in self.state_tk.items():
            inputs[k] = self.get_vars(v)

        # save
        self.save_data_via_filedialog(inputs, initialfilepath=fp, initialdir=self.runs_dir)

    def load_all_inputs(self, event=None):
        self.load_data(groupnames=list(self.state_tk.keys()), initialdir=self.runs_dir)
        self.update_all_plots()

    def save_shape(self, event=None):
        shape_data = self.read_shape_from_app_inputs()
        self.save_data_via_filedialog(shape_data, initialdir=self.shape_dir)    

    def read_shape_from_app_inputs(self):
        shape_data = {}
        for k in ["shape_params", "control_segments", "manual_points", "shape_plot_options"]:
            shape_data[k] = get_vars_format_type(self.state_tk[k])

        enough_inputs_set = None not in shape_data["shape_params"].values() and None not in shape_data["control_segments"].values()
        if enough_inputs_set:
            segs = get_segs(shape_data["control_segments"])
            [rbbbs, zbbbs] = shape_create_deadstart(shape_data["shape_params"])
            [cp_r, cp_z] = seg_intersections(segs, rbbbs, zbbbs)
        else:
            segs, rbbbs, zbbbs, cp_r, cp_z = None, None, None, None, None
        shape_data["derived_shape_data"] = {"segs": segs, "rbbbs": rbbbs, "zbbbs": zbbbs, "cp_r": cp_r, "cp_z": cp_z}

        return shape_data

    def update_text_widget(self, text_widget, var, *args):
        text_widget.configure(state="normal")
        text_widget.delete("1.0", tk.END)
        text_widget.insert("1.0", var.get())
        text_widget.configure(state="disabled")

    def ask_get_filename(self, filename_var, dirname_var, event=None):
        # query user for filenames
        filepaths = tk.filedialog.askopenfilenames(initialdir=str(self.shape_dir), filetypes=[("JSON files", "*.json")])
        if filepaths:
            folder = Path(filepaths[0]).parent
            dirname_var.set(str(folder.relative_to(gspulse_root)))
            filename_var.set(",\n".join([Path(fp).name for fp in filepaths]))

    def add_profile_editor(self, parent, fig_name, app_inputs_key, profile_from_coeffs_fun):
        self.add_plot_to_panel(parent, fig_name, figsize=(5, 3), row=len(self.state_tk[app_inputs_key]), column=0, columnspan=2)
        callback = partial(self.update_profile_editor_plot, fig_name, app_inputs_key, profile_from_coeffs_fun)
        self.plot_callbacks_list.append(callback)
        self.add_inputs_to_panel(parent, self.state_tk[app_inputs_key], callback=callback, width=12)    

    def update_profile_slider_plots(self, slider_value, event=None):
        def interpolate_dict_of_arrays(slider_value, slider_range, din):
            settings = get_vars_format_type(self.state_tk['settings'])
            times = replace_inf_times_with_edge_values(settings, din["time"])
            assert np.all(np.diff(times) > 0), "times must be sorted in ascending order"
            t = np.interp(slider_value, slider_range, times[[0, -1]])
            dout = {}
            for k, v in din.items():
                if v is None:
                    dout[k] = None
                else:
                    dout[k] = np.interp(t, times, v)
            dout["time"] = t
            return dout

        slider_value = float(slider_value)
        slider_range = np.array([0, 1])

        p_coeffs_array = get_vars_format_type(self.state_tk["pprime_profiles"])
        ff_coeffs_array = get_vars_format_type(self.state_tk["ffprime1_profiles"])

        try:
            p_coeffs = interpolate_dict_of_arrays(slider_value, slider_range, p_coeffs_array)
            if None in p_coeffs.values():
                self.display_invalid_plot("pprime_slider")
            else:
                x, y = pprime_from_coeffs(p_coeffs)
                title = f"Time {p_coeffs['time']:.3f}"
                self.update_xy_plot("pprime_slider", [x], [y], xlabel="psin", ylabel="pprime", title=title)
        except:
            self.display_invalid_plot("pprime_slider")

        try:
            ff_coeffs = interpolate_dict_of_arrays(slider_value, slider_range, ff_coeffs_array)
            if None in ff_coeffs.values():
                self.display_invalid_plot("ffprime1_slider")
            else:
                x, y = ffprime_from_coeffs(ff_coeffs)
                title = f"Time {ff_coeffs['time']:.3f}"
                self.update_xy_plot("ffprime1_slider", [x], [y], xlabel="psin", ylabel="ffprime", title=title)
        except:
            self.display_invalid_plot("ffprime1_slider")

    def update_profile_editor_plot(self, fig_name, app_inputs_key, profile_from_coeffs_fun, event=None):
        coeffs = get_vars_format_type(self.state_tk[app_inputs_key])
        if None in coeffs.values():
            self.display_invalid_plot(fig_name)
            return
        x, y = profile_from_coeffs_fun(coeffs)
        self.update_xy_plot(fig_name, [x], [y], xlabel="psin", ylabel="pprime")

    def display_invalid_plot(self, fig_name, msg=None):
        if msg is None:
            msg = "Invalid data"
        ax = self.axs[fig_name]
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()        
        ax.clear()
        # Wrap the text if it's too long
        max_chars = 30
        wrapped_msg = "\n".join(textwrap.wrap(str(msg), max_chars))
        ax.text(
            0.5,
            0.5,
            wrapped_msg,
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=12,
            wrap=True,
        )
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        self.figs[fig_name].tight_layout()
        self.canvases[fig_name].draw()

    

    def add_plot_to_panel(self, parent, name, figsize, **grid_kwargs):
        self.figs[name] = Figure(dpi=80, figsize=figsize)
        self.axs[name] = self.figs[name].add_subplot(1, 1, 1)
        self.canvases[name] = FigureCanvasTkAgg(self.figs[name], master=parent)
        self.canvases[name].get_tk_widget().grid(**grid_kwargs)
        self.canvases[name].draw()

    def get_vars(self, vars_dict):
        return {k: v.get() for k, v in vars_dict.items()}    

    def add_signal_panel(self, parent, fig_name, xvar, yvar, xlabel, ylabel):
        for i, (text, var) in enumerate(zip([xlabel, ylabel], [xvar, yvar])):
            label = tk.Label(parent, text=text)
            label.grid(row=i, column=0)
            entry = tk.Entry(parent, bd=0, width=35, textvariable=var)
            callback = partial(self.read_vars_update_plot, fig_name, [xvar], [yvar], xlabel=xlabel, ylabel=ylabel)
            entry.bind("<Return>", callback)
            entry.grid(row=i, column=1)

        self.plot_callbacks_list.append(callback)
        self.add_plot_to_panel(parent, fig_name, figsize=(5, 3), row=2, column=0, columnspan=2)

    def read_vars_update_plot(self, fig_name, xvars, yvars, event=None, xlabel=None, ylabel=None, legend_labels=None):
        multiple_signals = len(yvars) > 1
        if legend_labels is None:
            legend_labels = [None] * len(yvars)

        ax = self.axs[fig_name]
        ax.clear()
        ax.set_title(fig_name, fontsize=12)
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        for xvar, yvar, legend_label in zip(xvars, yvars, legend_labels):
            try:
                settings = get_vars_format_type(self.state_tk['settings'])
                x = np.array(eval(xvar.get()))                
                x = replace_inf_times_with_edge_values(settings, x)
                y = np.array(eval(yvar.get()))
                ax.plot(x.T, y.T, "-o", markersize=3, label=legend_label)
            except:
                if multiple_signals:
                    ax.plot([], [], label=f"Invalid {legend_label} data")
                else:
                    ax.text(
                        0.5,
                        0.5,
                        "Invalid data",
                        horizontalalignment="center",
                        verticalalignment="center",
                        transform=ax.transAxes,
                        fontsize=12,
                    )
            if multiple_signals:
                leg = ax.legend(fontsize=8)
                leg.set_draggable(True)

        self.figs[fig_name].tight_layout()
        self.canvases[fig_name].draw()

    def add_time_slider(self, parent, slider_callback=None, length=400, **grid_kwargs):
        slider = Scale(
            parent,
            from_=0,
            to=1,
            resolution=0.005,
            showvalue=0,
            orient="horizontal",
            width=16,
            length=length,
            label="Time",
            command=slider_callback,
        )
        slider.grid(**grid_kwargs)

    def add_settings_panel(self, parent):
        fig_name = "interval_t"
        panel = tk.LabelFrame(parent, text="Global settings", highlightbackground="gray", highlightthickness=2)
        panel.grid(row=0, column=0, padx=0, pady=0, sticky="nsew")
        input_frame = tk.Frame(panel, highlightbackground="gray", highlightthickness=0)
        input_frame.grid(row=0, column=0, rowspan=2, padx=0, pady=0, sticky="nsew")
        plot_frame = tk.Frame(panel, highlightbackground="gray", highlightthickness=0)
        plot_frame.grid(row=0, column=1, padx=0, pady=0, sticky="nsew")
        self.add_plot_to_panel(plot_frame, fig_name, figsize=(5, 3), row=0, column=0)
        self.add_toolbar(plot_frame, fig_name, row=1, column=0)
        callback = self.update_interval_t_plot
        self.plot_callbacks_list.append(callback)
        self.add_inputs_to_panel(input_frame, self.state_tk["settings"], callback=callback, width=15)

    def update_interval_t_plot(self, event=None):
        x = []
        y = []
        varnames = ["interval_t1", "interval_t2", "interval_t3", "interval_t4"]
        y_idx = 0
        for varname in varnames:
            try:
                xi = np.array(eval(self.state_tk["settings"][varname].get()))
                yi = y_idx + np.arange(len(xi))
                y_idx += len(xi)
                y.append(yi)
                x.append(xi)
            except:
                x.append(np.array([np.nan]))
                y.append(np.array([np.nan]))
                continue

        self.update_xy_plot("interval_t", x, y, xlabel="Time[s]", ylabel="Eq index", title="Global solution times", legend_labels=varnames)

    def add_fbt_params_panel(self, parent):
        panel = tk.LabelFrame(parent, text="FBT settings", highlightbackground="gray", highlightthickness=2)
        panel.grid(row=0, column=2, padx=0, pady=0, sticky="nsew")
        self.add_inputs_to_panel(panel, self.state_tk["fbt_params"], callback=self.placeholder_update)

    def add_inputs_to_panel(self, panel, varsdict, callback=None, width=12):
        for i, (name, var) in enumerate(varsdict.items()):
            label = tk.Label(panel, text=name)
            label.grid(row=i, column=0)

            if isinstance(var, tk.BooleanVar):
                B = tk.Checkbutton(panel, text="", variable=var, command=callback)
                B.grid(row=i, column=1)
            else:
                entry = tk.Entry(panel, bd=0, width=width, textvariable=var)
                entry.bind("<Return>", callback)
                entry.grid(row=i, column=1)

    def add_coil_limits_panel(self, parent):
        panel = tk.LabelFrame(parent, text="Coils (units are [kA] and [kV])", highlightbackground="gray", highlightthickness=2)
        panel.grid(row=0, column=1, padx=0, pady=0, sticky="nsew")

        # Column labels
        for j, text in enumerate(["maxI", "minI", "maxV", "minV", "initI", "initV"]):
            label = tk.Label(panel, text=text)
            label.grid(row=0, column=j + 2)

        for i, coil in enumerate(self.coils):
            label = tk.Label(panel, text=coil)
            label.grid(row=i + 1, column=1)

            for j, setting_type in enumerate(["ic_max", "ic_min", "vmax", "vmin", "ic_init", "vinit"]):
                entry = tk.Entry(panel, bd=0, width=4, textvariable=self.state_tk[setting_type][coil])
                entry.bind("<Return>", self.placeholder_update)
                entry.grid(row=i + 1, column=j + 2)

    def define_root_window(self):
        """METHOD: define_root_window
        DESCRIPTION:
        """
        self.root.title("GSPulse Editor")

        # center the window when opened
        w = 1600  # width for the root window
        h = 850  # height for the root window
        ws = self.root.winfo_screenwidth()  # width of the screen
        hs = self.root.winfo_screenheight()  # height of the screen
        x = ws / 2 - w / 2.1
        y = hs / 2.4 - h / 2
        self.root.geometry("%dx%d+%d+%d" % (w, h, x, y))

    def placeholder_update(self, event=None):
        print("Placeholder update")

    def get_default_state_tk(self):
        state_tk = {}

        # Global settings
        state_tk["settings"] = self.initialize_tk_vars(
            strvars=[
                "pulse_id",
                "tokamak",
                "interval_t1",
                "interval_t2",
                "interval_t3",
                "interval_t4",
                "picard_algo",
                "specify_psibry_mode",
                "qpsolver",
                "qpsolver_tol",
                "spline_basis_ratio",
                "initial_vess_currents",
                "fbt_Ip_threshold",
                "tol_ic_diff", 
                "Iy_relax_factor",
                "Iy_smooth_dt"
            ],
            intvars=["plotlevel", "niter", "verbose"],
            boolvars=["do_final_boundary_trace", "use_spline_basis", "calc_strike_pts", "inject_model_offset", "calc_post_prc_ext", "boundary_post_run"],
        )

        # Coil limits
        state_tk["ic_max"] = self.initialize_tk_vars(strvars=self.coils)
        state_tk["ic_min"] = self.initialize_tk_vars(strvars=self.coils)
        state_tk["vmax"] = self.initialize_tk_vars(strvars=self.coils)
        state_tk["vmin"] = self.initialize_tk_vars(strvars=self.coils)
        state_tk["ic_init"] = self.initialize_tk_vars(strvars=self.coils)
        state_tk["vinit"] = self.initialize_tk_vars(strvars=self.coils)

        # FBT parameters
        state_tk["fbt_params"] = self.initialize_tk_vars(
            strvars=["device", "zu", "zl", "ri", "ro", "nz", "nr", "nu", "ilim", "r0", "b0", "fbtagcon", "bfct"], boolvars=["icsint"]
        )

        # Plasma scalar parameters
        state_tk["plasma_scalars"] = self.initialize_tk_vars(
            strvars=["Ip_time", "Ip_data", "Wk_time", "Wk_data", "qA_time", "qA_data", "ag_time", "ag_data"]
        )

        # Plasma profiles
        state_tk["pprime_editor"] = self.initialize_tk_vars(strvars=["exp", "exp2", "ped_center", "ped_width", "ped_height"])
        state_tk["ffprime1_editor"] = self.initialize_tk_vars(
            strvars=["exp", "core_val", "flat_val", "ped_center", "ped_width", "ped_height"]
        )
        state_tk["ffprime2_editor"] = self.initialize_tk_vars(strvars=["psin", "ffprime2_expression"])
        state_tk["pprime_profiles"] = self.initialize_tk_vars(
            strvars=["time", "exp", "exp2", "ped_center", "ped_width", "ped_height"]
        )
        state_tk["ffprime1_profiles"] = self.initialize_tk_vars(
            strvars=["time", "exp", "core_val", "flat_val", "ped_center", "ped_width", "ped_height"]
        )

        # Shape editor
        state_tk["shape_params"] = self.initialize_tk_vars(
            strvars=["Zup", "Zlo", "Rout", "Rin", "triu", "tril", "sqou", "sqol", "sqiu", "sqil", "c_xplo", "c_xpup"]
        )
        state_tk["control_segments"] = self.initialize_tk_vars(
            strvars=["Num_segs", "Seg_length", "Ellipse_r0", "Ellipse_z0", "Ellipse_a", "Ellipse_b"]
        )
        self.manual_pt_names = (
            [f"control_pt_ref{i + 1}" for i in range(3)] + [f"x_pt{i + 1}" for i in range(4)] + [f"strike_pt{i + 1}" for i in range(8)]
        )
        state_tk["manual_points"] = self.initialize_tk_vars(
            strvars=[f"{prefix}_{name}" for prefix in ["R", "Z"] for name in self.manual_pt_names]
        )
        state_tk["shape_plot_options"] = self.initialize_tk_vars(
            boolvars=["plot_segments", "label_control_points", "label_cp_ref_pts", "label_xpts", "label_strike_points"]
        )

        # Multishape evolution
        state_tk["multishape"] = self.initialize_tk_vars(strvars=["shape_time", "multishape_dir", "shape_filenames"])
        state_tk["multishape"]["multishape_dir"].set("<shape dir>")
        state_tk["multishape"]["shape_filenames"].set("<shape files>")

        # Optimization signals
        # --------------------

        # control_point flux error
        for i in range(1, 4):
            state_tk[f"rel_flux{i}"] = self.initialize_tk_vars(
                boolvars=[
                    "use_in_optimization",
                    "r_is_manual_mode",
                    "z_is_manual_mode",
                    "r_multiref_is_manual_mode",
                    "z_multiref_is_manual_mode",
                    "wt_multiref_is_manual_mode",
                ],
                strvars=[
                    "r_channels",
                    "z_channels",
                    "r_multiref_channels",
                    "z_multiref_channels",
                    "r_time",
                    "r_data",
                    "z_time",
                    "z_data",
                    "r_multiref_time",
                    "r_multiref_data",
                    "z_multiref_time",
                    "z_multiref_data",
                    "wt_multiref_time",
                    "wt_multiref_data",
                    "wt_multiref_channels",
                    "targ_time",
                    "targ_time_val",
                    "targ_channel_val",
                    "wt_time",
                    "wt_time_val",
                    "wt_channel_val",
                    "dwt_time",
                    "dwt_time_val",
                    "dwt_channel_val",
                    "d2wt_time",
                    "d2wt_time_val",
                    "d2wt_channel_val",
                ],
            )

        for x in ["flux_abs_avg", "vac_flux_abs_avg", "Br", "Bz", "Br_vac", "Bz_vac"]:
            state_tk[x] = self.initialize_tk_vars(
                boolvars=["r_is_manual_mode", "z_is_manual_mode","use_in_optimization"],
                strvars=[
                    "r_channels",
                    "z_channels",
                    "r_time",
                    "r_data",
                    "z_time",
                    "z_data",
                    "targ_time",
                    "targ_time_val",
                    "targ_channel_val",
                    "wt_time",
                    "wt_time_val",
                    "wt_channel_val",
                    "dwt_time",
                    "dwt_time_val",
                    "dwt_channel_val",
                    "d2wt_time",
                    "d2wt_time_val",
                    "d2wt_channel_val",
                ],
            )

        for x in ["voltage", "current"]:
            state_tk[x] = self.initialize_tk_vars(
                boolvars=["use_in_optimization"],
                strvars=[
                    f"{coil}_{suffix}"
                    for coil in self.coils
                    for suffix in [
                        "targ_time",
                        "targ_data",
                        "constrain_time",
                        "constrain_data",
                        "wt_time",
                        "wt_data",
                        "dwt_time",
                        "dwt_data",
                        "d2wt_time",
                        "d2wt_data",
                    ]
                ]
            )

        combo_matrix_vars = [f"combo_{i}_{j}" for i in range(self.ncombos) for j in range(self.ncoils)]
        state_tk["current_combo_matrix"] = self.initialize_tk_vars(strvars=combo_matrix_vars)
        state_tk["current_combos"] = self.initialize_tk_vars(
            boolvars=["use_in_optimization"],
            strvars=[
                "targ_time",
                "targ_time_val",
                "targ_channel_val",
                "wt_time",
                "wt_time_val",
                "wt_channel_val",
                "dwt_time",
                "dwt_time_val",
                "dwt_channel_val",
                "d2wt_time",
                "d2wt_time_val",
                "d2wt_channel_val",
            ]
        )

        return state_tk

    @staticmethod
    def initialize_tk_vars(strvars=[], intvars=[], boolvars=[]):
        tk_vars = {}
        for varname in strvars:
            tk_vars[varname] = tk.StringVar()
        for varname in intvars:
            tk_vars[varname] = tk.IntVar()
        for varname in boolvars:
            tk_vars[varname] = tk.BooleanVar()
        return tk_vars

    def initialize_scrollable_frame(self):
        """This function written with AI help"""
        # Create a main frame to hold everything
        self.main_frame = tk.Frame(self.root, bg=self.root.cget("bg"), highlightthickness=0)
        self.main_frame.pack(fill="both", expand=True)

        # Create a canvas inside the main frame
        self.canvas = tk.Canvas(self.main_frame, highlightthickness=0, bg=self.root.cget("bg"))
        self.canvas.pack(side="left", fill="both", expand=True)

        # Add vertical scrollbar on the left
        self.v_scrollbar = tk.Scrollbar(self.main_frame, orient="vertical", command=self.canvas.yview)
        self.v_scrollbar.pack(side="left", fill="y")

        # Configure canvas to use vertical scrollbar
        self.canvas.configure(yscrollcommand=self.v_scrollbar.set)

        # Create a frame inside the canvas to hold all widgets
        self.scrollable_frame = tk.Frame(self.canvas, bg=self.root.cget("bg"), highlightthickness=0)
        self.scrollable_frame_id = self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        # Update scrollregion when the size of the frame changes
        def _on_frame_configure(event):
            self.canvas.configure(scrollregion=self.canvas.bbox("all"))

        self.scrollable_frame.bind("<Configure>", _on_frame_configure)

        # Allow scrolling with mouse wheel (cross-platform)
        def _on_mousewheel(event):
            # Windows and MacOS
            if hasattr(event, "delta"):
                if event.delta > 0:
                    self.canvas.yview_scroll(-1, "units")
                elif event.delta < 0:
                    self.canvas.yview_scroll(1, "units")
            # Linux (event.num is used)
            elif hasattr(event, "num"):
                if event.num == 4:
                    self.canvas.yview_scroll(-1, "units")
                elif event.num == 5:
                    self.canvas.yview_scroll(1, "units")

        # Bind mousewheel events for Windows, Mac, and Linux
        self.canvas.bind("<Enter>", lambda e: self.canvas.focus_set())
        # Windows
        self.canvas.bind_all("<MouseWheel>", _on_mousewheel)
        # MacOS (event.delta is smaller, but handled above)
        # Linux
        self.canvas.bind_all("<Button-4>", _on_mousewheel)
        self.canvas.bind_all("<Button-5>", _on_mousewheel)

        # Increase the width of the vertical scrollbar
        self.v_scrollbar.config(width=20)

        # All subsequent widgets should be added to self.scrollable_frame instead of self.root
        self.root_widget = self.scrollable_frame

def get_vars_format_type(vars_dict: dict) -> dict:
    """Convert tkinter variable dict to regular dict with proper data types."""
    vars_formatted = {}
    for k, v in vars_dict.items():
        val = v.get()
        if isinstance(v, tk.StringVar):
            if len(val) == 0:
                val = None
            else:
                try:
                    val = float(val)  # val is numeric scalar
                except ValueError:
                    try:  # noqa: SIM105
                        # val represents a list/array
                        val = np.array(eval(val)).astype(float)  # noqa: S307
                    except (NameError, ValueError, SyntaxError):
                        pass  # val is non-numeric string, keep as is

        elif isinstance(v, tk.IntVar):
            val = int(val)
        elif isinstance(v, tk.BooleanVar):
            val = bool(val)
        else:
            raise TypeError(f"Error parsing {k} of type {type(v)}")
        vars_formatted[k] = val

    return vars_formatted


def main():
    parser = argparse.ArgumentParser(description="GSPulse GUI")
    parser.add_argument("--tokamak", "-t", type=str, default="sparc", help="Tokamak name (default: sparc)")
    args = parser.parse_args()

    app = GspulseApp(tokamak=args.tokamak)
    app.root.protocol("WM_DELETE_WINDOW", app.root.quit)  # Ensure quit on window close
    app.root.mainloop()


if __name__ == "__main__":
    main()
