# Recent news

### Version 2.0.0 release
Version 2.0.0 introduces the GUI-based config editor as the standard format for defining a GSPulse run. Support for previous APIs is deprecated.

### GSPulse article published
GSPulse was recently published in *Nuclear Fusion* [here](https://iopscience.iop.org/article/10.1088/1741-4326/ae203c). Many thanks to all co-authors and developers!

# Installation:
To install GSPulse, run: 

```
bash scripts/install.sh
```

This script will fetch all dependent submodule repositories and compile the MEQ and SCS libraries. This script has not been extensively tested on all operating systems. Please consider contributing if you spot needed fixes. 

## Running a pulse

Launch the gspulse config editor GUI:

```
cd src/python
uv run gspulse_gui.py --tokamak anamak
```
Configure your pulse as desired (see below for more detail on each setting). To run the demo pulse: 

```
1. Click the "Load" button and load "demo_pulse_gui_state.json"
2. Click the "Run GSPulse" button    # Runs GSPulse via Octave. Outputs are saved to .mat file. 
3. (optional) visualise the saved .mat file outputs within Matlab (TODO: visualize via GUI): 
     - Open Matlab and cd to src/matlab. Within Matlab window: 
     - startup_gspulse('anamak')
     - plot_saved_run('anamak', 'demo_pulse')
```
 <img width="820" height="580" alt="image" src="https://github.com/user-attachments/assets/04cbf205-63c1-47d5-b04d-a55a9304c723"/>   

An alternative workflow is: 
```
1. Open Matlab and cd to src/matlab. Within Matlab window: 
2. startup_gspulse('anamak')
3. run_pulse('anamak', 'demo_pulse')  # Run and visualize GSPulse via Matlab
```

# Overview of GSPulse

GSPulse is a tool for designing time-dependent sequences of Grad-Shafranov plasma equilibria for a tokamak pulse. The equilibria satisfy both Grad-Shafranov force balance and axisymmetric circuit dynamics. 

<br />
<p>
 <img width="820" height="580" alt="image" src="https://github.com/user-attachments/assets/67783186-7f69-46a4-916a-3ff5a6bb55f2"/>   
</p>
<p>
     <em>GSPulse solves for time-dependent sequences of equilibria by performing a global optimization over the trajectory.</em>
</p>

<br />

To design a pulse, a user would configure a set of inputs that include: 

- Target shape evolution, including features such as x-points and strike points. 
- The power supply current and voltage limits for each of the tokamak shaping coils. 
- Plasma properties such as the evolution of stored thermal energy, and Grad-Shafranov profile basis function shapes. 
- Tokamak geometry data including mutual inductance matrices. 
- A set of optimizer weights used to tune the trajectory. 

Then GSPulse assembles and solves an optimization problem and converges towards an equilibrium trajectory. The output of GSPulse is the sequence of time-dependent equilibria as well as the coil currents and voltages required to achieve these equilibria. 

More detail on GSPulse is available in the publication [here](https://iopscience.iop.org/article/10.1088/1741-4326/ae203c). Please cite this article if you use GSPulse for any published work. 


# Description of configuration inputs

## Settings tab

| Variable| Type | Description |
|---|---|---|
| pulse_id | int or str | identifier for the pulse |
| tokamak | str | tokamak name |
| interval_t1 | array-like | global time basis for the solution, which can be specified in multiple intervals
| interval_t2 | " | "
| interval_t3 | " | " 
| interval_t4 | " | "
| picard_algo | str = fbt| currently only "fbt" supported
| specify_psibry_mode | str = direct | currently only "direct" supported
| qpsolver | str = {quadprog, scs, python_cvxopt} | Which QP solver program to use. If running in Matlab, all options are supported although quadprog requires the Matlab optimization toolbox. If running via Python through the GUI (run_pulse button) then python_cvxopt is likely to work best; scs may work but not robustly tested. 
| qpsolver_tol | float | tolerance for satisfying quadratic program
| spline_basis_ratio | float range (0-1) | How much to compress the global time basis during the solve using a spline basis. Only applies if use_spline_basis=true, and if the global time interval_t are uniform. 
| initial_vess_currents | float | initial vessel current values at first time
| fbt_Ip_threshold | float | When Ip is above this value, the plasma current density distribution is computed according to Grad-Shafranov force balance. When Ip is below this value, the distribution is computed via interpolation to the nearest Grad-Shafranov time. This parameter is used to improve convergence when simulating across the vacuum-plasma transition when Ip is low.  
| tol_ic_diff  | float | Tolerance for determining convergence. If >0, the solution is determined converged if the norm of difference of currents between iterations is below this value for multiple successive iterations. 
| Iy_relax_factor | float (0-1), recommend 1 unless having convergence difficulty | Relaxation factor on updating the plasma current density distribution across iterations. 0 = do not update the solution at all, 1 = full update of the solution. 
| Iy_smooth_dt | float, recommend 0 unless having convergence difficulty | Time constant for applying time-based smoothing to plasma current density distribution
| plotlevel | int = {0,1,2} | 0 = no plots, 1 = plot outputs, 2 = plot inputs and outputs. Only applicable when running via Matlab. 
| niter | int | Max number of Grad-Shafranov iterations to perform
| verbose | int = {0,1,2} | 0 = minimal, 1 = info, 2 = debug
| do_final_boundary_trace | bool | whether to trace the plasma boundary during post processing
| use_spline_basis | bool | whether to use a spline basis to compress the number of optimization moves
| calc_strike_pts | bool | Only supported for tokamak=sparc, whether to compute strike point locations
| inject_model_offset | bool | Currently not well-supported. Use false. 
| calc_post_prc_ext  | bool | Only supported for tokamak=sparc, whether to compute some additional derived shape parameters in post-processing
| boundary_post_run | bool | Only supported for tokamak=sparc, whether to compute boundary contour integrals. 
| maxI (per coil) | float [kA] or np.inf| Maximum allowed coil current
| minI (per coil) | float [kA] or -np.inf | Minimum allowed coil current
| maxV (per coil) | float [kV] or np.inf | Maximum allowed power supply voltage
| minV (per coil) | float [kV] or -np.inf | Minimum allowed power supply voltage
| initI (per coil) | float [kA] | Enforced initial coil currents at t=t0. Leave blank to specify no constraint. 
| initV (per coil) | float [kV] | Enforced initial power supply voltages at t=t0. Leave blank to specify no constraint.
| device | int or str | identifier for device version
| zu | float | Follows MEQ/FBT definition, upper z limit for grid. 
| zl | float | Follows MEQ/FBT definition, lower z limit for grid. 
| ri | float | Follows MEQ/FBT definition, inner r limit for grid. 
| ro | float | Follows MEQ/FBT definition, outer r limit for grid. 
| nz | int (must be power of 2, e.g. 32,64) | Follows MEQ/FBT definition, number of z grid points, final grid size is nz^2+1
| nr | int | Follows MEQ/FBT definition, number of r grid points, final number is nr+1
| nu | int | Follows MEQ/FBT definition, number of vessel eigenmodes to retain
| ilim | int (recommend 3) | Follows MEQ/FBT definition, method for computing limiter touch points
| r0 | float | radial position for defining vacuum toroidal field
| b0 | float | vacuum toroidal field
| fbtagcon | list of str | Follows MEQ/FBT definition, but use list of string syntax. E.g. ["Wk", "Ip", "qA"], or ["Wk", "Ip", "ag"]
| bfct | str = {bfabmex, bf3imex} | Follows MEQ/FBT definition. bfabmex=simple default parameterization for profiles. bf3imex = custom-defined basis functions. See also the profiles tab, which applies when bfct = bf3imex. 
| icsint | bool (recommend true) | Follows MEQ/FBT definition, method for computing touch points. 

## Plasma scalars tab
| Waveform | Description |
|---|---|
| Ip | Waveform for the plasma current in [MA] vs time in [s]. Only the waveforms that are specified in fbtagcon are used. For example, if fbtagcon=["Ip", "Wk", "qA"], those waveforms will be used and the "ag" waveform data is unused. 
| Wk | Follows MEQ/FBT definition, waveform for the plasma thermal energy in [MJ] vs time in [s]. 
| qA | Follows MEQ/FBT definition, waveform for the plasma on-axis q-value.
| ag | Follows MEQ/FBT definition, waveform for the profile basis function coefficient. 

## Profiles tab

The data in this tab is only used if bfct=bf3imex!

These inputs are used to define the profile basis functions in FBT. One basis function is used for pprime, and two basis functions are used for ffprime. 

For pprime and the first ffprime basis function, specify time-based arrays of coefficients. 

For pprime, the functional form is: 

$$P' = (1-\psi_N)^{exp1} * \psi_N^{exp2}$$


plus the addition of a pedestal defined by `ped_width`, `ped_height`, and `ped_center`.

For ffprime, the functional form is: 

$$FF' = (1-\psi_N)^{exp1}$$

with the core flattened to `core_val` where $\psi_N<$`flat value` and the addition of a pedestal defined by `ped_width`, `ped_height`, and `ped_center`.

The second basis function is constant across time and can be expressed as desired by user. 

## Shape Editor tab
The shaping inputs like `triangularity`, `squareness`, etc are used to parameterize and create shapes, as well as control segments used to define the actual isoflux shape control points. These shapes should then be saved to file. The time-dependent evolution is specified by selecting a set of shape files and corresponding time basis. 

## Optimization signals
These are the tuning parameters for the optimization, and are used to configure relative importance of aspects of the trajectory, such as smoothness of coil current evolution or importance of exactly achieving the desired shape. 

GSPulse uses isoflux-type signals within the optimization. Recall that GSPulse penalizes the error, derivative of error, and second derivative of error. For example, the coil currents signal is specified to include this data:


```
name        = 'current' 
target      = Timeseries(...) # target values
wt          = Timeseries(...) # error weights
dwt         = Timeseries(...) # error 1st-derivative weights
d2wt        = Timeseries(...) # error 2nd-derivative weights
```

The target, wt, dwt, and d2wt fields define the cost function associated with the coil currents. At each time step, the cost penalty is a quadratic function of the error and its time derivatives: 

```
e = Target - Actual   # coil current errors
cost = (e * wt * e)  +  (ė * dwt * ė) +  (ë * d2wt * ë)
```

Each optimization signal type must include waveforms for this data, as well as some signal-dependent additional specification. 

The error definitions for each signal are specified as follows:  

| Signal | Description | Math definition | Additional information
|---|---|---|---|
| Voltage | Power supply voltages | $e = Voltage$ | This signal supports the "constraint" option, where each coil can be constrained to the target value at specified times. 
| Current | Coil currents | $e = Target - CoilCurrents$ | This signal supports the "constraint" option, where each coil can be constrained to the target value at specified times. 
| Current combos | Linear combinations of coil currents | $e = Target - CoilCombinationsMatrix \times CoilCurrents$ | -- 
|rel_flux1,2,3| Weighted flux difference at multiple spatial locations | $e = Target - (Flux @ (R,Z) - WeightMultiref \times Flux @ (rMultiref, zMultiref))$ | This is a common signal for achieving a shaping objective. For example, if trying to create a double null plasma this can be configured for the flux at the control points minus an average of the flux at the x-points. 
| Flux absolute average | Average flux value at specified locations | $e = Target - Mean(Flux @ (r,z))$ | This is often used, for example, to specify the plasma surface voltage by specifying a ramp trajectory for the average flux at the plasma boundary. 
| Vac flux absolute average | Same as above, but only considering the vacuum contribution to the flux (not the plasma portion) | $e = Target - Mean(VacuumFlux @ (r,z))$ | -- 
| Br, Bz | Radial and vertical field | $e = Target - B_r @ (r,z),\\ e = Target - B_z @ (r,z)$ | Can be used to specify, for example, zero field at locations in order to create an x-point. 
| Br_vac, Bz_vac | Radial and vertical field only considering the vacuum contribution | $e = Target - VacuumB_r @ (r,z)\\ e = Target - VacuumB_z @ (r,z)$ | Useful for plasma breakdown scenarios, which require specified trajectories for the vacuum fields. 


# Adding another tokamak
Currently supported tokamaks in GSPulse V2 are SPARC, ARC, and the demo tokamak "Anamak". To add another tokamak, copy the `anamak` folder structure. The main body of work will be to define the `L` geometry data structure used by the MEQ/FBT codes under the hood. Reach out to the MEQ development team for help with this step if necessary (see https://gitlab.epfl.ch/spc/public/meq/meq). 

The GSPulse run function assumes that there is device-specific startup located in a `startup_<tokamak>.m` script and that the L geometry object can be loaded by calling `load_<tokamak>_L_fbt.m`. The GUI assumes that there is a `<tokamak>_gui_settings.json` file that contains some basic geometry information such as coil names and limiter definition for use by the GUI. After these steps are performed one just needs to develop and tune the pulse inputs using the GUI. 
