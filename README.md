*Note: This is the public version of the GSPulse repository, containing the core GSPulse algorithms without the MEQ submodule and only minimal examples. Additional examples pulses (including other tokamaks, scenarios, integration with TokSys, integration with MEQ) exist, even if not available here. Please be patient while more examples are created/approved for public release. Please feel free to reach out to Josiah Wai, jwai@cfs.energy for help in setting up your pulses or for submitting issues.*

# Grad-Shafranov-Pulse (GSPulse)

**Description:** GSPulse helps you design a time-dependent sequence of Grad-Shafranov plasma equilibria for a plasma pulse. Equilibria satisfy both force balance and circuit dynamics. The algorithm is described in detail in an arxiv [preprint](https://arxiv.org/abs/2306.13163) (note: updated description of the algorithm available in the `docs` folder), with peer-reviewed journal publication forthcoming. Please cite as appropriate if you use the code for any published work. 

<p align="center">
  <img src="https://github.com/cfs-energy-internal/GSPulse/assets/137820499/f1800ff1-8f17-4100-84c0-797c60775327" width="600" />
  <br />
   <em>Time-dependent equilibrium trajectories.</em>  
</p>


## Installation without MEQ
GSPulse has 2 different options for performing the plasma picard iteration: a built-in option, or one that relies on the FBT/MEQ code. This is the installation option if access to MEQ is not available or required. The FBT/MEQ option has advantages in robustness and performance, but additional geometry and solver specifications must be provided. 

`git clone https://github.com/cfs-energy-internal/GSPulse.git`


## Installation with MEQ
For access to MEQ, see the information page [here](https://www.epfl.ch/research/domains/swiss-plasma-center/meq/). First, clone the GSPulse repository including the MEQ submodule.  

`git clone --recurse-submodules https://github.com/cfs-energy-internal/GSPulse.git`

Compile the MEQ submodule. 

```
cd MEQ
make CC=gcc MATPATH=<path_to_matlab_root_directory>
```

On a MacOS with ARM processesor, one needs to switch the processor environment first
```
cd MEQ
env /usr/bin/arch -x86_64 /bin/zsh --login
make clean
make CC=gcc USE_OPENMP=no MATPATH=<path_to_matlab_root_directory>
```

If you run into issues, there is some more information on compiling MEQ [here](https://github.com/cfs-energy-internal/MEQ/blob/b44f2f091f9b1f18f90d54c98695a1ab15e1883c/README.md)

## Run a test case

NSTX-U pulse 100 and SPARC pulse 100 are good example test cases to run. NSTX-U pulse 100 does not use MEQ (i.e. `gsp_inputs.settings.picard_algo = 'gspulse'`), whereas SPARC pulse 100 requires the FBT code within MEQ (i.e. `gsp_inputs.settings.picard_algo = 'fbt'`)

```
# does not require MEQ submodule
# designs a 1-second lower single null NSTX-U pulse
startup_gspulse.m
[soln, gsp_inputs] = run_pulse('nstxu', 100);
```

```
# requires MEQ submodule to be compiled
# designs one cycle of a strikesweep for 8.7MA sparc equilibrium
startup_gspulse.m
[soln, gsp_inputs] = run_pulse('sparc', 100);
```


All of the equilibria and coil current and voltage trajectories were saved to the `soln` struct. Running the script should have produced some plots of the coil currents and voltages, as well as an interactive GUI plot of the equilibria. 

<p align="center">
  <img src="https://github.com/jwai-cfs/GSPulseDesign/assets/137820499/c3c45ba4-d66f-432a-98e9-6fe1be565b6b" width="300" />
  <br />
   <em>Sample equilibrium from SPARC.</em>  
</p>

## Designing your own equilibria trajectories

All of the pulse-specific configuration files are stored in the `tokamaks/<tokamak>/pulse<#>` folder. To design your own pulse, copy one of the other pulses and modify the contents of each of these scripts. Running the run_pulse.m just reads from each of these config scripts sequentially and then calls the equilibrium solver. 

```
settings         = define_general_settings();
tok              = define_tok();
shapes           = define_shapes();
plasma_params    = define_plasma_params();
targs            = define_optimization_targets();
weights          = define_optimization_weights();
init             = define_init();
```


Here's a look at each of these scripts. 

#### 1. define_general_settings.m

This script defines top-level solver parameters such as the time basis on which to solve equilibria and the number of solver iterations to perform. The code comments within the script itself describe each of these settings and their effects in more detail. 


#### 2. define_tok.m
This script defines geometry data for the tokamak, reading from stored device data. This is where grid definitions are defined. 


#### 3. define_shapes.m
The equilibrium shape targets are created using the shape GUI, and then saved to json files. This script defines which shape json files to use as target shapes, and also the time basis that these correspond to, via this line of the script:

```
shapes = shapefiles2shapes(shape_filenames, time);
```

For example if `time = [0 1]` and `shape_filenames = {'shapeA.json', 'shapeB.json'}`, the target shapes will be a linear interpolation between shapeA and shapeB over a 1 second period from 0 to 1. 


Creating shapes is done by using the shape GUI. To create and save a new shape

```
within a terminal, cd GSPulse/shape_gui
python shape_gui.py
```

This will open up the shape GUI. Edit the shaping parameters as desired, and then save the shape using the save button. 
<p align="center">
  <img src="https://github.com/cfs-energy-internal/GSPulse/assets/137820499/22ef9941-8e3f-4335-a663-9dfac06879e3" width="900" />   
  <br />
  <em>Shape GUI.</em>
</p>

#### 4. define_plasma_params.m

This script defines the trajectory for plasma parameters such as plasma current (Ip), thermal energy (wmhd), and on-axis safety factor (qA). It is also where time-dependent profile basis functions for P' and TT' can be specified. In general, each signal does not have to use the same time basis, whenever time bases differ they will be resolved by linear interpolation. These plasma parameters will be enforced by the equilibrium solver (the specific combination of parameters that is enforced is specified by L.agconc in define_tok.m). 

#### 5. define_optimization_targs.m

This script defines values for the isoflux targets used by the optimization algorithm, such as flux error values and coil current targets. Often, we want to control the isoflux error to be zero (i.e. specify that the boundary flux surface passes exactly through that shape control point), so many of these are set to zeros of the appropriate dimension. However, one could also specify offset flux values, to indicate that a point lies a certain distance away from the plasma boundary. 

#### 6. define_optimization_weights.m

This script defines the weights used by the optimizer. Be aware that good weight choices may span orders of magnitude, because the controlled physical parameters are not normalized. The weight trajectories are time dependent which is important for scenario design. For example, if designing a scenario a where the plasma starts limited and then diverts, this might be a good choice for the weight on the x-point flux gradient, since it starts at zero and then increases (i.e. turning on the x-point). 

<p align="center">
  <img src="https://github.com/cfs-energy-internal/GSPulse/assets/137820499/75f56de5-0a71-4ebb-803e-b56b5db9dd5c" width="650" />   
  <br />
  <em>Time-dependent X-point weight.</em>
</p>



#### 7. define_init.m

This script defines the initial condition for the solver, that is, the initial coil currents and voltages. Use nan to specify no constraint on the initial condition, in which case the optimizer will solve for the best initial condition. Values that are not nan will be enforced explicitly. 
