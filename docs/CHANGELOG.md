# Change Log

All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [2.0.0] - 2026-01-07
### Added 
- Config editor GUI. All user settings can be specified through GUI. 
- install.sh script to facilitate easier install and build
- scs-matlab as a submodule, to remove license requirement for matlab optimization_toolbox
- API to run GSPulse directly from GUI

### Changed
- Most pulses deprecated. SPARC 113 and ARC 109 updated to V2. 
- CI test updated for new API

## [1.0.1] - 2026-01-07
### Added 
- ARC pulse
- minor bugfixes

## [1.0.0] - 2025-08-13
 
### Added
- Refactor configuration inputs to use the "optimization signals" object with a set of supported calculation types.  
- Support optimization with non-uniform time. 
- Introduce semantic versioning. 
- Add QP solver tolerance as a user setting. 
- Add the shape evolution tab to shape editor GUI. 
- Improve Python API, out-of-the-box usage of run_pulse.py. 
- Added SPARC Pulse 113 which includes vacuum phases and plasma startup. 
 
### Changed
- Input config API has changed and is not backwards compatible. See README.md for the new API description. 
- Most pulses were deprecated. These pulses can be updated to the v1 API on an as-needed basis.  
- Cost function quantities associated with time 1st and 2nd derivatives have been modified. Previously, the cost function was based on first and second-order differences. However now they represent true time derivatives. This means that associated weights will need to be re-tuned to give the same result (To preserve the scaling, first derivative weights should be multiplied by $\Delta t^2$, and second derivative weights by $\Delta t^4$ where $\Delta t$ is the time step size.)
- Shape editor GUI was refactored and some of the variables re-named. Shapes were updated to the new GUI format, but external shape files that are referenced by GSPulse will likely need to be re-created.

 
### Fixed
- Bugfix on the how the initial plasma current distribution was estimated. (This bug did not affect fully converged solutions, but did increase the number of iterations to convergence.)
 
## [0.0.0] - 2025-08-12
   
### Added
- Semantic versioning was added with v1.0.0. V0 represents the state of the GSPulse repository main branch at the time of merging V1. 
