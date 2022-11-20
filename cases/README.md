Simulation cases
============

This folder contains all simulation setups from the thesis, divided by chapter:

- **Chapter 4**: Wave propagation (Stokes second-order) through a 2D wave flume.
- **Chapter 5**: Floating body:
	- Free heave-decay of a 2D-cylinder, based on [Ito's experiment](https://www-semanticscholar-org.tudelft.idm.oclc.org/paper/Study-of-the-transient-heave-oscillation-of-a-Ito/5274aa0d531b2b672c43d3a60f24b60e87eaebe9). 
	- Free and moored pitch-decay of a 3D-buoy, based on the work by [Paredes et al.](https://www-sciencedirect-com.tudelft.idm.oclc.org/science/article/pii/S2214166916300212)
- **Chapter 6**: ALM turbine with prescribed motion, based on the [OC6 Phase III campaign](https://wes.copernicus.org/preprints/wes-2022-74/).
- **Chapter 7**: Semi-submersible floating wind turbine under combined wind-wave conditions.
	- Coupled (turbine + floater + moorings).
	- Floater-only (prescribed thrust).
	- Turbine-only (prescribed motion).

The folder `Tests` contains basic _tutorials_ to test the new features, namely:
- Prescribed-motion turbine, `prescribed_ALM`
- Turbine coupling with rigid-body, `rigidBody_ALM`

Requisites
-----
All cases have been tested for OpenFOAM v2012 only.

- **Chapter 4**: requires `waves2Foam`.
- **Chapter 5**: requires `wavesDyMFoam`. For the moored-buoy case, also `floatingSixDoFRigidBodyMotion`.
- **Chapter 6**: requires `floatingTurbinesFoam`.
- **Chapter 7**: requires `wavesDyMFoam`, `floatingSixDoFRigidBodyMotion` and `floatingTurbinesFoam`.

Usage
-----
The simulations are run sequentially: 
1. Run pre-processing: `./preRun`
2. Run parallel decomposition: `./parRun`
3. Run the solver: `./solverRun`
4. Run post-processing: `./postRun`

Some cases might run with `./allRun` only.

The output from `./parRun` and `./solverRun` can be deleted with `./runClean`.

**Other tools**:
- For rigid-body simulations, extract the motion state with `./motionExtract.sh`
- For mooring line simulations, extract the loads with `./moorExtract.sh`
- For ALM simulations, extract the actuator lines' geometry with `./ALgeometry.sh`

**Note**: you might need to modify `./ALgeometry.sh` to specify the actuator lines' names.