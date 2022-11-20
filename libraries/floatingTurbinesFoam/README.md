floatingTurbinesFOAM
============
![OpenFOAM v2012](https://img.shields.io/badge/OpenFOAM-v2012-brightgreen.svg)

Library based on `turbinesFoam` by [P. Bachant](https://github.com/turbinesFoam/turbinesFoam) for the CFD Simulation of Floating Wind Turbines. 

Uploaded for educational purposes only. This is not an official release. 

Always cite the [original library](https://zenodo.org/record/3542301).

Features
-----
- Impose a prescribed motion (in any of the six DoFs) to `actuatorLineSource` and `axialFlowTurbineALSource`. 
	- So far, only multi-DoF harmonic motion is available.
	- To create new motions, modify the following functions from `actuatorLineSource.C` accordingly:
		- `prescribedMotionParameters`
		- `prescribedMotionInitialise`
- `axialFlowTurbineALSource` motion based on rigid-body dynamics.
	- Requires the `floatingSixDoFRigidBodyMotion` library.

Modifications
-----
Apart from the motion routines, some minor aspects were changed from the original `turbinesFoam`:
- **Blade pitch definition**. The changes are discussed [here](https://github.com/turbinesFoam/turbinesFoam/issues/350).
- **End-effects corrections**. In `axialFlowTurbineALSource::calcEndEffects()`, the rotor plane direction is defined as `-axis` instead of `freeStreamDirection`, since the rotor plane might move. 
- **Zero division**. If the turbine is moving upwind at exactly the freestream velocity, the relative velocity equals zero. Some modifications have been introduced to `actuatorLineElement.C` to avoid this situation. It affects:
	- The variables `relativeVelocity_`, `relativeVelocityGeom_` and `Re_`.
	- The function `inflowRefAngle()`.

Installation
-----
1. Copy the `floatingTurbinesFoam` folder to `$WM_PROJECT_USER_DIR`.
2. Run `./Allwmake`.
3. Include the following lines in your `controlDict` file:

```cpp
libs
(
    "libfloatingTurbinesFoam.so"
);
```
Usage
-----
Refer to the thesis document (sections 6.2.2 and 7.1.2) and tutorial cases from the `cases/Tests` folder.

Acknowledgements
----------------
__OpenFOAM__ is free, open source software for computational fluid dynamics (CFD),
developed primarily by [CFD Direct](http://cfd.direct), on behalf of the
[OpenFOAM](http://openfoam.org) Foundation.

__turbinesFoam__ is a library for simulating wind and marine hydrokinetic turbines
in OpenFOAM using the actuator line method, developed by [P. Bachant](https://github.com/turbinesFoam/turbinesFoam).
