floatingSixDoFRigidBodyMotion
============
![OpenFOAM v2012](https://img.shields.io/badge/OpenFOAM-v2012-brightgreen.svg)

Library based on OpenFOAM's `sixDoFRigidBodyMotion` (from v2012, you will find it in `$FOAM_SRC/sixDoFRigidBodyMotion`). 

Uploaded for educational purposes only. This is not an official release. 

Features
-----
- Coupling to an AL-turbine:
	- Turbine loads are applied onto the rigid body via the `turbineAL` restraint.
	- Requires the `floatingTurbinesFoam` library.
- New restraints:
	- `constantLoad`. Apply a constant force and torque.
	- `gyroscopicMoment`. Accounts for the gyroscopic effect imparted by a rotating body which 
    is attatched to the six-DoF rigid body, whose rotation is not modelled.
	- `mooringLine`. Quasi-steady catenary mooring model from [waves2Foam](https://www.researchgate.net/publication/319160515_waves2Foam_Manual).
	- `turbineAL`. This model accounts for the loads from a turbine from `floatingTurbinesFoam`. The `coupleLoads` 
    flag in _fvOptions_ must be set to true.

Modifications
-----
Apart from the above-mentioned features, the following has been modified:
- The `g` value defined in `constant/dynamicMeshDict` has priority over the value defined in `constant/g`.
- Added three new variables to the _status_ report: linear acceleration, angular momentum, and torque. 

Installation
-----
1. Create an empty folder called `src` within `$WM_PROJECT_USER_DIR` (only if it doesn't exist already).
2. Copy the `floatingSixDoFRigidBodyMotion` folder to `$WM_PROJECT_USER_DIR/src`.
3. Run `wmake`.
4. Include the following line in your `dynamicMeshDict` file:

```cpp
motionSolverLibs ("libfloatingSixDoFRigidBodyMotion.so");
```
Usage
-----
Refer to the thesis document (sections 7.1.2 and B.2) and tutorial cases from the `cases/Tests/RigidBody_ALM` folder.

Acknowledgements
----------------
__OpenFOAM__ is free, open source software for computational fluid dynamics (CFD),
developed primarily by [CFD Direct](http://cfd.direct), on behalf of the
[OpenFOAM](http://openfoam.org) Foundation.
