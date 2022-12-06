CFD Simulation of a Floating Wind Turbine with OpenFOAM
============
![OpenFOAM v2012](https://img.shields.io/badge/OpenFOAM-v2012-brightgreen.svg)

![simulationFloatingTurbine](https://user-images.githubusercontent.com/104892099/202909133-f9e5fe98-97e3-451a-807a-7d1902d009d9.png)

This page contains all simulation cases, modified libraries and post-processing tools from the Master Thesis
__"CFD Simulation of a Floating Wind Turbine with OpenFOAM: an FSI approach based on the actuator line and relaxation zone methods"__
by Pere Frontera.

Uploaded for educational purposes only.

Contents
-----
- **cases**: contains the necessary OpenFOAM files to reproduce the simulation cases from the thesis.
- **libraries**: contains the two modified libraries required for the simulations from chapters 6 and 7.
- **other**: contains post-processing and job submission (PBS) scripts.

Usage
-----
Refer to the `README.md` from within the different folders.

Publications
------------
Pere Frontera. _"CFD Simulation of a Floating Wind Turbine with OpenFOAM: an FSI approach based on the actuator line and relaxation zone methods"_. MSc Thesis. TU Delft, 2022. An electronic version of this thesis is available at [http://repository.tudelft.nl](http://resolver.tudelft.nl/uuid:1911a06e-da04-4047-b51c-8a4b563d43d2).

![floatingGIF](https://user-images.githubusercontent.com/104892099/205907005-3ed8ec98-5732-43bb-b2a0-f4c8d7cabe65.gif)

Acknowledgements
----------------
__OpenFOAM__ is free, open source software for computational fluid dynamics (CFD),
developed primarily by [CFD Direct](http://cfd.direct), on behalf of the
[OpenFOAM](http://openfoam.org) Foundation.

__turbinesFoam__ is a library for simulating wind and marine hydrokinetic turbines
in OpenFOAM using the actuator line method, developed by [P. Bachant](https://github.com/turbinesFoam/turbinesFoam).

__waves2Foam__  is a toolbox used to generate and absorb free surface water waves based on the relaxation zone method, developed by [Niels G. Jacobsen](https://www.researchgate.net/publication/319160515_waves2Foam_Manual) and mantained by Deltares.
