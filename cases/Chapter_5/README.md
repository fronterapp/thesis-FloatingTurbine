This file aims to describe the process of compiling a version of waves2Foam compatible with dynamic meshes: **waveDyMFoam**

# How OpenFOAM libraries work
OpenFOAM uses its own command,_ ./Allwmake_, to  compile libraries. To do so, in a folder called _make_, there must be two files:
1. _files_. Contains all the .C files from the library. They will be dynamically linked to a .so file. You specify if the output of the compilation should be an executable (like wavesFoam) or a library (like turbinesFoam).
2. _include_. Here you indicate all the executables and libraries upon which your library or executable also depend. 


# Dynamic mesh with waves2FOAM
waves2FOAM is compatible with dynamic mesh. However, it is not implemented by default. To do so, follow the steps below:
1. Create a folder in _\applications\utilities\waves2Foam\applications\solvers\solvers2012_PLUS_ with the name of the solver that we want to create, e.g. **waveDyMFoam**.
2. Copy the contents of _OpenFOAM\OpenFOAM-v2012\applications\solvers\multiphase\interFoam_ in the folder created above. If you won't be needing the _interMixingFoam_ and _overInterDyMFoam_ folders, you can delete them. 
3. Change the name of the **interFoam.C** to the solver name, e.g. **interDyMFoam.C**.
4. With a file comparison tool, copy to **interDyMFoam.C** the lines that are missing from _\applications\utilities\waves2Foam\applications\solvers\solvers2012_PLUS\waveFoam\waveFoam.C_. 
5. Modify the **createField.H** file:

```cpp
//#include "readGravitationalAcceleration.H"
//#include "readhRef.H"
//#include "gh.H"
volScalarField gh("gh", g & (mesh.C() - referencePoint));
surfaceScalarField ghf("ghf", g & (mesh.Cf() - referencePoint));
```
At the end of the file, add:
```cpp
relaxationZone relaxing(mesh, U, alpha1);
```

6. Go to *Make/files*, should look like this:
```bash
waveDyMFoam.C

EXE = $(FOAM_APPBIN)/waveDyMFoam
```

7. Now go to *Make/options*, should look like this:
```bash
EXE_INC = \
    -I../VoF \
    -I$(FOAM_SOLVERS)/multiphase/VoF \
    -I$(LIB_SRC)/phaseSystemModels/twoPhaseInter/incompressibleInterPhaseTransportModel/lnInclude \
    -I$(LIB_SRC)/phaseSystemModels/twoPhaseInter/VoFphaseIncompressibleTurbulenceModels/lnInclude \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/sampling/lnInclude \
    -I$(LIB_SRC)/dynamicFvMesh/lnInclude \
    -I$(LIB_SRC)/transportModels \
    -I$(LIB_SRC)/transportModels/incompressible/lnInclude \
    -I$(LIB_SRC)/transportModels/interfaceProperties/lnInclude \
    -I$(LIB_SRC)/transportModels/twoPhaseMixture/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/incompressible/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/phaseIncompressible/lnInclude \
    -I$(LIB_SRC)/transportModels/immiscibleIncompressibleTwoPhaseMixture/lnInclude \
    -DOFVERSION=2012 \
    -DEXTBRANCH=0 \
    -DOFPLUSBRANCH=1 \
    -DXVERSION=$(WAVES_XVERSION) \
    -I$(WAVES_SRC)/waves2Foam/lnInclude \
    -I$(WAVES_SRC)/waves2FoamSamplingNew/lnInclude \
    -I$(WAVES_GSL_INCLUDE)

EXE_LIBS = \
    -lfiniteVolume \
    -lfvOptions \
    -lmeshTools \
    -lsampling \
    -ldynamicFvMesh \
    -lincompressibleTransportModels \
    -linterfaceProperties \
    -limmiscibleIncompressibleTwoPhaseMixture \
    -lturbulenceModels \
    -lincompressibleTurbulenceModels \
    -lwaveModels \
    -lVoFphaseTurbulentTransportModels \
    -lincompressibleInterPhaseTransportModels\
    -L$(WAVES_LIBBIN) \
    -lwaves2Foam \
    -lwaves2FoamSampling \
    -L$(WAVES_GSL_LIB) \
    -lgsl \
    -lgslcblas
```

7. Source the waves2Foam bashrc before compiling. To do that, go with the terminal to _\applications\utilities\waves2Foam\bin_ and type:
```bash
$ . bashrc
```

8. Now go back to _\applications\utilities\waves2Foam\applications\solvers\solvers2012_PLUS\waveDyMFoam_ and type
```bash
$ wmake
```

9. To install in another computer, just copy the _waveDyMFoam_ folder in the appropriate directory and repeat steps 7 and 8.

