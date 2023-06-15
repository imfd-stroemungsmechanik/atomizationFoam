# atomizationFoam

## Features

- Volume-of-Fluid (VoF) solver for incompressible flows coupled with Lagrangian library
- Conversion of small VoF elements to Lagrangian parcels to reduce computational cost
- Support for
  - Particle-particle interaction (collision, coalescence) and secondary breakup
  - Adaptive mesh refinement
  - Fully parallelized
- Based on interIsoFoam in [OpenFOAM v2212](https://www.openfoam.com).

## Compilation

 1. Clone the `atomizationFoam` gitHub repository:
```
git clone https://github.com/imfd-stroemungsmechanik/atomizationFoam.git
```
 2. Use the `Allwmake` script to compile the `atomizationFoam` solver and its corresponding library:
```
./Allwmake
```

## Reference

The source code has been published in the following open-access research article:
```
Heinrich, M. and Schwarze, R. (2020).
3D-coupling of Volume-of-Fluid and Lagrangian particle tracking for spray atomization simulation in OpenFOAM
SoftwareX, 11 (2020)
DOI 10.1016/j.softx.2020.100483
```
