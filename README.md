# atomizationFoam

## Features

- Volume-of-Fluid (VoF) solver for incompressible flows coupled with Lagrangian library
- Conversion of small VoF elements to Lagrangian parcels to reduce computational cost
- Support for
  - Particle-particle interaction (collision, coalescence) and secondary breakup
  - Adaptive mesh refinement
  - Fully parallelized
- Based on interIsoFoam in OpenFOAM v1912

## Compilation

- First, compile the library libAtomization in the src folder.
- Afterwards, compile the solver atomizationFoam in the application folder.

## Reference

Heinrich, M. and Schwarze, R. (2020).
3D-coupling of Volume-of-Fluid and Lagrangian particle tracking
for spray atomization simulation in OpenFOAM
SoftwareX, 11 (2020)
DOI 10.1016/j.softx.2020.100483
