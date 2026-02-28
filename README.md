# thermalFOAM

`thermalFOAM` is an OpenFOAM-based solver for coupled heat transfer, phase change (enthalpy-porosity), and thaw-driven erosion (cell removal) for coastal permafrost / thermo-abrasion applications.

## Repository structure
- `thermalFoam/` : solver source code
- `tutorials/`   : example cases
- `Allwmake`     : build script

## Requirements
- OpenFOAM v2506 (or compatible)
- A C++ compiler supported by your OpenFOAM installation

## Build
From the repo root:
```bash
./Allwmake
