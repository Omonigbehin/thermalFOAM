# thermalFOAM

`thermalFOAM` is a standalone OpenFOAM application for coupled heat transfer, phase change (apparent heat capacity with a Gaussian SFCC), and thaw-driven erosion (cell removal) for coastal permafrost / thermo-abrasion applications. It is built on the `finiteVolume`, `solidThermo`, and `dynamicFvMesh` libraries.

thermalFOAM refers to the solver/project; the compiled executable is `thermalFoam`.

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
```

## How to cite / archived version

If you use thermalFOAM in your work, please cite the archived release:

> Omonigbehin, O., Stolle, J., Kurylyk, B., & Guimond, J. (2026).
> *thermalFOAM (v1.0): A parallel OpenFOAM solver for simulating ablative erosion of permafrost bluffs.*
> Zenodo. https://doi.org/10.5281/zenodo.19684137

The `v1.0.0` tag corresponds to the version archived on Zenodo and cited in the accompanying *Coastal Engineering* manuscript.

## License

GPL-3.0. See [LICENSE](LICENSE) for details.
