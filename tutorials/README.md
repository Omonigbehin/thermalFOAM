# thermalFOAM Tutorials

## test07 — Laboratory bluff erosion (Test 7 validation case)

This tutorial reproduces Test 7 from Omonigbehin et al. (2025): wave-driven thaw
erosion of a frozen sandy-soil specimen in a wave flume.

---

### Physical setup

| Parameter | Value |
|-----------|-------|
| Domain | 20 cm (x, into bluff) × 27 cm (y, span) × 27 cm (z, height) |
| Mesh | 60 × 27 × 86 = 139,320 hex cells, uniform 3.3 mm spacing |
| Initial temperature | −8 °C at wave face → −18 °C at back (from thermocouple data) |
| SFCC model | Gaussian, width W = 0.15 K, residual fraction = 0.0183 |
| Erosion zone | z ∈ (0.11, 0.27] m (`bluffBase` → `bluffCrest`) |
| Wave boundary | Time-varying h and T_water from lookup tables; free surface elevation from `waveRecord_WF.dat` |
| Run time | 0 – 1250 s (adaptive dt, max 5 s) |

---

### Prerequisites

1. OpenFOAM v2506 (or compatible ESI/OpenCFD release) sourced in your shell:
   ```bash
   source /usr/lib/openfoam/openfoam2506/etc/bashrc
   ```
2. thermalFOAM compiled and on your `$PATH`:
   ```bash
   # From the repository root:
   ./Allwmake
   which thermalFoam   # should return $FOAM_USER_APPBIN/thermalFoam
   ```

---

### Running the tutorial

```bash
cd tutorials/test07
./Allrun
```

`Allrun` performs these steps automatically:

| Step | Tool | What it does |
|------|------|-------------|
| 1 | `blockMesh` | Builds the 3-zone hex mesh |
| 2 | `topoSet` | Creates the empty `f_empty` face set |
| 3 | `createPatch -overwrite` | Adds the initially-empty `erodedWall` patch |
| 4 | `checkMesh` | Validates mesh quality |
| 5 | `setFields` | Patches the non-uniform initial temperature profile (40 × 5 mm x-slabs) |
| 6 | `thermalFoam` | Runs the thermal + erosion solver |

To clean and restart from scratch:
```bash
./Allclean
./Allrun
```

---

### Key input files

#### `constant/transportProperties`
All material and erosion parameters.  Important entries:

```c++
SFCCModel       gaussian;      // freezing curve model
SFCC_W          0.15;          // Gaussian width [K]
thetar_frac     0.018275;      // residual unfrozen fraction
L               334e6;         // volumetric latent heat [J/m³]

erosionControls
{
    enableErosion     true;
    thawThreshold     0.99;    // remove cell when 99% thawed
    bluffCrest        0.27;    // top of erosion zone [m]
    bluffBase         0.11;    // bottom of erosion zone [m]
    ...
}
```

#### `constant/waveRecord_WF.dat`
Two-column ASCII file `(time [s]  eta [m])` at 0.02 s intervals, covering
0–1314 s.  This drives the wave free surface elevation seen by the `waveFacing`
boundary condition.

#### `system/controlDict`
```c++
endTime         1250;     // [s]
deltaT          0.1;      // initial step
adjustTimeStep  yes;
writeInterval   20;       // write every 20 s
```

#### `system/fvSolution` — Picard subdict
```c++
nCorrectors     100;
tolerance       1e-6;
flRelaxation    0.9;
maxDeltaT       5;        // [s]
minDeltaT       0.01;     // [s]
```

---

### Expected output

- Write directories appear every 20 s: `20/`, `40/`, …, `1240/`
- `postProcessing/erosionHistory.dat` — time-series of eroded volume, cell count,
  and niche depth, appended whenever cells are removed
- `postProcessing/horizontalProbes/` — T, liquidFraction, Kth, CthEff at 11 points
  along the x-axis at z = 0.135 m (wave impact zone mid-height)

`liquidFraction ≥ 0.99`.  Probe signals go to −1e+300 as their host cells
are removed.

---

### Post-processing in ParaView

```bash
paraview test07.foam &   # if you kept the .foam file, otherwise:
touch test07.foam && paraview test07.foam &
```

Useful visualisations:
- **Temperature slice** (z-normal at z = 0.135 m) — shows thaw front advancing
- **liquidFraction** volume rendering — highlights the mushy zone
- **erosionMask** — marks cells that have been removed (value = 1 before removal)
- **T vs time** from `postProcessing/horizontalProbes/0/T` — directly comparable
  to experimental thermocouple data

---

