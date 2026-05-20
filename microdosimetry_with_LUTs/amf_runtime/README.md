# AMF TOPAS Runtime

This folder supports TOPAS runs that use the Analytical Microdosimetric
Function (AMF) extension. The AMF extension reconstructs microdosimetric
quantities from fitted analytical probability-density functions rather than
explicitly simulating nanometer-scale track structure.

## References

1. Hartzell S, Parisi A, Sato T, Beltran CJ, Furutani KM. Extending TOPAS with
   an analytical microdosimetric function: application and benchmarking with
   nBio track structure simulations. Phys. Med. Biol. 70 (2025), 105010.
   doi: 10.1088/1361-6560/adcfec.
2. Sato T et al. Improvement of the hybrid approach between Monte Carlo
   simulation and analytical function for calculating microdosimetric
   probability densities in macroscopic matter. Phys. Med. Biol. 68 (2023),
   155005. doi: 10.1088/1361-6560/ace14c.

## What Is Included

- `staged_runs/`: generated TOPAS run folders for phase-space replay.
- `reference_outputs/`: manually reviewed outputs that may be used for checks.
- `../amf_extension/`: AMF TOPAS scorer source files.
- `../lookup_tables/AMF/tsed.dat`: fitted AMF parameter database.

Generated staged runs are disposable. Keep only reviewed reference outputs or
comparison runs that you still need.

## AMF Quantities

The extension supports three primary scoring quantities:

- `AMF_yD`: dose-weighted mean lineal energy.
- `AMF_yS`: saturation-corrected dose mean lineal energy.
- `AMFSpectra`: dose-weighted microdosimetric spectra, `yd(y)`.

The extension also includes `AMFDose`, used internally for consistent
dose-weighting behavior.

## Workflows

AMF can be used in two ways.

**Phase-space replay:** use this repository's `AMF-stage` or `AMF-run` wrapper
to replay an existing TOPAS `.phsp/.header` pair through a selected detector.
This is useful for comparing detectors or re-scoring saved phase spaces without
rerunning the full beamline.

**Full simulation:** put the AMF detector and scorer block directly in a TOPAS
input file with the original source and geometry, then run TOPAS normally. This
is the cleanest workflow when you want the AMF output scored during the original
transport.

## Required Files

For phase-space replay, AMF needs:

- a TOPAS executable compiled with the AMF extension,
- a phase-space pair with matching base names, for example
  `input/PhaseSpace_curved_33mm.phsp` and
  `input/PhaseSpace_curved_33mm.header`,
- `lookup_tables/AMF/tsed.dat`,
- one AMF quantity: `AMFSpectra`, `AMF_yD`, or `AMF_yS`.

The `.header` file is required by TOPAS phase-space replay. It is metadata for
TOPAS, not a source of truth for scored particle counts.

For full simulation, place `tsed.dat` in the directory where TOPAS is launched,
or launch TOPAS from a directory that already contains it.

## Build

From `microdosimetry_with_LUTs`:

```bash
cmake --build build
```

If `build/` does not exist yet:

```bash
mkdir -p build
cd build
cmake ..
cmake --build .
```

## Phase-Space Replay

`AMF-stage` creates a self-contained TOPAS replay folder. It copies the phase
space pair and `tsed.dat`, writes `replay_amf.txt`, and writes
`amf_run_manifest.txt`.

```bash
./build/microdosimetry_with_LUTs AMF-stage \
  lookup_tables \
  input/PhaseSpace_curved_33mm \
  amf_runtime/staged_runs/PhaseSpace_curved_33mm_AMF_yD \
  AMF_yD \
  --detector water \
  --no-phase-space-precheck
```

`AMF-run` stages the same folder and then launches TOPAS from inside it:

```bash
./build/microdosimetry_with_LUTs AMF-run \
  /path/to/topas \
  lookup_tables \
  input/PhaseSpace_curved_33mm \
  amf_runtime/staged_runs/PhaseSpace_curved_33mm_AMF_yD \
  AMF_yD \
  --detector water \
  --no-phase-space-precheck
```

Return codes:

- `0`: TOPAS exited successfully and expected AMF output files were found.
- `2`: TOPAS returned a nonzero exit code.
- `3`: TOPAS exited successfully, but expected AMF output files were missing.

## Detector Presets

Choose the replay detector with:

```bash
--detector water
--detector silicon
--detector TEgas
```

Presets:

- `water`: `10 cm x 10 cm x 1 mm` `G4_WATER` slab named `AMFScoringVolume`.
- `silicon`: SOI active-layer approximation, `2.93 mm x 3.58 mm x 10 um`
  `G4_Si` slab named `SOISensitiveLayer`.
- `TEgas`: `6.35 mm` radius propane-gas sphere named `TEgasSV`.

The detector should match the surface where the phase space was scored. A phase
space scored on a TEPC gas sphere should be replayed with `--detector TEgas`.
A phase space scored on a slab face should use the corresponding slab detector.

Position the detector with:

```bash
--scoring-x-mm 0 --scoring-y-mm 0 --scoring-z-mm -1387.4
```

Resize slab detectors with:

```bash
--scoring-half-length-x-mm 50 \
--scoring-half-length-y-mm 50 \
--scoring-half-length-z-mm 0.5
```

Resize the TE-gas sphere with:

```bash
--scoring-radius-mm 6.35
```

## Mandatory AMF Parameters

All AMF quantities require a domain radius:

```text
d:Sc/AMF/DomainRadius = 0.28 um
```

In the wrapper, set this with:

```bash
--domain-radius 0.28
```

The allowed AMF development range is:

```text
0.0015 um <= DomainRadius <= 0.5 um
```

The CLI enforces this range. The default is `0.28 um`.

`AMF_yS` also requires:

```text
d:Sc/AMF/NucleusRadius = 3.9 um
d:Sc/AMF/BetaRef = 0.0615 /Gy2
```

In the wrapper:

```bash
--nucleus-radius 3.9 --beta-ref 0.0615
```

## Optional AMF Parameters

Stopping power mode:

```bash
--stopping-power Topas
--stopping-power ExternalTable
```

`Topas` is the recommended default and uses the stopping power calculated by
TOPAS. `ExternalTable` expects a `StoppingPower.txt` file in the TOPAS run
directory. The expected table format is kinetic energy in `MeV/u`, followed by
LET columns in `keV/um` for ions in increasing atomic number.

Step calculator:

```bash
--step-calculator MidStep
--step-calculator PreStep
```

`MidStep` is the default. It evaluates AMF parameters from the average of the
pre-step and post-step kinetic energies. `PreStep` uses the pre-step kinetic
energy.

Phase-space precheck:

```bash
--no-phase-space-precheck
```

Use this when a valid depth-specific `.phsp` file has stale count metadata in
its `.header`. The header is still required by TOPAS, but the precheck is not
allowed to reject the file based on mismatched counts.

## Electron Cut

The generated replay file sets:

```text
d:Ph/Default/CutForElectron = 1000 m
```

AMF includes secondary-electron and delta-ray contributions analytically. The
large electron production range cut suppresses explicit secondary-electron
transport so those contributions stay local to the ion step and are not counted
again as independent transported electrons. If the cut is too low, explicit
electron transport can double count part of the electron contribution and shift
the spectrum toward lower lineal energy.

## Full Simulation Example

For a full simulation, do not include the phase-space source block. Keep your
normal source and beamline, define the detector component, and attach the AMF
scorer to it:

```text
s:Sc/AMF/Quantity = "AMFSpectra"
s:Sc/AMF/Component = "AMFScoringVolume"
s:Sc/AMF/OutputType = "csv"
s:Sc/AMF/OutputFile = "AMF_Spectra"
s:Sc/AMF/IfOutputFileAlreadyExists = "Overwrite"
d:Sc/AMF/DomainRadius = 0.28 um
s:Sc/AMF/StoppingPowerCalculation = "Topas"
s:Sc/AMF/StepCalculator = "MidStep"

d:Ph/Default/CutForElectron = 1000 m
```

For `AMF_yD` or `AMF_yS`, replace the quantity and output file. For `AMF_yS`,
also include `NucleusRadius` and `BetaRef`.

## Outputs

Single-value scorers write one CSV:

```text
<OutputFile>.csv
```

The spectra scorer writes two CSVs:

```text
<OutputFile>.csv
<OutputFile>_MicrodosimetricSpectra.csv
```

The spectra file contains one row per scoring bin. The first three columns are
the scorer bin indices `x,y,z`. The remaining columns are normalized `yd(y)`
values for the lineal-energy bin centers listed in the header row. Plot those
values against the header bin centers to obtain the dose-weighted spectrum.

## Batch Replays

Use `macros/run_amf_depths.py` for multiple depth-specific phase-space pairs.

```bash
python3 macros/run_amf_depths.py \
  --topas /path/to/topas \
  --quantity AMF_yD \
  --input-dir input \
  --lookup-root lookup_tables \
  --staged-root amf_runtime/staged_runs \
  --auto-slab-geometry \
  -- \
  --detector water \
  --no-phase-space-precheck
```

`--auto-slab-geometry` centers a water slab on each phase-space file's z-plane
and defaults to a `10 cm x 10 cm x 1 mm` slab:

```text
X half-length = 50 mm
Y half-length = 50 mm
Z half-length = 0.5 mm
```

The script intentionally skips `.phsp` files that do not have same-base
`.header` files.

## Summarize Runs

After staging or running AMF cases:

```bash
python3 macros/summarize_amf_runs.py \
  --staged-root amf_runtime/staged_runs \
  --output-csv amf_runtime/amf_run_summary.csv \
  --print-missing
```

The summary includes expected output paths and whether each file exists.
