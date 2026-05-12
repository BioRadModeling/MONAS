# AMF Runtime

This folder holds TOPAS-backed AMF replay assets for the `microdosimetry_with_LUTs`
workflow.

- `templates/`: reusable TOPAS parameter-file templates.
- `staged_runs/`: generated per-case TOPAS run directories.
- `reference_outputs/`: manually validated AMF outputs used for regression checks.

Generated staged runs should be treated as disposable. Reference outputs should only
be added after a manual TOPAS AMF run has been reviewed.

## What AMF needs

For one AMF replay, provide:

- a TOPAS executable with the AMF extension compiled in,
- a phase-space pair with the same base name, for example
  `input/PhaseSpace_curved_33mm.phsp` and
  `input/PhaseSpace_curved_33mm.header`,
- `lookup_tables/AMF/tsed.dat`,
- the AMF quantity: `AMFSpectra`, `AMF_yD`, or `AMF_yS`.

The `.header` file is required by TOPAS phase-space replay. This workflow only
uses it as TOPAS input metadata; it does not trust particle counts in the header
for dose or spectra calculations.

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

## Stage only

Staging creates a self-contained TOPAS run directory. It copies the phase-space
pair and `tsed.dat`, writes `replay_amf.txt`, and writes an
`amf_run_manifest.txt` audit file.

From `microdosimetry_with_LUTs`:

```bash
./build/microdosimetry_with_LUTs AMF-stage \
  lookup_tables \
  input/PhaseSpace_curved_33mm \
  amf_runtime/staged_runs/PhaseSpace_curved_33mm_yD \
  AMF_yD
```

Use `AMFSpectra` for the full spectrum, `AMF_yD` for dose-weighted mean lineal
energy, or `AMF_yS` for saturation-corrected dose mean lineal energy.

## Run TOPAS

`AMF-run` stages the case and then runs TOPAS from inside the staged directory.
Pass the TOPAS executable as the first argument.

```bash
./build/microdosimetry_with_LUTs AMF-run \
  /path/to/topas \
  lookup_tables \
  input/PhaseSpace_curved_33mm \
  amf_runtime/staged_runs/PhaseSpace_curved_33mm_yD \
  AMF_yD
```

After TOPAS exits, the CLI reports the expected AMF output files as `FOUND` or
`MISSING`.

Return codes:

- `0`: TOPAS exited successfully and expected AMF output files were found.
- `2`: TOPAS returned a nonzero exit code.
- `3`: TOPAS exited successfully, but expected AMF output files were missing.

## Batch repeated depth runs

Use `macros/run_amf_depths.py` when you have multiple depth-specific
phase-space pairs in one input folder. The script discovers matching
`.phsp/.header` pairs and calls `AMF-stage` or `AMF-run` once per pair.

Stage every usable pair:

```bash
python3 macros/run_amf_depths.py \
  --stage-only \
  --quantity AMF_yD \
  --input-dir input \
  --lookup-root lookup_tables \
  --staged-root amf_runtime/staged_runs
```

Run every usable pair with TOPAS:

```bash
python3 macros/run_amf_depths.py \
  --topas /path/to/topas \
  --quantity AMF_yD \
  --input-dir input \
  --lookup-root lookup_tables \
  --staged-root amf_runtime/staged_runs
```

Limit the batch to a subset:

```bash
python3 macros/run_amf_depths.py \
  --stage-only \
  --quantity AMF_yD \
  --pattern "PhaseSpace_*mm.phsp"
```

Pass AMF scorer options after `--`:

```bash
python3 macros/run_amf_depths.py \
  --stage-only \
  --quantity AMF_yS \
  --pattern "PhaseSpace_curved_33mm.phsp" \
  -- \
  --domain-radius 0.28 \
  --nucleus-radius 3.9 \
  --beta-ref 0.0615
```

The script intentionally skips `.phsp` files that do not have a same-base
`.header` file, because TOPAS phase-space replay requires the pair.

### Recommended Phase-Space Slab Geometry

For depth-specific phase-space replay, use `--auto-slab-geometry`. This makes
each phase-space pair score in a macroscopic water slab centered on that file's
own phase-space z-plane:

```bash
python3 macros/run_amf_depths.py \
  --topas /Applications/TOPAS/OpenTOPAS-build/topas \
  --quantity AMF_yD \
  --input-dir input \
  --lookup-root lookup_tables \
  --staged-root amf_runtime/staged_runs \
  --pattern "PhaseSpace_curved_33mm.phsp" \
  --auto-slab-geometry \
  -- \
  --no-phase-space-precheck
```

Defaults for `--auto-slab-geometry` are:

```text
X half-length = 50 mm
Y half-length = 50 mm
Z half-length = 0.5 mm
```

That corresponds to a `10 cm x 10 cm x 1 mm` water slab, matching the slab
scale used in the AMF TOPAS paper for 1 mm depth resolution. The script derives
the slab z-position from the midpoint of column 3 in each `.phsp` file, whose
header defines the units as cm, and converts it to mm for TOPAS geometry.

For `PhaseSpace_curved_33mm.phsp`, this currently gives:

```text
Z center = -1387.4 mm
World half-length = 148.79 cm
```

The x/y slab size can be adjusted if the field of interest is smaller or larger:

```bash
--slab-half-x-mm 50 --slab-half-y-mm 50 --slab-half-z-mm 0.5
```

## Summarize completed runs

Use `macros/summarize_amf_runs.py` after staging or running AMF cases. It scans
`amf_run_manifest.txt` files and writes one CSV row per staged run, including
the expected output paths and whether each file exists.

```bash
python3 macros/summarize_amf_runs.py \
  --staged-root amf_runtime/staged_runs \
  --output-csv amf_runtime/amf_run_summary.csv \
  --print-missing
```

Summarize only one quantity:

```bash
python3 macros/summarize_amf_runs.py \
  --staged-root amf_runtime/staged_runs \
  --output-csv amf_runtime/amf_yD_summary.csv \
  --quantity AMF_yD
```

The script returns:

- `0` if every summarized run has its expected AMF output files.
- `1` if one or more summarized runs are missing expected AMF output files.
- `2` if the summary command itself fails.

## Accuracy defaults

The generated TOPAS file sets:

```text
d:Ph/Default/CutForElectron = 1000 m
d:Sc/AMF/DomainRadius = 0.28 um
s:Sc/AMF/StoppingPowerCalculation = "Topas"
s:Sc/AMF/StepCalculator = "MidStep"
```

The high electron range cut suppresses explicit secondary-electron transport so
the electron contribution already built into the AMF spectra is not counted a
second time.

For `AMF_yS`, the generated file also sets:

```text
d:Sc/AMF/NucleusRadius = 3.9 um
d:Sc/AMF/BetaRef = 0.0615 /Gy2
```

## Useful options

Override the domain radius:

```bash
--domain-radius 0.28
```

Override saturation parameters for `AMF_yS`:

```bash
--nucleus-radius 3.9 --beta-ref 0.0615
```

Move or resize the water scoring volume:

```bash
--scoring-half-length-mm 1.25 \
--scoring-x-mm 0.0 \
--scoring-y-mm 0.0 \
--scoring-z-mm 33.0
```

Use AMF's external stopping-power table path instead of TOPAS stopping power:

```bash
--stopping-power ExternalTable
```

Use pre-step parameters rather than mid-step parameters:

```bash
--step-calculator PreStep
```

Disable TOPAS phase-space precheck when the `.header` count metadata does not
match the actual depth-specific `.phsp` rows:

```bash
--no-phase-space-precheck
```

This still requires the `.header` file for TOPAS phase-space replay metadata,
but it prevents TOPAS from rejecting a valid phase-space file because stale
header particle counts disagree with the actual file.

## Expected outputs

For `AMF_yD` and `AMF_yS`, the expected scorer output is:

```text
<stagedRunDir>/<phaseSpaceBase>_AMF_yD.csv
<stagedRunDir>/<phaseSpaceBase>_AMF_yS.csv
```

For `AMFSpectra`, the extension is expected to produce:

```text
<stagedRunDir>/<phaseSpaceBase>_AMFSpectra.csv
<stagedRunDir>/<phaseSpaceBase>_AMFSpectra_MicrodosimetricSpectra.csv
```

The spectra CSV contains the normalized `yd(y)` bins per scoring voxel.
