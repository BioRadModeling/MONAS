# microdosimetry_with_LUTs

This module builds a **polyenergetic proton microdosimetric spectrum** from a TOPAS phase-space file and a library of precomputed monoenergetic lookup spectra.

It does **not** run the full MONAS survival/RBE workflow. Its job is:

1. read the phase-space file line by line,
2. keep only proton rows,
3. bracket each proton energy between two monoenergetic lookup spectra and linearly interpolate between them,
4. combine all matched spectra using the proton weights,
5. save the final polyenergetic spectrum as a CSV,
6. optional: plot selected quantities from that CSV as a JPEG.

---

## What this code expects

### Phase-space input

The phase-space file is expected to be a plain-text TOPAS phase-space file.

This code uses only these columns:

- **column 6** = proton kinetic energy in MeV
- **column 7** = particle weight
- **column 8** = particle type / PDG code

Only rows with **PDG = 2212** are kept, which means **protons only**.

### Lookup-table input

The lookup tables should already exist in this folder structure:

```text
lookup_tables/
  DeCunha/
    1mm_linear/
    1mm_logarithmic/
    5um_linear/
    5um_logarithmic/
```

The code only uses the `csv/` folders.

It automatically scans the selected folder and reads filenames like:

```text
Proton_71.499794_MeV.csv
```

The proton energy is parsed from the filename.

---

## How energy matching works

For each proton in the phase-space file, the code reads its kinetic energy
`E_12` from column 6 and searches the selected LUT library for the two
monoenergetic spectra with energies `E_1` and `E_2` such that:

```text
E_1 <= E_12 <= E_2
```

The monoenergetic spectrum used for that proton is then built by linear
interpolation:

```text
f_E12(y) = f_E2(y) + ((E_12 - E_2) / (E_1 - E_2)) * (f_E1(y) - f_E2(y))
```

This is equivalent to:

```text
f_E12(y) = w_1 * f_E1(y) + w_2 * f_E2(y)
```

with:

```text
w_1 = (E_2 - E_12) / (E_2 - E_1)
w_2 = (E_12 - E_1) / (E_2 - E_1)
```

The interpolation is applied after each monoenergetic LUT has already been
converted onto the common target `y` grid used for the final spectrum.

### Out-of-range energies

If a proton energy lies outside the available LUT range, the code keeps the
current endpoint behavior:

- if `E_12` is below the minimum LUT energy, it uses the lowest-energy LUT only
- if `E_12` is above the maximum LUT energy, it uses the highest-energy LUT only

That means the interpolation weights collapse to one endpoint:

- lower endpoint only: `w_1 = 1`, `w_2 = 0`
- upper endpoint only: `w_1 = 1`, `w_2 = 0` when both indices are the same

### Cartechini fallback behavior

When the selected family is `Cartechini`, the code behaves in two regions:

- for energies at or below the maximum Cartechini LUT energy, interpolation is
  done only within the Cartechini subset
- for energies above the maximum Cartechini LUT energy, the code switches
  entirely to the fallback `DeCunha/1mm_logarithmic` subset and performs the
  interpolation there

---

## How rebinning works

Before spectra from different proton energies can be combined, each
monoenergetic LUT is converted onto the same target `y` grid used for the final
polyenergetic spectrum. The code uses different rebinning strategies for the
two LUT families because the source data are stored differently.

### DeCunha LUTs: deterministic log-bin overlap

DeCunha LUTs store raw histogram counts `N(y)` together with `Ncpp`.

For these LUTs, the code:

1. interprets each source bin as carrying a fixed amount of spectral weight
2. compares that source bin to every target bin in `log10(y)` space
3. computes the fraction of overlap between the source bin and each target bin
4. transfers the source-bin weight to the target grid according to those
   overlap fractions

Because the rebinning is based on exact overlap in logarithmic bin space:

- no random sampling is used
- the result is noise-free
- the result is exactly reproducible from run to run

After that, the rebinned counts are converted into the internal monoenergetic
quantity used later in accumulation.

### Cartechini LUTs: stochastic sampling

Cartechini LUTs store tabulated `f(y)` values on a relatively coarse native
grid rather than raw histogram counts.

For these LUTs, the code first converts the source `f(y)` curve into interval
masses by integrating `f(y)` over each source interval. It then rebins those
interval masses onto the target grid by stochastic sampling:

1. build a cumulative distribution from the interval masses
2. repeatedly sample one source interval according to its weight
3. sample one `y` value inside that interval, uniformly in `log10(y)`
4. place that sampled `y` into the corresponding target bin

This Monte Carlo rebinning helps avoid comb-like empty bins that can appear if a
coarse Cartechini source grid is mapped directly onto a finer target grid.

Once the rebinned sample counts are obtained, they are normalized back into the
internal monoenergetic `f(y)`-like quantity used in the polyenergetic
accumulation.

If the stochastic rebinner produces no usable counts, the code falls back to a
deterministic overlap-based rebinning for robustness.

---

## What gets written to `output/`

After a successful run, you should see these files:

### 1. `proton_matches.csv`
A row-by-row audit file showing how each proton was matched.

Columns:

- `row_index`
- `energy_mev`
- `weight`
- `lower_matched_energy_mev`
- `upper_matched_energy_mev`
- `lower_interpolation_weight`
- `upper_interpolation_weight`
- `lower_matched_file`
- `upper_matched_file`
- `lower_matched_ncpp`
- `upper_matched_ncpp`
- `lower_matched_family`
- `upper_matched_family`

### 2. `poly_spectrum.csv`
The final polyenergetic spectrum.

Columns:

- `y_keV_per_um`
- `f_y`
- `yf_y`
- `d_y`
- `yd_y`

### 3. JPEG plot files
Examples:
- `yd_y_vs_y_keV_per_um.jpg`
- `d_y_vs_y_keV_per_um.jpg`
- `f_y_vs_y_keV_per_um.jpg`

---

## LET summary calculations

The `calculate-let` mode computes track-averaged LET and dose-averaged LET
from one phase-space file using the element-specific LET lookup tables in:

```text
lookup_tables/LET/
```

Expected LET files are:

```text
H_water.txt
He_water.txt
Li_water.txt
Be_water.txt
B_water.txt
C_water.txt
N_water.txt
O_water.txt
```

Each file is interpreted as:

```text
E(MeV) LET(keV/um)
```

The command is independent of the DeCunha and Cartechini spectrum modes:

```bash
./microdosimetry_with_LUTs calculate-let ../lookup_tables ../input/PhaseSpace_33mm.phsp ../output
```

For each element table, the same charged phase-space rows are processed against
that element's `LET(E)` curve. Neutral particles are ignored. The three output
groups are:

- `all_charged`: every charged phase-space row accepted by the LET reader
- `protons`: charged rows with PDG code `2212`
- `other_charged`: accepted charged rows with PDG code other than `2212`

The output directory receives one CSV per element:

```text
let_summary_H.csv
let_summary_He.csv
let_summary_Li.csv
let_summary_Be.csv
let_summary_B.csv
let_summary_C.csv
let_summary_N.csv
let_summary_O.csv
```

Each file has this format:

```csv
group,track_averaged_LET_keV_per_um,dose_averaged_LET_keV_per_um
all_charged,...
protons,...
other_charged,...
```

The averages are calculated as:

```text
track_averaged_LET = sum(w * LET) / sum(w)
dose_averaged_LET  = sum(w * LET * LET) / sum(w * LET)
```

where `w` is the phase-space weight from column 7 and `LET` is linearly
interpolated from the selected element table at the particle energy from column
6. Energies outside the LET table range are clamped to the nearest endpoint.
If a group has no contributing particles, its values are written as `nan`.

---

## Before you start

You need:

- a C++17 compiler
- CMake
- Python 3
- `matplotlib` for plotting

If `matplotlib` is missing, install it with:

```bash
python3 -m pip install matplotlib
```

---

## Step 1: build the executable

From the **MONAS repo root**:

```bash
mkdir -p microdosimetry_with_LUTs/build
cd microdosimetry_with_LUTs/build
cmake ..
make -j8
```

If the build succeeds, it creates an executable named:

```text
poly_micro_spectrum
```

---

## Step 2: test that the lookup library loads correctly

This is the safest first test. It does **not** read the phase-space file yet.

From `microdosimetry_with_LUTs/build`:

```bash
./microdosimetry_with_LUTs test-lookup ../lookup_tables DeCunha 1mm log 72.3
```

What this means:

- `test-lookup` = run only the lookup-table test
- `../lookup_tables` = where the lookup table library lives
- `DeCunha` = name of the LUT to be used
- `1mm` = use the 1 mm lookup library
- `log` = use the logarithmic energy library
- `72.3` = test proton energy in MeV

Expected output looks like this:

```text
Loaded library successfully.
Library folder:     "../lookup_tables/csv/1mm_logarithmic"
Number of tables:   300
Y bins:             3000
Requested energy:   72.3 MeV
Lower matched energy:   71.4998 MeV
Upper matched energy:   74.5801 MeV
Lower weight:           0.740218
Upper weight:           0.259782
Lower matched file:     ../lookup_tables/csv/1mm_logarithmic/Proton_71.499794_MeV.csv
Upper matched file:     ../lookup_tables/csv/1mm_logarithmic/Proton_74.580093_MeV.csv
Lower detected Ncpp:    917.274
Upper detected Ncpp:    910.492
```

If this works, your lookup tables are being found and parsed correctly.

---

## Step 3: audit the phase-space file

This reads the phase-space file, filters protons, brackets each proton energy
between two LUT energies, computes the interpolation weights, and writes
`proton_matches.csv`.

From `microdosimetry_with_LUTs/build`:

```bash
./microdosimetry_with_LUTs audit-phsp ../lookup_tables Decunha 1mm log ../input/PhaseSpace.phsp ../output
```

What this means:

- `audit-phsp` = read the phase-space file and create the proton match audit file
- `../lookup_tables` = lookup-table root folder
- `DeCunha` = name of the LUT to be used
- `1mm` = use the 1 mm library
- `log` = use the logarithmic library
- `../input/PhaseSpace.phsp` = input phase-space file
- `../output` = folder where output CSV files should be saved

Expected output looks like:

```text
Phase-space audit completed.
Library folder:     "../lookup_tables/csv/1mm_logarithmic"
Phase-space file:   "../input/PhaseSpace.phsp"
Protons found:      552711
Output CSV:         "../output/proton_matches.csv"
```

If this works, the proton filtering and interpolation matching are working.

---

## Step 4: build the polyenergetic spectrum

This is the main run.

From `microdosimetry_with_LUTs/build`:

```bash
./microdosimetry_with_LUTs build-spectrum ../lookup_tables Decunha 1mm log ../input/PhaseSpace.phsp ../output
```

This command:

- reads the phase-space file,
- keeps only protons,
- interpolates each proton between two monoenergetic lookup spectra,
- accumulates the weighted contributions,
- writes both:
  - `../output/proton_matches.csv`
  - `../output/poly_spectrum.csv`

Expected output:

```text
Phase-space processing completed.
Library folder:     "../lookup_tables/csv/1mm_logarithmic"
Phase-space file:   "../input/PhaseSpace.phsp"
Protons found:      552711
Match CSV:          "../output/proton_matches.csv"
Poly spectrum CSV:  "../output/poly_spectrum.csv"
```

---

## Step 5: plot the spectrum as a JPEG

The plotting script lives in:

```text
microdosimetry_with_LUTs/macros/plot_spectrum.py
```

It reads `poly_spectrum.csv` and saves a JPEG in the output folder.

### Example: plot `yd_y` vs `y_keV_per_um`

From `microdosimetry_with_LUTs/build`:

```bash
python3 ../macros/plot_spectrum.py \
  ../output/poly_spectrum.csv \
  --x y_keV_per_um \
  --y yd_y \
  --output-dir ../output \
  --title "yd(y) vs y" \
  --logx
```

This should create:

```text
../output/yd_y_vs_y_keV_per_um.jpg
```

### Example: Plot `d_y` vs `y_keV_per_um`

```bash
python3 ../macros/plot_spectrum.py \
  ../output/poly_spectrum.csv \
  --x y_keV_per_um \
  --y d_y \
  --output-dir ../output \
  --title "d(y) vs y" \
  --logx
```

---

## How to choose the lookup library

The command line lets you choose two things:

### 1. Voxel size
Use either:

- `1mm`
- `5um`

### 2. Energy spacing
Use either:

- `linear`
- `log`

Examples:

### 1 mm, logarithmic lookup library

```bash
./microdosimetry_with_LUTs build-spectrum ../lookup_tables Decunha 1mm log ../input/PhaseSpace.phsp ../output
```

### 1 mm, linear lookup library

```bash
./microdosimetry_with_LUTs build-spectrum ../lookup_tables 1mm Decunha linear ../input/PhaseSpace.phsp ../output
```

### 5 um, logarithmic lookup library

```bash
./microdosimetry_with_LUTs build-spectrum ../lookup_tables 5um Decunha log ../input/PhaseSpace.phsp ../output
```

### 5 um, linear lookup library

```bash
./microdosimetry_with_LUTs build-spectrum ../lookup_tables 5um Decunha linear ../input/PhaseSpace.phsp ../output
```

---

## Summary

If you are new and want the minimum set of commands, use these from `poly_micro_spectrum/build`:

```bash
cmake ..
make -j8
./microdosimetry_with_LUTs test-lookup ../lookup_tables Decunha 1mm log 72.3
./microdosimetry_with_LUTs build-spectrum ../lookup_tables Decunha 1mm log ../input/PhaseSpace.phsp ../output
./microdosimetry_with_LUTs calculate-let ../lookup_tables ../input/PhaseSpace.phsp ../output
python3 ../macros/plot_spectrum.py ../output/poly_spectrum.csv --x y_keV_per_um --y yd_y --output-dir ../output --title "yd(y) vs y" --logx
```

That is enough to:

- load the lookup library,
- process the phase-space file,
- create the polyenergetic spectrum,
- calculate LET summaries,
- save the JPEG plot.
