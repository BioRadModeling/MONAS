# poly_micro_spectrum

This module builds a **polyenergetic proton microdosimetric spectrum** from a TOPAS phase-space file and a library of precomputed monoenergetic lookup spectra.

It does **not** run the full MONAS survival/RBE workflow. Its job is:

1. read the phase-space file line by line,
2. keep only proton rows,
3. match each proton energy to the **nearest** monoenergetic lookup-table CSV,
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
  csv/
    1mm_linear/
    1mm_logarithmic/
    5um_linear/
    5um_logarithmic/
  les/
  root/
```

The code only uses the `csv/` folders.

It automatically scans the selected folder and reads filenames like:

```text
Proton_71.499794_MeV.csv
```

The proton energy is parsed from the filename.

---

## What gets written to `output/`

After a successful run, you should see these files:

### 1. `proton_matches.csv`
A row-by-row audit file showing how each proton was matched.

Columns:

- `row_index`
- `energy_mev`
- `weight`
- `matched_energy_mev`
- `matched_file`
- `matched_ncpp`

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
mkdir -p poly_micro_spectrum/build
cd poly_micro_spectrum/build
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

From `poly_micro_spectrum/build`:

```bash
./poly_micro_spectrum test-lookup ../lookup_tables 1mm log 72.3
```

What this means:

- `test-lookup` = run only the lookup-table test
- `../lookup_tables` = where the lookup table library lives
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
Matched energy:     71.4998 MeV
Matched CSV:        ../lookup_tables/csv/1mm_logarithmic/Proton_71.499794_MeV.csv
Detected Ncpp:      917.274
```

If this works, your lookup tables are being found and parsed correctly.

---

## Step 3: audit the phase-space file

This reads the phase-space file, filters protons, matches each proton to the nearest lookup spectrum, and writes `proton_matches.csv`.

From `poly_micro_spectrum/build`:

```bash
./poly_micro_spectrum audit-phsp ../lookup_tables 1mm log ../input/PhaseSpace.phsp ../output
```

What this means:

- `audit-phsp` = read the phase-space file and create the proton match audit file
- `../lookup_tables` = lookup-table root folder
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

If this works, the proton filtering and nearest-energy matching are working.

---

## Step 4: build the polyenergetic spectrum

This is the main run.

From `poly_micro_spectrum/build`:

```bash
./poly_micro_spectrum build-spectrum ../lookup_tables 1mm log ../input/PhaseSpace.phsp ../output
```

This command:

- reads the phase-space file,
- keeps only protons,
- matches each proton to the nearest monoenergetic lookup spectrum,
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
poly_micro_spectrum/macros/plot_spectrum.py
```

It reads `poly_spectrum.csv` and saves a JPEG in the output folder.

### Example: plot `yd_y` vs `y_keV_per_um`

From `poly_micro_spectrum/build`:

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
./poly_micro_spectrum build-spectrum ../lookup_tables 1mm log ../input/PhaseSpace.phsp ../output
```

### 1 mm, linear lookup library

```bash
./poly_micro_spectrum build-spectrum ../lookup_tables 1mm linear ../input/PhaseSpace.phsp ../output
```

### 5 um, logarithmic lookup library

```bash
./poly_micro_spectrum build-spectrum ../lookup_tables 5um log ../input/PhaseSpace.phsp ../output
```

### 5 um, linear lookup library

```bash
./poly_micro_spectrum build-spectrum ../lookup_tables 5um linear ../input/PhaseSpace.phsp ../output
```

---

## Sanity checks

After building `poly_spectrum.csv`, it is a good idea to verify normalization.

This small Python script checks whether:

- `sum(f_y) ≈ 1`
- `sum(d_y) ≈ 1`

Example script:

```python
import csv

f_sum = 0.0
d_sum = 0.0

with open("../output/poly_spectrum.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        f_sum += float(row["f_y"])
        d_sum += float(row["d_y"])

print("sum(f_y) =", f_sum)
print("sum(d_y) =", d_sum)
```

If everything is working, both values should be very close to 1.

---

## Summary

If you are new and want the minimum set of commands, use these from `poly_micro_spectrum/build`:

```bash
cmake ..
make -j8
./poly_micro_spectrum test-lookup ../lookup_tables 1mm log 72.3
./poly_micro_spectrum build-spectrum ../lookup_tables 1mm log ../input/PhaseSpace.phsp ../output
python3 ../macros/plot_spectrum.py ../output/poly_spectrum.csv --x y_keV_per_um --y yd_y --output-dir ../output --title "yd(y) vs y" --logx
```

That is enough to:

- load the lookup library,
- process the phase-space file,
- create the polyenergetic spectrum,
- save the JPEG plot.