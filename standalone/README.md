# MONAS Standalone
## Introduction
This is the standalone versione of the MONAS extension for TOPAS.\
It is a C++-based code which can be compiled and executed independentely from TOPAS.\
It takes as input the PhaseSpace produced by TOPAS microdosimetric extension which corresponds to a file where each raw is the event-by-event mircorodosimetric lineal energy **y**.

## Code structure
- **main**: Main function for Survival and RBE
	- **TsGetSurvival**: Class manager for computing MKM and GSM2 survival
		- **TsGSM2**: Class with GSM2 functions
			- **TsSpecificEnergy**: Class for specific energy spectra calculation
	- **TsLinealEnergy**: Class for microdosimetric spectra calculation

## LIBRARIES NEEDED
- GSL (Ubuntu): ```sudo apt-get install libgsl-dev```

## COMPILE
```
cmake .
make -j
```
## RUN
```
./bin/monas -parameters...
```
Example MKM (parameters for HSG cells):
```
./bin/monas -SMKM -DSMKM -fGetParticleContribution -fSetMultieventStatistic 10000 -topasScorer_08 $topas_08.phsp -topasScorer_8 $topas_8.phsp -MKM_rDomain 0.4571564 -MKM_rNucleus 8.00 -MKM_alpha 0.1626047 -MKM_beta 0.0789754 -MKM_y0 150 -MKM_alphaX 0.313 -MKM_betaX 0.030.6 -Doses 0 10 1
```
Example GSM2 (-CellLine requires "-H460" or "-H1437", otherwise no specific cell line is selected and default parameters are used, Nucleus Scorer is NOT mandatory):
```
./bin/monas -GSM2 -CellLine -fGetParticleContribution -topasScorerDomain $topas_domain_size.phsp -topasScorerNucleus $topas_nucleus_size.phsp -Doses 0 10 1
```

### Deterministic binned-spectrum MKM mode

Use `-MKMFromSpectra` when the input is already a binned lineal-energy spectrum instead of event-by-event phase-space samples. This mode reads each bin directly and can run `-MKMSatCorr`, `-MKMnonPoiss`, `-SMKM`, and/or `-GSM2` for one or more spectra in the same command.

Input spectrum formats:

- `name:poly:file.csv`
  - Generic CSV with columns `y_keV_per_um,f_y,yf_y,d_y,yd_y`.
- `name:topas:ySpecfile.txt`
  - TOPAS `ySpecfile.txt` style output.
- `name:amf:spectra.csv:moments.csv:x,y,z`
  - AMF spectrum table plus moments table for the voxel at `x,y,z`.

Outputs are written to `-outputDir` and include one file per source/model plus:

- `MKM_from_spectra_summary.csv`
- `MKM_from_spectra_manifest.txt`

Example TOPAS spectrum run:

```
./bin/monas -MKMFromSpectra -outputDir results/mkm_spectra \
  -spectrum TOPAS:topas:examples/ySpecfile.txt \
  -MKMSatCorr -MKMnonPoiss -SMKM \
  -MKM_rDomain 0.44 -MKM_rNucleus 8.0 \
  -MKM_alpha 0.19 -MKM_beta 0.07 -MKM_y0 150 \
  -MKM_alphaX 0.19 -MKM_betaX 0.05 \
  -Doses 0 10 1
```

Example multi-source comparison:

```
./bin/monas -MKMFromSpectra -outputDir results/mkm_compare \
  -spectrum DeCunha:poly:decunha_spectrum.csv \
  -spectrum AMF:amf:amf_spectra.csv:amf_moments.csv:0,0,0 \
  -spectrum TOPAS:topas:ySpecfile.txt \
  -MKMSatCorr -SMKM \
  -Doses 0 10 1
```

Example GSM2 binned-spectrum run:

```
./bin/monas -MKMFromSpectra -outputDir results/gsm2_spectra \
  -spectrum DeCunha:poly:decunha_spectrum.csv \
  -GSM2 \
  -fSetMultieventStatistic 10000 \
  -GSM2_rDomain 0.42 -GSM2_rNucleus 6.0 \
  -GSM2_alphaX 0.19 -GSM2_betaX 0.05 \
  -GSM2_a 0.1 -GSM2_b 0.1 -GSM2_r 0.1 \
  -Doses 0 10 1
```

In binned-spectrum GSM2 mode, MONAS converts the input frequency spectrum `f(y)` to a specific-energy spectrum using:

```
z [Gy] = 0.204 * y [keV/um] / (2*r_det [um])^2
```

This conversion is valid for spherical targets only. Different detector or target geometries require the corresponding geometry-specific conversion from lineal energy to specific energy. After changing variables from `y` to `z`, MONAS normalizes the converted `f(z)` with bin widths and samples the cumulative dose-weighted single-event distribution `z*f(z)` through the same GSM2 damage workflow used by event-derived spectra. Particle-contribution outputs are not available from binned spectra because the input histogram does not retain particle identity.

A minimal smoke example is available in `examples/run_binned_spectrum_smoke.sh`. It runs the DeCunha-style generic CSV, AMF, and TOPAS fixture spectra together and checks that the expected output files are created. A GSM2-specific binned-spectrum smoke example is available in `examples/run_binned_gsm2_smoke.sh`.

## HELP FOR INPUT PARAMETERS
```
./bin/monas -help
```
