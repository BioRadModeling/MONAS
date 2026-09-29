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

Use `-MKMFromSpectra` when the input is already a binned lineal-energy spectrum instead of event-by-event phase-space samples. This mode does not sample synthetic events. It reads each bin directly and can run `-MKMSatCorr`, `-MKMnonPoiss`, and/or `-SMKM` for one or more spectra in the same command.

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

A minimal smoke example is available in `examples/run_binned_spectrum_smoke.sh`. It runs the DeCunha-style generic CSV, AMF, and TOPAS fixture spectra together and checks that the expected output files are created.

## HELP FOR INPUT PARAMETERS
```
./bin/monas -help
```

## GSM2 WITH DETERMINISTIC BINNED SPECTRA

`-GSM2` can now be combined with `-MKMFromSpectra`, so GSM2 survival/RBE can be
computed directly from a binned microdosimetric spectrum (poly CSV, AMF or TOPAS
`ySpecfile.txt`) instead of from an event-by-event lineal-energy list.

```
./bin/monas -MKMFromSpectra -outputDir results/gsm2 \
  -spectrum TOPAS:topas:examples/ySpecfile.txt \
  -GSM2 -H460 \
  -fSetMultiEventStatistic 100000 \
  -Doses 0 10 1
```

Each source writes `<source>_GSM2.csv` and `<source>_yDvsRBE10.csv`, and GSM2 rows are
appended to `MKM_from_spectra_summary.csv` with the model name `GSM2`.

How it works: the binned f(y) is converted into the single-event specific-energy
spectrum f(z) on the same logarithmic z grid used by the event-by-event path, spreading
each lineal-energy bin over every z bin it overlaps. The domain and nucleus spectra
differ only through the y->z conversion radius (`-GSM2_rDomain` / `-GSM2_rNucleus`),
which mirrors the behaviour of the event-by-event path when no separate nucleus scorer
is supplied. Everything downstream (kappa/lambda, initial damage sampling, multi-event
nucleus spectrum, LQ fit, RBE10) is the unchanged GSM2 implementation.

Limitations:

* A binned spectrum carries no per-particle decomposition, so `-fGetParticleContribution`
  reports the total component only.
* `-DSMKM` is still not available with `-MKMFromSpectra`.

`-GSM2_nDomains N` overrides the number of domains per nucleus, which otherwise defaults
to `floor((rNucleus/rDomain)^2)` in both GSM2 modes.

A smoke example is available in `examples/run_binned_spectrum_gsm2_smoke.sh`.
