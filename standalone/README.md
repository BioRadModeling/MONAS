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
Example GSM2 (-CellLine requires "-H460" or "-H1437", otherwise no specific cell line is selected and default parameters are used):
```
./bin/monas -GSM2 -CellLine -fGetParticleContribution -topasScorer_08 $topas_08.phsp -topasScorer_8 $topas_8.phsp -Doses 0 10 1
```

## HELP FOR INPUT PARAMETERS
```
./bin/monas -help
```
