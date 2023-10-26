# MONAS Standalone
## Introduction
This is the standalone versione of the MONAS extension for TOPAS.\
It is a C++-based code which can be compiled and executed independentely from TOPAS.\
It takes as input the PhaseSpace produced by TOPAS microdosimetric extension which corresponds to a file where each raw is the event-by-event mircorodosimetric lineal energy **y**.

## Code structure
- **main**: Main function for Survival and RBE
	- **TsGetSurvival**: Class manager for computing MKM and GSM2 survival
		- **TsGSM2**: Class with GSM2 functions
			- **TsSpecificEnergy**: Class for sepcifi energy spectra calculation
	- **TsLinealEnergy**: Class for microdosimetric spectra calculation

## COMPILE
```
make all
```
## RUN
```
./monas -parameters...
```
Example MKM:
```
./monas -SMKM -DSMKM -fGetParticleContribution -fSetMultieventStatistic 10000 -topasScorer $topas -MKM_rDomain 0.4571564 -MKM_rNucleus 8.00 -MKM_alpha 0.1626047 -MKM_beta 0.0789754 -MKM_y0 150 -MKM_alphaX 0.313 -MKM_betaX 0.030.6 -Doses 0 10 1 #HSG
```
Example GSM2:
```
./monas  -GSM2 -fSetMultiEventStatistic 100000 -fGetParticleContribution -topasScorer $topas -GSM2_rDomain 0.4 -GSM2_rNucleus 8 -GSM2_alphaX 0.19 -GSM2_betaX 0.05 -GSM2_a 0.03691869  -GSM2_b 0.1822144  -GSM2_r 3.640363  -GSM2_kappa 0.0750 -GSM2_lambda 0.0000750 -Doses 0 10 1 #HSG
```

## HELP FOR INPUT PARAMETERS
```
./monas -help
```
