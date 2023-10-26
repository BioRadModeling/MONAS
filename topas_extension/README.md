# MONAS: MicrOdosimetry-based modeling for RBE Assessment

 ## Introduction
 MONAS expands TOPAS microdosimetric extension, by including novel specific energy scorers to calculate the single- and multi-event specific energy microdosimetric distributions at different micrometer scales. These spectra are used as physical input to three different formulations of the **Microdosimetric Kinetic Model (MKM)**, and to the **Generalized Stochastic Microdosimetric Model (GSM2)**, to predict dose-dependent cell survival fraction and RBE.

## Files

 - **TsYScorer**: It is the main file of the extension. This file mages the input parameters, the energy spectra, Survival and RBE calculation, and the output files. This file was already present in the original microdosimetric extension and it improved to include the MONAS calculations.
 - **TsGetSurvivalRBEQualityFactor**: Class that manages the cell survival and RBE calculation from 3 MKM formulations (MKM-z*, mSMKM, DSMKM) and GSM2.
 - **TsGSM2parallel**: Class used to manage the GSM2 cell survival fraction calculation.
 - **TsSpecificEnergy**: In this file it is defined the class for the storage of specific energy distributions: single- and multi-event.
 - **TsLinealEnergy**: Similar to TsSpecificEnergy, but for y-distributions.
 - **TsSOIMicrodosimeter**: In this file is defined the SOI detector geometry. 
 - **TsTrackerHit**:

## Compile
Just follow the topas README.txt file.
```
cd <path_to_topas>/
cmake -DTOPAS_EXTENSIONS_DIR=<path_to_topas_extensions>/MONAS/topas_extension
make -j
```
==NOTE:== You cannot compile the standard microdosimetric extension and monas extension together. They share same classes and the compile will give you an error.
==MONAS has the same features as the standare microdosimetric extensions. You can directly compile MONAS to have the same output==

## Input parameters
### Model independent parameters

### MKM models

 - `GetRBEWithMKModel`: boolean flag to activate MKM models
 `b:Sc/Scorer/GetRBEWithMKModel = "True" #True or False`

 - `MKMCalculation`: string vector with MKM formulations. Options: "SaturationCorrection", "nonPoisson", "DSMKM". **NOTE: when more formulations are selected, they share the same model parameters.**
 `sv:Sc/Scorer/MKMCalculation = 1 "SaturationCorrection"`
 
 - `MKModel_DomainRadius`: Cell domain radius [**um**]
 `u:Sc/Scorer/MKModel_DomainRadius = 0.44`
 
 - `MKModel_NucleusRadius`: Cell nucleus radius [**um**]
 `u:Sc/Scorer/MKModel_NucleusRadius = 3.9`
 
 - `u:Sc/Scorer/MKModel_alpha0 = 0.19`
 - `u:Sc/Scorer/MKModel_beta0 = 0.07`
 - `u:Sc/Scorer/MKModel_y0 = 150`
 
### GSM2 model parameters
 - `b:Sc/Scorer/GetRBEWithGSM2 = "True"`
 - `u:Sc/Scorer/GSM2_DomainRadius = 0.44`
 - `u:Sc/Scorer/GSM2_NucleusRadius = 3.9`
 - `u:Sc/Scorer/GSM2_a = 0.19`
 - `u:Sc/Scorer/GSM2_b = 0.07`
 - `u:Sc/Scorer/GSM2_r = 150`

## Output files

