# MONAS: MicrOdosimetry-based modeling for RBE Assessment

 ## Introduction
 MONAS expands TOPAS microdosimetric extension, by including novel specific energy scorers to calculate the single- and multi-event specific energy microdosimetric distributions at different micrometer scales. These spectra are used as physical input to three different formulations of the **Microdosimetric Kinetic Model (MKM)**, and to the **Generalized Stochastic Microdosimetric Model (GSM2)**, to predict dose-dependent cell survival fraction and RBE.

## Reference
- ```Cartechini, Giorgio, et al. "Integrating microdosimetric in vitro RBE models for particle therapy into TOPAS MC using the MicrOdosimetry-based modeliNg for RBE ASsessment (MONAS) tool." Physics in Medicine & Biology 69.4 (2024): 045005.```
- ```@article{cartechini2024integrating,
  title={Integrating microdosimetric in vitro RBE models for particle therapy into TOPAS MC using the MicrOdosimetry-based modeliNg for RBE ASsessment (MONAS) tool},
  author={Cartechini, Giorgio and Missiaggia, Marta and Scifoni, Emanuele and La Tessa, Chiara and Cordoni, Francesco G},
  journal={Physics in Medicine \& Biology},
  volume={69},
  number={4},
  pages={045005},
  year={2024},
  publisher={IOP Publishing}
}```

# Files

 - **TsYScorer**: It is the main file of the extension. This file mages the input parameters, the energy spectra, Survival and RBE calculation, and the output files. This file was already present in the original microdosimetric extension and it has been extended to include the MONAS calculations.
 - **TsGetSurvivalRBEQualityFactor**: Class that manages the cell survival and RBE calculation from 3 MKM formulations (MKM-z*, mSMKM, DSMKM) and GSM2.
 - **TsGSM2parallel**: Class used to manage the GSM2 cell survival fraction calculation. It uses multithread calculation.
 - **TsSpecificEnergy**: In this file it is defined the class for the storage of specific energy distributions: single- and multi-event.
 - **TsLinealEnergy**: Similar to TsSpecificEnergy, but for y-distributions.
 - **TsSOIMicrodosimeter**: In this file is defined the SOI detector geometry. 
 - **TsTrackerHit**: Class used to collect the particle information in the sensitive volume at each step 

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
 `u:Sc/Scorer/MKModel_NucleusRadius = 8.0`
 
 - `u:Sc/Scorer/MKModel_alpha0 = 0.19`
 - `u:Sc/Scorer/MKModel_beta0 = 0.07`
 - `u:Sc/Scorer/MKModel_y0 = 150`
 
### GSM2 model parameters
 - `b:Sc/Scorer/GetRBEWithGSM2 = "True"`
 - `u:Sc/Scorer/GSM2_DomainRadius = 0.8`
 - `u:Sc/Scorer/GSM2_NucleusRadius = 5.0`
 - `u:Sc/Scorer/GSM2_a = 0.037`
 - `u:Sc/Scorer/GSM2_b = 0.182`
 - `u:Sc/Scorer/GSM2_r = 3.641`

## Scorer Example Parameters
Example of MONAS Scorer using a spherical TEPC detector.

`s:Sc/Scorer/Quantity    = "TsYScorer"
s:Sc/Scorer/Component   = "TEgasSV"
#
# Mandatory parameters
i:Sc/Scorer/NumberOfHistoriesInRun    = So/Demo/NumberOfHistoriesInRun
i:Sc/Scorer/GeometryNumber            = 0
d:Sc/Scorer/SensitiveVolumeRadius     = 6.35 mm      # radius of sensitive volume (except for silicon microdosimeter)
d:Sc/Scorer/TissueEquivalentRadius    = 1.0 um      # radius of equivalent size of tissue equivalent volume (to calculate mean chord length)
d:Sc/Scorer/TransX                    = Ge/AlShellOut/TransX cm          # x position of sensitive volume from world center
d:Sc/Scorer/TransY                    = Ge/AlShellOut/TransY cm          # y position of sensitive volume from world center
d:Sc/Scorer/TransZ                    = Ge/AlShellOut/TransZ cm     # y position of sensitive volume from world center
s:Sc/Scorer/OutputType                = "ASCII"       # OutputType must be ASCII, Binary or ROOT (unit: keV/um)
s:Sc/Scorer/IfOutputFileAlreadyExists = "Overwrite"
#s:Sc/Scorer/IfOutputFileAlreadyExists = "Increment"
#
# Optional parameters
u:Sc/Scorer/LinealEnergyLowerlimit      = 0.1         # in unit of keV/um
u:Sc/Scorer/LinealEnergyUpperlimit     = 1000        # in unit of keV/um
b:Sc/Scorer/IncludeFrequencyMeanLinealEnergy   = "True"
b:Sc/Scorer/IncludeDoseMeanLinealEnergy   = "True"
u:Sc/Scorer/PileupProbability = 0.0
b:Sc/Scorer/GetStatisticInfo =  "False"
i:Sc/Scorer/SpectrumUpdateTimes = 1000
b:Sc/Scorer/SaveSpecificEnergySpecra = "True"

## MONAS PARAMETERS
#MKM
b:Sc/Scorer/GetRBEWithMKModel = "True"
sv:Sc/Scorer/MKMCalculation = 2 "SaturationCorrection" "SMKM"
u:Sc/Scorer/MKModel_alpha0 = 0.19
u:Sc/Scorer/MKModel_beta = 0.07
u:Sc/Scorer/MKModel_betaX = 0.05 #HSG 200 kV Kase et al. 2006
u:Sc/Scorer/MKModel_alphaX = 0.19 #HSG 200 kV Kase et al. 2006
u:Sc/Scorer/MKModel_DomainRadius = 0.44
u:Sc/Scorer/MKModel_NucleusRadius = 8.0
u:Sc/Scorer/MKModel_y0 = 150
u:Sc/Scorer/SetMultieventStatistic = 10000

#GSM2
b:Sc/Scorer/GetRBEWithGSM2 = "True"
u:Sc/Scorer/GSM2_DomainRadius = 0.8
u:Sc/Scorer/GSM2_NucleusRadius = 5.0
u:Sc/Scorer/GSM2_a = 0.037
u:Sc/Scorer/GSM2_b = 0.182
u:Sc/Scorer/GSM2_r = 3.641
u:Sc/Scorer/GSM2_alphaX = 0.19 #HSG 200 kV
u:Sc/Scorer/GSM2_betaX = 0.05 #HSG 200 kV

uv:Sc/Scorer/SurvivalDoses = 3 1 10 2`


