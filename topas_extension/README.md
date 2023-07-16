TOPAS MC Microdosimetric extension
========

### Reference paper:
G.Cartechini et al. Under Submission.

#### Extension Parameters

 - `b:Sc/.../GetRBEWithMKModel = "True"` Get RBE and Survival with one or more MKM formulations 
 - `sv:Sc/.../GetRBEWithMKModel = 4 "SaturationCorrection" "nonPoisson" ...`
      1. **SaturationCorrection**,  
      2. "meanValues" Then a file will be created, named "PROJECTNAME_survival_MKM.csv", containing the information about the parameters chosen for the simulation and the values of doses delivered and survival observed (a new line for each energy or dose evaluated).
   
      3. "cellValues" This kind of output is supported only by the MonteCarlo calculusType. It is a way to store the values of dose and survival obtained for each single cell irradiated during the monte carlo simulation. Then a directory will be created, named "PROJECTNAME_survival_data". In the directory the user will find a description file named "000_MonteCarlo_parameters.csv", listing the parameters used in the simulation, and a directory with the same name containing the corresponding data. In particular in the subdirectory some file will be created (a file for each level of dose imposed), each one containing two column with the dose delivered and the survival observed for each cell irradiated. When a new simulation is lauched with the same project name, the program will do a check over all the description files in the directory, if the parameters of the simulation are the same of another one already done, then it will enter the related subdirectory and append data there, if not a new description file (with progressive number) and corresponding subdirectory will be created.
      
   **Note**: the cell line is only determined by the chosen model parameters. `-cellType` is only a tag to indicate the cell, useful for bookkeeping.

### Bibliography


