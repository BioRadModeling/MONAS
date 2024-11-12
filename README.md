# MONAS
## Standalone
C++ toolkit that calculates cell survival and dose dependent RBE from micorodosimetric spectrum from TOPAS MC simulation. 
This is NOT an extension of TOPAS, therefore this should be compiled and run in a separate "build" folder by following the instructions on /Standalone/README.md.
This version was updated in July 2024 and contains the most recent implementation available for the GSM2 model.

## Topas Extension
Topas extension which is called at the end of the TOPAS MC simulation. It takes the TOPAS simulated y-spectrum to calculate cell survival and RBE.
