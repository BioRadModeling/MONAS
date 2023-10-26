# GSM2
## Class hierarchy
- **main**: Main function for Survival and RBE
	- **TsGetSurvival**: Class manager for computing MKM and GSM2 survival
		- **TsGSM2**: Class with GSM2 functions
			- **TsSpecificEnergy**: Class for sepcifi energy spectra calculation
	- **TsLinealEnergy**: Class for microdosimetric spectra calculation

## COMPILE
	make all
## RUN
`./monas -parameters...`
## EXAMPLE
./monas -Rd 0.1 -Rc 0.1 -k 0.1 -l 0.1 -a 0.1 -b 0.1 -r 0.1 -topasScorer ./Scorer.phsp
## HELP FOR INPUT PARAMETERS
./monas -help
