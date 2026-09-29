#ifndef AMFDose_hh
#define AMFDose_hh

// ***********************************************************************************************
//  This TOPAS extension was developed at Mayo Clinic Florida to score microdosimetric spectra 
//  and deterministic parameters using the analytical microdoimsetric function (AMF) introduced by 
//  Sato et al. in 2023, originally implemented in the Particle and Heavy Ion Transport System. 
//    
//  When using this extension, please cite the following publications:                    
//  [1]	S. Hartzell, A. Parisi, T. Sato, C. J. Beltran, and K. M. Furutani, 
//      "Extending TOPAS with an analytical microdosimetric function: application and 
//	benchmarking with nBio track structure simulations," Physics in Medicine & Biology, 
//	4/23/2025 2025, doi: doi.org/10.1088/1361-6560/adcfec.
//  [2]	T. Sato et al., "Improvement of the hybrid approach between Monte Carlo simulation 
//	and analytical function for calculating microdosimetric probability densities in 
//	macroscopic matter,"  vol. 68, ed. Phys. Med. Biol., 2023.
// 
//  Contacts: Shannon Hartzell, hartzell.shannon@mayo.edu
// ***********************************************************************************************

#include "TsVBinnedScorer.hh"

class AMFDose : public TsVBinnedScorer
{
public:
	AMFDose(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

	virtual ~ AMFDose();
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
};

#endif
