// Scorer for AMFDose

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

#include "AMFDose.hh"

AMFDose::AMFDose(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
								 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	SetUnit("Gy");
}

AMFDose::~AMFDose() {;}

G4bool AMFDose::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}
	G4double KE_prestep = aStep->GetPreStepPoint()->GetKineticEnergy();
	G4double KE_poststep = aStep->GetPostStepPoint()->GetKineticEnergy();
	G4double KE_mean = (KE_prestep + KE_poststep) / 2;
	G4double AtomicMass = aStep->GetTrack()->GetDefinition()->GetAtomicMass();
	G4double ene = KE_mean / AtomicMass;
	G4double edep = aStep->GetTotalEnergyDeposit();
	G4int AtomicNumber = aStep->GetTrack()->GetDefinition()->GetAtomicNumber();
	if ( edep > 0. ) {
		G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();
		ResolveSolid(aStep);
		if (AtomicNumber >= 1 && AtomicNumber <= 18 && ene >= 0.025) {
		G4double dose = edep / ( density * fSolid->GetCubicVolume() );
		dose *= aStep->GetPreStepPoint()->GetWeight();
		AccumulateHit(aStep, dose);

		return true;
		}
	}
	return false;
}
