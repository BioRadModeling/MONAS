#ifndef SCOREAMF_yS_HH
#define SCOREAMF_yS_HH

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
#include "G4ParticleDefinition.hh"

#include <vector>
#include <map>
#include <utility>
#include <string>

class LookupEntry_yS {
public:
    double energy;
    std::vector<double> stoppingPowers;

    LookupEntry_yS(double e, const std::vector<double>& sp)
        : energy(e), stoppingPowers(sp) {}
};

class ScoreAMF_yS : public TsVBinnedScorer {
public:
    ScoreAMF_yS(TsParameterManager* pM,
                TsMaterialManager* mM,
                TsGeometryManager* gM,
                TsScoringManager* scM,
                TsExtensionManager* eM,
                G4String scorerName,
                G4String quantity,
                G4String outFileName,
                G4bool isSubScorer = false);

    virtual ~ScoreAMF_yS();

    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* history) override;
    virtual G4int CombineSubScorers() override;

private:
    enum class StoppingPowerMode {
        Topas,
        ExternalTable
    };

    enum class StepCalculationMode {
        PreStep,
        MidStep
    };

    StoppingPowerMode fStoppingPowerMode = StoppingPowerMode::Topas;
    StepCalculationMode fStepCalculationMode = StepCalculationMode::MidStep;

    std::vector<double> yhig;
    std::vector<double> yfy;
    std::vector<double> ydy;
    std::vector<std::vector<double>> IonData;
    std::vector<LookupEntry_yS> fLookupTable;

    void LoadLookupTable(const std::string& filename);
    void loadIonData();

    double InterpolateStoppingPower(double izz, double energy);

    G4double GetStepKineticEnergy(G4Step* aStep);

    double CalculateStoppingPower(G4Step* aStep,
                                  G4ParticleDefinition* particleDef,
                                  double izz,
                                  double energyPerNucleon);

    double calculateMicrodosimetricSpectra(double izz,
                                           double iAA,
                                           double ene,
                                           double dEdx);

    void getAparaion(double& CelDiam,
                     double& ene,
                     double& izz,
                     double& ratioc,
                     double& ratioe,
                     double& ratiop,
                     int& ic1,
                     int& ie1,
                     int& ip1);

    double sedmean(double x,
                   double depev,
                   int ic1,
                   int ie1,
                   int ip1,
                   double ratioc,
                   double ratioe,
                   double ratiop,
                   double Apara[]);

    double sedfunc(double x, double depev, double Apara[]);

    static constexpr int nybin = 180;
    const int ROWS = 576;
    const int COLS = 9;
    static constexpr int mparased = 9;
    const int iunit = 2;

    double A9factor = 0;

    G4double CelDiam;
    G4double NucleusRadius;
    G4double DomainRadius;
    G4double BetaRef;

    std::vector<double> eincion = {
        1.0, 2.0, 3.0, 5.0, 7.0, 10.0,
        20.0, 30.0, 50.0, 100.0, 300.0, 999.0
    };

    std::vector<double> cdiamion = {
        0.003, 0.01, 0.03, 0.1,
        0.2, 0.3, 0.5, 1.0
    };

    std::vector<int> izion = {
        1, 2, 6, 10, 14, 26
    };
};

#endif
