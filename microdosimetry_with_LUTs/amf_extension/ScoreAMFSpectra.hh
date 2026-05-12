#ifndef ScoreAMFSpectra_hh
#define ScoreAMFSpectra_hh

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
#include <cmath>

class LookupEntryAMFSpectra {
public:
    double energy;
    std::vector<double> stoppingPowers;

    LookupEntryAMFSpectra(double e, const std::vector<double>& sp)
        : energy(e), stoppingPowers(sp) {}
};

class ScoreAMFSpectra : public TsVBinnedScorer {
public:
    ScoreAMFSpectra(TsParameterManager* pM,
                    TsMaterialManager* mM,
                    TsGeometryManager* gM,
                    TsScoringManager* scM,
                    TsExtensionManager* eM,
                    G4String scorerName,
                    G4String quantity,
                    G4String outFileName,
                    G4bool isSubScorer = false);

    virtual ~ScoreAMFSpectra();

    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) override;
    void UserHookForEndOfRun();
    void OutputFinalSpectra();

    double InterpolateStoppingPower(double izz, double ene);

    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* history) override;

    std::vector<std::pair<double, double>> calculateMicrodosimetricFunction(double izz,
                                                                            double iAA,
                                                                            double ene,
                                                                            double dEdx);

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

    static constexpr int nybin = 400;
    static constexpr int mparased = 9;
    static constexpr int iunit = 2;
    static constexpr int ROWS = 576;
    static constexpr int COLS = 9;
    static constexpr double yPowerMin = -3.0;
    static constexpr double yStep = 0.02;

    std::vector<double> yhig;
    std::vector<double> yfy;
    std::vector<double> ydy;

    std::vector<std::vector<double>> IonData;
    std::vector<LookupEntryAMFSpectra> fLookupTable;

    std::map<G4int, std::vector<G4double>> totalSpectra;
    std::map<G4int, G4double> cumulativeDose;

    const std::vector<double> eincion = {
        1.0, 2.0, 3.0, 5.0, 7.0, 10.0,
        20.0, 30.0, 50.0, 100.0, 300.0, 999.0
    };

    const std::vector<double> cdiamion = {
        0.003, 0.01, 0.03, 0.1,
        0.2, 0.3, 0.5, 1.0
    };

    const std::vector<int> izion = {
        1, 2, 6, 10, 14, 26
    };

    G4double CelDiam;
    G4double DomainRadius;

    double totalDose = 0.0;
    double A9factor = 0.0;

    void LoadLookupTable(const std::string& filename);
    void loadIonData();
    void initializeYGrid();

    G4double GetStepKineticEnergy(G4Step* aStep);

    double CalculateStoppingPower(G4Step* aStep,
                                  G4ParticleDefinition* particleDef,
                                  double izz,
                                  double energyPerNucleon);

    double sedmean(double x,
                   double depev,
                   int ic1,
                   int ie1,
                   int ip1,
                   double ratioc,
                   double ratioe,
                   double ratiop,
                   double Apara[]);

    double sedfunc(double x,
                   double depev,
                   const double Apara[],
                   size_t size);

    void getAparaion(const double& CelDiam,
                     const double& ene,
                     const int& iAA,
                     const int& izz,
                     double& ratioc,
                     double& ratioe,
                     double& ratiop,
                     int& ic1,
                     int& ie1,
                     int& ip1);
};

#endif
