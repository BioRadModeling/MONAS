// Scorer for AMF_yS

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

#include "ScoreAMF_yS.hh"
#include "TsParameterManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "TsVBinnedScorer.hh"
#include "G4EmCalculator.hh"
#include "G4Material.hh"
#include "G4Exception.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cstddef>
#include <string>
#include <utility>
#include <map>
#include <numeric>
#include <cctype>

ScoreAMF_yS::ScoreAMF_yS(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
    : TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    yhig.resize(nybin + 1, 0.0);
    yfy.resize(nybin, 0.0);
    ydy.resize(nybin, 0.0); // Initialize ydy
    SetUnit("");

    G4String stoppingPowerCalculation = "Topas";

    if (fPm->ParameterExists(GetFullParmName("StoppingPowerCalculation"))) {
        stoppingPowerCalculation = fPm->GetStringParameter(GetFullParmName("StoppingPowerCalculation"));
    }

    std::string stoppingPowerCalculationLower = stoppingPowerCalculation;
    std::transform(
        stoppingPowerCalculationLower.begin(),
        stoppingPowerCalculationLower.end(),
        stoppingPowerCalculationLower.begin(),
        [](unsigned char c) { return std::tolower(c); }
    );

    if (stoppingPowerCalculationLower == "topas") {
        fStoppingPowerMode = StoppingPowerMode::Topas;
        G4cout << "StoppingPowerCalculation: Topas" << G4endl;
    } else if (stoppingPowerCalculationLower == "externaltable") {
        fStoppingPowerMode = StoppingPowerMode::ExternalTable;
        G4cout << "StoppingPowerCalculation: ExternalTable" << G4endl;
        LoadLookupTable("StoppingPower.txt");
    } else {
        G4cerr << "Invalid StoppingPowerCalculation value: "
               << stoppingPowerCalculation
               << ". Valid options are Topas or ExternalTable. Defaulting to Topas."
               << G4endl;

        fStoppingPowerMode = StoppingPowerMode::Topas;
    }

    G4String stepCalculator = "MidStep";

    if (fPm->ParameterExists(GetFullParmName("StepCalculator"))) {
        stepCalculator = fPm->GetStringParameter(GetFullParmName("StepCalculator"));
    }

    std::string stepCalculatorLower = stepCalculator;
    std::transform(
        stepCalculatorLower.begin(),
        stepCalculatorLower.end(),
        stepCalculatorLower.begin(),
        [](unsigned char c) { return std::tolower(c); }
    );

    if (stepCalculatorLower == "prestep") {
        fStepCalculationMode = StepCalculationMode::PreStep;
        G4cout << "StepCalculator: PreStep" << G4endl;
    } else if (stepCalculatorLower == "midstep") {
        fStepCalculationMode = StepCalculationMode::MidStep;
        G4cout << "StepCalculator: MidStep" << G4endl;
    } else {
        G4cerr << "Invalid StepCalculator value: "
               << stepCalculator
               << ". Valid options are PreStep or MidStep. Defaulting to MidStep."
               << G4endl;

        fStepCalculationMode = StepCalculationMode::MidStep;
    }

    loadIonData();
    InstantiateSubScorer("AMFDose", outFileName, "AMFdose");
    
    DomainRadius = fPm->GetDoubleParameter(GetFullParmName("DomainRadius"), "Length");
    NucleusRadius = fPm->GetDoubleParameter(GetFullParmName("NucleusRadius"), "Length");
    BetaRef = fPm->GetDoubleParameter(GetFullParmName("BetaRef"), "perDoseSquare");
   
    G4cout << "NucleusRadius: " << NucleusRadius / um << " um" << G4endl;
    G4cout << "DomainRadius: " << DomainRadius / um << " um" << G4endl;
    G4cout << "BetaRef " << BetaRef/ (1. / (gray * gray)) << " /Gy^2" << G4endl;

    // Perform unit conversions 
    DomainRadius = DomainRadius / um;   // Convert DomainRadius from mm to um 
    NucleusRadius = NucleusRadius / um; // Convert NucleusRadius from mm to um 
    BetaRef = BetaRef / (1. / (gray * gray)) ;             // Convert BetaRef to /Gy2
    CelDiam = 2.0 * DomainRadius;
}

ScoreAMF_yS::~ScoreAMF_yS() {}

void ScoreAMF_yS::LoadLookupTable(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double energy;
        std::vector<double> stoppingPowers(18); // Stopping powers for Z = 1:18
        if (iss >> energy) {
            for (int i_spbin = 0; i_spbin < 18; ++i_spbin) iss >> stoppingPowers[i_spbin];
            fLookupTable.emplace_back(energy, stoppingPowers);
        }
    }
}

double ScoreAMF_yS::InterpolateStoppingPower(double izz, double energy) {
    int atomicNumberIndex = std::max(1, std::min(static_cast<int>(izz), 18)) - 1;

    if (izz != atomicNumberIndex + 1) {
        std::cerr << "Atomic Number out of valid range. Using closest valid value: " << atomicNumberIndex + 1 << std::endl;
    }

    if (fLookupTable.empty()) {
        std::cerr << "Lookup table is empty." << std::endl;
        return 0.0;
    }

    double minEnergy = fLookupTable.front().energy;
    double maxEnergy = fLookupTable.back().energy;
    double clampedEnergy = std::max(minEnergy, std::min(energy, maxEnergy));

    if (energy != clampedEnergy) {
        std::cerr << "Energy out of range. Using closest valid value: " << clampedEnergy << std::endl;
    }

    auto it = std::lower_bound(fLookupTable.begin(), fLookupTable.end(), clampedEnergy,
                               [](const LookupEntry_yS& entry, const double& energyValue) { return entry.energy < energyValue; });

    if (it != fLookupTable.end() && it->energy == clampedEnergy) {
        return it->stoppingPowers[atomicNumberIndex];
    }

    if (it == fLookupTable.end() || it == fLookupTable.begin()) {
        std::cerr << "Unexpected error: energy index calculation went out of bounds." << std::endl;
        return 0.0;
    }

    auto upper = it;
    auto lower = std::prev(it);

    double x0 = lower->energy, y0 = lower->stoppingPowers[atomicNumberIndex];
    double x1 = upper->energy, y1 = upper->stoppingPowers[atomicNumberIndex];
    double interpolatedStoppingPower = y0 + (y1 - y0) * (clampedEnergy - x0) / (x1 - x0);
    return interpolatedStoppingPower;
}


G4double ScoreAMF_yS::GetStepKineticEnergy(G4Step* aStep) {
    G4double KE_prestep = aStep->GetPreStepPoint()->GetKineticEnergy();

    if (fStepCalculationMode == StepCalculationMode::PreStep) {
        return KE_prestep;
    }

    G4double KE_poststep = aStep->GetPostStepPoint()->GetKineticEnergy();

    return (KE_prestep + KE_poststep) / 2;
}

double ScoreAMF_yS::CalculateStoppingPower(G4Step* aStep,
                                            G4ParticleDefinition* particleDef,
                                            double izz,
                                            double energyPerNucleon) {
    double dEdx = 0.0;

    if (fStoppingPowerMode == StoppingPowerMode::ExternalTable) {
        dEdx = InterpolateStoppingPower(izz, energyPerNucleon);
    } else {
        const G4Material* materialStep = aStep->GetPreStepPoint()->GetMaterial();

        G4double KE_step = GetStepKineticEnergy(aStep);

        G4EmCalculator emCal;

        G4double dEdxG4 = emCal.ComputeElectronicDEDX(
            KE_step,
            particleDef,
            materialStep
        );

        dEdx = dEdxG4 / (keV / um);
    }

    return dEdx;
}

G4bool ScoreAMF_yS::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    G4double edep = aStep->GetTotalEnergyDeposit();

    if (edep <= 0.) {
        return false;
    }

    G4ParticleDefinition* particleDef = aStep->GetTrack()->GetDefinition();

    G4double izz = particleDef->GetAtomicNumber();
    G4double iAA = particleDef->GetAtomicMass();

    if (izz < 1 || izz > 18) {
        return false;
    }

    G4double KE_step = GetStepKineticEnergy(aStep);
    G4double ene = KE_step / iAA;

    if (ene < 0.025) {
        return false;
    }

    G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();

    ResolveSolid(aStep);

    G4double dose = edep / (density * fSolid->GetCubicVolume());

    double dEdx = CalculateStoppingPower(
        aStep,
        particleDef,
        izz,
        ene
    );

    double yS_calc = calculateMicrodosimetricSpectra(izz, iAA, ene, dEdx);

    double yS = yS_calc * dose;
    yS *= aStep->GetPreStepPoint()->GetWeight();

    AccumulateHit(aStep, yS);

    return true;
}

void ScoreAMF_yS::loadIonData() {
    IonData.resize(ROWS, std::vector<double>(COLS));
    std::ifstream file("tsed.dat");
    if (!file.is_open()) {
        std::cerr << "Failed to open tsed.dat for reading." << std::endl;
        return;
    }

    std::string line;
    int row = 0;

    while (std::getline(file, line) && row < ROWS) {
        std::istringstream iss(line);
        double value;
        int col = 0;

        while (iss >> value && col < COLS) {
            IonData[row][col] = value;
            col++;
        }
        row++;
    }
    file.close();
}

void ScoreAMF_yS::getAparaion(double& CelDiam, double& ene, double& izz, double& ratioc, double& ratioe, double& ratiop, int& ic1, int& ie1, int& ip1) {
    double erg_AA = ene;
    double CD = std::abs(CelDiam);

   int ic = 0;
    for (ic = 1; ic <= 9; ++ic) {
        if (cdiamion[ic - 1] >= CD) {
            break;
        }
    }

    if (ic == 1) {
        ic1 = 1;
        ratioc = 0.0;
    } else {
        ic1 = ic - 1;
        ratioc = std::min(1.0, (std::log10(CD) - std::log10(cdiamion[ic1-1])) / (std::log10(cdiamion[ic1]) - std::log10(cdiamion[ic1-1])));
    }

    int ie = 0;
    for (ie = 1; ie <= 12; ++ie) {
        if (eincion[ie-1] >= erg_AA) {
            break;
        }
    }

    if (ie == 1) {
        ie1 = 1;
        ratioe = 0.0;
    } else {
        ie1 = ie - 1;
        ratioe = std::min(1.0, (std::log10(erg_AA) - std::log10(eincion[ie1-1])) / (std::log10(eincion[ie1]) - std::log10(eincion[ie1-1])));
    }

    int ip = 1;
    for (ip = 2; ip <= 6; ++ip) {
        if (izion[ip-1] >= izz) {
            break;
        }
    }
    ip1 = ip - 1;
    ratiop = std::min(1.0, static_cast<double>(izz - izion[ip1-1]) / (izion[ip1] - izion[ip1-1]));
}

double ScoreAMF_yS::sedmean(double x, double depev, int ic1, int ie1, int ip1, double ratioc, double ratioe, double ratiop, double Apara[]) {
    double sedmean = 0.0;
    double A9 = 0.0;

    for (int ip = ip1; ip <= ip1 + 1; ip++) {
        double Rp = (ip == ip1) ? (1.0 - ratiop) : ratiop;
        for (int ie = ie1; ie <= ie1 + 1; ie++) {
            double Re = (ie == ie1) ? (1.0 - ratioe) : ratioe;
            for (int ic = ic1; ic <= ic1 + 1; ic++) {
                double Rc = (ic == ic1) ? (1.0 - ratioc) : ratioc;
                int index = ((ip-1) * 96) + ((ie-1) * 8) + (ic-1);
                for (int i = 0; i < mparased; i++) {
                    Apara[i] = IonData[index][i];
 //                   std::cout << "ip1: "<< ip1 <<", ic1: " << ic1 <<", ie1: " << ie1 << std::endl;
                }

                double wei = Rp * Re * Rc;
                double sedfuncResult = sedfunc(x, depev, Apara);
           
 
                sedmean += sedfuncResult * wei;
                A9 += Apara[8] * wei;
            }
        }
    }
    Apara[8] = A9;
//    std::cout << "Apara0: " <<Apara[0] << ", Apara8: " << Apara[8] << std::endl;
    return sedmean;
}

double ScoreAMF_yS::sedfunc(double x, double depev, double Apara[]) {
    double getfirst = 0.0, getsecond = 0.0, getthird = 0.0;
    if (Apara[0] > 0.0) {
        double tmp;
        if (depev == 0.0) {
            tmp = std::min(50.0, std::pow(std::abs(x - Apara[1]), Apara[2]) / (2 * Apara[1]));
            getfirst = Apara[0] * std::exp(-tmp);
        } else {
            double cst1 = depev / Apara[8];
            tmp = std::min(50.0, Apara[1] * (x - cst1 * Apara[2]));
            getfirst = Apara[0] * x / (std::exp(tmp) + 1) * (2.0 / std::pow(cst1 * Apara[2], 2));
        }
    }

    if (Apara[3] > 0.0) {
        double tmp = std::min(50.0, std::pow(std::abs(x - Apara[4]), Apara[5]) / (2 * Apara[4]));
        getsecond = Apara[3] * std::exp(-tmp);
    }

    if (Apara[6] > 0.0) {
        getthird = Apara[6] / (Apara[7] - 1.0) * std::pow((Apara[7] - 1.0) / Apara[7], x);
    }

    double sedfunc = getfirst + getsecond + getthird;
    if (sedfunc < 1.0e-10) {
        sedfunc = 0.0;
    }

    return sedfunc;
}

double ScoreAMF_yS::calculateMicrodosimetricSpectra(double izz, double iAA, double ene, double dEdx) {
    double factor;
    double unitconv = 1.0;
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
    double Apara[mparased] = {0.0};

    yhig.resize(nybin + 1);
    yfy.resize(nybin);
    ydy.resize(nybin); // Resize ydy
    double ypower = -2.0;
    const double ystep = 0.05;

    for (size_t i = 0; i < yhig.size(); ++i) {
        yhig[i] = std::pow(10.0, ypower);
        ypower += ystep;
    }

    if (iunit <= 1) {
        unitconv = 1.0;
    } else if (iunit == 2) {
        unitconv = 1.0e-3 * (2.0 / 3.0 * CelDiam);
    } else if (iunit == 3) {
        unitconv = 4.0 / 3.0 * M_PI * std::pow(CelDiam / 2.0, 3) * 1.0e-15 / 1.602e-13;
    }

    int ic1, ie1, ip1;
    double ratioc, ratioe, ratiop;
    std::vector<std::pair<double, double>> microdosimetricSpectra(nybin);
    double erg = ene * iAA;
    double depev = std::min(dEdx * CelDiam * 1.0e3, erg * 1.0e6);
    getAparaion(CelDiam, ene, izz, ratioc, ratioe, ratiop, ic1, ie1, ip1);
    sedmean(1.0, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
    factor = (iunit == 0) ? 1.0 : 1.0e6 / Apara[8];

    for (size_t i = 0; i < nybin; ++i) {
        double ymid = (yhig[i] + yhig[i + 1]) / 2.0;
        double ywid = yhig[i + 1] - yhig[i];
        double eventmid = ymid * factor * unitconv;
        yfy[i] = ymid * sedmean(eventmid, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
        ydy[i] = yfy[i] * ymid; // Calculate ydy
        sum0 += yfy[i] * ywid / ymid;
        sum1 += yfy[i] * ywid;
        sum2 += yfy[i] * ywid * ymid;
    }

    // Calculate the bins per decade and normalization factor
    double binsperDecade = nybin / (std::log10(yhig.back() / yhig[0]));
    double normalization_factor = (binsperDecade / std::log(10)) / std::accumulate(ydy.begin(), ydy.end(), 0.0);

    // Apply normalization to ydy
    for (size_t i = 0; i < ydy.size(); ++i) {
        ydy[i] *= normalization_factor;
    }

//    double LinealEnergy_Dose = sum2 / sum1;
    double LinealEnergy_Freq = sum1 / sum0;
    // yS calculation
//    double y0 = std::pow(NucleusRadius / DomainRadius, 2) / std::sqrt(BetaRef * (1 + std::pow(NucleusRadius / DomainRadius, 2)));
    double y0 = (M_PI * DomainRadius * std::pow(NucleusRadius, 2)) / (std::sqrt(BetaRef * (std::pow(DomainRadius, 2) + std::pow(NucleusRadius, 2))) * 0.16022);
 //   G4cout << "y0: " << y0 << G4endl;
    std::vector<double> Z(nybin);
    for (size_t i = 0; i < nybin; ++i) {
        double F2 = (yhig[i] + yhig[i + 1]) / 2.0; // Bin center
        Z[i] = 1 - std::exp(-std::pow(F2, 2) / std::pow(y0, 2));
    }
    double sumNumerator = 0.0;
    double sumDenominator = 0.0;

    for (size_t i = 0; i < nybin; ++i) {
        sumNumerator += yfy[i] * Z[i];
        sumDenominator += yfy[i];
    }

    double LinealEnergyS = 0.0;
    if (sumDenominator > 0.0) {
        LinealEnergyS = ((sumNumerator / sumDenominator) / LinealEnergy_Freq ) * std::pow(y0, 2);
    }
//    double LinealEnergyS = LinealEnergyZ / ( M_PI * rho * std::pow(DomainRadius, 2) );

    // Output LinealEnergyZ
 /*   G4cout << "dEdx: " << dEdx << ", Z= " << izz << ", ene = " << ene << G4endl;
    G4cout << "Lineal Energy S (yS): " << LinealEnergyS << G4endl; */

    return LinealEnergyS;
}

G4int ScoreAMF_yS::CombineSubScorers()
{
	G4int counter = 0;
	TsVBinnedScorer* AMF_dose = dynamic_cast<TsVBinnedScorer*>(GetSubScorer("AMFdose"));
	for (G4int index=0; index < fNDivisions; index++) {
		if (AMF_dose->fFirstMomentMap[index] == 0.) {
			fFirstMomentMap[index] = 0;
			counter++;
		} else {
			fFirstMomentMap[index] = fFirstMomentMap[index] / AMF_dose->fFirstMomentMap[index]  ;
		}
	}

	return counter;
}
