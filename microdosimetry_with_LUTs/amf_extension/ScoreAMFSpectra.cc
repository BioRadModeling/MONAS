// Scorer for AMFSpectra

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

#include "ScoreAMFSpectra.hh"
#include "TsParameterManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "TsVBinnedScorer.hh"
#include "G4EmCalculator.hh"
#include "G4Material.hh"
#include "G4Exception.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
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
#include <limits>

ScoreAMFSpectra::ScoreAMFSpectra(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
    : TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    yhig.resize(nybin + 1, 0.0);
    yfy.resize(nybin, 0.0);
    ydy.resize(nybin, 0.0); // Initialize ydy
    initializeYGrid();
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
    
    DomainRadius = fPm->GetDoubleParameter(GetFullParmName("DomainRadius"), "Length");
    G4cout << "DomainRadius: " << DomainRadius / um << " um" << G4endl;
    DomainRadius = DomainRadius / um;   // Convert DomainRadius from mm to um
    CelDiam = 2.0 * DomainRadius;
}

ScoreAMFSpectra::~ScoreAMFSpectra() {}

void ScoreAMFSpectra::RatioMomentAccumulator::AddSample(
    G4double numerator,
    G4double denominator) {

    if (!std::isfinite(numerator) || !std::isfinite(denominator) ||
        denominator <= 0.0) {
        return;
    }

    ++count;
    numeratorSum += numerator;
    denominatorSum += denominator;
    numeratorSquaredSum += numerator * numerator;
    denominatorSquaredSum += denominator * denominator;
    numeratorDenominatorSum += numerator * denominator;
}

void ScoreAMFSpectra::RatioMomentAccumulator::Absorb(
    const RatioMomentAccumulator& other) {

    count += other.count;
    numeratorSum += other.numeratorSum;
    denominatorSum += other.denominatorSum;
    numeratorSquaredSum += other.numeratorSquaredSum;
    denominatorSquaredSum += other.denominatorSquaredSum;
    numeratorDenominatorSum += other.numeratorDenominatorSum;
}

G4double ScoreAMFSpectra::RatioMomentAccumulator::Mean() const {
    if (!(denominatorSum > 0.0)) {
        return std::numeric_limits<G4double>::quiet_NaN();
    }

    return numeratorSum / denominatorSum;
}

G4double ScoreAMFSpectra::RatioMomentAccumulator::StandardError() const {
    if (count < 2 || !(denominatorSum > 0.0)) {
        return std::numeric_limits<G4double>::quiet_NaN();
    }

    const G4double ratio = Mean();
    if (!std::isfinite(ratio)) {
        return std::numeric_limits<G4double>::quiet_NaN();
    }

    G4double residualSum =
        numeratorSquaredSum
        - 2.0 * ratio * numeratorDenominatorSum
        + ratio * ratio * denominatorSquaredSum;

    if (residualSum < 0.0 && std::abs(residualSum) < 1.0e-12) {
        residualSum = 0.0;
    }
    if (residualSum < 0.0) {
        return std::numeric_limits<G4double>::quiet_NaN();
    }

    const G4double variance =
        (static_cast<G4double>(count) / static_cast<G4double>(count - 1))
        * residualSum / (denominatorSum * denominatorSum);

    if (variance < 0.0) {
        return std::numeric_limits<G4double>::quiet_NaN();
    }

    return std::sqrt(variance);
}

ScoreAMFSpectra::SpectrumDistributionMoment
ScoreAMFSpectra::ComputeDistributionMoments(
    const std::vector<G4double>& yCenters,
    const std::vector<G4double>& densityValues,
    const std::string& distribution,
    G4double meanStandardError) const {

    SpectrumDistributionMoment moments;
    moments.distribution = distribution;
    moments.meanStandardErrorKeVPerUm = meanStandardError;

    if (yCenters.size() != densityValues.size() || yCenters.size() < 2) {
        moments.meanKeVPerUm = std::numeric_limits<G4double>::quiet_NaN();
        moments.varianceKeV2PerUm2 = std::numeric_limits<G4double>::quiet_NaN();
        moments.stdevKeVPerUm = std::numeric_limits<G4double>::quiet_NaN();
        moments.skewness = std::numeric_limits<G4double>::quiet_NaN();
        return moments;
    }

    const G4double normalizationFactor = std::log(10.0) * yStep;
    G4double totalProbability = 0.0;

    for (size_t i = 0; i < yCenters.size(); ++i) {
        const G4double y = yCenters[i];
        const G4double density = densityValues[i];

        if (!(y > 0.0) || !std::isfinite(density) || density < 0.0) {
            continue;
        }

        totalProbability += normalizationFactor * y * density;
    }

    if (!(totalProbability > 0.0)) {
        moments.meanKeVPerUm = std::numeric_limits<G4double>::quiet_NaN();
        moments.varianceKeV2PerUm2 = std::numeric_limits<G4double>::quiet_NaN();
        moments.stdevKeVPerUm = std::numeric_limits<G4double>::quiet_NaN();
        moments.skewness = std::numeric_limits<G4double>::quiet_NaN();
        return moments;
    }

    G4double meanNumerator = 0.0;
    for (size_t i = 0; i < yCenters.size(); ++i) {
        const G4double y = yCenters[i];
        const G4double density = densityValues[i];

        if (!(y > 0.0) || !std::isfinite(density) || density < 0.0) {
            continue;
        }

        const G4double probability = normalizationFactor * y * density;
        meanNumerator += probability * y;
    }
    moments.meanKeVPerUm = meanNumerator / totalProbability;

    G4double varianceNumerator = 0.0;
    for (size_t i = 0; i < yCenters.size(); ++i) {
        const G4double y = yCenters[i];
        const G4double density = densityValues[i];

        if (!(y > 0.0) || !std::isfinite(density) || density < 0.0) {
            continue;
        }

        const G4double probability = normalizationFactor * y * density;
        const G4double delta = y - moments.meanKeVPerUm;
        varianceNumerator += probability * delta * delta;
    }

    moments.varianceKeV2PerUm2 = varianceNumerator / totalProbability;
    if (moments.varianceKeV2PerUm2 < 0.0 &&
        std::abs(moments.varianceKeV2PerUm2) < 1.0e-12) {
        moments.varianceKeV2PerUm2 = 0.0;
    }

    if (moments.varianceKeV2PerUm2 < 0.0) {
        moments.stdevKeVPerUm = std::numeric_limits<G4double>::quiet_NaN();
        moments.skewness = std::numeric_limits<G4double>::quiet_NaN();
        return moments;
    }

    moments.stdevKeVPerUm = std::sqrt(moments.varianceKeV2PerUm2);
    if (moments.stdevKeVPerUm == 0.0) {
        moments.skewness = 0.0;
        return moments;
    }

    G4double skewnessNumerator = 0.0;
    for (size_t i = 0; i < yCenters.size(); ++i) {
        const G4double y = yCenters[i];
        const G4double density = densityValues[i];

        if (!(y > 0.0) || !std::isfinite(density) || density < 0.0) {
            continue;
        }

        const G4double probability = normalizationFactor * y * density;
        const G4double standardized =
            (y - moments.meanKeVPerUm) / moments.stdevKeVPerUm;
        skewnessNumerator +=
            probability * standardized * standardized * standardized;
    }

    moments.skewness = skewnessNumerator / totalProbability;
    return moments;
}

void ScoreAMFSpectra::initializeYGrid() {
    yhig.resize(nybin + 1);

    double ypower = yPowerMin;
    for (size_t i = 0; i < yhig.size(); ++i) {
        yhig[i] = std::pow(10.0, ypower);
        ypower += yStep;
    }
}

G4int ScoreAMFSpectra::GetCurrentEventId() const {
    G4RunManager* runManager = G4RunManager::GetRunManager();
    if (!runManager) {
        return -1;
    }

    const G4Event* event = runManager->GetCurrentEvent();
    if (!event) {
        return -1;
    }

    return event->GetEventID();
}

void ScoreAMFSpectra::FlushCurrentEventMoments() {
    for (const auto& eventPair : currentEventMoments) {
        const G4int binIndex = eventPair.first;
        const EventMomentContribution& contribution = eventPair.second;

        doseMeanStats[binIndex].AddSample(
            contribution.yDWeightedNumerator,
            contribution.doseDenominator);
        frequencyMeanStats[binIndex].AddSample(
            contribution.yFWeightedNumerator,
            contribution.yFWeightedDenominator);
    }

    currentEventMoments.clear();
}

void ScoreAMFSpectra::AccumulateCurrentEventMoments(
    G4int binIndex,
    G4double dose,
    const std::vector<std::pair<double, double>>& microdosimetricSpectra) {

    if (!(dose > 0.0)) {
        return;
    }

    const double normalizationFactor = std::log(10.0) * yStep;
    double yDStep = 0.0;
    double inverseYFStep = 0.0;

    for (const auto& bin : microdosimetricSpectra) {
        const double y = bin.first;
        const double yd = bin.second;

        if (!(y > 0.0) || !std::isfinite(yd) || yd < 0.0) {
            continue;
        }

        yDStep += normalizationFactor * y * yd;
        inverseYFStep += normalizationFactor * yd / y;
    }

    if (!(yDStep > 0.0) || !(inverseYFStep > 0.0)) {
        return;
    }

    EventMomentContribution& contribution = currentEventMoments[binIndex];
    contribution.yDWeightedNumerator += dose * yDStep;
    contribution.doseDenominator += dose;
    contribution.yFWeightedNumerator += dose;
    contribution.yFWeightedDenominator += dose * inverseYFStep;
}

void ScoreAMFSpectra::ComputeFinalSpectrumMoments() {
    initializeYGrid();
    spectrumMomentSummaries.clear();

    std::vector<G4double> yCenters(nybin, 0.0);
    for (size_t i = 0; i < yCenters.size(); ++i) {
        yCenters[i] = (yhig[i] + yhig[i + 1]) / 2.0;
    }

    const G4double normalizationFactor = std::log(10.0) * yStep;

    for (const auto& binPair : totalSpectra) {
        const G4int binIndex = binPair.first;
        const std::vector<G4double>& ydSpectrum = binPair.second;

        if (ydSpectrum.size() != yCenters.size()) {
            continue;
        }

        G4double doseProbabilityNorm = 0.0;
        for (size_t i = 0; i < ydSpectrum.size(); ++i) {
            const G4double yd = ydSpectrum[i];
            if (std::isfinite(yd) && yd >= 0.0) {
                doseProbabilityNorm += normalizationFactor * yd;
            }
        }

        if (!(doseProbabilityNorm > 0.0)) {
            continue;
        }

        std::vector<G4double> doseDensity(yCenters.size(), 0.0);
        G4double inverseYF = 0.0;

        for (size_t i = 0; i < yCenters.size(); ++i) {
            const G4double y = yCenters[i];
            const G4double yd = ydSpectrum[i] / doseProbabilityNorm;

            if (!(y > 0.0) || !std::isfinite(yd) || yd < 0.0) {
                continue;
            }

            doseDensity[i] = yd / y;
            inverseYF += normalizationFactor * yd / y;
        }

        if (!(inverseYF > 0.0)) {
            continue;
        }

        const G4double yF = 1.0 / inverseYF;
        std::vector<G4double> frequencyDensity(yCenters.size(), 0.0);

        for (size_t i = 0; i < yCenters.size(); ++i) {
            const G4double y = yCenters[i];
            if (y > 0.0) {
                frequencyDensity[i] = doseDensity[i] * yF / y;
            }
        }

        G4double frequencySem = std::numeric_limits<G4double>::quiet_NaN();
        auto frequencyStatsIt = frequencyMeanStats.find(binIndex);
        if (frequencyStatsIt != frequencyMeanStats.end()) {
            frequencySem = frequencyStatsIt->second.StandardError();
        }

        G4double doseSem = std::numeric_limits<G4double>::quiet_NaN();
        auto doseStatsIt = doseMeanStats.find(binIndex);
        if (doseStatsIt != doseMeanStats.end()) {
            doseSem = doseStatsIt->second.StandardError();
        }

        SpectrumMomentSummary summary;
        summary.frequency = ComputeDistributionMoments(
            yCenters,
            frequencyDensity,
            "frequency",
            frequencySem);
        summary.dose = ComputeDistributionMoments(
            yCenters,
            doseDensity,
            "dose",
            doseSem);
        spectrumMomentSummaries[binIndex] = summary;
    }
}

void ScoreAMFSpectra::LoadLookupTable(const std::string& filename) {
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

double ScoreAMFSpectra::InterpolateStoppingPower(double izz, double energy) {
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
                               [](const LookupEntryAMFSpectra& entry, const double& energyValue) { return entry.energy < energyValue; });

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


G4double ScoreAMFSpectra::GetStepKineticEnergy(G4Step* aStep) {
    G4double KE_prestep = aStep->GetPreStepPoint()->GetKineticEnergy();

    if (fStepCalculationMode == StepCalculationMode::PreStep) {
        return KE_prestep;
    }

    G4double KE_poststep = aStep->GetPostStepPoint()->GetKineticEnergy();

    return (KE_prestep + KE_poststep) / 2;
}

double ScoreAMFSpectra::CalculateStoppingPower(G4Step* aStep,
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

G4bool ScoreAMFSpectra::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
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
    G4int binIndex = GetIndex(aStep);

    const G4int eventId = GetCurrentEventId();
    if (eventId != currentEventId) {
        FlushCurrentEventMoments();
        currentEventId = eventId;
    }

    double dEdx = CalculateStoppingPower(
        aStep,
        particleDef,
        izz,
        ene
    );

    auto microdosimetricSpectra = calculateMicrodosimetricFunction(izz, iAA, ene, dEdx);

    if (totalSpectra.find(binIndex) == totalSpectra.end()) {
        totalSpectra[binIndex] = std::vector<G4double>(microdosimetricSpectra.size(), 0.0);
    } else if (totalSpectra[binIndex].size() < microdosimetricSpectra.size()) {
        totalSpectra[binIndex].resize(microdosimetricSpectra.size(), 0.0);
    }

    for (size_t i = 0; i < microdosimetricSpectra.size(); ++i) {
        totalSpectra[binIndex][i] += microdosimetricSpectra[i].second * dose;
    }

    cumulativeDose[binIndex] += dose;
    AccumulateCurrentEventMoments(binIndex, dose, microdosimetricSpectra);

    return true;
}



void ScoreAMFSpectra::loadIonData() {

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



// In getAparaion function

void ScoreAMFSpectra::getAparaion(const double& CelDiam, const double& ene, const int& iAA, const int& izz, double& ratioc, double& ratioe, double& ratiop, int& ic1, int& ie1, int& ip1) {

    double erg = ene * iAA;

    int modifiedIzz = (izz > 26) ? 26 : izz;

    double erg_AA = erg / iAA;

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





double ScoreAMFSpectra::sedfunc(double x, double depev, const double Apara[], size_t size) {

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



double ScoreAMFSpectra::sedmean(double x, double depev, int ic1, int ie1, int ip1, double ratioc, double ratioe, double ratiop, double Apara[]) {

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

                double sedfuncResult = sedfunc(x, depev, Apara, mparased);

 

                sedmean += sedfuncResult * wei;

                A9 += Apara[8] * wei;

            }

        }



    

    }

    Apara[8] = A9;
//    std::cout << "Apara0: " <<Apara[0] << ", Apara8: " << Apara[8] << std::endl;
    return sedmean;
}


std::vector<std::pair<double, double>> ScoreAMFSpectra::calculateMicrodosimetricFunction(double izz, double iAA, double ene, double dEdx) {
    double unitconv, factor;
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
    double Apara[mparased] = {0.0};

    yhig.resize(nybin + 1);
    yfy.resize(nybin);
    ydy.resize(nybin); // Resize ydy

    initializeYGrid();
  
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

    getAparaion(CelDiam, ene, iAA, izz, ratioc, ratioe, ratiop, ic1, ie1, ip1);
    sedmean(1.0, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
    factor = (iunit == 0) ? 1.0 : 1.0e6 / Apara[8];
//    double tmp = sedmean(1.0, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop);
//    std::cout << "factor: " << factor << std::endl;
//    std::cout << "Apara[8]: " << Apara[8] << std::endl;

    for (size_t i = 0; i < nybin; ++i) {
        double ymid = (yhig[i] + yhig[i + 1]) / 2.0;
//        std::cout << "ymid_bin: " << ymid << std::endl;
        double ywid = yhig[i + 1] - yhig[i];
        double eventmid = ymid * factor * unitconv;
//	std::cout << "eventmid: " << eventmid << std::endl;
//	std::cout << "unit conv: " << unitconv << std::endl;
//	std::cout << "factor: " << factor << std::endl;
//	std::cout << "ymid: " << ymid << std::endl;
        yfy[i] = ymid * sedmean(eventmid, depev, ic1, ie1, ip1, ratioc, ratioe, ratiop, Apara);
//        std::cout << "yfy_bin: " << ydy[i] << std::endl;
        ydy[i] = yfy[i] * ymid; // Calculate ydy
        sum0 += yfy[i] * ywid / ymid;
        sum1 += yfy[i] * ywid;
        sum2 += yfy[i] * ywid * ymid;
    }

    // Calculate the bins per decade and normalization factor
    double binsperDecade = nybin / (std::log10(yhig.back() / yhig[0]));
    double normalization_factor = (binsperDecade / std::log(10)) / std::accumulate(ydy.begin(), ydy.end(), 0.0);
//    std::cout << "yF: " << sum1/sum0 << std::endl;
//    std::cout << "yD: " << sum2/sum1 << std::endl;
    // Apply normalization to ydy
    for (size_t i = 0; i < ydy.size(); ++i) {
        ydy[i] *= normalization_factor;
    }

    for (size_t i = 0; i + 1 < yhig.size(); ++i) {
        double ymid = (yhig[i] + yhig[i + 1]) / 2.0;
        microdosimetricSpectra[i] = std::make_pair(ymid, ydy[i]); // Use ydy instead of yfy
    }
 /*   	double LinealEnergy_Freq = sum1/sum0;
	double LinealEnergy_Dose = sum2/sum1;
	    // Print the expectation value along with other relevant information before returning the spectra
    G4cout << "Atomic Number: " << izz
           << ", Kinetic Energy: " << ene
           << ", LET: " << dEdx
           << ", yF: " << LinealEnergy_Freq << ", yD: " << LinealEnergy_Dose << G4endl;
*/
    return microdosimetricSpectra;
}

void ScoreAMFSpectra::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
    ScoreAMFSpectra* worker = dynamic_cast<ScoreAMFSpectra*>(workerScorer);

    if (!worker) {
        G4cerr << "Error: Incorrect scorer type passed to AbsorbResultsFromWorkerScorer." << G4endl;
        return;
    }

    worker->FlushCurrentEventMoments();

    for (const auto& workerPair : worker->totalSpectra) {
        G4int binIndex = workerPair.first;
        const std::vector<G4double>& workerSpectra = workerPair.second;

        if (totalSpectra.find(binIndex) == totalSpectra.end()) {
            totalSpectra[binIndex] = std::vector<G4double>(workerSpectra.size(), 0.0);
        }

        for (size_t i = 0; i < workerSpectra.size(); ++i) {
            totalSpectra[binIndex][i] += workerSpectra[i];
        }
    }

    for (const auto& workerPair : worker->cumulativeDose) {

        G4int binIndex = workerPair.first;

        G4double workerDose = workerPair.second;


        if (cumulativeDose.find(binIndex) == cumulativeDose.end()) {
            cumulativeDose[binIndex] = 0.0;
        }
        cumulativeDose[binIndex] += workerDose;
    }

    for (const auto& workerPair : worker->frequencyMeanStats) {
        frequencyMeanStats[workerPair.first].Absorb(workerPair.second);
    }

    for (const auto& workerPair : worker->doseMeanStats) {
        doseMeanStats[workerPair.first].Absorb(workerPair.second);
    }
}

void ScoreAMFSpectra::UserHookForEndOfRun() {
    FlushCurrentEventMoments();

    for (auto& binSpectraPair : totalSpectra) {
        G4int binIndex = binSpectraPair.first;
        G4double binDose = cumulativeDose[binIndex];
        for (G4double& value : binSpectraPair.second) {
            value /= binDose;
        }
    }

    G4cout << "Cumulative Dose for each bin index:" << G4endl;
    for (const auto& dosePair : cumulativeDose) {
        G4int binIndex = dosePair.first;
        G4double binDose = dosePair.second;
        G4cout << "Bin Index " << binIndex << ": " << binDose << " Gy" << G4endl;
    }
    ComputeFinalSpectrumMoments();
    OutputFinalSpectra();
    OutputFinalSpectrumMoments();
}

void ScoreAMFSpectra::OutputFinalSpectra() {
    initializeYGrid();

    G4String out_file_fn = fPm->GetStringParameter(GetFullParmName("OutputFile")) + "_MicrodosimetricSpectra.csv";
    FILE* out_file = fopen(out_file_fn.c_str(), "w");

    if (!out_file) {
        std::cerr << "Failed to open output file" << std::endl;
        return;
    }

    // Write header
    fprintf(out_file, "# Microdosimetric Spectra Output\n");
    fprintf(out_file, "# Format: x, y, z, bin_1, bin_2, ..., bin_N\n");

    // Write bin center headers
    fprintf(out_file, "x,y,z");
    for (size_t i = 0; i < yhig.size() - 1; ++i) {
        double ymid = (yhig[i] + yhig[i + 1]) / 2.0;
        fprintf(out_file, ",%.5e", ymid);
    }
    fprintf(out_file, "\n");

    // Loop over all bins (voxels)
    for (const auto& binPair : totalSpectra) {
        G4int binIndex = binPair.first;
        const auto& spectra = binPair.second;

        // Get voxel coordinates (x, y, z)
        G4int ix = GetBin(binIndex, 0);
        G4int iy = GetBin(binIndex, 1);
        G4int iz = GetBin(binIndex, 2);

        fprintf(out_file, "%d,%d,%d", ix, iy, iz);

        for (const auto& value : spectra) {
            fprintf(out_file, ",%.5e", value);
        }

        fprintf(out_file, "\n");
    }

    fclose(out_file);
    G4cout << "Microdosimetric spectra output written to " << out_file_fn << G4endl;
}

void ScoreAMFSpectra::OutputFinalSpectrumMoments() {
    G4String out_file_fn =
        fPm->GetStringParameter(GetFullParmName("OutputFile")) +
        "_MicrodosimetricMoments.csv";
    FILE* out_file = fopen(out_file_fn.c_str(), "w");

    if (!out_file) {
        std::cerr << "Failed to open AMF spectrum moments output file" << std::endl;
        return;
    }

    fprintf(
        out_file,
        "x,y,z,distribution,mean_keV_per_um,variance_keV2_per_um2,"
        "stdev_keV_per_um,mean_standard_error_keV_per_um,skewness\n");

    const auto writeMomentRow =
        [this, out_file](G4int binIndex, const SpectrumDistributionMoment& moments) {
            G4int ix = GetBin(binIndex, 0);
            G4int iy = GetBin(binIndex, 1);
            G4int iz = GetBin(binIndex, 2);

            fprintf(
                out_file,
                "%d,%d,%d,%s,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                ix,
                iy,
                iz,
                moments.distribution.c_str(),
                moments.meanKeVPerUm,
                moments.varianceKeV2PerUm2,
                moments.stdevKeVPerUm,
                moments.meanStandardErrorKeVPerUm,
                moments.skewness);
        };

    for (const auto& summaryPair : spectrumMomentSummaries) {
        const G4int binIndex = summaryPair.first;
        const SpectrumMomentSummary& summary = summaryPair.second;
        writeMomentRow(binIndex, summary.frequency);
        writeMomentRow(binIndex, summary.dose);
    }

    fclose(out_file);
    G4cout << "Microdosimetric moments output written to " << out_file_fn << G4endl;
}
