#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "AtLookup.h"
#include "PhaseSpaceReader.h"

enum class AtParticle {
    Proton,
    Carbon
};

struct AtDiagnosticRecord {
    std::string energyBandMeVPerU;
    std::size_t rowCount{0};
    std::size_t matchedCount{0};
    std::size_t skippedBelowRangeCount{0};
    std::size_t skippedAboveRangeCount{0};
    double totalParticleWeight{0.0};
    double totalFrequencyWeight{0.0};
    double yFWeightedNumerator{0.0};
    double yDWeightedNumerator{0.0};
    double yStarWeightedNumerator{0.0};
    double zFWeightedNumerator{0.0};
    double zDWeightedNumerator{0.0};
    double zStarWeightedNumerator{0.0};
};

struct AtSummary {
    std::string selectedParticle;
    std::size_t chargedParticleCount{0};
    std::size_t selectedParticleCount{0};
    std::size_t matchedParticleCount{0};
    std::size_t skippedParticleCount{0};
    std::size_t skippedBelowRangeCount{0};
    std::size_t skippedAboveRangeCount{0};
    double totalMatchedWeight{0.0};
    double totalFrequencyWeight{0.0};
    double yFKeVPerUm{0.0};
    double yDKeVPerUm{0.0};
    double yStarKeVPerUm{0.0};
    double zFGy{0.0};
    double zDGy{0.0};
    double zStarGy{0.0};
    double yFWeightedNumerator{0.0};
    double yDWeightedNumerator{0.0};
    double yStarWeightedNumerator{0.0};
    double zFWeightedNumerator{0.0};
    double zDWeightedNumerator{0.0};
    double zStarWeightedNumerator{0.0};
    std::vector<AtDiagnosticRecord> diagnostics;
};

class AtCalculator {
public:
    // The AT totals are evaluated from charged phase-space rows whose decoded
    // ion identity matches the selected AT LUT. Each matched row queries the
    // LUT on the MeV/u axis using KE/A. The mixed-field totals are weighted by
    // the phase-space row weight, so each row contributes according to how
    // frequently that KE appears in the selected-particle phase-space sample.
    AtSummary calculate(
        const AtLookup& lookup,
        const std::vector<ChargedParticleRecord>& chargedParticles,
        AtParticle selectedParticle) const;
};
