#pragma once

#include <cstddef>
#include <vector>

#include "InaniwaLookup.h"
#include "PhaseSpaceReader.h"

struct InaniwaSummary {
    std::size_t matchedParticleCount{0};
    double totalEventEnergyMeV{0.0};
    double zdDMeanGy{0.0};
    double zdDStarMeanGy{0.0};
    double znDMeanGy{0.0};
};

class InaniwaCalculator {
public:
    // The Inaniwa totals are evaluated over all charged phase-space rows that
    // decode to supported ion families with Z=1..10. Each row selects its LUT
    // table from the decoded atomic number, uses e_k := KE from the phase
    // space, and queries the LUT on the MeV/u axis using KE/A.
    InaniwaSummary calculate(
        const InaniwaLookup& lookup,
        const std::vector<ChargedParticleRecord>& chargedParticles) const;
};
