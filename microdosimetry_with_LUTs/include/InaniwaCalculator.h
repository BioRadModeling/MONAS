#pragma once

#include <cstddef>
#include <vector>

#include "InaniwaLookup.h"
#include "PhaseSpaceReader.h"

struct InaniwaSummaryRecord {
    int atomicNumber{0};
    std::size_t particleCount{0};
    double totalEventEnergyMeV{0.0};
    double zdDMeanGy{0.0};
    double zdDStarMeanGy{0.0};
    double znDMeanGy{0.0};
};

class InaniwaCalculator {
public:
    // The Inaniwa totals are evaluated independently for each LUT atomic
    // number. The same charged phase-space rows are reused for every Z=1..10
    // table, with e_k := KE from the phase space and the LUT query energy
    // taken directly from the row kinetic energy.
    std::vector<InaniwaSummaryRecord> calculate(
        const InaniwaLookup& lookup,
        const std::vector<ChargedParticleRecord>& chargedParticles) const;
};
