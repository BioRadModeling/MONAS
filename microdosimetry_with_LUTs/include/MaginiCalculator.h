#pragma once

#include <cstddef>
#include <vector>

#include "MaginiLookup.h"
#include "PhaseSpaceReader.h"

struct MaginiSummary {
    std::size_t protonCount{0};
    double totalProtonWeight{0.0};
    double yFKeVPerUm{0.0};
    double yDKeVPerUm{0.0};
    double yStarKeVPerUm{0.0};
    double yDWeightedNumerator{0.0};
    double yStarWeightedNumerator{0.0};
};

class MaginiCalculator {
public:
    // The Magini totals are evaluated directly from weighted proton samples in
    // the phase-space input. Each proton contributes its interpolated LUT
    // values with its phase-space weight, which avoids introducing an
    // additional histogram/binning step on energy.
    MaginiSummary calculate(const MaginiLookup& lookup,
                            const std::vector<ProtonRecord>& protons) const;
};
