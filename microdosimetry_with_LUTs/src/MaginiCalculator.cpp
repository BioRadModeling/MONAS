#include "MaginiCalculator.h"

#include <cmath>
#include <stdexcept>

MaginiSummary MaginiCalculator::calculate(
    const MaginiLookup& lookup,
    const std::vector<ProtonRecord>& protons) const {

    MaginiSummary summary;
    summary.protonCount = protons.size();

    for (const auto& proton : protons) {
        if (proton.weight < 0.0) {
            throw std::runtime_error(
                "Magini calculation encountered a negative proton weight.");
        }

        summary.totalProtonWeight += proton.weight;
    }

    if (summary.protonCount == 0) {
        throw std::runtime_error(
            "Magini calculation requires at least one proton in the phase-space input.");
    }

    if (summary.totalProtonWeight <= 0.0) {
        throw std::runtime_error(
            "Magini calculation requires a positive total proton weight.");
    }

    // Interpret the proton phase-space rows as a discrete weighted sample of
    // the energy distribution. This makes the discrete sums in the Magini
    // formulas equivalent to weighted averages over proton records after LUT
    // interpolation, without introducing an arbitrary energy histogram.
    for (const auto& proton : protons) {
        const double normalizedWeight = proton.weight / summary.totalProtonWeight;
        const MaginiInterpolatedValues lut = lookup.interpolate(proton.energyMeV);

        summary.yFKeVPerUm += lut.yFLutKeVPerUm * normalizedWeight;
        summary.yFMaxErrorKeVPerUm +=
            lut.yFMaxErrorKeVPerUm * normalizedWeight;
        summary.yDWeightedNumerator +=
            lut.yDLutKeVPerUm * lut.yFLutKeVPerUm * normalizedWeight;
        summary.yStarWeightedNumerator +=
            lut.yStarLutKeVPerUm * lut.yFLutKeVPerUm * normalizedWeight;
    }

    if (summary.yFKeVPerUm <= 0.0) {
        throw std::runtime_error(
            "Magini calculation produced a non-positive y_F value.");
    }

    summary.yDKeVPerUm = summary.yDWeightedNumerator / summary.yFKeVPerUm;
    summary.yStarKeVPerUm = summary.yStarWeightedNumerator / summary.yFKeVPerUm;

    for (const auto& proton : protons) {
        const double normalizedWeight = proton.weight / summary.totalProtonWeight;
        const MaginiInterpolatedValues lut = lookup.interpolate(proton.energyMeV);

        summary.yDMaxErrorKeVPerUm +=
            std::abs(normalizedWeight * lut.yFLutKeVPerUm /
                     summary.yFKeVPerUm) *
                lut.yDMaxErrorKeVPerUm +
            std::abs(normalizedWeight *
                     (lut.yDLutKeVPerUm - summary.yDKeVPerUm) /
                     summary.yFKeVPerUm) *
                lut.yFMaxErrorKeVPerUm;
    }

    return summary;
}
