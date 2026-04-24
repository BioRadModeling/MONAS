#include "InaniwaCalculator.h"

#include <array>
#include <stdexcept>

namespace {

struct InaniwaAccumulator {
    std::size_t particleCount{0};
    double totalEventEnergyMeV{0.0};
    double zdDWeightedSum{0.0};
    double zdDStarWeightedSum{0.0};
    double znDWeightedSum{0.0};
};

}  // namespace

std::vector<InaniwaSummaryRecord> InaniwaCalculator::calculate(
    const InaniwaLookup& lookup,
    const std::vector<ChargedParticleRecord>& chargedParticles) const {

    std::array<InaniwaAccumulator, 11> accumulators;
    std::size_t contributingParticleCount = 0;

    for (const auto& particle : chargedParticles) {
        const double eventEnergyMeV = inaniwaEventEnergyMeV(particle);
        if (eventEnergyMeV < 0.0) {
            throw std::runtime_error(
                "Inaniwa calculation encountered a negative event energy.");
        }

        if (eventEnergyMeV == 0.0) {
            continue;
        }

        ++contributingParticleCount;

        for (int atomicNumber = 1; atomicNumber <= 10; ++atomicNumber) {
            const InaniwaInterpolatedValues lut =
                lookup.interpolate(atomicNumber, eventEnergyMeV);

            InaniwaAccumulator& accumulator = accumulators[atomicNumber];
            accumulator.particleCount += 1;
            accumulator.totalEventEnergyMeV += eventEnergyMeV;
            accumulator.zdDWeightedSum += eventEnergyMeV * lut.zdDMeanGy;
            accumulator.zdDStarWeightedSum += eventEnergyMeV * lut.zdDStarMeanGy;
            accumulator.znDWeightedSum += eventEnergyMeV * lut.znDMeanGy;
        }
    }

    if (contributingParticleCount == 0) {
        throw std::runtime_error(
            "Inaniwa calculation requires at least one charged particle "
            "with positive kinetic energy in the phase-space input.");
    }

    std::vector<InaniwaSummaryRecord> summaries;
    summaries.reserve(10);

    for (int atomicNumber = 1; atomicNumber <= 10; ++atomicNumber) {
        const InaniwaAccumulator& accumulator = accumulators[atomicNumber];
        if (accumulator.particleCount == 0) {
            continue;
        }

        if (accumulator.totalEventEnergyMeV <= 0.0) {
            throw std::runtime_error(
                "Inaniwa calculation produced a non-positive total event energy.");
        }

        summaries.push_back(InaniwaSummaryRecord{
            atomicNumber,
            accumulator.particleCount,
            accumulator.totalEventEnergyMeV,
            accumulator.zdDWeightedSum / accumulator.totalEventEnergyMeV,
            accumulator.zdDStarWeightedSum / accumulator.totalEventEnergyMeV,
            accumulator.znDWeightedSum / accumulator.totalEventEnergyMeV
        });
    }

    return summaries;
}
