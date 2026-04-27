#include "InaniwaCalculator.h"

#include <stdexcept>

InaniwaSummary InaniwaCalculator::calculate(
    const InaniwaLookup& lookup,
    const std::vector<ChargedParticleRecord>& chargedParticles) const {

    InaniwaSummary summary;
    double zdDWeightedSum = 0.0;
    double zdDStarWeightedSum = 0.0;
    double znDWeightedSum = 0.0;

    for (const auto& particle : chargedParticles) {
        if (!hasInaniwaLutIdentity(particle)) {
            continue;
        }

        if (!lookup.hasAtomicNumber(particle.atomicNumber)) {
            continue;
        }

        const double eventEnergyMeV = inaniwaEventEnergyMeV(particle);
        if (eventEnergyMeV < 0.0) {
            throw std::runtime_error(
                "Inaniwa calculation encountered a negative event energy.");
        }

        if (eventEnergyMeV == 0.0) {
            continue;
        }

        const double energyMeVPerU = inaniwaEnergyMeVPerU(particle);
        const InaniwaInterpolatedValues lut =
            lookup.interpolate(particle.atomicNumber, energyMeVPerU);

        summary.matchedParticleCount += 1;
        summary.totalEventEnergyMeV += eventEnergyMeV;
        zdDWeightedSum += eventEnergyMeV * lut.zdDMeanGy;
        zdDStarWeightedSum += eventEnergyMeV * lut.zdDStarMeanGy;
        znDWeightedSum += eventEnergyMeV * lut.znDMeanGy;
    }

    if (summary.matchedParticleCount == 0) {
        throw std::runtime_error(
            "Inaniwa calculation requires at least one supported charged ion "
            "with positive kinetic energy in the phase-space input.");
    }

    if (summary.totalEventEnergyMeV <= 0.0) {
        throw std::runtime_error(
            "Inaniwa calculation produced a non-positive total event energy.");
    }

    summary.zdDMeanGy = zdDWeightedSum / summary.totalEventEnergyMeV;
    summary.zdDStarMeanGy = zdDStarWeightedSum / summary.totalEventEnergyMeV;
    summary.znDMeanGy = znDWeightedSum / summary.totalEventEnergyMeV;

    return summary;
}
