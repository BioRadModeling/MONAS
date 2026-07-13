#include "AtCalculator.h"

#include <stdexcept>

namespace {

std::string selectedParticleName(AtParticle selectedParticle) {
    switch (selectedParticle) {
        case AtParticle::Proton:
            return "proton";
        case AtParticle::Carbon:
            return "carbon";
    }

    return "unknown";
}

bool matchesSelectedParticle(const ChargedParticleRecord& particle,
                             AtParticle selectedParticle) {
    switch (selectedParticle) {
        case AtParticle::Proton:
            return particle.atomicNumber == 1 && particle.massNumber == 1;
        case AtParticle::Carbon:
            return particle.atomicNumber == 6 && particle.massNumber == 12;
    }

    return false;
}

std::vector<AtDiagnosticRecord> makeDiagnosticRecords() {
    return std::vector<AtDiagnosticRecord>{
        AtDiagnosticRecord{"<10"},
        AtDiagnosticRecord{"10-20"},
        AtDiagnosticRecord{"20-50"},
        AtDiagnosticRecord{"50-100"},
        AtDiagnosticRecord{"100-250"},
        AtDiagnosticRecord{">250"}
    };
}

AtDiagnosticRecord& diagnosticForEnergy(std::vector<AtDiagnosticRecord>& records,
                                        double energyMeVPerU) {
    if (energyMeVPerU < 10.0) {
        return records[0];
    }
    if (energyMeVPerU < 20.0) {
        return records[1];
    }
    if (energyMeVPerU < 50.0) {
        return records[2];
    }
    if (energyMeVPerU < 100.0) {
        return records[3];
    }
    if (energyMeVPerU <= 250.0) {
        return records[4];
    }
    return records[5];
}

}  // namespace

AtSummary AtCalculator::calculate(
    const AtLookup& lookup,
    const std::vector<ChargedParticleRecord>& chargedParticles,
    AtParticle selectedParticle) const {

    AtSummary summary;
    summary.selectedParticle = selectedParticleName(selectedParticle);
    summary.chargedParticleCount = chargedParticles.size();
    summary.diagnostics = makeDiagnosticRecords();

    for (const auto& particle : chargedParticles) {
        if (!hasInaniwaLutIdentity(particle) ||
            !matchesSelectedParticle(particle, selectedParticle)) {
            summary.skippedParticleCount += 1;
            continue;
        }

        summary.selectedParticleCount += 1;

        if (particle.weight < 0.0) {
            throw std::runtime_error(
                "AT calculation encountered a negative particle weight.");
        }

        if (particle.energyMeV < 0.0) {
            throw std::runtime_error(
                "AT calculation encountered a negative particle energy.");
        }

        if (particle.weight == 0.0) {
            summary.skippedParticleCount += 1;
            continue;
        }

        if (!lookup.hasIon(particle.atomicNumber, particle.massNumber)) {
            summary.skippedParticleCount += 1;
            continue;
        }

        const double energyMeVPerU = inaniwaEnergyMeVPerU(particle);
        AtDiagnosticRecord& diagnostic =
            diagnosticForEnergy(summary.diagnostics, energyMeVPerU);
        diagnostic.rowCount += 1;
        diagnostic.totalParticleWeight += particle.weight;

        const double minEnergyMeVPerU =
            lookup.minEnergyMeVPerU(particle.atomicNumber, particle.massNumber);
        const double maxEnergyMeVPerU =
            lookup.maxEnergyMeVPerU(particle.atomicNumber, particle.massNumber);

        if (energyMeVPerU < minEnergyMeVPerU) {
            summary.skippedParticleCount += 1;
            summary.skippedBelowRangeCount += 1;
            diagnostic.skippedBelowRangeCount += 1;
            continue;
        }

        if (energyMeVPerU > maxEnergyMeVPerU) {
            summary.skippedParticleCount += 1;
            summary.skippedAboveRangeCount += 1;
            diagnostic.skippedAboveRangeCount += 1;
            continue;
        }

        const AtInterpolatedValues lut =
            lookup.interpolate(particle.atomicNumber,
                               particle.massNumber,
                               energyMeVPerU);

        summary.matchedParticleCount += 1;
        summary.totalMatchedWeight += particle.weight;
        diagnostic.matchedCount += 1;

        const double frequencyWeight = particle.weight;
        if (frequencyWeight <= 0.0) {
            summary.skippedParticleCount += 1;
            continue;
        }

        summary.totalFrequencyWeight += frequencyWeight;
        diagnostic.totalFrequencyWeight += frequencyWeight;

        summary.yFWeightedNumerator += frequencyWeight * lut.yFKeVPerUm;
        summary.yDWeightedNumerator += frequencyWeight * lut.yDKeVPerUm;
        summary.yStarWeightedNumerator += frequencyWeight * lut.yStarKeVPerUm;
        summary.zFWeightedNumerator += frequencyWeight * lut.zFGy;
        summary.zDWeightedNumerator += frequencyWeight * lut.zDGy;
        summary.zStarWeightedNumerator += frequencyWeight * lut.zStarGy;

        diagnostic.yFWeightedNumerator += frequencyWeight * lut.yFKeVPerUm;
        diagnostic.yDWeightedNumerator += frequencyWeight * lut.yDKeVPerUm;
        diagnostic.yStarWeightedNumerator += frequencyWeight * lut.yStarKeVPerUm;
        diagnostic.zFWeightedNumerator += frequencyWeight * lut.zFGy;
        diagnostic.zDWeightedNumerator += frequencyWeight * lut.zDGy;
        diagnostic.zStarWeightedNumerator += frequencyWeight * lut.zStarGy;
    }

    if (summary.matchedParticleCount == 0) {
        throw std::runtime_error(
            "AT calculation requires at least one charged particle matching "
            "the selected AT lookup table.");
    }

    if (summary.totalMatchedWeight <= 0.0) {
        throw std::runtime_error(
            "AT calculation produced a non-positive total matched weight.");
    }

    if (summary.totalFrequencyWeight <= 0.0) {
        throw std::runtime_error(
            "AT calculation produced a non-positive frequency-weighted total.");
    }

    summary.yFKeVPerUm = summary.yFWeightedNumerator / summary.totalFrequencyWeight;
    summary.yDKeVPerUm = summary.yDWeightedNumerator / summary.totalFrequencyWeight;
    summary.yStarKeVPerUm = summary.yStarWeightedNumerator / summary.totalFrequencyWeight;
    summary.zFGy = summary.zFWeightedNumerator / summary.totalFrequencyWeight;
    summary.zDGy = summary.zDWeightedNumerator / summary.totalFrequencyWeight;
    summary.zStarGy = summary.zStarWeightedNumerator / summary.totalFrequencyWeight;

    return summary;
}
