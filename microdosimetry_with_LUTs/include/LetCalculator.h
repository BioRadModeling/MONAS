#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "LetAccumulator.h"
#include "PhaseSpaceReader.h"

struct ElementLetSummary {
    std::string element;
    std::vector<LetSummaryRecord> records;
};

class LetCalculator {
public:
    static const std::vector<std::string>& supportedElements();

    std::vector<ElementLetSummary> calculateByElement(
        const std::filesystem::path& letDirectory,
        const std::vector<ChargedParticleRecord>& chargedParticles) const;

private:
    static void addParticleToAccumulator(LetAccumulator& accumulator,
                                         const ChargedParticleRecord& particle,
                                         double letKeVPerUm);
};
