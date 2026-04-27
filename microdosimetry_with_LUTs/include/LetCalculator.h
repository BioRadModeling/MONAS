#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "LetAccumulator.h"
#include "PhaseSpaceReader.h"

struct LetCalculationResult {
    std::vector<LetSummaryRecord> records;
    std::size_t matchedParticleCount{0};
};

class LetCalculator {
public:
    LetCalculationResult calculate(
        const std::filesystem::path& letDirectory,
        const std::vector<ChargedParticleRecord>& chargedParticles) const;

private:
    static std::optional<std::string> lutElementForPdg(int pdgCode);
    static void addParticleToAccumulator(LetAccumulator& accumulator,
                                         const ChargedParticleRecord& particle,
                                         double letKeVPerUm);
};
