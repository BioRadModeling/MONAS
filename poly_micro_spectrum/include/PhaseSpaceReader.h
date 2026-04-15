#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

struct ProtonRecord {
    std::size_t rowIndex;
    double energyMeV;
    double weight;
    int pdgCode;
};

class PhaseSpaceReader {
public:
    std::vector<ProtonRecord> readProtons(const std::filesystem::path& phspPath) const;
};