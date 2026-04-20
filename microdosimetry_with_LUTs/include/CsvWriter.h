#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include "SpectrumAccumulator.h"

struct ProtonMatchRecord {
    std::size_t rowIndex;
    double energyMeV;
    double weight;
    double matchedEnergyMeV;
    std::string matchedFile;
    double matchedNcpp;
    std::string matchedFamily;
};

class CsvWriter {
public:
    static void writeProtonMatches(const std::filesystem::path& outPath,
                                   const std::vector<ProtonMatchRecord>& matches);

    static void writePolySpectrum(const std::filesystem::path& outPath,
                                  const PolySpectrum& spectrum);
};