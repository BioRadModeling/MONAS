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
    double lowerMatchedEnergyMeV;
    double upperMatchedEnergyMeV;
    double lowerInterpolationWeight;
    double upperInterpolationWeight;
    std::string lowerMatchedFile;
    std::string upperMatchedFile;
    double lowerMatchedNcpp;
    double upperMatchedNcpp;
    std::string lowerMatchedFamily;
    std::string upperMatchedFamily;
};

class CsvWriter {
public:
    static void writeProtonMatches(const std::filesystem::path& outPath,
                                   const std::vector<ProtonMatchRecord>& matches);

    static void writePolySpectrum(const std::filesystem::path& outPath,
                                  const PolySpectrum& spectrum);
};
