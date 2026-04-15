#pragma once

#include <filesystem>
#include <vector>

#include "LookupTable.h"

class LookupLibrary {
public:
    void loadFromDirectory(const std::filesystem::path& libraryDir);

    const LookupTable& findNearest(double energyMeV) const;

    std::size_t size() const;
    const std::vector<double>& yReference() const;

private:
    std::vector<LookupTable> tables_;
    std::vector<double> yReference_;

    static bool isLookupCsvFile(const std::filesystem::path& filePath);
    static double parseEnergyFromFilename(const std::filesystem::path& filePath);

    void validateConsistentYGrid() const;
};