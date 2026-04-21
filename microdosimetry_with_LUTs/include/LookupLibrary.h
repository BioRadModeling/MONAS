#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include "LookupTable.h"

enum class LutFamily {
    DeCunha,
    Cartechini
};

class LookupLibrary {
public:
    static LutFamily inferFamily(const std::string& lutName);

    void loadFromDirectory(const std::filesystem::path& libraryDir,
                           LutFamily family);

    const LookupTable& findNearest(double energyMeV) const;
    std::size_t findNearestIndex(double energyMeV) const;

    std::size_t size() const;
    const std::vector<double>& yReference() const;
    const std::vector<LookupTable>& tables() const;

private:
    std::vector<LookupTable> tables_;
    std::vector<double> yReference_;

    LutFamily activeFamily_{LutFamily::DeCunha};

    std::vector<std::size_t> cartechiniIndices_;
    std::vector<std::size_t> fallbackDeCunhaIndices_;
    double maxCartechiniEnergyMeV_{-1.0};

    void loadDeCunhaDirectory(const std::filesystem::path& libraryDir);
    void loadCartechiniDirectory(const std::filesystem::path& libraryDir);

    static bool isLookupCsvFile(const std::filesystem::path& filePath);
    static double parseEnergyFromFilename(const std::filesystem::path& filePath);

    static bool isCartechiniEnergyFolder(const std::filesystem::path& folderPath);
    static double parseEnergyFromCartechiniFolderName(
        const std::filesystem::path& folderPath);

    void validateConsistentYGrid() const;
    void rebuildFamilyIndexCaches();

    std::size_t findNearestIndexInSubset(
        double energyMeV,
        const std::vector<std::size_t>& subset) const;
};