#include "LookupLibrary.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

LutFamily LookupLibrary::inferFamily(const std::string& lutName) {
    if (lutName == "Cartechini") {
        return LutFamily::Cartechini;
    }

    return LutFamily::DeCunha;
}

bool LookupLibrary::isLookupCsvFile(const fs::path& filePath) {
    if (!fs::is_regular_file(filePath)) {
        return false;
    }

    if (filePath.extension() != ".csv") {
        return false;
    }

    const std::string name = filePath.filename().string();

    return name.rfind("Proton_", 0) == 0 &&
           name.size() > std::string("Proton__MeV.csv").size() &&
           name.find("_MeV.csv") != std::string::npos;
}

double LookupLibrary::parseEnergyFromFilename(const fs::path& filePath) {
    const std::string name = filePath.filename().string();
    const std::string prefix = "Proton_";
    const std::string suffix = "_MeV.csv";

    if (name.rfind(prefix, 0) != 0) {
        throw std::runtime_error("Invalid lookup filename: " + name);
    }

    if (name.size() <= prefix.size() + suffix.size()) {
        throw std::runtime_error("Invalid lookup filename: " + name);
    }

    if (name.substr(name.size() - suffix.size()) != suffix) {
        throw std::runtime_error("Invalid lookup filename: " + name);
    }

    const std::string energyStr =
        name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());

    try {
        return std::stod(energyStr);
    } catch (...) {
        throw std::runtime_error("Failed to parse energy from filename: " + name);
    }
}

bool LookupLibrary::isCartechiniEnergyFolder(const fs::path& folderPath) {
    if (!fs::is_directory(folderPath)) {
        return false;
    }

    const std::string name = folderPath.filename().string();
    const std::string prefix = "H_E";
    const std::string suffix = "_R0.5";

    return name.rfind(prefix, 0) == 0 &&
           name.size() > prefix.size() + suffix.size() &&
           name.substr(name.size() - suffix.size()) == suffix;
}

double LookupLibrary::parseEnergyFromCartechiniFolderName(
    const fs::path& folderPath) {
    const std::string name = folderPath.filename().string();
    const std::string prefix = "H_E";
    const std::string suffix = "_R0.5";

    if (name.rfind(prefix, 0) != 0 ||
        name.size() <= prefix.size() + suffix.size() ||
        name.substr(name.size() - suffix.size()) != suffix) {
        throw std::runtime_error(
            "Invalid Cartechini folder name: " + name);
    }

    const std::string energyStr =
        name.substr(prefix.size(),
                    name.size() - prefix.size() - suffix.size());

    try {
        return std::stod(energyStr);
    } catch (...) {
        throw std::runtime_error(
            "Failed to parse energy from Cartechini folder name: " + name);
    }
}

void LookupLibrary::loadDeCunhaDirectory(const fs::path& libraryDir) {
    for (const auto& entry : fs::directory_iterator(libraryDir)) {
        const fs::path filePath = entry.path();

        if (!isLookupCsvFile(filePath)) {
            continue;
        }

        const double energyMeV = parseEnergyFromFilename(filePath);
        tables_.push_back(LookupTable::loadFromCsv(filePath, energyMeV));
    }
}

void LookupLibrary::loadCartechiniDirectory(const fs::path& libraryDir) {
    for (const auto& entry : fs::directory_iterator(libraryDir)) {
        const fs::path folderPath = entry.path();

        if (!isCartechiniEnergyFolder(folderPath)) {
            continue;
        }

        const double energyMeV =
            parseEnergyFromCartechiniFolderName(folderPath);

        const fs::path ySpecPath = folderPath / "ySpecfile.txt";
        if (!fs::exists(ySpecPath) || !fs::is_regular_file(ySpecPath)) {
            throw std::runtime_error(
                "Missing ySpecfile.txt in Cartechini folder: " +
                folderPath.string());
        }

        tables_.push_back(
            LookupTable::loadFromCartechiniYSpec(ySpecPath, energyMeV));
    }
}

void LookupLibrary::loadFromDirectory(const fs::path& libraryDir,
                                      LutFamily family) {
    tables_.clear();
    yReference_.clear();

    if (!fs::exists(libraryDir)) {
        throw std::runtime_error(
            "Lookup directory does not exist: " + libraryDir.string());
    }

    if (!fs::is_directory(libraryDir)) {
        throw std::runtime_error(
            "Lookup path is not a directory: " + libraryDir.string());
    }

    if (family == LutFamily::Cartechini) {
        loadCartechiniDirectory(libraryDir);
    } else {
        loadDeCunhaDirectory(libraryDir);
    }

    if (tables_.empty()) {
        throw std::runtime_error(
            "No lookup tables found in directory: " + libraryDir.string());
    }

    std::sort(tables_.begin(), tables_.end(),
              [](const LookupTable& a, const LookupTable& b) {
                  return a.monoEnergyMeV() < b.monoEnergyMeV();
              });

    yReference_ = tables_.front().yLowerEdges();
    validateConsistentYGrid();
}

std::size_t LookupLibrary::findNearestIndex(double energyMeV) const {
    if (tables_.empty()) {
        throw std::runtime_error(
            "LookupLibrary::findNearestIndex called on empty library.");
    }

    auto it = std::lower_bound(
        tables_.begin(), tables_.end(), energyMeV,
        [](const LookupTable& table, double value) {
            return table.monoEnergyMeV() < value;
        });

    if (it == tables_.begin()) {
        return 0;
    }

    if (it == tables_.end()) {
        return tables_.size() - 1;
    }

    const std::size_t upperIdx =
        static_cast<std::size_t>(std::distance(tables_.begin(), it));
    const std::size_t lowerIdx = upperIdx - 1;

    const double dLower = std::fabs(tables_[lowerIdx].monoEnergyMeV() - energyMeV);
    const double dUpper = std::fabs(tables_[upperIdx].monoEnergyMeV() - energyMeV);

    return (dLower <= dUpper) ? lowerIdx : upperIdx;
}

const LookupTable& LookupLibrary::findNearest(double energyMeV) const {
    return tables_[findNearestIndex(energyMeV)];
}

std::size_t LookupLibrary::size() const {
    return tables_.size();
}

const std::vector<double>& LookupLibrary::yReference() const {
    return yReference_;
}

const std::vector<LookupTable>& LookupLibrary::tables() const {
    return tables_;
}

void LookupLibrary::validateConsistentYGrid() const {
    if (tables_.empty()) {
        return;
    }

    const auto& ref = tables_.front().yLowerEdges();

    for (std::size_t i = 1; i < tables_.size(); ++i) {
        const auto& y = tables_[i].yLowerEdges();

        if (y.size() != ref.size()) {
            throw std::runtime_error(
                "Mismatch in y-grid size across lookup tables.");
        }

        for (std::size_t j = 0; j < ref.size(); ++j) {
            if (std::fabs(y[j] - ref[j]) > 1e-12) {
                throw std::runtime_error(
                    "Mismatch in y-grid values across lookup tables.");
            }
        }
    }
}