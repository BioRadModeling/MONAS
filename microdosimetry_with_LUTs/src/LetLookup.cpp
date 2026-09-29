#include "LetLookup.h"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

bool tryParseDouble(const std::string& s, double& value) {
    try {
        std::size_t idx = 0;
        value = std::stod(s, &idx);
        return idx > 0;
    } catch (...) {
        return false;
    }
}

}  // namespace

LetLookup::LetLookup(std::string element,
                     std::vector<double> energiesMeV,
                     std::vector<double> letKeVPerUm,
                     std::string sourceFile)
    : element_(std::move(element)),
      energiesMeV_(std::move(energiesMeV)),
      letKeVPerUm_(std::move(letKeVPerUm)),
      sourceFile_(std::move(sourceFile)) {}

LetLookup LetLookup::loadFromFile(const std::filesystem::path& tablePath,
                                  std::string element) {
    std::ifstream in(tablePath);
    if (!in) {
        throw std::runtime_error("Failed to open LET lookup table: " +
                                 tablePath.string());
    }

    std::vector<double> energiesMeV;
    std::vector<double> letKeVPerUm;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::string energyToken;
        std::string letToken;
        if (!(iss >> energyToken >> letToken)) {
            continue;
        }

        double energyMeV = 0.0;
        double letValue = 0.0;
        if (!tryParseDouble(energyToken, energyMeV) ||
            !tryParseDouble(letToken, letValue)) {
            continue;
        }

        energiesMeV.push_back(energyMeV);
        letKeVPerUm.push_back(letValue);
    }

    if (energiesMeV.empty() || letKeVPerUm.empty()) {
        throw std::runtime_error("No usable LET rows found in: " +
                                 tablePath.string());
    }

    if (energiesMeV.size() != letKeVPerUm.size()) {
        throw std::runtime_error("Mismatched LET table column sizes in: " +
                                 tablePath.string());
    }

    for (std::size_t i = 1; i < energiesMeV.size(); ++i) {
        if (energiesMeV[i] <= energiesMeV[i - 1]) {
            throw std::runtime_error(
                "LET table energies must be strictly increasing in: " +
                tablePath.string());
        }
    }

    return LetLookup(std::move(element),
                     std::move(energiesMeV),
                     std::move(letKeVPerUm),
                     tablePath.string());
}

double LetLookup::interpolate(double energyMeV) const {
    if (energiesMeV_.empty()) {
        throw std::runtime_error("Cannot interpolate an empty LET lookup table.");
    }

    if (energyMeV <= energiesMeV_.front()) {
        return letKeVPerUm_.front();
    }

    if (energyMeV >= energiesMeV_.back()) {
        return letKeVPerUm_.back();
    }

    const auto upperIt =
        std::upper_bound(energiesMeV_.begin(), energiesMeV_.end(), energyMeV);
    const std::size_t upperIndex =
        static_cast<std::size_t>(upperIt - energiesMeV_.begin());
    const std::size_t lowerIndex = upperIndex - 1;

    const double lowerEnergy = energiesMeV_[lowerIndex];
    const double upperEnergy = energiesMeV_[upperIndex];
    const double fraction =
        (energyMeV - lowerEnergy) / (upperEnergy - lowerEnergy);

    return letKeVPerUm_[lowerIndex] * (1.0 - fraction) +
           letKeVPerUm_[upperIndex] * fraction;
}

const std::string& LetLookup::element() const {
    return element_;
}

const std::vector<double>& LetLookup::energiesMeV() const {
    return energiesMeV_;
}

const std::vector<double>& LetLookup::letKeVPerUm() const {
    return letKeVPerUm_;
}

const std::string& LetLookup::sourceFile() const {
    return sourceFile_;
}
