#include "LookupTable.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

std::vector<std::string> splitCsvLine(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string item;

    while (std::getline(ss, item, ',')) {
        fields.push_back(item);
    }

    return fields;
}

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

LookupTable::LookupTable(double monoEnergyMeV,
                         LookupSpectrumKind spectrumKind,
                         std::vector<double> yLowerEdges,
                         std::vector<double> nY,
                         std::vector<double> fY,
                         double ncpp,
                         std::string sourceFile)
    : monoEnergyMeV_(monoEnergyMeV),
      spectrumKind_(spectrumKind),
      yLowerEdges_(std::move(yLowerEdges)),
      nY_(std::move(nY)),
      fY_(std::move(fY)),
      ncpp_(ncpp),
      sourceFile_(std::move(sourceFile)) {}

LookupTable LookupTable::loadFromCsv(const std::filesystem::path& csvPath,
                                     double monoEnergyMeV) {
    std::ifstream in(csvPath);
    if (!in) {
        throw std::runtime_error("Failed to open lookup CSV: " + csvPath.string());
    }

    std::vector<double> yLowerEdges;
    std::vector<double> nY;
    double detectedNcpp = -1.0;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        auto fields = splitCsvLine(line);
        if (fields.size() < 2) {
            continue;
        }

        double y = 0.0;
        double ny = 0.0;
        if (!tryParseDouble(fields[0], y) || !tryParseDouble(fields[1], ny)) {
            continue;
        }

        yLowerEdges.push_back(y);
        nY.push_back(ny);

        if (fields.size() >= 3) {
            double maybeNcpp = 0.0;
            if (tryParseDouble(fields[2], maybeNcpp) && maybeNcpp > 0.0) {
                if (detectedNcpp < 0.0) {
                    detectedNcpp = maybeNcpp;
                } else if (std::fabs(detectedNcpp - maybeNcpp) > 1e-12) {
                    throw std::runtime_error(
                        "Inconsistent Ncpp values found in: " + csvPath.string());
                }
            }
        }
    }

    if (yLowerEdges.empty() || nY.empty()) {
        throw std::runtime_error("No usable spectrum rows found in: " + csvPath.string());
    }

    if (detectedNcpp <= 0.0) {
        throw std::runtime_error("Failed to detect valid Ncpp in: " + csvPath.string());
    }

    return LookupTable(monoEnergyMeV,
                       LookupSpectrumKind::DeCunhaRawCounts,
                       std::move(yLowerEdges),
                       std::move(nY),
                       {},
                       detectedNcpp,
                       csvPath.string());
}

LookupTable LookupTable::loadFromCartechiniYSpec(
    const std::filesystem::path& ySpecPath,
    double monoEnergyMeV) {
    std::ifstream in(ySpecPath);
    if (!in) {
        throw std::runtime_error("Failed to open Cartechini ySpec file: " +
                                 ySpecPath.string());
    }

    std::vector<double> yLowerEdges;
    std::vector<double> fY;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        double y = 0.0;
        double fy = 0.0;

        // Numeric data rows begin with y then f(y).
        // Summary/header rows fail this parse and are skipped.
        if (!(iss >> y >> fy)) {
            continue;
        }

        yLowerEdges.push_back(y);
        fY.push_back(fy);
    }

    if (yLowerEdges.empty() || fY.empty()) {
        throw std::runtime_error(
            "No valid y / f(y) rows found in file: " + ySpecPath.string());
    }

    return LookupTable(monoEnergyMeV,
                       LookupSpectrumKind::CartechiniFy,
                       std::move(yLowerEdges),
                       {},
                       std::move(fY),
                       std::numeric_limits<double>::quiet_NaN(),
                       ySpecPath.string());
}

double LookupTable::monoEnergyMeV() const {
    return monoEnergyMeV_;
}

LookupSpectrumKind LookupTable::spectrumKind() const {
    return spectrumKind_;
}

double LookupTable::ncpp() const {
    return ncpp_;
}

const std::vector<double>& LookupTable::yLowerEdges() const {
    return yLowerEdges_;
}

const std::vector<double>& LookupTable::nY() const {
    return nY_;
}

const std::vector<double>& LookupTable::fY() const {
    return fY_;
}

const std::string& LookupTable::sourceFile() const {
    return sourceFile_;
}