#include "LookupTable.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

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

} // namespace

LookupTable::LookupTable(double monoEnergyMeV,
                         std::vector<double> yLowerEdges,
                         std::vector<double> nY,
                         double ncpp,
                         std::string sourceFile)
    : monoEnergyMeV_(monoEnergyMeV),
      yLowerEdges_(std::move(yLowerEdges)),
      nY_(std::move(nY)),
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
                       std::move(yLowerEdges),
                       std::move(nY),
                       detectedNcpp,
                       csvPath.string());
}

double LookupTable::monoEnergyMeV() const {
    return monoEnergyMeV_;
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

const std::string& LookupTable::sourceFile() const {
    return sourceFile_;
}