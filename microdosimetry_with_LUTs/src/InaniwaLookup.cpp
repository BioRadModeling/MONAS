#include "InaniwaLookup.h"

#include <algorithm>
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

double interpolateLinear(double lowerValue,
                         double upperValue,
                         double lowerWeight,
                         double upperWeight) {
    return lowerValue * lowerWeight + upperValue * upperWeight;
}

InaniwaTable loadSingleTable(const std::filesystem::path& csvPath,
                             int atomicNumber) {
    std::ifstream in(csvPath);
    if (!in) {
        throw std::runtime_error("Failed to open Inaniwa CSV: " + csvPath.string());
    }

    std::vector<InaniwaRow> rows;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        const auto fields = splitCsvLine(line);
        if (fields.size() < 4) {
            continue;
        }

        InaniwaRow row;
        if (!tryParseDouble(fields[0], row.energyMeVPerU) ||
            !tryParseDouble(fields[1], row.zdDMeanGy) ||
            !tryParseDouble(fields[2], row.zdDStarMeanGy) ||
            !tryParseDouble(fields[3], row.znDMeanGy)) {
            continue;
        }

        rows.push_back(row);
    }

    if (rows.empty()) {
        throw std::runtime_error(
            "No usable Inaniwa LUT rows found in: " + csvPath.string());
    }

    std::sort(rows.begin(), rows.end(),
              [](const InaniwaRow& a, const InaniwaRow& b) {
                  return a.energyMeVPerU < b.energyMeVPerU;
              });

    for (std::size_t i = 1; i < rows.size(); ++i) {
        if (rows[i - 1].energyMeVPerU == rows[i].energyMeVPerU) {
            throw std::runtime_error(
                "Duplicate Inaniwa LUT energy detected in: " + csvPath.string());
        }
    }

    return InaniwaTable{atomicNumber, std::move(rows), csvPath.string()};
}

}  // namespace

InaniwaLookup::InaniwaLookup(std::vector<InaniwaTable> tables)
    : tables_(std::move(tables)) {}

InaniwaLookup InaniwaLookup::loadFromDirectory(
    const std::filesystem::path& directory) {

    std::vector<InaniwaTable> tables;
    tables.reserve(10);

    for (int atomicNumber = 1; atomicNumber <= 10; ++atomicNumber) {
        const std::filesystem::path csvPath =
            directory / ("Zp_" + std::to_string(atomicNumber) + ".csv");

        tables.push_back(loadSingleTable(csvPath, atomicNumber));
    }

    return InaniwaLookup(std::move(tables));
}

bool InaniwaLookup::hasAtomicNumber(int atomicNumber) const {
    return std::any_of(
        tables_.begin(), tables_.end(),
        [atomicNumber](const InaniwaTable& table) {
            return table.atomicNumber == atomicNumber;
        });
}

InaniwaInterpolatedValues InaniwaLookup::interpolate(
    int atomicNumber,
    double energyMeVPerU) const {

    const InaniwaTable& table = tableForAtomicNumber(atomicNumber);
    const std::vector<InaniwaRow>& rows = table.rows;

    if (rows.empty()) {
        throw std::runtime_error(
            "Inaniwa lookup interpolation called on empty LUT table.");
    }

    auto it = std::lower_bound(
        rows.begin(), rows.end(), energyMeVPerU,
        [](const InaniwaRow& row, double value) {
            return row.energyMeVPerU < value;
        });

    if (it == rows.begin()) {
        const InaniwaRow& row = *it;
        return InaniwaInterpolatedValues{
            row.zdDMeanGy,
            row.zdDStarMeanGy,
            row.znDMeanGy,
            row.energyMeVPerU,
            row.energyMeVPerU,
            1.0,
            0.0
        };
    }

    if (it == rows.end()) {
        const InaniwaRow& row = rows.back();
        return InaniwaInterpolatedValues{
            row.zdDMeanGy,
            row.zdDStarMeanGy,
            row.znDMeanGy,
            row.energyMeVPerU,
            row.energyMeVPerU,
            1.0,
            0.0
        };
    }

    const InaniwaRow& upper = *it;
    const InaniwaRow& lower = *(it - 1);
    const double span = upper.energyMeVPerU - lower.energyMeVPerU;

    if (span <= 0.0) {
        throw std::runtime_error("Invalid Inaniwa LUT energy spacing encountered.");
    }

    const double upperWeight = (energyMeVPerU - lower.energyMeVPerU) / span;
    const double lowerWeight = 1.0 - upperWeight;

    return InaniwaInterpolatedValues{
        interpolateLinear(lower.zdDMeanGy, upper.zdDMeanGy,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.zdDStarMeanGy, upper.zdDStarMeanGy,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.znDMeanGy, upper.znDMeanGy,
                          lowerWeight, upperWeight),
        lower.energyMeVPerU,
        upper.energyMeVPerU,
        lowerWeight,
        upperWeight
    };
}

const std::vector<InaniwaTable>& InaniwaLookup::tables() const {
    return tables_;
}

const InaniwaTable& InaniwaLookup::tableForAtomicNumber(int atomicNumber) const {
    auto it = std::find_if(
        tables_.begin(), tables_.end(),
        [atomicNumber](const InaniwaTable& table) {
            return table.atomicNumber == atomicNumber;
        });

    if (it == tables_.end()) {
        throw std::runtime_error(
            "No Inaniwa LUT table loaded for atomic number Z=" +
            std::to_string(atomicNumber));
    }

    return *it;
}
