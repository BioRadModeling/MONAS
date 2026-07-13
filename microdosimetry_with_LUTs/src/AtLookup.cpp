#include "AtLookup.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <regex>
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

int atomicNumberForSymbol(const std::string& symbol) {
    if (symbol == "H") {
        return 1;
    }
    if (symbol == "He") {
        return 2;
    }
    if (symbol == "Li") {
        return 3;
    }
    if (symbol == "Be") {
        return 4;
    }
    if (symbol == "B") {
        return 5;
    }
    if (symbol == "C") {
        return 6;
    }
    if (symbol == "N") {
        return 7;
    }
    if (symbol == "O") {
        return 8;
    }
    if (symbol == "F") {
        return 9;
    }
    if (symbol == "Ne") {
        return 10;
    }

    return 0;
}

bool parseAtFilename(const std::filesystem::path& csvPath,
                     int& atomicNumber,
                     int& massNumber) {
    static const std::regex pattern(
        R"(^([0-9]+)([A-Z][a-z]?)_rd[^_]+_Rn[^_]+_LUT\.csv$)");

    std::smatch match;
    const std::string filename = csvPath.filename().string();
    if (!std::regex_match(filename, match, pattern)) {
        return false;
    }

    massNumber = std::stoi(match[1].str());
    atomicNumber = atomicNumberForSymbol(match[2].str());
    return atomicNumber > 0 && massNumber > 0;
}

AtTable loadSingleTable(const std::filesystem::path& csvPath,
                        int atomicNumber,
                        int massNumber) {
    std::ifstream in(csvPath);
    if (!in) {
        throw std::runtime_error("Failed to open AT CSV: " + csvPath.string());
    }

    std::vector<AtRow> rows;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        const auto fields = splitCsvLine(line);
        if (fields.size() < 8) {
            continue;
        }

        AtRow row;
        if (!tryParseDouble(fields[0], row.energyMeVPerU) ||
            !tryParseDouble(fields[1], row.letKeVPerUm) ||
            !tryParseDouble(fields[2], row.yFKeVPerUm) ||
            !tryParseDouble(fields[3], row.yDKeVPerUm) ||
            !tryParseDouble(fields[4], row.yStarKeVPerUm) ||
            !tryParseDouble(fields[5], row.zFGy) ||
            !tryParseDouble(fields[6], row.zDGy) ||
            !tryParseDouble(fields[7], row.zStarGy)) {
            continue;
        }

        rows.push_back(row);
    }

    if (rows.empty()) {
        throw std::runtime_error(
            "No usable AT LUT rows found in: " + csvPath.string());
    }

    std::sort(rows.begin(), rows.end(),
              [](const AtRow& a, const AtRow& b) {
                  return a.energyMeVPerU < b.energyMeVPerU;
              });

    for (std::size_t i = 1; i < rows.size(); ++i) {
        if (rows[i - 1].energyMeVPerU == rows[i].energyMeVPerU) {
            throw std::runtime_error(
                "Duplicate AT LUT energy detected in: " + csvPath.string());
        }
    }

    return AtTable{atomicNumber, massNumber, std::move(rows), csvPath.string()};
}

}  // namespace

AtLookup::AtLookup(std::vector<AtTable> tables)
    : tables_(std::move(tables)) {}

AtLookup AtLookup::loadFromDirectory(const std::filesystem::path& directory) {
    if (!std::filesystem::exists(directory) ||
        !std::filesystem::is_directory(directory)) {
        throw std::runtime_error("AT lookup folder not found: " + directory.string());
    }

    std::vector<AtTable> tables;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (!entry.is_regular_file() || entry.path().extension() != ".csv") {
            continue;
        }

        int atomicNumber = 0;
        int massNumber = 0;
        if (!parseAtFilename(entry.path(), atomicNumber, massNumber)) {
            continue;
        }

        tables.push_back(loadSingleTable(entry.path(), atomicNumber, massNumber));
    }

    if (tables.empty()) {
        throw std::runtime_error(
            "No AT lookup CSVs found in: " + directory.string());
    }

    std::sort(tables.begin(), tables.end(),
              [](const AtTable& a, const AtTable& b) {
                  if (a.atomicNumber != b.atomicNumber) {
                      return a.atomicNumber < b.atomicNumber;
                  }
                  return a.massNumber < b.massNumber;
              });

    for (std::size_t i = 1; i < tables.size(); ++i) {
        if (tables[i - 1].atomicNumber == tables[i].atomicNumber &&
            tables[i - 1].massNumber == tables[i].massNumber) {
            throw std::runtime_error(
                "Duplicate AT LUT ion identity detected in: " + directory.string());
        }
    }

    return AtLookup(std::move(tables));
}

bool AtLookup::hasIon(int atomicNumber, int massNumber) const {
    return std::any_of(
        tables_.begin(), tables_.end(),
        [atomicNumber, massNumber](const AtTable& table) {
            return table.atomicNumber == atomicNumber &&
                   table.massNumber == massNumber;
        });
}

double AtLookup::minEnergyMeVPerU(int atomicNumber, int massNumber) const {
    const AtTable& table = tableForIon(atomicNumber, massNumber);
    if (table.rows.empty()) {
        throw std::runtime_error("AT lookup min energy requested for empty LUT table.");
    }

    return table.rows.front().energyMeVPerU;
}

double AtLookup::maxEnergyMeVPerU(int atomicNumber, int massNumber) const {
    const AtTable& table = tableForIon(atomicNumber, massNumber);
    if (table.rows.empty()) {
        throw std::runtime_error("AT lookup max energy requested for empty LUT table.");
    }

    return table.rows.back().energyMeVPerU;
}

bool AtLookup::isEnergyInRange(int atomicNumber,
                               int massNumber,
                               double energyMeVPerU) const {
    const AtTable& table = tableForIon(atomicNumber, massNumber);
    if (table.rows.empty()) {
        return false;
    }

    return energyMeVPerU >= table.rows.front().energyMeVPerU &&
           energyMeVPerU <= table.rows.back().energyMeVPerU;
}

AtInterpolatedValues AtLookup::interpolate(int atomicNumber,
                                           int massNumber,
                                           double energyMeVPerU) const {
    const AtTable& table = tableForIon(atomicNumber, massNumber);
    const std::vector<AtRow>& rows = table.rows;

    if (rows.empty()) {
        throw std::runtime_error("AT lookup interpolation called on empty LUT table.");
    }

    auto it = std::lower_bound(
        rows.begin(), rows.end(), energyMeVPerU,
        [](const AtRow& row, double value) {
            return row.energyMeVPerU < value;
        });

    if (it == rows.begin()) {
        const AtRow& row = *it;
        return AtInterpolatedValues{
            row.letKeVPerUm,
            row.yFKeVPerUm,
            row.yDKeVPerUm,
            row.yStarKeVPerUm,
            row.zFGy,
            row.zDGy,
            row.zStarGy,
            row.energyMeVPerU,
            row.energyMeVPerU,
            1.0,
            0.0
        };
    }

    if (it == rows.end()) {
        const AtRow& row = rows.back();
        return AtInterpolatedValues{
            row.letKeVPerUm,
            row.yFKeVPerUm,
            row.yDKeVPerUm,
            row.yStarKeVPerUm,
            row.zFGy,
            row.zDGy,
            row.zStarGy,
            row.energyMeVPerU,
            row.energyMeVPerU,
            1.0,
            0.0
        };
    }

    const AtRow& upper = *it;
    const AtRow& lower = *(it - 1);
    const double span = upper.energyMeVPerU - lower.energyMeVPerU;

    if (span <= 0.0) {
        throw std::runtime_error("Invalid AT LUT energy spacing encountered.");
    }

    const double upperWeight = (energyMeVPerU - lower.energyMeVPerU) / span;
    const double lowerWeight = 1.0 - upperWeight;

    return AtInterpolatedValues{
        interpolateLinear(lower.letKeVPerUm, upper.letKeVPerUm,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.yFKeVPerUm, upper.yFKeVPerUm,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.yDKeVPerUm, upper.yDKeVPerUm,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.yStarKeVPerUm, upper.yStarKeVPerUm,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.zFGy, upper.zFGy,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.zDGy, upper.zDGy,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.zStarGy, upper.zStarGy,
                          lowerWeight, upperWeight),
        lower.energyMeVPerU,
        upper.energyMeVPerU,
        lowerWeight,
        upperWeight
    };
}

const std::vector<AtTable>& AtLookup::tables() const {
    return tables_;
}

const AtTable& AtLookup::tableForIon(int atomicNumber, int massNumber) const {
    auto it = std::find_if(
        tables_.begin(), tables_.end(),
        [atomicNumber, massNumber](const AtTable& table) {
            return table.atomicNumber == atomicNumber &&
                   table.massNumber == massNumber;
        });

    if (it == tables_.end()) {
        throw std::runtime_error(
            "No AT LUT table loaded for Z=" + std::to_string(atomicNumber) +
            ", A=" + std::to_string(massNumber));
    }

    return *it;
}
