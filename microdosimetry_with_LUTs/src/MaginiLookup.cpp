#include "MaginiLookup.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <initializer_list>
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

std::string trim(std::string value) {
    auto isNotSpace = [](unsigned char c) {
        return !std::isspace(c);
    };

    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(), isNotSpace));
    value.erase(std::find_if(value.rbegin(), value.rend(), isNotSpace).base(),
                value.end());
    return value;
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

std::size_t findHeaderIndex(const std::vector<std::string>& headers,
                            std::initializer_list<const char*> names) {
    for (std::size_t i = 0; i < headers.size(); ++i) {
        const std::string header = trim(headers[i]);
        for (const char* name : names) {
            if (header == name) {
                return i;
            }
        }
    }

    return headers.size();
}

struct MaginiColumnIndices {
    std::size_t energy{0};
    std::size_t yF{1};
    std::size_t yStar{2};
    std::size_t yD{3};
    std::size_t yFMaxError{static_cast<std::size_t>(-1)};
    std::size_t yDMaxError{static_cast<std::size_t>(-1)};
};

MaginiColumnIndices resolveHeaderColumns(const std::vector<std::string>& headers,
                                         const std::filesystem::path& csvPath) {
    MaginiColumnIndices indices;
    indices.energy = findHeaderIndex(headers, {"E", "E_i_MeV"});
    indices.yF = findHeaderIndex(headers, {"yF", "y_F_LUT_keV_per_um"});
    indices.yD = findHeaderIndex(headers, {"yD", "y_D_LUT_keV_per_um"});
    indices.yStar = findHeaderIndex(headers, {"yS", "y_star_LUT_keV_per_um"});
    indices.yFMaxError = findHeaderIndex(headers, {"y_fmax_err"});
    indices.yDMaxError = findHeaderIndex(headers, {"y_dmax_err"});

    if (indices.energy == headers.size() ||
        indices.yF == headers.size() ||
        indices.yD == headers.size() ||
        indices.yStar == headers.size()) {
        throw std::runtime_error(
            "Magini CSV missing one or more required columns "
            "(E/yF/yD/yS or legacy equivalents): " + csvPath.string());
    }

    return indices;
}

bool parseMaginiRow(const std::vector<std::string>& fields,
                    const MaginiColumnIndices& columns,
                    MaginiRow& row) {
    const std::size_t maxIndex =
        std::max({columns.energy, columns.yF, columns.yD, columns.yStar});
    if (fields.size() <= maxIndex) {
        return false;
    }

    if (!tryParseDouble(fields[columns.energy], row.energyMeV) ||
        !tryParseDouble(fields[columns.yF], row.yFLutKeVPerUm) ||
        !tryParseDouble(fields[columns.yD], row.yDLutKeVPerUm) ||
        !tryParseDouble(fields[columns.yStar], row.yStarLutKeVPerUm)) {
        return false;
    }

    if (columns.yFMaxError < fields.size() &&
        !tryParseDouble(fields[columns.yFMaxError], row.yFMaxErrorKeVPerUm)) {
        return false;
    }

    if (columns.yDMaxError < fields.size() &&
        !tryParseDouble(fields[columns.yDMaxError], row.yDMaxErrorKeVPerUm)) {
        return false;
    }

    return true;
}

}  // namespace

MaginiLookup::MaginiLookup(std::vector<MaginiRow> rows, std::string sourceFile)
    : rows_(std::move(rows)),
      sourceFile_(std::move(sourceFile)) {}

MaginiLookup MaginiLookup::loadFromCsv(const std::filesystem::path& csvPath) {
    std::ifstream in(csvPath);
    if (!in) {
        throw std::runtime_error("Failed to open Magini CSV: " + csvPath.string());
    }

    std::vector<MaginiRow> rows;
    std::string line;
    MaginiColumnIndices columns;
    bool sawData = false;

    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        const auto fields = splitCsvLine(line);
        if (fields.size() < 4) {
            continue;
        }

        MaginiRow row;
        if (!sawData && !parseMaginiRow(fields, columns, row)) {
            columns = resolveHeaderColumns(fields, csvPath);
            sawData = true;
            continue;
        }

        sawData = true;
        if (!parseMaginiRow(fields, columns, row)) {
            continue;
        }

        rows.push_back(row);
    }

    if (rows.empty()) {
        throw std::runtime_error(
            "No usable Magini LUT rows found in: " + csvPath.string());
    }

    std::sort(rows.begin(), rows.end(),
              [](const MaginiRow& a, const MaginiRow& b) {
                  return a.energyMeV < b.energyMeV;
              });

    for (std::size_t i = 1; i < rows.size(); ++i) {
        if (rows[i - 1].energyMeV == rows[i].energyMeV) {
            throw std::runtime_error(
                "Duplicate Magini LUT energy detected in: " + csvPath.string());
        }
    }

    return MaginiLookup(rows, csvPath.string());
}

MaginiInterpolatedValues MaginiLookup::interpolate(double energyMeV) const {
    if (rows_.empty()) {
        throw std::runtime_error("Magini lookup interpolation called on empty LUT.");
    }

    auto it = std::lower_bound(
        rows_.begin(), rows_.end(), energyMeV,
        [](const MaginiRow& row, double value) {
            return row.energyMeV < value;
        });

    if (it == rows_.begin()) {
        const MaginiRow& row = *it;
        return MaginiInterpolatedValues{
            row.yFLutKeVPerUm,
            row.yStarLutKeVPerUm,
            row.yDLutKeVPerUm,
            row.yFMaxErrorKeVPerUm,
            row.yDMaxErrorKeVPerUm,
            row.energyMeV,
            row.energyMeV,
            1.0,
            0.0
        };
    }

    if (it == rows_.end()) {
        const MaginiRow& row = rows_.back();
        return MaginiInterpolatedValues{
            row.yFLutKeVPerUm,
            row.yStarLutKeVPerUm,
            row.yDLutKeVPerUm,
            row.yFMaxErrorKeVPerUm,
            row.yDMaxErrorKeVPerUm,
            row.energyMeV,
            row.energyMeV,
            1.0,
            0.0
        };
    }

    const MaginiRow& upper = *it;
    const MaginiRow& lower = *(it - 1);
    const double span = upper.energyMeV - lower.energyMeV;

    if (span <= 0.0) {
        throw std::runtime_error("Invalid Magini LUT energy spacing encountered.");
    }

    const double upperWeight = (energyMeV - lower.energyMeV) / span;
    const double lowerWeight = 1.0 - upperWeight;

    return MaginiInterpolatedValues{
        interpolateLinear(lower.yFLutKeVPerUm, upper.yFLutKeVPerUm,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.yStarLutKeVPerUm, upper.yStarLutKeVPerUm,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.yDLutKeVPerUm, upper.yDLutKeVPerUm,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.yFMaxErrorKeVPerUm, upper.yFMaxErrorKeVPerUm,
                          lowerWeight, upperWeight),
        interpolateLinear(lower.yDMaxErrorKeVPerUm, upper.yDMaxErrorKeVPerUm,
                          lowerWeight, upperWeight),
        lower.energyMeV,
        upper.energyMeV,
        lowerWeight,
        upperWeight
    };
}

const std::vector<MaginiRow>& MaginiLookup::rows() const {
    return rows_;
}

const std::string& MaginiLookup::sourceFile() const {
    return sourceFile_;
}
