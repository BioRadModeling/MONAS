#include "AmfGeometryDeriver.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Bounds {
    double minX{0.0};
    double maxX{0.0};
    double minY{0.0};
    double maxY{0.0};
    double minZ{0.0};
    double maxZ{0.0};
    std::size_t rows{0};
};

struct Vec3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};
};

std::string trim(const std::string& value) {
    const std::size_t begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) {
        return "";
    }
    const std::size_t end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

std::string stripInlineComment(const std::string& line) {
    bool inQuote = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '"') {
            inQuote = !inQuote;
        } else if (line[i] == '#' && !inQuote) {
            return line.substr(0, i);
        }
    }
    return line;
}

std::string stripTypePrefix(std::string key) {
    key = trim(key);
    const std::size_t colon = key.find(':');
    if (colon != std::string::npos) {
        key = key.substr(colon + 1);
    }
    return trim(key);
}

std::string unquote(std::string value) {
    value = trim(value);
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
        return value.substr(1, value.size() - 2);
    }
    return value;
}

bool isNumber(const std::string& token) {
    char* end = nullptr;
    std::strtod(token.c_str(), &end);
    return end != token.c_str() && *end == '\0';
}

bool isUnit(const std::string& token) {
    return token == "mm" || token == "cm" || token == "m" || token == "um" ||
           token == "nm";
}

double unitToMm(const std::string& unit) {
    if (unit == "mm") {
        return 1.0;
    }
    if (unit == "cm") {
        return 10.0;
    }
    if (unit == "m") {
        return 1000.0;
    }
    if (unit == "um") {
        return 0.001;
    }
    if (unit == "nm") {
        return 0.000001;
    }
    throw std::runtime_error("Unsupported TOPAS length unit: " + unit);
}

std::vector<std::string> splitTokens(const std::string& value) {
    std::istringstream stream(value);
    std::vector<std::string> tokens;
    std::string token;
    while (stream >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

class TopasParameterFile {
public:
    explicit TopasParameterFile(const fs::path& path) {
        std::ifstream input(path);
        if (!input) {
            throw std::runtime_error("Failed to open TOPAS input file: " +
                                     path.string());
        }

        std::vector<std::string> unresolvedLines;
        std::string line;
        int lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;
            const std::string withoutComment = trim(stripInlineComment(line));
            if (withoutComment.empty()) {
                continue;
            }

            if (withoutComment.find("${") != std::string::npos) {
                unresolvedLines.push_back(
                    std::to_string(lineNumber) + ": " + withoutComment);
                continue;
            }

            const std::size_t equals = withoutComment.find('=');
            if (equals == std::string::npos) {
                continue;
            }

            const std::string key = stripTypePrefix(withoutComment.substr(0, equals));
            const std::string value = trim(withoutComment.substr(equals + 1));
            if (!key.empty()) {
                values_[key] = value;
            }
        }

        if (!unresolvedLines.empty()) {
            std::ostringstream message;
            message
                << "Unresolved TOPAS parameter handle(s) found in "
                << path
                << ". Replace placeholders such as ${DEPTH} or ${LATERAL} "
                << "with concrete values and try again.";
            for (const std::string& unresolvedLine : unresolvedLines) {
                message << "\n  " << unresolvedLine;
            }
            throw std::runtime_error(message.str());
        }
    }

    bool has(const std::string& key) const {
        return values_.find(key) != values_.end();
    }

    std::vector<std::string> keys() const {
        std::vector<std::string> result;
        result.reserve(values_.size());
        for (const auto& entry : values_) {
            result.push_back(entry.first);
        }
        return result;
    }

    std::string stringValue(const std::string& key) const {
        const auto it = values_.find(key);
        if (it == values_.end()) {
            throw std::runtime_error("Required TOPAS parameter missing: " + key);
        }
        return unquote(it->second);
    }

    std::string optionalStringValue(const std::string& key,
                                    const std::string& fallback) const {
        const auto it = values_.find(key);
        if (it == values_.end()) {
            return fallback;
        }
        return unquote(it->second);
    }

    double lengthMm(const std::string& key) {
        const std::string normalizedKey = stripTypePrefix(key);
        const auto cached = lengthCache_.find(normalizedKey);
        if (cached != lengthCache_.end()) {
            return cached->second;
        }

        if (activeKeys_.find(normalizedKey) != activeKeys_.end()) {
            throw std::runtime_error("Circular TOPAS parameter reference: " +
                                     normalizedKey);
        }

        const auto it = values_.find(normalizedKey);
        if (it == values_.end()) {
            throw std::runtime_error("Required TOPAS length parameter missing: " +
                                     normalizedKey);
        }

        activeKeys_.insert(normalizedKey);
        const double value = evaluateLengthExpression(it->second);
        activeKeys_.erase(normalizedKey);
        lengthCache_[normalizedKey] = value;
        return value;
    }

private:
    double evaluateLengthExpression(const std::string& expression) {
        const std::vector<std::string> tokens = splitTokens(expression);
        if (tokens.empty()) {
            throw std::runtime_error("Empty TOPAS length expression.");
        }

        double total = 0.0;
        int sign = 1;
        bool consumedValue = false;

        for (std::size_t i = 0; i < tokens.size(); ++i) {
            const std::string& token = tokens[i];
            if (token == "+") {
                sign = 1;
                consumedValue = false;
                continue;
            }
            if (token == "-") {
                sign = -1;
                consumedValue = false;
                continue;
            }
            if (isUnit(token)) {
                if (!consumedValue) {
                    throw std::runtime_error(
                        "Unit without preceding value in TOPAS expression: " +
                        expression);
                }
                continue;
            }

            double valueMm = 0.0;
            if (isNumber(token)) {
                if (i + 1 >= tokens.size() || !isUnit(tokens[i + 1])) {
                    if (std::abs(std::stod(token)) < 1e-15) {
                        valueMm = 0.0;
                    } else {
                        throw std::runtime_error(
                            "Numeric TOPAS length is missing a unit: " +
                            expression);
                    }
                } else {
                    valueMm = std::stod(token) * unitToMm(tokens[i + 1]);
                    ++i;
                }
            } else {
                valueMm = lengthMm(token);
            }

            total += sign * valueMm;
            sign = 1;
            consumedValue = true;
        }

        return total;
    }

    std::map<std::string, std::string> values_;
    std::map<std::string, double> lengthCache_;
    std::set<std::string> activeKeys_;
};

Bounds readPhaseSpaceBounds(const fs::path& phaseSpacePath) {
    std::ifstream input(phaseSpacePath);
    if (!input) {
        throw std::runtime_error("Failed to open phase-space file: " +
                                 phaseSpacePath.string());
    }

    Bounds bounds;
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream stream(line);
        double xCm = 0.0;
        double yCm = 0.0;
        double zCm = 0.0;
        if (!(stream >> xCm >> yCm >> zCm)) {
            continue;
        }

        const double xMm = xCm * 10.0;
        const double yMm = yCm * 10.0;
        const double zMm = zCm * 10.0;

        if (bounds.rows == 0) {
            bounds.minX = bounds.maxX = xMm;
            bounds.minY = bounds.maxY = yMm;
            bounds.minZ = bounds.maxZ = zMm;
        } else {
            bounds.minX = std::min(bounds.minX, xMm);
            bounds.maxX = std::max(bounds.maxX, xMm);
            bounds.minY = std::min(bounds.minY, yMm);
            bounds.maxY = std::max(bounds.maxY, yMm);
            bounds.minZ = std::min(bounds.minZ, zMm);
            bounds.maxZ = std::max(bounds.maxZ, zMm);
        }
        ++bounds.rows;
    }

    if (bounds.rows == 0) {
        throw std::runtime_error("No numeric rows found in phase-space file: " +
                                 phaseSpacePath.string());
    }

    return bounds;
}

std::vector<std::string> ancestorChain(TopasParameterFile& topas,
                                       const std::string& component) {
    std::vector<std::string> chain;
    std::string current = component;
    std::set<std::string> seen;

    while (!current.empty()) {
        if (seen.find(current) != seen.end()) {
            throw std::runtime_error("Cyclic TOPAS geometry parent chain at: " +
                                     current);
        }
        seen.insert(current);
        chain.push_back(current);

        if (current == "World") {
            break;
        }

        const std::string parentKey = "Ge/" + current + "/Parent";
        if (!topas.has(parentKey)) {
            break;
        }
        current = topas.stringValue(parentKey);
    }

    return chain;
}

double optionalLengthMm(TopasParameterFile& topas,
                        const std::string& key,
                        double fallback) {
    if (!topas.has(key)) {
        return fallback;
    }
    return topas.lengthMm(key);
}

void failOnNonZeroRotation(TopasParameterFile& topas,
                           const std::vector<std::string>& chain) {
    for (const std::string& component : chain) {
        for (const char axis : {'X', 'Y', 'Z'}) {
            const std::string key = "Ge/" + component + "/Rot" + axis;
            if (!topas.has(key)) {
                continue;
            }
            const std::vector<std::string> tokens =
                splitTokens(topas.optionalStringValue(key, "0 deg"));
            const double rotation = tokens.empty() ? 0.0 : std::stod(tokens.front());
            if (std::abs(rotation) > 1e-9) {
                throw std::runtime_error(
                    "AMF replay geometry derivation does not support non-zero "
                    "rotations in the scoring ancestry: " +
                    key);
            }
        }
    }
}

Vec3 worldTranslationMm(TopasParameterFile& topas,
                        const std::vector<std::string>& chain) {
    Vec3 translation;
    for (const std::string& component : chain) {
        translation.x += optionalLengthMm(topas, "Ge/" + component + "/TransX", 0.0);
        translation.y += optionalLengthMm(topas, "Ge/" + component + "/TransY", 0.0);
        translation.z += optionalLengthMm(topas, "Ge/" + component + "/TransZ", 0.0);
    }
    return translation;
}

std::string findPhaseSpaceScorer(TopasParameterFile& topas) {
    const std::string suffix = "/Quantity";
    for (const std::string& key : topas.keys()) {
        if (key.rfind("Sc/", 0) == 0 &&
            key.size() > suffix.size() &&
            key.compare(key.size() - suffix.size(), suffix.size(), suffix) == 0 &&
            topas.stringValue(key) == "PhaseSpace") {
            return key.substr(3, key.size() - 3 - suffix.size());
        }
    }

    throw std::runtime_error(
        "No phase-space scorer found. Expected a scorer with "
        "Quantity = \"PhaseSpace\".");
}

std::string findWaterPhantom(TopasParameterFile& topas,
                             const std::vector<std::string>& chain) {
    for (const std::string& component : chain) {
        if (component == "World") {
            continue;
        }
        const std::string material =
            topas.optionalStringValue("Ge/" + component + "/Material", "");
        if (material == "G4_WATER") {
            return component;
        }
    }

    throw std::runtime_error(
        "No G4_WATER ancestor found for the phase-space scoring component.");
}

void validatePhaseSpaceBounds(const AmfConfig& config, const Bounds& bounds) {
    const double midX = (bounds.minX + bounds.maxX) / 2.0;
    const double midY = (bounds.minY + bounds.maxY) / 2.0;
    const double midZ = (bounds.minZ + bounds.maxZ) / 2.0;

    const auto check = [&](const char* axis, double expected, double observed) {
        if (std::abs(expected - observed) > config.geometryToleranceMm) {
            std::ostringstream message;
            message << "Geometry mismatch: phase-space " << axis
                    << " midpoint is " << observed
                    << " mm, but the source TOPAS scorer center is "
                    << expected << " mm. The phase-space and source TOPAS "
                    << "file must match exactly.";
            throw std::runtime_error(message.str());
        }
    };

    check("X", config.scoringTransXmm, midX);
    check("Y", config.scoringTransYmm, midY);
    check("Z", config.scoringTransZmm, midZ);
}

}  // namespace

void AmfGeometryDeriver::deriveReplayGeometry(AmfConfig& config) {
    TopasParameterFile topas(config.sourceTopasPath);

    config.worldMaterial = topas.optionalStringValue("Ge/World/Material", "Air");
    config.worldHalfLengthXmm =
        optionalLengthMm(topas, "Ge/World/HLX", config.worldHalfLengthXmm);
    config.worldHalfLengthYmm =
        optionalLengthMm(topas, "Ge/World/HLY", config.worldHalfLengthYmm);
    config.worldHalfLengthZmm =
        optionalLengthMm(topas, "Ge/World/HLZ", config.worldHalfLengthZmm);

    const std::string scorerName = findPhaseSpaceScorer(topas);
    config.phaseSpaceScorerName = scorerName;
    config.phaseSpaceComponent =
        topas.stringValue("Sc/" + scorerName + "/Component");
    config.phaseSpaceSurface =
        topas.optionalStringValue("Sc/" + scorerName + "/Surface", "");

    const std::vector<std::string> scoringChain =
        ancestorChain(topas, config.phaseSpaceComponent);
    failOnNonZeroRotation(topas, scoringChain);
    const Vec3 sourceCenter = worldTranslationMm(topas, scoringChain);

    config.scoringComponent = "AMFScoringVolume";
    config.scoringMaterial = "G4_WATER";
    config.scoringTransXmm = sourceCenter.x;
    config.scoringTransYmm = sourceCenter.y;
    config.scoringTransZmm = sourceCenter.z;

    const std::string phantomComponent = findWaterPhantom(topas, scoringChain);
    const std::vector<std::string> phantomChain =
        ancestorChain(topas, phantomComponent);
    failOnNonZeroRotation(topas, phantomChain);
    const Vec3 phantomCenter = worldTranslationMm(topas, phantomChain);

    config.phantomComponent = phantomComponent;
    config.phantomMaterial = "G4_WATER";
    config.phantomHalfLengthXmm =
        topas.lengthMm("Ge/" + phantomComponent + "/HLX");
    config.phantomHalfLengthYmm =
        topas.lengthMm("Ge/" + phantomComponent + "/HLY");
    config.phantomHalfLengthZmm =
        topas.lengthMm("Ge/" + phantomComponent + "/HLZ");
    config.phantomTransXmm = phantomCenter.x;
    config.phantomTransYmm = phantomCenter.y;
    config.phantomTransZmm = phantomCenter.z;

    const Bounds bounds =
        readPhaseSpaceBounds(config.phaseSpaceBasePath.string() + ".phsp");
    config.hasPhaseSpaceBounds = true;
    config.phaseSpaceMinXmm = bounds.minX;
    config.phaseSpaceMaxXmm = bounds.maxX;
    config.phaseSpaceMinYmm = bounds.minY;
    config.phaseSpaceMaxYmm = bounds.maxY;
    config.phaseSpaceMinZmm = bounds.minZ;
    config.phaseSpaceMaxZmm = bounds.maxZ;

    validatePhaseSpaceBounds(config, bounds);
}
