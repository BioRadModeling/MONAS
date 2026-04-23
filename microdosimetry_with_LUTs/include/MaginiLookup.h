#pragma once

#include <filesystem>
#include <string>
#include <vector>

struct MaginiRow {
    double energyMeV{0.0};
    double yFLutKeVPerUm{0.0};
    double yStarLutKeVPerUm{0.0};
    double yDLutKeVPerUm{0.0};
};

struct MaginiInterpolatedValues {
    double yFLutKeVPerUm{0.0};
    double yStarLutKeVPerUm{0.0};
    double yDLutKeVPerUm{0.0};
    double lowerEnergyMeV{0.0};
    double upperEnergyMeV{0.0};
    double lowerWeight{1.0};
    double upperWeight{0.0};
};

class MaginiLookup {
public:
    static MaginiLookup loadFromCsv(const std::filesystem::path& csvPath);

    MaginiInterpolatedValues interpolate(double energyMeV) const;

    const std::vector<MaginiRow>& rows() const;
    const std::string& sourceFile() const;

private:
    explicit MaginiLookup(std::vector<MaginiRow> rows, std::string sourceFile);

    std::vector<MaginiRow> rows_;
    std::string sourceFile_;
};
