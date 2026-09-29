#pragma once

#include <filesystem>
#include <string>
#include <vector>

struct AtRow {
    double energyMeVPerU{0.0};
    double letKeVPerUm{0.0};
    double yFKeVPerUm{0.0};
    double yDKeVPerUm{0.0};
    double yStarKeVPerUm{0.0};
    double zFGy{0.0};
    double zDGy{0.0};
    double zStarGy{0.0};
};

struct AtInterpolatedValues {
    double letKeVPerUm{0.0};
    double yFKeVPerUm{0.0};
    double yDKeVPerUm{0.0};
    double yStarKeVPerUm{0.0};
    double zFGy{0.0};
    double zDGy{0.0};
    double zStarGy{0.0};
    double lowerEnergyMeVPerU{0.0};
    double upperEnergyMeVPerU{0.0};
    double lowerWeight{1.0};
    double upperWeight{0.0};
};

struct AtTable {
    int atomicNumber{0};
    int massNumber{0};
    std::vector<AtRow> rows;
    std::string sourceFile;
};

class AtLookup {
public:
    static AtLookup loadFromDirectory(const std::filesystem::path& directory);

    bool hasIon(int atomicNumber, int massNumber) const;
    double minEnergyMeVPerU(int atomicNumber, int massNumber) const;
    double maxEnergyMeVPerU(int atomicNumber, int massNumber) const;
    bool isEnergyInRange(int atomicNumber,
                         int massNumber,
                         double energyMeVPerU) const;
    AtInterpolatedValues interpolate(int atomicNumber,
                                     int massNumber,
                                     double energyMeVPerU) const;

    const std::vector<AtTable>& tables() const;

private:
    explicit AtLookup(std::vector<AtTable> tables);

    const AtTable& tableForIon(int atomicNumber, int massNumber) const;

    std::vector<AtTable> tables_;
};
