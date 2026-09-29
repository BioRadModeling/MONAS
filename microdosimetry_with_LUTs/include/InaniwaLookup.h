#pragma once

#include <filesystem>
#include <string>
#include <vector>

struct InaniwaRow {
    double energyMeVPerU{0.0};
    double zdDMeanGy{0.0};
    double zdDStarMeanGy{0.0};
    double znDMeanGy{0.0};
};

struct InaniwaInterpolatedValues {
    double zdDMeanGy{0.0};
    double zdDStarMeanGy{0.0};
    double znDMeanGy{0.0};
    double lowerEnergyMeVPerU{0.0};
    double upperEnergyMeVPerU{0.0};
    double lowerWeight{1.0};
    double upperWeight{0.0};
};

struct InaniwaTable {
    int atomicNumber{0};
    std::vector<InaniwaRow> rows;
    std::string sourceFile;
};

class InaniwaLookup {
public:
    static InaniwaLookup loadFromDirectory(const std::filesystem::path& directory);

    bool hasAtomicNumber(int atomicNumber) const;
    InaniwaInterpolatedValues interpolate(int atomicNumber,
                                          double energyMeVPerU) const;

    const std::vector<InaniwaTable>& tables() const;

private:
    explicit InaniwaLookup(std::vector<InaniwaTable> tables);

    const InaniwaTable& tableForAtomicNumber(int atomicNumber) const;

    std::vector<InaniwaTable> tables_;
};
