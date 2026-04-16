#pragma once

#include <filesystem>
#include <string>
#include <vector>

class LookupTable {
public:
    LookupTable() = default;

    LookupTable(double monoEnergyMeV,
                std::vector<double> yLowerEdges,
                std::vector<double> nY,
                double ncpp,
                std::string sourceFile);

    static LookupTable loadFromCsv(const std::filesystem::path& csvPath,
                                   double monoEnergyMeV);

    double monoEnergyMeV() const;
    double ncpp() const;

    const std::vector<double>& yLowerEdges() const;
    const std::vector<double>& nY() const;
    const std::string& sourceFile() const;

private:
    double monoEnergyMeV_{0.0};
    std::vector<double> yLowerEdges_;
    std::vector<double> nY_;
    double ncpp_{0.0};
    std::string sourceFile_;
};