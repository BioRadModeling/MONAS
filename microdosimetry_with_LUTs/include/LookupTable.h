#pragma once

#include <filesystem>
#include <string>
#include <vector>

enum class LookupSpectrumKind {
    DeCunhaRawCounts,
    CartechiniFy
};

class LookupTable {
public:
    LookupTable() = default;

    LookupTable(double monoEnergyMeV,
                LookupSpectrumKind spectrumKind,
                std::vector<double> yLowerEdges,
                std::vector<double> nY,
                std::vector<double> fY,
                double ncpp,
                std::string sourceFile);

    static LookupTable loadFromCsv(const std::filesystem::path& csvPath,
                                   double monoEnergyMeV);

    static LookupTable loadFromCartechiniYSpec(
        const std::filesystem::path& ySpecPath,
        double monoEnergyMeV);

    double monoEnergyMeV() const;
    LookupSpectrumKind spectrumKind() const;
    double ncpp() const;
    const std::vector<double>& yLowerEdges() const;
    const std::vector<double>& nY() const;
    const std::vector<double>& fY() const;
    const std::string& sourceFile() const;

private:
    double monoEnergyMeV_{0.0};
    LookupSpectrumKind spectrumKind_{LookupSpectrumKind::DeCunhaRawCounts};
    std::vector<double> yLowerEdges_;
    std::vector<double> nY_;
    std::vector<double> fY_;
    double ncpp_{0.0};
    std::string sourceFile_;
};