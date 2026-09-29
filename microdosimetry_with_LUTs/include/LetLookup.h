#pragma once

#include <filesystem>
#include <string>
#include <vector>

class LetLookup {
public:
    LetLookup() = default;

    LetLookup(std::string element,
              std::vector<double> energiesMeV,
              std::vector<double> letKeVPerUm,
              std::string sourceFile);

    static LetLookup loadFromFile(const std::filesystem::path& tablePath,
                                  std::string element);

    double interpolate(double energyMeV) const;

    const std::string& element() const;
    const std::vector<double>& energiesMeV() const;
    const std::vector<double>& letKeVPerUm() const;
    const std::string& sourceFile() const;

private:
    std::string element_;
    std::vector<double> energiesMeV_;
    std::vector<double> letKeVPerUm_;
    std::string sourceFile_;
};
