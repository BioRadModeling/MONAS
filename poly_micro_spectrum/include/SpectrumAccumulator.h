#pragma once

#include <vector>

#include "LookupTable.h"

struct PolySpectrum {
    std::vector<double> y;
    std::vector<double> countsPerPrimary;  // averaged counts-per-primary per bin
    std::vector<double> f;
    std::vector<double> yf;
    std::vector<double> d;
    std::vector<double> yd;
};

class SpectrumAccumulator {
public:
    explicit SpectrumAccumulator(std::vector<double> yRef);

    void addContribution(const LookupTable& table, double particleWeight);
    PolySpectrum finalize() const;

private:
    std::vector<double> yRef_;
    std::vector<double> weightedCountsPerPrimary_;
    double totalWeight_{0.0};
};