#pragma once

#include <vector>

#include "LookupTable.h"

struct PolySpectrum {
    std::vector<double> y;
    std::vector<double> f;
    std::vector<double> yf;
    std::vector<double> d;
    std::vector<double> yd;
};

class SpectrumAccumulator {
public:
    explicit SpectrumAccumulator(std::vector<double> yLowerEdges);

    void addContribution(const LookupTable& table, double particleWeight);
    PolySpectrum finalize() const;

private:
    // Native lookup-table y grid from the CSV first column
    std::vector<double> nativeYLowerEdges_;

    // R-style rebinned y grid
    std::vector<double> targetEdges_;
    std::vector<double> yCenters_;
    std::vector<double> binWidths_;

    // Weighted sum numerator on the rebinned grid
    std::vector<double> numerator_;
    double denominator_{0.0};

    void buildTargetYBins();
    std::vector<double> rebinCountsToTargetGrid(const LookupTable& table) const;
};