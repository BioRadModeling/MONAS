#pragma once

#include <string>
#include <unordered_map>
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
    SpectrumAccumulator();

    void addContribution(const LookupTable& table, double multiplicity);
    PolySpectrum finalize() const;

private:
    std::vector<double> targetEdges_;
    std::vector<double> yCenters_;
    std::vector<double> binWidths_;

    // Cache of precomputed monoenergetic contribution vectors,
    // keyed by unique source file path.
    std::unordered_map<std::string, std::vector<double>> precomputedRawFyCache_;

    // Accumulated numerator on the rebinned grid.
    std::vector<double> numerator_;
    double denominator_{0.0};

    void buildTargetYBins();

    std::vector<double> rebinDiscreteValuesToTargetGrid(
        const std::vector<double>& srcY,
        const std::vector<double>& srcValues) const;

    std::vector<double> buildPrecomputedRawFy(const LookupTable& table) const;
    std::vector<double> buildDeCunhaPrecomputedRawFy(
        const LookupTable& table) const;
    std::vector<double> buildCartechiniPrecomputedRawFy(
        const LookupTable& table) const;
};