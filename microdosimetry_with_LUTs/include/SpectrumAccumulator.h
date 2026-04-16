#pragma once

#include <cstddef>
#include <vector>

#include "LookupLibrary.h"
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
    explicit SpectrumAccumulator(const LookupLibrary& library);

    void addContributionByIndex(std::size_t tableIndex, double multiplicity);
    PolySpectrum finalize() const;

private:
    std::vector<double> targetEdges_;
    std::vector<double> yCenters_;
    std::vector<double> binWidths_;

    // Precomputed rebinned monoenergetic contribution for each lookup table
    std::vector<std::vector<double>> precomputedRawFyByTable_;

    // Accumulated numerator on the rebinned grid
    std::vector<double> numerator_;
    double denominator_{0.0};

    void buildTargetYBins();
    std::vector<double> rebinCountsToTargetGrid(const LookupTable& table) const;
    std::vector<double> buildPrecomputedRawFy(const LookupTable& table) const;
};