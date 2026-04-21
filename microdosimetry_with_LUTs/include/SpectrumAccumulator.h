#pragma once

#include <cstddef>
#include <cstdint>
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
    explicit SpectrumAccumulator(
        const LookupLibrary& library,
        std::size_t rebinSamples = 1000000,
        std::uint64_t rebinSeed = 0x5A17C3E4ULL);

    void addContributionByIndex(std::size_t tableIndex, double multiplicity);
    PolySpectrum finalize() const;

private:
    // Common interval representation used by both rebinners.
    // "masses" means spectral weight carried by each source interval.
    struct SourceIntervalData {
        std::vector<double> lower;
        std::vector<double> upper;
        std::vector<double> masses;
    };

    std::vector<double> targetEdges_;
    std::vector<double> yCenters_;
    std::vector<double> binWidths_;

    std::vector<std::vector<double>> precomputedRawFyByTable_;

    std::vector<double> numerator_;
    double denominator_{0.0};

    std::size_t rebinSamples_{1000000};
    std::uint64_t rebinSeed_{0x5A17C3E4ULL};

    void buildTargetYBins();

    SourceIntervalData buildSourceIntervals(const LookupTable& table) const;

    std::vector<double> sampleCountsToTargetGrid(
        const SourceIntervalData& intervals,
        std::size_t nSamples,
        std::uint64_t seed) const;

    // Exact log-bin overlap rebinning used mainly for DeCunha N(y) histograms.
    std::vector<double> rebinCountsToTargetGridStochastic(
        const LookupTable& table,
        std::size_t nSamples,
        std::uint64_t seed) const;

    // Monte Carlo rebinning used mainly for Cartechini f(y) LUTs.
    std::vector<double> rebinCountsToTargetGridDeterministic(
        const LookupTable& table) const;

    // Applies the family-specific rebinning policy and stores the rebinned
    // monoenergetic spectrum on the shared target grid.
    std::vector<double> buildPrecomputedRawFy(
        const LookupTable& table,
        std::size_t tableIndex) const;
};