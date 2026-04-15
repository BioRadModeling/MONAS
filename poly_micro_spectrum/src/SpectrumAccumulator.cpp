#include "SpectrumAccumulator.h"

#include <cmath>
#include <numeric>
#include <stdexcept>

SpectrumAccumulator::SpectrumAccumulator(std::vector<double> yRef)
    : yRef_(std::move(yRef)),
      weightedCountsPerPrimary_(yRef_.size(), 0.0) {}

void SpectrumAccumulator::addContribution(const LookupTable& table, double particleWeight) {
    if (table.yLowerEdges().size() != yRef_.size()) {
        throw std::runtime_error("Lookup table y-grid size does not match accumulator.");
    }

    if (particleWeight <= 0.0) {
        return;
    }

    const auto& ny = table.nY();
    const double totalNy = std::accumulate(ny.begin(), ny.end(), 0.0);

    if (totalNy <= 0.0) {
        throw std::runtime_error("Lookup table has non-positive total N(y).");
    }

    // monoCountsPerPrimary(bin) = f(bin) * Ncpp
    // where f(bin) = N(y)_bin / sum_bin N(y)
    for (std::size_t i = 0; i < ny.size(); ++i) {
        const double fBin = ny[i] / totalNy;
        const double countsPerPrimaryBin = fBin * table.ncpp();
        weightedCountsPerPrimary_[i] += particleWeight * countsPerPrimaryBin;
    }

    totalWeight_ += particleWeight;
}

PolySpectrum SpectrumAccumulator::finalize() const {
    if (totalWeight_ <= 0.0) {
        throw std::runtime_error("Cannot finalize spectrum with zero total particle weight.");
    }

    PolySpectrum out;
    out.y = yRef_;
    out.countsPerPrimary.resize(yRef_.size(), 0.0);
    out.f.resize(yRef_.size(), 0.0);
    out.yf.resize(yRef_.size(), 0.0);
    out.d.resize(yRef_.size(), 0.0);
    out.yd.resize(yRef_.size(), 0.0);

    for (std::size_t i = 0; i < yRef_.size(); ++i) {
        out.countsPerPrimary[i] = weightedCountsPerPrimary_[i] / totalWeight_;
    }

    const double ncppPoly =
        std::accumulate(out.countsPerPrimary.begin(), out.countsPerPrimary.end(), 0.0);

    if (ncppPoly <= 0.0) {
        throw std::runtime_error("Polyenergetic spectrum has non-positive total counts per primary.");
    }

    for (std::size_t i = 0; i < yRef_.size(); ++i) {
        out.f[i] = out.countsPerPrimary[i] / ncppPoly;
        out.yf[i] = out.y[i] * out.f[i];
    }

    const double sumYf = std::accumulate(out.yf.begin(), out.yf.end(), 0.0);

    if (sumYf <= 0.0) {
        throw std::runtime_error("Polyenergetic spectrum has non-positive sum(y*f(y)).");
    }

    for (std::size_t i = 0; i < yRef_.size(); ++i) {
        out.d[i] = out.yf[i] / sumYf;
        out.yd[i] = out.y[i] * out.d[i];
    }

    return out;
}