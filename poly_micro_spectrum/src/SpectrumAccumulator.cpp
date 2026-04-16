#include "SpectrumAccumulator.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

SpectrumAccumulator::SpectrumAccumulator(std::vector<double> yLowerEdges)
    : nativeYLowerEdges_(std::move(yLowerEdges)) {
    buildTargetYBins();
    numerator_.assign(yCenters_.size(), 0.0);
}

void SpectrumAccumulator::buildTargetYBins() {
    // Match the R script:
    // bin <- 10^(seq(log10(0.1), log10(1000), length.out = 100))
    const std::size_t nEdges = 100;
    const double yMin = 0.1;
    const double yMax = 1000.0;

    targetEdges_.resize(nEdges);
    yCenters_.resize(nEdges - 1);
    binWidths_.resize(nEdges - 1);

    const double logMin = std::log10(yMin);
    const double logMax = std::log10(yMax);

    for (std::size_t i = 0; i < nEdges; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(nEdges - 1);
        targetEdges_[i] = std::pow(10.0, logMin + t * (logMax - logMin));
    }

    for (std::size_t i = 0; i < nEdges - 1; ++i) {
        // Match the R code's arithmetic bin centers
        yCenters_[i] = 0.5 * (targetEdges_[i] + targetEdges_[i + 1]);
        binWidths_[i] = targetEdges_[i + 1] - targetEdges_[i];

        if (binWidths_[i] <= 0.0) {
            throw std::runtime_error("Non-positive rebinned y-bin width encountered.");
        }
    }
}

std::vector<double> SpectrumAccumulator::rebinCountsToTargetGrid(
    const LookupTable& table) const {

    const auto& srcY = table.yLowerEdges();
    const auto& srcN = table.nY();

    if (srcY.size() != srcN.size()) {
        throw std::runtime_error("Lookup table y-grid and N(y) size mismatch.");
    }

    std::vector<double> rebinnedCounts(yCenters_.size(), 0.0);

    // Mimic the R code's hist(..., breaks = bin) behavior on the lookup y values.
    // Keep only values strictly inside (0.1, 1000), as in the R script.
    for (std::size_t i = 0; i < srcY.size(); ++i) {
        const double yVal = srcY[i];
        const double count = srcN[i];

        if (count <= 0.0) {
            continue;
        }

        if (!(yVal > targetEdges_.front() && yVal < targetEdges_.back())) {
            continue;
        }

        auto it = std::upper_bound(targetEdges_.begin(), targetEdges_.end(), yVal);
        if (it == targetEdges_.begin() || it == targetEdges_.end()) {
            continue;
        }

        const std::size_t idx =
            static_cast<std::size_t>(std::distance(targetEdges_.begin(), it) - 1);

        if (idx < rebinnedCounts.size()) {
            rebinnedCounts[idx] += count;
        }
    }

    return rebinnedCounts;
}

void SpectrumAccumulator::addContribution(const LookupTable& table, double particleWeight) {
    if (particleWeight <= 0.0) {
        return;
    }

    const std::vector<double> rebinnedCounts = rebinCountsToTargetGrid(table);
    const double totalRebinnedCounts =
        std::accumulate(rebinnedCounts.begin(), rebinnedCounts.end(), 0.0);

    if (totalRebinnedCounts <= 0.0) {
        return;
    }

    // Build a monoenergetic spectrum on the rebinned grid as a density per unit y.
    // This keeps the "sum monoenergetic f(y) to get polyenergetic f(y)" interpretation,
    // but now on the same y-grid used by the R script.
    for (std::size_t i = 0; i < rebinnedCounts.size(); ++i) {
        const double fMono = (rebinnedCounts[i] / totalRebinnedCounts) / binWidths_[i];
        numerator_[i] += particleWeight * fMono * table.ncpp();
    }

    denominator_ += particleWeight;
}

PolySpectrum SpectrumAccumulator::finalize() const {
    if (denominator_ <= 0.0) {
        throw std::runtime_error("Cannot finalize spectrum with zero total weight.");
    }

    PolySpectrum out;
    out.y = yCenters_;
    out.f.resize(yCenters_.size(), 0.0);
    out.yf.resize(yCenters_.size(), 0.0);
    out.d.resize(yCenters_.size(), 0.0);
    out.yd.resize(yCenters_.size(), 0.0);

    // Step 1: provisional polyenergetic f(y) from the weighted sum
    std::vector<double> rawFy(yCenters_.size(), 0.0);
    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        rawFy[i] = numerator_[i] / denominator_;
    }

    // Step 2: R-style post-summation normalization
    //
    // Match the alternative normalization branch in the R code:
    // C <- log(10) * diff(log10(h$BinWidth))[1]
    // fy_norm <- fy / (C * sum(y * fy))
    //
    if (binWidths_.size() < 2) {
        throw std::runtime_error("Need at least two rebinned y bins.");
    }

    const double C =
        std::log(10.0) * (std::log10(binWidths_[1]) - std::log10(binWidths_[0]));

    if (C <= 0.0) {
        throw std::runtime_error("Invalid logarithmic-bin normalization factor C.");
    }

    double fyNormDen = 0.0;
    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        fyNormDen += out.y[i] * rawFy[i];
    }
    fyNormDen *= C;

    if (fyNormDen <= 0.0) {
        throw std::runtime_error("Polyenergetic f(y) has non-positive normalization.");
    }

    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        out.f[i] = rawFy[i] / fyNormDen;
        out.yf[i] = out.y[i] * out.f[i];
    }

    // Step 3: R-style d(y) normalization
    //
    // Match:
    // yfy_norm <- yfy / (C * sum(y * yfy))
    // ydy <- y * yfy_norm
    //
    double dNormDen = 0.0;
    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        dNormDen += out.y[i] * out.yf[i];
    }
    dNormDen *= C;

    if (dNormDen <= 0.0) {
        throw std::runtime_error("Polyenergetic d(y) has non-positive normalization.");
    }

    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        out.d[i] = out.yf[i] / dNormDen;
        out.yd[i] = out.y[i] * out.d[i];
    }

    return out;
}