#include "SpectrumAccumulator.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

SpectrumAccumulator::SpectrumAccumulator() {
    buildTargetYBins();
    numerator_.assign(yCenters_.size(), 0.0);
}

void SpectrumAccumulator::buildTargetYBins() {
    // Keep the current rebinned-y behavior:
    // 100 logarithmically spaced edges from 0.1 to 1000 keV/um
    const std::size_t nEdges = 100;
    const double yMin = 0.1;
    const double yMax = 1000.0;

    targetEdges_.resize(nEdges);
    yCenters_.resize(nEdges - 1);
    binWidths_.resize(nEdges - 1);

    const double logMin = std::log10(yMin);
    const double logMax = std::log10(yMax);

    for (std::size_t i = 0; i < nEdges; ++i) {
        const double t =
            static_cast<double>(i) / static_cast<double>(nEdges - 1);
        targetEdges_[i] = std::pow(10.0, logMin + t * (logMax - logMin));
    }

    for (std::size_t i = 0; i < nEdges - 1; ++i) {
        yCenters_[i] = 0.5 * (targetEdges_[i] + targetEdges_[i + 1]);
        binWidths_[i] = targetEdges_[i + 1] - targetEdges_[i];

        if (binWidths_[i] <= 0.0) {
            throw std::runtime_error(
                "Non-positive rebinned y-bin width encountered.");
        }
    }
}

std::vector<double> SpectrumAccumulator::rebinDiscreteValuesToTargetGrid(
    const std::vector<double>& srcY,
    const std::vector<double>& srcValues) const {
    if (srcY.size() != srcValues.size()) {
        throw std::runtime_error(
            "Source y-grid and source value vector size mismatch.");
    }

    std::vector<double> rebinnedValues(yCenters_.size(), 0.0);

    for (std::size_t i = 0; i < srcY.size(); ++i) {
        const double yVal = srcY[i];
        const double value = srcValues[i];

        if (value <= 0.0) {
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

        if (idx < rebinnedValues.size()) {
            rebinnedValues[idx] += value;
        }
    }

    return rebinnedValues;
}

std::vector<double> SpectrumAccumulator::buildDeCunhaPrecomputedRawFy(
    const LookupTable& table) const {
    const auto& srcY = table.yLowerEdges();
    const auto& srcN = table.nY();

    if (srcY.size() != srcN.size()) {
        throw std::runtime_error(
            "DeCunha lookup table y-grid and N(y) size mismatch.");
    }

    const std::vector<double> rebinnedCounts =
        rebinDiscreteValuesToTargetGrid(srcY, srcN);

    const double totalRebinnedCounts =
        std::accumulate(rebinnedCounts.begin(), rebinnedCounts.end(), 0.0);

    std::vector<double> rawFy(yCenters_.size(), 0.0);

    if (totalRebinnedCounts <= 0.0) {
        return rawFy;
    }

    const double ncpp = table.ncpp();
    if (!std::isfinite(ncpp) || ncpp <= 0.0) {
        throw std::runtime_error(
            "Invalid Ncpp encountered for DeCunha lookup table: " +
            table.sourceFile());
    }

    for (std::size_t i = 0; i < rebinnedCounts.size(); ++i) {
        const double fMono =
            (rebinnedCounts[i] / totalRebinnedCounts) / binWidths_[i];
        rawFy[i] = fMono * ncpp;
    }

    return rawFy;
}

std::vector<double> SpectrumAccumulator::buildCartechiniPrecomputedRawFy(
    const LookupTable& table) const {
    const auto& srcY = table.yLowerEdges();
    const auto& srcFy = table.fY();

    if (srcY.size() != srcFy.size()) {
        throw std::runtime_error(
            "Cartechini lookup table y-grid and f(y) size mismatch.");
    }

    // Cartechini provides monoenergetic f(y) directly.
    // Rebin it onto the target grid and use it directly.
    return rebinDiscreteValuesToTargetGrid(srcY, srcFy);
}

std::vector<double> SpectrumAccumulator::buildPrecomputedRawFy(
    const LookupTable& table) const {
    switch (table.spectrumKind()) {
        case LookupSpectrumKind::DeCunhaRawCounts:
            return buildDeCunhaPrecomputedRawFy(table);

        case LookupSpectrumKind::CartechiniFy:
            return buildCartechiniPrecomputedRawFy(table);

        default:
            throw std::runtime_error(
                "Unsupported lookup spectrum kind encountered.");
    }
}

void SpectrumAccumulator::addContribution(const LookupTable& table,
                                          double multiplicity) {
    if (multiplicity <= 0.0) {
        return;
    }

    const std::string cacheKey = table.sourceFile();
    auto it = precomputedRawFyCache_.find(cacheKey);

    if (it == precomputedRawFyCache_.end()) {
        const std::vector<double> rawFy = buildPrecomputedRawFy(table);
        it = precomputedRawFyCache_.emplace(cacheKey, rawFy).first;
    }

    const auto& rawFy = it->second;
    if (rawFy.size() != numerator_.size()) {
        throw std::runtime_error(
            "Cached monoenergetic contribution size mismatch.");
    }

    for (std::size_t i = 0; i < rawFy.size(); ++i) {
        numerator_[i] += multiplicity * rawFy[i];
    }

    denominator_ += multiplicity;
}

PolySpectrum SpectrumAccumulator::finalize() const {
    if (denominator_ <= 0.0) {
        throw std::runtime_error(
            "Cannot finalize spectrum with zero total weight.");
    }

    PolySpectrum out;
    out.y = yCenters_;
    out.f.resize(yCenters_.size(), 0.0);
    out.yf.resize(yCenters_.size(), 0.0);
    out.d.resize(yCenters_.size(), 0.0);
    out.yd.resize(yCenters_.size(), 0.0);

    std::vector<double> rawFy(yCenters_.size(), 0.0);
    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        rawFy[i] = numerator_[i] / denominator_;
    }

    const double C =
        std::log(10.0) *
        (std::log10(binWidths_[1]) - std::log10(binWidths_[0]));

    if (C <= 0.0) {
        throw std::runtime_error(
            "Invalid logarithmic-bin normalization factor C.");
    }

    double fyNormDen = 0.0;
    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        fyNormDen += out.y[i] * rawFy[i];
    }
    fyNormDen *= C;

    if (fyNormDen <= 0.0) {
        throw std::runtime_error(
            "Polyenergetic f(y) has non-positive normalization.");
    }

    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        out.f[i] = rawFy[i] / fyNormDen;
        out.yf[i] = out.y[i] * out.f[i];
    }

    double dNormDen = 0.0;
    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        dNormDen += out.y[i] * out.yf[i];
    }
    dNormDen *= C;

    if (dNormDen <= 0.0) {
        throw std::runtime_error(
            "Polyenergetic d(y) has non-positive normalization.");
    }

    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        out.d[i] = out.yf[i] / dNormDen;
        out.yd[i] = out.y[i] * out.d[i];
    }

    return out;
}