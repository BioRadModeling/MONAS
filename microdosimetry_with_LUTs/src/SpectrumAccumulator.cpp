#include "SpectrumAccumulator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace {

bool isFinitePositive(double x) {
    return std::isfinite(x) && x > 0.0;
}

void validateSpectrumShape(const PolySpectrum& spectrum) {
    const std::size_t nBins = spectrum.y.size();
    if (nBins < 2) {
        throw std::runtime_error("Insufficient rebinned y-grid size.");
    }

    if (spectrum.f.size() != nBins || spectrum.yf.size() != nBins ||
        spectrum.d.size() != nBins || spectrum.yd.size() != nBins) {
        throw std::runtime_error(
            "PolySpectrum vectors must all have the same size.");
    }
}

SpectrumDistributionMoments computeDistributionMoments(
    const std::vector<double>& yValues,
    const std::vector<double>& densityValues,
    const std::string& distribution,
    double normalizationFactor) {

    if (yValues.size() != densityValues.size()) {
        throw std::runtime_error(
            "Cannot compute distribution moments with mismatched vector sizes.");
    }

    SpectrumDistributionMoments moments;
    moments.distribution = distribution;

    double totalProbability = 0.0;
    for (std::size_t i = 0; i < yValues.size(); ++i) {
        const double y = yValues[i];
        const double density = densityValues[i];

        if (!isFinitePositive(y)) {
            throw std::runtime_error(
                "Cannot compute distribution moments with a non-positive y value.");
        }
        if (!std::isfinite(density) || density < 0.0) {
            throw std::runtime_error(
                "Cannot compute distribution moments with an invalid density value.");
        }

        totalProbability += normalizationFactor * y * density;
    }

    if (!(totalProbability > 0.0)) {
        throw std::runtime_error(
            "Cannot compute distribution moments with non-positive total probability.");
    }

    double meanNumerator = 0.0;
    for (std::size_t i = 0; i < yValues.size(); ++i) {
        const double probability = normalizationFactor * yValues[i] * densityValues[i];
        meanNumerator += probability * yValues[i];
    }
    moments.meanKeVPerUm = meanNumerator / totalProbability;

    double varianceNumerator = 0.0;
    for (std::size_t i = 0; i < yValues.size(); ++i) {
        const double probability = normalizationFactor * yValues[i] * densityValues[i];
        const double delta = yValues[i] - moments.meanKeVPerUm;
        varianceNumerator += probability * delta * delta;
    }

    moments.varianceKeV2PerUm2 = varianceNumerator / totalProbability;
    if (moments.varianceKeV2PerUm2 < 0.0 &&
        std::abs(moments.varianceKeV2PerUm2) < 1.0e-12) {
        moments.varianceKeV2PerUm2 = 0.0;
    }
    if (moments.varianceKeV2PerUm2 < 0.0) {
        throw std::runtime_error(
            "Cannot compute distribution moments with a negative variance.");
    }

    moments.stdevKeVPerUm = std::sqrt(moments.varianceKeV2PerUm2);

    if (moments.stdevKeVPerUm == 0.0) {
        moments.skewness = 0.0;
        return moments;
    }

    double skewnessNumerator = 0.0;
    for (std::size_t i = 0; i < yValues.size(); ++i) {
        const double probability = normalizationFactor * yValues[i] * densityValues[i];
        const double standardized =
            (yValues[i] - moments.meanKeVPerUm) / moments.stdevKeVPerUm;
        skewnessNumerator += probability * standardized * standardized * standardized;
    }

    moments.skewness = skewnessNumerator / totalProbability;
    return moments;
}

}  // namespace

SpectrumAccumulator::SpectrumAccumulator(
    const LookupLibrary& library,
    std::size_t rebinSamples,
    std::uint64_t rebinSeed)
    : rebinSamples_(rebinSamples), rebinSeed_(rebinSeed) {

    if (rebinSamples_ == 0) {
        throw std::runtime_error("Rebin sample count must be greater than zero.");
    }

    buildTargetYBins();
    numerator_.assign(yCenters_.size(), 0.0);

    const auto& tables = library.tables();
    precomputedRawFyByTable_.reserve(tables.size());

    for (std::size_t tableIndex = 0; tableIndex < tables.size(); ++tableIndex) {
        precomputedRawFyByTable_.push_back(
            buildPrecomputedRawFy(tables[tableIndex], tableIndex));
    }
}

void SpectrumAccumulator::buildTargetYBins() {
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
        yCenters_[i] = 0.5 * (targetEdges_[i] + targetEdges_[i + 1]);
        binWidths_[i] = targetEdges_[i + 1] - targetEdges_[i];
        if (binWidths_[i] <= 0.0) {
            throw std::runtime_error("Non-positive rebinned y-bin width encountered.");
        }
    }
}

// Build a common interval representation of the source LUT before rebinning.
//
// For Cartechini, the source data are tabulated f(y) values on a coarse native
// y grid. In that case, each interval mass is estimated by integrating f(y)
// over [y_i, y_{i+1}] using the trapezoidal rule.
//
// For DeCunha, the source data are raw histogram counts N(y). Each native bin
// already carries a fixed amount of spectral weight, so the interval mass is
// taken directly from the stored bin content rather than multiplying by bin
// width again.
//
// The returned SourceIntervalData stores:
//   - lower edge of each source interval
//   - upper edge of each source interval
//   - spectral weight ("mass") carried by that interval
//
// This interval representation is then consumed by either the stochastic
// rebinner or the deterministic overlap rebinner.
SpectrumAccumulator::SourceIntervalData
SpectrumAccumulator::buildSourceIntervals(const LookupTable& table) const {
    const auto& srcY = table.yLowerEdges();

    SourceIntervalData out;
    if (srcY.size() < 2) {
        return out;
    }

    if (table.spectrumKind() == LookupSpectrumKind::CartechiniFy) {
        const auto& srcF = table.fY();
        if (srcF.size() != srcY.size()) {
            throw std::runtime_error(
                "Cartechini LUT has mismatched y and f(y) sizes: y size = " +
                std::to_string(srcY.size()) + ", f(y) size = " +
                std::to_string(srcF.size()));
        }

        out.lower.reserve(srcY.size() - 1);
        out.upper.reserve(srcY.size() - 1);
        out.masses.reserve(srcY.size() - 1);

        for (std::size_t i = 0; i + 1 < srcY.size(); ++i) {
            const double y0 = srcY[i];
            const double y1 = srcY[i + 1];
            const double f0 = std::max(0.0, srcF[i]);
            const double f1 = std::max(0.0, srcF[i + 1]);

            if (!(y0 > 0.0) || !(y1 > y0)) {
                continue;
            }

            const double mass = 0.5 * (f0 + f1) * (y1 - y0);
            if (mass > 0.0) {
                out.lower.push_back(y0);
                out.upper.push_back(y1);
                out.masses.push_back(mass);
            }
        }

        return out;
    }

    const auto& srcN = table.nY();

    // DeCunha raw N(y) histogram:
    // Case 1: explicit final upper edge is present
    if (srcY.size() == srcN.size() + 1) {
        out.lower.reserve(srcN.size());
        out.upper.reserve(srcN.size());
        out.masses.reserve(srcN.size());

        for (std::size_t i = 0; i < srcN.size(); ++i) {
            const double y0 = srcY[i];
            const double y1 = srcY[i + 1];
            const double n = std::max(0.0, srcN[i]);

            if (!(y0 > 0.0) || !(y1 > y0) || !(n > 0.0)) {
                continue;
            }

            out.lower.push_back(y0);
            out.upper.push_back(y1);
            out.masses.push_back(n);
        }

        return out;
    }

    // DeCunha raw N(y) histogram:
    // Case 2: only lower edges are present; infer the last upper edge
    if (srcY.size() == srcN.size()) {
        out.lower.reserve(srcN.size());
        out.upper.reserve(srcN.size());
        out.masses.reserve(srcN.size());

        for (std::size_t i = 0; i < srcN.size(); ++i) {
            const double y0 = srcY[i];
            const double n = std::max(0.0, srcN[i]);

            double y1 = 0.0;
            if (i + 1 < srcY.size()) {
                y1 = srcY[i + 1];
            } else {
                // Infer final upper edge from the last two lower edges.
                // These LUTs are logarithmically binned, so extrapolate by ratio.
                const double yPrev = srcY[i - 1];
                if (!(yPrev > 0.0) || !(y0 > yPrev)) {
                    continue;
                }
                const double ratio = y0 / yPrev;
                y1 = y0 * ratio;
            }

            if (!(y0 > 0.0) || !(y1 > y0) || !(n > 0.0)) {
                continue;
            }

            out.lower.push_back(y0);
            out.upper.push_back(y1);
            out.masses.push_back(n);
        }

        return out;
    }

    throw std::runtime_error(
        "Unsupported lookup-table y-grid layout for source file " +
        table.sourceFile() + ": y size = " + std::to_string(srcY.size()) +
        ", nY size = " + std::to_string(srcN.size()) +
        ", fY size = " + std::to_string(table.fY().size()));
}

// Rebin a source spectrum onto the common target y grid by stochastic sampling.
//
// The input intervals define a piecewise distribution over y. We first build a
// CDF from the interval spectral weights, then repeatedly:
//
//   1. sample a source interval according to its weight
//   2. sample one y value inside that interval
//   3. place the sampled y into the corresponding target bin
//
// Sampling inside the chosen interval is uniform in log10(y), which is natural
// for microdosimetric spectra spanning multiple orders of magnitude and helps
// avoid artificial empty bins when the native source grid is coarse.
//
// This method is primarily used for Cartechini f(y) LUTs, whose native spacing
// is sparse enough that direct point-to-bin remapping can produce comb-like
// zero bins on the target grid.
std::vector<double> SpectrumAccumulator::sampleCountsToTargetGrid(
    const SourceIntervalData& intervals,
    std::size_t nSamples,
    std::uint64_t seed) const {

    std::vector<double> rebinnedCounts(yCenters_.size(), 0.0);

    if (intervals.lower.empty() || intervals.upper.empty() ||
        intervals.masses.empty() || nSamples == 0) {
        return rebinnedCounts;
    }

    const double totalMass = std::accumulate(
        intervals.masses.begin(), intervals.masses.end(), 0.0);

    if (!(totalMass > 0.0)) {
        return rebinnedCounts;
    }

    std::vector<double> cdf(intervals.masses.size(), 0.0);
    double running = 0.0;
    for (std::size_t i = 0; i < intervals.masses.size(); ++i) {
        running += intervals.masses[i] / totalMass;
        cdf[i] = running;
    }
    cdf.back() = 1.0;

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> unit01(0.0, 1.0);

    for (std::size_t sample = 0; sample < nSamples; ++sample) {
        const double uInterval = unit01(rng);
        auto itInterval = std::lower_bound(cdf.begin(), cdf.end(), uInterval);
        std::size_t idx = static_cast<std::size_t>(
            std::distance(cdf.begin(), itInterval));

        if (idx >= intervals.lower.size()) {
            idx = intervals.lower.size() - 1;
        }

        const double y0 = intervals.lower[idx];
        const double y1 = intervals.upper[idx];

        if (!(y0 > 0.0) || !(y1 > y0)) {
            continue;
        }

        const double logY0 = std::log10(y0);
        const double logY1 = std::log10(y1);
        const double uWithin = unit01(rng);
        const double sampledY = std::pow(10.0, logY0 + uWithin * (logY1 - logY0));

        if (!(sampledY >= targetEdges_.front() && sampledY < targetEdges_.back())) {
            continue;
        }

        auto itTarget = std::upper_bound(
            targetEdges_.begin(), targetEdges_.end(), sampledY);

        if (itTarget == targetEdges_.begin()) {
            continue;
        }

        if (itTarget == targetEdges_.end()) {
            rebinnedCounts.back() += 1.0;
            continue;
        }

        const std::size_t targetIndex =
            static_cast<std::size_t>(std::distance(targetEdges_.begin(), itTarget) - 1);

        if (targetIndex < rebinnedCounts.size()) {
            rebinnedCounts[targetIndex] += 1.0;
        }
    }

    return rebinnedCounts;
}

// Wrapper for the stochastic rebinner.
//
// This converts the LUT into interval form and then performs Monte Carlo
// transport of spectral weight onto the target grid. It is kept as a separate
// function so the family-specific rebinning policy can be expressed clearly
// in buildPrecomputedRawFy().
std::vector<double> SpectrumAccumulator::rebinCountsToTargetGridStochastic(
    const LookupTable& table,
    std::size_t nSamples,
    std::uint64_t seed) const {

    const SourceIntervalData intervals = buildSourceIntervals(table);
    return sampleCountsToTargetGrid(intervals, nSamples, seed);
}

// Rebin a source spectrum onto the common target y grid by deterministic
// log-bin overlap.
//
// Each source interval is treated as carrying a fixed amount of spectral weight.
// That weight is redistributed across target bins according to the fraction of
// overlap in log10(y) space between the source interval and each target bin.
//
// No random sampling is performed, so this method is noise-free and exactly
// reproducible. It is the preferred method for DeCunha LUTs because those
// files store native histogram bin contents N(y), meaning each source bin
// already has a well-defined amount of spectral weight that should be
// transported directly to the target grid.
std::vector<double> SpectrumAccumulator::rebinCountsToTargetGridDeterministic(
    const LookupTable& table) const {

    const SourceIntervalData intervals = buildSourceIntervals(table);
    std::vector<double> rebinnedCounts(yCenters_.size(), 0.0);

    if (intervals.lower.empty() || intervals.upper.empty() || intervals.masses.empty()) {
        return rebinnedCounts;
    }

    std::vector<double> logTargetEdges(targetEdges_.size(), 0.0);
    for (std::size_t i = 0; i < targetEdges_.size(); ++i) {
        logTargetEdges[i] = std::log10(targetEdges_[i]);
    }

    for (std::size_t i = 0; i < intervals.masses.size(); ++i) {
        const double y0 = intervals.lower[i];
        const double y1 = intervals.upper[i];
        const double mass = intervals.masses[i];

        if (!(y0 > 0.0) || !(y1 > y0) || !(mass > 0.0)) {
            continue;
        }

        const double lo = std::max(y0, targetEdges_.front());
        const double hi = std::min(y1, targetEdges_.back());

        if (!(hi > lo)) {
            continue;
        }

        const double logY0 = std::log10(lo);
        const double logY1 = std::log10(hi);
        const double logWidth = logY1 - logY0;

        if (!(logWidth > 0.0)) {
            continue;
        }

        for (std::size_t j = 0; j < yCenters_.size(); ++j) {
            const double a = std::max(logY0, logTargetEdges[j]);
            const double b = std::min(logY1, logTargetEdges[j + 1]);

            if (b > a) {
                const double frac = (b - a) / logWidth;
                rebinnedCounts[j] += mass * frac;
            }
        }
    }

    return rebinnedCounts;
}

// Precompute the monoenergetic contribution on the shared target y grid.
//
// This is where the family-specific rebinning policy is applied:
//
//   - Cartechini:
//       The source LUT stores f(y) directly on a coarse native grid.
//       We rebin it stochastically to avoid artificial zero bins, then
//       normalize the rebinned result back to a density on the target grid.
//
//   - DeCunha:
//       The source LUT stores raw histogram counts N(y) together with Ncpp.
//       We rebin by deterministic log-bin overlap so native bin content is
//       preserved without Monte Carlo noise, then convert the rebinned counts
//       into the raw f(y)-like quantity used later in polyenergetic summation.
//
// The output of this function is stored once per LUT and reused during
// accumulation to avoid repeating the rebinning work for every proton.
std::vector<double> SpectrumAccumulator::buildPrecomputedRawFy(
    const LookupTable& table,
    std::size_t tableIndex) const {

    const std::uint64_t energyKey =
        static_cast<std::uint64_t>(std::llround(table.monoEnergyMeV() * 1.0e6));

    const std::uint64_t tableSeed =
        rebinSeed_
        ^ (energyKey + 0x9E3779B97F4A7C15ULL)
        ^ (static_cast<std::uint64_t>(tableIndex) * 0xBF58476D1CE4E5B9ULL);

    std::vector<double> rebinnedCounts;

    if (table.spectrumKind() == LookupSpectrumKind::CartechiniFy) {
        rebinnedCounts = rebinCountsToTargetGridStochastic(table, rebinSamples_, tableSeed);

        double totalRebinnedCounts =
            std::accumulate(rebinnedCounts.begin(), rebinnedCounts.end(), 0.0);

        if (!(totalRebinnedCounts > 0.0)) {
            rebinnedCounts = rebinCountsToTargetGridDeterministic(table);
            totalRebinnedCounts =
                std::accumulate(rebinnedCounts.begin(), rebinnedCounts.end(), 0.0);
        }

        std::vector<double> rawFy(yCenters_.size(), 0.0);
        if (!(totalRebinnedCounts > 0.0)) {
            return rawFy;
        }

        for (std::size_t i = 0; i < rebinnedCounts.size(); ++i) {
            rawFy[i] = (rebinnedCounts[i] / totalRebinnedCounts) / binWidths_[i];
        }
        return rawFy;
    }

    // For DeCunha raw N(y) histograms, use deterministic log-bin overlap only.
    rebinnedCounts = rebinCountsToTargetGridDeterministic(table);

    const double totalRebinnedCounts =
        std::accumulate(rebinnedCounts.begin(), rebinnedCounts.end(), 0.0);

    std::vector<double> rawFy(yCenters_.size(), 0.0);
    if (!(totalRebinnedCounts > 0.0)) {
        return rawFy;
    }

    const double ncpp = table.ncpp();
    if (!std::isfinite(ncpp) || ncpp <= 0.0) {
        throw std::runtime_error(
            "Invalid Ncpp for DeCunha table: " + table.sourceFile());
    }

    for (std::size_t i = 0; i < rebinnedCounts.size(); ++i) {
        const double fMono = (rebinnedCounts[i] / totalRebinnedCounts) / binWidths_[i];
        rawFy[i] = fMono * ncpp;
    }

    return rawFy;
}

void SpectrumAccumulator::addContributionByIndex(
    std::size_t tableIndex,
    double multiplicity) {

    addInterpolatedContribution(tableIndex, tableIndex, 1.0, 0.0, multiplicity);
}

// Interpolate between the two precomputed monoenergetic spectra selected for a
// proton, then scale that interpolated spectrum by the proton multiplicity.
void SpectrumAccumulator::addInterpolatedContribution(
    std::size_t lowerTableIndex,
    std::size_t upperTableIndex,
    double lowerWeight,
    double upperWeight,
    double multiplicity) {

    if (multiplicity <= 0.0) {
        return;
    }

    if (lowerTableIndex >= precomputedRawFyByTable_.size() ||
        upperTableIndex >= precomputedRawFyByTable_.size()) {
        throw std::runtime_error(
            "SpectrumAccumulator received invalid interpolation table index.");
    }

    if (lowerWeight < 0.0 || upperWeight < 0.0) {
        throw std::runtime_error(
            "SpectrumAccumulator received negative interpolation weights.");
    }

    const double weightSum = lowerWeight + upperWeight;
    if (!(weightSum > 0.0)) {
        return;
    }

    const auto& lowerRawFy = precomputedRawFyByTable_[lowerTableIndex];
    const auto& upperRawFy = precomputedRawFyByTable_[upperTableIndex];

    for (std::size_t i = 0; i < lowerRawFy.size(); ++i) {
        const double interpolatedRawFy =
            lowerWeight * lowerRawFy[i] + upperWeight * upperRawFy[i];
        numerator_[i] += multiplicity * interpolatedRawFy;
    }
    denominator_ += multiplicity;
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

    std::vector<double> rawFy(yCenters_.size(), 0.0);
    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        rawFy[i] = numerator_[i] / denominator_;
    }

    const double C = logarithmicBinNormalizationFactor();

    double fyNormDen = 0.0;
    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        fyNormDen += out.y[i] * rawFy[i];
    }
    fyNormDen *= C;

    if (!(fyNormDen > 0.0)) {
        throw std::runtime_error("Polyenergetic f(y) has non-positive normalization.");
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

    if (!(dNormDen > 0.0)) {
        throw std::runtime_error("Polyenergetic d(y) has non-positive normalization.");
    }

    for (std::size_t i = 0; i < yCenters_.size(); ++i) {
        out.d[i] = out.yf[i] / dNormDen;
        out.yd[i] = out.y[i] * out.d[i];
    }

    return out;
}

PolySpectrumMomentsSummary SpectrumAccumulator::summarizeMoments(
    const PolySpectrum& spectrum) const {

    validateSpectrumShape(spectrum);

    if (spectrum.y.size() != yCenters_.size()) {
        throw std::runtime_error(
            "Cannot summarize moments for a spectrum built on a different y grid.");
    }

    const double normalizationFactor = logarithmicBinNormalizationFactor();

    PolySpectrumMomentsSummary summary;
    summary.frequency = computeDistributionMoments(
        spectrum.y, spectrum.f, "frequency", normalizationFactor);
    summary.dose = computeDistributionMoments(
        spectrum.y, spectrum.d, "dose", normalizationFactor);
    return summary;
}

double SpectrumAccumulator::logarithmicBinNormalizationFactor() const {
    if (yCenters_.size() < 2) {
        throw std::runtime_error("Insufficient rebinned y-grid size.");
    }

    const double deltaLogY = std::log10(yCenters_[1]) - std::log10(yCenters_[0]);
    const double normalizationFactor = std::log(10.0) * deltaLogY;

    if (!(normalizationFactor > 0.0)) {
        throw std::runtime_error("Invalid logarithmic-bin normalization factor C.");
    }

    return normalizationFactor;
}
