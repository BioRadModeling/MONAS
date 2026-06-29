#include "SpectrumMoments.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace {

bool isFinitePositive(double x) {
    return std::isfinite(x) && x > 0.0;
}

}  // namespace

void RatioMeanStandardErrorAccumulator::addSample(
    double numerator,
    double denominator) {

    if (!std::isfinite(numerator) || !std::isfinite(denominator) ||
        denominator <= 0.0) {
        return;
    }

    ++count_;
    numeratorSum_ += numerator;
    denominatorSum_ += denominator;
    numeratorSquaredSum_ += numerator * numerator;
    denominatorSquaredSum_ += denominator * denominator;
    numeratorDenominatorSum_ += numerator * denominator;
}

double RatioMeanStandardErrorAccumulator::mean() const {
    if (!(denominatorSum_ > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return numeratorSum_ / denominatorSum_;
}

double RatioMeanStandardErrorAccumulator::standardError() const {
    if (count_ < 2 || !(denominatorSum_ > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double ratio = mean();
    if (!std::isfinite(ratio)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double residualSum =
        numeratorSquaredSum_
        - 2.0 * ratio * numeratorDenominatorSum_
        + ratio * ratio * denominatorSquaredSum_;

    if (residualSum < 0.0 && std::abs(residualSum) < 1.0e-12) {
        residualSum = 0.0;
    }
    if (residualSum < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double variance =
        (static_cast<double>(count_) / static_cast<double>(count_ - 1))
        * residualSum / (denominatorSum_ * denominatorSum_);

    if (variance < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::sqrt(variance);
}

double logarithmicBinNormalizationFactorForCenters(
    const std::vector<double>& yValues) {

    if (yValues.size() < 2) {
        throw std::runtime_error("Insufficient y-grid size.");
    }

    const double deltaLogY = std::log10(yValues[1]) - std::log10(yValues[0]);
    const double normalizationFactor = std::log(10.0) * deltaLogY;

    if (!(normalizationFactor > 0.0)) {
        throw std::runtime_error("Invalid logarithmic-bin normalization factor C.");
    }

    return normalizationFactor;
}

SpectrumDistributionMoments computeSpectrumDistributionMoments(
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
