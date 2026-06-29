#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

struct SpectrumDistributionMoments {
    std::string distribution;
    double meanKeVPerUm{0.0};
    double varianceKeV2PerUm2{0.0};
    double stdevKeVPerUm{0.0};
    double meanStandardErrorKeVPerUm{
        std::numeric_limits<double>::quiet_NaN()};
    double skewness{0.0};
};

class RatioMeanStandardErrorAccumulator {
public:
    void addSample(double numerator, double denominator);

    std::size_t count() const { return count_; }
    double numeratorSum() const { return numeratorSum_; }
    double denominatorSum() const { return denominatorSum_; }
    double mean() const;
    double standardError() const;

private:
    std::size_t count_{0};
    double numeratorSum_{0.0};
    double denominatorSum_{0.0};
    double numeratorSquaredSum_{0.0};
    double denominatorSquaredSum_{0.0};
    double numeratorDenominatorSum_{0.0};
};

double logarithmicBinNormalizationFactorForCenters(
    const std::vector<double>& yValues);

SpectrumDistributionMoments computeSpectrumDistributionMoments(
    const std::vector<double>& yValues,
    const std::vector<double>& densityValues,
    const std::string& distribution,
    double normalizationFactor);
