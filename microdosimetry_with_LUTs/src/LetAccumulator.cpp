#include "LetAccumulator.h"

#include <limits>
#include <stdexcept>

void LetAccumulator::add(LetGroup group,
                         double weight,
                         double letKeVPerUm) {
    Bucket* bucket = nullptr;

    switch (group) {
        case LetGroup::Protons:
            bucket = &protons_;
            break;
        case LetGroup::AllCharged:
            bucket = &allCharged_;
            break;
    }

    if (bucket == nullptr) {
        throw std::runtime_error("Unknown LET accumulator group.");
    }

    bucket->trackNumerator += weight * letKeVPerUm;
    bucket->trackDenominator += weight;
    bucket->doseNumerator += weight * letKeVPerUm * letKeVPerUm;
    bucket->doseDenominator += weight * letKeVPerUm;
}

std::vector<LetSummaryRecord> LetAccumulator::summarize() const {
    return {
        summarizeBucket("protons", protons_),
        summarizeBucket("all_charged", allCharged_)
    };
}

LetSummaryRecord LetAccumulator::summarizeBucket(const std::string& name,
                                                 const Bucket& bucket) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double trackAverage =
        bucket.trackDenominator > 0.0 ? bucket.trackNumerator /
                                            bucket.trackDenominator
                                      : nan;
    const double doseAverage =
        bucket.doseDenominator > 0.0 ? bucket.doseNumerator /
                                           bucket.doseDenominator
                                     : nan;

    return LetSummaryRecord{name, trackAverage, doseAverage};
}
