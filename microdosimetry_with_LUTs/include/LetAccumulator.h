#pragma once

#include <string>
#include <vector>

enum class LetGroup {
    AllCharged,
    Protons,
    OtherCharged
};

struct LetSummaryRecord {
    std::string group;
    double trackAveragedLetKeVPerUm;
    double doseAveragedLetKeVPerUm;
};

class LetAccumulator {
public:
    void add(LetGroup group, double weight, double letKeVPerUm);
    std::vector<LetSummaryRecord> summarize() const;

private:
    struct Bucket {
        double trackNumerator{0.0};
        double trackDenominator{0.0};
        double doseNumerator{0.0};
        double doseDenominator{0.0};
    };

    Bucket allCharged_;
    Bucket protons_;
    Bucket otherCharged_;

    static LetSummaryRecord summarizeBucket(const std::string& name,
                                            const Bucket& bucket);
};
