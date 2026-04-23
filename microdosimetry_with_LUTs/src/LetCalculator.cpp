#include "LetCalculator.h"

#include <filesystem>
#include <string>
#include <vector>

#include "LetLookup.h"

namespace {

constexpr int protonPdgCode = 2212;

std::filesystem::path letTablePath(const std::filesystem::path& letDirectory,
                                   const std::string& element) {
    return letDirectory / (element + "_water.txt");
}

}  // namespace

const std::vector<std::string>& LetCalculator::supportedElements() {
    static const std::vector<std::string> elements{
        "H", "He", "Li", "Be", "B", "C", "N", "O"
    };
    return elements;
}

std::vector<ElementLetSummary> LetCalculator::calculateByElement(
    const std::filesystem::path& letDirectory,
    const std::vector<ChargedParticleRecord>& chargedParticles) const {

    std::vector<ElementLetSummary> summaries;
    summaries.reserve(supportedElements().size());

    for (const std::string& element : supportedElements()) {
        const LetLookup lookup =
            LetLookup::loadFromFile(letTablePath(letDirectory, element), element);
        LetAccumulator accumulator;

        for (const auto& particle : chargedParticles) {
            const double letKeVPerUm = lookup.interpolate(particle.energyMeV);
            addParticleToAccumulator(accumulator, particle, letKeVPerUm);
        }

        summaries.push_back(ElementLetSummary{
            element,
            accumulator.summarize()
        });
    }

    return summaries;
}

void LetCalculator::addParticleToAccumulator(
    LetAccumulator& accumulator,
    const ChargedParticleRecord& particle,
    double letKeVPerUm) {

    accumulator.add(LetGroup::AllCharged, particle.weight, letKeVPerUm);

    if (particle.pdgCode == protonPdgCode) {
        accumulator.add(LetGroup::Protons, particle.weight, letKeVPerUm);
    } else {
        accumulator.add(LetGroup::OtherCharged, particle.weight, letKeVPerUm);
    }
}
