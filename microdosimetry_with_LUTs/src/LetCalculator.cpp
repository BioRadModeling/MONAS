#include "LetCalculator.h"

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "LetLookup.h"

namespace {

constexpr int protonPdgCode = 2212;
constexpr int protonIonPdgCode = 1000010010;

std::filesystem::path letTablePath(const std::filesystem::path& letDirectory,
                                   const std::string& element) {
    return letDirectory / (element + "_water.txt");
}

bool isProtonPdg(int pdgCode) {
    return pdgCode == protonPdgCode || pdgCode == protonIonPdgCode;
}

}  // namespace

LetCalculationResult LetCalculator::calculate(
    const std::filesystem::path& letDirectory,
    const std::vector<ChargedParticleRecord>& chargedParticles) const {

    static const std::array<std::string, 8> elements{
        "H", "He", "Li", "Be", "B", "C", "N", "O"
    };

    std::unordered_map<std::string, LetLookup> lookups;
    lookups.reserve(elements.size());
    for (const auto& element : elements) {
        lookups.emplace(element,
                        LetLookup::loadFromFile(letTablePath(letDirectory, element),
                                                element));
    }

    LetAccumulator accumulator;
    std::size_t matchedParticleCount = 0;

    for (const auto& particle : chargedParticles) {
        const std::optional<std::string> element = lutElementForPdg(particle.pdgCode);
        if (!element.has_value()) {
            continue;
        }

        const auto lookupIt = lookups.find(*element);
        if (lookupIt == lookups.end()) {
            continue;
        }

        const double letKeVPerUm =
            lookupIt->second.interpolate(particle.energyMeV);
        addParticleToAccumulator(accumulator, particle, letKeVPerUm);
        ++matchedParticleCount;
    }

    return LetCalculationResult{accumulator.summarize(), matchedParticleCount};
}

std::optional<std::string> LetCalculator::lutElementForPdg(int pdgCode) {
    switch (pdgCode) {
        case protonPdgCode:
        case protonIonPdgCode:
        case 1000010020:
        case 1000010030:
            return std::string("H");
        case 1000020030:
        case 1000020040:
            return std::string("He");
        case 1000030060:
        case 1000030070:
            return std::string("Li");
        case 1000040070:
        case 1000040090:
            return std::string("Be");
        case 1000050100:
        case 1000050110:
            return std::string("B");
        case 1000060110:
        case 1000060120:
        case 1000060130:
            return std::string("C");
        case 1000070130:
        case 1000070140:
        case 1000070150:
            return std::string("N");
        case 1000080150:
        case 1000080160:
            return std::string("O");
        default:
            return std::nullopt;
    }
}

void LetCalculator::addParticleToAccumulator(
    LetAccumulator& accumulator,
    const ChargedParticleRecord& particle,
    double letKeVPerUm) {

    if (isProtonPdg(particle.pdgCode)) {
        accumulator.add(LetGroup::Protons, particle.weight, letKeVPerUm);
    }

    accumulator.add(LetGroup::AllCharged, particle.weight, letKeVPerUm);
}
