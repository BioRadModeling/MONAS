#pragma once

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <vector>

struct ProtonRecord {
    std::size_t rowIndex;
    double energyMeV;
    double weight;
    int pdgCode;
    int atomicNumber;
    int massNumber;
};

struct ChargedParticleRecord {
    std::size_t rowIndex;
    double energyMeV;
    double weight;
    int pdgCode;
    int atomicNumber;
    int massNumber;
};

template <typename ParticleRecord>
inline bool hasInaniwaLutIdentity(const ParticleRecord& particle) {
    return particle.atomicNumber > 0 && particle.massNumber > 0;
}

template <typename ParticleRecord>
inline double inaniwaEventEnergyMeV(const ParticleRecord& particle) {
    return particle.energyMeV;
}

template <typename ParticleRecord>
inline double inaniwaEnergyMeVPerU(const ParticleRecord& particle) {
    if (!hasInaniwaLutIdentity(particle)) {
        throw std::invalid_argument(
            "Inaniwa LUT energy conversion requires a particle with valid "
            "atomic and mass numbers.");
    }

    return particle.energyMeV / static_cast<double>(particle.massNumber);
}

class PhaseSpaceReader {
public:
    static bool isChargedParticlePdg(int pdgCode);
    static int atomicNumberFromPdg(int pdgCode);
    static int massNumberFromPdg(int pdgCode);

    std::vector<ProtonRecord> readProtons(const std::filesystem::path& phspPath) const;
    std::vector<ChargedParticleRecord> readChargedParticles(
        const std::filesystem::path& phspPath) const;
};
