#include "PhaseSpaceReader.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct ParticleIdentity {
    bool isCharged{false};
    int atomicNumber{0};
    int massNumber{0};
};

ParticleIdentity decodeParticleIdentity(int pdgCode) {
    const int absPdg = std::abs(pdgCode);

    // Charged leptons.
    if (absPdg == 11 || absPdg == 13 || absPdg == 15) {
        return ParticleIdentity{true, 0, 0};
    }

    // Common charged hadrons found in phase-space files.
    if (absPdg == 2212) {
        return ParticleIdentity{true, 1, 1};
    }

    if (absPdg == 211 || absPdg == 321) {
        return ParticleIdentity{true, 0, 0};
    }

    // Ions use the PDG nuclear code 10LZZZAAAI. A nonzero Z is charged.
    if (absPdg >= 1000000000) {
        const int z = (absPdg / 10000) % 1000;
        const int a = (absPdg / 10) % 1000;
        if (z > 0 && a > 0) {
            return ParticleIdentity{true, z, a};
        }
    }

    return ParticleIdentity{};
}

}  // namespace

bool PhaseSpaceReader::isChargedParticlePdg(int pdgCode) {
    return decodeParticleIdentity(pdgCode).isCharged;
}

int PhaseSpaceReader::atomicNumberFromPdg(int pdgCode) {
    return decodeParticleIdentity(pdgCode).atomicNumber;
}

int PhaseSpaceReader::massNumberFromPdg(int pdgCode) {
    return decodeParticleIdentity(pdgCode).massNumber;
}

std::vector<ChargedParticleRecord> PhaseSpaceReader::readChargedParticles(
    const std::filesystem::path& phspPath) const {

    std::ifstream in(phspPath);
    if (!in) {
        throw std::runtime_error("Failed to open phase space file: " + phspPath.string());
    }

    std::vector<ChargedParticleRecord> chargedParticles;
    std::string line;
    std::size_t rowIndex = 0;

    while (std::getline(in, line)) {
        ++rowIndex;

        if (line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::vector<std::string> fields;
        std::string token;

        while (iss >> token) {
            fields.push_back(token);
        }

        // Need at least columns 1..8
        if (fields.size() < 8) {
            continue;
        }

        try {
            const double energyMeV = std::stod(fields[5]);   // column 6
            const double weight    = std::stod(fields[6]);   // column 7
            const int pdgCode      = static_cast<int>(std::stod(fields[7])); // column 8
            const ParticleIdentity identity = decodeParticleIdentity(pdgCode);

            if (identity.isCharged) {
                chargedParticles.push_back(ChargedParticleRecord{
                    rowIndex,
                    energyMeV,
                    weight,
                    pdgCode,
                    identity.atomicNumber,
                    identity.massNumber
                });
            }
        } catch (...) {
            // skip malformed rows
            continue;
        }
    }

    return chargedParticles;
}

std::vector<ProtonRecord> PhaseSpaceReader::readProtons(
    const std::filesystem::path& phspPath) const {

    const std::vector<ChargedParticleRecord> chargedParticles =
        readChargedParticles(phspPath);

    std::vector<ProtonRecord> protons;
    protons.reserve(chargedParticles.size());

    for (const auto& particle : chargedParticles) {
        if (particle.pdgCode == 2212) {
            protons.push_back(ProtonRecord{
                particle.rowIndex,
                particle.energyMeV,
                particle.weight,
                particle.pdgCode,
                particle.atomicNumber,
                particle.massNumber
            });
        }
    }

    return protons;
}
