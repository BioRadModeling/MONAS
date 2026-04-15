#include "PhaseSpaceReader.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<ProtonRecord> PhaseSpaceReader::readProtons(
    const std::filesystem::path& phspPath) const {

    std::ifstream in(phspPath);
    if (!in) {
        throw std::runtime_error("Failed to open phase space file: " + phspPath.string());
    }

    std::vector<ProtonRecord> protons;
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

            if (pdgCode == 2212) {
                protons.push_back(ProtonRecord{
                    rowIndex,
                    energyMeV,
                    weight,
                    pdgCode
                });
            }
        } catch (...) {
            // skip malformed rows
            continue;
        }
    }

    return protons;
}