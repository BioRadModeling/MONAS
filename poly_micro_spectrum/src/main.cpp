#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CsvWriter.h"
#include "LookupLibrary.h"
#include "PhaseSpaceReader.h"
#include "SpectrumAccumulator.h"

namespace fs = std::filesystem;

static std::string resolveFolderName(const std::string& voxelSize,
                                     const std::string& energyGrid) {
    if (voxelSize == "1mm" && energyGrid == "linear") {
        return "1mm_linear";
    }
    if (voxelSize == "1mm" && energyGrid == "log") {
        return "1mm_logarithmic";
    }
    if (voxelSize == "5um" && energyGrid == "linear") {
        return "5um_linear";
    }
    if (voxelSize == "5um" && energyGrid == "log") {
        return "5um_logarithmic";
    }

    throw std::runtime_error("Invalid voxelSize / energyGrid combination.");
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            std::cerr
                << "Usage:\n"
                << "  " << argv[0]
                << " test-lookup <lookup_root> <1mm|5um> <linear|log> <test_energy_mev>\n"
                << "  " << argv[0]
                << " audit-phsp <lookup_root> <1mm|5um> <linear|log> <phsp_file> <output_dir>\n"
                << "  " << argv[0]
                << " build-spectrum <lookup_root> <1mm|5um> <linear|log> <phsp_file> <output_dir>\n";
            return 1;
        }

        const std::string mode = argv[1];

        if (mode == "test-lookup") {
            if (argc != 6) {
                throw std::runtime_error(
                    "test-lookup requires: <lookup_root> <1mm|5um> <linear|log> <test_energy_mev>");
            }

            const fs::path lookupRoot = argv[2];
            const std::string voxelSize = argv[3];
            const std::string energyGrid = argv[4];
            const double testEnergyMeV = std::stod(argv[5]);

            const std::string folderName = resolveFolderName(voxelSize, energyGrid);
            const fs::path libraryDir = lookupRoot / "csv" / folderName;

            LookupLibrary library;
            library.loadFromDirectory(libraryDir);

            const LookupTable& match = library.findNearest(testEnergyMeV);

            std::cout << "Loaded library successfully.\n";
            std::cout << "Library folder:     " << libraryDir << "\n";
            std::cout << "Number of tables:   " << library.size() << "\n";
            std::cout << "Y bins:             " << library.yReference().size() << "\n";
            std::cout << "Requested energy:   " << testEnergyMeV << " MeV\n";
            std::cout << "Matched energy:     " << match.monoEnergyMeV() << " MeV\n";
            std::cout << "Matched CSV:        " << match.sourceFile() << "\n";
            std::cout << "Detected Ncpp:      " << match.ncpp() << "\n";

            return 0;
        }

        if (mode == "audit-phsp" || mode == "build-spectrum") {
            if (argc != 7) {
                throw std::runtime_error(
                    std::string(mode) +
                    " requires: <lookup_root> <1mm|5um> <linear|log> <phsp_file> <output_dir>");
            }

            const fs::path lookupRoot = argv[2];
            const std::string voxelSize = argv[3];
            const std::string energyGrid = argv[4];
            const fs::path phspFile = argv[5];
            const fs::path outputDir = argv[6];

            const std::string folderName = resolveFolderName(voxelSize, energyGrid);
            const fs::path libraryDir = lookupRoot / "csv" / folderName;

            LookupLibrary library;
            library.loadFromDirectory(libraryDir);

            PhaseSpaceReader reader;
            const std::vector<ProtonRecord> protons = reader.readProtons(phspFile);

            fs::create_directories(outputDir);

            std::vector<ProtonMatchRecord> matches;
            matches.reserve(protons.size());

            SpectrumAccumulator accumulator(library.yReference());

            for (const auto& proton : protons) {
                const LookupTable& match = library.findNearest(proton.energyMeV);

                matches.push_back(ProtonMatchRecord{
                    proton.rowIndex,
                    proton.energyMeV,
                    proton.weight,
                    match.monoEnergyMeV(),
                    match.sourceFile(),
                    match.ncpp()
                });

                if (mode == "build-spectrum") {
                    accumulator.addContribution(match, proton.weight);
                }
            }

            const fs::path matchCsv = outputDir / "proton_matches.csv";
            CsvWriter::writeProtonMatches(matchCsv, matches);

            std::cout << "Phase-space processing completed.\n";
            std::cout << "Library folder:     " << libraryDir << "\n";
            std::cout << "Phase-space file:   " << phspFile << "\n";
            std::cout << "Protons found:      " << protons.size() << "\n";
            std::cout << "Match CSV:          " << matchCsv << "\n";

            if (mode == "build-spectrum") {
                const PolySpectrum spectrum = accumulator.finalize();
                const fs::path polyCsv = outputDir / "poly_spectrum.csv";
                CsvWriter::writePolySpectrum(polyCsv, spectrum);

                std::cout << "Poly spectrum CSV:  " << polyCsv << "\n";
            }

            return 0;
        }

        throw std::runtime_error("Unknown mode: " + mode);

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}