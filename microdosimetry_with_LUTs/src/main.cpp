#include <cstddef>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CsvWriter.h"
#include "LookupLibrary.h"
#include "LookupTable.h"
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

    throw std::runtime_error(
        "Invalid voxelSize / energyGrid combination. "
        "Use voxelSize = 1mm or 5um, and energyGrid = linear or log.");
}

static fs::path resolveLibraryDir(const fs::path& lookupRoot,
                                  const std::string& lutName,
                                  LutFamily family,
                                  const std::string& voxelSize = "",
                                  const std::string& energyGrid = "") {
    if (family == LutFamily::Cartechini) {
        return lookupRoot / lutName;
    }

    const std::string folderName = resolveFolderName(voxelSize, energyGrid);
    return lookupRoot / lutName / folderName;
}

static void printUsage(const char* programName) {
    std::cerr
        << "Usage:\n"
        << "  DeCunha:\n"
        << "    " << programName
        << " test-lookup <lookupRoot> DeCunha <1mm|5um> <linear|log> <energyMeV>\n"
        << "    " << programName
        << " audit-phsp <lookupRoot> DeCunha <1mm|5um> <linear|log> <phspFile> <outputDir>\n"
        << "    " << programName
        << " build-spectrum <lookupRoot> DeCunha <1mm|5um> <linear|log> <phspFile> <outputDir>\n"
        << "\n"
        << "  Cartechini:\n"
        << "    " << programName
        << " test-lookup <lookupRoot> Cartechini <energyMeV>\n"
        << "    " << programName
        << " audit-phsp <lookupRoot> Cartechini <phspFile> <outputDir>\n"
        << "    " << programName
        << " build-spectrum <lookupRoot> Cartechini <phspFile> <outputDir>\n";
}

static std::string lutFamilyToString(LutFamily family) {
    switch (family) {
        case LutFamily::DeCunha:
            return "DeCunha";
        case LutFamily::Cartechini:
            return "Cartechini";
        default:
            return "Unknown";
    }
}

static std::string spectrumKindToString(LookupSpectrumKind kind) {
    switch (kind) {
        case LookupSpectrumKind::DeCunhaRawCounts:
            return "DeCunhaRawCounts";
        case LookupSpectrumKind::CartechiniFy:
            return "CartechiniFy";
        default:
            return "Unknown";
    }
}

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            printUsage(argv[0]);
            return 1;
        }

        const std::string mode = argv[1];

        if (mode == "test-lookup") {
            if (argc < 5) {
                printUsage(argv[0]);
                throw std::runtime_error("test-lookup requires more arguments.");
            }

            const fs::path lookupRoot = argv[2];
            const std::string lutName = argv[3];
            const LutFamily family = LookupLibrary::inferFamily(lutName);

            fs::path libraryDir;
            double testEnergyMeV = 0.0;

            if (family == LutFamily::DeCunha) {
                if (argc != 7) {
                    printUsage(argv[0]);
                    throw std::runtime_error(
                        "test-lookup for DeCunha requires: "
                        "<lookupRoot> DeCunha <1mm|5um> <linear|log> <energyMeV>");
                }

                const std::string voxelSize = argv[4];
                const std::string energyGrid = argv[5];
                testEnergyMeV = std::stod(argv[6]);

                libraryDir = resolveLibraryDir(
                    lookupRoot, lutName, family, voxelSize, energyGrid);
            } else {
                if (argc != 5) {
                    printUsage(argv[0]);
                    throw std::runtime_error(
                        "test-lookup for Cartechini requires: "
                        "<lookupRoot> Cartechini <energyMeV>");
                }

                testEnergyMeV = std::stod(argv[4]);
                libraryDir = resolveLibraryDir(lookupRoot, lutName, family);
            }

            if (!fs::exists(libraryDir) || !fs::is_directory(libraryDir)) {
                throw std::runtime_error(
                    "Lookup library folder not found: " + libraryDir.string());
            }

            LookupLibrary library;
            library.loadFromDirectory(libraryDir, family);

            const LookupTable& match = library.findNearest(testEnergyMeV);

            std::cout << "Loaded library successfully.\n";
            std::cout << "LUT family: " << lutFamilyToString(family) << "\n";
            std::cout << "Library folder: " << libraryDir << "\n";
            std::cout << "Number of tables: " << library.size() << "\n";
            std::cout << "Y bins: " << library.yReference().size() << "\n";
            std::cout << "Requested energy: " << testEnergyMeV << " MeV\n";
            std::cout << "Matched energy: " << match.monoEnergyMeV() << " MeV\n";
            std::cout << "Matched file: " << match.sourceFile() << "\n";
            std::cout << "Spectrum kind: "
                      << spectrumKindToString(match.spectrumKind()) << "\n";

            if (match.spectrumKind() == LookupSpectrumKind::DeCunhaRawCounts) {
                std::cout << "Detected Ncpp: " << match.ncpp() << "\n";
            } else if (match.spectrumKind() == LookupSpectrumKind::CartechiniFy) {
                std::cout << "Detected Ncpp: N/A\n";
            }

            return 0;
        }

        if (mode == "audit-phsp" || mode == "build-spectrum") {
            if (argc < 6) {
                printUsage(argv[0]);
                throw std::runtime_error(mode + " requires more arguments.");
            }

            const fs::path lookupRoot = argv[2];
            const std::string lutName = argv[3];
            const LutFamily family = LookupLibrary::inferFamily(lutName);

            fs::path primaryLibraryDir;
            fs::path phspFile;
            fs::path outputDir;

            if (family == LutFamily::DeCunha) {
                if (argc != 8) {
                    printUsage(argv[0]);
                    throw std::runtime_error(
                        mode + " for DeCunha requires: "
                        "<lookupRoot> DeCunha <1mm|5um> <linear|log> <phspFile> <outputDir>");
                }

                const std::string voxelSize = argv[4];
                const std::string energyGrid = argv[5];
                phspFile = argv[6];
                outputDir = argv[7];

                primaryLibraryDir = resolveLibraryDir(
                    lookupRoot, lutName, family, voxelSize, energyGrid);
            } else {
                if (argc != 6) {
                    printUsage(argv[0]);
                    throw std::runtime_error(
                        mode + " for Cartechini requires: "
                        "<lookupRoot> Cartechini <phspFile> <outputDir>");
                }

                phspFile = argv[4];
                outputDir = argv[5];
                primaryLibraryDir = resolveLibraryDir(lookupRoot, lutName, family);
            }

            if (!fs::exists(primaryLibraryDir) || !fs::is_directory(primaryLibraryDir)) {
                throw std::runtime_error(
                    "Lookup library folder not found: " + primaryLibraryDir.string());
            }

            LookupLibrary primaryLibrary;
            primaryLibrary.loadFromDirectory(primaryLibraryDir, family);

            LookupLibrary fallbackLibrary;
            fs::path fallbackLibraryDir;
            bool useFallback = false;
            double cartechiniMaxEnergy = -1.0;

            if (family == LutFamily::Cartechini) {
                fallbackLibraryDir = resolveLibraryDir(
                    lookupRoot, "DeCunha", LutFamily::DeCunha, "1mm", "log");

                if (!fs::exists(fallbackLibraryDir) ||
                    !fs::is_directory(fallbackLibraryDir)) {
                    throw std::runtime_error(
                        "Fallback DeCunha library folder not found: " +
                        fallbackLibraryDir.string());
                }

                fallbackLibrary.loadFromDirectory(
                    fallbackLibraryDir, LutFamily::DeCunha);

                if (primaryLibrary.tables().empty()) {
                    throw std::runtime_error(
                        "Cartechini library loaded but contains no tables.");
                }

                cartechiniMaxEnergy =
                    primaryLibrary.tables().back().monoEnergyMeV();
                useFallback = true;
            }

            PhaseSpaceReader reader;
            const std::vector<ProtonRecord> protons = reader.readProtons(phspFile);

            fs::create_directories(outputDir);

            std::vector<ProtonMatchRecord> matches;
            matches.reserve(protons.size());

            SpectrumAccumulator accumulator;

            std::cout << "Phase-space processing started.\n";
            std::cout << "LUT family: " << lutFamilyToString(family) << "\n";
            std::cout << "Library folder: " << primaryLibraryDir << "\n";

            if (useFallback) {
                std::cout
                    << "WARNING: Cartechini LUT is only available up to "
                    << cartechiniMaxEnergy
                    << " MeV. For protons with kinetic energy above this range, "
                    << "DeCunha 1mm logarithmic LUT will be used as fallback.\n";
            }

            for (const auto& proton : protons) {
                const LookupTable* matchPtr = nullptr;
                std::string matchedFamily;

                if (family == LutFamily::Cartechini &&
                    proton.energyMeV > cartechiniMaxEnergy) {
                    matchPtr = &fallbackLibrary.findNearest(proton.energyMeV);
                    matchedFamily = "DeCunha";
                } else {
                    matchPtr = &primaryLibrary.findNearest(proton.energyMeV);
                    matchedFamily = lutFamilyToString(family);
                }

                const LookupTable& match = *matchPtr;

                matches.push_back(ProtonMatchRecord{
                    proton.rowIndex,
                    proton.energyMeV,
                    proton.weight,
                    match.monoEnergyMeV(),
                    match.sourceFile(),
                    match.ncpp(),
                    matchedFamily
                });

                if (mode == "build-spectrum") {
                    accumulator.addContribution(match, 1.0);
                }
            }

            const fs::path matchCsv = outputDir / "proton_matches.csv";
            CsvWriter::writeProtonMatches(matchCsv, matches);

            std::cout << "Phase-space processing completed.\n";
            std::cout << "Phase-space file: " << phspFile << "\n";
            std::cout << "Protons found: " << protons.size() << "\n";
            std::cout << "Match CSV: " << matchCsv << "\n";

            if (mode == "build-spectrum") {
                const PolySpectrum spectrum = accumulator.finalize();
                const fs::path polyCsv = outputDir / "poly_spectrum.csv";
                CsvWriter::writePolySpectrum(polyCsv, spectrum);

                std::cout << "Poly spectrum CSV: " << polyCsv << "\n";
            }

            return 0;
        }

        printUsage(argv[0]);
        throw std::runtime_error("Unknown mode: " + mode);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}