#include <cctype>
#include <cstdint>
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

namespace {

std::string toLower(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

LutFamily parseLutFamily(const std::string& lutName) {
    const std::string lower = toLower(lutName);

    if (lower == "cartechini") {
        return LutFamily::Cartechini;
    }
    if (lower == "decunha" || lower == "de_cunha" || lower == "de-cunha") {
        return LutFamily::DeCunha;
    }

    throw std::runtime_error(
        "Unknown LUT family '" + lutName +
        "'. Expected Cartechini or DeCunha.");
}

std::string resolveFolderName(const std::string& voxelSize,
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

fs::path resolveLibraryDir(const fs::path& lookupRoot,
                           const std::string& lutName,
                           LutFamily family,
                           const std::string& voxelSize,
                           const std::string& energyGrid) {
    if (family == LutFamily::Cartechini) {
        return lookupRoot / lutName;
    }

    return lookupRoot / lutName / resolveFolderName(voxelSize, energyGrid);
}

void parseOptionalBuildArgs(int argc,
                            char* argv[],
                            int startIndex,
                            std::size_t& rebinSamples,
                            std::uint64_t& rebinSeed) {
    for (int i = startIndex; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "--rebin-samples") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --rebin-samples");
            }
            rebinSamples = static_cast<std::size_t>(std::stoull(argv[++i]));
            if (rebinSamples == 0) {
                throw std::runtime_error("--rebin-samples must be > 0");
            }
        } else if (arg == "--rebin-seed") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value after --rebin-seed");
            }
            rebinSeed = static_cast<std::uint64_t>(std::stoull(argv[++i]));
        } else {
            throw std::runtime_error("Unknown option: " + arg);
        }
    }
}

void printUsage(const char* programName) {
    std::cerr
        << "Usage:\n"
        << "  Cartechini:\n"
        << "    " << programName
        << " test-lookup <lookupRoot> Cartechini <testEnergyMeV>\n"
        << "    " << programName
        << " audit-phsp <lookupRoot> Cartechini <phspFile> <outputDir>\n"
        << "    " << programName
        << " build-spectrum <lookupRoot> Cartechini <phspFile> <outputDir>"
        << " [--rebin-samples N] [--rebin-seed S]\n"
        << "\n"
        << "  DeCunha:\n"
        << "    " << programName
        << " test-lookup <lookupRoot> DeCunha <1mm|5um> <linear|log> <testEnergyMeV>\n"
        << "    " << programName
        << " audit-phsp <lookupRoot> DeCunha <1mm|5um> <linear|log> <phspFile> <outputDir>\n"
        << "    " << programName
        << " build-spectrum <lookupRoot> DeCunha <1mm|5um> <linear|log> <phspFile> <outputDir>"
        << " [--rebin-samples N] [--rebin-seed S]\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            printUsage(argv[0]);
            return 1;
        }

        const std::string mode = argv[1];

        if (mode == "test-lookup") {
            if (argc < 5) {
                throw std::runtime_error("Insufficient arguments for test-lookup.");
            }

            const fs::path lookupRoot = argv[2];
            const std::string lutName = argv[3];
            const LutFamily family = parseLutFamily(lutName);

            fs::path libraryDir;
            double testEnergyMeV = 0.0;

            if (family == LutFamily::Cartechini) {
                if (argc != 5) {
                    throw std::runtime_error(
                        "Cartechini test-lookup requires: <lookupRoot> Cartechini <testEnergyMeV>");
                }

                testEnergyMeV = std::stod(argv[4]);
                libraryDir = resolveLibraryDir(lookupRoot, lutName, family, "", "");
            } else {
                if (argc != 7) {
                    throw std::runtime_error(
                        "DeCunha test-lookup requires: <lookupRoot> DeCunha <1mm|5um> <linear|log> <testEnergyMeV>");
                }

                const std::string voxelSize = argv[4];
                const std::string energyGrid = argv[5];
                testEnergyMeV = std::stod(argv[6]);
                libraryDir = resolveLibraryDir(lookupRoot, lutName, family, voxelSize, energyGrid);
            }

            if (!fs::exists(libraryDir) || !fs::is_directory(libraryDir)) {
                throw std::runtime_error(
                    "Lookup library folder not found: " + libraryDir.string());
            }

            LookupLibrary library;
            library.loadFromDirectory(libraryDir, family);

            const LookupTable& match = library.findNearest(testEnergyMeV);

            std::cout << "Loaded library successfully.\n";
            std::cout << "Library folder: " << libraryDir << "\n";
            std::cout << "Number of tables: " << library.size() << "\n";
            std::cout << "Y bins: " << library.yReference().size() << "\n";
            std::cout << "Requested energy: " << testEnergyMeV << " MeV\n";
            std::cout << "Matched energy: " << match.monoEnergyMeV() << " MeV\n";
            std::cout << "Matched CSV: " << match.sourceFile() << "\n";
            std::cout << "Detected Ncpp: " << match.ncpp() << "\n";
            return 0;
        }

        if (mode == "audit-phsp" || mode == "build-spectrum") {
            if (argc < 6) {
                throw std::runtime_error("Insufficient arguments.");
            }

            const fs::path lookupRoot = argv[2];
            const std::string lutName = argv[3];
            const LutFamily family = parseLutFamily(lutName);

            fs::path libraryDir;
            fs::path phspFile;
            fs::path outputDir;
            int firstOptionalArgIndex = argc;

            if (family == LutFamily::Cartechini) {
                if (mode == "audit-phsp") {
                    if (argc != 6) {
                        throw std::runtime_error(
                            "Cartechini audit-phsp requires: <lookupRoot> Cartechini <phspFile> <outputDir>");
                    }
                } else {
                    if (argc < 6) {
                        throw std::runtime_error(
                            "Cartechini build-spectrum requires: <lookupRoot> Cartechini <phspFile> <outputDir> [--rebin-samples N] [--rebin-seed S]");
                    }
                }

                libraryDir = resolveLibraryDir(lookupRoot, lutName, family, "", "");
                phspFile = argv[4];
                outputDir = argv[5];
                firstOptionalArgIndex = 6;
            } else {
                if (mode == "audit-phsp") {
                    if (argc != 8) {
                        throw std::runtime_error(
                            "DeCunha audit-phsp requires: <lookupRoot> DeCunha <1mm|5um> <linear|log> <phspFile> <outputDir>");
                    }
                } else {
                    if (argc < 8) {
                        throw std::runtime_error(
                            "DeCunha build-spectrum requires: <lookupRoot> DeCunha <1mm|5um> <linear|log> <phspFile> <outputDir> [--rebin-samples N] [--rebin-seed S]");
                    }
                }

                const std::string voxelSize = argv[4];
                const std::string energyGrid = argv[5];
                libraryDir = resolveLibraryDir(lookupRoot, lutName, family, voxelSize, energyGrid);
                phspFile = argv[6];
                outputDir = argv[7];
                firstOptionalArgIndex = 8;
            }

            std::size_t rebinSamples = 1000000;
            std::uint64_t rebinSeed = 0x5A17C3E4ULL;

            if (mode == "build-spectrum") {
                parseOptionalBuildArgs(
                    argc, argv, firstOptionalArgIndex, rebinSamples, rebinSeed);
            }

            if (!fs::exists(libraryDir) || !fs::is_directory(libraryDir)) {
                throw std::runtime_error(
                    "Lookup library folder not found: " + libraryDir.string());
            }

            LookupLibrary library;
            library.loadFromDirectory(libraryDir, family);

            PhaseSpaceReader reader;
            const std::vector<ProtonRecord> protons = reader.readProtons(phspFile);

            fs::create_directories(outputDir);

            std::vector<ProtonMatchRecord> matches;
            matches.reserve(protons.size());

            SpectrumAccumulator accumulator(library, rebinSamples, rebinSeed);
            std::vector<std::size_t> matchCounts(library.size(), 0);

            for (const auto& proton : protons) {
                const std::size_t matchIdx = library.findNearestIndex(proton.energyMeV);
                const LookupTable& match = library.tables()[matchIdx];

                const std::string matchedFamily =
                    (match.spectrumKind() == LookupSpectrumKind::CartechiniFy)
                        ? "Cartechini"
                        : "DeCunha";

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
                    matchCounts[matchIdx] += 1;
                }
            }

            const fs::path matchCsv = outputDir / "proton_matches.csv";
            CsvWriter::writeProtonMatches(matchCsv, matches);

            std::cout << "Phase-space processing completed.\n";
            std::cout << "Library folder: " << libraryDir << "\n";
            std::cout << "Phase-space file: " << phspFile << "\n";
            std::cout << "Protons found: " << protons.size() << "\n";
            std::cout << "Match CSV: " << matchCsv << "\n";

            if (mode == "build-spectrum") {
                std::cout << "Rebin samples per LUT: " << rebinSamples << "\n";
                std::cout << "Rebin seed: " << rebinSeed << "\n";

                for (std::size_t i = 0; i < matchCounts.size(); ++i) {
                    if (matchCounts[i] == 0) {
                        continue;
                    }

                    accumulator.addContributionByIndex(
                        i, static_cast<double>(matchCounts[i]));
                }

                const PolySpectrum spectrum = accumulator.finalize();
                const fs::path polyCsv = outputDir / "poly_spectrum.csv";
                CsvWriter::writePolySpectrum(polyCsv, spectrum);

                std::cout << "Poly spectrum CSV: " << polyCsv << "\n";
            }

            return 0;
        }

        throw std::runtime_error("Unknown mode: " + mode);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}