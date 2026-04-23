#include <cctype>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CsvWriter.h"
#include "LetCalculator.h"
#include "MaginiCalculator.h"
#include "MaginiLookup.h"
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

std::string resolveCartechiniRadiusFolder(const std::string& radius) {
    const std::string normalized = toLower(radius);

    if (normalized == "0.5" || normalized == "0.5um" ||
        normalized == "0.5_um" || normalized == "0.5-um" ||
        normalized == "r0.5") {
        return "R0.5";
    }

    if (normalized == "8" || normalized == "8.0" ||
        normalized == "8um" || normalized == "8.0um" ||
        normalized == "8_um" || normalized == "8.0_um" ||
        normalized == "8-um" || normalized == "8.0-um" ||
        normalized == "r8" || normalized == "r8.0") {
        return "R8.0";
    }

    throw std::runtime_error(
        "Invalid Cartechini scoring radius '" + radius +
        "'. Expected 0.5um or 8um.\n"
        "Usage: ./microdosimetry_with_LUTs build-spectrum <lookupRoot> "
        "Cartechini <0.5um|8um> <phspFile> <outputDir> [options]");
}

fs::path resolveLibraryDir(const fs::path& lookupRoot,
                           const std::string& lutName,
                           LutFamily family,
                           const std::string& voxelSize,
                           const std::string& energyGrid) {
    if (family == LutFamily::Cartechini) {
        return lookupRoot / lutName / resolveCartechiniRadiusFolder(voxelSize);
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
        << " test-lookup <lookupRoot> Cartechini <0.5um|8um> <testEnergyMeV>\n"
        << "    " << programName
        << " audit-phsp <lookupRoot> Cartechini <0.5um|8um> <phspFile> <outputDir>\n"
        << "    " << programName
        << " build-spectrum <lookupRoot> Cartechini <0.5um|8um> <phspFile> <outputDir>"
        << " [--rebin-samples N] [--rebin-seed S]\n"
        << "\n"
        << "  DeCunha:\n"
        << "    " << programName
        << " test-lookup <lookupRoot> DeCunha <1mm|5um> <linear|log> <testEnergyMeV>\n"
        << "    " << programName
        << " audit-phsp <lookupRoot> DeCunha <1mm|5um> <linear|log> <phspFile> <outputDir>\n"
        << "    " << programName
        << " build-spectrum <lookupRoot> DeCunha <1mm|5um> <linear|log> <phspFile> <outputDir>"
        << " [--rebin-samples N] [--rebin-seed S]\n"
        << "\n"
        << "  LET:\n"
        << "    " << programName
        << " LET <lookupRoot> <phspFile> <outputDir>\n"
        << "\n"
        << "  Magini:\n"
        << "    " << programName
        << " Magini <lookupRoot> <phspFile> <outputDir>\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            printUsage(argv[0]);
            return 1;
        }

        const std::string mode = argv[1];

        if (mode == "LET") {
            if (argc != 5) {
                throw std::runtime_error(
                    "LET requires: <lookupRoot> <phspFile> <outputDir>");
            }

            const fs::path lookupRoot = argv[2];
            const fs::path letDirectory = lookupRoot / "LET";
            const fs::path phspFile = argv[3];
            const fs::path outputDir = argv[4];

            if (!fs::exists(letDirectory) || !fs::is_directory(letDirectory)) {
                throw std::runtime_error(
                    "LET lookup folder not found: " + letDirectory.string());
            }

            PhaseSpaceReader reader;
            const std::vector<ChargedParticleRecord> chargedParticles =
                reader.readChargedParticles(phspFile);

            LetCalculator calculator;
            const std::vector<ElementLetSummary> summaries =
                calculator.calculateByElement(letDirectory, chargedParticles);

            CsvWriter::writeElementLetSummaries(outputDir, summaries);

            std::cout << "LET calculation completed.\n";
            std::cout << "LET folder: " << letDirectory << "\n";
            std::cout << "Phase-space file: " << phspFile << "\n";
            std::cout << "Charged particles found: " << chargedParticles.size() << "\n";
            std::cout << "Output directory: " << outputDir << "\n";
            std::cout << "Element summaries written: " << summaries.size() << "\n";
            return 0;
        }

        if (mode == "Magini") {
            if (argc != 5) {
                throw std::runtime_error(
                    "Magini requires: <lookupRoot> <phspFile> <outputDir>");
            }

            const fs::path lookupRoot = argv[2];
            const fs::path maginiCsv = lookupRoot / "Magini" / "Magini.csv";
            const fs::path phspFile = argv[3];
            const fs::path outputDir = argv[4];

            if (!fs::exists(maginiCsv) || !fs::is_regular_file(maginiCsv)) {
                throw std::runtime_error(
                    "Magini lookup CSV not found: " + maginiCsv.string());
            }

            PhaseSpaceReader reader;
            const std::vector<ProtonRecord> protons = reader.readProtons(phspFile);

            const MaginiLookup lookup = MaginiLookup::loadFromCsv(maginiCsv);

            MaginiCalculator calculator;
            const MaginiSummary summary = calculator.calculate(lookup, protons);
            fs::create_directories(outputDir);
            const fs::path summaryCsv = outputDir / "magini_summary.csv";
            CsvWriter::writeMaginiSummary(summaryCsv, summary);

            std::cout << "Magini calculation completed.\n";
            std::cout << "Phase-space file: " << phspFile << "\n";
            std::cout << "Output directory: " << outputDir << "\n";
            std::cout << "Protons found: " << summary.protonCount << "\n";
            std::cout << "Total proton weight: " << summary.totalProtonWeight << "\n";
            std::cout << "y_F [keV/um]: " << summary.yFKeVPerUm << "\n";
            std::cout << "y_D [keV/um]: " << summary.yDKeVPerUm << "\n";
            std::cout << "y* [keV/um]: " << summary.yStarKeVPerUm << "\n";
            std::cout << "Magini summary CSV: " << summaryCsv << "\n";
            return 0;
        }

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
                if (argc != 6) {
                    throw std::runtime_error(
                        "Cartechini test-lookup requires: <lookupRoot> Cartechini <0.5um|8um> <testEnergyMeV>");
                }

                const std::string scoringRadius = argv[4];
                testEnergyMeV = std::stod(argv[5]);
                libraryDir = resolveLibraryDir(lookupRoot, lutName, family,
                                              scoringRadius, "");
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

            const InterpolationMatch match =
                library.findInterpolationMatch(testEnergyMeV);
            const LookupTable& lowerMatch = library.tables()[match.lowerIndex];
            const LookupTable& upperMatch = library.tables()[match.upperIndex];

            std::cout << "Loaded library successfully.\n";
            std::cout << "Library folder: " << libraryDir << "\n";
            std::cout << "Number of tables: " << library.size() << "\n";
            std::cout << "Y bins: " << library.yReference().size() << "\n";
            std::cout << "Requested energy: " << testEnergyMeV << " MeV\n";
            std::cout << "Lower matched energy: " << lowerMatch.monoEnergyMeV() << " MeV\n";
            std::cout << "Upper matched energy: " << upperMatch.monoEnergyMeV() << " MeV\n";
            std::cout << "Lower weight: " << match.lowerWeight << "\n";
            std::cout << "Upper weight: " << match.upperWeight << "\n";
            std::cout << "Lower matched file: " << lowerMatch.sourceFile() << "\n";
            std::cout << "Upper matched file: " << upperMatch.sourceFile() << "\n";
            std::cout << "Lower detected Ncpp: " << lowerMatch.ncpp() << "\n";
            std::cout << "Upper detected Ncpp: " << upperMatch.ncpp() << "\n";
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
                    if (argc != 7) {
                        throw std::runtime_error(
                            "Cartechini audit-phsp requires: <lookupRoot> Cartechini <0.5um|8um> <phspFile> <outputDir>");
                    }
                } else {
                    if (argc < 7) {
                        throw std::runtime_error(
                            "Cartechini build-spectrum requires: <lookupRoot> Cartechini <0.5um|8um> <phspFile> <outputDir> [--rebin-samples N] [--rebin-seed S]");
                    }
                }

                const std::string scoringRadius = argv[4];
                libraryDir = resolveLibraryDir(lookupRoot, lutName, family,
                                              scoringRadius, "");
                phspFile = argv[5];
                outputDir = argv[6];
                firstOptionalArgIndex = 7;
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

            for (const auto& proton : protons) {
                const InterpolationMatch match =
                    library.findInterpolationMatch(proton.energyMeV);
                const LookupTable& lowerMatch = library.tables()[match.lowerIndex];
                const LookupTable& upperMatch = library.tables()[match.upperIndex];

                const std::string lowerMatchedFamily =
                    (lowerMatch.spectrumKind() == LookupSpectrumKind::CartechiniFy)
                        ? "Cartechini"
                        : "DeCunha";
                const std::string upperMatchedFamily =
                    (upperMatch.spectrumKind() == LookupSpectrumKind::CartechiniFy)
                        ? "Cartechini"
                        : "DeCunha";

                matches.push_back(ProtonMatchRecord{
                    proton.rowIndex,
                    proton.energyMeV,
                    proton.weight,
                    lowerMatch.monoEnergyMeV(),
                    upperMatch.monoEnergyMeV(),
                    match.lowerWeight,
                    match.upperWeight,
                    lowerMatch.sourceFile(),
                    upperMatch.sourceFile(),
                    lowerMatch.ncpp(),
                    upperMatch.ncpp(),
                    lowerMatchedFamily,
                    upperMatchedFamily
                });

                if (mode == "build-spectrum") {
                    accumulator.addInterpolatedContribution(
                        match.lowerIndex,
                        match.upperIndex,
                        match.lowerWeight,
                        match.upperWeight,
                        proton.weight);
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
