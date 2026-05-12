#include "AmfTemplateWriter.h"

#include <fstream>
#include <stdexcept>
#include <string>

namespace {

std::string boolToTopas(bool value) {
    return value ? "True" : "False";
}

std::string phaseSpaceBasename(const std::filesystem::path& phaseSpaceBasePath) {
    return phaseSpaceBasePath.filename().string();
}

void writeCommonSource(std::ostream& out, const AmfConfig& config) {
    out << "# Phase-space source\n";
    out << "s:So/Replay/Type = \"PhaseSpace\"\n";
    out << "s:So/Replay/Component = \"World\"\n";
    out << "s:So/Replay/PhaseSpaceFileName = \""
        << phaseSpaceBasename(config.phaseSpaceBasePath) << "\"\n";
    out << "b:So/Replay/PhaseSpacePreCheck = \""
        << boolToTopas(config.phaseSpacePreCheck) << "\"\n\n";
}

void writeGeometryPlaceholder(std::ostream& out, const AmfConfig& config) {
    out << "# Geometry\n";
    out << "d:Ge/World/HLX = " << config.worldHalfLengthCm << " cm\n";
    out << "d:Ge/World/HLY = " << config.worldHalfLengthCm << " cm\n";
    out << "d:Ge/World/HLZ = " << config.worldHalfLengthCm << " cm\n";
    out << "s:Ge/World/Material = \"Air\"\n\n";

    out << "s:Ge/" << config.scoringComponent << "/Parent = \"World\"\n";
    out << "s:Ge/" << config.scoringComponent << "/Type = \"TsBox\"\n";
    out << "s:Ge/" << config.scoringComponent << "/Material = \""
        << config.scoringMaterial << "\"\n";
    out << "d:Ge/" << config.scoringComponent << "/HLX = "
        << config.scoringHalfLengthXmm << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/HLY = "
        << config.scoringHalfLengthYmm << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/HLZ = "
        << config.scoringHalfLengthZmm << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/TransX = "
        << config.scoringTransXmm << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/TransY = "
        << config.scoringTransYmm << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/TransZ = "
        << config.scoringTransZmm << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/RotX = 0.0 deg\n";
    out << "d:Ge/" << config.scoringComponent << "/RotY = 0.0 deg\n";
    out << "d:Ge/" << config.scoringComponent << "/RotZ = 0.0 deg\n\n";
}

void writeElectronCutNote(std::ostream& out, const AmfConfig& config) {
    out << "# Electron production cut\n";
    out << "# AMF spectra include secondary-electron/delta-ray contributions\n";
    out << "# analytically. Use a high electron range cut to suppress explicit\n";
    out << "# secondary-electron transport and avoid double counting.\n";
    out << "d:Ph/Default/CutForElectron = "
        << config.electronRangeCutM << " m\n\n";
}

void writeScorer(std::ostream& out, const AmfConfig& config) {
    const char* scorerName = "AMF";

    out << "# AMF scorer\n";
    out << "s:Sc/" << scorerName << "/Quantity = \""
        << toTopasQuantityName(config.quantity) << "\"\n";
    out << "s:Sc/" << scorerName << "/Component = \""
        << config.scoringComponent << "\"\n";
    out << "s:Sc/" << scorerName << "/OutputFile = \""
        << config.outputFile << "\"\n";
    out << "d:Sc/" << scorerName << "/DomainRadius = "
        << config.domainRadiusUm << " um\n";

    if (config.quantity == AmfQuantity::YS) {
        out << "d:Sc/" << scorerName << "/NucleusRadius = "
            << config.nucleusRadiusUm << " um\n";
        out << "d:Sc/" << scorerName << "/BetaRef = "
            << config.betaRefPerGy2 << " /Gy2\n";
    }

    out << "s:Sc/" << scorerName << "/StoppingPowerCalculation = \""
        << toTopasStoppingPowerModeName(config.stoppingPowerMode) << "\"\n";
    out << "s:Sc/" << scorerName << "/StepCalculator = \""
        << toTopasStepCalculatorModeName(config.stepCalculatorMode) << "\"\n";
}

}  // namespace

void AmfTemplateWriter::writeReplayParameterFile(
    const AmfConfig& config,
    const std::filesystem::path& outPath) {

    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open AMF TOPAS parameter file for writing: " +
            outPath.string());
    }

    out << "# Generated TOPAS AMF phase-space replay file\n";
    out << "# Keep tsed.dat in the same directory where this file is launched.\n\n";

    writeCommonSource(out, config);
    writeGeometryPlaceholder(out, config);
    writeElectronCutNote(out, config);
    writeScorer(out, config);
}
