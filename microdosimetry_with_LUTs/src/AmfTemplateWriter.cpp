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

double localScoringXmm(const AmfConfig& config) {
    return config.scoringTransXmm - config.phantomTransXmm;
}

double localScoringYmm(const AmfConfig& config) {
    return config.scoringTransYmm - config.phantomTransYmm;
}

double localScoringZmm(const AmfConfig& config) {
    return config.scoringTransZmm - config.phantomTransZmm;
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

void writeGeometry(std::ostream& out, const AmfConfig& config) {
    out << "# Geometry derived from source TOPAS input\n";
    out << "# Source file: " << config.sourceTopasPath << "\n";
    out << "# Phase-space scorer: Sc/" << config.phaseSpaceScorerName << "\n";
    out << "# Phase-space component: " << config.phaseSpaceComponent << "\n";
    out << "# Phase-space surface: " << config.phaseSpaceSurface << "\n";
    out << "d:Ge/World/HLX = " << config.worldHalfLengthXmm << " mm\n";
    out << "d:Ge/World/HLY = " << config.worldHalfLengthYmm << " mm\n";
    out << "d:Ge/World/HLZ = " << config.worldHalfLengthZmm << " mm\n";
    out << "s:Ge/World/Material = \"" << config.worldMaterial << "\"\n\n";

    out << "# Water phantom copied from the source simulation placement\n";
    out << "s:Ge/" << config.phantomComponent << "/Parent = \"World\"\n";
    out << "s:Ge/" << config.phantomComponent << "/Type = \"TsBox\"\n";
    out << "s:Ge/" << config.phantomComponent << "/Material = \""
        << config.phantomMaterial << "\"\n";
    out << "d:Ge/" << config.phantomComponent << "/HLX = "
        << config.phantomHalfLengthXmm << " mm\n";
    out << "d:Ge/" << config.phantomComponent << "/HLY = "
        << config.phantomHalfLengthYmm << " mm\n";
    out << "d:Ge/" << config.phantomComponent << "/HLZ = "
        << config.phantomHalfLengthZmm << " mm\n";
    out << "d:Ge/" << config.phantomComponent << "/TransX = "
        << config.phantomTransXmm << " mm\n";
    out << "d:Ge/" << config.phantomComponent << "/TransY = "
        << config.phantomTransYmm << " mm\n";
    out << "d:Ge/" << config.phantomComponent << "/TransZ = "
        << config.phantomTransZmm << " mm\n";
    out << "d:Ge/" << config.phantomComponent << "/RotX = 0.0 deg\n";
    out << "d:Ge/" << config.phantomComponent << "/RotY = 0.0 deg\n";
    out << "d:Ge/" << config.phantomComponent << "/RotZ = 0.0 deg\n\n";

    out << "# AMF detector is always water and is centered on the source "
           "phase-space scoring region\n";
    out << "s:Ge/" << config.scoringComponent << "/Parent = \""
        << config.phantomComponent << "\"\n";
    out << "s:Ge/" << config.scoringComponent << "/Type = \"TsSphere\"\n";
    out << "s:Ge/" << config.scoringComponent << "/Material = \"G4_WATER\"\n";
    out << "d:Ge/" << config.scoringComponent << "/RMin = 0 mm\n";
    out << "d:Ge/" << config.scoringComponent << "/RMax = "
        << config.scoringRadiusMm << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/SPhi = 0 deg\n";
    out << "d:Ge/" << config.scoringComponent << "/DPhi = 360 deg\n";
    out << "d:Ge/" << config.scoringComponent << "/STheta = 0 deg\n";
    out << "d:Ge/" << config.scoringComponent << "/DTheta = 180 deg\n";
    out << "d:Ge/" << config.scoringComponent << "/TransX = "
        << localScoringXmm(config) << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/TransY = "
        << localScoringYmm(config) << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/TransZ = "
        << localScoringZmm(config) << " mm\n";
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
    out << "s:Sc/" << scorerName << "/OutputType = \"csv\"\n";
    out << "s:Sc/" << scorerName << "/OutputFile = \""
        << config.outputFile << "\"\n";
    out << "s:Sc/" << scorerName
        << "/IfOutputFileAlreadyExists = \"Overwrite\"\n";
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
    out << "# TOPAS phase-space replay can encounter rare excited-state ion PDG codes.\n";
    out << "b:Ts/TreatExcitedIonsAsGroundState = \"True\"\n\n";

    writeCommonSource(out, config);
    writeGeometry(out, config);
    writeElectronCutNote(out, config);
    writeScorer(out, config);
}
