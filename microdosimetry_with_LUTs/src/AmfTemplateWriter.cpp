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

void writeBoxDetector(std::ostream& out, const AmfConfig& config) {
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

void writeTEGasMaterial(std::ostream& out, const AmfConfig& config) {
    if (config.scoringMaterial != "PropaneGas") {
        return;
    }

    out << "# Tissue-equivalent propane gas\n";
    out << "sv:Ma/PropaneGas/Components = 2 \"Hydrogen\" \"Carbon\"\n";
    out << "uv:Ma/PropaneGas/Fractions = 2 0.182864 0.817136\n";
    out << "d:Ma/PropaneGas/Density = 0.108 mg/cm3\n\n";
}

void writeSphereDetector(std::ostream& out, const AmfConfig& config) {
    out << "s:Ge/" << config.scoringComponent << "/Parent = \"World\"\n";
    out << "s:Ge/" << config.scoringComponent << "/Type = \"TsSphere\"\n";
    out << "s:Ge/" << config.scoringComponent << "/Material = \""
        << config.scoringMaterial << "\"\n";
    out << "d:Ge/" << config.scoringComponent << "/RMin = 0 mm\n";
    out << "d:Ge/" << config.scoringComponent << "/RMax = "
        << config.scoringRadiusMm << " mm\n";
    out << "d:Ge/" << config.scoringComponent << "/SPhi = 0 deg\n";
    out << "d:Ge/" << config.scoringComponent << "/DPhi = 360 deg\n";
    out << "d:Ge/" << config.scoringComponent << "/STheta = 0 deg\n";
    out << "d:Ge/" << config.scoringComponent << "/DTheta = 180 deg\n";
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

void writeGeometryPlaceholder(std::ostream& out, const AmfConfig& config) {
    out << "# Geometry\n";
    out << "# Detector preset: "
        << toAmfDetectorTypeName(config.detectorType) << "\n";
    out << "d:Ge/World/HLX = " << config.worldHalfLengthCm << " cm\n";
    out << "d:Ge/World/HLY = " << config.worldHalfLengthCm << " cm\n";
    out << "d:Ge/World/HLZ = " << config.worldHalfLengthCm << " cm\n";
    out << "s:Ge/World/Material = \"Air\"\n\n";

    if (config.detectorType == AmfDetectorType::TEGas) {
        writeTEGasMaterial(out, config);
        writeSphereDetector(out, config);
    } else {
        writeBoxDetector(out, config);
    }
}

void writeDetectorNotes(std::ostream& out, const AmfConfig& config) {
    out << "# Detector notes\n";

    switch (config.detectorType) {
        case AmfDetectorType::Water:
            out << "# water: 10 cm x 10 cm x 1 mm G4_WATER slab by default.\n";
            break;
        case AmfDetectorType::Silicon:
            out << "# silicon: SOI active-layer approximation, "
                << "2.93 mm x 3.58 mm x 10 um G4_Si by default.\n";
            break;
        case AmfDetectorType::TEGas:
            out << "# TEgas: 6.35 mm radius propane-gas sphere by default.\n";
            break;
    }

    out << "\n";
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
    writeDetectorNotes(out, config);
    writeElectronCutNote(out, config);
    writeScorer(out, config);
}
