#pragma once

#include <filesystem>
#include <string>

enum class AmfQuantity {
    Spectra,
    YD,
    YS
};

enum class AmfStoppingPowerMode {
    Topas,
    ExternalTable
};

enum class AmfStepCalculatorMode {
    MidStep,
    PreStep
};

enum class AmfDetectorType {
    Water,
    Silicon,
    TEGas
};

inline constexpr double kMinAmfDomainRadiusUm = 0.0015;
inline constexpr double kMaxAmfDomainRadiusUm = 0.5;

struct AmfConfig {
    AmfQuantity quantity{AmfQuantity::YD};
    AmfStoppingPowerMode stoppingPowerMode{AmfStoppingPowerMode::Topas};
    AmfStepCalculatorMode stepCalculatorMode{AmfStepCalculatorMode::MidStep};
    AmfDetectorType detectorType{AmfDetectorType::Water};

    std::filesystem::path topasExecutable{"topas"};
    std::filesystem::path phaseSpaceBasePath;
    std::filesystem::path tsedPath;
    std::filesystem::path outputDir;
    std::filesystem::path stagedRunDir;

    std::string scoringComponent{"AMFScoringVolume"};
    std::string scoringMaterial{"G4_WATER"};
    std::string outputFile{"amf_output"};

    double worldHalfLengthCm{20.0};
    double scoringHalfLengthXmm{50.0};
    double scoringHalfLengthYmm{50.0};
    double scoringHalfLengthZmm{0.5};
    double scoringRadiusMm{6.35};
    double scoringTransXmm{0.0};
    double scoringTransYmm{0.0};
    double scoringTransZmm{0.0};

    double domainRadiusUm{0.28};
    double nucleusRadiusUm{3.9};
    double betaRefPerGy2{0.0615};
    double electronRangeCutM{1000.0};
    bool phaseSpacePreCheck{true};
};

inline const char* toTopasQuantityName(AmfQuantity quantity) {
    switch (quantity) {
        case AmfQuantity::Spectra:
            return "AMFSpectra";
        case AmfQuantity::YD:
            return "AMF_yD";
        case AmfQuantity::YS:
            return "AMF_yS";
    }

    return "AMF_yD";
}

inline const char* toTopasStoppingPowerModeName(AmfStoppingPowerMode mode) {
    switch (mode) {
        case AmfStoppingPowerMode::Topas:
            return "Topas";
        case AmfStoppingPowerMode::ExternalTable:
            return "ExternalTable";
    }

    return "Topas";
}

inline const char* toTopasStepCalculatorModeName(AmfStepCalculatorMode mode) {
    switch (mode) {
        case AmfStepCalculatorMode::MidStep:
            return "MidStep";
        case AmfStepCalculatorMode::PreStep:
            return "PreStep";
    }

    return "MidStep";
}

inline const char* toAmfDetectorTypeName(AmfDetectorType detectorType) {
    switch (detectorType) {
        case AmfDetectorType::Water:
            return "water";
        case AmfDetectorType::Silicon:
            return "silicon";
        case AmfDetectorType::TEGas:
            return "TEgas";
    }

    return "water";
}

inline bool isValidAmfDomainRadiusUm(double radiusUm) {
    return radiusUm >= kMinAmfDomainRadiusUm &&
           radiusUm <= kMaxAmfDomainRadiusUm;
}

inline void applyAmfDetectorPreset(AmfConfig& config,
                                   AmfDetectorType detectorType) {
    config.detectorType = detectorType;

    switch (detectorType) {
        case AmfDetectorType::Water:
            config.scoringComponent = "AMFScoringVolume";
            config.scoringMaterial = "G4_WATER";
            config.scoringHalfLengthXmm = 50.0;
            config.scoringHalfLengthYmm = 50.0;
            config.scoringHalfLengthZmm = 0.5;
            break;
        case AmfDetectorType::Silicon:
            config.scoringComponent = "SOISensitiveLayer";
            config.scoringMaterial = "G4_Si";
            config.scoringHalfLengthXmm = 1.465;
            config.scoringHalfLengthYmm = 1.790;
            config.scoringHalfLengthZmm = 0.005;
            break;
        case AmfDetectorType::TEGas:
            config.scoringComponent = "TEgasSV";
            config.scoringMaterial = "PropaneGas";
            config.scoringRadiusMm = 6.35;
            break;
    }
}
