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

inline constexpr double kMinAmfDomainRadiusUm = 0.0015;
inline constexpr double kMaxAmfDomainRadiusUm = 0.5;

struct AmfConfig {
    AmfQuantity quantity{AmfQuantity::YD};
    AmfStoppingPowerMode stoppingPowerMode{AmfStoppingPowerMode::Topas};
    AmfStepCalculatorMode stepCalculatorMode{AmfStepCalculatorMode::MidStep};

    std::filesystem::path topasExecutable{"topas"};
    std::filesystem::path phaseSpaceBasePath;
    std::filesystem::path sourceTopasPath;
    std::filesystem::path tsedPath;
    std::filesystem::path stagedRunDir;

    std::string outputFile{"amf_output"};
    std::string phaseSpaceScorerName;
    std::string phaseSpaceComponent;
    std::string phaseSpaceSurface;

    std::string worldMaterial{"Air"};
    double worldHalfLengthXmm{3000.0};
    double worldHalfLengthYmm{3000.0};
    double worldHalfLengthZmm{3000.0};

    std::string phantomComponent{"Phantom"};
    std::string phantomMaterial{"G4_WATER"};
    double phantomHalfLengthXmm{150.0};
    double phantomHalfLengthYmm{150.0};
    double phantomHalfLengthZmm{150.0};
    double phantomTransXmm{0.0};
    double phantomTransYmm{0.0};
    double phantomTransZmm{-1500.0};

    std::string scoringComponent{"AMFScoringVolume"};
    std::string scoringMaterial{"G4_WATER"};
    double scoringRadiusMm{6.35};
    double scoringTransXmm{0.0};
    double scoringTransYmm{0.0};
    double scoringTransZmm{0.0};

    bool hasPhaseSpaceBounds{false};
    double phaseSpaceMinXmm{0.0};
    double phaseSpaceMaxXmm{0.0};
    double phaseSpaceMinYmm{0.0};
    double phaseSpaceMaxYmm{0.0};
    double phaseSpaceMinZmm{0.0};
    double phaseSpaceMaxZmm{0.0};
    double geometryToleranceMm{1.0};

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

inline bool isValidAmfDomainRadiusUm(double radiusUm) {
    return radiusUm >= kMinAmfDomainRadiusUm &&
           radiusUm <= kMaxAmfDomainRadiusUm;
}
