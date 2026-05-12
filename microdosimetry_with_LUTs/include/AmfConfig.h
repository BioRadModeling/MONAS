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

struct AmfConfig {
    AmfQuantity quantity{AmfQuantity::YD};
    AmfStoppingPowerMode stoppingPowerMode{AmfStoppingPowerMode::Topas};
    AmfStepCalculatorMode stepCalculatorMode{AmfStepCalculatorMode::MidStep};

    std::filesystem::path topasExecutable{"topas"};
    std::filesystem::path phaseSpaceBasePath;
    std::filesystem::path tsedPath;
    std::filesystem::path outputDir;
    std::filesystem::path stagedRunDir;

    std::string scoringComponent{"AMFScoringVolume"};
    std::string scoringMaterial{"G4_WATER"};
    std::string outputFile{"amf_output"};

    double worldHalfLengthCm{20.0};
    double scoringHalfLengthXmm{0.5};
    double scoringHalfLengthYmm{0.5};
    double scoringHalfLengthZmm{0.5};
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
