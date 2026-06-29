#pragma once

#include <filesystem>

#include "AmfConfig.h"

struct AmfStagedRun {
    std::filesystem::path runDirectory;
    std::filesystem::path parameterFile;
    std::filesystem::path manifestFile;
    std::filesystem::path stagedTsedPath;
    std::filesystem::path stagedPhaseSpacePath;
    std::filesystem::path stagedHeaderPath;
    std::size_t phaseSpaceEnergyFlooredRows{0};
};

struct AmfRunResult {
    AmfStagedRun stagedRun;
    int exitCode{0};
    std::filesystem::path stdoutLog;
    std::filesystem::path stderrLog;
};

class AmfRunner {
public:
    static AmfStagedRun stageRun(const AmfConfig& config);
    static AmfRunResult runTopas(const AmfConfig& config);
};
