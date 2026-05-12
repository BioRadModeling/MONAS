#pragma once

#include <filesystem>
#include <optional>

#include "AmfConfig.h"
#include "AmfRunner.h"

struct AmfResultFiles {
    std::filesystem::path scorerOutputFile;
    std::optional<std::filesystem::path> spectraCsvFile;
};

class AmfResultParser {
public:
    static AmfResultFiles expectedResultFiles(const AmfConfig& config,
                                              const AmfStagedRun& stagedRun);

    static bool hasExpectedResults(const AmfConfig& config,
                                   const AmfStagedRun& stagedRun);
};
