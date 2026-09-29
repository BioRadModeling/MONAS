#pragma once

#include <filesystem>

#include "AmfConfig.h"

class AmfTemplateWriter {
public:
    static void writeReplayParameterFile(const AmfConfig& config,
                                         const std::filesystem::path& outPath);
};
