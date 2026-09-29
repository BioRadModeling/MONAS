#pragma once

#include <filesystem>
#include <string>

struct ProgramOptions {
    std::filesystem::path phspPath;
    std::filesystem::path headerPath;
    std::filesystem::path lookupRoot;
    std::filesystem::path outputDir;
    std::string voxelSize;   // "1mm" or "5um"
    std::string energyGrid;  // "linear" or "log"
};