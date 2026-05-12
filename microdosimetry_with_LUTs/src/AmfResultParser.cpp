#include "AmfResultParser.h"

#include <filesystem>
#include <string>

namespace fs = std::filesystem;

namespace {

fs::path scorerOutputPath(const AmfConfig& config,
                          const AmfStagedRun& stagedRun) {
    return stagedRun.runDirectory / (config.outputFile + ".csv");
}

fs::path spectraOutputPath(const AmfConfig& config,
                           const AmfStagedRun& stagedRun) {
    return stagedRun.runDirectory /
           (config.outputFile + "_MicrodosimetricSpectra.csv");
}

}  // namespace

AmfResultFiles AmfResultParser::expectedResultFiles(
    const AmfConfig& config,
    const AmfStagedRun& stagedRun) {

    AmfResultFiles files;
    files.scorerOutputFile = scorerOutputPath(config, stagedRun);

    if (config.quantity == AmfQuantity::Spectra) {
        files.spectraCsvFile = spectraOutputPath(config, stagedRun);
    }

    return files;
}

bool AmfResultParser::hasExpectedResults(const AmfConfig& config,
                                         const AmfStagedRun& stagedRun) {
    const AmfResultFiles files = expectedResultFiles(config, stagedRun);

    if (!fs::exists(files.scorerOutputFile)) {
        return false;
    }

    if (files.spectraCsvFile.has_value() &&
        !fs::exists(*files.spectraCsvFile)) {
        return false;
    }

    return true;
}
