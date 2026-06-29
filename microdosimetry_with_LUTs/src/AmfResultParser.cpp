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

fs::path spectraMomentsOutputPath(const AmfConfig& config,
                                  const AmfStagedRun& stagedRun) {
    return stagedRun.runDirectory /
           (config.outputFile + "_MicrodosimetricMoments.csv");
}

}  // namespace

AmfResultFiles AmfResultParser::expectedResultFiles(
    const AmfConfig& config,
    const AmfStagedRun& stagedRun) {

    AmfResultFiles files;
    files.scorerOutputFile = scorerOutputPath(config, stagedRun);

    if (config.quantity == AmfQuantity::Spectra) {
        files.spectraCsvFile = spectraOutputPath(config, stagedRun);
        files.spectraMomentsCsvFile = spectraMomentsOutputPath(config, stagedRun);
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

    if (files.spectraMomentsCsvFile.has_value() &&
        !fs::exists(*files.spectraMomentsCsvFile)) {
        return false;
    }

    return true;
}
