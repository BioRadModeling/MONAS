#include "AmfRunner.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sys/wait.h>

#include "AmfTemplateWriter.h"

namespace fs = std::filesystem;

namespace {

void requireRegularFile(const fs::path& path, const char* label) {
    if (!fs::exists(path) || !fs::is_regular_file(path)) {
        throw std::runtime_error(
            std::string(label) + " not found: " + path.string());
    }
}

fs::path requireStagedRunDirectory(const AmfConfig& config) {
    if (config.stagedRunDir.empty()) {
        throw std::runtime_error("AMF staged run directory is required.");
    }

    fs::create_directories(config.stagedRunDir);
    return config.stagedRunDir;
}

struct PhaseSpacePair {
    fs::path phaseSpacePath;
    fs::path headerPath;
};

PhaseSpacePair validatePhaseSpacePair(const fs::path& phaseSpaceBasePath) {
    if (phaseSpaceBasePath.empty()) {
        throw std::runtime_error("AMF phase-space base path is required.");
    }

    PhaseSpacePair pair{
        phaseSpaceBasePath.string() + ".phsp",
        phaseSpaceBasePath.string() + ".header"
    };

    requireRegularFile(pair.phaseSpacePath, "AMF phase-space file");
    requireRegularFile(pair.headerPath, "AMF phase-space header");

    return pair;
}

void writeManifest(const AmfConfig& config,
                   const fs::path& manifestFile,
                   const AmfStagedRun& stagedRun) {
    std::ofstream out(manifestFile);
    if (!out) {
        throw std::runtime_error(
            "Failed to open AMF manifest for writing: " +
            manifestFile.string());
    }

    out << "AMF staged run manifest\n";
    out << "run_directory = " << stagedRun.runDirectory << "\n";
    out << "parameter_file = " << stagedRun.parameterFile << "\n";
    out << "manifest_file = " << stagedRun.manifestFile << "\n";
    out << "tsed_dat = " << stagedRun.stagedTsedPath << "\n";
    out << "phase_space = " << stagedRun.stagedPhaseSpacePath << "\n";
    out << "phase_space_header = " << stagedRun.stagedHeaderPath << "\n";
    out << "quantity = " << toTopasQuantityName(config.quantity) << "\n";
    out << "output_file = " << config.outputFile << "\n";
    out << "detector = " << toAmfDetectorTypeName(config.detectorType) << "\n";
    out << "scoring_component = " << config.scoringComponent << "\n";
    out << "scoring_material = " << config.scoringMaterial << "\n";
    out << "world_half_length_cm = " << config.worldHalfLengthCm << "\n";
    out << "scoring_half_length_x_mm = " << config.scoringHalfLengthXmm << "\n";
    out << "scoring_half_length_y_mm = " << config.scoringHalfLengthYmm << "\n";
    out << "scoring_half_length_z_mm = " << config.scoringHalfLengthZmm << "\n";
    out << "scoring_radius_mm = " << config.scoringRadiusMm << "\n";
    out << "scoring_trans_x_mm = " << config.scoringTransXmm << "\n";
    out << "scoring_trans_y_mm = " << config.scoringTransYmm << "\n";
    out << "scoring_trans_z_mm = " << config.scoringTransZmm << "\n";
    out << "domain_radius_um = " << config.domainRadiusUm << "\n";
    out << "nucleus_radius_um = " << config.nucleusRadiusUm << "\n";
    out << "beta_ref_per_gy2 = " << config.betaRefPerGy2 << "\n";
    out << "electron_cut_m = " << config.electronRangeCutM << "\n";
    out << "stopping_power = "
        << toTopasStoppingPowerModeName(config.stoppingPowerMode) << "\n";
    out << "step_calculator = "
        << toTopasStepCalculatorModeName(config.stepCalculatorMode) << "\n";
    out << "phase_space_pre_check = "
        << (config.phaseSpacePreCheck ? "True" : "False") << "\n";
}

std::string shellQuote(const std::string& value) {
    std::string quoted = "'";
    for (const char c : value) {
        if (c == '\'') {
            quoted += "'\\''";
        } else {
            quoted += c;
        }
    }
    quoted += "'";
    return quoted;
}

std::string buildTopasCommand(const AmfConfig& config,
                              const AmfStagedRun& stagedRun,
                              const fs::path& stdoutLog,
                              const fs::path& stderrLog) {
    return "cd " + shellQuote(stagedRun.runDirectory.string()) + " && " +
           shellQuote(config.topasExecutable.string()) + " " +
           shellQuote(stagedRun.parameterFile.filename().string()) +
           " > " + shellQuote(stdoutLog.filename().string()) +
           " 2> " + shellQuote(stderrLog.filename().string());
}

int normalizeSystemExitCode(int systemStatus) {
    if (systemStatus == -1) {
        return -1;
    }

    if (WIFEXITED(systemStatus)) {
        return WEXITSTATUS(systemStatus);
    }

    return systemStatus;
}

}  // namespace

AmfStagedRun AmfRunner::stageRun(const AmfConfig& config) {
    const PhaseSpacePair phaseSpacePair =
        validatePhaseSpacePair(config.phaseSpaceBasePath);
    requireRegularFile(config.tsedPath, "AMF tsed.dat");

    const fs::path runDirectory = requireStagedRunDirectory(config);
    const fs::path stagedTsedPath = runDirectory / "tsed.dat";
    const fs::path stagedPhaseSpacePath =
        runDirectory / phaseSpacePair.phaseSpacePath.filename();
    const fs::path stagedHeaderPath =
        runDirectory / phaseSpacePair.headerPath.filename();
    const fs::path parameterFile = runDirectory / "replay_amf.txt";
    const fs::path manifestFile = runDirectory / "amf_run_manifest.txt";

    fs::copy_file(config.tsedPath,
                  stagedTsedPath,
                  fs::copy_options::overwrite_existing);
    fs::copy_file(phaseSpacePair.phaseSpacePath,
                  stagedPhaseSpacePath,
                  fs::copy_options::overwrite_existing);
    fs::copy_file(phaseSpacePair.headerPath,
                  stagedHeaderPath,
                  fs::copy_options::overwrite_existing);

    AmfTemplateWriter::writeReplayParameterFile(config, parameterFile);

    AmfStagedRun stagedRun{
        runDirectory,
        parameterFile,
        manifestFile,
        stagedTsedPath,
        stagedPhaseSpacePath,
        stagedHeaderPath
    };

    writeManifest(config, manifestFile, stagedRun);

    return stagedRun;
}

AmfRunResult AmfRunner::runTopas(const AmfConfig& config) {
    const AmfStagedRun stagedRun = stageRun(config);
    const fs::path stdoutLog = stagedRun.runDirectory / "topas_stdout.log";
    const fs::path stderrLog = stagedRun.runDirectory / "topas_stderr.log";

    const std::string command =
        buildTopasCommand(config, stagedRun, stdoutLog, stderrLog);
    const int exitCode = normalizeSystemExitCode(std::system(command.c_str()));

    return AmfRunResult{
        stagedRun,
        exitCode,
        stdoutLog,
        stderrLog
    };
}
