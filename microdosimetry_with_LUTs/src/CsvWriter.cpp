#include "CsvWriter.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace {

std::filesystem::path letSummaryPath(const std::filesystem::path& outputDir,
                                     const std::string& element) {
    return outputDir / ("let_summary_" + element + ".csv");
}

}  // namespace

void CsvWriter::writeProtonMatches(
    const std::filesystem::path& outPath,
    const std::vector<ProtonMatchRecord>& matches) {
    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open output CSV for writing: " + outPath.string());
    }

    out << "row_index,energy_mev,weight,"
        << "lower_matched_energy_mev,upper_matched_energy_mev,"
        << "lower_interpolation_weight,upper_interpolation_weight,"
        << "lower_matched_file,upper_matched_file,"
        << "lower_matched_ncpp,upper_matched_ncpp,"
        << "lower_matched_family,upper_matched_family\n";

    for (const auto& m : matches) {
        out << m.rowIndex << ','
            << m.energyMeV << ','
            << m.weight << ','
            << m.lowerMatchedEnergyMeV << ','
            << m.upperMatchedEnergyMeV << ','
            << m.lowerInterpolationWeight << ','
            << m.upperInterpolationWeight << ','
            << '"' << m.lowerMatchedFile << '"' << ','
            << '"' << m.upperMatchedFile << '"' << ','
            << m.lowerMatchedNcpp << ','
            << m.upperMatchedNcpp << ','
            << '"' << m.lowerMatchedFamily << '"' << ','
            << '"' << m.upperMatchedFamily << '"'
            << '\n';
    }
}

void CsvWriter::writePolySpectrum(const std::filesystem::path& outPath,
                                  const PolySpectrum& spectrum) {
    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open output CSV for writing: " + outPath.string());
    }

    out << "y_keV_per_um,f_y,yf_y,d_y,yd_y\n";
    for (std::size_t i = 0; i < spectrum.y.size(); ++i) {
        out << spectrum.y[i] << ','
            << spectrum.f[i] << ','
            << spectrum.yf[i] << ','
            << spectrum.d[i] << ','
            << spectrum.yd[i] << '\n';
    }
}

void CsvWriter::writeLetSummary(
    const std::filesystem::path& outPath,
    const std::vector<LetSummaryRecord>& records) {

    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open LET summary CSV for writing: " + outPath.string());
    }

    out << "group,track_averaged_LET_keV_per_um,"
        << "dose_averaged_LET_keV_per_um\n";

    for (const auto& record : records) {
        out << record.group << ','
            << record.trackAveragedLetKeVPerUm << ','
            << record.doseAveragedLetKeVPerUm << '\n';
    }
}

void CsvWriter::writeElementLetSummaries(
    const std::filesystem::path& outputDir,
    const std::vector<ElementLetSummary>& summaries) {

    std::filesystem::create_directories(outputDir);

    for (const auto& summary : summaries) {
        writeLetSummary(letSummaryPath(outputDir, summary.element),
                        summary.records);
    }
}

void CsvWriter::writeInaniwaSummary(
    const std::filesystem::path& outPath,
    const std::vector<InaniwaSummaryRecord>& summaries) {

    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open Inaniwa summary CSV for writing: " + outPath.string());
    }

    out << "atomic_number,z_d_D_mean_Gy,z_d_D_star_mean_Gy,z_n_D_mean_Gy\n";

    for (const auto& summary : summaries) {
        out << summary.atomicNumber << ','
            << summary.zdDMeanGy << ','
            << summary.zdDStarMeanGy << ','
            << summary.znDMeanGy << '\n';
    }
}

void CsvWriter::writeMaginiSummary(const std::filesystem::path& outPath,
                                   const MaginiSummary& summary) {
    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open Magini summary CSV for writing: " + outPath.string());
    }

    out << "proton_count,total_proton_weight,"
        << "y_F_keV_per_um,y_D_keV_per_um,y_star_keV_per_um,"
        << "y_D_weighted_numerator,y_star_weighted_numerator\n";

    out << summary.protonCount << ','
        << summary.totalProtonWeight << ','
        << summary.yFKeVPerUm << ','
        << summary.yDKeVPerUm << ','
        << summary.yStarKeVPerUm << ','
        << summary.yDWeightedNumerator << ','
        << summary.yStarWeightedNumerator << '\n';
}
