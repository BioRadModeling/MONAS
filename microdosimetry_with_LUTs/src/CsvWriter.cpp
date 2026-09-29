#include "CsvWriter.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>

void CsvWriter::writeParticleMatches(
    const std::filesystem::path& outPath,
    const std::vector<ParticleMatchRecord>& matches) {
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

void CsvWriter::writePolySpectrumMoments(
    const std::filesystem::path& outPath,
    const PolySpectrumMomentsSummary& summary) {

    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open output CSV for writing: " + outPath.string());
    }

    out << "distribution,mean_keV_per_um,variance_keV2_per_um2,"
        << "stdev_keV_per_um,mean_standard_error_keV_per_um,skewness\n";

    const auto writeRow = [&out](const SpectrumDistributionMoments& moments) {
        out << moments.distribution << ','
            << moments.meanKeVPerUm << ','
            << moments.varianceKeV2PerUm2 << ','
            << moments.stdevKeVPerUm << ','
            << moments.meanStandardErrorKeVPerUm << ','
            << moments.skewness << '\n';
    };

    writeRow(summary.frequency);
    writeRow(summary.dose);
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

void CsvWriter::writeInaniwaSummary(
    const std::filesystem::path& outPath,
    const InaniwaSummary& summary) {

    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open Inaniwa summary CSV for writing: " + outPath.string());
    }

    out << "z_d_D_mean_Gy,z_d_D_star_mean_Gy,z_n_D_mean_Gy\n";
    out << summary.zdDMeanGy << ','
        << summary.zdDStarMeanGy << ','
        << summary.znDMeanGy << '\n';
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
        << "y_F_lut_error_max_keV_per_um,"
        << "y_D_lut_error_max_keV_per_um,"
        << "y_D_weighted_numerator,y_star_weighted_numerator\n";

    out << summary.protonCount << ','
        << summary.totalProtonWeight << ','
        << summary.yFKeVPerUm << ','
        << summary.yDKeVPerUm << ','
        << summary.yStarKeVPerUm << ','
        << summary.yFMaxErrorKeVPerUm << ','
        << summary.yDMaxErrorKeVPerUm << ','
        << summary.yDWeightedNumerator << ','
        << summary.yStarWeightedNumerator << '\n';
}

void CsvWriter::writeAtSummary(const std::filesystem::path& outPath,
                               const AtSummary& summary) {
    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open AT summary CSV for writing: " + outPath.string());
    }

    out << "selected_particle,"
        << "charged_particle_count,matched_particle_count,skipped_particle_count,"
        << "selected_particle_count,skipped_below_range_count,"
        << "skipped_above_range_count,total_matched_weight,total_frequency_weight,"
        << "y_F_keV_per_um,y_D_keV_per_um,y_star_keV_per_um,"
        << "z_F_Gy,z_D_Gy,z_star_Gy,"
        << "y_F_weighted_numerator,"
        << "y_D_weighted_numerator,y_star_weighted_numerator,"
        << "z_F_weighted_numerator,"
        << "z_D_weighted_numerator,z_star_weighted_numerator\n";

    out << summary.selectedParticle << ','
        << summary.chargedParticleCount << ','
        << summary.matchedParticleCount << ','
        << summary.skippedParticleCount << ','
        << summary.selectedParticleCount << ','
        << summary.skippedBelowRangeCount << ','
        << summary.skippedAboveRangeCount << ','
        << summary.totalMatchedWeight << ','
        << summary.totalFrequencyWeight << ','
        << summary.yFKeVPerUm << ','
        << summary.yDKeVPerUm << ','
        << summary.yStarKeVPerUm << ','
        << summary.zFGy << ','
        << summary.zDGy << ','
        << summary.zStarGy << ','
        << summary.yFWeightedNumerator << ','
        << summary.yDWeightedNumerator << ','
        << summary.yStarWeightedNumerator << ','
        << summary.zFWeightedNumerator << ','
        << summary.zDWeightedNumerator << ','
        << summary.zStarWeightedNumerator << '\n';
}

void CsvWriter::writeAtDiagnostics(const std::filesystem::path& outPath,
                                   const AtSummary& summary) {
    std::ofstream out(outPath);
    if (!out) {
        throw std::runtime_error(
            "Failed to open AT diagnostics CSV for writing: " + outPath.string());
    }

    out << "selected_particle,energy_band_MeV_per_u,row_count,matched_count,"
        << "skipped_below_range_count,skipped_above_range_count,"
        << "total_particle_weight,total_frequency_weight,frequency_weight_fraction,"
        << "y_F_contribution_keV_per_um,y_D_contribution_keV_per_um,"
        << "y_star_contribution_keV_per_um,z_F_contribution_Gy,"
        << "z_D_contribution_Gy,z_star_contribution_Gy\n";

    for (const auto& record : summary.diagnostics) {
        const double frequencyWeightFraction =
            summary.totalFrequencyWeight > 0.0
                ? record.totalFrequencyWeight / summary.totalFrequencyWeight
                : 0.0;
        const double yFContribution =
            summary.totalFrequencyWeight > 0.0
                ? record.yFWeightedNumerator / summary.totalFrequencyWeight
                : 0.0;
        const double yDContribution =
            summary.totalFrequencyWeight > 0.0
                ? record.yDWeightedNumerator / summary.totalFrequencyWeight
                : 0.0;
        const double yStarContribution =
            summary.totalFrequencyWeight > 0.0
                ? record.yStarWeightedNumerator / summary.totalFrequencyWeight
                : 0.0;
        const double zFContribution =
            summary.totalFrequencyWeight > 0.0
                ? record.zFWeightedNumerator / summary.totalFrequencyWeight
                : 0.0;
        const double zDContribution =
            summary.totalFrequencyWeight > 0.0
                ? record.zDWeightedNumerator / summary.totalFrequencyWeight
                : 0.0;
        const double zStarContribution =
            summary.totalFrequencyWeight > 0.0
                ? record.zStarWeightedNumerator / summary.totalFrequencyWeight
                : 0.0;

        out << summary.selectedParticle << ','
            << record.energyBandMeVPerU << ','
            << record.rowCount << ','
            << record.matchedCount << ','
            << record.skippedBelowRangeCount << ','
            << record.skippedAboveRangeCount << ','
            << record.totalParticleWeight << ','
            << record.totalFrequencyWeight << ','
            << frequencyWeightFraction << ','
            << yFContribution << ','
            << yDContribution << ','
            << yStarContribution << ','
            << zFContribution << ','
            << zDContribution << ','
            << zStarContribution << '\n';
    }
}
