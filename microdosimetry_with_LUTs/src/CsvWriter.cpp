#include "CsvWriter.h"

#include <fstream>
#include <stdexcept>

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
