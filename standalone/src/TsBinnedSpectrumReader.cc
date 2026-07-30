// *************************************************************************
// * MONAS is a C++ package that calculates cell survival curves and        *
// * dose dependent RBE from microdosimetric spectra.                       *
// *************************************************************************

#include "TsBinnedSpectrumReader.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

std::string Trim(const std::string& input)
{
	const std::string whitespace = " \t\r\n";
	const std::string::size_type first = input.find_first_not_of(whitespace);
	if(first == std::string::npos)
		return "";
	const std::string::size_type last = input.find_last_not_of(whitespace);
	return input.substr(first, last - first + 1);
}

std::vector<std::string> SplitCsvLine(const std::string& line)
{
	std::vector<std::string> fields;
	std::string field;
	bool inQuotes = false;

	for(std::size_t i = 0; i < line.size(); ++i)
	{
		const char c = line[i];
		if(c == '"')
		{
			if(inQuotes && i + 1 < line.size() && line[i + 1] == '"')
			{
				field.push_back('"');
				++i;
			}
			else
			{
				inQuotes = !inQuotes;
			}
		}
		else if(c == ',' && !inQuotes)
		{
			fields.push_back(Trim(field));
			field.clear();
		}
		else
		{
			field.push_back(c);
		}
	}

	fields.push_back(Trim(field));
	return fields;
}

double ParseDouble(const std::string& text, const std::string& context)
{
	try
	{
		std::size_t index = 0;
		const double value = std::stod(Trim(text), &index);
		if(index != Trim(text).size())
			throw std::invalid_argument("trailing text");
		return value;
	}
	catch(const std::exception&)
	{
		throw std::runtime_error("Cannot parse numeric value '" + text + "' in " + context + ".");
	}
}

std::vector<double> InferBinWidthsFromCenters(const std::vector<double>& centers)
{
	if(centers.empty())
		throw std::runtime_error("Cannot infer bin widths from an empty y grid.");

	std::vector<double> widths(centers.size(), 0.0);
	if(centers.size() == 1)
	{
		widths[0] = centers[0];
		return widths;
	}

	for(std::size_t i = 0; i < centers.size(); ++i)
	{
		if(!(centers[i] > 0.0))
			throw std::runtime_error("All y-bin centers must be positive.");
	}

	std::vector<double> edges(centers.size() + 1, 0.0);
	edges[0] = centers[0] / std::sqrt(centers[1] / centers[0]);
	for(std::size_t i = 1; i < centers.size(); ++i)
		edges[i] = std::sqrt(centers[i - 1] * centers[i]);
	edges[centers.size()] = centers.back() * std::sqrt(centers.back() / centers[centers.size() - 2]);

	for(std::size_t i = 0; i < centers.size(); ++i)
		widths[i] = edges[i + 1] - edges[i];

	return widths;
}

void ValidateSameSize(const TsBinnedSpectrum& spectrum)
{
	const std::size_t size = spectrum.YCenter.size();
	if(size == 0)
		throw std::runtime_error("Spectrum '" + spectrum.SourceName + "' has no y bins.");
	if(spectrum.BinWidth.size() != size ||
	   spectrum.FrequencyDensity.size() != size ||
	   spectrum.DoseDensity.size() != size)
	{
		throw std::runtime_error("Spectrum '" + spectrum.SourceName + "' has mismatched vector sizes.");
	}
}

double IntegrateDensity(const std::vector<double>& density,
                        const std::vector<double>& binWidth)
{
	double integral = 0.0;
	for(std::size_t i = 0; i < density.size(); ++i)
		integral += density[i] * binWidth[i];
	return integral;
}

void NormalizeDensity(std::vector<double>& density,
                      const std::vector<double>& binWidth,
                      const std::string& sourceName,
                      const std::string& densityName)
{
	for(std::size_t i = 0; i < density.size(); ++i)
	{
		if(!std::isfinite(density[i]) || density[i] < 0.0)
			throw std::runtime_error("Spectrum '" + sourceName + "' has invalid " + densityName + " values.");
	}

	const double integral = IntegrateDensity(density, binWidth);
	if(!(integral > 0.0) || !std::isfinite(integral))
		throw std::runtime_error("Spectrum '" + sourceName + "' has non-positive " + densityName + " integral.");

	for(std::size_t i = 0; i < density.size(); ++i)
		density[i] /= integral;
}

std::vector<std::string> SplitVoxelId(const std::string& voxelId)
{
	std::vector<std::string> parts;
	std::stringstream stream(voxelId);
	std::string part;
	while(std::getline(stream, part, ','))
		parts.push_back(Trim(part));
	if(parts.size() != 3)
		throw std::runtime_error("AMF voxel id must be formatted as x,y,z.");
	return parts;
}

} // namespace

TsBinnedSpectrum TsBinnedSpectrumReader::ReadPolySpectrumCsv(const std::string& fileName,
                                                             const std::string& sourceName)
{
	std::ifstream input(fileName.c_str());
	if(!input)
		throw std::runtime_error("Cannot open poly spectrum CSV: " + fileName);

	std::string line;
	if(!std::getline(input, line))
		throw std::runtime_error("Poly spectrum CSV is empty: " + fileName);

	const std::vector<std::string> header = SplitCsvLine(line);
	if(header.size() < 5 ||
	   header[0] != "y_keV_per_um" ||
	   header[1] != "f_y" ||
	   header[3] != "d_y")
	{
		throw std::runtime_error("Unexpected poly spectrum CSV header in: " + fileName);
	}

	TsBinnedSpectrum spectrum;
	spectrum.SourceName = sourceName;

	while(std::getline(input, line))
	{
		if(Trim(line).empty())
			continue;
		const std::vector<std::string> fields = SplitCsvLine(line);
		if(fields.size() < 5)
			throw std::runtime_error("Malformed poly spectrum row in: " + fileName);

		spectrum.YCenter.push_back(ParseDouble(fields[0], fileName));
		spectrum.FrequencyDensity.push_back(ParseDouble(fields[1], fileName));
		spectrum.DoseDensity.push_back(ParseDouble(fields[3], fileName));
	}

	FinalizeSpectrum(spectrum);
	return spectrum;
}

TsBinnedSpectrum TsBinnedSpectrumReader::ReadAmfSpectrumCsv(const std::string& spectraFileName,
                                                            const std::string& momentsFileName,
                                                            const std::string& voxelId,
                                                            const std::string& sourceName)
{
	std::ifstream spectraInput(spectraFileName.c_str());
	if(!spectraInput)
		throw std::runtime_error("Cannot open AMF spectra CSV: " + spectraFileName);

	std::string line;
	while(std::getline(spectraInput, line))
	{
		if(!Trim(line).empty() && Trim(line)[0] != '#')
			break;
	}
	if(Trim(line).empty())
		throw std::runtime_error("AMF spectra CSV has no header row: " + spectraFileName);

	const std::vector<std::string> header = SplitCsvLine(line);
	if(header.size() < 4 || header[0] != "x" || header[1] != "y" || header[2] != "z")
		throw std::runtime_error("Unexpected AMF spectra CSV header in: " + spectraFileName);

	TsBinnedSpectrum spectrum;
	spectrum.SourceName = sourceName;
	for(std::size_t i = 3; i < header.size(); ++i)
		spectrum.YCenter.push_back(ParseDouble(header[i], spectraFileName));

	const std::vector<std::string> voxel = SplitVoxelId(voxelId);
	bool foundVoxel = false;
	while(std::getline(spectraInput, line))
	{
		if(Trim(line).empty() || Trim(line)[0] == '#')
			continue;
		const std::vector<std::string> fields = SplitCsvLine(line);
		if(fields.size() != header.size())
			throw std::runtime_error("Malformed AMF spectra row in: " + spectraFileName);

		if(fields[0] == voxel[0] && fields[1] == voxel[1] && fields[2] == voxel[2])
		{
			for(std::size_t i = 3; i < fields.size(); ++i)
			{
				const double y = spectrum.YCenter[i - 3];
				const double yd = ParseDouble(fields[i], spectraFileName);
				spectrum.DoseDensity.push_back(y > 0.0 ? yd / y : 0.0);
			}
			foundVoxel = true;
			break;
		}
	}

	if(!foundVoxel)
		throw std::runtime_error("Voxel " + voxelId + " not found in AMF spectra CSV: " + spectraFileName);

	std::ifstream momentsInput(momentsFileName.c_str());
	if(!momentsInput)
		throw std::runtime_error("Cannot open AMF moments CSV: " + momentsFileName);

	if(!std::getline(momentsInput, line))
		throw std::runtime_error("AMF moments CSV is empty: " + momentsFileName);

	bool foundFrequency = false;
	bool foundDose = false;
	while(std::getline(momentsInput, line))
	{
		if(Trim(line).empty())
			continue;
		const std::vector<std::string> fields = SplitCsvLine(line);
		if(fields.size() < 9)
			throw std::runtime_error("Malformed AMF moments row in: " + momentsFileName);

		if(fields[0] == voxel[0] && fields[1] == voxel[1] && fields[2] == voxel[2])
		{
			if(fields[3] == "frequency")
			{
				spectrum.yF = ParseDouble(fields[4], momentsFileName);
				foundFrequency = true;
			}
			else if(fields[3] == "dose")
			{
				spectrum.yD = ParseDouble(fields[4], momentsFileName);
				foundDose = true;
			}
		}
	}

	if(!foundFrequency || !foundDose)
		throw std::runtime_error("AMF moments for voxel " + voxelId + " not found in: " + momentsFileName);

	spectrum.FrequencyDensity.resize(spectrum.DoseDensity.size(), 0.0);
	for(std::size_t i = 0; i < spectrum.DoseDensity.size(); ++i)
	{
		const double y = spectrum.YCenter[i];
		spectrum.FrequencyDensity[i] = y > 0.0 ? spectrum.DoseDensity[i] * spectrum.yF / y : 0.0;
	}

	FinalizeSpectrum(spectrum);
	return spectrum;
}

TsBinnedSpectrum TsBinnedSpectrumReader::ReadTopasYSpecfile(const std::string& fileName,
                                                            const std::string& sourceName)
{
	std::ifstream input(fileName.c_str());
	if(!input)
		throw std::runtime_error("Cannot open TOPAS y spectrum file: " + fileName);

	TsBinnedSpectrum spectrum;
	spectrum.SourceName = sourceName;

	std::string line;
	while(std::getline(input, line))
	{
		const std::string trimmed = Trim(line);
		if(trimmed.empty())
			continue;
		if(trimmed.find("yF") == 0 || trimmed.find("yD") == 0 ||
		   trimmed.find("y (") == 0 || trimmed[0] == '#')
			continue;

		std::stringstream stream(trimmed);
		double y = 0.0;
		double f = 0.0;
		double fStd = 0.0;
		double yf = 0.0;
		double yfStd = 0.0;
		double d = 0.0;
		double dStd = 0.0;
		double yd = 0.0;
		double ydStd = 0.0;

		if(stream >> y >> f >> fStd >> yf >> yfStd >> d >> dStd >> yd >> ydStd)
		{
			spectrum.YCenter.push_back(y);
			spectrum.FrequencyDensity.push_back(f);
			spectrum.DoseDensity.push_back(d);
		}
	}

	FinalizeSpectrum(spectrum);
	return spectrum;
}

void TsBinnedSpectrumReader::FinalizeSpectrum(TsBinnedSpectrum& spectrum)
{
	if(spectrum.BinWidth.empty())
		spectrum.BinWidth = InferBinWidthsFromCenters(spectrum.YCenter);

	ValidateSameSize(spectrum);
	NormalizeDensity(spectrum.FrequencyDensity, spectrum.BinWidth, spectrum.SourceName, "f(y)");
	NormalizeDensity(spectrum.DoseDensity, spectrum.BinWidth, spectrum.SourceName, "d(y)");

	spectrum.yF = 0.0;
	spectrum.yD = 0.0;
	for(std::size_t i = 0; i < spectrum.YCenter.size(); ++i)
	{
		spectrum.yF += spectrum.YCenter[i] * spectrum.FrequencyDensity[i] * spectrum.BinWidth[i];
		spectrum.yD += spectrum.YCenter[i] * spectrum.DoseDensity[i] * spectrum.BinWidth[i];
	}
}
