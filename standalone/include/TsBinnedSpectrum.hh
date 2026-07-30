// *************************************************************************
// * MONAS is a C++ package that calculates cell survival curves and        *
// * dose dependent RBE from microdosimetric spectra.                       *
// *************************************************************************

#ifndef TsBinnedSpectrum_hh
#define TsBinnedSpectrum_hh

#include <string>
#include <vector>

struct TsBinnedSpectrum {
	std::string SourceName;
	std::vector<double> YCenter;
	std::vector<double> BinWidth;
	std::vector<double> FrequencyDensity;
	std::vector<double> DoseDensity;
	double yF = 0.0;
	double yD = 0.0;
};

#endif
