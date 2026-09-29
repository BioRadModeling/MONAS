// *************************************************************************
// * MONAS is a C++ package that calculates cell survival curves and        *
// * dose dependent RBE from microdosimetric spectra.                       *
// *************************************************************************

#ifndef TsBinnedSpectrumReader_hh
#define TsBinnedSpectrumReader_hh

#include <string>

#include "TsBinnedSpectrum.hh"

class TsBinnedSpectrumReader {
	public:
		static TsBinnedSpectrum ReadPolySpectrumCsv(const std::string& fileName,
		                                            const std::string& sourceName);
		static TsBinnedSpectrum ReadAmfSpectrumCsv(const std::string& spectraFileName,
		                                          const std::string& momentsFileName,
		                                          const std::string& voxelId,
		                                          const std::string& sourceName);
		static TsBinnedSpectrum ReadTopasYSpecfile(const std::string& fileName,
		                                          const std::string& sourceName);

	private:
		static void FinalizeSpectrum(TsBinnedSpectrum& spectrum);
};

#endif
