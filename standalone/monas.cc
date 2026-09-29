// *************************************************************************
// * MONAS is a C++ package that calculates cell surviavl curvs and        *
// * dose dependednt RBE from microdosimetric spectra.			   *
// *									   *
// * Copyright © 2023 Giorgio Cartechini <giorgio.cartechini@maastro.nl>	   *
// * 									   *
// * This program is free software: you can redistribute it and/or modify  *
// * it under the terms of the GNU General Public License as published by  *
// * the Free Software Foundation, either version 3 of the License, or     *
// * (at your option) any later version.			           *
// * 									   *
// * This program is distributed in the hope that it will be useful,       *
// * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
// * GNU General Public License for more details.			   *
// * 									   *
// * You should have received a copy of the GNU General Public license     *
// * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
// **************************************************************************
#include<iostream>
#include<vector>
#include<cstring>
#include<cmath>
#include<fstream>
#include<algorithm>
#include<random>
#include<sstream>
#include<stdexcept>
#include<cerrno>
#include<sys/stat.h>
#include<sys/types.h>
//#include<filesystem>
#include <chrono>
#include <ctime>
#include "TsGetSurvivalRBEQualityFactor.hh"
#include "TsBinnedSpectrumReader.hh"
#include "TsLinealEnergy.hh"
#include "TsSpecificEnergy.hh"
using namespace std;

namespace {

struct SpectrumInputSpec {
	string SourceName;
	string Type;
	vector<string> Arguments;
};

vector<string> SplitString(const string& text, char delimiter)
{
	vector<string> parts;
	string part;
	stringstream stream(text);
	while(getline(stream, part, delimiter))
		parts.push_back(part);
	return parts;
}

string SanitizeForFilename(const string& text)
{
	string output;
	for(char c:text)
	{
		if((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '-')
			output.push_back(c);
		else
			output.push_back('_');
	}
	return output.empty() ? "Spectrum" : output;
}

bool DirectoryExists(const string& path)
{
	struct stat info;
	return stat(path.c_str(), &info) == 0 && (info.st_mode & S_IFDIR);
}

void EnsureDirectory(const string& path)
{
	if(path.empty() || path == ".")
		return;
	if(DirectoryExists(path))
		return;

	string current;
	for(std::size_t i = 0; i < path.size(); ++i)
	{
		current.push_back(path[i]);
		if(path[i] != '/' && i + 1 != path.size())
			continue;
		if(current.empty() || current == "/")
			continue;
		if(!DirectoryExists(current))
		{
			if(mkdir(current.c_str(), 0755) != 0 && errno != EEXIST)
				throw runtime_error("Cannot create output directory: " + current);
		}
	}
}

string JoinPath(const string& directory, const string& filename)
{
	if(directory.empty() || directory == ".")
		return filename;
	const char lastChar = directory[directory.size() - 1];
	if(lastChar == '/' || lastChar == '\\')
		return directory + filename;
	return directory + "/" + filename;
}

SpectrumInputSpec ParseSpectrumSpec(const string& spec)
{
	const vector<string> fields = SplitString(spec, ':');
	if(fields.size() < 3)
		throw runtime_error("Invalid -spectrum spec '" + spec + "'. Expected name:type:paths...");

	SpectrumInputSpec parsed;
	parsed.SourceName = fields[0];
	parsed.Type = fields[1];
	for(std::size_t i = 2; i < fields.size(); ++i)
		parsed.Arguments.push_back(fields[i]);
	return parsed;
}

TsBinnedSpectrum ReadSpectrumFromSpec(const SpectrumInputSpec& spec)
{
	if(spec.Type == "poly")
	{
		if(spec.Arguments.size() != 1)
			throw runtime_error("Poly spectrum spec requires name:poly:poly_spectrum.csv");
		return TsBinnedSpectrumReader::ReadPolySpectrumCsv(spec.Arguments[0], spec.SourceName);
	}
	if(spec.Type == "amf")
	{
		if(spec.Arguments.size() != 3)
			throw runtime_error("AMF spectrum spec requires name:amf:spectra.csv:moments.csv:x,y,z");
		return TsBinnedSpectrumReader::ReadAmfSpectrumCsv(spec.Arguments[0], spec.Arguments[1], spec.Arguments[2], spec.SourceName);
	}
	if(spec.Type == "topas")
	{
		if(spec.Arguments.size() != 1)
			throw runtime_error("TOPAS spectrum spec requires name:topas:ySpecfile.txt");
		return TsBinnedSpectrumReader::ReadTopasYSpecfile(spec.Arguments[0], spec.SourceName);
	}

	throw runtime_error("Unknown spectrum type '" + spec.Type + "' for source '" + spec.SourceName + "'.");
}

double FiniteStdDev(double variance)
{
	const double stddev = sqrt(variance);
	return std::isfinite(stddev) ? stddev : 0.0;
}

void ApplyMKMParameters(TsGetSurvivalRBEQualityFactor& calculator,
                        vector<double>& Doses,
                        double MKModel_alpha0,
                        double MKModel_beta,
                        double MKModel_alphaX,
                        double MKModel_betaX,
                        double MKModel_rd,
                        double MKModel_Rn,
                        double MKModel_y0,
                        int fSetMultiEventStatistic)
{
	calculator.SetDosesMacro(&Doses[0]);
	calculator.SetMCMultieventIterations(fSetMultiEventStatistic);
	calculator.SetMKModel_alpha0(MKModel_alpha0);
	calculator.SetMKModel_beta(MKModel_beta);
	calculator.SetMKModel_alphaX(MKModel_alphaX);
	calculator.SetMKModel_betaX(MKModel_betaX);
	calculator.SetMKModel_rd(MKModel_rd);
	calculator.SetMKModel_Rn(MKModel_Rn);
	calculator.SetMKModel_y0(MKModel_y0);
}

void AppendSummaryRows(std::ofstream& summary,
                       const TsBinnedSpectrum& spectrum,
                       const TsGetSurvivalRBEQualityFactor& calculator)
{
	const vector<double>& doses = calculator.GetLastDoses();
	const vector<double>& survival = calculator.GetLastSurvival();
	const vector<double>& survivalVariance = calculator.GetLastSurvivalVariance();
	const vector<double>& rbe = calculator.GetLastRBE();
	const vector<double>& rbeVariance = calculator.GetLastRBEVariance();

	for(std::size_t i = 0; i < doses.size(); ++i)
	{
		summary << spectrum.SourceName << ','
		        << calculator.GetLastModelName() << ','
		        << doses[i] << ','
		        << survival[i] << ','
		        << FiniteStdDev(survivalVariance[i]) << ','
		        << rbe[i] << ','
		        << FiniteStdDev(rbeVariance[i]) << ','
		        << spectrum.yF << ','
		        << spectrum.yD << '\n';
	}
}

void WriteBinnedModeManifest(const string& outputDirectory,
                             int argc,
                             char* argv[],
                             const vector<string>& spectrumSpecs,
                             bool MKMSatCorrFlag,
                             bool MKMnonPoissFlag,
                             bool SMKMFlag,
                             bool GSM2Flag,
                             const vector<double>& Doses,
                             double MKModel_alpha0,
                             double MKModel_beta,
                             double MKModel_alphaX,
                             double MKModel_betaX,
                             double MKModel_rd,
                             double MKModel_Rn,
                             double MKModel_y0)
{
	const string manifestPath = JoinPath(outputDirectory, "MKM_from_spectra_manifest.txt");
	ofstream manifest(manifestPath.c_str());
	if(!manifest)
		throw runtime_error("Cannot write manifest file: " + manifestPath);

	const time_t now = time(0);
	manifest << "MONAS deterministic binned-spectrum MKM run\n";
	manifest << "Run timestamp: " << ctime(&now);
	manifest << "Command:";
	for(int i = 0; i < argc; ++i)
		manifest << ' ' << argv[i];
	manifest << "\n\n";

	manifest << "Output directory: " << outputDirectory << "\n\n";
	manifest << "Input spectra:\n";
	for(const string& spec:spectrumSpecs)
		manifest << "- " << spec << "\n";

	manifest << "\nModels:\n";
	if(MKMSatCorrFlag)
		manifest << "- MKM_SaturationCorrected\n";
	if(MKMnonPoissFlag)
		manifest << "- MKM_nonPoisson\n";
	if(SMKMFlag)
		manifest << "- SMKM\n";
	if(GSM2Flag)
		manifest << "- GSM2\n";

	manifest << "\nDose macro: " << Doses[0] << " " << Doses[1] << " " << Doses[2] << "\n";
	manifest << "\nMKM parameters:\n";
	manifest << "MKM_alpha0 = " << MKModel_alpha0 << " Gy-1\n";
	manifest << "MKM_beta = " << MKModel_beta << " Gy-2\n";
	manifest << "MKM_alphaX = " << MKModel_alphaX << " Gy-1\n";
	manifest << "MKM_betaX = " << MKModel_betaX << " Gy-2\n";
	manifest << "MKM_rDomain = " << MKModel_rd << " um\n";
	manifest << "MKM_rNucleus = " << MKModel_Rn << " um\n";
	manifest << "MKM_y0 = " << MKModel_y0 << " keV/um\n";
}

} // namespace

int main(int argc, char *argv[])
{

	/////////////////////////////////////////////////////////////	
	// TIME VARIABLES	////////////////////////////////////////
	/////////////////////////////////////////////////////////////
	auto start = std::chrono::high_resolution_clock::now();
	/////////////////////////////////////////////////////////////
	//
	// INITIALIZE VARIABLES
	//
	////////////////////////////////////////////////////////////
	bool fGetStatitisticInfo = false;
	int fSpectrumUpdateTimes = 100;
	int fSetMultiEventStatistic = 100000;
	bool fGetParticleContribution = false;
	//Initialize Microdosimetric variables
	double yF=0., yD=0.;
	double yF_var=0., yF_std=0.;
	double yD_var=0., yD_std=0.;
	vector<double> hfy, hdy, hyfy, hydy, BinLimit, BinWidth;
	std::vector<double> fy_var, ydy_var, yfy_var, dy_var;
	double MKModel_alpha0 = 0.13; // Unit:Gy-1
	double MKModel_beta   = 0.05; // Unit:Gy-2,
	double MKModel_alphaX = 0.19; // Unit:Gy-1
	double MKModel_betaX = 0.05;
	double MKModel_rd     = 0.42; // Unit:um
	double MKModel_Rn     = 6; //Unit:um
	double MKModel_rho    = 1.;    // Unit:g/cm3
	double MKModel_y0     = 150;  // Unit:keV/um

	//Split Dose MKM
	double MKModel_D1 = 1; //Unit:Gy
	double MKModel_D2 = 1; //Unit:Gy
	double MKModel_ac = 2.187; //Unit:h-1
	double MKModel_tr = 2.284; //Unit:h

	//GSM2 default parameters
	double GSM2_kappa      = 0.5;
	double GSM2_lambda     = 0.5;
	double GSM2_rd         = 0.42; //Unit: um
	double GSM2_Rn         = 6;    //Unit: um
	double GSM2_zStart     = 0.01; //Gy
	double GSM2_zEnd       = 100; //Gy
	double GSM2_zBins      = 100; //Bins of spectra
	double GSM2_a          = 0.1;
	double GSM2_b          = 0.1;
	double GSM2_r          = 0.1;
	double GSM2_alphaX 	= 0.19;
	double GSM2_betaX 	= 0.05;
	string GSM2_ion        = "H";
	double GSM2_LET        = 10.0; //keV/um
	//Macroscopic Doses
	vector<double> Doses = {0,10,0.5}; //Unit:Gy
	bool MKMSatCorrFlag = 0, MKMnonPoissFlag = 0, SMKMFlag = 0, DSMKMFlag = 0, GSM2Flag = 0, RBEWeightingFlag = 0, QfICRUFlag = 0, QfKellFlag = 0, UseTwoSpecraFlag =0;
	bool H460Flag = 0, H1437Flag = 0;
	bool MKMFromSpectraFlag = false;
	string OutputDirectory = ".";
	vector<string> SpectrumSpecs;
	/////////////////////////////////////////////////////////////
	//
	// READ INPUT PARAMETERS
	//
	////////////////////////////////////////////////////////////

	string TopasScorerFileDomain;
	string TopasScorerFileNucleus; // new line
	string BioWeightFunctionDataFile = "BioWeightFuncData_interpolation.txt";
	//Loop on inputs
	for(int i=0; i<argc; i++)
	{
		if(strcmp(argv[i],"-MKM_rDomain") == 0) {MKModel_rd = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_rNucleus") == 0) {MKModel_Rn = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_alpha") == 0) {MKModel_alpha0 = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_beta") == 0) {MKModel_beta = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_alphaX") == 0) {MKModel_alphaX = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_betaX") == 0) {MKModel_betaX = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_y0") == 0) {MKModel_y0 = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_D1") == 0) {MKModel_D1 = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_D2") == 0) {MKModel_D2 = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_ac") == 0) {MKModel_ac = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_tr") == 0) {MKModel_tr= stod(argv[i+1]);}

		// PARAMETERS FROM CELL LINE FLAG (TO ADD, ALSO SAME alphaX, betaX FOR BOTH GSM2 AND MKM)
		if(strcmp(argv[i],"-H460") == 0) {H460Flag = 1;}
		if(strcmp(argv[i],"-H1437") == 0) {H1437Flag = 1;}
		
		if(strcmp(argv[i],"-GSM2_rDomain") == 0) {GSM2_rd = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_rNucleus") == 0) {GSM2_Rn = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_kappa") == 0) {GSM2_kappa = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_a") == 0) {GSM2_a = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_b") == 0) {GSM2_b = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_r") == 0) {GSM2_r = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_alphaX") == 0) {GSM2_alphaX = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_betaX") == 0) {GSM2_betaX = stod(argv[i+1]);}		
		if(strcmp(argv[i],"-GSM2_ion") == 0) {GSM2_ion = argv[i+1];}
		if(strcmp(argv[i],"-GSM2_LET") == 0) {GSM2_LET = stod(argv[i+1]);}

		if(strcmp(argv[i],"-MKMSatCorr") == 0) {MKMSatCorrFlag = 1;}
		if(strcmp(argv[i],"-MKMnonPoiss") == 0) {MKMnonPoissFlag = 1;}
		if(strcmp(argv[i],"-SMKM") == 0) {SMKMFlag = 1;}
		if(strcmp(argv[i],"-DSMKM") == 0) {DSMKMFlag = 1;}
		if(strcmp(argv[i],"-GSM2") == 0) {GSM2Flag = 1;}
		if(strcmp(argv[i],"-MKMFromSpectra") == 0) {MKMFromSpectraFlag = true;}

		if(strcmp(argv[i],"-fGetStatitisticInfo") == 0) {fGetStatitisticInfo = true;}
		if(strcmp(argv[i],"-fSpectrumUpdateTimes") == 0) {fSpectrumUpdateTimes = stoi(argv[i+1]);}
		if(strcmp(argv[i],"-fSetMultiEventStatistic") == 0) {fSetMultiEventStatistic = stoi(argv[i+1]);}
		if(strcmp(argv[i],"-fGetParticleContribution") == 0) {fGetParticleContribution = true;}
		if(strcmp(argv[i], "-RBEWeighting") == 0) 
		{	
			BioWeightFunctionDataFile = argv[i+1];
			RBEWeightingFlag = 1;
		}
		if(strcmp(argv[i], "-QfICRU") == 0) {QfICRUFlag = 1;}
		if(strcmp(argv[i], "-QfKeller") == 0) {QfKellFlag = 1;}
		if(strcmp(argv[i],"-Doses") == 0) {
			Doses.clear();
			Doses = {stod(argv[i+1]), stod(argv[i+2]), stod(argv[i+3])};
		}
		if(strcmp(argv[i],"-outputDir") == 0) {OutputDirectory = argv[i+1];}
		if(strcmp(argv[i],"-spectrum") == 0) {SpectrumSpecs.push_back(argv[i+1]);}

		if(strcmp(argv[i],"-topasScorerDomain") == 0) {TopasScorerFileDomain = argv[i+1];} //input file
		if(strcmp(argv[i],"-topasScorerNucleus") == 0) 
		{
			UseTwoSpecraFlag = true;
			TopasScorerFileNucleus = argv[i+1];
		} //input file // new line
		if(strcmp(argv[i],"-help") == 0) 
		{
			cout 	<<"-Rd: Domain radius [um]" <<endl
				<<"-Rc: Cell Nucleus radius [um]" <<endl
				<<"-CellLine: replace with cell line name (-H460 or -H1437)" <<endl
				<<"-topasScorer: path/file.phsp with y values" <<endl
				<<"-MKMFromSpectra: run deterministic binned-spectrum MKM mode" <<endl
				<<"-spectrum: repeatable deterministic spectrum input:" <<endl
				<<"           name:poly:file.csv where CSV columns are y_keV_per_um,f_y,yf_y,d_y,yd_y" <<endl
				<<"           name:topas:ySpecfile.txt" <<endl
				<<"           name:amf:spectra.csv:moments.csv:x,y,z" <<endl
				<<"-outputDir: directory for output files; deterministic-spectrum mode also writes MKM_from_spectra_summary.csv and MKM_from_spectra_manifest.txt" <<endl
				<<"-MKMSatCorr/-MKMnonPoiss/-SMKM/-GSM2: deterministic-spectrum radiobiology models; defaults to -MKMSatCorr when none is selected" <<endl
				<<"-DSMKM: not available with -MKMFromSpectra" <<endl
				<<"-Doses: (initial dose value) (final dose value) (step value)" <<endl
				<<"-help: list of definitions and inputs" <<endl;
			return 0;			     
		}
	}
	
	// Parameters for each cell line (H460, H1437)
	if(H460Flag)
	{
		cout 	<<"\nH460 CELL LINE IS SELECTED \n" <<endl;
		GSM2_alphaX = 0.29;
		GSM2_betaX = 0.083;
		MKModel_alphaX = 0.29;
		MKModel_betaX = 0.083;
		GSM2_rd = 0.8;
		GSM2_Rn = 6;
		GSM2_a = 0.000899;
		GSM2_b = 0.0642;
		GSM2_r = 2.71;
	}else if(H1437Flag){
		cout 	<<"\nH1437 CELL LINE IS SELECTED \n" <<endl;
		GSM2_alphaX = 0.05;
		GSM2_betaX = 0.041;
		MKModel_alphaX = 0.05;
		MKModel_betaX = 0.041;
		GSM2_rd = 0.6;
		GSM2_Rn = 8;
		GSM2_a = 0.0148;
		GSM2_b = 0.0149;
		GSM2_r = 2.70;
	}else{
		// Default values
		cout << "\033[1;33m\nWARNING: NO SPECIFIC CELL LINE IS SELECTED (USING DEFAULT VALUES)\033[0m\n" << endl;
	}

	if(MKMFromSpectraFlag)
	{
		try
		{
			if(SpectrumSpecs.empty())
				throw runtime_error("-MKMFromSpectra requires at least one -spectrum argument.");
			if(DSMKMFlag)
				throw runtime_error("DSMKM is not available in deterministic binned-spectrum mode yet.");
			if(!MKMSatCorrFlag && !MKMnonPoissFlag && !SMKMFlag && !GSM2Flag)
				MKMSatCorrFlag = true;

			EnsureDirectory(OutputDirectory);
			const string summaryPath = JoinPath(OutputDirectory, "MKM_from_spectra_summary.csv");
			ofstream summary(summaryPath.c_str());
			if(!summary)
				throw runtime_error("Cannot write summary file: " + summaryPath);
			summary << "source,model,dose_Gy,survival,survival_std,RBE,RBE_std,yF_keV_per_um,yD_keV_per_um\n";
			WriteBinnedModeManifest(OutputDirectory, argc, argv, SpectrumSpecs, MKMSatCorrFlag, MKMnonPoissFlag, SMKMFlag, GSM2Flag, Doses, MKModel_alpha0, MKModel_beta, MKModel_alphaX, MKModel_betaX, MKModel_rd, MKModel_Rn, MKModel_y0);

			cout << "\nRunning deterministic binned-spectrum MKM mode\n";
			cout << "Output directory: " << OutputDirectory << endl;
			cout << "Summary file: " << summaryPath << endl;

			for(const string& specText:SpectrumSpecs)
			{
				const SpectrumInputSpec spec = ParseSpectrumSpec(specText);
				TsBinnedSpectrum spectrum = ReadSpectrumFromSpec(spec);

				cout << "\nSpectrum source: " << spectrum.SourceName << endl;
				cout << "yF: " << spectrum.yF << " keV/um" << endl;
				cout << "yD: " << spectrum.yD << " keV/um" << endl;

				TsGetSurvivalRBEQualityFactor calculator(spectrum);
				ApplyMKMParameters(calculator, Doses, MKModel_alpha0, MKModel_beta, MKModel_alphaX, MKModel_betaX, MKModel_rd, MKModel_Rn, MKModel_y0, fSetMultiEventStatistic);
				calculator.SetOutputDirectory(OutputDirectory);
				calculator.SetOutputPrefix(SanitizeForFilename(spectrum.SourceName) + "_");

				if(MKMSatCorrFlag)
				{
					calculator.GetSurvWithMKModel_SaturationCorr();
					AppendSummaryRows(summary, spectrum, calculator);
				}
				if(MKMnonPoissFlag)
				{
					calculator.GetSurvWithMKModel_nonPoissonCorr();
					AppendSummaryRows(summary, spectrum, calculator);
				}
				if(SMKMFlag)
				{
					calculator.GetSurvWithSMKModel();
					AppendSummaryRows(summary, spectrum, calculator);
				}
				if(GSM2Flag)
				{
					calculator.SetGSM2_alphaX(GSM2_alphaX);
					calculator.SetGSM2_betaX(GSM2_betaX);
					calculator.SetGSM2_rd(GSM2_rd);
					calculator.SetGSM2_Rn(GSM2_Rn);
					calculator.SetGSM2_a(GSM2_a);
					calculator.SetGSM2_b(GSM2_b);
					calculator.SetGSM2_r(GSM2_r);
					calculator.SetGSM2_ion(GSM2_ion);
					calculator.SetGSM2_LET(GSM2_LET);
					calculator.GetSurvWithGSM2();
					AppendSummaryRows(summary, spectrum, calculator);
				}
			}
		}
		catch(const exception& error)
		{
			cerr << "\033[1;31mERROR:: " << error.what() << "\033[0m" << endl;
			return -1;
		}

		return 0;
	}

	/////////////////////////////////////////////////////////////
	//
	// READ TOPAS SCORER DOMAIN: MANDATORY!!!
	//
	////////////////////////////////////////////////////////////

	vector<vector<double>> yVector_Particle;
	vector<double>  yVector;

	ifstream infile(&TopasScorerFileDomain[0]);
	if(!infile) // checks to see if file opened
	{
		cout << "\033[1;31mERROR:: INPUT FILE (DOMAIN) " << TopasScorerFileNucleus << " NOT FOUND!!!\033[0m" << endl;
		return -1; // no point continuing if the file didn't open...
	}
	while(!infile.eof()) // reads file to end of *file*, not line
	{ 

		double y_total, y_z0, y_z1_prim, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_;
		//double y_total, y_z0, y_z1_prim, y_z1_seco, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_;
		infile >> y_total ;
			// >> y_z0
			// >> y_z1_prim
			// >> y_z1_seco
			// >> y_z2
			// >> y_z3
			// >> y_z4
			// >> y_z5
			// >> y_z6
			// >> y_z_;

		// vector<double> yParticle {y_z0, y_z1_prim, 0, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_, y_total};
		vector<double> yParticle {0, 0, 0, 0, 0, 0, 0, 0, 0, y_total};
		
		//vector<double> yParticle {y_z0, y_z1_prim, y_z1_seco, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_, y_total};
		
		yVector_Particle.push_back(yParticle);
		yVector.push_back(y_total);  
	}
	infile.close();
	
	
	/////////////////////////////////////////////////////////////
	//
	// READ TOPAS SCORER FOR NUCLEUS // new line
	//
	////////////////////////////////////////////////////////////

	vector<vector<double>> yVector_Particle_Nucleus;
	vector<double>  yVector_Nucleus;

	if (UseTwoSpecraFlag)  // if the second spectrum is given use this
	{
		ifstream infile_Nucleus(&TopasScorerFileNucleus[0]);
		if(!infile_Nucleus) // checks to see if file opended 
		{
			cout << "\033[1;31mERROR:: INPUT FILE (NUCLEUS) " << TopasScorerFileNucleus << " NOT FOUND!!!\033[0m" << endl;
			return -1; // no point continuing if the file didn't open...
		}
		while(!infile_Nucleus.eof()) // reads file to end of *file*, not line
		{ 

			double y_total, y_z0, y_z1_prim, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_;
			//double y_total, y_z0, y_z1_prim, y_z1_seco, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_;
			infile_Nucleus >> y_total 
				>> y_z0
				>> y_z1_prim
				//>> y_z1_seco
				>> y_z2
				>> y_z3
				>> y_z4
				>> y_z5
				>> y_z6
				>> y_z_;

			vector<double> yParticle {y_z0, y_z1_prim, 0, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_, y_total};

			//vector<double> yParticle {y_z0, y_z1_prim, y_z1_seco, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_, y_total};

			yVector_Particle_Nucleus.push_back(yParticle);
			yVector_Nucleus.push_back(y_total);  
		}
		infile_Nucleus.close();
	} 
	else // use the domain one and it will be rescaled afterwords
	{
		yVector_Particle_Nucleus = yVector_Particle;
		yVector_Nucleus = yVector;
		cout << "\033[1;33mWARNING: Using domain spectrum for nucleus calculations. Results may not be accurate.\033[0m" << endl;
	}
    
	auto stop_input = std::chrono::high_resolution_clock::now();


	/////////////////////////////////////////////////////////////
	//
	// GET MICRODOSIMETRIC SPECTRA AND QUANTITIES
	//
	////////////////////////////////////////////////////////////

	//TsLinealEnergy* aLinealEnergy = new TsLinealEnergy(yVector, yVector_Particle);
	TsLinealEnergy* aLinealEnergy = new TsLinealEnergy(yVector_Nucleus, yVector_Particle_Nucleus);
	BinLimit = aLinealEnergy->GetyBinLimit();
	BinWidth = aLinealEnergy->GetyBinWidth();
	hfy = aLinealEnergy -> Getfy();
	hdy = aLinealEnergy -> Getdy();
	yF  = aLinealEnergy -> GetyF();
	yD  = aLinealEnergy -> GetyD();
	fy_var = aLinealEnergy -> GetfyVariance();
	dy_var = aLinealEnergy -> GetdyVariance();
	yF_var = aLinealEnergy -> GetyFvar();
	yD_var = aLinealEnergy -> GetyDvar();
	int yBinNum = BinWidth.size();
	std::cout << "\033[1;32mINPUT DOMAIN SPECTRUM PARAMETERS:\033[0m\n";
	std::cout << "\033[1;32m" << "yF: " << yF << " -+ " << sqrt(yF_var) << "\033[0m" << endl
			  << "\033[1;32m" << "yD: " << yD << " -+ " << sqrt(yD_var) << "\033[0m" << endl;

	vector<vector<double>> contribution = aLinealEnergy -> GetParticleContribution();
	hydy = aLinealEnergy->Getydy();

	TsGetSurvivalRBEQualityFactor *aSurvRBEQf = new TsGetSurvivalRBEQualityFactor(contribution, yVector, yVector_Particle, yVector_Nucleus, yVector_Particle_Nucleus, &BinLimit[0], &BinWidth[0], &hfy[0], &hdy[0], yF, yD, yF_var, yD_var, fy_var, dy_var, yBinNum, fGetStatitisticInfo, fSpectrumUpdateTimes, fGetParticleContribution);
	
	aSurvRBEQf->SetDosesMacro(&Doses[0]);
	aSurvRBEQf -> SetMCMultieventIterations(fSetMultiEventStatistic);
	aSurvRBEQf->SetMKModel_alpha0(MKModel_alpha0);
	aSurvRBEQf->SetMKModel_beta(MKModel_beta);
	aSurvRBEQf->SetMKModel_alphaX(MKModel_alphaX);
	aSurvRBEQf->SetMKModel_betaX(MKModel_betaX);
	aSurvRBEQf->SetMKModel_rd(MKModel_rd);
	aSurvRBEQf->SetMKModel_Rn(MKModel_Rn);

	TsSpecificEnergy* aSpecEne = new TsSpecificEnergy(yVector_Particle, MKModel_rd, fGetStatitisticInfo, fSpectrumUpdateTimes);
	std::vector<double> fz = aSpecEne -> GetHfz();

	//std::vector<std::vector<double>> fzParticle = aSpecEne -> GetParticleContribution();
	//std::vector<double> fz1 = aSpecEne ->GetHfzMultiEvent(0.5, 1e6);
	//std::vector<double> fz10 = aSpecEne ->GetHfzMultiEvent(1, 1e6);
	//std::vector<std::vector<double>> fz10Particle = aSpecEne -> GetMultiEventParticleContribution();
	//std::vector<double> fz100 = aSpecEne ->GetHfzMultiEvent(15, 1e6);
	
	/*for(int i=0; i<hfy.size(); i++)
	{
		cout << zbin[i] <<',' <<hfy[i]; //<<',' <<fz1[i]; //<<','<<fz10[i]<<','<<fz100[i];
		//for(int j=0; j<10; j++)
		//	cout <<','<<fz10Particle[i][j]*fz10[i]; //<<',' <<fz1[i]<<','<<fz10[i]<<','<<fz100[i] <<endl;
		cout << endl;
	}
	std::vector<double> zbin = aSpecEne -> GetBinCenter();
	for(int i=0; i<fz.size(); i++)
	{
		cout << zbin[i] <<',' <<fz[i]; //<<',' <<fz1[i]; //<<','<<fz10[i]<<','<<fz100[i];
		//for(int j=0; j<10; j++)
		//	cout <<','<<fz10Particle[i][j]*fz10[i]; //<<',' <<fz1[i]<<','<<fz10[i]<<','<<fz100[i] <<endl;
		cout << endl;
	}
	*/
	
	auto GSM2_time = std::chrono::seconds(0);

	if(MKMSatCorrFlag)
		aSurvRBEQf->GetSurvWithMKModel_SaturationCorr();
	if(MKMnonPoissFlag)
		aSurvRBEQf->GetSurvWithMKModel_nonPoissonCorr();
	if(SMKMFlag)
		aSurvRBEQf->GetSurvWithSMKModel();
	if(DSMKMFlag)
		aSurvRBEQf->GetSurvWithDSMKModel();

	if(RBEWeightingFlag)
	{
		aSurvRBEQf -> SetBioWeightFunctionDataFile(BioWeightFunctionDataFile);
		aSurvRBEQf -> GetRBEWithBioWeightFunction();
	}
	if(QfICRUFlag)
		aSurvRBEQf -> GetQualityFactorWithICRU40();
	if(QfKellFlag)
		aSurvRBEQf -> GetQualityFactorWithKellereHahn();
	
	//aSurvRBEQf->GetSurvWithMKModel_SplitDoseIrradiation();
	if(GSM2Flag)
	{
		auto start_surv = std::chrono::high_resolution_clock::now();
		aSurvRBEQf->SetGSM2_alphaX(GSM2_alphaX);
		aSurvRBEQf->SetGSM2_betaX(GSM2_betaX);
		aSurvRBEQf->SetGSM2_rd(GSM2_rd);
		aSurvRBEQf->SetGSM2_Rn(GSM2_Rn);
		aSurvRBEQf->SetGSM2_a(GSM2_a);
		aSurvRBEQf->SetGSM2_b(GSM2_b);
		aSurvRBEQf->SetGSM2_r(GSM2_r);
		aSurvRBEQf->SetGSM2_ion(GSM2_ion);
		aSurvRBEQf->SetGSM2_LET(GSM2_LET);
		
		aSurvRBEQf->GetSurvWithGSM2();
		auto stop = std::chrono::high_resolution_clock::now();
		GSM2_time = std::chrono::duration_cast<std::chrono::seconds>(stop - start_surv);
	}	
	

	auto stop = std::chrono::high_resolution_clock::now();
	auto totaltime = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

	auto input_time = std::chrono::duration_cast<std::chrono::seconds>(stop_input - start);
	auto calculation_time = std::chrono::duration_cast<std::chrono::seconds>(stop - stop_input);

	double input_percentage = (double)input_time.count() / totaltime.count() * 100;
	double calculation_percentage = (double)calculation_time.count() / totaltime.count() * 100;
	double GSM2_percentage = (double)GSM2_time.count() / totaltime.count() * 100;

	cout.precision(0);
	cout << "\nTime Table:\n";
	cout << "-----------------------------------\n";
	cout << "Total Time: " << totaltime.count() << " seconds (" << fixed << totaltime.count() / 60.0 << " minutes)\n";
	cout << "Input Parsing Time: " << input_time.count() << " seconds (" << fixed << input_percentage << "%)\n";
	cout << "Total Calculation Time: " << calculation_time.count() << " seconds (" << fixed << calculation_percentage << "%)\n";
	cout << "GSM2 Calculation Time: " << GSM2_time.count() << " seconds (" << fixed << GSM2_percentage << "%)\n";
	cout << "-----------------------------------\n";

//	aSurvRBEQf->GetSurvWithGSM2();	
	return 0;
}
