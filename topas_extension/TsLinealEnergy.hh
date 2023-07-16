// Claculate GSM2
// Author: Giorgio Cartechini
// Date: 13/11/2021

#ifndef TsLinealEnergy_hh
#define TsLinealEnergy_hh

#include <stdint.h>
#include <vector>
#include <string>
#include<random>

//#include "globals.hh"
//#include "G4RandomDirection.hh"
//#include "G4SystemOfUnits.hh"
//#include "g4root.hh"
using namespace std;

class TsLinealEnergy 
{
	public:
		TsLinealEnergy(std::vector<double> yVector, std::vector<std::vector<double>>);
		~TsLinealEnergy();
	
		void InitializeMicrodosimetricSpectrum();
		void InitializeStatistic();
		void GetSpectrum();
		void GetErrorPropagation();
		void Calculatefy(std::vector<double> dataVector);
		void GetStatisticInfo(int Binindex, double variable);

		double GetyF() {return yF;};
		double GetyD() {return yD;};
		double GetyFvar() {return yF_var;};
		double GetyDvar() {return yD_var;};
		
		std::vector<double> GetfyVariance() {return fVariance;};
		std::vector<double> Getfy() {return hfy;};
		std::vector<double> Getyfy() {return hyfy;};
		std::vector<double> GetyfyVariance() {return yfy_var;};	
		std::vector<double> Getdy() {return hdy;};
		std::vector<double> GetdyVariance() {return dy_var;};
		std::vector<double> Getydy() {return hydy;};
		std::vector<double> GetydyVariance() {return ydy_var;};

		std::vector<double> GetyBinWidth() {return BinWidth;};
		std::vector<double> GetyBinLimit() {return BinLimit;};

		std::vector<std::vector<double>> GetParticleContribution() {return yParticleContibution;};
	private:

	std::vector<double> fyVector;
	std::vector<std::vector<double>> fyVector_Particle;
	const int    yBinNum = 100; // yBinNum == yBinMagnitude*yBinMagnitudeInterval
	const double yBinMagnitudeInterval = 20.;
	const double yBinMagnitude = 5.;
	
	// statisic
	std::vector<double> fFirstMomentMap;
	std::vector<double> fSecondMomentMap;
	std::vector<double> fCountMap;
	std::vector<double> fVariance;
	std::vector<double> fStandardDeviation;
	int fSpectrumUpdateTimes;
	double yF, yD;
	double yF_var, yF_std;
	double yD_var, yD_std;
	std::vector<double> ydy_var, ydy_std;
	std::vector<double> yfy_var, yfy_std;
	std::vector<double>  dy_var,  dy_std;

	//std::vector<std::vector<double>> hfy_particle;
	std::vector<std::vector<double>> yParticleContibution;
	double **hfy_particle;
	std::vector<double> hfy, hdy, hyfy, hydy, BinLimit, BinWidth;
	

};
#endif
