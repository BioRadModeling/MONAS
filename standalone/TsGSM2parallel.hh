// Claculate GSM2
// Author: Giorgio Cartechini
// Date: 13/11/2021

#ifndef TsGSM2_hh
#define TsGSM2_hh

#include <stdint.h>
#include <vector>
#include <string>
#include <random>
#include <iostream>

#include "TsSpecificEnergy.hh"

//#include "globals.hh"
//#include "G4RandomDirection.hh"
//#include "G4SystemOfUnits.hh"
//#include "g4root.hh"



using namespace std;

//CONTROLLARE GLI HISO (BIN START, BIN END)
class TsGSM2 
{
	public:
		TsGSM2(double yF, double Rd, double Rc, double kinA, double kinB, double kinR, std::vector<std::vector<double>> yVector_Particle, bool GetStatisticInfo, int SpectrumUpdateTimes);
		~TsGSM2();

		//input solo parametri biologici
		void InitializeHistograms(int bins, double start, double end);
		void SetSpecificEnergySpectraCellNucleus();
		void ParallelGetInitialLethalNonLethalDamages(vector<double> &p0x, vector<double> &p0xy, double zn, int NumberOfSamples);
		vector<vector<double>> GetInitialLethalNonLethalDamages(double zn, int NumberOfSamples);
		std::vector<double> GetSurvivalDomain(double a, double b, double r, std::vector<double>p0x, std::vector<double> p0y,  std::vector<double>p0x_var, std::vector<double> p0y_var);
		std::vector<double> GetMultieventNucleus(double D, int NumberOfSamples) { return fSpecificEnergy_C->GetHfzMultiEvent(D,NumberOfSamples); };
		std::vector<double> GetMultieventNucleusVariance() { return fSpecificEnergy_C->GetHfzMultievent_var(); };
		std::vector<std::vector<double>> GetMultieventNucleusParticleContribution() {return fSpecificEnergy_C -> GetMultiEventParticleContribution();};
		
		vector<vector<double>> GSM2StochasticEvolution(vector<double> p0x, vector<double> p0y, double Tmax);
		//void GSM2EvolutionDoseRate();
		void CalculateInitialDamageDistribution(std::vector<double> NonLethalDamage, std::vector<double> LethalDamage);
		void GetStatisticInfo(int Binindex, double variable,std::vector<double> &CountMap, std::vector<double> &FirstMomentMap, std::vector<double> &SecondMomentMap, std::vector<double> &Variance, std::vector<double> &StandardDeviation);
		void InitializeStatistic();


		vector<int> SampleDistribution(vector<double> p, int Niter); //obtain a vector of indeces distributed as p


		vector<double> GetZn() {return zBinCenter;};
		vector<double> GetzBinWidth() {return zBinWidth;};
	private:
		//Histograms settings
		int fzBins; 
		double fzStart; 
		double fzEnd;
		std::vector<double> zBinLimit;
		std::vector<double> zBinWidth;
		std::vector<double> zBinCenter;


		std::vector<double> fCountMap_p0x, fFirstMomentMap_p0x, fSecondMomentMap_p0x, fVariance_p0x, fStandardDeviation_p0x;
		std::vector<double> fCountMap_p0y, fFirstMomentMap_p0y, fSecondMomentMap_p0y, fVariance_p0y, fStandardDeviation_p0y;
		//y events per particle matrix
		std::vector<std::vector<double>> fyVector_Particle;

		//microdosimetric z on cell domain
		std::vector<double> hfz_D;
		std::vector<double> hzfz_D;
		std::vector<double> hzfz_cumulative_D;
		double zF_D;

		//microdosimetric z on cell nucleus
		std::vector<double> hfz_C;
		std::vector<double> hzfz_C;
		std::vector<double> hzfz_cumulative_C;
		double zF_C;

		//multievent on cell population for each zn
		std::vector<vector<double>> hzfz_mevent_cell;


		double GSM2Model_yF, GSM2Model_kappa, GSM2Model_lambda, GSM2Model_rd, GSM2Model_rc, GSM2_a, GSM2_r, GSM2_b;

		vector<vector<double>> p0x_zn, p0y_zn;

		bool fGetStatisticInfo;
		int fSpectrumUpdateTimes;
		TsSpecificEnergy *fSpecificEnergy_D, *fSpecificEnergy_C;
};



#endif
