// *************************************************************************
// * MONAS is a C++ package that calculates cell surviavl curvs and        *
// * dose dependednt RBE from microdosimetric spectra.			   *
// *									   *
// * Copyright © 2023 Giorgio Cartechini <giorgio.cartechini@miami.edu>	   *
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
#include "TsLinealEnergy.hh"

//#include "globals.hh"
//#include "G4RandomDirection.hh"
//#include "G4SystemOfUnits.hh"
//#include "g4root.hh"


using namespace std;

//CONTROLLARE GLI HISO (BIN START, BIN END)
class TsGSM2 
{
	public:
		TsGSM2(double yF, 
			double Rd, 
			double Rc, 
			double kinA, 
			double kinB, 
			double kinR,
			const string& ion,
			double LET, 
			std::vector<double> yVector, 
			std::vector<std::vector<double>> yVector_Particle, 
			std::vector<double> yVector_Nucleus, 
			std::vector<std::vector<double>> yVector_Particle_Nucleus, 
			bool GetStatisticInfo, 
			int SpectrumUpdateTimes);
		~TsGSM2();

		//input solo parametri biologici
		double CalculateKappaFromLET(const string& ion, double LET); // New Kappa formulation
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
		std::vector<double> fyVector;
		std::vector<std::vector<double>> fyVector_Particle;
		std::vector<double> fyVector_Nucleus;
		std::vector<std::vector<double>> fyVector_Particle_Nucleus;

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
		
		// Variable for Kappa from spectra calculation
		const int yBinNum = 100; // yBinNum == yBinMagnitude*yBinMagnitudeInterval
		const double yBinMagnitudeInterval = 20.;
		const double yBinMagnitude = 5.;
		double **hfy_particle;
		std::vector<double> hfy, hdy, hyfy, hydy, BinLimit, BinWidth;
		
		const string& GSM2_ion = "H";
		const double GSM2_LET = 10.0; // keV/um
		
		TsLinealEnergy* ySpectra_F;
};



#endif
