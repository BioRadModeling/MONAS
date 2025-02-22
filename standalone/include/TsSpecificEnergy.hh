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
// Claculate GSM2
// Author: Giorgio Cartechini
// Date: 13/11/2021

#ifndef TsSpecificEnergy_hh
#define TsSpecificEnergy_hh

#include <stdint.h>
#include <vector>
#include <string>
#include<random>
#include<iostream>

//#include "globals.hh"
//#include "G4RandomDirection.hh"
//#include "G4SystemOfUnits.hh"
//#include "g4root.hh"
using namespace std;

class TsSpecificEnergy 
{
	public:
		TsSpecificEnergy(std::vector<std::vector<double>> yVector_Particle, double radius, bool GetStatisticInfo, int SpectrumUpdateTimes);
		~TsSpecificEnergy();

		void InizializeHistograms();
		void InitializeStatistic();
		void InitializeStatisticMultievent();
		void Calculatefz(std::vector<double> dataVector);
		void CalculatefzMultievent(std::vector<double> dataVector);
		void GetStatisticInfo(int Binindex, double variable,std::vector<double> &CountMap, std::vector<double> &FirstMomentMap, std::vector<double> &SecondMomentMap, std::vector<double> &Variance, std::vector<double> &StandardDeviation);
		void GetErrorPropagation();
		void GetErrorPropagationMultievent();
		void SetSpecificEnergySpectra();


		std::vector<double> GetBinCenter(){return zBinCenter;};
		std::vector<double> GetBinWidth() {return zBinWidth;};
		std::vector<double> GetBinLimit() {return zBinLimit;};

		std::vector<double> GetHfz(){return hfz;};
		std::vector<double> GetHfz_std() {return fStandardDeviation;};
		std::vector<double> GetHzfz(){return hzfz;};
		std::vector<double> GetHzfz_std() {return hzfz_std;};
		std::vector<double> GetHzfzCumulative() {return hzfz_cumulative;};
		std::vector<std::vector<double>> GetParticleContribution() {return zParticleContibution;}
		double GetzF() {return zF;};
		double GetzF_std() {return zF_std;};

		void ParallelGetHfzMultiEvent(std::vector<std::vector<double>> &zMultievent_Particle, double dose, int NumberOfSamples);
		std::vector<double> GetHfzMultiEvent(double dose, int NumberOfSamples);
		std::vector<std::vector<double>> GetMultiEventParticleContribution() {return zMultieventParticleContibution;};
		void CalculateMultieventStatisticUncertainty();
		std::vector<double> GetHfzMultievent_var() {return fVarianceMultievent;};
	private:


		int zBins;
		double zStart;
		double zEnd;

		double y2z_factor;

		std::vector<std::vector<double>> fyVector_Particle;
		double fRadius;

		std::vector<double> zBinLimit, zBinCenter, zBinWidth;
		std::vector<double> hfz;
		std::vector<double> hzfz, hzfz_var, hzfz_std, hzfzMultievent_var, hzfzMultievent_std;
		std::vector<double> hzfz_cumulative;
		std::vector<double> fFirstMomentMap;
		std::vector<double> fSecondMomentMap;
		std::vector<double> fCountMap;
		std::vector<double> fVariance;
		std::vector<double> fStandardDeviation;
		
		std::vector<double> fFirstMomentMapMultievent;
		std::vector<double> fSecondMomentMapMultievent;
		std::vector<double> fCountMapMultievent;
		std::vector<double> fVarianceMultievent;
		std::vector<double> fStandardDeviationMultievent;
		
		int fSpectrumUpdateTimes;
		bool fGetStatisticInfo;
		
		std::vector<std::vector<double>> zParticleContibution;
		std::vector<std::vector<double>> zMultieventParticleContibution;
		std::vector<std::vector<double>> zMultieventVector_Particle;
		double **hfz_particle;

		double zF, zF_var, zF_std;
};



#endif
