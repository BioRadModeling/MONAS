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
		TsGSM2(double kappa, double lambda, double Rd, double Rc, double kinA, double kinB, double kinR, std::vector<double*> yVector_Particle);
		~TsGSM2();

		//input solo parametri biologici
		void InitializeHistograms(int bins, double start, double end);
		void SetSpecificEnergySpectraCellNucleus();
		vector<vector<double>> GetInitialLethalNonLethalDamages(double zn, int NumberOfSamples);
		double GetSurvivalDomain(double a, double b, double r, vector<double>p0x, vector<double> p0y);
		vector<double> GetMultieventNucleus(double D, int NumberOfSamples) { return fSpecificEnergy_C->GetHfzMultiEvent(D, NumberOfSamples); };

		vector<vector<double>> GSM2StochasticEvolution(vector<double> p0x, vector<double> p0y, double Tmax);
		//void GSM2EvolutionDoseRate();

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

		//y events per particle matrix
		std::vector<double* > fyVector_Particle;

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


		double GSM2Model_kappa, GSM2Model_lambda, GSM2Model_rd, GSM2Model_rc, GSM2_a, GSM2_r, GSM2_b;

		vector<vector<double>> p0x_zn, p0y_zn;

		TsSpecificEnergy *fSpecificEnergy_D, *fSpecificEnergy_C;
};



#endif
