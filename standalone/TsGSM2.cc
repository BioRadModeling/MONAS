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
// Extra Class for TsYScorer

// Claculate RBE 
// Author: Hongyu Zhu
// Date: 04/19/2019

#include "TsGSM2.hh"
#include "TsSpecificEnergy.hh"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>


//#include "globals.hh"
//#include "G4RandomDirection.hh"
//#include "G4SystemOfUnits.hh"
//#include "g4root.hh"


using namespace std;

TsGSM2::TsGSM2(double kappa, double lambda, double Rd, double Rc, double kinA, double kinB, double kinR, std::vector<double*> yVector_Particle)
	:GSM2Model_kappa(kappa), GSM2Model_lambda(lambda), GSM2Model_rd(Rd), GSM2Model_rc(Rc), GSM2_a(kinA), GSM2_b(kinB), GSM2_r(kinR), fyVector_Particle(yVector_Particle)
{
	cout << "************** GSM2 **************\n"
		<< "Kappa: " << GSM2Model_kappa <<endl
		<< "Lambda: " << GSM2Model_lambda << endl
		<< "Rd[um]: " <<GSM2Model_rd << endl
		<< "Rc[um]: " <<GSM2Model_rc << endl
		<< "kinetic a: " << GSM2_a << endl
		<< "kinetic b: " << GSM2_b << endl
		<< "kinetic r: " << GSM2_r <<endl
		<< "**********************************\n";


	TsSpecificEnergy* zSpectra_D = new TsSpecificEnergy(fyVector_Particle, GSM2Model_rd);
	fSpecificEnergy_D = zSpectra_D;

	//GetHistoInfo
	zBinCenter = fSpecificEnergy_D->GetBinCenter();
	zBinLimit = fSpecificEnergy_D-> GetBinLimit();
	zBinWidth = fSpecificEnergy_D->GetBinWidth();
	fzBins = zBinCenter.size();

	zF_D = fSpecificEnergy_D->GetzF();
	hzfz_cumulative_D = fSpecificEnergy_D->GetHzfzCumulative();
	
	TsSpecificEnergy* zSpectra_C = new TsSpecificEnergy(fyVector_Particle, GSM2Model_rc);
	fSpecificEnergy_C = zSpectra_C;
	zF_C = fSpecificEnergy_C->GetzF();
	hzfz_cumulative_C = fSpecificEnergy_C -> GetHzfzCumulative();
};

TsGSM2::~TsGSM2()
{};

vector<vector<double>> TsGSM2::GetInitialLethalNonLethalDamages(double zn, int NumberOfSamples)
{
	vector<double> p0x;
	vector<double> p0y;

	std::default_random_engine generator;
	std::uniform_real_distribution<double> uniform(0.0,1.);
	double rU;

	for(int k=0; k<NumberOfSamples; k++)
	{
		std::poisson_distribution<int> poisson(zn/zF_D); //N tracce
		double rP = poisson(generator); //NU


		double z_tot = 0.;
		if(rP>0)
		{
			for(int nu=0; nu<rP; ++nu)
			{
				//cout << nu << '\t' << poiss_rnd << endl;
				//campiono f(z) e sommo per ogni traccia nu
				rU = uniform(generator);
				for (int j = 0; j < fzBins; j++)
				{
					if (rU <= hzfz_cumulative_D[j])
					{
						z_tot += zBinCenter[j];
						break;
					}
				}

			} //chiudo su poisson nu
		}//chiudo if()
		else
			z_tot = zBinLimit[0]; //CHIEDERE

		//PASSAGGIO 3
		//ESTRAGGP I DANNI DA TUTTE LE TRACCE
		std::poisson_distribution<int> poissonX(GSM2Model_kappa*z_tot);
		double x0 = poissonX(generator); //numero di danni


		std::poisson_distribution<int> poissonY(GSM2Model_lambda*z_tot);
		double y0 = poissonY(generator); //numero di danni

		p0x.push_back(x0);
		p0y.push_back(y0);

	} //chiudo samples


	int Bins_danno_x =  500;//*std::max_element(p0x.begin(), p0x.end()) + 1;
	int Bins_danno_y =  500;//*std::max_element(p0y.begin(), p0y.end()) + 1;


	std::vector<double> p0x_dis(Bins_danno_x, 0.);
	std::vector<double> p0y_dis(Bins_danno_y, 0.);

	for(int bin = 0; bin<Bins_danno_x; bin++)
		p0x_dis[bin] = 1.*std::count(p0x.begin(), p0x.end(), bin)/p0x.size();

	for(int bin = 0; bin<Bins_danno_y; bin++)
		p0y_dis[bin] = 1.*std::count(p0y.begin(), p0y.end(), bin)/p0y.size();

	return {p0x_dis, p0y_dis};

	//Prendo la f(y). Faccio la f(z) in cui z è sul dominio grande 5um
	//calcolcol la zF(5um)
	//calcolo la fn(z) che è z_tot faccio lo spettro micro di questa z_tot
	//Poi fare l'intergrale f(z_tot)*Sn

}

double TsGSM2::GetSurvivalDomain(double a, double b, double r, vector<double> p0x, vector<double> p0y)
{
	double sn = p0x[0]*p0y[0];
	for(int x=1; x<p0x.size(); ++x)
	{
		double c = 1.;
		for(int x0=1; x0<=x; ++x0)
		{
			c *= r*x0 / ((a+r)*x0 + b*x0*(x0-1));
		}
		sn += c*p0x[x]*p0y[0];
	}

	return sn;

}

vector<vector<double>> TsGSM2::GSM2StochasticEvolution(vector<double> p0x, vector<double> p0y, double Tmax)
{
	double tt=0;
	double a = GSM2_a; //lethal damage
	double r = GSM2_r; //repair rate
	double b = GSM2_b; //

	vector<double> X, Y, Time;
	vector<int> tmp;
	double x0, y0;

	//simulate initial damages to evolve
	tmp =  SampleDistribution(p0x,1);
	x0 = tmp[0];

	tmp =  SampleDistribution(p0x,1);
	y0 = tmp[0];

	X.push_back(1.*x0);
	Y.push_back(1.*y0);
	Time.push_back(tt);

	int cntr = 0;
	while(tt<Tmax)
	{
		double h = r*X[cntr] + a*X[cntr] + b*X[cntr]*(X[cntr]-1);

		if(h==0)
		{
			cout <<"h = 0!";
			break;
		}

		std::default_random_engine generator;
		std::exponential_distribution<double> rexp(h);

		tt += rexp(generator);
		Time.push_back(tt);


		vector<double> v;
		v.push_back(r*X[cntr]/h);
		v.push_back(a*X[cntr]/h);
		v.push_back(( b*X[cntr]*(X[cntr]-1) )/h);

		tmp = SampleDistribution(v, 1);
		int rate = tmp[0];

		if(rate == 0) 
		{
			X.push_back(X[cntr]-1);
			Y.push_back(Y[cntr]);
		}
		else if(rate == 1)
		{
			X.push_back(X[cntr]-1);
			Y.push_back(Y[cntr]+1);
		}
		else
		{
			X.push_back(X[cntr]-2);
			Y.push_back(Y[cntr]+1);
		}

		cntr += 1;
	}

	vector<vector<double>> output = {Time, X, Y};
	return output;
}

/*
   void TsGSM2::GSM2EvolutionDoseRate()
   {
   double tt=0;
   double a = GSM2_a; //lethal damage
   double r = GSM2_r; //repair rate
   double b = GSM2_b; //
   double DoseRate = 0.;

   vector<int> tmp;
   vector<double> X, Y, Time;

   X.push_back(0.);
   Y.push_back(0.);

   Time.push_back(tt);
   int cntr = 0;
   while(tt<Tmax)
   {
   double h = r*X[cntr] + a*X[cntr] + b*X[cntr]*(X[cntr]-1) + DoseRate;

   if(h==0)
   {
   cout <<"h = 0!";
   break;
   }

   std::default_random_engine generator;
   std::exponential_distribution<double> rexp(h);

   tt += rexp(generator);
   Time.push_back(tt);


   vector<double> v;
   v.push_back(r*X[cntr]/h);
   v.push_back(a*X[cntr]/h);
   v.push_back(( b*X[cntr]*(X[cntr]-1) )/h);

   tmp = SampleDistribution(v,1);
   int rate = tmp[0];

   if(rate == 0) 
   {
   X.push_back(X[cntr]-1);
   Y.push_back(Y[cntr]);
   }
   else if(rate == 1)
   {
   X.push_back(X[cntr]-1);
   Y.push_back(Y[cntr]+1);
   }
   else if(rate==2)
   {
   X.push_back(X[cntr]-2);
   Y.push_back(Y[cntr]+1);
   }
   else 
   {
   std::default_random_engine generator;
   std::uniform_real_distribution<double> uniform(0.0,1.);
   double z_sampled;

   double rU = uniform(generator); //G4UniformRand();
   for (int j = 0; j < yBinNum-1; j++) 
   {
   if (rU <= hzfz_cumulative[j])
   {
   z_sampled = zBinCenter[j];
   break;
   }
}


std::poisson_distribution<int> poisson1(z_sampled*kappa); //N tracce
std::poisson_distribution<int> poisson2(z_sampled*lambda); //N tracce
double poiss_x = poisson1(generator); //NU
double poiss_y = poisson2(generator); //NU


X.push_back(X[cntr]+poiss_x);
Y.push_back(Y[cntr]+poiss_y);
} //close else

cntr++;
} // close while

vector<vector<double>> outpuut = {Time, X, Y};
}
*/

//Sample a discrete distribution p Niter times. The output is the index of histgram/distribution p
std::vector<int> TsGSM2::SampleDistribution(vector<double> p, int Niter)
{

	int N = p.size();
	std::vector<double> p_cumulative(N, 0.);
	p_cumulative[0] = p[0];
	// compute cumulative probabilities
	double sum_cumulative = p_cumulative[0];
	for (int i = 1; i < N; i++)
	{
		p_cumulative[i] = p_cumulative[i - 1] + p[i];
		sum_cumulative += p[i];
	}

	for (int i = 0; i < N; i++)
		p_cumulative[i] /= sum_cumulative;

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.);

	vector<int> output;
	for(int i=0; i<Niter; i++)
	{
		double rU = distribution(generator);
		for (int index = 0; index < N; index++)
		{
			if (rU <= p_cumulative[index])
			{
				output.push_back(index);
				break;
			}
		}
	}

	return output;
}
