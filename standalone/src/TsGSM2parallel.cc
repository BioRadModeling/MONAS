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

#include "TsGSM2parallel.hh"
#include "TsSpecificEnergy.hh"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <thread>
#include <algorithm>
#include <numeric> // Addition to sum a vector


//#include "globals.hh"
//#include "G4RandomDirection.hh"
//#include "G4SystemOfUnits.hh"
//#include "g4root.hh"


using namespace std;

TsGSM2::TsGSM2(double yF, double Rd, double Rc, double kinA, double kinB, double kinR, std::vector<double> yVector, std::vector<std::vector<double>> yVector_Particle, bool GetStatisticInfo, int SpectrumUpdateTimes)
	:GSM2Model_yF(yF), GSM2Model_rd(Rd), GSM2Model_rc(Rc), GSM2_a(kinA), GSM2_b(kinB), GSM2_r(kinR), fyVector(yVector), fyVector_Particle(yVector_Particle), fGetStatisticInfo(GetStatisticInfo), fSpectrumUpdateTimes(SpectrumUpdateTimes)
{
	// Old Kappa and lambda
	//double nDBS = 139.6*exp(0.0002568*GSM2Model_yF) -92.28*exp(-0.01855*GSM2Model_yF);
	
	// New Kappa formulation
	
	double nDBS = CalculateKappaFromSpectra();
	
	GSM2Model_kappa = nDBS*pow(GSM2Model_rd/GSM2Model_rc,3);
	GSM2Model_lambda = GSM2Model_kappa*1e-3;

	cout << "************** GSM2 **************\n"
		<< "Kappa: " << GSM2Model_kappa <<endl
		<< "Lambda: " << GSM2Model_lambda << endl
		<< "Rd[um]: " <<GSM2Model_rd << endl
		<< "Rc[um]: " <<GSM2Model_rc << endl
		<< "kinetic a: " << GSM2_a << endl
		<< "kinetic b: " << GSM2_b << endl
		<< "kinetic r: " << GSM2_r <<endl
		<< "**********************************\n";

	TsSpecificEnergy* zSpectra_D = new TsSpecificEnergy(fyVector_Particle, GSM2Model_rd, fGetStatisticInfo, fSpectrumUpdateTimes);
	fSpecificEnergy_D = zSpectra_D;

	//GetHistoInfo
	zBinCenter = fSpecificEnergy_D->GetBinCenter();
	zBinLimit = fSpecificEnergy_D-> GetBinLimit();
	zBinWidth = fSpecificEnergy_D->GetBinWidth();
	fzBins = zBinCenter.size();

	zF_D = fSpecificEnergy_D->GetzF();
	hzfz_cumulative_D = fSpecificEnergy_D->GetHzfzCumulative();
	
	TsSpecificEnergy* zSpectra_C = new TsSpecificEnergy(fyVector_Particle, GSM2Model_rc, fGetStatisticInfo, fSpectrumUpdateTimes);
	fSpecificEnergy_C = zSpectra_C;
	zF_C = fSpecificEnergy_C->GetzF();
	hzfz_cumulative_C = fSpecificEnergy_C -> GetHzfzCumulative();
};

TsGSM2::~TsGSM2()
{};

// New method for Kappa
double TsGSM2::CalculateKappaFromSpectra()
{
	// New formulation of Kappa and Lambda from PARTRAC simulations on DSBsites (Kundrát, Baiocco et al.)	
	// "Total DSBsites yield" parameters (H,H_sec,He,Li,Be,B,C) e- and other missing (first and ninth column)
	
	vector<double> p1{6.8,6.8,6.8,6.8,6.8,6.8,6.8};  	                  
	vector<double> p2{0.1773,0.1773,0.1471,0.1653,0.1425,0.1587,0.156};
	vector<double> p3{0.9314,0.9314,1.038,0.8782,0.95,0.8714,0.9214};
	vector<double> p4{0,0,0.006239,0.004284,0.005151,0.004345,0.005245};
	vector<double> p5{0,0,1.582,1.406,1.407,1.389,1.395};
	std::cout<<"DEBUG 3\n";
	// Calculation of particle contribution to y
	hfy.resize(yBinNum);
	hdy.resize(yBinNum);
	hyfy.resize(yBinNum);
	hydy.resize(yBinNum);
	hfy_particle = new double *[yBinNum];
	for (int i=0; i<yBinNum; i++)
		hfy_particle[i] =new double [10];	

	BinLimit.resize(yBinNum+1);
	BinWidth.resize(yBinNum);

	BinLimit[0]=0.1;
	for (int i=0;i<yBinNum;i++)
	{
		double aa = (double)((i+1)/yBinMagnitudeInterval);
		BinLimit[i+1] = pow(10,(aa -1.0));
		BinWidth[i] = BinLimit[i+1]-BinLimit[i];
	}
	
	
	int nnum=0;
	int index=0;
	for (std::vector<double>::const_iterator i = fyVector.begin(); i != fyVector.end(); ++i){
		for (int n=0;n<yBinNum;n++){
			if(*i<=BinLimit[n+1]){
				hfy[n] = hfy[n]+1;
				for(int particle = 0; particle<10; particle++)
				{
					hfy_particle[n][particle] += fyVector_Particle[index][particle];
				}
				break;
			}
		}
		nnum=nnum+1;
		index++;
	}
	std::cout<<"DEBUG 4\n";
	for (int i=0;i<yBinNum;i++)
		{
		hfy[i] = hfy[i]/(BinWidth[i]*nnum); 	              // normalization. divide by bin width * number of entries
		hyfy[i] = (BinLimit[i]+BinLimit[i+1])/2*hfy[i];        // calculate y*f(y) = BinCenter * BinContent
		
		//Calculate particle contibution
		std::vector<double> contribution;
		for(int particle = 0; particle<10; particle++)
		{
			if(hfy_particle[i][9] == 0) 
				contribution.push_back(0.);
			else
				contribution.push_back(hfy_particle[i][particle]/hfy_particle[i][9]);
			
		}
		yParticleContibution.push_back(contribution);
	}
	
	std::cout<<"DEBUG 5\n";
	// Calculate yF for each particle (multiply by bin width)
	for (int particle = 0; particle<9; particle++)
	{std::cout<<"DEBUG 6\n";
		for (int i=0;i<yBinNum;i++){
			yF[particle] += hyfy[i]*BinWidth[i]*yParticleContibution[i][particle]; 
		}      
	}
	
	// Evaluate total contribution (std::accumulate)
	sum_BinWidth = accumulate(BinWidth.begin(), BinWidth.end(), 0);
	for (int particle = 0; particle<9; particle++){
		TotalContributionParticle[particle] = 0;
		for (int i=0;i<yBinNum;i++){
			TotalContributionParticle[particle] = TotalContributionParticle[particle] + BinWidth[i]*yParticleContibution[i][particle];    
		} 
		TotalContributionParticle[particle] = TotalContributionParticle[particle]/sum_BinWidth;     
	}
	
	// Calculate Kappa for each particle
	for (int particle = 0; particle<9; particle++){
		if (particle == 0)
			KappaParticle[particle] = 9*(p1[particle]+pow((p2[particle]*yF[particle]),p3[particle]))/(1+pow((p4[particle]*yF[particle]),p5[particle]));
		else if (particle == 8)
			KappaParticle[particle] = 9*(p1[particle-2]+pow((p2[particle-2]*yF[particle]),p3[particle-2]))/(1+pow((p4[particle-2]*yF[particle]),p5[particle-2]));
		else
			KappaParticle[particle] = 9*(p1[particle-1]+pow((p2[particle-1]*yF[particle]),p3[particle-1]))/(1+pow((p4[particle-1]*yF[particle]),p5[particle-1]));
	}
	
	// Return Kappa weighted with particle contributions
	sum_TotalContributionParticle = accumulate(TotalContributionParticle.begin(), TotalContributionParticle.end(), 0);
	KappaValue = 0;
	for(int particle = 0; particle<9; particle++){
		KappaValue = KappaValue + KappaParticle[particle]*TotalContributionParticle[particle];
	}
	return KappaValue/sum_TotalContributionParticle;
}

void TsGSM2::ParallelGetInitialLethalNonLethalDamages(vector<double> &p0x, vector<double> &p0y, double zn, int NumberOfSamples)
{
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
			z_tot = zBinCenter[0]; //CHIEDERE

		//PASSAGGIO 3
		//ESTRAGGP I DANNI DA TUTTE LE TRACCE
		std::poisson_distribution<int> poissonX(GSM2Model_kappa*z_tot);
		double x0 = poissonX(generator); //numero di danni

		std::poisson_distribution<int> poissonY(GSM2Model_lambda*z_tot);
		double y0 = poissonY(generator); //numero di danni

		p0x.push_back(x0);
		p0y.push_back(y0);

	} //chiudo samples

}

vector<vector<double>> TsGSM2::GetInitialLethalNonLethalDamages(double zn, int NumberOfSamples)
{
	vector<double> p0x;
	vector<double> p0y;

	//Parallel computation of damages
	unsigned int NThreads = std::thread::hardware_concurrency();
	NThreads--;
	if(NThreads == 0)
		NThreads = 1;

	thread t[NThreads];
	vector<double> p0x_worker[NThreads];
	vector<double> p0y_worker[NThreads];
	int NumberOfSamples_worker = ceil(NumberOfSamples/NThreads);

	for(int i=0; i < NThreads; i++)
		t[i] = thread(&TsGSM2::ParallelGetInitialLethalNonLethalDamages, this,  ref(p0x_worker[i]), ref(p0y_worker[i]), zn, NumberOfSamples_worker);

	for(int i=0; i < NThreads; i++)
		t[i].join(); 

	for(int i=0; i < NThreads; i++)
	{
		for(int j=0; j<p0x_worker[i].size(); j++)
		{
			p0x.push_back(p0x_worker[i][j]);
			p0y.push_back(p0y_worker[i][j]);
		} 
	}


	int Bins_danno_x =  500;//*std::max_element(p0x.begin(), p0x.end()) + 1;
	int Bins_danno_y =  500;//*std::max_element(p0y.begin(), p0y.end()) + 1;


	std::vector<double> p0x_dis(Bins_danno_x, 0.);
	std::vector<double> p0y_dis(Bins_danno_y, 0.);

	for(int bin = 0; bin<Bins_danno_x; bin++)
		p0x_dis[bin] = 1.*std::count(p0x.begin(), p0x.end(), bin)/p0x.size();

	for(int bin = 0; bin<Bins_danno_y; bin++)
		p0y_dis[bin] = 1.*std::count(p0y.begin(), p0y.end(), bin)/p0y.size();

	//CALCULATE STATISTICAL UNCERTAINTY
	InitializeStatistic();
	if(fGetStatisticInfo)
	{
		for (int i=1; i<=fSpectrumUpdateTimes; i++)
		{
			int dataSize = ceil(p0x.size()/fSpectrumUpdateTimes);
			std::vector<double> dataVector_p0x, dataVector_p0y;
			for(int j=0; j<dataSize*i && j<p0x.size(); j++)
			{
				dataVector_p0x.push_back(p0x[j]);
				dataVector_p0y.push_back(p0y[j]);
			}
			CalculateInitialDamageDistribution(dataVector_p0x, dataVector_p0y);
			dataVector_p0x.clear();
			dataVector_p0y.clear();

		}
	}

	return {p0x_dis, p0y_dis, fVariance_p0x, fVariance_p0y};

	

	//Prendo la f(y). Faccio la f(z) in cui z è sul dominio grande 5um
	//calcolcol la zF(5um)
	//calcolo la fn(z) che è z_tot faccio lo spettro micro di questa z_tot
	//Poi fare l'intergrale f(z_tot)*Sn

}

void TsGSM2::CalculateInitialDamageDistribution(std::vector<double> NonLethalDamage, std::vector<double> LethalDamage)
{
	int Bins_danno_x =  500;//*std::max_element(p0x.begin(), p0x.end()) + 1;
	int Bins_danno_y =  500;//*std::max_element(p0y.begin(), p0y.end()) + 1;


	std::vector<double> p0x_dis(Bins_danno_x, 0.);
	std::vector<double> p0y_dis(Bins_danno_y, 0.);

	for(int bin = 0; bin<Bins_danno_x; bin++)
	{
		p0x_dis[bin] = 1.*std::count(NonLethalDamage.begin(), NonLethalDamage.end(), bin)/NonLethalDamage.size();
		GetStatisticInfo(bin, p0x_dis[bin], fCountMap_p0x, fFirstMomentMap_p0x, fSecondMomentMap_p0x, fVariance_p0x, fStandardDeviation_p0x);
	}

	for(int bin = 0; bin<Bins_danno_y; bin++)
	{
		p0y_dis[bin] = 1.*std::count(LethalDamage.begin(), LethalDamage.end(), bin)/LethalDamage.size();
		GetStatisticInfo(bin, p0y_dis[bin], fCountMap_p0y, fFirstMomentMap_p0y, fSecondMomentMap_p0y, fVariance_p0y, fStandardDeviation_p0y);
	}

}

void TsGSM2::InitializeStatistic()
{
	fCountMap_p0x.resize(500,0.);
	fCountMap_p0y.resize(500, 0.);
	
	fFirstMomentMap_p0x.resize(500,0.);
	fFirstMomentMap_p0y.resize(500, 0.);

	fSecondMomentMap_p0x.resize(500, 0.);
	fSecondMomentMap_p0y.resize(500, 0.);

	fVariance_p0x.resize(500, 0.);
	fVariance_p0y.resize(500, 0.);

	fStandardDeviation_p0x.resize(500, 0.);
	fStandardDeviation_p0y.resize(500, 0.);
	
}

void TsGSM2::GetStatisticInfo(int Binindex, double variable,std::vector<double> &CountMap, std::vector<double> &FirstMomentMap, std::vector<double> &SecondMomentMap, std::vector<double> &Variance, std::vector<double> &StandardDeviation)
{

	// Bin index
	int index;

	// Value from one specific bin
	double x;
	double mean;
	double delta;
	double mom2;
	double recorededHistories;

	// set value
	index  = Binindex;
	x =  variable;
	CountMap[index] ++;
	recorededHistories = CountMap[index];

	// Use numerically stable algoritm from Donald E. Knuth (1998).
	// The Art of Computer Programming, volume 2: Seminumerical Algorithms,
	// 3rd edn., p. 232. Boston: Addison-Wesley.
	// for x in data:
	//   n = n + 1
	//   delta = x - mean
	//   mean = mean + delta/n
	//   mom2 = mom2 + delta*(x - mean)
	//   variance = mom2/(n - 1)

	if ( CountMap[index]==1){
		// Initialize values to account for all previous histories having zero value
		// If we want Mean but don't want SecondMoment, can use a faster method at end of scoring.
		mean = x/recorededHistories;
		FirstMomentMap[index] = mean;
		mom2 = (recorededHistories-1)*mean*mean + (x - mean)*(x - mean);
		SecondMomentMap[index] = mom2;
	}
	else
	{
		mean = FirstMomentMap[index];
		delta = x - mean;

		mean += delta/recorededHistories;
		mom2 = SecondMomentMap[index];
		mom2 += delta*(x-mean);

		SecondMomentMap[index] = mom2;
		FirstMomentMap[index] = mean;
		Variance[index] = SecondMomentMap[index]/(recorededHistories-1);
		StandardDeviation[index] = sqrt(Variance[index]);
	}

	// if (x>0){
	// cout << "index="<<index<<" x="<<x<<" recorededHistories="<<recorededHistories<<" fCountMap[index]"<<fCountMap[index];
	// cout<<" FirstMoment="<<fFirstMomentMap[index]<<" SecondMoment=" <<fSecondMomentMap[index]
	// << " var="<<fVariance[index]<<" std="<<fStandardDeviation[index]<<endl;
	// }
}

std::vector<double> TsGSM2::GetSurvivalDomain(double a, double b, double r, std::vector<double> p0x, std::vector<double> p0y, std::vector<double> p0x_var, std::vector<double> p0y_var)
{
	double sn = p0x[0]*p0y[0];
	double sn_var = p0x[0]*p0x[0]*p0y_var[0] + p0y[0]*p0y[0]*p0x_var[0];

	for(int x=1; x<p0x.size(); ++x)
	{
		double c = 1.;
		for(int x0=1; x0<=x; ++x0)
		{
			c *= r*x0 / ((a+r)*x0 + b*x0*(x0-1));
		}
		sn += c*p0x[x]*p0y[0];
		sn_var += c*c*(p0x[x]*p0x[x]*p0y_var[0] + p0y[0]*p0y[0]*p0x_var[x]);
	}


	return {sn, sn_var};

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
