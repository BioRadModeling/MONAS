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

//#include "globals.hh"
//#include "G4RandomDirection.hh"
//#include "G4SystemOfUnits.hh"
//#include "g4root.hh"


using namespace std;

TsGSM2::TsGSM2(double yF, double Rd, double Rc, double kinA, double kinB, double kinR, std::vector<std::vector<double>> yVector_Particle, bool GetStatisticInfo, int SpectrumUpdateTimes)
	:GSM2Model_yF(yF), GSM2Model_rd(Rd), GSM2Model_rc(Rc), GSM2_a(kinA), GSM2_b(kinB), GSM2_r(kinR), fyVector_Particle(yVector_Particle), fGetStatisticInfo(GetStatisticInfo), fSpectrumUpdateTimes(SpectrumUpdateTimes)
{


	//Kappa and lambda
	double nDBS = 139.6*exp(0.0002568*GSM2Model_yF) -92.28*exp(-0.01855*GSM2Model_yF);
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
	unsigned int availThread = std::thread::hardware_concurrency();
	int NThreads = (int)availThread-1;
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
