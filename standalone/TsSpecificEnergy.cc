// Extra Class for TsYScorer

// Claculate RBE 
// Author: Hongyu Zhu
// Date: 04/19/2019

#include "TsSpecificEnergy.hh"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <thread>
#include <chrono>

//#include "globals.hh"
//#include "G4RandomDirection.hh"
//#include "G4SystemOfUnits.hh"
//#include "g4root.hh"


using namespace std;

TsSpecificEnergy::TsSpecificEnergy(std::vector<std::vector<double>> yVector_Particle, double radius, bool GetStatisticInfo, int SpectrumUpdateTimes)
	:fyVector_Particle(yVector_Particle), fRadius(radius), fGetStatisticInfo(GetStatisticInfo), fSpectrumUpdateTimes(SpectrumUpdateTimes)
{
	//Default parameters
	zBins = 100;
	zStart = 0.001;
	zEnd = 100;

	double pi = 3.1415;
	double rho = 1; //g/cm3
	y2z_factor = 0.16/(pi*rho*fRadius*fRadius);


	InizializeHistograms();
	// calculate statistic information
	InitializeStatistic();
	InitializeStatisticMultievent();

	if(fGetStatisticInfo)
	{


		for (int i=1; i<=fSpectrumUpdateTimes; i++){
			int dataSize = ceil(fyVector_Particle.size()/fSpectrumUpdateTimes);
			std::vector<double> dataVector;
			for(int j=0; j<dataSize*i && j<fyVector_Particle.size(); j++)
				dataVector.push_back(fyVector_Particle[j][9]);
			Calculatefz(dataVector);
			dataVector.clear();

		}
	}
	SetSpecificEnergySpectra();
	GetErrorPropagation();

};

TsSpecificEnergy::~TsSpecificEnergy()
{};

void TsSpecificEnergy::InizializeHistograms()
{

	zBinLimit.resize(zBins+1, 0.);
	zBinWidth.resize(zBins, 0.);
	zBinCenter.resize(zBins, 0.);

	hfz.resize(zBins, 0.);
	hzfz.resize(zBins, 0.);
	hzfz_cumulative.resize(zBins,0.);

	hfz_particle = new double *[zBins];
	for (int i=0; i<zBins; i++)
		hfz_particle[i] =new double [10];	


	zMultieventParticleContibution.resize(zBins, std::vector<double>(10,0.));

	zBinLimit[0] = zStart; //lowest z value
	double zMax = zEnd; //highest z value
	double step = (log10(zMax) - log10(zStart))/zBins;

	double Binlog = log10(zStart);
	for (int i=0; i<zBins; i++)
	{
		Binlog += step;
		zBinLimit[i+1] = pow(10, Binlog);
		zBinWidth[i] = zBinLimit[i+1]-zBinLimit[i];
		zBinCenter[i] = (zBinLimit[i+1]+zBinLimit[i])/2.;
	}

}

void TsSpecificEnergy::InitializeStatistic()
{
	fFirstMomentMap.resize(zBins);
	fSecondMomentMap.resize(zBins);
	fCountMap.resize(zBins);
	fVariance.resize(zBins);
	fStandardDeviation.resize(zBins);

	fFirstMomentMapMultievent.resize(zBins);
	fSecondMomentMapMultievent.resize(zBins);
	fCountMapMultievent.resize(zBins);
	fVarianceMultievent.resize(zBins);
	fStandardDeviationMultievent.resize(zBins);

	zF_var=0.;
	zF_std=0.;
	hzfz_var.resize(zBins);
	hzfz_std.resize(zBins);

	hzfzMultievent_var.resize(zBins);
	hzfzMultievent_std.resize(zBins);
}

void TsSpecificEnergy::InitializeStatisticMultievent()
{

	fFirstMomentMapMultievent.resize(zBins);
	fSecondMomentMapMultievent.resize(zBins);
	fCountMapMultievent.resize(zBins);
	fVarianceMultievent.resize(zBins);
	fStandardDeviationMultievent.resize(zBins);

	hzfzMultievent_var.resize(zBins);
	hzfzMultievent_std.resize(zBins);
}

void TsSpecificEnergy::Calculatefz(std::vector<double> dataVector)
{	

	//CALCULTAE HISTO f(z_domain)
	int ScorerEvents = dataVector.size();

	double * histfz = new double [zBins] {0};
	int nnum=0;
	for (int i = 0; i < ScorerEvents; ++i){
		for (int n=0;n<zBins;n++){
			if(dataVector[i]*y2z_factor<=zBinLimit[n+1]){
				histfz[n] = histfz[n]+1;
				break;
			}
		}
		nnum=nnum+1;
	}

	for (int i=0;i<zBins;i++){
		histfz[i] = histfz[i]/(zBinWidth[i]*nnum);    // normalization. divide by bin width * number of entries
		GetStatisticInfo(i, histfz[i], fCountMap, fFirstMomentMap, fSecondMomentMap, fVariance, fStandardDeviation);              // Calculate statistic error
	}
}

void TsSpecificEnergy::CalculatefzMultievent(std::vector<double> dataVector)
{	

	//CALCULTAE HISTO f(z_domain)
	int ScorerEvents = dataVector.size();

	double * histfz = new double [zBins] {0};
	int nnum=0;
	for (int i = 0; i < ScorerEvents; ++i){
		for (int n=0;n<zBins;n++){
			if(dataVector[i]<=zBinLimit[n+1]){
				histfz[n] = histfz[n]+1;
				break;
			}
		}
		nnum=nnum+1;
	}

	for (int i=0;i<zBins;i++){
		histfz[i] = histfz[i]/(zBinWidth[i]*nnum);    // normalization. divide by bin width * number of entries
		GetStatisticInfo(i, histfz[i], fCountMapMultievent, fFirstMomentMapMultievent, fSecondMomentMapMultievent, fVarianceMultievent, fStandardDeviationMultievent);              // Calculate statistic error
	}
}

void TsSpecificEnergy::GetStatisticInfo(int Binindex, double variable,std::vector<double> &CountMap, std::vector<double> &FirstMomentMap, std::vector<double> &SecondMomentMap, std::vector<double> &Variance, std::vector<double> &StandardDeviation)
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

void TsSpecificEnergy::GetErrorPropagation()
{
	// calculate statistic error for yF
	for (int i = 0; i < zBins; i++) 
		zF_var += pow(zBinWidth[i],2)*pow(zBinCenter[i], 2)*fVariance[i];
	zF_std = sqrt(zF_var);


	// calculate statistic error for yd(y), yf(y), d(y)
	for (int i = 0; i < zBins; i++)
	{
		double zi2 = pow(zBinCenter[i], 2);

		hzfz_var[i] =  zi2*fVariance[i];
		hzfz_std[i] = sqrt(hzfz_var[i]);
	}
}

void TsSpecificEnergy::GetErrorPropagationMultievent()
{
	// calculate statistic error for yF
	//	for (int i = 0; i < zBins; i++) 
	//		zF_var += pow(zBinWidth[i],2)*pow(zBinCenter[i], 2)*fVarianceMultievent[i];
	//	zF_std = sqrt(zF_var);


	// calculate statistic error for yd(y), yf(y), d(y)
	for (int i = 0; i < zBins; i++)
	{
		double zi2 = pow(zBinCenter[i], 2);

		hzfzMultievent_var[i] =  zi2*fVariance[i];
		hzfzMultievent_std[i] = sqrt(hzfzMultievent_var[i]);
	}
}

void TsSpecificEnergy::SetSpecificEnergySpectra()
{

	//CALCULTAE HISTO f(z_domain)
	int ScorerEvents = fyVector_Particle.size();
	int index=0;   
	for (int i=0; i < ScorerEvents; ++i) //per ora prendo solo y_total
	{
		for (int n=0; n<zBins; ++n)
		{
			if(fyVector_Particle[i][9]*y2z_factor<=zBinLimit[n+1])
			{
				hfz[n] = hfz[n]+1;
				for(int particle = 0; particle<10; particle++)
				{
					hfz_particle[n][particle] += fyVector_Particle[index][particle]*y2z_factor;
				}
				break;
			}
		}
		index++;
	}

	//hfz[0] = 0.;
	//hfz_cell[0] = 0.;

	for (int i=0; i<zBins; ++i)
	{
		hfz[i] = hfz[i]/(zBinWidth[i]*ScorerEvents);                    // normalization. divide by bin width * number of entries
		hzfz[i] = zBinCenter[i]*hfz[i];        // calculate z*f(z) = BinCenter * BinContent

		//Calculate particle contibution
		std::vector<double> contribution;
		for(int particle = 0; particle<10; particle++)
		{
			if(hfz_particle[i][9] == 0) 
				contribution.push_back(0.);
			else
				contribution.push_back(hfz_particle[i][particle]/hfz_particle[i][9]);

		}
		zParticleContibution.push_back(contribution);

	}

	zF = 0.;
	for (int i=0;i<zBins;i++)
		zF += hzfz[i]*zBinWidth[i];          // multiply by bin width

	//CALCOLO LA CUMULATIVA PER IL SAMPLING
	hzfz_cumulative[0] = hzfz[0];

	// compute cumulative probabilities
	double sum_cumulative = hzfz_cumulative[0];
	for (int i = 1; i <zBins; i++)
	{
		hzfz_cumulative[i] = hzfz_cumulative[i-1] + hzfz[i];
		sum_cumulative += hzfz[i];
	}

	for (int i = 0; i < zBins; ++i)
		hzfz_cumulative[i] /= sum_cumulative;
};


void TsSpecificEnergy::ParallelGetHfzMultiEvent(std::vector<std::vector<double>> &zVectorPart, double dose, int NumberOfSamples)
{

	std::default_random_engine generator;
	std::uniform_real_distribution<double> uniform(0.0,1.);
	double rU;

	zVectorPart.resize(NumberOfSamples, std::vector<double>(10, 0.));

	for(int k=0; k<NumberOfSamples; k++)
	{       
		std::poisson_distribution<int> poisson_cell(dose/zF); //N tracce
		double poiss_cell = poisson_cell(generator); //NU

		std::vector<double> zMultievent_Particle(10,0.);
		double zMultievent = 0.;
		if(poiss_cell>0)
		{       
			for(int nu=0; nu<poiss_cell; ++nu)
			{       
				//cout << nu << '\t' << poiss_rnd << endl;
				//campiono f(z) e sommo per ogni traccia nu
				rU = uniform(generator);
				for (int j = 0; j < zBins-1; j++)
				{       
					if (rU <= hzfz_cumulative[j]) {
						zMultievent += zBinCenter[j];
						for(int particle=0; particle<10; particle++)
						{
							zMultievent_Particle[particle] += zBinCenter[j]*zParticleContibution[j][particle];
						}
						break;
					}

				}


			} //chiudo su poisson nu
		}//chiudo if()
		else   
		{
			for(int particle=0; particle<10; particle++)
			{
				zMultievent_Particle[particle] += zBinCenter[0]*zParticleContibution[0][particle];
			}

		}

		zVectorPart[k]= zMultievent_Particle;
	} //chiudo samples

}

std::vector<double> TsSpecificEnergy::GetHfzMultiEvent(double dose, int NumberOfSamples)
{

	zMultieventVector_Particle.clear();
	std::vector<double> hfzMultievent(zBins, 0);
	std::vector<double> hzfzMultievent(zBins, 0);
	std::vector<std::vector<double>> hfzMultievent_Particle(zBins, std::vector<double>(10,0.));
	//hfzMultievent_Particle = new double *[zBins];
	//for (int i=0; i<zBins; i++)
	//hfzMultievent_Particle[i] =new double [10];	


	unsigned int NThreads = std::thread::hardware_concurrency();
	NThreads = 1;
	if(NThreads == 0)
		NThreads = 1;
	thread t[NThreads];
	std::vector<std::vector<double>> zMultieventVector_Particle_worker[NThreads];

	int NumberOfSamples_worker = ceil(NumberOfSamples/NThreads);

	for(int i=0; i < NThreads; i++)
		t[i] = thread(&TsSpecificEnergy::ParallelGetHfzMultiEvent, this,  ref(zMultieventVector_Particle_worker[i]), dose, NumberOfSamples_worker);

	for(int i=0; i < NThreads; i++)
		t[i].join(); 

	//ABSORB VECTORS FROM WORKWERS
	for(int i=0; i<NThreads; i++)
	{
		for(int j=0; j<zMultieventVector_Particle_worker[i].size(); j++)
			zMultieventVector_Particle.push_back(zMultieventVector_Particle_worker[i][j]);
	}



	//CALCULATE HISTOGRAMS
	int zVectorSize = zMultieventVector_Particle.size();
	int nnum = 0;
	for(int i=0; i<zVectorSize; i++)
	{
		if(zMultieventVector_Particle[i][9]>zBinLimit[zBinLimit.size()-1])
		{
				nnum++;
		}

		for (int n=0;n<zBins;n++)
		{       
			double zMultieventTot = zMultieventVector_Particle[i][9];
			if(zMultieventTot<=zBinLimit[n+1])
			{       
				hfzMultievent[n] = hfzMultievent[n]+1;
				for(int particle=0; particle<10; particle++)
				{
					hfzMultievent_Particle[n][particle] += zMultieventVector_Particle[i][particle]; 
				}
				break;
			}
			
		}
	}

	for (int i=0;i<zBins;i++)
	{       
		hfzMultievent[i] = hfzMultievent[i]/(zBinWidth[i]*(zVectorSize-nnum)); //*NumberOfSamples_worker*NThreads);                    // normalization. divide by bin width * number of entries
		hzfzMultievent[i] = zBinCenter[i]*hfzMultievent[i];        // calculate z*f(z) = BinCenter * BinContent
		//Calculate particle contibution
		std::vector<double> contribution;
		for(int particle = 0; particle<10; particle++)
		{
			if(hfzMultievent_Particle[i][9] == 0) 
				zMultieventParticleContibution[i][particle] = 0.;
			else
			{

				zMultieventParticleContibution[i][particle] = hfzMultievent_Particle[i][particle]/hfzMultievent_Particle[i][9];
			}
		}

	}

	

	if(fGetStatisticInfo)
	{	
		CalculateMultieventStatisticUncertainty();
		GetErrorPropagation();
	}

	//	double aaa = 0;
	//	for(int i=0; i<zBins; i++)
	//		aaa += hfzMultievent[i]*zBinWidth[i];
	//	cout << "probability: " << dose << '\t' << aaa << endl; 
	return hfzMultievent;

}

void TsSpecificEnergy::CalculateMultieventStatisticUncertainty()
{
	InitializeStatisticMultievent();
	int size = zMultieventVector_Particle.size();	

	for (int i=1; i<=fSpectrumUpdateTimes; i++){
		int dataSize = ceil(size/fSpectrumUpdateTimes);
		std::vector<double> dataVector;
		for(int j=0; j<dataSize*i && j<size; j++)
			dataVector.push_back(zMultieventVector_Particle[j][9]);
		CalculatefzMultievent(dataVector);
		dataVector.clear();
	}

}

