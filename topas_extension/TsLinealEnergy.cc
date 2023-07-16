// Extra Class for TsYScorer

// Claculate RBE 
// Author: Hongyu Zhu
// Date: 04/19/2019

#include "TsLinealEnergy.hh"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <thread>
#include <chrono>
#include<vector>

using namespace std;

TsLinealEnergy::TsLinealEnergy(std::vector<double> yVector, std::vector<std::vector<double>> yVector_Particle)
	:fyVector(yVector), fyVector_Particle(yVector_Particle)
{
	fSpectrumUpdateTimes = 100;
	InitializeMicrodosimetricSpectrum();
	// calculate statistic information
	InitializeStatistic();
	//GetSpectrum();
	cout <<"Start statistic error calcuation ..."<<endl;
	clock_t start,end;
	start = clock();
	if (fSpectrumUpdateTimes>fyVector.size()){
		int size = fyVector.size();
		fSpectrumUpdateTimes = min(size,1000);
		cout << "Statistic Update Frequency is larger than scored events("<<fyVector.size() <<")!"<<endl;
		cout << "Statistic Update Frequency is reset as "<<fSpectrumUpdateTimes<<endl;
	}

	for (int i=1; i<=fSpectrumUpdateTimes; i++){
		int dataSize = ceil(fyVector.size()/fSpectrumUpdateTimes);
		std::vector<double> dataVector;
		for(int j=0; j<dataSize*i && j<fyVector.size(); j++)
			dataVector.push_back(fyVector[j]);

		Calculatefy(dataVector);
		dataVector.clear();

		if(fSpectrumUpdateTimes>=100 )
		{       
			int tenPercent = ceil(fSpectrumUpdateTimes/10);
			if (i%tenPercent==0){
				end = clock(); 
				float duration = (float) (end - start)/ CLOCKS_PER_SEC;
				cout << (i/tenPercent)*10<< " % of statistic information update finished, total time used :"<<duration<<" sec; ("<<duration/60<<" min)"<<endl;
			}
		}
	}
	cout<<"\n"<<endl;

	GetSpectrum();
	GetErrorPropagation();
};

TsLinealEnergy::~TsLinealEnergy()
{};

void TsLinealEnergy::InitializeMicrodosimetricSpectrum()
{

	yF = 0.;
	yD = 0.;
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
	
	

}

void TsLinealEnergy::InitializeStatistic()
{
	fFirstMomentMap.resize(yBinNum);
	fSecondMomentMap.resize(yBinNum);
	fCountMap.resize(yBinNum);
	fVariance.resize(yBinNum);
	fStandardDeviation.resize(yBinNum);

	yF_var=0;
	yF_std=0;
	yD_var=0;
	yD_std=0;
	ydy_var.resize(yBinNum);
	ydy_std.resize(yBinNum);
	yfy_var.resize(yBinNum);
	yfy_std.resize(yBinNum);
	dy_var.resize(yBinNum);
	dy_std.resize(yBinNum);
}

void TsLinealEnergy::GetSpectrum()
{
	cout<<"yVector size="<<fyVector.size()<< endl;
	cout<<"yVector_Particle size" << fyVector_Particle.size();
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


	
	for (int i=0;i<yBinNum;i++)
	{
		hfy[i] = hfy[i]/(BinWidth[i]*nnum);                    // normalization. divide by bin width * number of entries
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

	//******************************************************************
	//                 Validate f(y) & calculate yF
	//******************************************************************
	double Probability_fy =0;
	for(int i=0;i<yBinNum;i++)
		Probability_fy += hfy[i]*BinWidth[i];
//	cout << "sum of f(y)*delta_y ="<< Probability_fy<<endl;    

	//calculate yF
	yF=0;
	for (int i=0;i<yBinNum;i++){
		yF = yF + hyfy[i]*BinWidth[i];          // multiply by bin width
	}

	for (int i=0;i<yBinNum;i++){
		hdy[i] = hyfy[i]/yF;                                    //calculate d(y) = y*f(y)/yF (cf. Burigo et al., NIMB 320 (2014))
		hydy[i] = (BinLimit[i]+BinLimit[i+1])/2*hdy[i];         // calculate y*d(y) = BinCenter * d(y)
	}


	//******************************************************************
	//               Validate d(y) & calculate yF
	//******************************************************************
	double Probability_dy =0;
	for(int i=0;i<yBinNum;i++)
		Probability_dy += hdy[i]*BinWidth[i];
//	cout << "sum of d(y)*delta_y ="<< Probability_dy<<endl;

	//calculate yD
	yD=0;
	for (int i=0;i<yBinNum;i++){
		yD = yD + hydy[i]*BinWidth[i];          // multiply by bin width
	}
}

void TsLinealEnergy::Calculatefy(std::vector<double> dataVector)
{
	double * histfy = new double [yBinNum] {0};
	int nnum=0;
	for (std::vector<double>::const_iterator i = dataVector.begin(); i != dataVector.end(); ++i){
		for (int n=0;n<yBinNum;n++){
			if(*i<=BinLimit[n+1]){
				histfy[n] = histfy[n]+1;
				break;
			}
		}
		nnum=nnum+1;
	}

	for (int i=0;i<yBinNum;i++){
		histfy[i] = histfy[i]/(BinWidth[i]*nnum);    // normalization. divide by bin width * number of entries
		GetStatisticInfo(i, histfy[i]);              // Calculate statistic error
	}
}

void TsLinealEnergy::GetStatisticInfo(int Binindex, double variable)
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
	fCountMap[index] ++;
	recorededHistories = fCountMap[index];

	// Use numerically stable algoritm from Donald E. Knuth (1998).
	// The Art of Computer Programming, volume 2: Seminumerical Algorithms,
	// 3rd edn., p. 232. Boston: Addison-Wesley.
	// for x in data:
	//   n = n + 1
	//   delta = x - mean
	//   mean = mean + delta/n
	//   mom2 = mom2 + delta*(x - mean)
	//   variance = mom2/(n - 1)

	if ( fCountMap[index]==1){
		// Initialize values to account for all previous histories having zero value
		// If we want Mean but don't want SecondMoment, can use a faster method at end of scoring.
		mean = x/recorededHistories;
		fFirstMomentMap[index] = mean;
		mom2 = (recorededHistories-1)*mean*mean + (x - mean)*(x - mean);
		fSecondMomentMap[index] = mom2;
	}
	else
	{
		mean = fFirstMomentMap[index];
		delta = x - mean;

		mean += delta/recorededHistories;
		mom2 = fSecondMomentMap[index];
		mom2 += delta*(x-mean);

		fSecondMomentMap[index] = mom2;
		fFirstMomentMap[index] = mean;
		fVariance[index] = fSecondMomentMap[index]/(recorededHistories-1);
		fStandardDeviation[index] = sqrt(fVariance[index]);
	}

	// if (x>0){
	// cout << "index="<<index<<" x="<<x<<" recorededHistories="<<recorededHistories<<" fCountMap[index]"<<fCountMap[index];
	// cout<<" FirstMoment="<<fFirstMomentMap[index]<<" SecondMoment=" <<fSecondMomentMap[index]
	// << " var="<<fVariance[index]<<" std="<<fStandardDeviation[index]<<endl;
	// }
}

void TsLinealEnergy::GetErrorPropagation()
{
	// calculate statistic error for yF
	for (int i = 0; i < yBinNum; i++) 
		yF_var += pow(BinWidth[i],2)*pow((BinLimit[i]+BinLimit[i+1])/2, 2)*fVariance[i];
	yF_std = sqrt(yF_var);

	// calculate statistic error for yD
	for (int i = 0; i < yBinNum; i++)
	{
		double aa = BinWidth[i]*pow((BinLimit[i]+BinLimit[i+1])/2, 2);
		yD_var += pow(aa/yF,2)*fVariance[i] + pow(aa*hfy[i]/(yF*yF),2)*yF_var;
	}
	yD_std = sqrt(yD_var);

	// calculate statistic error for yd(y), yf(y), d(y)
	for (int i = 0; i < yBinNum; i++)
	{
		double yi2 = pow((BinLimit[i]+BinLimit[i+1])/2, 2);
		ydy_var[i] = pow(yi2/yF, 2)*fVariance[i]+pow(yi2*hfy[i]/(yF*yF), 2)*yF_var;
		ydy_std[i] = sqrt(ydy_var[i]);   

		yfy_var[i] =  yi2*fVariance[i];
		yfy_std[i] = sqrt(yfy_var[i]);

		double yi = (BinLimit[i]+BinLimit[i+1])/2;
		dy_var[i] = pow(yi/yF,2)*fVariance[i] + pow(yi*hfy[i]/(yF*yF), 2)*yF_var;
		dy_std[i] = sqrt(dy_var[i]);
	}
}
