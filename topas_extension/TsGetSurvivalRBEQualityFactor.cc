// Extra Class for TsYScorer

#include "TsGetSurvivalRBEQualityFactor.hh"
#include "TsGSM2parallel.hh"
#include "TsSpecificEnergy.hh"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <thread>
#include <chrono>

using namespace std;

TsGetSurvivalRBEQualityFactor::TsGetSurvivalRBEQualityFactor(std::vector<std::vector<double>> yParticleContribution, std::vector<std::vector<double>> yVector_Particle, double*hBinLimit, double* hBinWidth,  double* hfy,double* hdy, double hyF, double hyD, double hyF_var, double hyD_var, std::vector<double> hfy_var, std::vector<double> hdy_var, int SpecLength, bool GetStatisticInfo, int SpectrumUpdateTimes, bool GetParticleContribution)
	:fyParticleContribution(yParticleContribution), fyVector_Particle(yVector_Particle),fBinLimit(hBinLimit), fBinWidth(hBinWidth), fhfy(hfy), fhdy(hdy), yF(hyF), yD(hyD), yF_var(hyF_var), yD_var(hyD_var), fy_var(hfy_var), dy_var(hdy_var),fSpecLength(SpecLength), fGetStatitisticInfo(GetStatisticInfo), fSpectrumUpdateTimes(SpectrumUpdateTimes), fGetParticleContribution(GetParticleContribution)
{
	// default values of MK model
	// 10% cell survival relative to 200 kVp X-rays for HSG cells
	MKModel_alpha0 = 0.13; // Unit:Gy-1
	MKModel_beta   = 0.05; // Unit:Gy-2,
	MKModel_alphaX = 0.19; // Unit:Gy-1
	MKModel_betaX  = 0.05;
	MKModel_rd     = 0.42; // Unit:um
	MKModel_Rn     = 6; //Unit:um
	MKModel_rho    = 1.;    // Unit:g/cm3
	MKModel_y0     = 150;  // Unit:keV/um

	//Split Dose MKM
	MKModel_D1 = 1; //Unit:Gy
	MKModel_D2 = 1; //Unit:Gy
	MKModel_ac = 2.187; //Unit:h-1
	MKModel_tr = 2.284; //Unit:h

	//GSM2 default parameters
	GSM2_kappa 	= 0.5;
	GSM2_lambda 	= 0.5;
	GSM2_rd 	= 0.42; //Unit: um
	GSM2_Rn 	= 6; 	//Unit: um
	GSM2_zStart 	= 0.01; //Gy
	GSM2_zEnd 	= 100; //Gy
	GSM2_zBins 	= 100; //Bins of spectra
	GSM2_a 		= 0.1;
	GSM2_b 		= 0.1;
	GSM2_r 		= 0.1;
	GSM2_alphaX 	= 0.19;
	GSM2_betaX 	= 0.05;
	//Macroscopic Doses
	Doses = {0, 1,2,3,4,5,6,7,8,9,10}; //Unit:Gy

	//MultiEventIterations
	MCMultieventIterations = 1e4;
	BioWeightFunctionDataFile = "BioWeightFuncData_interpolation.txt";

};

TsGetSurvivalRBEQualityFactor::~TsGetSurvivalRBEQualityFactor()
{};

// See reference: Y. Kase, et al.(2006). 
// Microdosimetric measurements and estimation of human cell survival for heavy-ion beams. 
// Radiat Res, 166, 629-38.
void TsGetSurvivalRBEQualityFactor::GetSurvWithMKModel_SaturationCorr()
{
	//Set parameters
	double alpha0 = MKModel_alpha0; // Unit:Gy-1
	double alphaX = MKModel_alphaX; // Unit:Gy-1
	double beta   = MKModel_beta;   // Unit:Gy-2
	double betaX  = MKModel_betaX;
	double rd     = MKModel_rd;     // Unit:um
	double y0     = MKModel_y0;     // Unit:keV/um
	double rho    = MKModel_rho;    // Unit:g/cm3
	double pi     = 3.1415;

	// *************************** Calculate y_star *****************************
	double Integrate_y = 0;
	int Ncomponents = fyParticleContribution[0].size();
	vector<double> ystar_Particle(Ncomponents,0.);

	for(int i=0; i<fSpecLength-1;i++) 
	{
		double yi = (fBinLimit[i]+fBinLimit[i+1])/2;  
		Integrate_y += (1-exp(-pow(yi,2)/pow(y0,2)))*fhfy[i]*fBinWidth[i];

		for(int comp=0; comp<Ncomponents; comp++)
			ystar_Particle[comp] += y0*y0*((1-exp(-pow(yi,2)/pow(y0,2)))*fhfy[i]*fBinWidth[i])*fyParticleContribution[i][comp]/yF;
	}
	double y_star = y0*y0*Integrate_y/yF;  // Unit:keV/um

	// ***************************** Calculate alphta *****************************
	double zstar = y_star/(rho*pi*pow(rd,2)); // unit:keV*cm3/g/um3
	zstar = zstar*0.16;                        // unit:J/kg = Gy
	vector<double> zstar_Particle(Ncomponents,0.);
	for(int comp=0; comp<Ncomponents; comp++)
	{
		zstar_Particle[comp] = ystar_Particle[comp]*0.16/(rho*pi*pow(rd,2));
	}
	double alpha = alpha0 + beta*zstar;       // unit:Gy-1

	// ***************************** Calculate Survival and RBE *****************************

	std::vector<double> S, RBE;
	std::vector<vector<double>> S_Particle, RBE_Particle;
	for(double D:Doses)
	{
		double s = exp( -(alpha*D) - (beta*D*D) );
		double rbe = 0;
		if(D>0)	
			rbe = (sqrt( (alphaX*alphaX) - (4*betaX*log(s)) ) - alphaX)/(2*betaX*D); 
		S.push_back(s);
		RBE.push_back(rbe);	
	}

	// ************************* Calculate statistic error ************************
	double ystar_var =0.;
	double zstar_var =0.;
	double alpha_var=0.;
	std::vector<double> S_var, RBE_var;

	for(int i=0; i<fSpecLength-1;i++) 
	{    
		double yi = (fBinLimit[i]+fBinLimit[i+1])/2;     
		double aa = (1-exp(-pow(yi,2)/pow(y0,2)));
		ystar_var += pow(y0*y0*aa*fBinWidth[i]/yF,2)*fy_var[i] +pow(y0*y0*aa*fhfy[i]*fBinWidth[i]/(yF*yF),2)*yF_var;
	}

	zstar_var = pow(0.16/(rho*pi*rd*rd),2)*ystar_var;
	alpha_var = beta*beta*zstar_var; // 0.16 convert keV*cm3/g/um3 into Gy
	for(int i=0; i<Doses.size(); ++i)
	{
		double s_var = Doses[i]*Doses[i]*alpha_var;
		double rbe_var = (1./(pow(Doses[i]*S[i],2)*(alphaX*alphaX-4*betaX*log(S[i]))))*s_var;
		S_var.push_back(s_var);
		RBE_var.push_back(rbe_var);

	}

	WriteMKMSurvival("MKModel_StaurationCorrected.csv", Doses, S, S_var, RBE, RBE_var);

	std::cout<<"************************************ Get RBE with MK method - Saturation Corrected ***************************"<<std::endl;
	std::cout<<"Doses\tSurvival\tRBE\n";
	for(int i=0; i<Doses.size(); i++)
	{
		std::cout<<Doses[i]<<'\t';
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<S[i];
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(S_var[i])<<")\t" <<'\t';
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<RBE[i];
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(RBE_var[i])<<")" << '\t' <<std::endl;
	}

	std::cout<<"\nSpectrum Parameters:"<<endl;
	std::cout<<"yF: " << yF  <<" ( " <<sqrt(yF_var) << " )" << endl
		<<"yD: " << yD  <<" ( " <<sqrt(yD_var) << " )"<< endl
		<<"y*: " << y_star <<" ( " <<sqrt(ystar_var) << " )\t" << endl
		<<"z*: " << zstar <<" ( " <<sqrt(zstar_var) << " )\t"  << endl;
	std::cout<<"***************************************************************************************\n"<<std::endl;
}

//
//RBE at 10% survival 
void TsGetSurvivalRBEQualityFactor::GetSurvWithMKModel_nonPoissonCorr()
{
	//Set parameters
	double alpha0 = MKModel_alpha0; // Unit:Gy-1
	double alphaX = MKModel_alphaX; // Unit:Gy-1
	double beta   = MKModel_beta;   // Unit:Gy-2
	double betaX  = MKModel_betaX;
	double rd     = MKModel_rd;     // Unit:um
	double Rn     = MKModel_rd;     //unit:um
	double y0     = MKModel_y0;     // Unit:keV/um
	double rho    = MKModel_rho;    // Unit:g/cm3
	double pi     = 3.1415;

	// ***************************** Calculate alphtaNP *****************************
	double zD = yD/(rho*pi*pow(rd,2)); // unit:keV*cm3/g/um3
	zD = zD*0.16;                        // unit:J/kg = Gy

	double znD = yD/(rho*pi*pow(Rn,2)); // unit:keV*cm3/g/um3
	znD = znD*0.16;

	double alphaP = alpha0 + beta*zD;       // unit:Gy-1
	double alphaNP = (1 - exp(-alphaP*znD))/znD;
	//cout<<"alphta = "<<alpha<<" beta = "<<beta<<endl;

	// ***************************** Calculate Survival and RBE *****************************
	std::vector<double> S, RBE;
	for(double D:Doses)
	{
		double s = exp( - (alphaNP*D) - (beta*D*D) );
		double rbe = 0;
		if(D>0)
			rbe = (sqrt( (alphaX*alphaX) - (4*betaX*log(s)) ) - alphaX)/(2*betaX*D); 
		S.push_back(s);
		RBE.push_back(rbe);
	}

	// ************************* Calculate statistic error ************************ 
	double alphaP_var = pow(0.16*beta/(rho*pi*pow(rd,2)),2)*yD_var; // 0.16 convert keV*cm3/g/um3 into Gy
	double znD_var = pow(0.16/(rho*pi*pow(Rn,2)),2)*yD_var;
	double alphaNP_var_a = pow(alphaP*exp(-alphaP*znD)/znD,2);
	double alphaNP_var_b = pow( (1-(alphaP*znD)*exp(-alphaP*znD))/(znD*znD) ,2);
	double alphaNP_var = (alphaNP_var_a*alphaP_var) + (alphaNP_var_b*znD_var);
	std::vector<double> S_var, RBE_var;

	for(int i=0; i<Doses.size(); ++i)
	{
		double s_var = Doses[i]*Doses[i]*alphaNP_var;
		double rbe_var = (1./(pow(Doses[i]*S[i],2)*(alphaX*alphaX-4*betaX*log(S[i]))))*s_var;
		S_var.push_back(s_var);
		RBE_var.push_back(rbe_var);

	}

	WriteMKMSurvival("MKModel_nonPoisson.csv", Doses, S, S_var, RBE, RBE_var);

	std::cout<<"************************************ Get RBE with MK method - non-Poisson ***************************"<<std::endl;
	std::cout<<"Doses\tSurvival\tRBE\n";
	for(int i=0; i<Doses.size(); i++)
	{
		std::cout<<Doses[i]<<'\t';
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<S[i];
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(S_var[i])<<")\t";
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<RBE[i];
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(RBE_var[i])<<" )"<<std::endl;
	}

	std::cout<<"\nSpectrum Parameters:"<<endl;
	std::cout<<"yF: " << yF  <<" ( " <<sqrt(yF_var) << " )" << endl
		<<"yD: " << yD  <<" ( " <<sqrt(yD_var) << " )"<< endl
		<<"znD: " << znD <<" ( " <<sqrt(znD_var) << " )"  << endl;
	std::cout<<"***************************************************************************************\n"<<std::endl;



} 

void TsGetSurvivalRBEQualityFactor::GetSurvWithMKModel_SplitDoseIrradiation()
{
	//Set parameters
	double alpha0 = MKModel_alpha0; // Unit:Gy-1
	double alphaX = MKModel_alphaX; // Unit:Gy-1
	double beta   = MKModel_beta;   // Unit:Gy-2
	double rd     = MKModel_rd;     // Unit:um
	double Rn     = MKModel_rd;     //unit:um
	double y0     = MKModel_y0;     // Unit:keV/um
	double rho    = MKModel_rho;    // Unit:g/cm3
	double pi     = 3.1415;

	//Split DoseIrradaiation parameters
	double D1 = MKModel_D1; //Unit:Gy
	double D2 = MKModel_D2; //Unit:Gy
	double ac = MKModel_ac;  //Unit:h-1
	double tr = MKModel_tr; //Unit:h

	// *************************** Calculate y_star *****************************
	double Integrate_y = 0; 
	for(int i=0; i<fSpecLength-1;i++) 
	{
		double yi = (fBinLimit[i]+fBinLimit[i+1])/2;  
		Integrate_y += (1-exp(-pow(yi,2)/pow(y0,2)))*fhfy[i]*fBinWidth[i];
	}
	double y_star = y0*y0*Integrate_y/yF;  // Unit:keV/um


	// ***************************** Calculate z_star *****************************
	double z_star = y_star/(rho*pi*pow(rd,2)); // unit:keV*cm3/g/um3
	z_star = z_star*0.16;

	vector<double> TimeInteravl;
	vector<double> S;
	cout << "Parameters:\n"
		<<"D1, D2, (a+c), beta, tr:\n"
		<<D1 << ' ' << D2 << ' ' << ac << ' ' << beta << ' ' << tr << endl;
	for(double t=0; t<10; t+=0.5)
	{
		TimeInteravl.push_back(t);
		if(t<tr)
		{	
			//double linearquad = (-alpha0 - beta*z_star)*(D1+D2) -beta*pow(D1+D2,2);
			double linearquad = -0.237*(D1+D2) -beta*pow(D1+D2,2); //FOR DEBUGGING
			double mixed = 2.*beta*D1*D2*( 1. -  exp(-ac*t)*(1-exp(-ac*(tr-t)))/(1-exp(-2*ac*tr)));
			double s = exp(linearquad+mixed);
			S.push_back(s);
		}
		else
		{
			//double s = exp((-alpha0 - beta*z_star)*(D1+D2) -beta*(D1*D1+D2*D2));
			double s = exp(-0.237*(D1+D2) -beta*(D1*D1+D2*D2)); //FOR DEBUGGING
			S.push_back(s);
		}
	}	
	std::cout<<"************************************ Get RBE with MK method - SplitDose ***************************"<<std::endl;
	std::cout<<"Time\tSurvival\tRBE\n";
	for(int i=0; i<TimeInteravl.size(); i++)
	{
		std::cout<<TimeInteravl[i]<<'\t';
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<S[i] << endl; 
		//		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(S_var[i])<<")\t";
		//		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<RBE[i];
		//		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(RBE_var[i])<<")"<<std::endl;
	}

	cout<<"Default parameters:"<<endl;
	cout<<"1.The reference radiation is X-ray(200 kVp) with alpha = 0.19 Gy-1 and beta = 0.05 Gy-2"<<std::endl;
	cout<<"2.THe bilogical end point is 10% survival of he human salivary gland (HSG) tumor cells."<<std::endl;
	cout<<"Parameter used in this calculation:"     <<std::endl;
	cout<<"alpha0="<<MKModel_alpha0<<" Gy-1; "<<"beta="<<MKModel_beta<<" Gy-2; "
		<<"rd="<<MKModel_rd<<" um; "<< "Rn="<< MKModel_Rn <<" um; "<<"y0="<<MKModel_y0<<" keV/um; "
		<<"rho="<<MKModel_rho<<" g/cm3"  <<std::endl;
	cout<<"***************************************************************************************\n"<<std::endl;


}
//
//void TsGetSurvivalRBEQualityFactor::GetSurvWithMKModel_ContinousIrradiation()
//{
//	//Set parameters
//	double alpha0 = MKModel_alpha0; // Unit:Gy-1
//	double alphaX = MKModel_alphaX; // Unit:Gy-1
//	double beta   = MKModel_beta;   // Unit:Gy-2
//	double rd     = MKModel_rd;     // Unit:um
//	double Rn     = MKModel_rd;     //unit:um
//	double y0     = MKModel_y0;     // Unit:keV/um
//	double rho    = MKModel_rho;    // Unit:g/cm3
//	double pi     = 3.1415;
//	
//	//Split DoseIrradaiation parameters
//	double D = MKM_D; //Unit:Gy
//	double ac = MKM_ac;  //Unit:h-1
//	double IrradiationTime = MKM_T; //Unit:h
//
//	int steps = 2000.;
//	double DoseRate = D/IrradiationTime;
//	
//	// *************************** Calculate y_star *****************************
//	double Integrate_y = 0; 
//	for(int i=0; i<fSpecLength-1;i++) 
//	{
//		double yi = (fBinLimit[i]+fBinLimit[i+1])/2;  
//		Integrate_y += (1-exp(-pow(yi,2)/pow(y0,2)))*fhfy[i]*fBinWidth[i];
//	}
//	double y_star = y0*y0*Integrate_y/yF;  // Unit:keV/um
//
//
//	// ***************************** Calculate z_star *****************************
//	double z_star = y_star/(rho*pi*pow(rd,2)); // unit:keV*cm3/g/um3
//	z_star = z_star*0.16;
//
//	double alpha = alpha0+z_star*beta;
//	double aa = 0.;
//	for(int n=1; n<steps; n++)
//	{
//		double MicroIrradiationTime = n*IrradiationTime/(steps-n);
//		aa += (steps-n)*( 1 - exp(-ac*MicroIrradiationTime)*(1-exp(-2*ac*(tr-MicroIrradiationTime)))/(1-exp(-2*ac*tr))  );
//	}	
//	
//	double betaPrime = beta*(1-2.*aa/(steps*steps));
//
//	
//}



void TsGetSurvivalRBEQualityFactor::GetSurvWithDSMKModel()
{
	//Set parameters
	double alpha0 = MKModel_alpha0; // Unit:Gy-1
	double alphaX = MKModel_alphaX; // Unit:Gy-1
	double beta   = MKModel_beta;   // Unit:Gy-2
	double betaX  = MKModel_betaX;
	double rd     = MKModel_rd;     // Unit:um
	double Rn     = MKModel_Rn;     //unit:um
	double y0     = MKModel_y0;     // Unit:keV/um
	double rho    = MKModel_rho;    // Unit:g/cm3
	double pi     = 3.1415;

	TsSpecificEnergy* aSpecificEnergy_D = new TsSpecificEnergy(fyVector_Particle, rd, fGetStatitisticInfo, fSpectrumUpdateTimes);
	vector<double> zn = aSpecificEnergy_D->GetBinCenter();
	vector<double> hfz = aSpecificEnergy_D->GetHfz();
	vector<double> zBinWidth = aSpecificEnergy_D->GetBinWidth(); 

	TsSpecificEnergy* aSpecificEnergy_C = new TsSpecificEnergy(fyVector_Particle, Rn, fGetStatitisticInfo, fSpectrumUpdateTimes);

	std::vector<double> Szn, S, S_var, RBE, RBE_var;
	std::vector<std::vector<double>> S_Particle, RBE_Particle;

	//PRINT CALCULATION UPDATES
	int tenPercent = ceil((zn.size()+Doses.size())/10.);
	int update = 1;
	clock_t start,end;
	start = clock();

	std::cout<<"************************************ Get RBE with MK method - DSMKModel ***************************\n"<<std::endl;

	//CALCULATE SURVIVAL AND RBE
	double z0 = y0*0.16/(rho*pi*pow(rd,2));
	for(int i=0; i<zn.size(); i++)
	{
		vector<double> multieventDomain = aSpecificEnergy_D->GetHfzMultiEvent(zn[i],MCMultieventIterations);
		double log_szn = 0.;
		for(int j=0; j<zn.size(); j++)
		{	
			double zj = zn[j];
			double zDPrime = z0*sqrt(1-exp(-pow(zj/z0,2)));
			double zDPrime2 = z0*z0*(1-exp(-pow(zj/z0,2)));

			log_szn += -(alpha0*zDPrime*multieventDomain[j]*zBinWidth[j]) - (beta*zDPrime2*multieventDomain[j]*zBinWidth[j]);
		}
		Szn.push_back(exp(log_szn));

		if (update%tenPercent==0)
		{
			end = clock();
			float duration = (float) (end - start)/ CLOCKS_PER_SEC;
			std::cout << (update/tenPercent)*10<< " % of calculation finished, total time used :"<<duration<<" sec; ("<<duration/60<<" min)"<<std::endl;

		}
		update++;
	}


	for(double D:Doses)
	{
		std::vector<double>  multieventNucleus = aSpecificEnergy_C -> GetHfzMultiEvent(D,MCMultieventIterations); //Array di zfz Vs z //DEVE AVERE IL BINCENTER DELLE Zn
		std::vector<std::vector<double>>  multieventNucleusParticleContribution = aSpecificEnergy_C -> GetMultiEventParticleContribution();

		double Survival=0.;
		std::vector<double> scomponent(10, 0.);
		std::vector<double> rbecomponent(10,0.);

		for(int i=0; i<multieventNucleus.size(); ++i) 
		{
			Survival += multieventNucleus[i]*Szn[i]*zBinWidth[i]; //stesso del DSMKM
			for(int comp = 0; comp<10; comp++ )
			{
				scomponent[comp] += multieventNucleusParticleContribution[i][comp]*multieventNucleus[i]*Szn[i]*zBinWidth[i];
			}
		}


		//METTERE RBE = (aX*D+bX*D2)/SurvModello;	
		double rbe =0., rbe_var = 0.;
		if(D>0)
		{ 
			rbe = (sqrt( (alphaX*alphaX) - (4*betaX*log(Survival)) ) - alphaX)/(2*betaX*D);
			//rbe_var = (1./(pow(D*s,2)*(alphaX*alphaX-4*betaX*log(s))))*s_var;
			for(int comp=0; comp<10; comp++)
				rbecomponent.push_back( (sqrt( (alphaX*alphaX) - (4*betaX*log(scomponent[comp])) ) - alphaX)/(2*betaX*D)  );

		}

		S.push_back(Survival);
		RBE.push_back(rbe);

		S_Particle.push_back(scomponent);
		RBE_Particle.push_back(rbecomponent);

		S_var.push_back(0);
		RBE_var.push_back(0);
		if (update%tenPercent==0)
		{
			end = clock();
			float duration = (float) (end - start)/ CLOCKS_PER_SEC;
			std::cout << (update/tenPercent)*10<< " % of Radiobiological model update finished, total time used :"<<duration<<" sec; ("<<duration/60<<" min)"<<std::endl;

		}
		update++;
	}


	WriteMKMSurvival("MKModel_DSMKM.csv", Doses, S, S_var, RBE, RBE_var);

	if(fGetParticleContribution)
	{
		WriteSurvivivalRBEParticleContribution("MKModel_DSMKM_SurvParticle.csv", Doses, S_Particle);
		WriteSurvivivalRBEParticleContribution("MKModel_DSMKM_RBEParticle.csv", Doses, RBE_Particle);
	}

	std::cout<<"\nDoses\tSurvival\tRBE\n";
	for(int i=0; i<Doses.size(); i++)
	{
		std::cout<<Doses[i]<<'\t';
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<S[i]; 
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(S_var[i])<<")\t";
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<RBE[i];
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(RBE_var[i])<<")"<<std::endl;
	}

	cout<<"Default parameters:"<<endl;
	cout<<"1.The reference radiation is X-ray(200 kVp) with alpha = 0.19 Gy-1 and beta = 0.05 Gy-2"<<std::endl;
	cout<<"2.THe bilogical end point is 10% survival of he human salivary gland (HSG) tumor cells."<<std::endl;
	cout<<"Parameter used in this calculation:"     <<std::endl;
	cout<<"alpha0="<<MKModel_alpha0<<" Gy-1; "<<"beta="<<MKModel_beta<<" Gy-2; "
		<<"rd="<<MKModel_rd<<" um; "<< "Rn="<< MKModel_Rn <<" um; "<<"y0="<<MKModel_y0<<" keV/um; "
		<<"rho="<<MKModel_rho<<" g/cm3"  <<std::endl;
	cout<<"***************************************************************************************\n"<<std::endl;



}

void TsGetSurvivalRBEQualityFactor::GetSurvWithSMKModel()
{
	//Set parameters
	double alpha0 = MKModel_alpha0; // Unit:Gy-1
	double alphaX = MKModel_alphaX; // Unit:Gy-1
	double beta   = MKModel_beta;   // Unit:Gy-2
	double betaX  = MKModel_betaX;
	double rd     = MKModel_rd;     // Unit:um
	double Rn     = MKModel_Rn;     //unit:um
	double y0     = MKModel_y0;     // Unit:keV/um
	double rho    = MKModel_rho;    // Unit:g/cm3
	double pi     = 3.1415;

	// *************************** Calculate y_star *****************************
	double Integrate_y = 0;
	for(int i=0; i<fSpecLength-1;i++) 
	{
		double yi = (fBinLimit[i]+fBinLimit[i+1])/2;  
		Integrate_y += (1-exp(-pow(yi,2)/pow(y0,2)))*fhfy[i]*fBinWidth[i];	
	}
	double y_star = y0*y0*Integrate_y/yF;  // Unit:keV/um


	// ***************************** Calculate alphta *****************************
	double zstar = y_star/(rho*pi*pow(rd,2)); // unit:keV*cm3/g/um3
	zstar = zstar*0.16;                        // unit:J/kg = Gy

	double zD = yD/(rho*pi*pow(rd,2));
	zD *= 0.16;

	double znD = yD/(rho*pi*pow(Rn,2));
	znD *= 0.16;

	double alphaSMK = alpha0 + beta*zstar;       // unit:Gy-1
	double betaSMK = beta*(zstar/zD);  	       // Unit:Gy-2

	// ***************************** Calculate Survival and RBE *****************************
	std::vector<double> S, RBE;
	for(double D:Doses)
	{
		double s =exp( -(alphaSMK*D) - (betaSMK*D*D))*( 1 + znD*D*(-betaSMK + 0.5*pow(alphaSMK + (2*betaSMK*D) ,2)) );
		double rbe = 0;
		if(D>0) 
			rbe = (sqrt( (alphaX*alphaX) - (4*betaX*log(s)) ) - alphaX)/(2*betaX*D); 


		cout << alphaX << ' ' <<betaX << ' ' << s << ' ' <<rbe << endl;

		S.push_back(s);
		RBE.push_back(rbe);

	}

	// ************************* Calculate statistic error ************************
	double ystar_var =0.;
	double zstar_var =0.;
	double zD_var = pow(0.16/(rho*pi*pow(rd,2)),2)*yD_var;
	double znD_var = pow(0.16/(rho*pi*pow(Rn,2)),2)*yD_var;
	double alphaSMK_var = 0.;
	double betaSMK_var = 0.;
	std::vector<double> S_var, RBE_var;

	for(int i=0; i<fSpecLength-1;i++) 
	{    
		double yi = (fBinLimit[i]+fBinLimit[i+1])/2;     
		double aa = (1-exp(-pow(yi,2)/pow(y0,2)));
		ystar_var += pow(y0*y0*aa*fBinWidth[i]/yF,2)*fy_var[i] +pow(y0*y0*aa*fhfy[i]*fBinWidth[i]/(yF*yF),2)*yF_var;
	}

	zstar_var = pow(0.16/(rho*pi*pow(rd,2)),2)*ystar_var;
	alphaSMK_var = beta*beta*zstar_var;
	betaSMK_var = beta*beta*(znD_var/(zD*zD) + pow(zstar/(zD*zD),2)*zD_var);

	for(int i=0; i<Doses.size(); ++i)
	{
		double dSdznD = S[i]*( (Doses[i]*(-betaSMK + 0.5*pow(alphaSMK + (2*betaSMK*Doses[i]) ,2)))/(1 + znD*Doses[i]*(-betaSMK + 0.5*pow(alphaSMK + (2*betaSMK*Doses[i]) ,2))));
		double dSdAlpha = S[i]*(-Doses[i] + (znD*Doses[i]*(alphaSMK + (2*betaSMK*Doses[i])))/(1 + znD*Doses[i]*(-betaSMK + 0.5*pow(alphaSMK + (2*betaSMK*Doses[i]) ,2))) );
		double dSdBeta = S[i]*( -pow(Doses[i],2) + (znD*Doses[i]*(-1+2*Doses[i]*(alphaSMK+2*betaSMK*Doses[i])))/(1 + znD*Doses[i]*(-betaSMK + 0.5*pow(alphaSMK + (2*betaSMK*Doses[i]) ,2))) );
		double s_var = dSdznD*dSdznD*znD_var + dSdAlpha*dSdAlpha*alphaSMK_var + dSdBeta*dSdBeta*betaSMK_var;
		double rbe_var = (1./(pow(Doses[i]*S[i],2)*(alphaX*alphaX-4*betaX*log(S[i]))))*s_var;
		S_var.push_back(s_var);
		RBE_var.push_back(rbe_var);

	}

	WriteMKMSurvival("MKModel_SMKM.csv", Doses, S, S_var, RBE, RBE_var);
	std::cout<<"************************************ Get RBE with MK method - SMKModel ***************************"<<std::endl;
	std::cout<<"Doses\tSurvival\tRBE\n";
	for(int i=0; i<Doses.size(); i++)
	{
		std::cout<<Doses[i]<<'\t';
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<S[i];
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(S_var[i])<<")\t";
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<RBE[i];
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(RBE_var[i])<<")\t"<<std::endl;
	}

	std::cout<<"\nSpectrum Parameters:"<<endl;
	std::cout<<"yF: " << yF  <<" ( " <<sqrt(yF_var) << " )" << endl
		<<"yD: " << yD  <<" ( " <<sqrt(yD_var) << " )"<< endl
		<<"y*: " << y_star <<" ( " <<sqrt(ystar_var) << " )"  << endl
		<<"z*: " << zstar <<" ( " <<sqrt(zstar_var) << " )" << endl;
	std::cout<<"***************************************************************************************\n"<<std::endl;



}

void TsGetSurvivalRBEQualityFactor::GetRBEWithBioWeightFunction()
{
	vector< vector<double> > BioWeightFuncData;

	//******************************************************************
	//                            ReadData
	//******************************************************************
	ifstream infile;
	infile.open(BioWeightFunctionDataFile);
	double x=0; // y (kev/um)
	double y=0; //r(y)


	if (!infile)    {
		std::cout<<"Error in open "<<BioWeightFunctionDataFile<<" !"<<std::endl;
		return;
	}


	while (!infile.eof())
	{
		infile >> x >>y;

		vector<double> singleData;
		singleData.push_back(x);
		singleData.push_back(y);

		BioWeightFuncData.push_back(singleData);
	}


	//******************************************************************
	//                           Calculate RBE
	//******************************************************************
	double RBE = 0;
	double RBE_var =0;
	double RBE_std =0;
	double ry  = 1;
	int Ncomponents = fyParticleContribution[0].size();
	std::vector<double> RBEComponents(Ncomponents, 0.);
	std::vector<double> RBEComponents_var(Ncomponents, 0.);
	for(int i=0; i<fSpecLength; i++)
	{
		if(fBinLimit[i]>1000 )                // if y(keV/um) set r(y)=0
			ry=0;   
		else
		{
			for(int j=0; j<BioWeightFuncData.size()-1; j++)
			{ 
				if( fBinLimit[i]>=BioWeightFuncData[j][0] && fBinLimit[i]<BioWeightFuncData[j+1][0])
					ry = BioWeightFuncData[j][1];
			}           
		}
		RBE += ry*fhdy[i]*fBinWidth[i];
		RBE_var += pow(ry*fBinWidth[i], 2)*dy_var[i];
		for(int comp=0; comp<Ncomponents; comp++)
		{
			RBEComponents[comp] += ry*fhdy[i]*fBinWidth[i]*fyParticleContribution[i][comp];
			RBEComponents_var[comp] += pow(ry*fBinWidth[i]*fyParticleContribution[i][comp], 2)*dy_var[i];
		}

	}
	RBE_std = sqrt(RBE_var);

	std::cout<<"******************** Get RBE with biological weight function **************************"<<std::endl;
	std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<"RBE = "<<RBE;
	std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< RBE_std<<")"<<std::endl;
	std::cout<<"Default parameters:"<<endl;
	std::cout<<"1. RBE is calculated with endpoint of intestinal tolerance in mice."<<endl;
	std::cout<<"2. The value of biological weight function, r(y), was set as 0 when y> 1000 keV/um."<<endl;
	std::cout<<"***************************************************************************************\n"<<std::endl;
}

void TsGetSurvivalRBEQualityFactor::GetQualityFactorWithICRU40()
{
	//******************************************************************
	//                           Calculate Q(y)
	//******************************************************************
	double Q = 0;
	double Q_var =0;
	double Q_std =0;
	double qy  = 1;

	int Ncomponents = fyParticleContribution[0].size();
	std::vector<double> QComponents(Ncomponents, 0.);
	std::vector<double> QComponents_var(Ncomponents, 0.);
	for(int i=0; i<fSpecLength-1; i++)
	{
		double yBinCenter = (fBinLimit[i+1] + fBinLimit[i])/2.;
		double qy = ( 5510. / yBinCenter )*( 1 - exp(-(5e-5*yBinCenter*yBinCenter) - (2e-7*yBinCenter*yBinCenter*yBinCenter) ) );
		Q += qy*fhdy[i]*fBinWidth[i];
		Q_var += pow(qy*fBinWidth[i], 2)*dy_var[i];
		for(int comp=0; comp<Ncomponents; comp++)
		{
			QComponents[comp] += qy*fhdy[i]*fBinWidth[i]*fyParticleContribution[i][comp];
			QComponents_var[comp] += pow(qy*fBinWidth[i]*fyParticleContribution[i][comp], 2)*dy_var[i];
		}


	}

	Q_std = sqrt(Q_var);

	if(fGetParticleContribution)
		WriteQParticleContribution("QICRU40_Particle.csv", QComponents);

	std::cout<<"******************** Get Quality Factor with ICRU40 **************************"<<std::endl;
	std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<"Q = "<<Q;
	std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< Q_std<<")"<<std::endl;
	std::cout<<"***************************************************************************************\n"<<std::endl;
}

void TsGetSurvivalRBEQualityFactor::GetQualityFactorWithKellereHahn()
{
	//******************************************************************
	//                           Calculate Q(y)
	//******************************************************************
	double Q = 0;
	double Q_var =0;
	double Q_std =0;
	double qy  = 1;

	int Ncomponents = fyParticleContribution[0].size();
	std::vector<double> QComponents(Ncomponents, 0.);
	std::vector<double> QComponents_var(Ncomponents, 0.);
	for(int i=0; i<fSpecLength-1; i++)
	{
		double yBinCenter = (fBinLimit[i+1] + fBinLimit[i])/2.;
		double qy = 0.3* yBinCenter*pow(1+pow(yBinCenter/137.,5),-0.4);

		Q += qy*fhdy[i]*fBinWidth[i];
		Q_var += pow(qy*fBinWidth[i], 2)*dy_var[i];
		for(int comp=0; comp<Ncomponents; comp++)
		{
			QComponents[comp] += qy*fhdy[i]*fBinWidth[i]*fyParticleContribution[i][comp];
			QComponents_var[comp] += pow(qy*fBinWidth[i]*fyParticleContribution[i][comp], 2)*dy_var[i];
		}
	}
	Q_std = sqrt(Q_var);

	if(fGetParticleContribution)
	{
		cout << "MERDA1"<< endl;
		WriteQParticleContribution("QKellerer_Particle.csv", QComponents);
	}
	std::cout<<"******************** Get Quality Factor with Kellerer-Hahn approximation **************************"<<std::endl;
	std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<"Q = "<<Q;
	std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< Q_std<<")"<<std::endl;
	std::cout<<"***************************************************************************************\n"<<std::endl;
}

void TsGetSurvivalRBEQualityFactor::GetSurvWithGSM2()
{

	double alphaX = GSM2_alphaX;
	double betaX = GSM2_betaX;
	TsGSM2* aGSM2 = new TsGSM2(yF,  GSM2_rd, GSM2_Rn, GSM2_a, GSM2_b, GSM2_r, fyVector_Particle, fGetStatitisticInfo, fSpectrumUpdateTimes);
	cout << MCMultieventIterations << endl;
	vector<double> zBinCenter = aGSM2->GetZn();
	vector<double> zBinWidth = aGSM2->GetzBinWidth();
	double Nbins = zBinWidth.size();

	vector<double> Szn, Szn_var;

	//PRINT CALCULATION UPDATES
	int tenPercent = ceil((zBinCenter.size()+Doses.size())/10.);
	int update = 1;
	clock_t start,end;
	start = clock();

	//ofstream p0X("p0X.csv");
	//p0X<<"zn,0,1,2,3,..."<<endl; //DEBUG
	
	//ofstream p0Y("p0Y.csv");
	//p0Y<<"zn,0,1,2,3,..."<<endl; //DEBUG
	
	std::cout<<"************************************ Get RBE with GSM2 method ***************************\n"<<std::endl;
	for(double zn:zBinCenter)
	{
		vector<vector<double>> p0xy = aGSM2->GetInitialLethalNonLethalDamages(zn, MCMultieventIterations);
		std::vector<double> sn_snVar = aGSM2-> GetSurvivalDomain(GSM2_a, GSM2_b, GSM2_r, p0xy[0], p0xy[1], p0xy[2], p0xy[3]);
		Szn.push_back( sn_snVar[0] ); //array di S(Zn) per ogni zN
		Szn_var.push_back(sn_snVar[1]);

	/*	p0X << zn;
		p0Y << zn;
		for(int x=0; x<p0xy[0].size(); x++)
		{
			p0X<<','<<p0xy[0][x];
			p0Y<<','<<p0xy[1][x];

		}
		p0X << endl;
		p0Y << endl; */
		//cout << zn <<'\t' << sn_snVar[0] << '\t' << sqrt(sn_snVar[1]) << endl;
		if (update%tenPercent==0)
		{
			end = clock();
			float duration = (float) (end - start)/ CLOCKS_PER_SEC;
			std::cout << (update/tenPercent)*10<< " % of Radiobiological model update finished, total time used :"<<duration<<" sec; ("<<duration/60<<" min)"<<std::endl;

		}	
		update++;
	}

//	p0X.close();
//	p0Y.close();


//	ofstream fn_Nucleus("fn_Nucleus.csv"); //DEBUG
//	fn_Nucleus<<"Dose,z,fn"<<endl; //DEBUG

	int Ndomains = std::floor( pow(GSM2_Rn/GSM2_rd,3) );
	//integrro S(zn)
	//Numero domini = rapporto dei raggi al quadrato.
	vector<double> S, S_var, RBE, RBE_var;
	std::vector<std::vector<double>> S_Particle, RBE_Particle;
	for(double D:Doses)
	{
		std::vector<double>  multieventNucleus = aGSM2 -> GetMultieventNucleus(D, MCMultieventIterations); //Array di zfz Vs z //DEVE AVERE IL BINCENTER DELLE Zn
		std::vector<double> multieventNucleus_var = aGSM2 -> GetMultieventNucleusVariance();
		std::vector<std::vector<double>> multieventNucleusParticleContribution = aGSM2 -> GetMultieventNucleusParticleContribution();

		double s=0., s_var =0., tmp = 0., tmp1 = 0.;	
		std::vector<double> scomponent(10, 0.);
		std::vector<double> rbecomponent(10,0.);
		for(int i=0; i<Nbins; ++i) 
		{
//			fn_Nucleus <<D<<','<< zBinCenter[i] << ',' << multieventNucleus[i]*zBinWidth[i] <<  endl; //DEBUG
			s += multieventNucleus[i]*pow(Szn[i],Ndomains)*zBinWidth[i]; //stesso del DSMKM
			s_var += pow(Ndomains*pow(Szn[i], Ndomains-1)*multieventNucleus[i]*zBinWidth[i],2)*Szn_var[i] + pow(pow(Szn[i],Ndomains)*zBinWidth[i],2)*multieventNucleus_var[i];
			tmp += pow(pow(Szn[i], Ndomains-1)*multieventNucleus[i]*zBinWidth[i],2)*Szn_var[i];
			tmp1 += pow(pow(Szn[i],Ndomains)*zBinWidth[i],2)*multieventNucleus_var[i];

			for(int comp = 0; comp<10; comp++ )
			{
				scomponent[comp] += multieventNucleusParticleContribution[i][comp]*multieventNucleus[i]*pow(Szn[i],Ndomains)*zBinWidth[i];
			}
		}

		double rbe = 0., rbe_var = 0.;
		if(D>0)
		{ 
			rbe = (sqrt( (alphaX*alphaX) - (4*betaX*log(s)) ) - alphaX)/(2*betaX*D);
			rbe_var = (1./(pow(D*s,2)*(alphaX*alphaX-4*betaX*log(s))))*s_var;

			for(int comp=0; comp<10; comp++)
				rbecomponent.push_back( (sqrt( (alphaX*alphaX) - (4*betaX*log(scomponent[comp])) ) - alphaX)/(2*betaX*D)  );

		}


		S.push_back(s);
		RBE.push_back(rbe);


		S_Particle.push_back(scomponent);
		RBE_Particle.push_back(rbecomponent);

		S_var.push_back(s_var);
		RBE_var.push_back(rbe_var);

		if (update%tenPercent==0)
		{
			end = clock();
			float duration = (float) (end - start)/ CLOCKS_PER_SEC;
			std::cout << (update/tenPercent)*10<< " % of Radiobiological model update finished, total time used :"<<duration<<" sec; ("<<duration/60<<" min)"<<std::endl;

		}	
		update++;
	}

//	fn_Nucleus.close(); //DEBUG


	WriteGSM2Survival("GSM2.csv", Doses, S, S_var, RBE, RBE_var);

	if(fGetParticleContribution)
	{
		WriteSurvivivalRBEParticleContribution("GSM2_SurvParticle.csv", Doses, S_Particle);
		WriteSurvivivalRBEParticleContribution("GSM2_RBEParticle.csv", Doses, RBE_Particle);
	}

	std::cout<<"Doses\tSurvival\tRBE\n";
	for(int i=0; i<Doses.size(); i++)
	{
		std::cout<<Doses[i]<<'\t';
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<S[i];
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(S_var[i])<<")\t" << S_Particle[i][9] << '\t';
		std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<RBE[i];
		std::cout<<setiosflags(ios::fixed)<<setprecision(6)<<" ("<< sqrt(RBE_var[i])<<")\t" << RBE_Particle[i][9]<<std::endl;
	}
	std::cout<<"Default parameters:"<<endl;
	std::cout<<"1.The reference radiation is X-ray(200 kVp) with alpha = 0.19 Gy-1 and beta = 0.05 Gy-2"<<std::endl;
	std::cout<<"2.THe bilogical end point is 10% survival of he human salivary gland (HSG) tumor cells."<<std::endl;
	std::cout<<"Parameter used in this calculation:"     <<std::endl;
	// std::cout<<"alpha0="<<MKModel_alpha0<<" Gy-1; "<<"beta="<<MKModel_beta<<" Gy-2; "
	//         <<"rd="<<MKModel_rd<<" um; "<< "Rn="<< MKModel_Rn <<" um; "<<"y0="<<MKModel_y0<<" keV/um; "
	//         <<"rho="<<MKModel_rho<<" g/cm3"  <<std::endl;
	std::cout<<"***************************************************************************************\n"<<std::endl;




}

void TsGetSurvivalRBEQualityFactor::WriteMKMSurvival(string filename, std::vector<double> D, std::vector<double> S, std::vector<double> Svar, std::vector<double> RBE, std::vector<double> RBEvar)
{
	std::ofstream output(filename);
	output << "# MKM Parameters\n#\n";
	output << "# Alpha0 = " << MKModel_alpha0 << " Gy-1\n"
		<< "# Beta = " << MKModel_beta << " Gy-2\n"
		<< "# AlphaX = " << MKModel_alphaX << " Gy-1 Reference radiation\n"
		<< "# Domain Radius = " << MKModel_rd << " um\n"
		<< "# Nucleus Radius = " << MKModel_Rn << " um\n"
		<< "#\n";

	output << "# Dose[Gy], Survival, SurvivalStd, RBE, RBEStd\n";
	for(int i=0; i<D.size(); i++)
		output << std::fixed << std::setprecision(7) << D[i] << ", " << S[i] << ", " << sqrt(Svar[i]) << ", " << RBE[i] << ", " << sqrt(RBEvar[i]) << endl;

	output.close();
	std::cout << "Output file " << filename << " written!" <<std::endl;

}

void TsGetSurvivalRBEQualityFactor::WriteGSM2Survival(string filename, std::vector<double> D, std::vector<double> S, std::vector<double> Svar, std::vector<double> RBE, std::vector<double> RBEvar)
{
	std::ofstream output(filename);
	output << "# GSM2 Parameters\n#\n";
	output << "# kappa = " << GSM2_kappa << " Gy-1\n"
		<< "# lambda = " << GSM2_lambda << " Gy-1\n"
		<< "# a = " << GSM2_a << "??\n"
		<< "# b = " << GSM2_b << "??\n"
		<< "# r = " << GSM2_r << "??\n"
		<< "# AlphaX = " << GSM2_alphaX << " Gy-1 Reference radiation\n"
		<< "# BetaX = " << GSM2_betaX << "Gy-2 Reference radiation\n"
		<< "# Domain Radius = " << GSM2_rd << " um\n"
		<< "# Nucleus Radius = " << GSM2_Rn << " um\n"
		<< "#\n";

	output << "# Dose[Gy], Survival, SurvivalStd, RBE, RBEStd\n";
	for(int i=0; i<D.size(); i++)
		output << std::fixed << std::setprecision(7) << D[i] << ", " << S[i] << ", " << sqrt(Svar[i]) << ", " << RBE[i] << ", " << sqrt(RBEvar[i]) << endl;

	output.close();
	std::cout << "Output file " << filename << " written!" <<std::endl;




}

void TsGetSurvivalRBEQualityFactor::WriteSurvivivalRBEParticleContribution(string filename, std::vector<double> D, std::vector<std::vector<double>> Vector_Particle)
{
	std::ofstream outputParticle(filename);

	outputParticle << "D(Gy)    e-       Hprim    Hsec    He      Li      Be      B        C       Other    Total[Survival]\n";

	for (int i=0;i<D.size();i++)
	{       
		outputParticle << std::fixed << std::setprecision(7) << D[i] <<"  ";

		for(int loop=0; loop<10; loop++)
			outputParticle << std::fixed << std::setprecision(7) << Vector_Particle[i][loop] <<"  ";

		outputParticle << std::endl;

	}
	outputParticle.close();
}

void TsGetSurvivalRBEQualityFactor::WriteQParticleContribution(string filename, std::vector<double> Vector_Particle)
{

	cout <<"MERDA"<< endl;
	std::ofstream outputParticle(filename);

	outputParticle << "e-       Hprim    Hsec    He      Li      Be      B        C       Other    Total[Q]\n"; 
	for(int loop=0; loop<10; loop++)
		outputParticle << std::fixed << std::setprecision(7) << Vector_Particle[loop] <<"  ";

	outputParticle << std::endl;

	outputParticle.close();
}
