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
#include<iostream>
#include<vector>
#include<cstring>
#include<cmath>
#include<fstream>
#include<algorithm>
#include<random>
//#include<filesystem>
#include <chrono>
#include "TsGetSurvivalRBEQualityFactor.hh"
#include "TsLinealEnergy.hh"
#include "TsSpecificEnergy.hh"
using namespace std;

int main(int argc, char *argv[])
{

	/////////////////////////////////////////////////////////////
	//
	// INITIALIZE VARIABLES
	//
	////////////////////////////////////////////////////////////
	bool fGetStatitisticInfo = false;
	int fSpectrumUpdateTimes = 100;
	int fSetMultiEventStatistic = 100000;
	bool fGetParticleContribution = false;
	//Initialize Microdosimetric variables
	double yF=0., yD=0.;
	double yF_var=0., yF_std=0.;
	double yD_var=0., yD_std=0.;
	vector<double> hfy, hdy, hyfy, hydy, BinLimit, BinWidth;
	std::vector<double> fy_var, ydy_var, yfy_var, dy_var;
	double MKModel_alpha0 = 0.13; // Unit:Gy-1
	double MKModel_beta   = 0.05; // Unit:Gy-2,
	double MKModel_alphaX = 0.19; // Unit:Gy-1
	double MKModel_betaX = 0.05;
	double MKModel_rd     = 0.42; // Unit:um
	double MKModel_Rn     = 6; //Unit:um
	double MKModel_rho    = 1.;    // Unit:g/cm3
	double MKModel_y0     = 150;  // Unit:keV/um

	//Split Dose MKM
	double MKModel_D1 = 1; //Unit:Gy
	double MKModel_D2 = 1; //Unit:Gy
	double MKModel_ac = 2.187; //Unit:h-1
	double MKModel_tr = 2.284; //Unit:h

	//GSM2 default parameters
	double GSM2_kappa      = 0.5;
	double GSM2_lambda     = 0.005;
	double GSM2_rd         = 0.42; //Unit: um
	double GSM2_Rn         = 6;    //Unit: um
	double GSM2_zStart     = 0.01; //Gy
	double GSM2_zEnd       = 100; //Gy
	double GSM2_zBins      = 100; //Bins of spectra
	double GSM2_a          = 0.1;
	double GSM2_b          = 0.1;
	double GSM2_r          = 0.1;
	double GSM2_alphaX 	= 0.19;
	double GSM2_betaX 	= 0.05;
	//Macroscopic Doses
	vector<double> Doses = {0,10,0.5}; //Unit:Gy
	bool MKMSatCorrFlag = 0, MKMnonPoissFlag = 0, SMKMFlag = 0, DSMKMFlag = 0, GSM2Flag = 0, RBEWeightingFlag = 0, QfICRUFlag = 0, QfKellFlag = 0;
	bool H460Flag = 0, H1437Flag = 0;
	/////////////////////////////////////////////////////////////
	//
	// READ INPUT PARAMETERS
	//
	////////////////////////////////////////////////////////////

	string TopasScorerFile;
	string TopasScorerFile_Nucleus; // new line
	string BioWeightFunctionDataFile = "BioWeightFuncData_interpolation.txt";
	//Loop on inputs
	for(int i=0; i<argc; i++)
	{
		if(strcmp(argv[i],"-MKM_rDomain") == 0) {MKModel_rd = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_rNucleus") == 0) {MKModel_Rn = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_alpha") == 0) {MKModel_alpha0 = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_beta") == 0) {MKModel_beta = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_alphaX") == 0) {MKModel_alphaX = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_betaX") == 0) {MKModel_betaX = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_y0") == 0) {MKModel_y0 = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_D1") == 0) {MKModel_D1 = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_D2") == 0) {MKModel_D2 = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_ac") == 0) {MKModel_ac = stod(argv[i+1]);}
		if(strcmp(argv[i],"-MKM_tr") == 0) {MKModel_tr= stod(argv[i+1]);}

		// PARAMETERS FROM CELL LINE FLAG (TO ADD, ALSO SAME alphaX, betaX FOR BOTH GSM2 AND MKM)
		if(strcmp(argv[i],"-H460") == 0) {H460Flag = 1;}
		if(strcmp(argv[i],"-H1437") == 0) {H1437Flag = 1;}
		
		if(strcmp(argv[i],"-GSM2_rDomain") == 0) {GSM2_rd = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_rNucleus") == 0) {GSM2_Rn = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_kappa") == 0) {GSM2_kappa = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_a") == 0) {GSM2_a = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_b") == 0) {GSM2_b = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_r") == 0) {GSM2_r = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_alphaX") == 0) {GSM2_alphaX = stod(argv[i+1]);}
		if(strcmp(argv[i],"-GSM2_betaX") == 0) {GSM2_betaX = stod(argv[i+1]);}		

		if(strcmp(argv[i],"-MKMSatCorr") == 0) {MKMSatCorrFlag = 1;}
		if(strcmp(argv[i],"-MKMnonPoiss") == 0) {MKMnonPoissFlag = 1;}
		if(strcmp(argv[i],"-SMKM") == 0) {SMKMFlag = 1;}
		if(strcmp(argv[i],"-DSMKM") == 0) {DSMKMFlag = 1;}
		if(strcmp(argv[i],"-GSM2") == 0) {GSM2Flag = 1;}

		if(strcmp(argv[i],"-fGetStatitisticInfo") == 0) {fGetStatitisticInfo = true;}
		if(strcmp(argv[i],"-fSpectrumUpdateTimes") == 0) {fSpectrumUpdateTimes = stoi(argv[i+1]);}
		if(strcmp(argv[i],"-fSetMultiEventStatistic") == 0) {fSetMultiEventStatistic = stoi(argv[i+1]);}
		if(strcmp(argv[i],"-fGetParticleContribution") == 0) {fGetParticleContribution = true;}
		if(strcmp(argv[i], "-RBEWeighting") == 0) 
		{	
			BioWeightFunctionDataFile = argv[i+1];
			RBEWeightingFlag = 1;
		}
		if(strcmp(argv[i], "-QfICRU") == 0) {QfICRUFlag = 1;}
		if(strcmp(argv[i], "-QfKeller") == 0) {QfKellFlag = 1;}
		if(strcmp(argv[i],"-Doses") == 0) {
			Doses.clear();
			Doses = {stod(argv[i+1]), stod(argv[i+2]), stod(argv[i+3])};
		}

		if(strcmp(argv[i],"-topasScorer_08") == 0) {TopasScorerFile = argv[i+1];} //input file
		if(strcmp(argv[i],"-topasScorer_8") == 0) {TopasScorerFile_Nucleus = argv[i+1];} //input file // new line
		if(strcmp(argv[i],"-help") == 0) 
		{
			cout 	<<"-Rd: Domain radius [um]" <<endl
				<<"-Rc: Cell Nucleus radius [um]" <<endl
				<<"-CellLine: replace with cell line name (-H460 or -H1437)" <<endl
				<<"-topasScorer: path/file.phsp with y values" <<endl
				<<"-Doses: (initial dose value) (final dose value) (step value)" <<endl
				<<"-help: list of definitions and inputs" <<endl;
			return 0;			     
		}
	}
	
	// Parameters for each cell line (H460, H1437)
	if(H460Flag)
	{
		cout 	<<"\nH460 CELL LINE IS SELECTED \n" <<endl;
		GSM2_alphaX = 0.29;
		GSM2_betaX = 0.083;
		MKModel_alphaX = 0.29;
		MKModel_betaX = 0.083;
		GSM2_rd = 0.8;
		GSM2_Rn = 6;
		GSM2_a = 0.000899;
		GSM2_b = 0.0642;
		GSM2_r = 2.71;
	}else if(H1437Flag){
		cout 	<<"\nH1437 CELL LINE IS SELECTED \n" <<endl;
		GSM2_alphaX = 0.05;
		GSM2_betaX = 0.041;
		MKModel_alphaX = 0.05;
		MKModel_betaX = 0.041;
		GSM2_rd = 0.6;
		GSM2_Rn = 8;
		GSM2_a = 0.0148;
		GSM2_b = 0.0149;
		GSM2_r = 2.70;
	}else{
		// Default values
		cout 	<<"\nNO SPECIFIC CELL LINE IS SELECTED (USING DEFAULT VALUES) \n" <<endl;
	}

	/////////////////////////////////////////////////////////////
	//
	// READ TOPAS SCORER
	//
	////////////////////////////////////////////////////////////

	vector<vector<double>> yVector_Particle;
	vector<double>  yVector;

    	ifstream infile(&TopasScorerFile[0]);
	if(infile.fail()) // checks to see if file opended 
	{
		cout << "ERROR::" <<TopasScorerFile <<" NOT FOUND!!!" << endl;
		return -1; // no point continuing if the file didn't open...
	}
	while(!infile.eof()) // reads file to end of *file*, not line
	{ 

		double y_total, y_z0, y_z1_prim, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_;
		//double y_total, y_z0, y_z1_prim, y_z1_seco, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_;
		infile >> y_total 
			>> y_z0
			>> y_z1_prim
			//>> y_z1_seco
			>> y_z2
			>> y_z3
			>> y_z4
			>> y_z5
			>> y_z6
			>> y_z_;

		vector<double> yParticle {y_z0, y_z1_prim, 0, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_, y_total};
		
		//vector<double> yParticle {y_z0, y_z1_prim, y_z1_seco, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_, y_total};
		
		yVector_Particle.push_back(yParticle);
		yVector.push_back(y_total);  
	}
	infile.close();
	
	
	/////////////////////////////////////////////////////////////
	//
	// READ TOPAS SCORER FOR NUCLEUS // new line
	//
	////////////////////////////////////////////////////////////

	vector<vector<double>> yVector_Particle_Nucleus;
	vector<double>  yVector_Nucleus;

    	ifstream infile_Nucleus(&TopasScorerFile_Nucleus[0]);
	if(infile_Nucleus.fail()) // checks to see if file opended 
	{
		cout << "ERROR::" <<TopasScorerFile_Nucleus <<" NOT FOUND!!!" << endl;
		return -1; // no point continuing if the file didn't open...
	}
	while(!infile_Nucleus.eof()) // reads file to end of *file*, not line
	{ 

		double y_total, y_z0, y_z1_prim, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_;
		//double y_total, y_z0, y_z1_prim, y_z1_seco, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_;
		infile_Nucleus >> y_total 
			>> y_z0
			>> y_z1_prim
			//>> y_z1_seco
			>> y_z2
			>> y_z3
			>> y_z4
			>> y_z5
			>> y_z6
			>> y_z_;

		vector<double> yParticle {y_z0, y_z1_prim, 0, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_, y_total};
		
		//vector<double> yParticle {y_z0, y_z1_prim, y_z1_seco, y_z2, y_z3, y_z4, y_z5, y_z6, y_z_, y_total};
		
		yVector_Particle_Nucleus.push_back(yParticle);
		yVector_Nucleus.push_back(y_total);  
	}
	infile_Nucleus.close();


	/////////////////////////////////////////////////////////////
	//
	// GET MICRODOSIMETRIC SPECTRA AND QUANTITIES
	//
	////////////////////////////////////////////////////////////

	//TsLinealEnergy* aLinealEnergy = new TsLinealEnergy(yVector, yVector_Particle);
	TsLinealEnergy* aLinealEnergy = new TsLinealEnergy(yVector_Nucleus, yVector_Particle_Nucleus);
	BinLimit = aLinealEnergy->GetyBinLimit();
	BinWidth = aLinealEnergy->GetyBinWidth();
	hfy = aLinealEnergy -> Getfy();
	hdy = aLinealEnergy -> Getdy();
	yF  = aLinealEnergy -> GetyF();
	yD  = aLinealEnergy -> GetyD();
	fy_var = aLinealEnergy -> GetfyVariance();
	dy_var = aLinealEnergy -> GetdyVariance();
	yF_var = aLinealEnergy -> GetyFvar();
	yD_var = aLinealEnergy -> GetyDvar();
	int yBinNum = BinWidth.size();
	std::cout<<"Parameters:\n";
	std::cout<<"yF: " <<yF << "-+ " <<sqrt(yF_var) <<endl
		<<"yD: " << yD << "-+ " <<sqrt(yD_var) << endl;

	vector<vector<double>> contribution = aLinealEnergy -> GetParticleContribution();
	hydy = aLinealEnergy->Getydy();

	TsGetSurvivalRBEQualityFactor *aSurvRBEQf = new TsGetSurvivalRBEQualityFactor(contribution, yVector, yVector_Particle, yVector_Nucleus, yVector_Particle_Nucleus, &BinLimit[0], &BinWidth[0], &hfy[0], &hdy[0], yF, yD, yF_var, yD_var, fy_var, dy_var, yBinNum, fGetStatitisticInfo, fSpectrumUpdateTimes, fGetParticleContribution);
	
	aSurvRBEQf->SetDosesMacro(&Doses[0]);
	aSurvRBEQf -> SetMCMultieventIterations(fSetMultiEventStatistic);
	aSurvRBEQf->SetMKModel_alpha0(MKModel_alpha0);
	aSurvRBEQf->SetMKModel_beta(MKModel_beta);
	aSurvRBEQf->SetMKModel_alphaX(MKModel_alphaX);
	aSurvRBEQf->SetMKModel_betaX(MKModel_betaX);
	aSurvRBEQf->SetMKModel_rd(MKModel_rd);
	aSurvRBEQf->SetMKModel_Rn(MKModel_Rn);

	TsSpecificEnergy* aSpecEne = new TsSpecificEnergy(yVector_Particle, MKModel_rd, fGetStatitisticInfo, fSpectrumUpdateTimes);
	std::vector<double> fz = aSpecEne -> GetHfz();

	//std::vector<std::vector<double>> fzParticle = aSpecEne -> GetParticleContribution();
	//std::vector<double> fz1 = aSpecEne ->GetHfzMultiEvent(0.5, 1e6);
	//std::vector<double> fz10 = aSpecEne ->GetHfzMultiEvent(1, 1e6);
	//std::vector<std::vector<double>> fz10Particle = aSpecEne -> GetMultiEventParticleContribution();
	//std::vector<double> fz100 = aSpecEne ->GetHfzMultiEvent(15, 1e6);
	
	/*for(int i=0; i<hfy.size(); i++)
	{
		cout << zbin[i] <<',' <<hfy[i]; //<<',' <<fz1[i]; //<<','<<fz10[i]<<','<<fz100[i];
		//for(int j=0; j<10; j++)
		//	cout <<','<<fz10Particle[i][j]*fz10[i]; //<<',' <<fz1[i]<<','<<fz10[i]<<','<<fz100[i] <<endl;
		cout << endl;
	}
	std::vector<double> zbin = aSpecEne -> GetBinCenter();
	for(int i=0; i<fz.size(); i++)
	{
		cout << zbin[i] <<',' <<fz[i]; //<<',' <<fz1[i]; //<<','<<fz10[i]<<','<<fz100[i];
		//for(int j=0; j<10; j++)
		//	cout <<','<<fz10Particle[i][j]*fz10[i]; //<<',' <<fz1[i]<<','<<fz10[i]<<','<<fz100[i] <<endl;
		cout << endl;
	}
	*/
	
	if(MKMSatCorrFlag)
		aSurvRBEQf->GetSurvWithMKModel_SaturationCorr();
	if(MKMnonPoissFlag)
		aSurvRBEQf->GetSurvWithMKModel_nonPoissonCorr();
	if(SMKMFlag)
		aSurvRBEQf->GetSurvWithSMKModel();
	if(DSMKMFlag)
		aSurvRBEQf->GetSurvWithDSMKModel();

	if(RBEWeightingFlag)
	{
		aSurvRBEQf -> SetBioWeightFunctionDataFile(BioWeightFunctionDataFile);
		aSurvRBEQf -> GetRBEWithBioWeightFunction();
	}
	if(QfICRUFlag)
		aSurvRBEQf -> GetQualityFactorWithICRU40();
	if(QfKellFlag)
		aSurvRBEQf -> GetQualityFactorWithKellereHahn();
	
	//aSurvRBEQf->GetSurvWithMKModel_SplitDoseIrradiation();
	if(GSM2Flag)
	{
		aSurvRBEQf->SetGSM2_alphaX(GSM2_alphaX);
		aSurvRBEQf->SetGSM2_betaX(GSM2_betaX);
		aSurvRBEQf->SetGSM2_rd(GSM2_rd);
		aSurvRBEQf->SetGSM2_Rn(GSM2_Rn);
		aSurvRBEQf->SetGSM2_a(GSM2_a);
		aSurvRBEQf->SetGSM2_b(GSM2_b);
		aSurvRBEQf->SetGSM2_r(GSM2_r);
		
		aSurvRBEQf->GetSurvWithGSM2();
	}	
	
//	aSurvRBEQf->GetSurvWithGSM2();	
	return 0;
}
