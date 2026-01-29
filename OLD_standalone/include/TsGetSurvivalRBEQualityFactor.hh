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
// Claculate RBE 
// Author: Giorgio Cartechini
// Date:

#ifndef TsGetSurvivalRBEQualityFactor_hh
#define TsGetSurvivalRBEQualityFactor_hh

#include <stdint.h>
#include <vector>
#include <string>


using namespace std;

class TsGetSurvivalRBEQualityFactor 
{
	public:
		TsGetSurvivalRBEQualityFactor(std::vector<std::vector<double>> yParticleContribution, std::vector<double> yVector, std::vector<std::vector<double>> yVector_Particle, std::vector<double> yVector_Nucleus, std::vector<std::vector<double>> yVector_Particle_Nucleus, double*BinLimit, double* hBinWidth, double* hfy, double* hdy, double hyF, double hyD, double hyF_var, double hyD_var, std::vector<double> hfy_var, std::vector<double>hdy_var, int SpecLength, bool GetStatisticInfo, int SpectrumUpdateTimesi, bool GetParticleContribution);
		~TsGetSurvivalRBEQualityFactor();

		void GetSurvWithMKModel_SaturationCorr();
		void GetSurvWithMKModel_nonPoissonCorr();
	        void GetSurvWithMKModel_SplitDoseIrradiation();
		void GetSurvWithDSMKModel();
		void GetSurvWithSMKModel();

		void GetRBEWithBioWeightFunction();
		void GetQualityFactorWithICRU40();
		void GetQualityFactorWithKellereHahn();
	
		void GetSurvWithGSM2();

		void WriteMKMSurvival(string filename, std::vector<double> D, std::vector<double> S, std::vector<double> Svar, std::vector<double> RBE, std::vector<double> RBEvar);
		void WriteGSM2Survival(string filename, std::vector<double> D, std::vector<double> S, std::vector<double> Svar, std::vector<double> RBE, std::vector<double> RBEvar);
		void WriteSurvivivalRBEParticleContribution(string filename, std::vector<double> D, std::vector<std::vector<double>> Vector_Particle);
		void WriteQParticleContribution(string filename, std::vector<double> Vector_Particle);
		void Write_yD_RBE10(string filename, double yD, double Dose10, double RBE10);
		
		// TO DO: LQ/linear fit functions
		vector<double> logTransform(const vector<double>& S);
		void quadraticFit(const vector<double>& doses, const vector<double>& logS, double& alpha, double& beta, double& error);
		// Quadratic fit with GSL library
		void fit_quadratic_GSL(const std::vector<double>& doses, const std::vector<double>& logS, double& alpha, double& beta, double& error);
		void linearFit(const vector<double>& doses, const vector<double>& logS, double& alpha, double& error);
		void linearFitDoseSquare(const vector<double>& doses, const vector<double>& logS, double& beta, double& error);
		double calculateDose(double alpha, double beta, double targetS);

		void SetDosesMacro(double* vec) 
		{ 
			Doses.clear();
			for(double d = vec[0]; d <= vec[1]; d+=vec[2])
				Doses.push_back(d);
 		}
		
		//Set MKM parameters
		void SetMKModel_alpha0 	(double num ){ MKModel_alpha0 = num; }
		void SetMKModel_alphaX 	(double num ){ MKModel_alphaX = num; }
		void SetMKModel_betaX   (double num ){ MKModel_betaX  = num; }
		void SetMKModel_beta   	(double num ){ MKModel_beta   = num; }
		void SetMKModel_rho    	(double num ){ MKModel_rho    = num; }
		void SetMKModel_rd    	(double num ){ MKModel_rd     = num; }
		void SetMKModel_Rn    	(double num ){ MKModel_Rn     = num; }
		void SetMKModel_y0      (double num ){ MKModel_y0     = num; } 
		void SetMKModel_D1    	(double num ){ MKModel_D1     = num; }
		void SetMKModel_D2    	(double num ){ MKModel_D2     = num; }
		void SetMKModel_ac    	(double num ){ MKModel_ac     = num; }
		void SetMKModel_tr    	(double num ){ MKModel_tr     = num; }

		//Set GSM2 parameters
		void SetGSM2_rd 	(double num) { GSM2_rd 		= num; }
		void SetGSM2_Rn 	(double num) { GSM2_Rn 		= num; }
		void SetGSM2_a 		(double num) { GSM2_a 		= num; }
		void SetGSM2_b 		(double num) { GSM2_b 		= num; }
		void SetGSM2_r          (double num) { GSM2_r 		= num; }
		void SetGSM2_alphaX     (double num) { GSM2_alphaX 	= num; }
		void SetGSM2_betaX 	(double num) { GSM2_betaX 	= num; }
		void SetBioWeightFunctionDataFile(string fileName ){ BioWeightFunctionDataFile= fileName; }

		//MultieventIterations
		void SetMCMultieventIterations (int num) { MCMultieventIterations = num; }
	private:
		
		bool fGetStatitisticInfo;
		int  fSpectrumUpdateTimes;
		bool fGetParticleContribution;
		int MCMultieventIterations;
		std::vector<double> fyVector;
		std::vector<std::vector<double>> fyVector_Particle;
                std::vector<std::vector<double>> fyParticleContribution;
                std::vector<double> fyVector_Nucleus;
		std::vector<std::vector<double>> fyVector_Particle_Nucleus;
		double* fBinLimit;
		double* fBinWidth;
		double* fhy;
		double* fhfy;
		double* fhdy;
		double* fhydy;
		double  yF;
		double  yF_var;
		double  yD;
		double  yD_var;
		std::vector<double> fy_var;
		std::vector<double> dy_var;

		int fSpecLength;
		std::vector<double> Doses;

		double MKModel_alpha0,MKModel_alphaX, MKModel_beta, MKModel_betaX, MKModel_rho, MKModel_rd, MKModel_Rn, MKModel_y0, MKModel_D1, MKModel_D2, MKModel_ac, MKModel_tr;
		double GSM2_kappa, GSM2_lambda, GSM2_rd, GSM2_Rn, GSM2_a, GSM2_b, GSM2_r, GSM2_zStart, GSM2_zEnd, GSM2_zBins, GSM2_alphaX, GSM2_betaX;

		string BioWeightFunctionDataFile;

};



#endif
