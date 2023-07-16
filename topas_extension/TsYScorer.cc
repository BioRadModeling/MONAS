// Scorer for TsYScorer

// input
//
// mandatory parameters
// 1. GeometryNumber: 0 for spherical TEPC, 1 for cylindrical mini-TEPC and 2 for silicon microdosimeter
// 2. TransX: x position of sensitive volume
// 3. TransY: y position of sensitive volume
// 4. TransZ: z position of sensitive volume
// 5. SensitiveVolumeRadius: radius of sensitive volume (invalid for silicon microdosimeter)
// 6.1. TissueEquivalentRadius: radius of tissue equivalent volume, for spherical and cylindrical TEPC
// 6.2. MeanPathLength : Mean path lengt for SOD detecor
//
// optional parameters
// 7. LinealEnergyLowerlimit: lower threshold of y scorer
// 8. LinealEnergyUpperlimit: upper threshold of y scorer
// 9. IncludeFrequencyMeanLinealEnergy: (boolean) whether to output yF
// 10. IncludeDoseMeanLinealEnergy: (boolean) whether to output yD
// 11. GetRBEWithBioWeightFunction: true or false(default),to calculate RBE with biological weight function with endpoint of intestinal tolerance in mice
// 12. GetRBEWithMKModel : true or false(default), calculate RBE with MKM model, Kase, Y., et al. (2006).
//
// output
// y: lineal energy (keV/um) saved in ntuple
// y-yd(y): microdimetric spectrum saved in ydy.txt
// ySpecfile: all data of y, f(y), yf(y), d(y), yd(y) saved in ySpecfile.txt
//
// Last edit: 2019.04.03
//GIORGIO EDIT: 2021.03.08
//INSERTED A NEW PARTICLE: DIVISION BETWEEN PRIM AND SECO PROTONS

#include "TsYScorer.hh"
#include "TsTrackerHit.hh"
#include "TsGetSurvivalRBEQualityFactor.hh"
#include "TsSpecificEnergy.hh"

#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "g4root.hh"

#include "G4GeometryTolerance.hh"
#include "G4SDManager.hh"

#include "string.h"
#include "math.h"
#include <vector>


TsYScorer::TsYScorer(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
		G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
	: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)//,n_count(0)
{
	InitializeMicrodosimetricSpectrum();
	yVector.clear();
	yVector_Particle.clear();

	fHitsCollection = new TsTrackerHitsCollection();
	fNtuple->RegisterColumnD(&fy,    "y_total (keV/um)","");
	fNtuple->RegisterColumnD(&fy_z0, "y_z0 (keV/um)","");
	fNtuple->RegisterColumnD(&fy_z_Prim, "y_z_Prim (keV/um)",""); //DEBUG
	fNtuple->RegisterColumnD(&fy_z1_Seco, "y_z1_Seco (keV/um)",""); //DEBUG
	fNtuple->RegisterColumnD(&fy_z2, "y_z2 (keV/um)","");
	fNtuple->RegisterColumnD(&fy_z3, "y_z3 (keV/um)","");
	fNtuple->RegisterColumnD(&fy_z4, "y_z4 (keV/um)","");
	fNtuple->RegisterColumnD(&fy_z5, "y_z5 (keV/um)","");
	fNtuple->RegisterColumnD(&fy_z6, "y_z6 (keV/um)","");
	fNtuple->RegisterColumnD(&fy_z_, "y_z_ (keV/um)","");

	G4cout << "Computed tolerance = "
		<< G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/nm
		<< " nm" << G4endl;

	// *************************** mandatory parameters ***************************
	GeoNo   = fPm->GetIntegerParameter(GetFullParmName("GeometryNumber"));
	CenterX = fPm->GetDoubleParameter(GetFullParmName("TransX"),"Length");
	CenterY = fPm->GetDoubleParameter(GetFullParmName("TransY"),"Length");
	CenterZ = fPm->GetDoubleParameter(GetFullParmName("TransZ"),"Length");
	SVradius = 	fPm->GetDoubleParameter(GetFullParmName("SensitiveVolumeRadius"),"Length");
	TEradius =  fPm->GetDoubleParameter(GetFullParmName("TissueEquivalentRadius"),"Length"); 
	NumberOfHistoriesInRun =  fPm->GetIntegerParameter(GetFullParmName("NumberOfHistoriesInRun"));

	if(GeoNo == 1)
	{
		SVheight = SVradius;
		if(fPm->ParameterExists(GetFullParmName("SensitiveVolumeHalfLength")))
			SVheight=fPm->GetDoubleParameter(GetFullParmName("SensitiveVolumeHalfLength"),"Length");
	}


	// *************************** optional parameters ***************************
	ylimitL = 0 ;  // in unit of keV/um
	if ( fPm->ParameterExists(GetFullParmName("LinealEnergyLowerlimit")) ) {
		ylimitL  =  fPm->GetUnitlessParameter(GetFullParmName("LinealEnergyLowerlimit"));
	}

	ylimitU = 1E4 ;  // in unit of keV/um
	if ( fPm->ParameterExists(GetFullParmName("LinealEnergyUpperlimit")) ) {
		ylimitU  =  fPm->GetUnitlessParameter(GetFullParmName("LinealEnergyUpperlimit"));
	}

	IncludeYF = true;
	if ( fPm->ParameterExists(GetFullParmName("IncludeFrequencyMeanLinealEnergy")) ) {
		IncludeYF = fPm->GetBooleanParameter(GetFullParmName("IncludeFrequencyMeanLinealEnergy"));
	}

	IncludeYD = true;
	if ( fPm->ParameterExists(GetFullParmName("IncludeDoseMeanLinealEnergy")) ) {
		IncludeYD = fPm->GetBooleanParameter(GetFullParmName("IncludeDoseMeanLinealEnergy"));
	}

	fGetSecondariesContribution = true;
	if ( fPm->ParameterExists(GetFullParmName("GetContributionOfSecondaries")) ) {
		fGetSecondariesContribution = fPm->GetBooleanParameter(GetFullParmName("GetContributionOfSecondaries"));
	}

	fGetStatisticInfo = false;
	if ( fPm->ParameterExists(GetFullParmName("GetStatisticInfo")) ) {
		fGetStatisticInfo = fPm->GetBooleanParameter(GetFullParmName("GetStatisticInfo"));
	}

	fSpectrumUpdateTimes = 1000;
	if ( fPm->ParameterExists(GetFullParmName("SpectrumUpdateTimes")) ) {
		fSpectrumUpdateTimes = fPm->GetIntegerParameter(GetFullParmName("SpectrumUpdateTimes"));
	}


	P_pileup=0;
	if ( fPm->ParameterExists(GetFullParmName("PileupProbability")) )
		P_pileup=fPm->GetUnitlessParameter(GetFullParmName("PileupProbability"));
	G4cout<<"pile-up probability = "<< P_pileup<<G4endl;



	if(GeoNo == 0 || GeoNo == 1)
		fMeanChordLength = 4.*TEradius/3;
	else if (GeoNo ==2)
	{
		MeanPathLength = fPm->GetDoubleParameter(GetFullParmName("MeanPathLength"),"Length");
		fMeanChordLength = MeanPathLength;
	}   
	else
	{
		G4cout<<"Error: The detector is not supported yet."<< G4endl;  
		exit(0);
	}


	//**************************************************  Survival, RBE, QUality Factor *******************************
	fGetRBEWithBioWeightFunction = false;
	if ( fPm->ParameterExists(GetFullParmName("GetRBEWithBiologicalWeightFunction")) ) {
		fGetRBEWithBioWeightFunction = fPm->GetBooleanParameter(GetFullParmName("GetRBEWithBiologicalWeightFunction"));
	}
	fGetRBEWithMKModel = false;
	if ( fPm->ParameterExists(GetFullParmName("GetRBEWithMKModel")) ) {
		fGetRBEWithMKModel = fPm->GetBooleanParameter(GetFullParmName("GetRBEWithMKModel"));
	} 
	fGetRBEWithGSM2 = false;
	if ( fPm->ParameterExists(GetFullParmName("GetRBEWithGSM2")) ) {
		fGetRBEWithGSM2 = fPm->GetBooleanParameter(GetFullParmName("GetRBEWithGSM2"));
	}
	fGetQualityFactor = false;
	if ( fPm->ParameterExists(GetFullParmName("GetQualityFactor")) ) {
		fGetQualityFactor = fPm->GetBooleanParameter(GetFullParmName("GetQualityFactor"));
	}
	fSaveSpecificEnergySpectra = false;
	if ( fPm->ParameterExists(GetFullParmName("SaveSpecificEnergySpecra")) ) {
		fSaveSpecificEnergySpectra = fPm->GetBooleanParameter(GetFullParmName("SaveSpecificEnergySpecra"));
	}
	fSetMultiEventStatistic = 1e4;
	if ( fPm->ParameterExists(GetFullParmName("SetMultieventStatistic")) ) {
		fSetMultiEventStatistic = fPm->GetUnitlessParameter(GetFullParmName("SetMultieventStatistic"));
	}



	G4cout << "********************** Simulation info. ********************************"<< G4endl; 
	G4cout << "Detector type : ";
	if (GeoNo ==0 ) G4cout << "spherical TEPC "<< G4endl; 
	if (GeoNo ==1 ) G4cout << "cylindrical TEPC "<< G4endl; 
	if (GeoNo ==2 ) G4cout << "silicon microdosimeter "<< G4endl; 
	G4cout << "Radius of sensitive volume =" << SVradius/um << " um" << G4endl;   
	if(GeoNo ==1 )
		G4cout << "Height of sensitive volume =" << SVheight/um << " um" << G4endl; 
	G4cout << "Radius of tissue equivalent volume = "<< TEradius/um<< " um" << G4endl;
	G4cout << "Mean chord lengh of sensitive volume = "<< fMeanChordLength/um << " um" << G4endl;
	G4cout << "Sensitive volume center = ("<< CenterX/mm <<", "<<CenterY/mm<<", "<<CenterZ/mm<<")*mm" << G4endl;
	G4cout << "************************************************************************"<< G4endl; 

	for(int i=0; i<10;i++){
		EnergyDepInE_Particle_pu[i] = 0;
	}
	NumofPileupEvent =0;
	SVtotalEdep = 0;
}


TsYScorer::~TsYScorer() {;}

G4bool TsYScorer::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

	// energy deposit
	G4double edep = aStep->GetTotalEnergyDeposit();
	if (edep==0.) return false;

	TsTrackerHit* newHit = new TsTrackerHit();
	newHit->SetTrackID (aStep->GetTrack()->GetTrackID());
	newHit->SetEdep(edep);
	newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());

	G4String ParticleName =aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
	newHit->SetParticleName(ParticleName);

	G4int flagParticle = -1;
	if(aStep->GetTrack()->GetTrackID() == 1 && aStep->GetTrack()->GetParentID() == 0)
		flagParticle = 1; //DEBUG
	else
	{
		if (ParticleName == "e-") flagParticle = 0;
		else if (ParticleName == "proton")   flagParticle = 2; //DEBUG
		else if (ParticleName == "deuteron") flagParticle = 2;
		else if (ParticleName == "triton")   flagParticle = 2;
		else if (ParticleName == "He3")      flagParticle = 3;
		else if (ParticleName == "alpha")    flagParticle = 3;
		else if (ParticleName == "He6")      flagParticle = 3;
		else if (ParticleName == "Li7")      flagParticle = 4;
		else if (ParticleName == "Be7")      flagParticle = 5;
		else if (ParticleName == "Be9")      flagParticle = 5;
		else if (ParticleName == "B10")      flagParticle = 6;
		else if (ParticleName == "B11")      flagParticle = 6;
		else if (ParticleName == "C11")      flagParticle = 7;
		else if (ParticleName == "C12")      flagParticle = 7;
		else                                 flagParticle = 8;
		newHit->SetParticleFlag(flagParticle);

	}

	
	if (aStep->GetTrack()->GetTrackID() == 1 && aStep->GetTrack()->GetParentID() == 0)
		newHit->SetIncidentEnergy(aStep->GetTrack()->GetVertexKineticEnergy());

	fHitsCollection->insert(newHit);

	return true;
}

void TsYScorer::AccumulateEvent()
{
	G4int nofHits = fHitsCollection->entries();
	G4double Einc=0;
	G4double ylowerlimit = ylimitL;
	G4double yupperlimit = ylimitU;
	const G4int    ParticleType    = 10; //DEBUG
	G4double * EdepInEvent_Particle = new G4double[ParticleType] {0};
	G4double * yVectorInEvent = new G4double[ParticleType] {0};

	// ************************************
	if (nofHits != 0)
	{
		G4double SphereCenterX =  CenterX;
		G4double SphereCenterY =  CenterY;
		G4double SphereCenterZ =  CenterZ;
		//******************************************************************
		//                    Spherical TEPC  
		//******************************************************************
		if(GeoNo==0)    
		{	
			G4ThreeVector CenterPos(SphereCenterX,SphereCenterY,SphereCenterZ);
			G4double epsilon = 0;

			for ( G4int i=0; i<nofHits; i++ )
			{
				if ((*fHitsCollection)[i]->GetIncidentEnergy()>0)
					Einc = (*fHitsCollection)[i]->GetIncidentEnergy();

				G4ThreeVector localPos = (*fHitsCollection)[i]->GetPos();
				if (
						(localPos.x()-CenterPos.x()) * (localPos.x()-CenterPos.x()) +
						(localPos.y()-CenterPos.y()) * (localPos.y()-CenterPos.y()) +
						(localPos.z()-CenterPos.z()) * (localPos.z()-CenterPos.z())
						< 1.00001*SVradius*SVradius
				   )
				{
					G4double edep = (*fHitsCollection)[i]->GetEdep();
					G4int flag = (*fHitsCollection)[i]->GetParticleFlag();
					EdepInEvent_Particle[flag] += edep;  // save different particle edep
					EdepInEvent_Particle[9] += edep;     // save total particle edep

				}
			}

			if (EdepInEvent_Particle[9]!=0) //ene tot !=0
			{
				G4double pp=G4UniformRand();
				if (pp<=P_pileup) 
				{
					// save energy deposition info into pile up 
					NumofPileupEvent ++;
					for(int i=0; i<10;i++)
						EnergyDepInE_Particle_pu[i] += EdepInEvent_Particle[i];	
				}
				else 
				{
					// read energy deposition info from pile up  && clear pile up info
					for(int i=0; i<10;i++){
						EdepInEvent_Particle[i] += EnergyDepInE_Particle_pu[i];
						EnergyDepInE_Particle_pu[i] = 0;
					}
					PushBackData(EdepInEvent_Particle);

				}
			}
			//PushBackData(EdepInEvent_Particle);

		}       //  End of spherical TEPC  


		//******************************************************************
		//                    Cylindrical TEPC  
		//   The axis of the cylinder should along the y axis of the world
		//******************************************************************
		else if(GeoNo==1)       
		{    
			G4ThreeVector CenterPos(SphereCenterX,SphereCenterY,SphereCenterZ);
			G4double epsilon = 0;

			for ( G4int i=0; i<nofHits; i++ )
			{
				if ((*fHitsCollection)[i]->GetIncidentEnergy()>0)
					Einc = (*fHitsCollection)[i]->GetIncidentEnergy();

				G4ThreeVector localPos = (*fHitsCollection)[i]->GetPos();                                  
				if (
						((localPos.x()-CenterPos.x()) * (localPos.x()-CenterPos.x()) +
						 (localPos.z()-CenterPos.z()) * (localPos.z()-CenterPos.z()) 
						 < 1.00001*SVradius*SVradius)
						&& ((localPos.y()-CenterPos.y()) * (localPos.y()-CenterPos.y()) < 1.00001*SVheight*SVheight)
				   )
				{
					G4double edep = (*fHitsCollection)[i]->GetEdep();
					G4int flag = (*fHitsCollection)[i]->GetParticleFlag();
					EdepInEvent_Particle[flag] += edep;  // save different particle edep
					EdepInEvent_Particle[9] += edep;     // save total particle edep
				}  
			}

			if ( EdepInEvent_Particle[9]!=0)
				PushBackData(EdepInEvent_Particle);

		}       // End of cylindrical mini-TEPC  

		//******************************************************************
		//                Silicon microdosimeter array  
		//******************************************************************
		else if(GeoNo==2)  
		{

			G4double x1 = 15.0*um;
			G4double y1 = 5.0*um;
			G4double CenterXminSV = -1450.;
			G4double CenterYminSV = -1775.;
			G4double CenterXminBridge = -1457.5;
			G4double CenterYminBridge = -1750.;
			G4double EdepinChannel[59]={0};		

			for (G4int ii=0;ii<59;ii++)
			{            
				for (G4int jj=0;jj<72;jj++)
				{
					G4double epsilon = 0;

					SphereCenterX = (ii*50.0+CenterXminSV)*um + CenterX;
					SphereCenterY = (jj*50.0+CenterYminSV)*um + CenterY;
					SphereCenterZ = -148.935*um + CenterZ;                    
					G4ThreeVector CenterPos(SphereCenterX,SphereCenterY,SphereCenterZ);

					for ( G4int i=0; i<nofHits; i++ )
					{
						if ((*fHitsCollection)[i]->GetIncidentEnergy()>0)
							Einc = (*fHitsCollection)[i]->GetIncidentEnergy();

						G4ThreeVector localPos = (*fHitsCollection)[i]->GetPos();
						if (						
								((localPos.x()-CenterPos.x()) * (localPos.x()-CenterPos.x()) <1.0001*x1*x1) &&
								((localPos.y()-CenterPos.y()) * (localPos.y()-CenterPos.y()) <1.0001*x1*x1) &&
								((localPos.z()-CenterPos.z()) * (localPos.z()-CenterPos.z()) <1.0001*y1*y1)
						   )
						{
							G4double TEFactor = 0.58;  // tissue equivalent conversion factor
							G4double edep = (*fHitsCollection)[i]->GetEdep();
							edep = edep*TEFactor;
							G4int flag = (*fHitsCollection)[i]->GetParticleFlag();
							EdepInEvent_Particle[flag] += edep;  // save different particle edep
							EdepInEvent_Particle[9] += edep;     // save total particle edep   

							EdepinChannel[ii] = EdepinChannel[ii] + edep ;
						}
					}
				}

				if (EdepinChannel[ii]!=0)
					PushBackData(EdepInEvent_Particle);
			}

		}  // End of silicon microdosimeter array   

	} // MUST BE HITS
	// initialization of fHitsCollection for every single event

	delete fHitsCollection;
	fHitsCollection = new TsTrackerHitsCollection();

}

void TsYScorer::PushBackData(G4double* Edep)
{
	const G4int    Type_of_Particle   = 10; //DEBUG
	G4double * yVectorInEvent = new G4double[Type_of_Particle] {0};

	fy_z0 = (Edep[0]/keV)/(fMeanChordLength/um);
	fy_z_Prim = (Edep[1]/keV)/(fMeanChordLength/um); //DEBUG
	fy_z1_Seco = (Edep[2]/keV)/(fMeanChordLength/um); //DEBUG
	fy_z2 = (Edep[3]/keV)/(fMeanChordLength/um);
	fy_z3 = (Edep[4]/keV)/(fMeanChordLength/um);
	fy_z4 = (Edep[5]/keV)/(fMeanChordLength/um);
	fy_z5 = (Edep[6]/keV)/(fMeanChordLength/um);
	fy_z6 = (Edep[7]/keV)/(fMeanChordLength/um);
	fy_z_ = (Edep[8]/keV)/(fMeanChordLength/um);
	fy    = (Edep[9]/keV)/(fMeanChordLength/um);

	for(G4int i=0; i<Type_of_Particle; i++)
		yVectorInEvent[i] = (Edep[i]/keV)/(fMeanChordLength/um);

	if (fy > ylimitL && fy < ylimitU)
	{
		fNtuple->Fill();
		yVector.push_back(fy);
		yVector_Particle.push_back(yVectorInEvent);  
	}

	SVtotalEdep += Edep[9]; //DEBUG

}

void TsYScorer::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

	TsYScorer* workerMTScorer = dynamic_cast<TsYScorer*>(workerScorer);
	std::vector<G4double> yVector_worker = workerMTScorer->yVector;
	std::vector<G4double * > yVector_Particle_worker= workerMTScorer->yVector_Particle;

	SVtotalEdep += workerMTScorer->SVtotalEdep;
	for (G4int i=0; i<yVector_Particle_worker.size(); i++){
		G4double * aVector = yVector_Particle_worker[i];
		yVector_Particle.push_back(aVector);
		yVector.push_back(yVector_worker[i]);
	}

	workerMTScorer->yVector.clear();
	workerMTScorer->yVector_Particle.clear();
}

void TsYScorer::UserHookForEndOfRun()
{
	G4cout<<"Total energy deposition in sensitive volume is: "<<SVtotalEdep/MeV<<" MeV"<<G4endl;
	G4cout<<"Total number of pile-up events is: "<<NumofPileupEvent<<G4endl;

	// calculate statistic information
	InitializeStatistic();
	//GetSpectrum();
	if (fGetStatisticInfo){
		G4cout <<"Start statistic error calcuation ..."<<G4endl;
		clock_t start,end;
		start = clock();
		if (fSpectrumUpdateTimes>yVector.size()){
			G4int size = yVector.size();
			fSpectrumUpdateTimes = min(size,1000);
			G4cout << "Statistic Update Frequency is larger than scored events("<<yVector.size() <<")!"<<G4endl;
			G4cout << "Statistic Update Frequency is reset as "<<fSpectrumUpdateTimes<<G4endl;
		}

		for (G4int i=1; i<=fSpectrumUpdateTimes; i++){
			G4int dataSize = ceil(yVector.size()/fSpectrumUpdateTimes);
			std::vector<G4double> dataVector;
			for(int j=0; j<dataSize*i && j<yVector.size(); j++)
				dataVector.push_back(yVector[j]);

			Calculatefy(dataVector);
			dataVector.clear();

			if(fSpectrumUpdateTimes>=100 )
			{
				G4int tenPercent = ceil(fSpectrumUpdateTimes/10);
				if (i%tenPercent==0){
					end = clock();
					float duration = (float) (end - start)/ CLOCKS_PER_SEC;
					G4cout << (i/tenPercent)*10<< " % of statistic information update finished, total time used :"<<duration<<" sec; ("<<duration/60<<" min)"<<G4endl;
				}
			}
		}
		G4cout<<"\n"<<G4endl;
	}

	GetSpectrum();
	GetErrorPropagation();
	GetRBE();
	OutputResult();

}

void TsYScorer::Calculatefy( std::vector<G4double> dataVector)
{
	G4double * histfy = new G4double [yBinNum] {0};
	G4int nnum=0;
	for (std::vector<G4double>::const_iterator i = dataVector.begin(); i != dataVector.end(); ++i){
		for (G4int n=0;n<yBinNum;n++){
			if(*i<=BinLimit[n+1]){
				histfy[n] = histfy[n]+1;
				break;
			}
		}
		nnum=nnum+1;
	}

	for (G4int i=0;i<yBinNum;i++){
		histfy[i] = histfy[i]/(BinWidth[i]*nnum);    // normalization. divide by bin width * number of entries
		GetStatisticInfo(i, histfy[i]);              // Calculate statistic error
	}
}


void TsYScorer::InitializeMicrodosimetricSpectrum()
{
	hfy      = new G4double [yBinNum] {0};
	hdy      = new G4double [yBinNum] {0};
	hyfy     = new G4double [yBinNum] {0};
	hydy     = new G4double [yBinNum] {0};
	hfy_particle = new G4double *[yBinNum];
	for (G4int i=0; i<yBinNum; i++)
		hfy_particle[i] =new G4double [10];

	yParticleContribution.resize(yBinNum, std::vector<double>(10,0.)); //DEBUG

	BinLimit = new G4double [yBinNum] {0};
	BinWidth = new G4double [yBinNum] {0};

	BinLimit[0]=0.1;
	for (G4int i=0;i<yBinNum;i++){
		G4double aa = (double)((i+1)/yBinMagnitudeInterval);
		BinLimit[i+1] = pow(10,(aa -1.0));  
		BinWidth[i] = BinLimit[i+1]-BinLimit[i];   
	}

}

void TsYScorer::InitializeStatistic()
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


void TsYScorer::GetSpectrum()
{
	G4cout<<"yVector_Particle size="<<yVector_Particle.size()<< G4endl;
	G4cout<<"yVector size="<<yVector.size()<< G4endl;

	G4int nnum=0;
	G4int index=0;
	for (std::vector<G4double>::const_iterator i = yVector.begin(); i != yVector.end(); ++i){
		for (G4int n=0;n<yBinNum;n++){
			if(*i<=BinLimit[n+1]){
				hfy[n] = hfy[n]+1;
				for(G4int loop = 0; loop<10; loop++)
					hfy_particle[n][loop] += yVector_Particle[index][loop];
				break;
			}
		}
		nnum=nnum+1;
		index++;
	}


	for (G4int i=0;i<yBinNum;i++){
		hfy[i] = hfy[i]/(BinWidth[i]*nnum);                    // normalization. divide by bin width * number of entries
		hyfy[i] = (BinLimit[i]+BinLimit[i+1])/2*hfy[i];        // calculate y*f(y) = BinCenter * BinContent

		for(int loop=0; loop<10; loop++) //DEBUG
		{
			if(hfy_particle[i][9] != 0) //DEBUG
				yParticleContribution[i][loop] = hfy_particle[i][loop]/hfy_particle[i][9]; //DEBUG
		}
	}

	//******************************************************************
	//                 Validate f(y) & calculate yF
	//******************************************************************
	G4double Probability_fy =0;
	for(G4int i=0;i<yBinNum;i++)
		Probability_fy += hfy[i]*BinWidth[i];
	//G4cout << "sum of f(y)*delta_y ="<< Probability_fy<<G4endl;    

	//calculate yF
	yF=0;
	for (G4int i=0;i<yBinNum;i++){
		yF = yF + hyfy[i]*BinWidth[i];          // multiply by bin width
	}

	for (G4int i=0;i<yBinNum;i++){
		hdy[i] = hyfy[i]/yF;                                    //calculate d(y) = y*f(y)/yF (cf. Burigo et al., NIMB 320 (2014))
		hydy[i] = (BinLimit[i]+BinLimit[i+1])/2*hdy[i];         // calculate y*d(y) = BinCenter * d(y)
	}


	//******************************************************************
	//               Validate d(y) & calculate yF
	//******************************************************************
	G4double Probability_dy =0;
	for(G4int i=0;i<yBinNum;i++)
		Probability_dy += hdy[i]*BinWidth[i];
	//G4cout << "sum of d(y)*delta_y ="<< Probability_dy<<G4endl;

	//calculate yD
	yD=0;
	for (G4int i=0;i<yBinNum;i++){
		yD = yD + hydy[i]*BinWidth[i];          // multiply by bin width
	}
}

void TsYScorer::GetRBE()
{
	//******************************************************************
	//                           Calculate S, RBE, Qf
	//******************************************************************
	if(fGetRBEWithBioWeightFunction == true || fGetRBEWithMKModel == true || fGetRBEWithGSM2 == true || fGetQualityFactor == true) //DEBUG
	{

		//Create temp matrix defined as vector x vector instead of vector x double*
		std::vector<std::vector<double>> tmp_yVector_Particle(yVector_Particle.size(), std::vector<double>(10,0.)); //DEBUG 
		for(int i=0; i<yVector_Particle.size(); i++) 	//DEBUG
		{
			for(int loop=0; loop<10; loop++) 	//DEBUG
			{
				tmp_yVector_Particle[i][loop] = yVector_Particle[i][loop]; //DEBUG
			}
		}

		TsGetSurvivalRBEQualityFactor *aSurvRBEQf = new TsGetSurvivalRBEQualityFactor(yParticleContribution, tmp_yVector_Particle, BinLimit, BinWidth, hfy, hdy, yF, yD, yF_var, yD_var, fVariance, dy_var, yBinNum, fGetStatisticInfo, fSpectrumUpdateTimes, fGetSecondariesContribution); //DEBUG

		aSurvRBEQf -> SetMCMultieventIterations(fSetMultiEventStatistic);

		if(fGetQualityFactor == true)
		{	
			if ( fPm->ParameterExists(GetFullParmName("QualityFactorCalculation")) ) 
			{
				const G4int Ncalc = fPm -> GetVectorLength(GetFullParmName("QualityFactorCalculation"));
				G4String* MKMCalculation_array = fPm -> GetStringVector(GetFullParmName("QualityFactorCalculation"));

				for(G4int i=0; i<Ncalc; i++)
				{
					if(strcmp(MKMCalculation_array[i], "ICRU40") == 0)
						aSurvRBEQf->GetQualityFactorWithICRU40();

					if(strcmp(MKMCalculation_array[i], "KellererHahn") == 0)
						aSurvRBEQf -> GetQualityFactorWithKellereHahn();

				}

			}
			else //Default Saturation Correction is used
				aSurvRBEQf->GetQualityFactorWithICRU40();			//DEBUG	
		}

		if(fGetRBEWithBioWeightFunction == true)
		{
			if(fPm->ParameterExists(GetFullParmName("BiologicalWeightFunctionDataFile")))           
			{
				G4String DataFilename =fPm->GetStringParameter("BiologicalWeightFunctionDataFile");
				aSurvRBEQf->SetBioWeightFunctionDataFile(DataFilename);   
			}
			aSurvRBEQf->GetRBEWithBioWeightFunction();
		}


		//Set Macroscopic Doses in GSM2 and MKM for survival calculation
		std::vector<double> doses;
		if( fPm->ParameterExists(GetFullParmName("SurvivalDoses")) )
		{
			G4double* parameter = fPm -> GetUnitlessVector(GetFullParmName("SurvivalDoses"));
			aSurvRBEQf -> SetDosesMacro(parameter);
			for(double d=parameter[0]; d<=parameter[1]; d+=parameter[2])
				doses.push_back(d);
		}


		//******************************************************************
		//                           MKM model parameters
		//******************************************************************

		if(fGetRBEWithMKModel==true)
		{
			if ( fPm->ParameterExists(GetFullParmName("MKModel_alpha0")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("MKModel_alpha0"));
				aSurvRBEQf->SetMKModel_alpha0(parameter);
			}
			if ( fPm->ParameterExists(GetFullParmName("MKModel_beta")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("MKModel_beta"));
				aSurvRBEQf->SetMKModel_beta(parameter);
			}
			if ( fPm->ParameterExists(GetFullParmName("MKModel_alphaX")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("MKModel_alphaX"));
				aSurvRBEQf->SetMKModel_alphaX(parameter);
			}
			if ( fPm->ParameterExists(GetFullParmName("MKModel_betaX")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("MKModel_betaX"));
				aSurvRBEQf->SetMKModel_betaX(parameter);
			}	
			if ( fPm->ParameterExists(GetFullParmName("MKModel_rho")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("MKModel_rho"));
				aSurvRBEQf->SetMKModel_rho(parameter);
			}
			if ( fPm->ParameterExists(GetFullParmName("MKModel_DomainRadius")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("MKModel_DomainRadius"));
				aSurvRBEQf->SetMKModel_rd(parameter);
			}
			G4double MKMNucleusRadius = 0.;
			if ( fPm->ParameterExists(GetFullParmName("MKModel_NucleusRadius")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("MKModel_NucleusRadius"));
				aSurvRBEQf->SetMKModel_Rn(parameter);
				MKMNucleusRadius = parameter;
			}
			if ( fPm->ParameterExists(GetFullParmName("MKModel_y0")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("MKModel_y0"));
				aSurvRBEQf->SetMKModel_y0(parameter);
			}

			if(fSaveSpecificEnergySpectra==true) {
				OutputResultSpecificEnergy(MKMNucleusRadius, doses);
			}

			if ( fPm->ParameterExists(GetFullParmName("MKMCalculation")) ) 
			{
				const G4int Ncalc = fPm -> GetVectorLength(GetFullParmName("MKMCalculation"));
				G4String* MKMCalculation_array = fPm -> GetStringVector(GetFullParmName("MKMCalculation"));

				for(G4int i=0; i<Ncalc; i++)
				{
					if(strcmp(MKMCalculation_array[i], "SaturationCorrection") == 0)
						aSurvRBEQf->GetSurvWithMKModel_SaturationCorr();

					if(strcmp(MKMCalculation_array[i], "nonPoisson") == 0)
						aSurvRBEQf -> GetSurvWithMKModel_nonPoissonCorr();

					if(strcmp(MKMCalculation_array[i], "SMKM") == 0)
						aSurvRBEQf -> GetSurvWithSMKModel();

					if(strcmp(MKMCalculation_array[i], "DSMKM") == 0)
						aSurvRBEQf -> GetSurvWithDSMKModel();
				}

			}
			else //Default Saturation Correction is used
				aSurvRBEQf->GetSurvWithMKModel_SaturationCorr();
		}

		//******************************************************************
		//                           GSM2 model parameters
		//******************************************************************

		if(fGetRBEWithGSM2==true)
		{	
			if ( fPm->ParameterExists(GetFullParmName("GSM2_DomainRadius")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("GSM2_DomainRadius"));
				aSurvRBEQf->SetGSM2_rd(parameter);
			}
			G4double GSM2NucleusRadius = 0.;
			if ( fPm->ParameterExists(GetFullParmName("GSM2_NucleusRadius")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("GSM2_NucleusRadius"));
				aSurvRBEQf->SetGSM2_Rn(parameter);
				GSM2NucleusRadius = parameter;
			}
			if ( fPm->ParameterExists(GetFullParmName("GSM2_a")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("GSM2_a"));
				aSurvRBEQf->SetGSM2_a(parameter);
			}
			if ( fPm->ParameterExists(GetFullParmName("GSM2_b")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("GSM2_b"));
				aSurvRBEQf->SetGSM2_b(parameter);
			}
			if ( fPm->ParameterExists(GetFullParmName("GSM2_r")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("GSM2_r"));
				aSurvRBEQf->SetGSM2_r(parameter);
			}
			if ( fPm->ParameterExists(GetFullParmName("GSM2_alphaX")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("GSM2_alphaX"));
				aSurvRBEQf->SetGSM2_alphaX(parameter);
			}
			if ( fPm->ParameterExists(GetFullParmName("GSM2_betaX")) ) {
				G4double parameter = fPm->GetUnitlessParameter(GetFullParmName("GSM2_betaX"));
				aSurvRBEQf->SetGSM2_betaX(parameter);
			}
			if(fSaveSpecificEnergySpectra==true)
				OutputResultSpecificEnergy(GSM2NucleusRadius, doses);


			aSurvRBEQf -> GetSurvWithGSM2();
		}

	}
}

void TsYScorer::OutputResultSpecificEnergy(G4double NucleusRadius, std::vector<double> MacroscopicDoses)
{
	//Create temp matrix defined as vector x vector instead of vector x double*
	std::vector<std::vector<double>> tmp_yVector_Particle(yVector_Particle.size(), std::vector<double>(10,0.)); //DEBUG 
	for(int i=0; i<yVector_Particle.size(); i++) 	//DEBUG
	{
		for(int loop=0; loop<10; loop++) 	//DEBUG
		{
			tmp_yVector_Particle[i][loop] = yVector_Particle[i][loop]; //DEBUG
		}
	}



	TsSpecificEnergy* aSpecificEnergy = new TsSpecificEnergy(tmp_yVector_Particle, NucleusRadius, fGetStatisticInfo, fSpectrumUpdateTimes);
	std::vector<double> zBinCenter = aSpecificEnergy -> GetBinCenter();
	std::vector<double> zBinWidth = aSpecificEnergy -> GetBinWidth();
	double zF = aSpecificEnergy -> GetzF();
	double zF_std = aSpecificEnergy -> GetzF_std();
	std::vector<std::vector<double>> zParticleContribution = aSpecificEnergy -> GetParticleContribution();
	std::vector<double> hfz = aSpecificEnergy -> GetHfz();
	std::vector<double> hfz_std = aSpecificEnergy -> GetHfz_std();
	std::vector<double> hzfz = aSpecificEnergy -> GetHzfz();
	std::vector<double> hzfz_std = aSpecificEnergy -> GetHzfz_std();
	std::vector<std::vector<double>> hfzMultievent, hfzMultievent_std;
	for(double D:MacroscopicDoses)
	{
		hfzMultievent.push_back(aSpecificEnergy->GetHfzMultiEvent(D, fSetMultiEventStatistic));
		hfzMultievent_std.push_back(aSpecificEnergy->GetHfzMultievent_var());
	}
	//******************************************************************
	//                          OutputResult
	//******************************************************************
	FILE* zSpecfile;
	zSpecfile=fopen("zSpecfile.txt","w");

	G4cout <<"\n***************************** SPECIFIC ENERGY **********************\n"; 
	G4cout << "zF = " << zF <<" ( std: "<< zF_std <<")"<< " keV/um" << G4endl;
	G4cout <<"**********************************************************************\n";

	fprintf(zSpecfile,"zF = %.4e (std: %.4e) keV/um\n",zF, zF_std);
	fprintf(zSpecfile,"zn (keV/um)     fn,1(zn)(std)                  znfn,1(zn)(std)          \n");
	for (G4int i=0;i<zBinCenter.size();i++)
	{       
		fprintf(zSpecfile,"%.4e ",zBinCenter[i]);
		fprintf(zSpecfile,"    ");
		fprintf(zSpecfile,"%.4e  %.4e ",hfz[i], hfz_std[i]);
		fprintf(zSpecfile,"    ");
		fprintf(zSpecfile,"%.4e  %.4e ",hzfz[i], hzfz_std[i]);
		fprintf(zSpecfile,"    ");
		fprintf(zSpecfile,"\n");

	}
	fclose(zSpecfile);

	FILE* zSpecfileMultievent;
	zSpecfileMultievent=fopen("zSpecfileMultievent.txt","w");
	//HEADER FILE
	fprintf(zSpecfileMultievent,"zn [keV/um]");
	for(double D:MacroscopicDoses) 
		fprintf(zSpecfileMultievent,"      fn,n(zn,%.1f)(std)", D);
	fprintf(zSpecfileMultievent, "\n");
	//BODY FILE
	for (G4int i=0;i<zBinCenter.size();i++)
	{       
		fprintf(zSpecfileMultievent,"%.4e ",zBinCenter[i]);
		fprintf(zSpecfileMultievent,"    ");
		for(int d=0; d<MacroscopicDoses.size(); d++)
		{
			fprintf(zSpecfileMultievent,"%.4e  %.4e ", hfzMultievent[d][i], hfzMultievent_std[d][i]);
			fprintf(zSpecfileMultievent,"    ");
		}
		fprintf(zSpecfileMultievent,"\n");

	}
	fclose(zSpecfileMultievent);

	//******************************************************************
	//                          Output particle contribution
	//******************************************************************
	if(fGetSecondariesContribution)
	{
		FILE* zSpec_particle;
		zSpec_particle=fopen("zSpecfile_Particle.txt","w"); 
		fprintf(zSpec_particle,"z(Gy)    e-           Zprim            Hsec            He           Li           Be           B            C            Other        Total[zn,1f(zn,1)]\n");

		for (G4int i=0;i<zBinCenter.size();i++)
		{       
			fprintf(zSpec_particle,"%.4e ",zBinCenter[i]);
			fprintf(zSpec_particle,"  ");

			for(G4int loop=0; loop<10; loop++){
				fprintf(zSpec_particle,"%.4e ",hfz[i]*zParticleContribution[i][loop]);
				fprintf(zSpec_particle,"  ");
			}

			fprintf(zSpec_particle,"\n");

		}
		fclose(zSpec_particle);
	}

}

void TsYScorer::OutputResult()
{
	//******************************************************************
	//                          OutputResult
	//******************************************************************
	FILE* ydyfile;
	ydyfile=fopen("ydy.txt","w");  
	for (G4int i=0;i<yBinNum;i++){
		fprintf(ydyfile,"%.4e keV/um",(BinLimit[i]+BinLimit[i+1])/2);
		fprintf(ydyfile,"   ");
		fprintf(ydyfile,"%.4e ",hydy[i]);
		fprintf(ydyfile,"\n");
	}
	fclose(ydyfile);

	FILE* ySpecfile;
	ySpecfile=fopen("ySpecfile.txt","w"); 
	if(IncludeYF){
		G4cout << "yF = " << yF <<" ( std: "<<yF_std<<")"<< " keV/um" << G4endl;
		fprintf(ySpecfile,"yF = %.4e (std: %.4e) keV/um\n",yF, yF_std);
	}
	if(IncludeYD){
		G4cout << "yD = " << yD <<" ( std: "<<yD_std<<")"<< " keV/um" << G4endl;  
		fprintf(ySpecfile,"yD = %.4e (std: %.4e) keV/um\n",yD, yD_std);
	}
	fprintf(ySpecfile,"y (keV/um)     f(y)(std)                  yf(y)(std)                 d(y)(std)                  yd(y)(std)  \n");
	for (G4int i=0;i<yBinNum;i++){       
		fprintf(ySpecfile,"%.4e ",(BinLimit[i]+BinLimit[i+1])/2);
		fprintf(ySpecfile,"    ");
		fprintf(ySpecfile,"%.4e  %.4e ",hfy[i], fStandardDeviation[i]);
		fprintf(ySpecfile,"    ");
		fprintf(ySpecfile,"%.4e  %.4e ",hyfy[i], yfy_std[i]);
		fprintf(ySpecfile,"    ");
		fprintf(ySpecfile,"%.4e  %.4e ",hdy[i], dy_std[i] );
		fprintf(ySpecfile,"    ");
		fprintf(ySpecfile,"%.4e  %.4e ",hydy[i], ydy_std[i]);
		fprintf(ySpecfile,"\n");
	}
	fclose(ySpecfile);

	//******************************************************************
	//                          Output particle contribution
	//******************************************************************
	if(fGetSecondariesContribution)
	{
		FILE* ySpec_particle;
		ySpec_particle=fopen("ySpecfile_Particle.txt","w"); 
		fprintf(ySpec_particle,"y(keV/um)    e-           Zprim            Hsec            He           Li           Be           B            C            Other        Total[yd(y)]\n");

		for (G4int i=0;i<yBinNum;i++)
		{       
			fprintf(ySpec_particle,"%.4e ",(BinLimit[i]+BinLimit[i+1])/2);
			fprintf(ySpec_particle,"  ");

			for(G4int loop=0; loop<10; loop++){
				G4double ydyParticleContribution = (hfy_particle[i][loop]/hfy_particle[i][9])*hydy[i];
				if(hfy_particle[i][9]==0)                 ydyParticleContribution=0;
				fprintf(ySpec_particle,"%.4e ",ydyParticleContribution);
				fprintf(ySpec_particle,"  ");
			}

			fprintf(ySpec_particle,"\n");

		}
		fclose(ySpec_particle);
	}

}


void TsYScorer::GetStatisticInfo(G4int Binindex , G4double variable /*  hfy[i]*/)
{
	// Bin index
	G4int index;

	// Value from one specific bin
	G4double x;
	G4double mean;
	G4double delta;
	G4double mom2;
	G4double recorededHistories;

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
	// G4cout << "index="<<index<<" x="<<x<<" recorededHistories="<<recorededHistories<<" fCountMap[index]"<<fCountMap[index];
	// G4cout<<" FirstMoment="<<fFirstMomentMap[index]<<" SecondMoment=" <<fSecondMomentMap[index]
	// << " var="<<fVariance[index]<<" std="<<fStandardDeviation[index]<<G4endl;
	// }

}

void TsYScorer::GetErrorPropagation()
{
	// calculate statistic error for yF
	for (G4int i = 0; i < yBinNum; i++) 
		yF_var += pow(BinWidth[i],2)*pow((BinLimit[i]+BinLimit[i+1])/2, 2)*fVariance[i];
	yF_std = sqrt(yF_var);

	// calculate statistic error for yD
	for (G4int i = 0; i < yBinNum; i++)
	{
		G4double aa = BinWidth[i]*pow((BinLimit[i]+BinLimit[i+1])/2, 2);
		yD_var += pow(aa/yF,2)*fVariance[i] + pow(aa*hfy[i]/(yF*yF),2)*yF_var;
	}
	yD_std = sqrt(yD_var);

	// calculate statistic error for yd(y), yf(y), d(y)
	for (G4int i = 0; i < yBinNum; i++)
	{
		G4double yi2 = pow((BinLimit[i]+BinLimit[i+1])/2, 2);
		ydy_var[i] = pow(yi2/yF, 2)*fVariance[i]+pow(yi2*hfy[i]/(yF*yF), 2)*yF_var;
		ydy_std[i] = sqrt(ydy_var[i]);   

		yfy_var[i] =  yi2*fVariance[i];
		yfy_std[i] = sqrt(yfy_var[i]);

		G4double yi = (BinLimit[i]+BinLimit[i+1])/2;
		dy_var[i] = pow(yi/yF,2)*fVariance[i] + pow(yi*hfy[i]/(yF*yF), 2)*yF_var;
		dy_std[i] = sqrt(dy_var[i]);
	}
}
