
//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//	This code implementation is the intellectual property of the         //
//	MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//	intellectual property of the RD44 GEANT4 collaboration.              //
//                                                                           //
//                        *********************                              //
//                                                                           //
//    Neither the authors of this software system, nor their employing       //
//    institutes, nor the agencies providing financial support for this      //
//    work  make  any representation or  warranty, express or implied,       //
//    regarding this software system or assume any liability for its use.    //
//    By copying, distributing or modifying the Program (or any work based   //
//    on on the Program) you indicate your acceptance of this statement,     //
//    and all its terms.                                                     //
//                                                                           //
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//---------------------------------------------------------------------------//
/**
 * CLASS IMPLEMENTATION:
 *   GEADLDetectorTrace
 *
 * AUTHOR:  Y. Kermaidic
 * CONTACT: yoann.kermaidic@mpi-hd.pmg.de
 *
 *
 */

#include "GEADLDetectorTrace.hh"
#include "TRandom.h"
#include <fstream>

using namespace std;

//---------------------------------------------------------------------------//
GEADLDetectorTrace::GEADLDetectorTrace(int debug, int channel){

  debugADL = debug;
  SetADLDebug(0);
  detector_channel = channel;

}

//---------------------------------------------------------------------------//

GEADLDetectorTrace::~GEADLDetectorTrace() {
}

//---------------------------------------------------------------------------//

void GEADLDetectorTrace::SetSetupFile(int channel){

  detector_channel = channel;

  std::ostringstream oss;
  if(channel <10) oss << 0;
  else oss.str("");
  oss << channel;

  if(channel == 666) detector_setupfile = "configfiles/Det_HADES/ICOAX.txt";
  else detector_setupfile = "configfiles/Det_" + oss.str() + "/Det.txt";
}

std::string GEADLDetectorTrace::GetSetupFile() {return detector_setupfile;}

void GEADLDetectorTrace::ConfigureADL(int isalreadyread,std::string setupfile)
{
  ADLSetup(const_cast<char*>(setupfile.c_str()));
  if(!isalreadyread){
    SetupFields(const_cast<char*>(setupfile.c_str()));

    ADL_Epot[detector_channel] = GetADLEpot();
    ADL_Wpot[detector_channel] = GetADLWpot();
    ADL_Stru[detector_channel] = GetADLStru();

    GridSize[detector_channel] = GetSimionGridSize();
    Center[detector_channel]   = GetSimionCenter();
    Height[detector_channel]   = GetSimionHeight();
  }
  else{
    SetADLEpot(ADL_Epot[detector_channel]);
    SetADLWpot(ADL_Wpot[detector_channel]);
    SetADLStru(ADL_Stru[detector_channel]);

    SetSimionGridSize(GridSize[detector_channel]);
    SetSimionCenter(Center[detector_channel]);
    SetSimionHeight(Height[detector_channel]);
  }
  
  if(debugADL) 
    std::cout << detector_channel << "       " << 
    ADL_Epot[detector_channel]->h.nx  << " " << 
    ADL_Epot[detector_channel]->h.ny  << " " << 
    ADL_Epot[detector_channel]->h.nz  << " " <<  
    ADL_Wpot[detector_channel]->h.nx  << " " << 
    ADL_Wpot[detector_channel]->h.ny  << " " << 
    ADL_Wpot[detector_channel]->h.nz  << " " <<
    ADL_Stru[detector_channel]->h.nx  << " " << 
    ADL_Stru[detector_channel]->h.ny  << " " << 
    ADL_Stru[detector_channel]->h.nz  << std::endl;
  
  if(debugADL) StatusFields();
}

void GEADLDetectorTrace::SetPositionOffset(double inv, double x, double y, double z){

  inverted = inv;

  xoffset = x;
  yoffset = y;
  zoffset = z;

  gridsize = GetSimionGridSize();
  xcenter  = GetSimionCenter()    - xoffset/10.;           //center of detector in cm
  ycenter  = GetSimionCenter()    - yoffset/10.;           //center of detector in cm
  height   = inverted*GetSimionHeight()/2. - zoffset/10.;           //bottom of detector in cm
  
  if(debugADL){
    std::cout << "DEBUG: Set position offset : " << x << " " << y << " " << z/2. << std::endl;  
    std::cout << "DEBUG: ADL detector gridsize : " << gridsize          << std::endl;
    std::cout << "DEBUG: ADL detector center   : " << GetSimionCenter() << std::endl;
    std::cout << "DEBUG: ADL detector height   : " << GetSimionHeight() << std::endl;
    std::cout << "DEBUG: Position offset set " << std::endl;
    std::cout << "DEBUG: ADL to MaGe xcenter   : " << xcenter << std::endl;
    std::cout << "DEBUG: ADL to MaGe ycenter   : " << ycenter << std::endl;
    std::cout << "DEBUG: ADL to MaGe height    : " << height  << std::endl;
  }
}

void GEADLDetectorTrace::CreateADLevent(){
  ADL_evt = new_event();
  ADL_evt->HP.T0= 1.000;               //Time the interaction occurs in the trace
}

void GEADLDetectorTrace::DeleteADLevent(){
  delete ADL_evt;
}

void GEADLDetectorTrace::SetWaveformAttribute(double wfpretrigger, double baseline, double amplitude, double rms_noise){
  wfPreTrigger = wfpretrigger;
  Baseline = baseline;
  Amplitude = amplitude;
  RMS_noise = rms_noise;
}

void GEADLDetectorTrace::SetAuxWaveformAttribute(double wfpretrigger, double baseline, double amplitude, double rms_noise){
  AuxwfPreTrigger = wfpretrigger;
  AuxBaseline = baseline;
  AuxAmplitude = amplitude;
  AuxRMS_noise = rms_noise;
}

double GEADLDetectorTrace::SetADLhits(int hits_totnum, Float_t*hits_edep, Float_t*hits_xpos, Float_t*hits_ypos, Float_t*hits_zpos, Int_t*hits_iddet){
  
  double ETotDet = 0;
  int j = 0;
  
  for(Int_t i = 0;i<hits_totnum;i++){
    
    if(hits_iddet[i] == detector_channel){ // Consider only hits in the given detector.
      //Fill in the Hit Pattern (HP):
      double P0[4]={0,hits_xpos[i] + xcenter,hits_ypos[i] + ycenter,inverted*(hits_zpos[i] + height)};
      
      if(IsInDetector(P0)){
	if(debugADL) std::cout << "Hits in detector " << detector_channel << std::endl;
	ETotDet += hits_edep[i];
	ADL_evt->HP.Eint[j]  =hits_edep[i];             //Energy of interaction
	ADL_evt->HP.Pos[j][0]=hits_xpos[i] + xcenter;	  //Position where this interaction occures in the ADL referential
	ADL_evt->HP.Pos[j][1]=hits_ypos[i] + ycenter;
	ADL_evt->HP.Pos[j][2]=inverted*(hits_zpos[i] + height);
	
	j++; // Only iterate on points in the detector
      }
      else if(debugADL) std::cout << "Hits not in detector " << detector_channel << std::endl;
      
      if(debugADL) std::cout << "    MaGe ref. hits position : " << hits_xpos[i] << " " << hits_ypos[i] << " " << hits_zpos[i] << std::endl;
      if(debugADL) std::cout << "    ADL  ref. hits position : " << P0[1] << " " << P0[2] << " " << P0[3] << std::endl;
      if(debugADL) std::cout << "    Cen. ref. hits position : " << P0[1]-GetSimionCenter() << " " << P0[2]-GetSimionCenter() << " " << P0[3] << std::endl;
    }
  }
  return ETotDet;
}

double GEADLDetectorTrace::SetADLhits(int hits_totnum, Float_t*hits_edep, Float_t*hits_xpos, Float_t*hits_ypos, Float_t*hits_zpos, Int_t*hits_iddet, Float_t*hits_ADLpos, Int_t*hits_isOut){

  double ETotDet = 0;
  int j = 0;

  for(Int_t i = 0;i<hits_totnum;i++){
    if(hits_iddet[i] == detector_channel){
      //Fill in the Hit Pattern (HP):
      double P0[4]={0,hits_xpos[i] + xcenter,hits_ypos[i] + ycenter,inverted*(hits_zpos[i] + height)};
      hits_ADLpos[i] = GetDetectorPos(P0);
      if(IsInDetector(P0)){
	hits_isOut[i] = 0;
	ETotDet += hits_edep[i];
	ADL_evt->HP.Eint[j]  =hits_edep[i];             //Energy of interaction
	ADL_evt->HP.Pos[j][0]=hits_xpos[i] + xcenter;	  //Position where this interaction occures in the ADL referential
	ADL_evt->HP.Pos[j][1]=hits_ypos[i] + ycenter;
	ADL_evt->HP.Pos[j][2]=inverted*(hits_zpos[i] + height);
	
	j++; // Only iterate on points in the detector
      }
      else{
	// Hit not in detector. Set outer position to hit position 
	hits_isOut[i] = 1;
	if(debugADL) std::cout << "Hits not in detector " << detector_channel << std::endl;
	if(debugADL) std::cout << "    Hits position " << hits_ADLpos[i] << std::endl;
	if(debugADL) std::cout << "    MaGe ref. hits position : " << hits_xpos[i] << " " << hits_ypos[i] << " " << hits_zpos[i] << std::endl;
	if(debugADL) std::cout << "    ADL ref. hits position  : " << P0[1] << " " << P0[2] << " " << P0[3] << std::endl;
	if(debugADL) std::cout << "    Cen. ref. hits position : " << P0[1]-GetSimionCenter() << " " << P0[2]-GetSimionCenter() << " " << P0[3] << std::endl;
      }
    }
    else if(hits_isOut[i] < 0) hits_isOut[i] = 2; // If hits_isOut has not been set yet and is not related to detector channel
  }
  return ETotDet;
}

int GEADLDetectorTrace::CalculateTrace(std::string setupfile){
  if(debugADL) StatusTraces(ADL_evt);
  CalculateTraces(ADL_evt);

  return 0;
}

int GEADLDetectorTrace::SetADLWaveform(MGTWaveform* waveform) {

  int TraceLength = GetDIMT();
  int adl_iter = 0;
  
  if(TraceLength>0){
    for (size_t i=0;i<waveform->GetLength();i++){
      if(i<wfPreTrigger) (*waveform)[i] = 4. * Baseline + 2. * gRandom->Gaus(0,RMS_noise);
      else if(i<wfPreTrigger + TraceLength/4 - 5){ 
	(*waveform)[i] = 4. * Baseline + 2. * gRandom->Gaus(0,RMS_noise) + Amplitude * ((ADL_evt->TD).Tr[0][adl_iter]
											+ (ADL_evt->TD).Tr[0][adl_iter+1] 
											+ (ADL_evt->TD).Tr[0][adl_iter+2] 
											+ (ADL_evt->TD).Tr[0][adl_iter+3]);
	adl_iter += 4;
      }
      else (*waveform)[i] = 4. * Baseline + 2. * gRandom->Gaus(0,RMS_noise) + 4. * Amplitude * (ADL_evt->TD).Tr[0][TraceLength];
    }
  }
  else return 1;
  return 0;
}

int GEADLDetectorTrace::SetADLauxWaveform(MGTWaveform* waveform) {
  
  int TraceLength = GetDIMT();
  int adl_iter = 0;
  
  if(TraceLength>0){
    for (size_t i=0;i<waveform->GetLength();i++){
      if(i<AuxwfPreTrigger) (*waveform)[i] = AuxBaseline + gRandom->Gaus(0,AuxRMS_noise);
      else if(i<AuxwfPreTrigger + TraceLength){ 
	if(debugADL) std::cout << adl_iter 
				    << " " << (ADL_evt->TD).Tr[0][adl_iter] 
				    << " " << GetNUMRES_XYZh()[adl_iter][1] 
				    << " " << GetNUMRES_XYZh()[adl_iter][2] 
				    << " " << GetNUMRES_XYZh()[adl_iter][3] 
				    << std::endl;
	(*waveform)[i] = AuxBaseline + gRandom->Gaus(0,AuxRMS_noise) + AuxAmplitude * (ADL_evt->TD).Tr[0][adl_iter];
	adl_iter ++;
      }
      else (*waveform)[i] = AuxBaseline + gRandom->Gaus(0,AuxRMS_noise) + AuxAmplitude * (ADL_evt->TD).Tr[0][TraceLength];					 
    }
  }
  else return 1;
  return 0;
}

int GEADLDetectorTrace::GetTraceDim(){return GetDIMT();}
int GEADLDetectorTrace::GetCenter(){return GetSimionCenter();}
double** GEADLDetectorTrace::GetElectronPath(){return GetNUMRES_XYZe();}
double** GEADLDetectorTrace::GetHolePath(){return GetNUMRES_XYZh();}
