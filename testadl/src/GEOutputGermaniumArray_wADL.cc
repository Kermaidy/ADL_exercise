
//---------------------------------------------------------------------------//
//bb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nubb0nu//
//                                                                           //
//                                                                           //
//                         MaGe Simulation                                   //
//                                                                           //
//      This code implementation is the intellectual property of the         //
//      MAJORANA and Gerda Collaborations. It is based on Geant4, an         //
//      intellectual property of the RD44 GEANT4 collaboration.              //
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
//
// include files for ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TObject.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector2.h"

#include <sstream>
#include <string>
#include <typeinfo>

//---------------------------------------------------------------------------//

#include "GEOutputGermaniumArray_wADL.hh"
using namespace std;

//---------------------------------------------------------------------------//
GEOutputGermaniumArray_wADL::GEOutputGermaniumArray_wADL():
      fSimulateTraces(true),
      fRecordADLOutPos(true),
      fRecordADLTraces(true)
{
  MGTree = 0;
  theRun = 0;
  event = 0;
  waveform = 0;
  auxwaveform = 0;
  digiData = 0;
}

GEOutputGermaniumArray_wADL::~GEOutputGermaniumArray_wADL()
{
 
}

void GEOutputGermaniumArray_wADL::DefineSchema()
{
    nT = new TTree("fTree","fTree");
    
   // hits in sensitive detectors
    nT->Branch("hits_totnum",&hits_totnum,"hits_totnum/I");
    nT->Branch("hits_tote",&hits_tote,"hits_tote/F");
    nT->Branch("hits_edep",hits_edep,"hits_edep[hits_totnum]/F");
    nT->Branch("hits_xpos",hits_xpos,"hits_xpos[hits_totnum]/F");
    nT->Branch("hits_ypos",hits_ypos,"hits_ypos[hits_totnum]/F");
    nT->Branch("hits_zpos",hits_zpos,"hits_zpos[hits_totnum]/F");
    if( fRecordADLOutPos ) nT->Branch("hits_ADLpos", hits_ADLpos, "hits_ADLpos[hits_totnum]/F");
    if( fRecordADLOutPos ) nT->Branch("hits_isOut", hits_isOut, "hits_isOut[hits_totnum]/I");
 
    if( fRecordADLTraces ) nT->Branch("trace_totnum", &trace_totnum, "trace_totnum/I");
    if( fRecordADLTraces ) nT->Branch("trace_iddet", &trace_iddet, "trace_iddet/I");
    if( fRecordADLTraces ) nT->Branch("trace_xposE", trace_xposE, "trace_xposE[trace_totnum]/F");
    if( fRecordADLTraces ) nT->Branch("trace_yposE", trace_yposE, "trace_yposE[trace_totnum]/F");
    if( fRecordADLTraces ) nT->Branch("trace_zposE", trace_zposE, "trace_zposE[trace_totnum]/F");
    if( fRecordADLTraces ) nT->Branch("trace_xposH", trace_xposH, "trace_xposH[trace_totnum]/F");
    if( fRecordADLTraces ) nT->Branch("trace_yposH", trace_yposH, "trace_yposH[trace_totnum]/F");
    if( fRecordADLTraces ) nT->Branch("trace_zposH", trace_zposH, "trace_zposH[trace_totnum]/F");

    /// Optional optimization for long running jobs generating large files.
    /// By default ROOT flushes the buffers and writes the header of the TTree
    /// to the file when the amount of data written to the file exceeds 300 MB.
    /// For long simulations with small amount of data per event this can cause
    /// a lot of events to be lost in case of forced job termination.
    /// By changing this value, the buffers and the header are flushed for different sizes

    if(fSimulateTraces){
      for(int i=0; i < NDET; i++) isAlreadyRead[i] = 0;
      setupADLdetPos = 0;
      
      // Set waveform amplitudes/noise
      FEP_kev = 2.6145;
      Baseline = 5604.0;
      FEP_ADC = 6046;
      RMS_noise = 5.3;

      // ADL detector initialization
      ADLDetector = new  GEADLDetectorTrace(0,0);

      // Tier1 tree definition
      MGTree = new TTree("MGTree","MGTree");
      
      // Waveform branch
      MGTree->Branch("event",&event, 32000,-99);
      
      // Run info
      MGRunType::RunType fRunType;
      fRunType = MGRunType::kData;
      std::string theDAQ = "Struck";
      int fRunNumber = 1;
      time_t fStartTime = 1479394213;
      time_t fStopTime = 1479394213.000160000;
      std::string theRunDescription = "No description";
      
      theRun = new MGTRun();
      theRun->SetRunType(fRunType);
      theRun->SetRunNumber(fRunNumber);
      theRun->SetRunDescription(theRunDescription);
      theRun->SetStartTime(fStartTime);
      theRun->SetStopTime(fStopTime);
      //The Stop time is available at the end of run
      theRun->SetParentDAQLabel(theDAQ);
      
      // Set MGDO digitizer
      fPreTrigger = 100;
      fTimestamp = 1479394213;
      fDecimalTimestamp = 0;
      IsPulseInverted = false;
      TriggerNumber = 1;
      fEventNumber = 1;
      fMuVetoed = false;
      fMuVetoSample = 0;
      wfLength = 4000;
      wfPreTrigger = 1900;
      auxwfLength = 1000;
      auxwfPreTrigger = 300;
      
      theTime = 1479394213;
      SamplingFrequency = 0.025; // GHz
      AuxSamplingFrequency = 0.1; // GHz
      iChannel = 0;
      validChannelCounter = 0;
  }
}

void GEOutputGermaniumArray_wADL::RunSimulation(int EventNumber){

    double xDet[NDET];
    double yDet[NDET];
    double zDet[NDET];
    double emean = 2035;
    double esigma = 5;
    
    eventnumber = EventNumber;
    hits_tote = 0.0;
 
    for(int i=0; i < NDET; i++){
      traceCalculated[i] = 0;
      xDet[i] = 0;
      yDet[i] = 0;
      zDet[i] = 0;
    }

    hits_totnum = 1;
    
    TRandom3* rand = new TRandom3(0);
    
    //Set hits informations
    for (int i=0; i<hits_totnum; i++) {
      rand->SetSeed(0);

      hits_edep[i] = rand->Gaus(emean,esigma);
      hits_tote += hits_edep[i];
      hits_xpos[i] = rand->Uniform(-4,4);
      hits_ypos[i] = rand->Uniform(-4,4);
      hits_zpos[i] = rand->Uniform(-2,2);
      hits_iddet[i] = roundf(rand->Uniform(0,NDET-1));
      
      if(fSimulateTraces){ // Get detector position
        hits_isOut[i] = -1;
        xDet[hits_iddet[i]] = 0;
        yDet[hits_iddet[i]] = 0;
        zDet[hits_iddet[i]] = 0;
//          xDet[hits_iddet[i]] = theDetDB->GetCrystalPosition(fOuterCounter,fInnerCounter).x();
//          yDet[hits_iddet[i]] = theDetDB->GetCrystalPosition(fOuterCounter,fInnerCounter).y();
//          zDet[hits_iddet[i]] = theDetDB->GetCrystalPosition(fOuterCounter,fInnerCounter).z();
      }
    }

    /*
      std::cout << "E   : " << hits_edep[0]  << std::endl;
      std::cout << "x   : " << hits_xpos[0]  << std::endl;
      std::cout << "y   : " << hits_ypos[0]  << std::endl;
      std::cout << "z   : " << hits_zpos[0]  << std::endl;
      std::cout << "Det : " << hits_iddet[0] << std::endl;
    */

    if(fSimulateTraces){
      // Set MGDO event structure
      event = new MGTEvent();
      MGEventType::EventType fEventType;
      fEventType = MGEventType::kBaseline;
      event->SetAuxWaveformArrayStatus(true);
      event->SetTotalNumberOfIDs(NDET);
      event->SetAllActiveIDs(true);
      event->InitializeArrays("GETGERDADigitizerData",1);
      event->SetBypassStreamerArray(kFALSE); //! Event type info
      event->SetEventType(fEventType);
      event->SetTime(theTime);
      event->SetWFEncScheme(MGTWaveform::kDiffVarInt);
      event->SetAuxWFEncScheme(MGTWaveform::kDiffVarInt);
      event->SetETotal(hits_tote);
      validChannelCounter = 0;
    }

    //Simulate pulses with ADL once per event
    for (int i=0; i<hits_totnum; i++) {
      if(fSimulateTraces && traceCalculated[hits_iddet[i]] == 0){
          SimulatePulse(isAlreadyRead[hits_iddet[i]],hits_iddet[i],xDet[hits_iddet[i]],yDet[hits_iddet[i]],zDet[hits_iddet[i]]);
          traceCalculated[hits_iddet[i]] = 1;
          isAlreadyRead[hits_iddet[i]] = 1;
      }
    }
    nT->Fill();
}

void GEOutputGermaniumArray_wADL::SimulatePulse(int isalreadyread,int channel, double x,double y,double z){
 
  int debugADL = 0;

  double ETotDet = 0;

  ADLDetector->SetSetupFile(channel);
  ADLDetector->ConfigureADL(isalreadyread,ADLDetector->GetSetupFile()); // ADL-4.2
  if(channel == 9 || channel == 14 ||channel == 16 ||channel == 20 ||channel == 22 ||channel == 33) // Consider inverted BEGe
    ADLDetector->SetPositionOffset(-1,x,y,z);
  else ADLDetector->SetPositionOffset(1,x,y,z);

  ADLDetector->CreateADLevent();
  if(debugADL) std::cout << "DEBUG: ADL event created" << std::endl;

  if(fRecordADLOutPos)
    ETotDet = ADLDetector->SetADLhits(hits_totnum,hits_edep,hits_xpos,hits_ypos,hits_zpos,hits_iddet,hits_ADLpos,hits_isOut);
  else
    ETotDet = ADLDetector->SetADLhits(hits_totnum,hits_edep,hits_xpos,hits_ypos,hits_zpos,hits_iddet);
  if(debugADL) std::cout << "DEBUG: ADL hits set" << std::endl;
  if(debugADL) std::cout << "DEBUG: Deposited energy in channel " << channel << " is " << ETotDet << "/" << hits_tote << " MeV" << std::endl;

  if(hits_tote > 0.){

    if(ADLDetector->CalculateTrace(ADLDetector->GetSetupFile())) std::cerr<< "Failed to calculate trace" <<std::endl;
    if(debugADL) std::cout << "DEBUG: ADL calculate trace" << std::endl;
    
    MGWaveformTag::EWaveformTag fWaveformTag = MGWaveformTag::kNormal;
    
    //Fill digiData. 
    digiData = new ((*(event->GetDigitizerData()))[validChannelCounter]) GETGERDADigitizerData(); 
    digiData->SetClockFrequency(SamplingFrequency);
    
    digiData->SetPretrigger(fPreTrigger);
    digiData->SetTimeStamp(fTimestamp);
    digiData->SetDecimalTimeStamp(fDecimalTimestamp);
    digiData->SetIsInverted(IsPulseInverted);
    digiData->SetTriggerNumber(TriggerNumber);
    digiData->SetID(channel);
    digiData->SetEnergy(ETotDet);
    digiData->SetEventNumber(eventnumber);
    digiData->SetIsMuVetoed(fMuVetoed);
    digiData->SetMuVetoSample(fMuVetoSample);
    digiData->SetWaveformTag(fWaveformTag);
    if(debugADL) std::cout << "DEBUG: ADL digitizer set" << std::endl;
    
    // Set MGDO waveform 
    waveform = new ((*(event->GetWaveforms()))[validChannelCounter]) MGTWaveform(NULL,0,SamplingFrequency,0.0,MGWaveform::kADC,0);
    waveform->SetLength(wfLength);
    waveform->SetTOffset(0.);
    waveform->SetID(channel);
    
    ADLDetector->SetWaveformAttribute(wfPreTrigger,Baseline,FEP_ADC/FEP_kev,RMS_noise);
    if(debugADL) std::cout << "DEBUG: Waveform attribute set" << std::endl;
    
    if(ADLDetector->SetADLWaveform(waveform)) std::cerr<< "Failed to set waveform" <<std::endl;
    if(debugADL) std::cout << "DEBUG: Waveform set" << std::endl;
    
    if (event->GetAuxWaveformArrayStatus()){
      auxwaveform = new ((*(event->GetAuxWaveforms()))[validChannelCounter])
	MGTWaveform(NULL,0,AuxSamplingFrequency,0.0,MGWaveform::kADC,0);     
      auxwaveform->SetSamplingFrequency(AuxSamplingFrequency);
      auxwaveform->SetID(channel);
      auxwaveform->SetTOffset(0.);      
      auxwaveform->SetLength(auxwfLength);
      
      ADLDetector->SetAuxWaveformAttribute(auxwfPreTrigger,Baseline,FEP_ADC/FEP_kev,RMS_noise);
      if(ADLDetector->SetADLauxWaveform(auxwaveform)) std::cerr<< "Failed to set aux waveform" <<std::endl;
    }

    validChannelCounter++;

    if(fRecordADLTraces){
      double** ePath = ADLDetector->GetElectronPath();  // Get matrix containing e- path in (x,y,z) coord.
      double** hPath = ADLDetector->GetHolePath();      // Get matrix containing h  path in (x,y,z) coord.
      trace_totnum = ADLDetector->GetTraceDim();
      trace_iddet = channel;
      for(int i = 0;i<trace_totnum;i++){
          trace_xposE[i] = ePath[i][1] - ADLDetector->GetCenter(); // center e/h traces around 0 in the (x,y) plane
          trace_yposE[i] = ePath[i][2] - ADLDetector->GetCenter();
          trace_zposE[i] = ePath[i][3];
          trace_xposH[i] = hPath[i][1] - ADLDetector->GetCenter();
          trace_yposH[i] = hPath[i][2] - ADLDetector->GetCenter();
          trace_zposH[i] = hPath[i][3];
      }
    }

    MGTree->Fill();
    if(debugADL) std::cout << "Print MGTree : " << MGTree->GetEntries() << " events recorded" << std::endl;
  }
  ADLDetector->DeleteADLevent();
}

int GEOutputGermaniumArray_wADL::WriteOutputs(TFile* fOutputFile)
{
  fOutputFile->cd();
  if(fOutputFile->IsOpen()){
    MGTree->Write();
    nT->Write();

    fOutputFile->Close();
    return 1;
  }
  else
    return 0;
}
