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


#ifndef _GEOUTPUTGERMANIUMARRAY_WADL_HH
#define _GEOUTPUTGERMANIUMARRAY_WADL_HH

// MGDO includes
#include "MGTWaveform.hh"
#include "MGTEvent.hh"
#include "MGTRun.hh"

// GELATIO includes
#include "GETGERDADigitizerData.hh"

// ADL includes
#include "GEADLDetectorTrace.hh"

//---------------------------------------------------------------------------//
// ROOT declarations
class TFile;
class Th1D;
class TTree;

//---------------------------------------------------------------------------//


class GEOutputGermaniumArray_wADL
{
public:
  //default constructor
  GEOutputGermaniumArray_wADL();

  //destructor
  ~GEOutputGermaniumArray_wADL();

  void DefineSchema();
  void RunSimulation(int);
  int WriteOutputs(TFile*);

  //private  members
private:

  void SimulatePulse(int,int,double,double,double);   // ADL-4.2

  static const int MAX_NTRACE=2000;
  static const int MAX_NHITS=10;
  static const int NDET=2;

  // hits
  Int_t    hits_totnum;
  Float_t  hits_tote;
  Float_t  hits_edep[MAX_NHITS];
  Float_t  hits_xpos[MAX_NHITS];
  Float_t  hits_ypos[MAX_NHITS];
  Float_t  hits_zpos[MAX_NHITS];
  Int_t    hits_iddet[MAX_NHITS];  // which Ge detector this hit is in

  // hits : ADL informations
  Float_t  hits_ADLpos[MAX_NHITS];
  Int_t  hits_isOut[MAX_NHITS];

  Int_t    trace_totnum;
  Int_t    trace_iddet;
  Int_t    eventnumber;
  Float_t  trace_xposE[MAX_NTRACE];
  Float_t  trace_yposE[MAX_NTRACE];
  Float_t  trace_zposE[MAX_NTRACE];
  Float_t  trace_xposH[MAX_NTRACE];
  Float_t  trace_yposH[MAX_NTRACE];
  Float_t  trace_zposH[MAX_NTRACE];

  bool fSimulateTraces;
  bool fRecordADLOutPos;
  bool fRecordADLTraces;

    TTree* nT;

  /***********************/
  /*** ADL simulation  ***/
  /***********************/

  GEADLDetectorTrace* ADLDetector;

  int isAlreadyRead[NDET];   // Flag to read only once E/W/Stru potential from files (then stored into memory)
  int traceCalculated[NDET]; // Flag to calculate trace in Ge det only once per event
  int setupADLdetPos;    // Used to define the detectors positon offset in ADL only for the 1st event

    TTree* MGTree;
  MGTRun* theRun;
  MGTEvent* event;
  MGTWaveform* waveform;
  MGTWaveform* auxwaveform;
  GETGERDADigitizerData* digiData;

  size_t fPreTrigger;
  unsigned long long fTimestamp;
  unsigned long long fDecimalTimestamp;
  bool IsPulseInverted;
  unsigned int TriggerNumber;
  int fEventNumber;
  bool fMuVetoed;
  unsigned int fMuVetoSample;
  int wfLength;
  int wfPreTrigger;
  int auxwfLength;
  int auxwfPreTrigger;

  Double_t theTime;
  double SamplingFrequency; // GHz
  double AuxSamplingFrequency; // GHz
  size_t iChannel;
  int validChannelCounter;

  double FEP_kev;
  double FEP_ADC;
  double Baseline;
  double RMS_noise;
  
  /***********************/

};

#endif

