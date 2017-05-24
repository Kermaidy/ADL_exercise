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
 *
 * DESCRIPTION:
 *
 *  A class to calculate traces via ADL in Ge detectors for each energy deposition
 *  
 *
 * AUTHOR:  Y. Kermaidic
 *
 * REVISION: MM-DD-YYYY
 *
 *
 *
 */

#ifndef _GEADLDETECTORTRACE_HH
#define _GEADLDETECTORTRACE_HH

#include <sstream>      // std::ostringstream

// MGDO includes
#include "MGTWaveform.hh"

// ADL includes
#include "ADL.h"

class GEADLDetectorTrace
{
public:
  //default constructor
  GEADLDetectorTrace(int,int);

  //copy constructor
  GEADLDetectorTrace(const GEADLDetectorTrace &);

  //destructor
  ~GEADLDetectorTrace();

  // public functions
  std::string GetSetupFile();
  void SetSetupFile(int);
  void ConfigureADL(int,std::string);
  void SetPositionOffset(double, double, double, double);
  void CreateADLevent();
  void DeleteADLevent();
  void SetWaveformAttribute(double,double,double,double);
  void SetAuxWaveformAttribute(double,double,double,double);
  double SetADLhits(int, Float_t*, Float_t*, Float_t*, Float_t*, Int_t*, Float_t*, Int_t*);
  double SetADLhits(int, Float_t*, Float_t*, Float_t*, Float_t*, Int_t*);
//  double SetADLhits(int, Float_t&, Float_t&, Float_t&, Float_t&, Int_t&);
  int CalculateTrace(std::string);
  int SetADLWaveform(MGTWaveform*);
  int SetADLauxWaveform(MGTWaveform*);
  int GetTraceDim();
  int GetCenter();
  double** GetElectronPath();
  double** GetHolePath();

private:

  struct SIMION_PA *ADL_Epot[40];
  struct SIMION_PA *ADL_Wpot[40];
  struct SIMION_PA *ADL_Stru[40];

  double GridSize[40];
  double Center[40];
  double Height[40];

  std::string detector_setupfile; // ADL configuration file
  int detector_channel; // Detector channel
  int debugADL;

  double inverted; // Correct ADL hits position for inverted detectors
  double xoffset;  // Transpose MaGe coordinate
  double yoffset;  //      into the ADL referential
  double zoffset;

  double xcenter;   // Detector center in cm in ADL coordinate 
  double ycenter;   // Detector center in cm in ADL coordinate 
  double height;   // Detector height in cm 
  double gridsize; // grid size used for det. potentials computation

  struct ADL_EVENT *ADL_evt; // Define ADL event to set hits position and calculate carriers path
  
  double Amplitude;     // Signal amplitude in ADC
  double Baseline;      // Signal baseline
  double RMS_noise;     // Baseline noise amplitude (assumed to be gaussian)
  double wfPreTrigger;  // Time before trigger (usually half of the signal length)

  double AuxAmplitude;
  double AuxBaseline;
  double AuxRMS_noise;
  double AuxwfPreTrigger;
};

#endif
