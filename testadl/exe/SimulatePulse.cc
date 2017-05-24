#include "GEOutputGermaniumArray_wADL.hh"
#include "GEADLDetectorTrace.hh"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

using namespace std;

int main(int argc, const char* argv[])
{

 cout << " " << endl;
 cout << " " << endl;
 cout << " THE ADL TESTER " << endl;
 cout << " " << endl;
 cout << " " << endl;

 int Nevts = 1000;
 string rootfilename = "./ADLTest-Tier1.root";
 string simNumber;
 if(argc>1){
	simNumber = argv[1];
	rootfilename = "./ADLTest_" + simNumber + "-Tier1.root";
 }

 GEOutputGermaniumArray_wADL ADLtest;

 cout << "Initialize the simulation " << endl;
 ADLtest.DefineSchema();

 for(int i = 0; i<Nevts;i++){
	cout << "\rSimulating event : " << i+1 << "/" << Nevts << flush;	
	ADLtest.RunSimulation(i);
 }

 std::cout << " " << std::endl;
 std::cout << "Open outputs ROOT file" << endl;

 TFile* fOutputFile = new TFile(rootfilename.c_str(),"RECREATE");

 if(ADLtest.WriteOutputs(fOutputFile)){
   cout << "Simulation ended properly " << endl;
   cout << " " << std::endl;
 }
 else{  
   cout << "Data not written. Exit. " << endl;
   cout << " " << std::endl;
 }
}
