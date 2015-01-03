
#include "FileList.h"
#include <fstream>

void Run()
{
  FileStat_t x;
  TString libName = "libDmpBT2014AMS";
  if(gSystem->GetPathInfo(libName,x)){
    //gROOT->ProcessLine(".include $DMPSWSYS/include"); // same as the next line
    gInterpreter->AddIncludePath("$DMPSWSYS/include");
    gSystem->Load("$DMPSWSYS/lib/libDmpEvent.so");
    //gSystem->Load("$DMPSWSYS/lib/libDmpKernel.so");
    //gSystem->Load("$DMPSWWORK/lib/libDmpEvtAms.so");
    //gSystem->Load("$DMPSWWORK/lib/libDmpEventRaw.so");
    gSystem->CompileMacro("./FileList.h","k",libName);
  }else{
    gSystem->Load(libName);
  }
}

//-------------------------------------------------------------------
enum FileType{
  electron = 0,
  muon = 1,
  photon = 2,
};

//-------------------------------------------------------------------
int Load(short type, TString InFN="NO")
{
  TString f = "NO";
  switch (type){
    case electron:
      f = "./AMS/SPS/Combine_run_1416070809_ANC_366.root";
      break;
    case muon:
      f = "./AMS/SPS/Combine_run_1416338929_ANC_476.root";
      break;
    case photon:
      f = "./AMS/SPS/Combine_run_1416276173_ANC_451.root";
      break;
  }
  return Conf::LoadInput(f);
}

//-------------------------------------------------------------------
void test()
{
  AMS::Performance::Clusters();
  AMS::Alignment::SingleStrack_S_Side();
  //Tracking::Plots(nEvt);
}

