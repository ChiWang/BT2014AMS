/*
 *  $Id: Debug.C, 2015-01-03 15:38:47 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 03/01/2015
//-------------------------------------------------------------------
 *  Usage:
 *      root -l Debug.C
 *      root[1] compile(aFile)
//-------------------------------------------------------------------
*/

void Debug(){
    gInterpreter->AddIncludePath("$DMPSWSYS/include");
    //gInterpreter->AddIncludePath("./");
    gSystem->Load("$DMPSWSYS/lib/libDmpEvent.so");
}

bool compile(TString file){
    TString libName = "lib";
    libName += file;
    libName.Remove(libName.Length()-2);
    gSystem->CompileMacro(file,"k",libName);
}

void Test(){
  AMS::Common::LoadInput("./Rec0/ALL/DAMPE_AMS_ANC_20141118_154848_Rec0.root");
  AMS::DataQuality();
}

