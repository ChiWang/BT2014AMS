/*
 *  $Id: DataQuality.C, 2015-01-03 15:42:05 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 30/11/2014
*/

#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "TTree.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "AMSSetup.h"

#ifndef AMS_DataQuality_C
#define AMS_DataQuality_C

namespace AMS{
using namespace Common;
  void DataQuality(){  // without any cuts
    vector<TH1F*>   h_clsSeed(NLadder*2);  // cluster seed
    vector<TH1F*>   h_clsNB(NLadder*2);  // cluster numbers
    vector<TH2F*>   h_COG_SNR(NLadder*2);  // COG, center of geometry. GetCoG(). SNR: signal inoise ratio of cluster. GetTotSN()
    vector<TH2F*>   h_COG_Length(NLadder*2);  // length of this cluster: GetLength()
    vector<TH2F*>   h_Seed_Value_SN(NLadder*2);  // GetSeedVal(), GetSeedSN()
    for(short i =0;i<NLadder; ++i){
      for(short j =0;j<2;++j){
        h_clsSeed[i*2+j] = new TH1F(Form("L%d_S%d--cluster seed",i,j),Form("L%d_S%d cluster seed",i,j),1024,0,1024);
        h_clsSeed[i*2+j]->SetLabelSize(0.12);
        h_clsSeed[i*2+j]->SetLabelSize(0.08,"Y");
        h_clsSeed[i*2+j]->SetLineColor(j+3);
        h_clsNB[i*2+j] = new TH1F(Form("L%d_S%d--cluster number",i,j),Form("L%d_S%d cluster number",i,j),8,0,8);
        h_clsNB[i*2+j]->SetLabelSize(0.12);
        h_clsNB[i*2+j]->SetLabelSize(0.08,"Y");
        h_COG_SNR[i*2+j] = new TH2F(Form("L%d_S%d--COG VS SNR",i,j),Form("L%d_S%d CoG VS SNR",i,j),1024*2,j*640,640+j*384,200,0,40);
        h_COG_SNR[i*2+j]->SetLabelSize(0.12);
        h_COG_SNR[i*2+j]->SetLabelSize(0.08,"Y");
        h_COG_Length[i*2+j]=new TH2F(Form("L%d_S%d--COG VS Length",i,j),Form("L%d_S%d CoG VS Length",i,j),1024*2,j*640,640+j*384,30,0,30);
        h_COG_Length[i*2+j]->SetLabelSize(0.12);
        h_COG_Length[i*2+j]->SetLabelSize(0.08,"Y");
        h_Seed_Value_SN[i*2+j]=new TH2F(Form("L%d_S%d--Seed ADC VS Seed SNR",i,j),Form("L%d_S%d Seed ADC VS Seed SNR",i,j),200,0,200,30,0,30);
        h_Seed_Value_SN[i*2+j]->SetLabelSize(0.12);
        h_Seed_Value_SN[i*2+j]->SetLabelSize(0.08,"Y");
      }
    }

    // event loop
    for(_evtID=0;_evtID<_entries;++_evtID){
      _tree->GetEntry(_evtID);
      vector<int> clusBN = ClusterNumberInLadder();
      for(short id=0;id<NLadder;++id){
        for(short s=0;s<2;++s){
          h_clsNB[id*2+s]->Fill(clusBN[id*2+s]);
        }
      }
      int n_cls = _evt_AMS->fAMSCls->GetEntriesFast();
      for(short ic=0;ic<n_cls;++ic){
        Cluster *aCluster = dynamic_cast<Cluster*>(_evt_AMS->fAMSCls->At(ic));
        short order = aCluster->ladder;
        short side = aCluster->side;
        float CoG = aCluster->GetCoG();
        h_COG_SNR[order*2+side]->Fill(CoG,aCluster->GetTotSN());
        h_COG_Length[order*2+side]->Fill(CoG,aCluster->GetLength());
        h_Seed_Value_SN[order*2+side]->Fill(aCluster->GetSeedVal(),aCluster->GetSeedSN());
        h_clsSeed[order*2+side]->Fill(aCluster->GetSeedAdd());
      }
    }

    // draw
    TCanvas *c0 = new TCanvas(_file+"  Cluster seed",_file+"  Cluster seed");
    c0->Divide(1,5,0,0);
    TCanvas *c1 = new TCanvas(_file+"  Cluster number",_file+"  Clusters number");
    c1->Divide(2,5,0,0);
    TCanvas *c2 = new TCanvas(_file+"  Cluster CoG_SNR",_file+"  Cluster CoG_SNR");
    c2->Divide(2,5,0,0);
    TCanvas *c3 = new TCanvas(_file+"  Cluster CoG_Length",_file+"  Cluster CoG_Length");
    c3->Divide(2,5,0,0);
    TCanvas *c4 = new TCanvas(_file+"  Cluster SNR",_file+"  Cluster SNR");
    c4->Divide(2,5,0,0);
    TCanvas *c5 = new TCanvas(_file+"  Cluster Length",_file+"  Cluster Length");
    c5->Divide(2,5,0,0);
    //TCanvas *c6 = new TCanvas(_file+"  Seed ADC and SN",_file+"  Seed ADC and SN");
    //c6->Divide(2,5,0,0);
    TCanvas *c6 = new TCanvas(_file+"  Seed SN",_file+"  Seed SN");
    c6->Divide(2,5,0,0);
    TCanvas *c7 = new TCanvas(_file+"  Seed ADC",_file+"  Seed ADC");
    c7->Divide(2,5,0,0);
    for(short id=0;id<NLadder;++id){
        c0->cd(id+1);
        h_clsSeed[id*2+0]->Draw();
        h_clsSeed[id*2+1]->Draw("same");
      for(short s=0;s<2;++s){
        c1->cd(id*2+s+1);
        h_clsNB[id*2+s]->Draw();
        c2->cd(id*2+s+1);
        h_COG_SNR[id*2+s]->Draw("colz");
        //h_COG_SNR[id][s]->ProfileX()->Draw("same");
        c3->cd(id*2+s+1);
        h_COG_Length[id*2+s]->Draw("colz");
        c4->cd(id*2+s+1);
        TH1D *h_SNR = h_COG_SNR[id*2+s]->ProjectionY();
        h_SNR->SetTitle(Form("L%d_S%d Cluster SNR",id,s));
        h_SNR->SetLabelSize(0.12);
        h_SNR->SetLabelSize(0.08,"Y");
        gPad->SetLogy();
        h_SNR->Draw();
        c5->cd(id*2+s+1);
        TH1D *h_Length = h_COG_Length[id*2+s]->ProjectionY();
        h_Length->SetTitle(Form("L%d_S%d Cluster Length",id,s));
        h_Length->SetLabelSize(0.12);
        h_Length->SetLabelSize(0.08,"Y");
        gPad->SetLogy();
        h_Length->Draw();
        c6->cd(id*2+s+1);
        TH1D *h_Seed_SN = h_Seed_Value_SN[id*2+s]->ProjectionY();
        h_Seed_SN->SetTitle(Form("L%d_S%d Seed SN",id,s));
        h_Seed_SN->SetLabelSize(0.12);
        h_Seed_SN->SetLabelSize(0.08,"Y");
        gPad->SetLogy();
        h_Seed_SN->Draw();
        c7->cd(id*2+s+1);
        TH1D *h_Seed_ADC = h_Seed_Value_SN[id*2+s]->ProjectionX();
        h_Seed_ADC->SetTitle(Form("L%d_S%d Seed ADC",id,s));
        h_Seed_ADC->SetLabelSize(0.12);
        h_Seed_ADC->SetLabelSize(0.08,"Y");
        gPad->SetLogy();
        h_Seed_ADC->Draw();
      }
    }
  }
};

#endif

