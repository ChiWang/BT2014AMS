/*
 *  $Id: DataQuality.C, 2015-01-01 16:59:38 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 30/11/2014
*/

/*
 * refer to Readme.md
 *
 * // unit for length: cm 
 * // strip 0,0 as reference.
 *
 */

#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "TFile.h"
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

#ifndef AMS_Alignment_C
#define AMS_Alignment_C

namespace AMS{
using namespace Common;
namespace Alignment{
// *
// *  TODO: ios has not be declared
// *
  TString outFilename = "./Calibration/AMS/"+Conf::File+"_align.txt";
  //ofstream output(outFilename,ios::out | iso::app);
  //output<<"ladder id\tmean\terror";

  void SingleStrack_S_Side(bool reload=true){    // for each ladder
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(100111);

    TH1F *h_align_CoG[NLadder][3] = {{0}};  // cluster numbers, 0: s-side, 1: k-side-sensor0, 2: k-side-sensor-1
    h_align_CoG[0][0] = new TH1F("L0_S0 Reference Position","L0_S0 Reference Position",1500,0,10);
    h_align_CoG[0][1] = new TH1F("L0_S1_0 Reference Position","L0_S1 Reference Position",1500,0,10);
    h_align_CoG[0][2] = new TH1F("L0_S1_1 Reference Position","L0_S1 Reference Position",1500,0,10);
    for(short i =1;i<NLadder; ++i){
      h_align_CoG[i][0] = new TH1F(Form("L%d_S0 Offset",i),Form("L%d_S0 Offset",i),1000,-2.5,1.5);
      h_align_CoG[i][1] = new TH1F(Form("L%d_S1_0 Offset",i),Form("L%d_S0 Offset",i),1000,-2.5,1.5);
      h_align_CoG[i][2] = new TH1F(Form("L%d_S1_1 Offset",i),Form("L%d_S0 Offset",i),1000,-2.5,1.5);
    }

    TH2F *h_offset_SeedAdd[NLadder][3] = {{0}};  // cluster numbers, 0: s-side, 1: k-side-sensor0, 2: k-side-sensor-1
    for(short i =0;i<NLadder-1; ++i){
      h_offset_SeedAdd[i][0] = new TH2F(Form("L%d_S0 Offset_SeedAdd",i+1),  Form("L%d_S0 Offset_SeedAdd",i+1),640,0,640,1000,-2.5,1.5);
      h_offset_SeedAdd[i][1] = new TH2F(Form("L%d_S1_0 Offset_SeedAdd",i+1),Form("L%d_S1 Offset_SeedAdd",i+1),384,640,1024,1000,-2.5,1.5);
      h_offset_SeedAdd[i][2] = new TH2F(Form("L%d_S1_1 Offset_SeedAdd",i+1),Form("L%d_S1 Offset_SeedAdd",i+1),384,640,1024,1000,-2.5,1.5);
    }

    for(Conf::evtID =0;Conf::evtID<Conf::entries;++Conf::evtID){
      LoadEvent();
      // one track event
      if(! N_ClustersInLadder_I(1,0)){   // both sides
        continue;
      }
      if(! ClusterNumberLessThan2_forAllS_Side()){
        continue;
      }
      // update reference
      float Posi_Ref_ladder0[2]={0.0,0.};   // ladder 0, side 0, 1
      short k_Ref_SensorID = 0;       // for k-side, short(long) ladder has 2(6) group, one group contains 2 silicon sensors. While alignmenting, we'd use vertical tracks (means, clusters of ladder 1~4 must in the same sensor(0 or 1) as the cluster which belongs to ladder 0)
      int n_cls = Conf::AMS_Evt->Cls->GetEntriesFast();
      for(short ic=0;ic<n_cls;++ic){
        Cluster *aCluster = Conf::AMS_Evt->GetCluster(ic);
        short order = Setup::LadderInOrder[Conf::ExHall][aCluster->ladder];
        if(order != 0){
          continue;
        }
        short side = aCluster->side;
        Posi_Ref_ladder0[side] = GetPosition(aCluster,false);
        if(side == 1 && aCluster->GetSeedAdd()>(640+192)){  // 640: number of s-side readout strips(one sensor), 192: number of k-side readout strips (one sensor)
          k_Ref_SensorID = 1;
        }
        h_align_CoG[0][side+(side ==1 ? k_Ref_SensorID : 0)]->Fill(Posi_Ref_ladder0[side]);
      }
      // alignment
      for(short ic=0;ic<n_cls;++ic){
        Cluster *aCluster = Conf::AMS_Evt->GetCluster(ic);
        short order = Setup::LadderInOrder[Conf::ExHall][aCluster->ladder];
        if(order == 0){
          continue;
        }
        short side = aCluster->side;
        float offV = GetPosition(aCluster,false) - Posi_Ref_ladder0[side];
        if(side == 1){
          short k_SensorID = (aCluster->GetSeedAdd()<(640+192)) ? 0 : 1;
          if(k_SensorID == k_Ref_SensorID){
            h_align_CoG[order][side + k_SensorID]->Fill(offV);
            h_offset_SeedAdd[order-1][side+k_SensorID]->Fill(aCluster->GetSeedAdd(),offV);
          }
        }else{
          h_align_CoG[order][side]->Fill(offV);
          h_offset_SeedAdd[order-1][side]->Fill(aCluster->GetSeedAdd(),offV);
        }
      }
    }

    TCanvas *c0 = new TCanvas(Conf::File+"  Alignment Ref.",Conf::File+"  Alignment Ref.");
    c0->Divide(1,2);
    h_align_CoG[0][0]->SetXTitle("X / cm");
    h_align_CoG[0][1]->SetXTitle("Y / cm");
    for(short s=0;s<2;++s){
      c0->cd(s+1);
      gPad->SetLogy();
      h_align_CoG[0][s]->SetLabelSize(0.06);
      h_align_CoG[0][s]->SetLabelSize(0.04,"Y");
      h_align_CoG[0][s]->SetTitleSize(0.05,"X");
      h_align_CoG[0][s]->Draw();
      float mean = h_align_CoG[0][s]->GetMean(), rms = h_align_CoG[0][s]->GetRMS();
      Conf::gausFit->SetRange(mean-rms,mean+rms);
      h_align_CoG[0][s]->Fit(Conf::gausFit,"R0Q");
      Conf::gausFit->DrawCopy("lsame");
      if(s ==1){
        h_align_CoG[0][2]->SetLineColor(3);
        h_align_CoG[0][2]->Draw("same");
      }
    }

    TCanvas *c1 = new TCanvas(Conf::File+"  Alignment",Conf::File+"  Alignment");
    c1->Divide(2,4,0.,0.0);
    for(short id=1;id<NLadder;++id){
      for(short s=0;s<2;++s){
        c1->cd((id-1)*2+s+1);
        gPad->SetLogy();
        if(s==0){
          h_align_CoG[id][s]->SetXTitle("X / cm");
        }else{
          h_align_CoG[id][s]->SetXTitle("Y / cm");
        }
        h_align_CoG[id][s]->SetLabelSize(0.12);
        h_align_CoG[id][s]->SetLabelSize(0.08,"Y");
        h_align_CoG[id][s]->SetTitleSize(0.04,"X");
        h_align_CoG[id][s]->Draw();
        float mean = h_align_CoG[id][s]->GetMean(), rms = h_align_CoG[id][s]->GetRMS();
        Conf::gausFit->SetRange(mean-rms,mean+rms);
        h_align_CoG[id][s]->Fit(Conf::gausFit,"R0Q");
        Conf::gausFit->DrawCopy("lsame");
        //output<<
      }
      c1->cd((id-1)*2+2);
      h_align_CoG[id][2]->SetLineColor(3);
      h_align_CoG[id][2]->Draw("same");
    }
    TCanvas *c2 = new TCanvas(Conf::File+"  Offset_SeedAdd",Conf::File+"  Offset_SeedAdd");
    c2->Divide(2,4,0.,0.0);
    for(short id=0;id<NLadder-1;++id){
      for(short s=0;s<2;++s){
        c2->cd(id*2+s+1);
        gStyle->SetOptStat(11111111);
        gStyle->SetOptFit(111111111);
        h_offset_SeedAdd[id][s]->SetXTitle("Seed ID");
        h_offset_SeedAdd[id][s]->SetYTitle("Offset / cm");
        h_offset_SeedAdd[id][s]->SetLabelSize(0.12);
        h_offset_SeedAdd[id][s]->SetLabelSize(0.08,"Y");
        //h_offset_SeedAdd[id][s]->SetTitleSize(0.04,"Y");
        h_offset_SeedAdd[id][s]->Draw("colz");
        //h_offset_SeedAdd[id][s]->ProfileX()->Fit(Conf::linearFit,"0Q");//QF
        //Conf::linearFit->DrawCopy("same");
        //output<<
      }
      c2->cd(id*2+2);
      h_offset_SeedAdd[id][2]->Draw("same");
    }
    if(reload) Conf::UpdateOffset(outFilename); 
  }
}
};

#endif


