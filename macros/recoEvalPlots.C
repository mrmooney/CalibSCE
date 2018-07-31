#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"

const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;
//const Double_t Lx = 2.25;
//const Double_t Ly = 2.0;
//const Double_t Lz = 1.6;

const Int_t numIter = 100;

void makeRecoEvalPlot(Int_t numDispDivisions, Double_t zBin, Int_t logScale, Int_t whichPlot, Int_t dispComp);

void recoEvalPlots(Int_t numDispDivisions, Double_t zBin, Int_t logScale, Int_t whichPlot, Int_t dispComp)
{
  makeRecoEvalPlot(numDispDivisions,zBin,logScale,whichPlot,dispComp);
}

void makeRecoEvalPlot(Int_t numDispDivisions, Double_t zBin, Int_t logScale, Int_t whichPlot, Int_t dispComp)
{
  Double_t numDispDivisions_x = numDispDivisions;
  Double_t numDispDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)numDispDivisions));
  Double_t numDispDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)numDispDivisions));

  Double_t logFactor = 1.0;
  if(logScale == 2)
  {
    logScale = 1;
    logFactor = -1.0;
  }

  Double_t zpointVal = (((Double_t) zBin)/((Double_t) numDispDivisions_z))*Lz;

  Char_t *dispCompString = (Char_t*) "";
  Char_t *dispCompString2 = (Char_t*) "";
  Char_t *dispCompString3 = (Char_t*) "";
  if(dispComp == 1)
  {
    dispCompString = (Char_t*) "Dx";
    dispCompString2 = (Char_t*) "{}^{}X_{true} - X_{reco}";
    dispCompString3 = (Char_t*) "X";
  }
  else if(dispComp == 2)
  {
    dispCompString = (Char_t*) "Dy";
    dispCompString2 = (Char_t*) "{}^{}Y_{true} - Y_{reco}";
    dispCompString3 = (Char_t*) "Y";
  }
  else if(dispComp == 3)
  {
    dispCompString = (Char_t*) "Dz";
    dispCompString2 = (Char_t*) "{}^{}Z_{true} - Z_{reco}";
    dispCompString3 = (Char_t*) "Z";
  }

  TCanvas *c = new TCanvas(Form("c_%d_%d_%d_%d_%d",numDispDivisions,(Int_t)zBin,logScale,whichPlot,dispComp),"",600,600);
  c->cd();

  TGaxis::SetMaxDigits(2);

  Double_t stops[5] = {0.00,0.34,0.61,0.84,1.00};
  Double_t red[5] = {0.00,0.00,0.87,1.00,0.51};
  Double_t green[5] = {0.00,0.81,1.00,0.20,0.00};
  Double_t blue[5] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);

  TFile* fileActual = new TFile("data/dispOutput_MicroBooNE_E273.root");
  TTree* treeActual = (TTree*)fileActual->Get("SpaCEtree_bkwdDisp");
  //TProfile2D* profileActual = new TProfile2D("profileActual","",numDispDivisions_x+1,-1.0*Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),Ly+Ly/(2.0*((Double_t) numDispDivisions_y)));
  TProfile2D* profileActual = new TProfile2D("profileActual","",numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  //treeActual->Draw(Form("100.0*%s:y_reco:x_reco>>profileActual",dispCompString),Form("z_reco < %f && z_reco > %f && elecFate == 1",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
  if(dispComp == 1) {
    treeActual->Draw(Form("-100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileActual",dispCompString),Form("z_reco < %f && z_reco > %f && elecFate == 1",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
  }
  else {
    treeActual->Draw(Form("100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileActual",dispCompString),Form("z_reco < %f && z_reco > %f && elecFate == 1",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
  }
  profileActual->SetStats(0);
  
  //TFile* fileCalib = new TFile("results/output_toyMC_isotropic_10000cosmics.root");
  //TFile* fileCalib = new TFile("results/output_toyMC_isotropic_10000cosmics_withInterp.root");
  //TFile* fileCalib = new TFile("results/output_toyMC_isotropic_10000cosmics_withMCS_withInterp.root");
  //TFile* fileCalib = new TFile("results/output_toyMC_laserscan.root");
  //TFile* fileCalib = new TFile("results/output_toyMC_laserscan_withInterp.root");
  //TFile* fileCalib = new TFile("results/output_actualMC_6p0AnodeMCS_2p5CathodeMCS.root");
  //TFile* fileCalib = new TFile("results/output_actualMC_6p0AnodeMCS_2p5CathodeMCS_withInterp.root");
  //TFile* fileCalib = new TFile("results/CollabMeetingResults_Data_20kTracks.root");
  //TFile* fileCalib = new TFile("results/CollabMeetingResults_Data_20kTracks_withInterp.root");
  //TFile* fileCalib = new TFile("results/CollabMeetingResults_MC_20kTracks.root");
  //TFile* fileCalib = new TFile("results/PostCollabMeetingResults_MC_20kTracks.root");
  //TFile* fileCalib = new TFile("results/CollabMeetingResults_MC_20kTracks_withInterp.root");
  TFile* fileCalib = new TFile("output.root");
  //TFile* fileCalib = new TFile("results/new/cosmics_MC_large.root");
  //TFile* fileCalib = new TFile("results/dispOutput_MicroBooNE_calib_2cm_3.root");
  //TFile* fileCalib = new TFile("results/dispOutput_MicroBooNE_calib_2cm_3_lasersOnly.root");
  TTree* treeCalib = (TTree*)fileCalib->Get("SpaCEtree_calib");
  //TProfile2D* profileCalib = new TProfile2D("profileCalib","",numDispDivisions_x+1,-1.0*Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),Ly+Ly/(2.0*((Double_t) numDispDivisions_y)));
  TProfile2D* profileCalib = new TProfile2D("profileCalib","",numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  //treeCalib->Draw(Form("100.0*%s:y_reco:x_reco>>profileCalib",dispCompString),Form("z_reco < %f && z_reco > %f && ((elecFate == 1) || (x_reco==%f))",zpointVal+0.025,zpointVal-0.025,Lx),"profCOLZ");
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileCalib",dispCompString),Form("z_reco < %f && z_reco > %f && ((elecFate == 1) || (x_reco==%f))",zpointVal+0.025,zpointVal-0.025,Lx),"profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileCalib",dispCompString),Form("z_reco < %f && z_reco > %f && ((elecFate == 1) || (x_reco==%f))",zpointVal+0.025,zpointVal-0.025,Lx),"profCOLZ");
  }
  profileCalib->SetStats(0);

  //TFile* fileCalib = new TFile("RecoField5.root");
  //TH3F* histCalib = (TH3F*)fileCalib->Get(Form("Reco_Field_%s",dispCompString3));
  //TProfile2D* profileCalib = new TProfile2D("profileCalib","",numDispDivisions_x+1,-1.0*Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),Ly+Ly/(2.0*((Double_t) numDispDivisions_y)));
  //Double_t recoBinsX = histCalib->GetNbinsX();
  //Double_t recoBinsY = histCalib->GetNbinsY();
  //Double_t recoBinsZ = histCalib->GetNbinsZ();

  TProfile2D* profileDiff = new TProfile2D("profileDiff","",numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  for(Int_t i = 1; i <= numDispDivisions_x+1; i++)
  {
    for(Int_t j = 1; j <= numDispDivisions_y+1; j++)
    {
      if(fabs(profileActual->GetBinContent(i,j)) == 0.0)
        profileDiff->SetBinContent(i,j,0.0);
      else
        profileDiff->SetBinContent(i,j,logFactor*(1.0*profileCalib->GetBinContent(i,j)-1.0*profileActual->GetBinContent(i,j)));

      profileDiff->SetBinEntries(profileDiff->GetBin(i,j),1);
    }
  }
  profileDiff->SetStats(0);  

  //TProfile2D* profileDiff = new TProfile2D("profileDiff","",numDispDivisions_x+1,-1.0*Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),Ly+Ly/(2.0*((Double_t) numDispDivisions_y)));
  //profileDiff->SetStats(0);  
  //for(Int_t i = 1; i <= histCalib->GetNbinsX(); i++)
  //{
  //  for(Int_t j = 1; j <= histCalib->GetNbinsY(); j++)
  //  {
  //    for(Int_t k = 1; k <= histCalib->GetNbinsZ(); k++)
  //    {
  //      if(((((Double_t) k-1.0)/(recoBinsZ-1.0))*Lz > zpointVal-0.025) && ((((Double_t) k-1.0)/(recoBinsZ-1.0))*Lz < zpointVal+0.025))
  //      {
  //        //profileCalib->Fill((((Double_t) i-1.0)/(recoBinsX-1.0))*Lx,(((Double_t) j-1.0)/(recoBinsY-1.0))*Ly,histCalib->GetBinContent(i,j,k)); // W/O Coordinate Shift
  //        //profileCalib->Fill(Lx-(((Double_t) i-1.0)/(recoBinsX-1.0))*Lx,(((Double_t) j-1.0)/(recoBinsY-1.0))*Ly,histCalib->GetBinContent(i,j,k));
  //
  //        if(dispComp == 1)          
  //          profileCalib->Fill(Lx-(((Double_t) i-1.0)/(recoBinsX-1.0))*Lx,(((Double_t) j-1.0)/(recoBinsY-1.0))*Ly,-100.0*histCalib->GetBinContent(i,j,k));
  //	  else
  //          profileCalib->Fill(Lx-(((Double_t) i-1.0)/(recoBinsX-1.0))*Lx,(((Double_t) j-1.0)/(recoBinsY-1.0))*Ly,100.0*histCalib->GetBinContent(i,j,k));
  //      }
  //    }
  //  }
  //}
  //for(Int_t i = 1; i <= numDispDivisions_x+1; i++)
  //{
  //  for(Int_t j = 1; j <= numDispDivisions_y+1; j++)
  //  {
  //    Double_t normFactor = profileCalib->GetBinEntries(profileCalib->GetBin(i,j));
  //    profileCalib->SetBinContent(i,j,profileCalib->GetBinContent(i,j)/normFactor);
  //
  //    profileDiff->SetBinContent(i,j,logFactor*(1.0*profileCalib->GetBinContent(i,j)-1.0*profileActual->GetBinContent(i,j)));
  //    profileDiff->SetBinEntries(profileDiff->GetBin(i,j),1);
  //  }
  //}
  //profileCalib->SetStats(0);
    
  TProfile2D* profileRatio = (TProfile2D*)profileDiff->Clone("profileRatio");
  profileRatio->Divide(profileActual);
  for(Int_t i = 1; i <= numDispDivisions_x+1; i++)
  {
    for(Int_t j = 1; j <= numDispDivisions_y+1; j++)
    {
      if(fabs(profileActual->GetBinContent(i,j)) < 0.001)
        profileRatio->SetBinContent(i,j,0.0);

      if(fabs(profileActual->GetBinContent(i,j)) < 0.00001)
        profileActual->SetBinContent(i,j,0.0);
    }
  }
  profileRatio->SetStats(0);  

  gPad->Clear();
  gStyle->SetTitleW(0.9);

  TProfile2D* profilePlot = new TProfile2D;
  
  if(whichPlot == 1)
    profilePlot = (TProfile2D*)profileActual->Clone("profilePlot");
  else if(whichPlot == 2)
    profilePlot = (TProfile2D*)profileCalib->Clone("profilePlot");
  else if(whichPlot == 3)
    profilePlot = (TProfile2D*)profileDiff->Clone("profilePlot");
  else if(whichPlot == 4)
    profilePlot = (TProfile2D*)profileRatio->Clone("profilePlot");

  if(whichPlot == 1)
    profilePlot->SetTitle(Form("Actual %s [cm]:  Z = %.2f m",dispCompString2,zpointVal*(5.18/5.0)));
  else if(whichPlot == 2)
    profilePlot->SetTitle(Form("Calibrated %s [cm]:  Z = %.2f m",dispCompString2,zpointVal*(5.18/5.0)));
  else if(whichPlot == 3)
    profilePlot->SetTitle(Form("Calibrated-Actual %s [cm]:  Z = %.2f m",dispCompString2,zpointVal*(5.18/5.0)));
  else if(whichPlot == 4)
    profilePlot->SetTitle(Form("(Calibrated-Actual)/Actual %s:  Z = %.2f m",dispCompString2,zpointVal*(5.18/5.0)));

  profilePlot->GetXaxis()->SetTitle("X_{reco} [m]");
  profilePlot->GetXaxis()->SetTitleOffset(1.0);
  profilePlot->GetXaxis()->SetTitleSize(0.04);
  profilePlot->GetYaxis()->SetTitle("Y_{reco} [m]");
  profilePlot->GetYaxis()->SetTitleOffset(1.2);
  profilePlot->GetYaxis()->SetTitleSize(0.04);
  profilePlot->GetZaxis()->SetNoExponent(kTRUE);
  profilePlot->GetZaxis()->SetLabelSize(0.025);
  profilePlot->SetStats(0);
  profilePlot->Draw("COLZ");
  
  gPad->SetLogz(logScale);  
}
