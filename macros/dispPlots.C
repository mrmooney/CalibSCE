#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TMath.h"

const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;
//const Double_t Lx = 2.25;
//const Double_t Ly = 2.0;
//const Double_t Lz = 1.6;
//const Double_t Lx = 2.4; // was 3.6
//const Double_t Ly = 6.0;
//const Double_t Lz = 7.2;
//const Double_t Lx = 2.0;
//const Double_t Ly = 4.0;
//const Double_t Lz = 5.0;
//const Double_t Lx = 6.0;
//const Double_t Ly = 6.0;
//const Double_t Lz = 6.0;

const Int_t numIter = 100;

void makeDispPlot(Int_t numDivisions, Double_t zBin, Int_t logScale, Int_t whichPlot, Int_t dispComp, Double_t Efield);

void dispPlots(Int_t numDivisions, Double_t zBin, Int_t logScale, Int_t whichPlot, Int_t dispComp, Double_t Efield = 500.0)
{
  makeDispPlot(numDivisions,zBin,logScale,whichPlot,dispComp,Efield);
}

void makeDispPlot(Int_t numDivisions, Double_t zBin, Int_t logScale, Int_t whichPlot, Int_t dispComp, Double_t Efield)
{
  Double_t numDivisions_x = numDivisions;
  Double_t numDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)numDivisions));
  Double_t numDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)numDivisions));

  Double_t logFactor = 1.0;
  if(logScale == 2)
  {
    logScale = 1;
    logFactor = -1.0;
  }

  Double_t zpointVal = (((Double_t) zBin)/((Double_t) numDivisions_z))*Lz;

  Char_t *dispCompString = (Char_t*) "";
  Char_t *dispCompString2 = (Char_t*) "";
  
  if(whichPlot == 1)
  {
    if(dispComp == 1)
    {
      dispCompString = (Char_t*) "Dx";
      dispCompString2 = (Char_t*) "X_{reco} #minus X_{true}";
    }
    else if(dispComp == 2)
    {
      dispCompString = (Char_t*) "Dy";
      dispCompString2 = (Char_t*) "Y_{reco} #minus Y_{true}";
    }
    else if(dispComp == 3)
    {
      dispCompString = (Char_t*) "Dz";
      dispCompString2 = (Char_t*) "Z_{reco} #minus Z_{true}";
    }
  }
  else if(whichPlot == 2)
  {
    if(dispComp == 1)
    {
      dispCompString = (Char_t*) "Dx";
      dispCompString2 = (Char_t*) "X_{true} #minus X_{reco}";
    }
    else if(dispComp == 2)
    {
      dispCompString = (Char_t*) "Dy";
      dispCompString2 = (Char_t*) "Y_{true} #minus Y_{reco}";
    }
    else if(dispComp == 3)
    {
      dispCompString = (Char_t*) "Dz";
      dispCompString2 = (Char_t*) "Z_{true} #minus Z_{reco}";
    }
  }

  TCanvas *c = new TCanvas(Form("c_%d_%d_%d_%d_%d_%d",numDivisions,(Int_t)zBin,logScale,whichPlot,dispComp,(Int_t)Efield),"",600,600);
  c->cd();

  TGaxis::SetMaxDigits(2);

  Double_t stops[5] = {0.00,0.34,0.61,0.84,1.00};
  Double_t red[5] = {0.00,0.00,0.87,1.00,0.51};
  Double_t green[5] = {0.00,0.81,1.00,0.20,0.00};
  Double_t blue[5] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);

  TFile* fileSCE = new TFile(Form("dispOutput_MicroBooNE_E%d.root",(Int_t)Efield));
  TTree* treeSCE;
  if(whichPlot == 1)
    treeSCE = (TTree*)fileSCE->Get("SpaCEtree_fwdDisp");
    //treeSCE = (TTree*)fileSCE->Get("SpaCEtree_disp");
  else if(whichPlot == 2)
    treeSCE = (TTree*)fileSCE->Get("SpaCEtree_bkwdDisp");

  TProfile2D* profileSCE = new TProfile2D("profileSCE","",numDivisions_x+1,-Lx/(2.0*((Double_t) numDivisions_x)),Lx+Lx/(2.0*((Double_t) numDivisions_x)),numDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDivisions_y)));

  if(whichPlot == 1)
  {
    if(dispComp == 1) {
      treeSCE->Draw(Form("-100.0*%s:(y_true-1.25):(2.5-x_true)>>profileSCE",dispCompString),Form("z_true < %f && z_true > %f",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
    }
    else {
      treeSCE->Draw(Form("100.0*%s:(y_true-1.25):(2.5-x_true)>>profileSCE",dispCompString),Form("z_true < %f && z_true > %f",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
    }
  }
  else if(whichPlot == 2)
  {
    if(dispComp == 1) {
      treeSCE->Draw(Form("-100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileSCE",dispCompString),Form("z_reco < %f && z_reco > %f && elecFate == 1",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
    }
    else {
      treeSCE->Draw(Form("100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileSCE",dispCompString),Form("z_reco < %f && z_reco > %f && elecFate == 1",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
    }
  }

  profileSCE->SetStats(0);

  //gPad->Clear();
  gStyle->SetTitleW(0.9);

  TProfile2D* profilePlot = new TProfile2D;
  
  profilePlot = (TProfile2D*)profileSCE->Clone("profilePlot");

  //profilePlot->SetTitle(Form("%s [cm]:  Z = %.2f m",dispCompString2,zpointVal));
  profilePlot->SetTitle(Form("%s [cm]:  Z = %.2f m",dispCompString2,5.18)); // for public note (central Z slice only)
  profilePlot->GetXaxis()->SetTitle("X [m]");
  profilePlot->GetXaxis()->SetTitleOffset(1.0);
  profilePlot->GetXaxis()->SetTitleSize(0.04);
  profilePlot->GetYaxis()->SetTitle("Y [m]");
  profilePlot->GetYaxis()->SetTitleOffset(1.2);
  profilePlot->GetYaxis()->SetTitleSize(0.04);
  //((TGaxis*) profilePlot->GetZaxis())->SetMaxDigits(1);
  profilePlot->GetZaxis()->SetNoExponent(kTRUE);
  profilePlot->GetZaxis()->SetLabelSize(0.025);
  profilePlot->GetYaxis()->SetRange(2,profilePlot->GetNbinsY()-1); // for public note
  //profilePlot->GetZaxis()->SetRangeUser(-1.0,1.0); // for public note (central Z slice, dZ only)
  profilePlot->SetStats(0);
  profilePlot->Draw("COLZ");

  gPad->SetLogz(logScale);

  c->SaveAs(Form("spatial_%d_%d_%d_%d_%d_%d.png",numDivisions,(Int_t)zBin,logScale,whichPlot,dispComp,(Int_t)Efield));
  c->SaveAs(Form("spatial_%d_%d_%d_%d_%d_%d.pdf",numDivisions,(Int_t)zBin,logScale,whichPlot,dispComp,(Int_t)Efield));
}
