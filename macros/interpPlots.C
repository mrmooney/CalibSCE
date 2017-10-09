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

void makeInterpPlot(Int_t numDivisions, Double_t zBin, Int_t logScale, Int_t whichPlot, Int_t fieldComp, Double_t Efield);

void interpPlots(Int_t numDivisions, Double_t zBin, Int_t logScale, Int_t whichPlot, Int_t fieldComp, Double_t Efield = 500.0)
{
  makeInterpPlot(numDivisions,zBin,logScale,whichPlot,fieldComp,Efield);
}

void makeInterpPlot(Int_t numDivisions, Double_t zBin, Int_t logScale, Int_t whichPlot, Int_t fieldComp, Double_t Efield)
{
  const Char_t *interpExt = "_5OuterLayers";

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

  Char_t *fieldCompString = (Char_t*) "";
  Char_t *fieldCompString2 = (Char_t*) "";
  if(fieldComp == 1)
  {
    fieldCompString = (Char_t*) "Ex";
    fieldCompString2 = (Char_t*) "(^{}E_{x} #minus E_{0})";
  }
  else if(fieldComp == 2)
  {
    fieldCompString = (Char_t*) "Ey";
    fieldCompString2 = (Char_t*) "E_{y}";
  }
  else if(fieldComp == 3)
  {
    fieldCompString = (Char_t*) "Ez";
    fieldCompString2 = (Char_t*) "E_{z}";
  }

  TCanvas *c = new TCanvas(Form("c_%d_%d_%d_%d_%d_%d",numDivisions,(Int_t)zBin,logScale,whichPlot,fieldComp,(Int_t)Efield),"",600,600);
  c->cd();

  TGaxis::SetMaxDigits(2);

  Double_t stops[5] = {0.00,0.34,0.61,0.84,1.00};
  Double_t red[5] = {0.00,0.00,0.87,1.00,0.51};
  Double_t green[5] = {0.00,0.81,1.00,0.20,0.00};
  Double_t blue[5] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);

  TFile* fileSCE = new TFile(Form("dispOutput_MicroBooNE_E%d.root",(Int_t)Efield));
  TTree* treeSCE = (TTree*)fileSCE->Get("SpaCEtree");
  TProfile2D* profileSCE = new TProfile2D("profileSCE","",numDivisions_x+1,-Lx/(2.0*((Double_t) numDivisions_x)),Lx+Lx/(2.0*((Double_t) numDivisions_x)),numDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDivisions_y)));
  //TProfile2D* profileSCE = new TProfile2D("profileSCE","",numDivisions_x+1,-1.0*Lx/(2.0*((Double_t) numDivisions_x)),Lx+Lx/(2.0*((Double_t) numDivisions_x)),numDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDivisions_y)));

  treeSCE->Draw(Form("(-1./%f)*%s:(ypoint-1.25):(2.5-xpoint)>>profileSCE",Efield,fieldCompString),Form("zpoint < %f && zpoint > %f",zpointVal+0.025,zpointVal-0.025),"profCOLZ");

  profileSCE->SetStats(0);  
  
  //gPad->Clear();
  gStyle->SetTitleW(0.9);

  TProfile2D* profilePlot = new TProfile2D;
  
  profilePlot = (TProfile2D*)profileSCE->Clone("profilePlot");
  //profilePlot->SetTitle(Form("Simulated %s / ^{}E_{0} [%%]:   Z = %.2f m",fieldCompString2,zpointVal));
  profilePlot->SetTitle(Form("Simulated %s / ^{}E_{0} [%%]:   Z = %.2f m",fieldCompString2,5.18)); // for public note (central Z slice only)
  profilePlot->GetXaxis()->SetTitle("X [m]");
  profilePlot->GetXaxis()->SetTitleOffset(1.0);
  profilePlot->GetXaxis()->SetTitleSize(0.04);
  profilePlot->GetYaxis()->SetTitle("Y [m]");
  profilePlot->GetYaxis()->SetTitleOffset(1.2);
  profilePlot->GetYaxis()->SetTitleSize(0.04);
  profilePlot->SetStats(0);
  profilePlot->GetZaxis()->SetNoExponent(kTRUE);
  profilePlot->GetZaxis()->SetLabelSize(0.025);
  profilePlot->GetYaxis()->SetRange(2,profilePlot->GetNbinsY()-1); // for public note
  //profilePlot->GetZaxis()->SetRangeUser(-1.0,1.0); // for public note (central Z slice, dZ only)
  profilePlot->Draw("COLZ");

  gPad->SetLogz(logScale);

  c->SaveAs(Form("Efield_%d_%d_%d_%d_%d_%d.png",numDivisions,(Int_t)zBin,logScale,whichPlot,fieldComp,(Int_t)Efield));
  c->SaveAs(Form("Efield_%d_%d_%d_%d_%d_%d.pdf",numDivisions,(Int_t)zBin,logScale,whichPlot,fieldComp,(Int_t)Efield));
}
