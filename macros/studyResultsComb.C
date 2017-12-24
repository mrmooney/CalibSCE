#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TTreeReader.h"

const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;
//const Double_t Lx = 3.6
//const Double_t Ly = 6.0;
//const Double_t Lz = 7.2;

const Int_t numDivisions = 25;

Double_t comb_Dx[26][26][101];
Double_t comb_Dy[26][26][101];
Double_t comb_Dz[26][26][101];
Double_t comb_N[26][26][101];

void studyResultsComb()
{
  TGaxis::SetMaxDigits(2);

  Double_t numDivisions_x = numDivisions;
  Double_t numDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)numDivisions));
  Double_t numDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)numDivisions));

  Double_t stops[5] = {0.00,0.34,0.61,0.84,1.00};
  Double_t red[5] = {0.00,0.00,0.87,1.00,0.51};
  Double_t green[5] = {0.00,0.81,1.00,0.20,0.00};
  Double_t blue[5] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);

  TH1F *DxHist = new TH1F("DxHist","",100,-5.0,5.0);
  TH1F *DyHist = new TH1F("DyHist","",100,-5.0,5.0);
  TH1F *DzHist = new TH1F("DzHist","",100,-5.0,5.0);
  
  TFile* fileActual = new TFile("data/dispOutput_MicroBooNE_E273_7Sets.root");
  TTreeReader readerActual("SpaCEtree_bkwdDisp", fileActual);
  TTreeReaderValue<Double_t> Dx_actual(readerActual, "Dx.data_bkwdDisp");
  TTreeReaderValue<Double_t> Dy_actual(readerActual, "Dy.data_bkwdDisp");
  TTreeReaderValue<Double_t> Dz_actual(readerActual, "Dz.data_bkwdDisp");
  TTreeReaderValue<Double_t> xReco_actual(readerActual, "x_reco.data_bkwdDisp");
  TTreeReaderValue<Double_t> yReco_actual(readerActual, "y_reco.data_bkwdDisp");
  TTreeReaderValue<Double_t> zReco_actual(readerActual, "z_reco.data_bkwdDisp");
  TTreeReaderValue<Double_t> xTrue_actual(readerActual, "x_true.data_bkwdDisp");
  TTreeReaderValue<Double_t> yTrue_actual(readerActual, "y_true.data_bkwdDisp");
  TTreeReaderValue<Double_t> zTrue_actual(readerActual, "z_true.data_bkwdDisp");
  TTreeReaderValue<Int_t> elecFate_actual(readerActual, "elecFate.data_bkwdDisp");
  
  TFile* fileCalib = new TFile("output.root");
  TTreeReader readerCalib("SpaCEtree_calib", fileCalib);
  TTreeReaderValue<Double_t> Dx_calib(readerCalib, "Dx.data_calib");
  TTreeReaderValue<Double_t> Dy_calib(readerCalib, "Dy.data_calib");
  TTreeReaderValue<Double_t> Dz_calib(readerCalib, "Dz.data_calib");
  TTreeReaderValue<Double_t> xReco_calib(readerCalib, "x_reco.data_calib");
  TTreeReaderValue<Double_t> yReco_calib(readerCalib, "y_reco.data_calib");
  TTreeReaderValue<Double_t> zReco_calib(readerCalib, "z_reco.data_calib");
  TTreeReaderValue<Double_t> xTrue_calib(readerCalib, "x_true.data_calib");
  TTreeReaderValue<Double_t> yTrue_calib(readerCalib, "y_true.data_calib");
  TTreeReaderValue<Double_t> zTrue_calib(readerCalib, "z_true.data_calib");
  TTreeReaderValue<Int_t> elecFate_calib(readerCalib, "elecFate.data_calib");

  for (int x = 0; x <= numDivisions_x; x++ ) {
    for (int y = 0; y <= numDivisions_y; y++ ) {
      for (int z = 0; z <= numDivisions_z; z++ ) {
        comb_Dx[x][y][z] = 0.0;
        comb_Dy[x][y][z] = 0.0;
        comb_Dz[x][y][z] = 0.0;
	comb_N[x][y][z] = 0.0;
      }
    }
  }
  
  while (readerActual.Next())
  {
    readerCalib.Next();

    if((*elecFate_actual == 1) && (*elecFate_calib == 1)) {
    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*xReco_calib > 0.5)) {
    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*xReco_calib > 2.0)) {
      comb_Dx[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += 100.0*(*Dx_calib - *Dx_actual);
      comb_Dy[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += 100.0*(*Dy_calib - *Dy_actual);
      comb_Dz[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += 100.0*(*Dz_calib - *Dz_actual);
      comb_N[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += 1.0;
    }
  }

  for (int x = 0; x <= numDivisions_x; x++ ) {
    for (int y = 0; y <= numDivisions_y; y++ ) {
      for (int z = 0; z <= numDivisions_z; z++ ) {
	if (comb_N[x][y][z] > 0.0) {
	  DxHist->Fill(comb_Dx[x][y][z]/comb_N[x][y][z]);
          DyHist->Fill(comb_Dy[x][y][z]/comb_N[x][y][z]);
          DzHist->Fill(comb_Dz[x][y][z]/comb_N[x][y][z]);	
	}
      }
    }
  }
  DxHist->Scale(1.0/DxHist->Integral());
  DyHist->Scale(1.0/DyHist->Integral());
  DzHist->Scale(1.0/DzHist->Integral());

  gStyle->SetTitleW(0.9);
  gStyle->SetOptStat(0);

  TCanvas *c_Dx = new TCanvas();
  c_Dx->cd();
  DxHist->GetXaxis()->SetTitle("#Deltax_{calibrated} #minus #Deltax_{actual} [cm]");
  DxHist->GetXaxis()->SetTitleOffset(0.95);
  DxHist->GetXaxis()->SetTitleSize(0.045);
  DxHist->GetYaxis()->SetTitle("Arb. Units");
  DxHist->GetYaxis()->SetTitleOffset(0.95);
  DxHist->GetYaxis()->SetTitleSize(0.05);
  DxHist->GetYaxis()->SetNoExponent(kTRUE);
  DxHist->SetStats(0);
  DxHist->SetLineWidth(2.0);
  DxHist->SetLineColor(kRed);
  DxHist->Draw("HIST");
  DxHist->Draw("AXISsame");
  DxHist->SetMinimum(0.0001);
  c_Dx->SaveAs("DxHist.png");

  TCanvas *c_Dy = new TCanvas();
  c_Dy->cd();
  DyHist->GetXaxis()->SetTitle("#Deltay_{calibrated} #minus #Deltay_{actual} [cm]");
  DyHist->GetXaxis()->SetTitleOffset(0.95);
  DyHist->GetXaxis()->SetTitleSize(0.045);
  DyHist->GetYaxis()->SetTitle("Arb. Units");
  DyHist->GetYaxis()->SetTitleOffset(0.95);
  DyHist->GetYaxis()->SetTitleSize(0.05);
  DyHist->GetYaxis()->SetNoExponent(kTRUE);
  DyHist->SetStats(0);
  DyHist->SetLineWidth(2.0);
  DyHist->SetLineColor(kBlue);
  DyHist->Draw("HIST");
  DyHist->Draw("AXISsame");
  DyHist->SetMinimum(0.0001);
  c_Dy->SaveAs("DyHist.png");

  TCanvas *c_Dz = new TCanvas();
  c_Dz->cd();
  DzHist->GetXaxis()->SetTitle("#Deltaz_{calibrated} #minus #Deltaz_{actual} [cm]");
  DzHist->GetXaxis()->SetTitleOffset(0.95);
  DzHist->GetXaxis()->SetTitleSize(0.045);
  DzHist->GetYaxis()->SetTitle("Arb. Units");
  DzHist->GetYaxis()->SetTitleOffset(0.95);
  DzHist->GetYaxis()->SetTitleSize(0.05);
  DzHist->GetYaxis()->SetNoExponent(kTRUE);
  DzHist->SetStats(0);
  DzHist->SetLineWidth(2.0);
  DzHist->SetLineColor(kGreen+2);
  DzHist->Draw("HIST");
  DzHist->Draw("AXISsame");
  DzHist->SetMinimum(0.0001);
  c_Dz->SaveAs("DzHist.png");

  TCanvas *c_Dxyz = new TCanvas();
  c_Dxyz->cd();
  TLegend *leg_Dxyz = new TLegend(0.73,0.65,0.88,0.85);
  leg_Dxyz->SetLineColor(kWhite);
  leg_Dxyz->AddEntry(DxHist,"#Deltax","L");
  leg_Dxyz->AddEntry(DyHist,"#Deltay","L");
  leg_Dxyz->AddEntry(DzHist,"#Deltaz","L");
  DxHist->GetXaxis()->SetTitle("#Delta_{calibrated} #minus #Delta_{actual} [cm]");
  DxHist->Draw("HIST");
  DyHist->Draw("HISTsame");
  DzHist->Draw("HISTsame");
  leg_Dxyz->Draw("same");
  DxHist->Draw("AXISsame");
  DxHist->SetMaximum(1.1*max(max(DyHist->GetMaximum(),DzHist->GetMaximum()),DxHist->GetMaximum()));
  DxHist->SetMinimum(0.0001);
  c_Dxyz->SaveAs("DxyzHist.png");
}
