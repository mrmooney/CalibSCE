#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TGaxis.h"

using namespace std;

Double_t Lx = 2.5;
Double_t Ly = 2.5;
Double_t Lz = 10.0;

void studyResults()
{
  //double xScale = 1.0;
  //double yScale = 1.0; 
  //double zScale = 1.0;
  double xScale = 2.56/Lx;
  double yScale = 2.33/Ly;
  double zScale = 10.37/Lz;
  
  TGaxis::SetMaxDigits(2);
  Double_t stops[5] = {0.00,0.34,0.61,0.84,1.00};
  Double_t red[5] = {0.00,0.00,0.87,1.00,0.51};
  Double_t green[5] = {0.00,0.81,1.00,0.20,0.00};
  Double_t blue[5] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);

  TH1F *DxHist = new TH1F("DxHist","",100,-5.0,5.0);
  TH1F *DxHist_calib = new TH1F("DxHist_calib","",100,-5.0,5.0);
  TH1F *DxHist_actual = new TH1F("DxHist_actual","",100,-5.0,5.0);
  TH1F *DyHist = new TH1F("DyHist","",100,-5.0,5.0);
  TH1F *DyHist_calib = new TH1F("DyHist_calib","",100,-20.0,20.0);
  TH1F *DyHist_actual = new TH1F("DyHist_actual","",100,-20.0,20.0);
  TH1F *DzHist = new TH1F("DzHist","",100,-5.0,5.0);

  TH2F *DyHist2D = new TH2F("DyHist2D","",26,-0.05,2.55,26,-0.05,2.55);
  TH2F *DyHist2D_Num = new TH2F("DyHist2D_Num","",26,-0.05,2.55,26,-0.05,2.55);
  
  TH2F *DxHist2D = new TH2F("DxHist2D","",26,-0.05,2.55,26,-0.05,2.55);
  TH2F *DxHist2D_Num = new TH2F("DxHist2D_Num","",26,-0.05,2.55,26,-0.05,2.55);

  TFile* fileActual = new TFile("data/dispOutput_MicroBooNE_E273.root");
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
  TTreeReader readerFwd("SpaCEtree_fwdDisp", fileActual);
  TTreeReaderValue<Double_t> Dx_fwd(readerFwd, "Dx.data_fwdDisp");
  TTreeReaderValue<Double_t> Dy_fwd(readerFwd, "Dy.data_fwdDisp");
  TTreeReaderValue<Double_t> Dz_fwd(readerFwd, "Dz.data_fwdDisp");
  TTreeReaderValue<Double_t> xReco_fwd(readerFwd, "x_reco.data_fwdDisp");
  TTreeReaderValue<Double_t> yReco_fwd(readerFwd, "y_reco.data_fwdDisp");
  TTreeReaderValue<Double_t> zReco_fwd(readerFwd, "z_reco.data_fwdDisp");
  TTreeReaderValue<Double_t> xTrue_fwd(readerFwd, "x_true.data_fwdDisp");
  TTreeReaderValue<Double_t> yTrue_fwd(readerFwd, "y_true.data_fwdDisp");
  TTreeReaderValue<Double_t> zTrue_fwd(readerFwd, "z_true.data_fwdDisp");
  TTreeReaderValue<Int_t> elecFate_fwd(readerFwd, "elecFate.data_fwdDisp");
  
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

  while (readerActual.Next())
  {
    readerCalib.Next();
    readerFwd.Next();

    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*zReco_actual > 4.55) && (*zReco_actual < 5.45)) {
    if((*elecFate_actual == 1) && (*elecFate_calib == 1)) {
      
    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*xReco_calib > 0.1)) {
    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*xReco_calib < 0.5)) {
    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*xReco_calib > 0.2)) {
    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*yReco_calib > 1.0) && (*yReco_calib < 1.5) && (*xReco_calib < 1.0)) {
    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*yReco_calib < 0.5) && (*xReco_calib < 1.0) && (*xReco_calib > 0.1)) {
    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*yReco_calib > 2.0) && (*xReco_calib < 1.0) && (*xReco_calib > 0.1)) {
    //if((*elecFate_actual == 1) && (*elecFate_calib == 1) && (*zReco_calib > 9.5) && (*xReco_calib < 1.0) && (*xReco_calib > 0.1)) {

      
      DxHist->Fill(-100.0*(xScale*(*Dx_calib) - *Dx_actual)); // negative sign in for LArSoft coordinate system (and below for dX)
      DyHist->Fill(100.0*(yScale*(*Dy_calib) - *Dy_actual));
      DzHist->Fill(100.0*(zScale*(*Dz_calib) - *Dz_actual));
      
      DyHist2D->Fill(*xReco_calib,*yReco_calib,100.0*(yScale*(*Dy_calib) - *Dy_actual));
      DyHist2D_Num->Fill(*xReco_calib,*yReco_calib,1.0);
      
      DxHist2D->Fill(*xReco_calib,*yReco_calib,-100.0*(xScale*(*Dx_calib) - *Dx_actual));
      DxHist2D_Num->Fill(*xReco_calib,*yReco_calib,1.0);
      
      DxHist_actual->Fill(-100.0*(*Dx_actual));
      DxHist_calib->Fill(-100.0*xScale*(*Dx_calib));
      
      DyHist_actual->Fill(100.0*(*Dy_actual));
      DyHist_calib->Fill(100.0*yScale*(*Dy_calib));


      //DxHist->Fill(-100.0*(xScale*(*Dx_calib) - 0.0));
      //DyHist->Fill(100.0*(yScale*(*Dy_calib) - 0.0));
      //DzHist->Fill(100.0*(zScale*(*Dz_calib) - 0.0));
      //
      //DyHist2D->Fill(*xReco_calib,*yReco_calib,100.0*(yScale*(*Dy_calib) - 0.0));
      //DyHist2D_Num->Fill(*xReco_calib,*yReco_calib,1.0);
      //
      //DxHist2D->Fill(*xReco_calib,*yReco_calib,-100.0*(xScale*(*Dx_calib) - 0.0));
      //DxHist2D_Num->Fill(*xReco_calib,*yReco_calib,1.0);
      //
      //DxHist_actual->Fill(100.0*(0.0));
      //DxHist_calib->Fill(100.0*xScale*(*Dx_calib));
      //
      //DyHist_actual->Fill(100.0*(0.0));
      //DyHist_calib->Fill(100.0*yScale*(*Dy_calib));
    }
  }
  DxHist->Scale(1.0/DxHist->Integral());
  DyHist->Scale(1.0/DyHist->Integral());
  DzHist->Scale(1.0/DzHist->Integral());
  DxHist_calib->Scale(1.0/DxHist_calib->Integral());
  DxHist_actual->Scale(1.0/DxHist_actual->Integral());
  DyHist_calib->Scale(1.0/DyHist_calib->Integral());
  DyHist_actual->Scale(1.0/DyHist_actual->Integral());

  for (int i = 1; i <= DyHist2D->GetNbinsX(); i++) {
    for (int j = 1; j <= DyHist2D->GetNbinsY(); j++) {
      if(DyHist2D_Num->GetBinContent(i,j) > 0.0) {
        DyHist2D->SetBinContent(i,j,DyHist2D->GetBinContent(i,j)/DyHist2D_Num->GetBinContent(i,j));
      }
      if(DxHist2D_Num->GetBinContent(i,j) > 0.0) {
        DxHist2D->SetBinContent(i,j,DxHist2D->GetBinContent(i,j)/DxHist2D_Num->GetBinContent(i,j));
      }
    }
  }

  cout << endl << DxHist->GetRMS() << " " << DxHist->GetMean() << endl;
  cout << endl << DyHist->GetRMS() << " " << DyHist->GetMean() << endl;
  cout << endl << DzHist->GetRMS() << " " << DzHist->GetMean() << endl;
  
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

  DxHist_calib->SetLineWidth(2.0);
  DxHist_calib->SetLineColor(kBlue);

  DxHist_actual->SetLineWidth(2.0);
  DxHist_actual->SetLineColor(kGreen+2);

  DxHist->Draw("HIST");
  DxHist_actual->Draw("HISTsame");
  DxHist_calib->Draw("HISTsame");
  DxHist->Draw("AXISsame");
  
  DxHist->SetMinimum(0.0001);
  c_Dx->SaveAs("studyResults_DxHist.png");

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

  DyHist_calib->SetLineWidth(2.0);
  DyHist_calib->SetLineColor(kRed);

  DyHist_actual->SetLineWidth(2.0);
  DyHist_actual->SetLineColor(kGreen+2);

  DyHist->Draw("HIST");
  DyHist_actual->Draw("HISTsame");
  DyHist_calib->Draw("HISTsame");
  DyHist->Draw("AXISsame");

  DyHist->SetMinimum(0.0001);
  c_Dy->SaveAs("studyResults_DyHist.png");

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
  c_Dz->SaveAs("studyResults_DzHist.png");

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
  c_Dxyz->SaveAs("studyResults_DxyzHist.png");

  TCanvas *c_Dx2D = new TCanvas();
  c_Dx2D->cd();
  DxHist2D->GetXaxis()->SetRangeUser(0.0,2.5);
  DxHist2D->GetYaxis()->SetRangeUser(0.0,2.5);
  DxHist2D->Draw("COLZ");
  c_Dx2D->SaveAs("studyResults_DxHist2D.png");

  TCanvas *c_Dy2D = new TCanvas();
  c_Dy2D->cd();
  DyHist2D->GetXaxis()->SetRangeUser(0.0,2.5);
  DyHist2D->GetYaxis()->SetRangeUser(0.0,2.5);
  DyHist2D->Draw("COLZ");
  c_Dy2D->SaveAs("studyResults_DyHist2D.png");

  return 0;
}
