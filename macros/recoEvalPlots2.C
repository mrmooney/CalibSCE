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

void makeRecoEvalPlot(Int_t numDispDivisions, Int_t faceNum, Int_t whichPlot, Int_t dispComp);

void recoEvalPlots2(Int_t numDispDivisions, Int_t faceNum, Int_t whichPlot, Int_t dispComp)
{
  makeRecoEvalPlot(numDispDivisions,faceNum,whichPlot,dispComp);
}

void makeRecoEvalPlot(Int_t numDispDivisions, Int_t faceNum, Int_t whichPlot, Int_t dispComp)
{
  Double_t numDispDivisions_x = numDispDivisions;
  Double_t numDispDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)numDispDivisions));
  Double_t numDispDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)numDispDivisions));

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

  TCanvas *c;

  if((faceNum == 2) || (faceNum == 3)) {
    c = new TCanvas(Form("c_%d_%d_%d_%d",numDispDivisions,faceNum,whichPlot,dispComp),"",600,600);
  }
  else {
    c = new TCanvas(Form("c_%d_%d_%d_%d",numDispDivisions,faceNum,whichPlot,dispComp),"",2400,600);
  }
  c->cd();

  TGaxis::SetMaxDigits(2);

  Double_t stops[5] = {0.00,0.34,0.61,0.84,1.00};
  Double_t red[5] = {0.00,0.00,0.87,1.00,0.51};
  Double_t green[5] = {0.00,0.81,1.00,0.20,0.00};
  Double_t blue[5] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);
  
//  TFile* fileActual = new TFile("data/dispOutput_MicroBooNE_E273.root");
//  TTree* treeActual = (TTree*)fileActual->Get("SpaCEtree_bkwdDisp");
//  //TProfile2D* profileActual = new TProfile2D("profileActual","",numDispDivisions_x+1,-1.0*Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),Ly+Ly/(2.0*((Double_t) numDispDivisions_y)));
//  TProfile2D* profileActual = new TProfile2D("profileActual","",numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
//  //treeActual->Draw(Form("100.0*%s:y_reco:x_reco>>profileActual",dispCompString),Form("z_reco < %f && z_reco > %f && elecFate == 1",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
//  if(dispComp == 1) {
//    treeActual->Draw(Form("-100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileActual",dispCompString),Form("z_reco < %f && z_reco > %f && elecFate == 1",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
//  }
//  else {
//    treeActual->Draw(Form("100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileActual",dispCompString),Form("z_reco < %f && z_reco > %f && elecFate == 1",zpointVal+0.025,zpointVal-0.025),"profCOLZ");
//  }
//  profileActual->SetStats(0);
  
  TFile* fileCalib = new TFile("output.root");
  TTree* treeCalib = (TTree*)fileCalib->Get("SpaCEtree_calibFaces");

  TProfile2D* profileCalibTop = new TProfile2D("profileCalibTop","",numDispDivisions_z+1,-Lz/(2.0*((Double_t) numDispDivisions_z)),Lz+Lz/(2.0*((Double_t) numDispDivisions_z)),numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)));
  TProfile2D* profileActualTop = (TProfile2D*)profileCalibTop->Clone("profileActualTop");
  TProfile2D* profileTrackerTop = (TProfile2D*)profileCalibTop->Clone("profileTrackerTop");
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s:(2.5-x_reco):z_reco>>profileCalibTop",dispCompString),"faceNum==0 && elecFate>=0","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s:(2.5-x_reco):z_reco>>profileCalibTop",dispCompString),"faceNum==0 && elecFate>=0","profCOLZ");
  }
  profileCalibTop->SetStats(0);
  
  TProfile2D* profileCalibBottom = new TProfile2D("profileCalibBottom","",numDispDivisions_z+1,-Lz/(2.0*((Double_t) numDispDivisions_z)),Lz+Lz/(2.0*((Double_t) numDispDivisions_z)),numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)));
  TProfile2D* profileActualBottom = (TProfile2D*)profileCalibBottom->Clone("profileActualBottom");
  TProfile2D* profileTrackerBottom = (TProfile2D*)profileCalibBottom->Clone("profileTrackerBottom");
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s:(2.5-x_reco):z_reco>>profileCalibBottom",dispCompString),"faceNum==1 && elecFate>=0","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s:(2.5-x_reco):z_reco>>profileCalibBottom",dispCompString),"faceNum==1 && elecFate>=0","profCOLZ");
  }
  profileCalibBottom->SetStats(0);

  TProfile2D* profileCalibUpstream = new TProfile2D("profileCalibUpstream","",numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  TProfile2D* profileActualUpstream = (TProfile2D*)profileCalibUpstream->Clone("profileActualUpstream");
  TProfile2D* profileTrackerUpstream = (TProfile2D*)profileCalibUpstream->Clone("profileTrackerUpstream");
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileCalibUpstream",dispCompString),"faceNum==2 && elecFate>=0","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileCalibUpstream",dispCompString),"faceNum==2 && elecFate>=0","profCOLZ");
  }
  profileCalibUpstream->SetStats(0);

  TProfile2D* profileCalibDownstream = new TProfile2D("profileCalibDownstream","",numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  TProfile2D* profileActualDownstream = (TProfile2D*)profileCalibDownstream->Clone("profileActualDownstream");
  TProfile2D* profileTrackerDownstream = (TProfile2D*)profileCalibDownstream->Clone("profileTrackerDownstream");
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileCalibDownstream",dispCompString),"faceNum==3 && elecFate>=0","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s:(y_reco-1.25):(2.5-x_reco)>>profileCalibDownstream",dispCompString),"faceNum==3 && elecFate>=0","profCOLZ");
  }
  profileCalibDownstream->SetStats(0);

  TProfile2D* profileCalibCathode = new TProfile2D("profileCalibCathode","",numDispDivisions_z+1,-Lz/(2.0*((Double_t) numDispDivisions_z)),Lz+Lz/(2.0*((Double_t) numDispDivisions_z)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  TProfile2D* profileActualCathode = (TProfile2D*)profileCalibCathode->Clone("profileActualCathode");
  TProfile2D* profileTrackerCathode = (TProfile2D*)profileCalibCathode->Clone("profileTrackerCathode");
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s:(y_reco-1.25):z_reco>>profileCalibCathode",dispCompString),"faceNum==4 && elecFate>=0","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s:(y_reco-1.25):z_reco>>profileCalibCathode",dispCompString),"faceNum==4 && elecFate>=0","profCOLZ");
  }
  profileCalibCathode->SetStats(0);

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

  while (readerActual.Next())
  {
    if(*elecFate_actual == 1) {
      if(*yReco_actual > profileTrackerTop->GetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1))
      {
	profileTrackerTop->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,*yReco_actual);
	profileTrackerTop->SetBinEntries(profileTrackerTop->GetBin(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1),1);
	      
	if (dispComp == 1)
	{
	  profileActualTop->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,-100.0*(*Dx_actual));
	}
	else if (dispComp == 2)
	{
	  profileActualTop->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,100.0*(*Dy_actual));
	}
	else if (dispComp == 3)
	{
	  profileActualTop->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,100.0*(*Dz_actual));
	}
	profileActualTop->SetBinEntries(profileActualTop->GetBin(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1),1);
      }

      if((Ly-*yReco_actual) > profileTrackerBottom->GetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1))
      {
	profileTrackerBottom->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,Ly-*yReco_actual);
	profileTrackerBottom->SetBinEntries(profileTrackerBottom->GetBin(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1),1);
	      
	if (dispComp == 1)
	{
	  profileActualBottom->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,-100.0*(*Dx_actual));
	}
	else if (dispComp == 2)
	{
	  profileActualBottom->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,100.0*(*Dy_actual));
	}
	else if (dispComp == 3)
	{
	  profileActualBottom->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,100.0*(*Dz_actual));
	}
	profileActualBottom->SetBinEntries(profileActualBottom->GetBin(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1),1);
      }
      
      if((Lz-*zReco_actual) > profileTrackerUpstream->GetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1))
      {
	profileTrackerUpstream->SetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,Lz-*zReco_actual);
	profileTrackerUpstream->SetBinEntries(profileTrackerUpstream->GetBin(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1),1);
	      
	if (dispComp == 1)
	{
	  profileActualUpstream->SetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,-100.0*(*Dx_actual));
	}
	else if (dispComp == 2)
	{
	  profileActualUpstream->SetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,100.0*(*Dy_actual));
	}
	else if (dispComp == 3)
	{
	  profileActualUpstream->SetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,100.0*(*Dz_actual));
	}
	profileActualUpstream->SetBinEntries(profileActualUpstream->GetBin(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1),1);
      }

      if(*zReco_actual > profileTrackerDownstream->GetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1))
      {
	profileTrackerDownstream->SetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,*zReco_actual);
	profileTrackerDownstream->SetBinEntries(profileTrackerDownstream->GetBin(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1),1);
	      
	if (dispComp == 1)
	{
	  profileActualDownstream->SetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,-100.0*(*Dx_actual));
	}
	else if (dispComp == 2)
	{
	  profileActualDownstream->SetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,100.0*(*Dy_actual));
	}
	else if (dispComp == 3)
	{
	  profileActualDownstream->SetBinContent(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,100.0*(*Dz_actual));
	}
	profileActualDownstream->SetBinEntries(profileActualDownstream->GetBin(TMath::Nint(numDispDivisions_x*(Lx-*xReco_actual)/Lx)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1),1);
      }

      if((Lx-*xReco_actual) > profileTrackerCathode->GetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1))
      {
	profileTrackerCathode->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,Lx-*xReco_actual);
	profileTrackerCathode->SetBinEntries(profileTrackerCathode->GetBin(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1),1);
	      
	if (dispComp == 1)
	{
	  profileActualCathode->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,-100.0*(*Dx_actual));
	}
	else if (dispComp == 2)
	{
	  profileActualCathode->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,100.0*(*Dy_actual));
	}
	else if (dispComp == 3)
	{
	  profileActualCathode->SetBinContent(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1,100.0*(*Dz_actual));
	}
	profileActualCathode->SetBinEntries(profileActualCathode->GetBin(TMath::Nint(numDispDivisions_z*(*zReco_actual)/Lz)+1,TMath::Nint(numDispDivisions_y*(*yReco_actual)/Ly)+1),1);
      }
    }
  }
  
  //TProfile2D* profileDiff = new TProfile2D("profileDiff","",numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  //for(Int_t i = 1; i <= numDispDivisions_x+1; i++)
  //{
  //  for(Int_t j = 1; j <= numDispDivisions_y+1; j++)
  //  {
  //    if(fabs(profileActual->GetBinContent(i,j)) == 0.0)
  //      profileDiff->SetBinContent(i,j,0.0);
  //    else
  //      profileDiff->SetBinContent(i,j,logFactor*(1.0*profileCalib->GetBinContent(i,j)-1.0*profileActual->GetBinContent(i,j)));
  //
  //    profileDiff->SetBinEntries(profileDiff->GetBin(i,j),1);
  //  }
  //}
  //profileDiff->SetStats(0);  
    
  //TProfile2D* profileRatio = (TProfile2D*)profileDiff->Clone("profileRatio");
  //profileRatio->Divide(profileActual);
  //for(Int_t i = 1; i <= numDispDivisions_x+1; i++)
  //{
  //  for(Int_t j = 1; j <= numDispDivisions_y+1; j++)
  //  {
  //    if(fabs(profileActual->GetBinContent(i,j)) < 0.001)
  //      profileRatio->SetBinContent(i,j,0.0);
  //
  //    if(fabs(profileActual->GetBinContent(i,j)) < 0.00001)
  //      profileActual->SetBinContent(i,j,0.0);
  //  }
  //}
  //profileRatio->SetStats(0);  

  gPad->Clear();
  gStyle->SetTitleW(0.9);

  //if(whichPlot == 1)
  //  profilePlot = (TProfile2D*)profileActual->Clone("profilePlot");
  //else if(whichPlot == 2)
  //  profilePlot = (TProfile2D*)profileCalib->Clone("profilePlot");
  //else if(whichPlot == 3)
  //  profilePlot = (TProfile2D*)profileDiff->Clone("profilePlot");
  //else if(whichPlot == 4)
  //  profilePlot = (TProfile2D*)profileRatio->Clone("profilePlot");
  //
  //if(whichPlot == 1)
  //  profilePlot->SetTitle(Form("Actual %s [cm]:  Z = %.2f m",dispCompString2,zpointVal*(5.18/5.0)));
  //else if(whichPlot == 2)
  //  profilePlot->SetTitle(Form("Calibrated %s [cm]:  Z = %.2f m",dispCompString2,zpointVal*(5.18/5.0)));
  //else if(whichPlot == 3)
  //  profilePlot->SetTitle(Form("Calibrated-Actual %s [cm]:  Z = %.2f m",dispCompString2,zpointVal*(5.18/5.0)));
  //else if(whichPlot == 4)
  //  profilePlot->SetTitle(Form("(Calibrated-Actual)/Actual %s:  Z = %.2f m",dispCompString2,zpointVal*(5.18/5.0)));

  TProfile2D* profilePlot = new TProfile2D;
  if(faceNum == 0)
  {
    if(whichPlot == 1) {
      profilePlot = (TProfile2D*)profileActualTop->Clone("profilePlot");
    }
    else if(whichPlot == 2) {
      profilePlot = (TProfile2D*)profileCalibTop->Clone("profilePlot");
    }
    profilePlot->GetXaxis()->SetTitle("Z_{reco} [m]");
    profilePlot->GetXaxis()->SetTitleOffset(1.0);
    profilePlot->GetYaxis()->SetTitle("X_{reco} [m]");
    profilePlot->GetYaxis()->SetTitleOffset(0.4);
  }
  else if(faceNum == 1)
  {
    if(whichPlot == 1) {
      profilePlot = (TProfile2D*)profileActualBottom->Clone("profilePlot");
    }
    else if(whichPlot == 2) {
      profilePlot = (TProfile2D*)profileCalibBottom->Clone("profilePlot");
    }
    profilePlot->GetXaxis()->SetTitle("Z_{reco} [m]");
    profilePlot->GetXaxis()->SetTitleOffset(1.0);
    profilePlot->GetYaxis()->SetTitle("X_{reco} [m]");
    profilePlot->GetYaxis()->SetTitleOffset(0.4);
  }
  else if(faceNum == 2)
  {
    if(whichPlot == 1) {
      profilePlot = (TProfile2D*)profileActualUpstream->Clone("profilePlot");
    }
    else if(whichPlot == 2) {
      profilePlot = (TProfile2D*)profileCalibUpstream->Clone("profilePlot");
    }
    profilePlot->GetXaxis()->SetTitle("X_{reco} [m]");
    profilePlot->GetXaxis()->SetTitleOffset(1.0);
    profilePlot->GetYaxis()->SetTitle("Y_{reco} [m]");
    profilePlot->GetYaxis()->SetTitleOffset(1.2);
  }
  else if(faceNum == 3)
  {
    if(whichPlot == 1) {
      profilePlot = (TProfile2D*)profileActualDownstream->Clone("profilePlot");
    }
    else if(whichPlot == 2) {
      profilePlot = (TProfile2D*)profileCalibDownstream->Clone("profilePlot");
    }
    profilePlot->GetXaxis()->SetTitle("X_{reco} [m]");
    profilePlot->GetXaxis()->SetTitleOffset(1.0);
    profilePlot->GetYaxis()->SetTitle("Y_{reco} [m]");
    profilePlot->GetYaxis()->SetTitleOffset(1.2);
  }
  else if(faceNum == 4)
  {
    if(whichPlot == 1) {
      profilePlot = (TProfile2D*)profileActualCathode->Clone("profilePlot");
    }
    else if(whichPlot == 2) {
      profilePlot = (TProfile2D*)profileCalibCathode->Clone("profilePlot");
    }
    profilePlot->GetXaxis()->SetTitle("Z_{reco} [m]");
    profilePlot->GetXaxis()->SetTitleOffset(1.0);
    profilePlot->GetYaxis()->SetTitle("Y_{reco} [m]");
    profilePlot->GetYaxis()->SetTitleOffset(0.4);
  }

  profilePlot->GetXaxis()->SetTitleSize(0.04);
  profilePlot->GetYaxis()->SetTitleSize(0.04);
  profilePlot->GetZaxis()->SetNoExponent(kTRUE);
  profilePlot->GetZaxis()->SetLabelSize(0.03);
  profilePlot->SetStats(0);
  profilePlot->Draw("COLZ");

  c->Update();
  TPaletteAxis *palette = (TPaletteAxis*)profilePlot->GetListOfFunctions()->FindObject("palette");
  if((faceNum < 2) || (faceNum > 3)) {
    palette->SetX1NDC(0.9025);
    palette->SetX2NDC(0.9125);
    profilePlot->GetZaxis()->SetTickLength(0.012);
  }
  else {
    palette->SetX1NDC(0.905);
    palette->SetX2NDC(0.93);
    profilePlot->GetZaxis()->SetTickLength(0.03);
  }
  c->Modified();
}
