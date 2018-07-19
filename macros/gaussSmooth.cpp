#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TProfile2D.h>
#include <TGraph.h>
#include <TView.h>
#include <TView3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TVector3.h>
#include <TVirtualFFT.h>
#include <TSystem.h>
#include <TGraph2D.h>

#include "spline.h"

using namespace std;

const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;

const Int_t minEntries = 0;

TH2F* getFaceHist(Int_t numDispDivisions, Int_t faceNum, Int_t dispComp);
void extrapolate(TH2F *interpHist, Int_t edgeXbins, Int_t edgeYbins, Int_t faceNum);

int main(int argc, char** argv)
{
  Int_t edgeXbins = 0;
  Int_t edgeYbins = 0;
  Int_t medianFiltXbins = 1;
  Int_t medianFiltYbins = 1;
  Double_t gaussFilterWidth = 0.10;
  Int_t faceNum = 2;
  Int_t dispComp = 1;
  
  gErrorIgnoreLevel = kError;

  TH2F *tempHist = getFaceHist(25,faceNum,dispComp);
  Int_t numXbins = tempHist->GetNbinsX()-2*edgeXbins;
  Int_t numYbins = tempHist->GetNbinsY()-2*edgeYbins;
  
  TH2F *origHist = new TH2F("origHist","",numXbins,0,numXbins,numYbins,0,numYbins);
  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      origHist->SetBinContent(i+1,j+1,tempHist->GetBinContent(i+1+edgeXbins,j+1+edgeYbins));
    }
  }

  Int_t numFiltXbins = origHist->GetNbinsX()-2*medianFiltXbins;
  Int_t numFiltYbins = origHist->GetNbinsY()-2*medianFiltYbins;
  TH2F *medianHist = new TH2F("medianHist","",numFiltXbins,0,numFiltXbins,numFiltYbins,0,numFiltYbins);
  for(Int_t i = medianFiltXbins; i < numXbins-medianFiltXbins; i++)
  {
    for(Int_t j = medianFiltYbins; j < numYbins-medianFiltYbins; j++)
    {
      Double_t vals[(2*medianFiltXbins+1)*(2*medianFiltYbins+1)];
      Double_t weights[(2*medianFiltXbins+1)*(2*medianFiltYbins+1)];

      Int_t numCounts = 0;
      for(Int_t h = -1*medianFiltXbins; h <= medianFiltXbins; h++)
      {
        for(Int_t k = -1*medianFiltYbins; k <= medianFiltYbins; k++)
        {
          vals[numCounts] = origHist->GetBinContent(i+h+1,j+k+1);
          if(origHist->GetBinContent(i+h+1,j+k+1) > -1000.0)
          {
            weights[numCounts] = 1.0;
          }
	  else
	  {
            weights[numCounts] = 0.0;
	  }
	  numCounts++;
        }
      }

      if(numCounts >= 5)
      {
        medianHist->SetBinContent(i-medianFiltXbins+1,j-medianFiltYbins+1,TMath::Median(numCounts,&vals[0],&weights[0]));
      }
      else
      {
        medianHist->SetBinContent(i-medianFiltXbins+1,j-medianFiltYbins+1,-10000.0);
      }
    }
  }

  TH2F *filledHist = new TH2F("filledHist","",numFiltXbins,0,numFiltXbins,numFiltYbins,0,numFiltYbins);
  TGraph2D *interpGraph = new TGraph2D();
  Int_t pointNum = 0;
  for(Int_t i = 0; i < numFiltXbins; i++)
  {
    for(Int_t j = 0; j < numFiltYbins; j++)
    {
      if(medianHist->GetBinContent(i+1,j+1) > -1000.0)
      {
        interpGraph->SetPoint(pointNum,medianHist->GetXaxis()->GetBinCenter(i+1),medianHist->GetYaxis()->GetBinCenter(j+1),medianHist->GetBinContent(i+1,j+1));
	pointNum++;
      }
    }
  }
  for(Int_t i = 0; i < numFiltXbins; i++)
  {
    for(Int_t j = 0; j < numFiltYbins; j++)
    {
      filledHist->SetBinContent(i+1,j+1,interpGraph->Interpolate(filledHist->GetXaxis()->GetBinCenter(i+1),filledHist->GetYaxis()->GetBinCenter(j+1)));
    }
  }
  
  TH2F *expandHist = new TH2F("expandHist","",3*numFiltXbins,0,3*numFiltXbins,3*numFiltYbins,0,3*numFiltYbins);
  for(Int_t i = 0; i < 3*numFiltXbins; i++)
  {
    for(Int_t j = 0; j < 3*numFiltYbins; j++)
    {
      expandHist->SetBinContent(i+1,j+1,filledHist->GetBinContent(fabs(numFiltXbins-0.5-(i%(2*numFiltXbins)))+0.5,fabs(numFiltYbins-0.5-(j%(2*numFiltYbins)))+0.5));
    }
  }
  
  TH1F **inputHist = new TH1F*[3*numFiltYbins];
  for(Int_t j = 0; j < 3*numFiltYbins; j++)
  {
    inputHist[j] = new TH1F(Form("inputHist_%d",j),"",3*numFiltXbins,0,3*numFiltXbins);
    
    for(Int_t i = 0; i < 3*numFiltXbins; i++)
    {
      inputHist[j]->SetBinContent(i+1,expandHist->GetBinContent(i+1,j+1));
    }
  }

  TH2F *resultHist = new TH2F("resultHist","",3*numFiltXbins,0,3*numFiltXbins,3*numFiltYbins,0,3*numFiltYbins);
  TH2F *cropHist = new TH2F("cropHist","",numFiltXbins,0,numFiltXbins,numFiltYbins,0,numFiltYbins);

  // Do preliminary setup
  TVirtualFFT::SetTransform(0);
  Int_t numTotalX = 3*numFiltXbins;
  Int_t numTotalY = 3*numFiltYbins;

  // Carry out FFT for convoluted signal
  auto rho = new Double_t[500][500];
  auto ph = new Double_t[500][500];
  for(Int_t i = 0; i < numTotalY; i++)
  {
    TH1 *hm1 = 0;
    hm1 = inputHist[i]->FFT(hm1,"MAG");
    TH1 *hp1 = 0;
    hp1 = inputHist[i]->FFT(hp1,"PH");
    for(Int_t j = 0; j < numTotalX; j++)
    {
      rho[i][j] = hm1->GetBinContent(j+1);
      ph[i][j] = hp1->GetBinContent(j+1);
    }
    delete hm1;
    delete hp1;
  }

  // Do deconvolution for given frequency
  auto result_re = new Double_t[500][500];
  auto result_im = new Double_t[500][500];
  Double_t value_re[numTotalY];
  Double_t value_im[numTotalY];
  Double_t temp_re[numTotalY];
  Double_t temp_im[numTotalY];
  Double_t temp2_re[numTotalY];
  Double_t temp2_im[numTotalY];
  TVirtualFFT *ifft;
  TVirtualFFT *ifft2;
  TF1 *gaussFunc = new TF1("gaussFunc","gaus");
  gaussFunc->SetParameter(0,1.0);
  gaussFunc->SetParameter(1,0.0);
  gaussFunc->SetParameter(2,gaussFilterWidth);
  Double_t filterFactor;
  for(Int_t i = 0; i < numTotalX; i++)
  {
    // Build vectors in wire-number-space for given frequency
    for(Int_t j = 0; j < numTotalY; j++)
    {
      value_re[j] = rho[j][i]*cos(ph[j][i]);
      value_im[j] = rho[j][i]*sin(ph[j][i]);
    }
  
    // Do secondary FFT over wire number for convoluted signal
    ifft = TVirtualFFT::FFT(1,&numTotalY,"C2CFORWARD M K");
    ifft->SetPointsComplex(value_re,value_im);
    ifft->Transform();
    ifft->GetPointsComplex(temp_re,temp_im);

    for(Int_t j = 0; j < numTotalY; j++)
    {
      filterFactor = gaussFunc->Eval(TMath::Min(((Double_t) i),((Double_t) numTotalX-i))/((Double_t) numTotalX))*gaussFunc->Eval(TMath::Min(((Double_t) j),((Double_t) numTotalY-j))/((Double_t) numTotalY)); // define filter here
      
      temp_re[j] = filterFactor*temp_re[j]/numTotalY;
      temp_im[j] = filterFactor*temp_im[j]/numTotalY;
    }

    // Do anti-FFT in wire-number space
    ifft2 = TVirtualFFT::FFT(1,&numTotalY,"C2CBACKWARD M K");
    ifft2->SetPointsComplex(temp_re,temp_im);
    ifft2->Transform();
    ifft2->GetPointsComplex(temp2_re,temp2_im);
  
    // Prepare deconvoluted results for given frequency (on all wires)
    for(Int_t j = 0; j < numTotalY; j++)
    {
      result_re[j][i] = temp2_re[j]/numTotalX;
      result_im[j][i] = temp2_im[j]/numTotalX;
    }
  }
  delete ifft;
  delete ifft2;

  // Final anti-FFT to go to time domain and store final results
  TVirtualFFT *ifft3;
  for(Int_t i = 0; i < numTotalY; i++)
  {
    ifft3 = TVirtualFFT::FFT(1,&numTotalX,"C2R M K");
    ifft3->SetPointsComplex(result_re[i],result_im[i]);
    ifft3->Transform();
  
    TH1 *fb = 0;
    fb = TH1::TransformHisto(ifft3,fb,"Re");
  
    for(Int_t j = 0; j < numTotalX; j++)
    {
      Double_t content = fb->GetBinContent(j+1);
      resultHist->SetBinContent(j+1,i+1,content); // fill results here
    }
    delete fb;
  }
  delete ifft3;
  delete[] rho;
  delete[] ph;
  delete[] result_re;
  delete[] result_im;

  for(Int_t i = 0; i < numFiltXbins; i++)
  {
    for(Int_t j = 0; j < numFiltYbins; j++)
    {
      cropHist->SetBinContent(i+1,j+1,resultHist->GetBinContent(i+1+numFiltXbins,j+1+numFiltYbins));
    }
  }

  // Interpolate to fix boundaries of map
  TH2F *interpHist = new TH2F("interpHist","",tempHist->GetNbinsX(),0,tempHist->GetNbinsX(),tempHist->GetNbinsY(),0,tempHist->GetNbinsY());
  for(Int_t i = (edgeXbins+medianFiltXbins); i < numFiltXbins+(edgeXbins+medianFiltXbins); i++)
  {
    for(Int_t j = (edgeYbins+medianFiltYbins); j < numFiltYbins+(edgeYbins+medianFiltYbins); j++)
    {
      interpHist->SetBinContent(i+1,j+1,cropHist->GetBinContent(i+1-(edgeXbins+medianFiltXbins),j+1-(edgeYbins+medianFiltYbins)));
    }
  }

  extrapolate(interpHist,(edgeXbins+medianFiltXbins),(edgeYbins+medianFiltYbins),faceNum);

  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      if(origHist->GetBinContent(i+1,j+1) <= -1000.0)
      {
  	origHist->SetBinContent(i+1,j+1,0.0);
      }
    }
  }

  for(Int_t i = 0; i < numFiltXbins; i++)
  {
    for(Int_t j = 0; j < numFiltYbins; j++)
    {
      if(medianHist->GetBinContent(i+1,j+1) <= -1000.0)
      {
	medianHist->SetBinContent(i+1,j+1,0.0);
      }
    }
  }

  TFile outputfile("output_smooth.root","RECREATE");
  outputfile.cd();
  origHist->Write();
  medianHist->Write();
  filledHist->Write();
  expandHist->Write();
  resultHist->Write();
  cropHist->Write();
  interpHist->Write();
  outputfile.Close();
  
  return 0;
}

TH2F* getFaceHist(Int_t numDispDivisions, Int_t faceNum, Int_t dispComp)
{
  Double_t numDispDivisions_x = numDispDivisions;
  Double_t numDispDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)numDispDivisions));
  Double_t numDispDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)numDispDivisions));

  Char_t *dispCompString = (Char_t*) "";
  if(dispComp == 1)
  {
    dispCompString = (Char_t*) "Dx";
  }
  else if(dispComp == 2)
  {
    dispCompString = (Char_t*) "Dy";
  }
  else if(dispComp == 3)
  {
    dispCompString = (Char_t*) "Dz";
  }

  TFile* fileCalib = new TFile("output.root","READ");
  TTree* treeCalib = (TTree*)fileCalib->Get("SpaCEtree_calibFaces");

  TProfile2D* profileCalibTop = new TProfile2D("profileCalibTop","",numDispDivisions_z+1,-Lz/(2.0*((Double_t) numDispDivisions_z)),Lz+Lz/(2.0*((Double_t) numDispDivisions_z)),numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)));
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s-(numEntries<%d)*10000.0:(2.5-x_reco):z_reco>>profileCalibTop",dispCompString,minEntries),"faceNum==0","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s-(numEntries<%d)*10000.0:(2.5-x_reco):z_reco>>profileCalibTop",dispCompString,minEntries),"faceNum==0","profCOLZ");
  }
  profileCalibTop->SetStats(0);

  TProfile2D* profileCalibBottom = new TProfile2D("profileCalibBottom","",numDispDivisions_z+1,-Lz/(2.0*((Double_t) numDispDivisions_z)),Lz+Lz/(2.0*((Double_t) numDispDivisions_z)),numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)));
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s-(numEntries<%d)*10000.0:(2.5-x_reco):z_reco>>profileCalibBottom",dispCompString,minEntries),"faceNum==1","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s-(numEntries<%d)*10000.0:(2.5-x_reco):z_reco>>profileCalibBottom",dispCompString,minEntries),"faceNum==1","profCOLZ");
  }
  profileCalibBottom->SetStats(0);

  TProfile2D* profileCalibUpstream = new TProfile2D("profileCalibUpstream","",numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s-(numEntries<%d)*10000.0:(y_reco-1.25):(2.5-x_reco)>>profileCalibUpstream",dispCompString,minEntries),"faceNum==2","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s-(numEntries<%d)*10000.0:(y_reco-1.25):(2.5-x_reco)>>profileCalibUpstream",dispCompString,minEntries),"faceNum==2","profCOLZ");
  }
  profileCalibUpstream->SetStats(0);

  TProfile2D* profileCalibDownstream = new TProfile2D("profileCalibDownstream","",numDispDivisions_x+1,-Lx/(2.0*((Double_t) numDispDivisions_x)),Lx+Lx/(2.0*((Double_t) numDispDivisions_x)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s-(numEntries<%d)*10000.0:(y_reco-1.25):(2.5-x_reco)>>profileCalibDownstream",dispCompString,minEntries),"faceNum==3","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s-(numEntries<%d)*10000.0:(y_reco-1.25):(2.5-x_reco)>>profileCalibDownstream",dispCompString,minEntries),"faceNum==3","profCOLZ");
  }
  profileCalibDownstream->SetStats(0);

  TProfile2D* profileCalibCathode = new TProfile2D("profileCalibCathode","",numDispDivisions_z+1,-Lz/(2.0*((Double_t) numDispDivisions_z)),Lz+Lz/(2.0*((Double_t) numDispDivisions_z)),numDispDivisions_y+1,-1.0*(Ly/2.0)-1.0*Ly/(2.0*((Double_t) numDispDivisions_y)),(Ly/2.0)+Ly/(2.0*((Double_t) numDispDivisions_y)));
  if(dispComp == 1) {
    treeCalib->Draw(Form("-100.0*%s-(numEntries<%d)*10000.0:(y_reco-1.25):z_reco>>profileCalibCathode",dispCompString,minEntries),"faceNum==4","profCOLZ");
  }
  else {
    treeCalib->Draw(Form("100.0*%s-(numEntries<%d)*10000.0:(y_reco-1.25):z_reco>>profileCalibCathode",dispCompString,minEntries),"faceNum==4","profCOLZ");
  }
  profileCalibCathode->SetStats(0);

  TProfile2D* profilePlot = new TProfile2D;
  if(faceNum == 0)
  {
    profilePlot = (TProfile2D*)profileCalibTop->Clone("profilePlot");
  }
  else if(faceNum == 1)
  {
    profilePlot = (TProfile2D*)profileCalibBottom->Clone("profilePlot");
  }
  else if(faceNum == 2)
  {
    profilePlot = (TProfile2D*)profileCalibUpstream->Clone("profilePlot");
  }
  else if(faceNum == 3)
  {
    profilePlot = (TProfile2D*)profileCalibDownstream->Clone("profilePlot");
  }
  else if(faceNum == 4)
  {
    profilePlot = (TProfile2D*)profileCalibCathode->Clone("profilePlot");
  }

  return (TH2F*) profilePlot->ProjectionXY();
}

void extrapolate(TH2F *interpHist, Int_t edgeXbins, Int_t edgeYbins, Int_t faceNum)
{
  for(Int_t i = edgeXbins; i < interpHist->GetNbinsX()-edgeXbins; i++)
  {
    std::vector<Double_t> X, Y;
    tk::spline s;

    if(((faceNum == 0) || (faceNum == 1)) && (edgeYbins != 0))
    {
      X.push_back(interpHist->GetYaxis()->GetBinCenter(1));
      Y.push_back(0.0);
    }
    
    for(Int_t j = edgeYbins; j < interpHist->GetNbinsY()-edgeYbins; j++)
    {
      X.push_back(interpHist->GetYaxis()->GetBinCenter(j+1));
      Y.push_back(interpHist->GetBinContent(i+1,j+1));
    }

    s.set_points(X,Y);

    for(Int_t k = 0; k < edgeYbins; k++)
    {
      interpHist->SetBinContent(i+1,k+1,s(interpHist->GetYaxis()->GetBinCenter(k+1)));
      interpHist->SetBinContent(i+1,interpHist->GetNbinsY()-k,s(interpHist->GetYaxis()->GetBinCenter(interpHist->GetNbinsY()-k)));
    }
  }

  for(Int_t j = edgeYbins; j < interpHist->GetNbinsY()-edgeYbins; j++)
  {
    std::vector<Double_t> X, Y;
    tk::spline s;

    if(((faceNum == 2) || (faceNum == 3)) && (edgeXbins != 0))
    {
      X.push_back(interpHist->GetXaxis()->GetBinCenter(1));
      Y.push_back(0.0);
    }
    
    for(Int_t i = edgeXbins; i < interpHist->GetNbinsX()-edgeXbins; i++)
    {
      X.push_back(interpHist->GetXaxis()->GetBinCenter(i+1));
      Y.push_back(interpHist->GetBinContent(i+1,j+1));
    }

    s.set_points(X,Y);

    for(Int_t k = 0; k < edgeXbins; k++)
    {
      interpHist->SetBinContent(k+1,j+1,s(interpHist->GetXaxis()->GetBinCenter(k+1)));
      interpHist->SetBinContent(interpHist->GetNbinsX()-k,j+1,s(interpHist->GetXaxis()->GetBinCenter(interpHist->GetNbinsX()-k)));
    }
  }

  for(Int_t m = 0; m < edgeXbins; m++)
  {
    for(Int_t n = 0; n < edgeYbins; n++)
    {
      std::vector<Double_t> X1, Y1;
      tk::spline s1;

      std::vector<Double_t> X2, Y2;
      tk::spline s2;

      std::vector<Double_t> X3, Y3;
      tk::spline s3;

      std::vector<Double_t> X4, Y4;
      tk::spline s4;

      if(((faceNum == 0) || (faceNum == 1)) && (edgeYbins != 0))
      {
        X1.push_back(interpHist->GetYaxis()->GetBinCenter(1));
        Y1.push_back(0.0);
      }

      for(Int_t j = edgeYbins; j < interpHist->GetNbinsY()-edgeYbins; j++)
      {
        X1.push_back(interpHist->GetYaxis()->GetBinCenter(j+1));
        Y1.push_back(interpHist->GetBinContent(m+1,j+1));
      }

      s1.set_points(X1,Y1);

      if(((faceNum == 2) || (faceNum == 3)) && (edgeXbins != 0))
      {
        X2.push_back(interpHist->GetXaxis()->GetBinCenter(1));
        Y2.push_back(0.0);
      }

      for(Int_t i = edgeXbins; i < interpHist->GetNbinsX()-edgeXbins; i++)
      {
        X2.push_back(interpHist->GetXaxis()->GetBinCenter(i+1));
        Y2.push_back(interpHist->GetBinContent(i+1,n+1));	
      }
      
      s2.set_points(X2,Y2);

      if(((faceNum == 0) || (faceNum == 1)) && (edgeYbins != 0))
      {
        X3.push_back(interpHist->GetYaxis()->GetBinCenter(1));
        Y3.push_back(0.0);
      }

      for(Int_t j = edgeYbins; j < interpHist->GetNbinsY()-edgeYbins; j++)
      {
        X3.push_back(interpHist->GetYaxis()->GetBinCenter(j+1));
        Y3.push_back(interpHist->GetBinContent(interpHist->GetNbinsX()-m,j+1));
      }

      s3.set_points(X3,Y3);

      if(((faceNum == 2) || (faceNum == 3)) && (edgeXbins != 0))
      {
        X4.push_back(interpHist->GetXaxis()->GetBinCenter(1));
        Y4.push_back(0.0);
      }

      for(Int_t i = edgeXbins; i < interpHist->GetNbinsX()-edgeXbins; i++)
      {
        X4.push_back(interpHist->GetXaxis()->GetBinCenter(i+1));
        Y4.push_back(interpHist->GetBinContent(i+1,interpHist->GetNbinsY()-n));	
      }
      
      s4.set_points(X4,Y4);

      interpHist->SetBinContent(m+1,n+1,(s1(interpHist->GetYaxis()->GetBinCenter(n+1))+s2(interpHist->GetXaxis()->GetBinCenter(m+1)))/2.0);
      interpHist->SetBinContent(m+1,interpHist->GetNbinsY()-n,(s1(interpHist->GetYaxis()->GetBinCenter(interpHist->GetNbinsY()-n))+s4(interpHist->GetXaxis()->GetBinCenter(m+1)))/2.0);
      interpHist->SetBinContent(interpHist->GetNbinsX()-m,n+1,(s3(interpHist->GetYaxis()->GetBinCenter(n+1))+s2(interpHist->GetXaxis()->GetBinCenter(interpHist->GetNbinsX()-m)))/2.0);
      interpHist->SetBinContent(interpHist->GetNbinsX()-m,interpHist->GetNbinsY()-n,(s3(interpHist->GetYaxis()->GetBinCenter(interpHist->GetNbinsY()-n))+s4(interpHist->GetXaxis()->GetBinCenter(interpHist->GetNbinsX()-m)))/2.0);
    }
  }

  return;
}
