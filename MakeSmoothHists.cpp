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
#include "TTreeReader.h"

#include "spline.h"
#include <Fade_3D.h>

using namespace std;
using namespace FADE3D;

const Char_t *inputFileName = "output.root"; // NEED TO RUN CALIBSCE FIRST

const Char_t *simInterpFileName = "data/output_siminterp_MicroBooNE_4p5_gap.root";

const Bool_t isMC = true;
//const Bool_t isMC = false;

const Double_t Efield = 0.2739;
const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;

const Int_t numDispDivisions = 25;
const Int_t minMedianPoints = 23;
const Int_t interpMode = 1;
const Int_t numInterpBufferBins = 3;
const Int_t testBin = 50;

const Double_t WF_top = 117.1;
const Double_t WF_bottom = -115.1;
const Double_t WF_upstream = 0.2;
const Double_t WF_downstream = 1036.6;
const Double_t WF_cathode = 254.4;

TFile outputfile("output_smooth.root","RECREATE");

Int_t numDispDivisions_x = numDispDivisions;
Int_t numDispDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)numDispDivisions));
Int_t numDispDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)numDispDivisions));

vector<TH3F*> getCalibHists();
vector<TH3F*> getSimInterpHists();
void smoothOneDim(TH3F *inputHist, TH3F* simInterpHist, Int_t dispComp);
void conditionHist(TH3F* &fixHist, Int_t driftSign, Bool_t doCornerFix, Int_t component);
vector<TH3F*> makeEfieldHists(const vector<TH3F*> &spatialHists, const Int_t driftSign);
void doMedianFiltering(TH3F* &inputHist);
Double_t findElecField(Double_t input);
Point3 getInterpVal(Point3 point, pair<Point3,Point3> vtx1, pair<Point3,Point3> vtx2, pair<Point3,Point3> vtx3, pair<Point3,Point3> vtx4);

Double_t doInvCoordTransformX(const Double_t inputX)
{
  Double_t outputX;
  outputX = 100.0*(2.56/Lx)*(Lx - inputX);

  return outputX;
}

Double_t doInvCoordTransformY(const Double_t inputY)
{
  Double_t outputY;
  outputY = 100.0*(2.33/Ly)*inputY - 116.5;

  return outputY;
}

Double_t doInvCoordTransformZ(const Double_t inputZ)
{
  Double_t outputZ;
  outputZ = 100.0*(10.37/Lz)*inputZ;

  return outputZ;
}

int main(int argc, char** argv)
{
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(2);
  Double_t stops[5] = {0.00,0.34,0.61,0.84,1.00};
  Double_t red[5] = {0.00,0.00,0.87,1.00,0.51};
  Double_t green[5] = {0.00,0.81,1.00,0.20,0.00};
  Double_t blue[5] = {0.51,1.00,0.12,0.00,0.00};
  TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
  gStyle->SetNumberContours(255);
  gStyle->SetOptStat(0);
   
  vector<TH3F*> RecoBkwdHists = getCalibHists();
  vector<TH3F*> TempSimInterpHists = getSimInterpHists();
  
  for(Int_t k = 0; k < 3; k++)
  {
    smoothOneDim(RecoBkwdHists.at(k),TempSimInterpHists.at(k),k+1);
  }

  // STILL NEED TO PRODUCE FORWARD MAPS!!!!

  vector<TH3F*> RecoEfieldHists = makeEfieldHists(RecoBkwdHists,-1);
  
  outputfile.cd();

  for(Int_t k = 0; k < 3; k++)
  {
    (RecoBkwdHists.at(k))->Write();
  }
  for(Int_t k = 0; k < 4; k++)
  {
    (RecoEfieldHists.at(k))->Write();
  }

  outputfile.Close();
    
  return 0;
}

void smoothOneDim(TH3F *inputHist, TH3F* simInterpHist, Int_t dispComp)
{
  const Int_t medianFiltXbins = 1;
  const Int_t medianFiltYbins = 1;
  const Int_t medianFiltZbins = 1;

  const Int_t padXbins = 1;
  const Int_t padYbins = 1;
  const Int_t padZbins = 1;

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

  Int_t numXbins = inputHist->GetNbinsX();
  Int_t numYbins = inputHist->GetNbinsY();
  Int_t numZbins = inputHist->GetNbinsZ();
  
  TH3F *origHist = new TH3F(Form("%s_origHist",dispCompString),"",numXbins,0,numXbins,numYbins,0,numYbins,numZbins,0,numZbins);
  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      for(Int_t k = 0; k < numZbins; k++)
      {
        origHist->SetBinContent(i+1,j+1,k+1,inputHist->GetBinContent(i+1,j+1,k+1));
      }
    }
  }

  TH2F *origHistCenter = new TH2F(Form("%s_origHistCenter",dispCompString),"",numXbins,0,numXbins,numYbins,0,numYbins);
  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      origHistCenter->SetBinContent(i+1,j+1,origHist->GetBinContent(i+1,j+1,testBin));
    }
  }

  TH3F *medianHist = new TH3F(Form("%s_medianHist",dispCompString),"",numXbins,0,numXbins,numYbins,0,numYbins,numZbins,0,numZbins);
  for(Int_t i = medianFiltXbins; i < numXbins-medianFiltXbins; i++)
  {
    for(Int_t j = medianFiltYbins; j < numYbins-medianFiltYbins; j++)
    {
      for(Int_t k = medianFiltZbins; k < numZbins-medianFiltZbins; k++)
      {
        Double_t vals[(2*medianFiltXbins+1)*(2*medianFiltYbins+1)*(2*medianFiltZbins+1)];
        Double_t weights[(2*medianFiltXbins+1)*(2*medianFiltYbins+1)*(2*medianFiltZbins+1)];
      
        Int_t numCounts = 0;
        Int_t numFilled = 0;
        for(Int_t p = -1*medianFiltXbins; p <= medianFiltXbins; p++)
        {
          for(Int_t q = -1*medianFiltYbins; q <= medianFiltYbins; q++)
          {
            for(Int_t r = -1*medianFiltZbins; r <= medianFiltZbins; r++)
            {
              vals[numCounts] = origHist->GetBinContent(i+p+1,j+q+1,k+r+1);
              if(origHist->GetBinContent(i+p+1,j+q+1,k+r+1) > -1000.0)
              {
                weights[numCounts] = 1.0;
		numFilled++;
              }
      	      else
      	      {
                weights[numCounts] = 0.0;
      	      }
      	      numCounts++;
	    }
          }
        }
	
        if(numFilled >= minMedianPoints)
        {
          medianHist->SetBinContent(i+1,j+1,k+1,TMath::Median(numCounts,&vals[0],&weights[0]));
        }
        else
        {
          medianHist->SetBinContent(i+1,j+1,k+1,-10000.0);
        }
      }
    }
  }

  TH2F *medianHistCenter = new TH2F(Form("%s_medianHistCenter",dispCompString),"",numXbins,0,numXbins,numYbins,0,numYbins);
  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      medianHistCenter->SetBinContent(i+1,j+1,medianHist->GetBinContent(i+1,j+1,testBin));
    }
  }

  TH3F *paddedHist = new TH3F(Form("%s_paddedHist",dispCompString),"",numXbins,0,numXbins,numYbins,0,numYbins,numZbins,0,numZbins);
  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      for(Int_t k = 0; k < numZbins; k++)
      {
        if((i < padXbins) || (i > numXbins-padXbins-1) || (j < padYbins) || (j > numYbins-padYbins-1) || (k < padZbins) || (k > numZbins-padZbins-1)) {
          paddedHist->SetBinContent(i+1,j+1,k+1,simInterpHist->GetBinContent(i+1,j+1,k+1));
	}
	else {
          paddedHist->SetBinContent(i+1,j+1,k+1,medianHist->GetBinContent(i+1,j+1,k+1));
	}
      }
    }
  }
  
  TH2F *paddedHistCenter = new TH2F(Form("%s_paddedHistCenter",dispCompString),"",numXbins,0,numXbins,numYbins,0,numYbins);
  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      paddedHistCenter->SetBinContent(i+1,j+1,paddedHist->GetBinContent(i+1,j+1,testBin));
    }
  }

  TH3F *filledHist = new TH3F(Form("%s_filledHist",dispCompString),"",numXbins,0,numXbins,numYbins,0,numYbins,numZbins,0,numZbins);
  TGraph2D **interpGraphXY = new TGraph2D*[numZbins];
  TGraph2D **interpGraphXZ = new TGraph2D*[numYbins];
  TGraph2D **interpGraphYZ = new TGraph2D*[numXbins];

  for(Int_t k = 0; k < numZbins; k++)
  {
    interpGraphXY[k] = new TGraph2D();
    
    Int_t pointNum = 0;
    for(Int_t i = 0; i < numXbins; i++)
    {
      for(Int_t j = 0; j < numYbins; j++)
      {
        if(paddedHist->GetBinContent(i+1,j+1,k+1) > -1000.0)
        {
          interpGraphXY[k]->SetPoint(pointNum,paddedHist->GetXaxis()->GetBinCenter(i+1),paddedHist->GetYaxis()->GetBinCenter(j+1),paddedHist->GetBinContent(i+1,j+1,k+1));
    	  pointNum++;
        }
      }
    }
  }
  
  for(Int_t k = numInterpBufferBins; k < numZbins-numInterpBufferBins; k++)
  {
    for(Int_t i = 0; i < numXbins; i++)
    {
      for(Int_t j = 0; j < numYbins; j++)
      {
        if(paddedHist->GetBinContent(i+1,j+1,k+1) > -1000.0)
        {
          filledHist->SetBinContent(i+1,j+1,k+1,paddedHist->GetBinContent(i+1,j+1,k+1));
	}
	else
	{
          filledHist->SetBinContent(i+1,j+1,k+1,interpGraphXY[k]->Interpolate(filledHist->GetXaxis()->GetBinCenter(i+1),filledHist->GetYaxis()->GetBinCenter(j+1)));
	}
      }
    }
  }

  for(Int_t i = 0; i < numXbins; i++)
  {
    interpGraphYZ[i] = new TGraph2D();
    
    Int_t pointNum = 0;
    for(Int_t j = 0; j < numYbins; j++)
    {
      for(Int_t k = 0; k < numZbins; k++)
      {
        if(paddedHist->GetBinContent(i+1,j+1,k+1) > -1000.0)
        {
          interpGraphYZ[i]->SetPoint(pointNum,paddedHist->GetYaxis()->GetBinCenter(j+1),paddedHist->GetZaxis()->GetBinCenter(k+1),paddedHist->GetBinContent(i+1,j+1,k+1));
    	  pointNum++;
        }
        else if(filledHist->GetBinContent(i+1,j+1,k+1) != 0.0)
        {
          interpGraphYZ[i]->SetPoint(pointNum,filledHist->GetYaxis()->GetBinCenter(j+1),filledHist->GetZaxis()->GetBinCenter(k+1),filledHist->GetBinContent(i+1,j+1,k+1));
    	  pointNum++;	  
        }
      }
    }
  }

  for(Int_t j = 0; j < numYbins; j++)
  {
    interpGraphXZ[j] = new TGraph2D();
    
    Int_t pointNum = 0;
    for(Int_t i = 0; i < numXbins; i++)
    {
      for(Int_t k = 0; k < numZbins; k++)
      {
        if(paddedHist->GetBinContent(i+1,j+1,k+1) > -1000.0)
        {
          interpGraphXZ[j]->SetPoint(pointNum,paddedHist->GetXaxis()->GetBinCenter(i+1),paddedHist->GetZaxis()->GetBinCenter(k+1),paddedHist->GetBinContent(i+1,j+1,k+1));
    	  pointNum++;
        }
	else if(filledHist->GetBinContent(i+1,j+1,k+1) != 0.0)
	{
          interpGraphXZ[j]->SetPoint(pointNum,filledHist->GetXaxis()->GetBinCenter(i+1),filledHist->GetZaxis()->GetBinCenter(k+1),filledHist->GetBinContent(i+1,j+1,k+1));
    	  pointNum++;
	}
      }
    }
  }

  for(Int_t k = 0; k < numInterpBufferBins; k++)
  {
    for(Int_t i = 0; i < numXbins; i++)
    {
      for(Int_t j = 0; j < numYbins; j++)
      {
        if(paddedHist->GetBinContent(i+1,j+1,k+1) > -1000.0)
        {
          filledHist->SetBinContent(i+1,j+1,k+1,paddedHist->GetBinContent(i+1,j+1,k+1));
	}
	else
	{
	  if(interpMode == 1)
	  {
            filledHist->SetBinContent(i+1,j+1,k+1,interpGraphYZ[i]->Interpolate(filledHist->GetYaxis()->GetBinCenter(j+1),filledHist->GetZaxis()->GetBinCenter(k+1)));
	  }
	  else if(interpMode == 2)
	  {
            filledHist->SetBinContent(i+1,j+1,k+1,interpGraphXZ[j]->Interpolate(filledHist->GetXaxis()->GetBinCenter(i+1),filledHist->GetZaxis()->GetBinCenter(k+1)));
	  }
	}
      }
    }
  }
  for(Int_t k = numZbins-numInterpBufferBins; k < numZbins; k++)
  {
    for(Int_t i = 0; i < numXbins; i++)
    {
      for(Int_t j = 0; j < numYbins; j++)
      {
        if(paddedHist->GetBinContent(i+1,j+1,k+1) > -1000.0)
        {
          filledHist->SetBinContent(i+1,j+1,k+1,paddedHist->GetBinContent(i+1,j+1,k+1));
	}
	else
	{
	  if(interpMode == 1)
	  {
            filledHist->SetBinContent(i+1,j+1,k+1,interpGraphYZ[i]->Interpolate(filledHist->GetYaxis()->GetBinCenter(j+1),filledHist->GetZaxis()->GetBinCenter(k+1)));
	  }
	  else if(interpMode == 2)
	  {
            filledHist->SetBinContent(i+1,j+1,k+1,interpGraphXZ[j]->Interpolate(filledHist->GetXaxis()->GetBinCenter(i+1),filledHist->GetZaxis()->GetBinCenter(k+1)));
	  }
	}
      }
    }
  }
  
  TH2F *filledHistCenter = new TH2F(Form("%s_filledHistCenter",dispCompString),"",numXbins,0,numXbins,numYbins,0,numYbins);
  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      filledHistCenter->SetBinContent(i+1,j+1,filledHist->GetBinContent(i+1,j+1,testBin));
    }
  }

  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      for(Int_t k = 0; k < numZbins; k++)
      {
        inputHist->SetBinContent(i+1,j+1,k+1,filledHist->GetBinContent(i+1,j+1,k+1));
      }
    }
  }

  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      for(Int_t k = 0; k < numZbins; k++)
      {
        if(origHist->GetBinContent(i+1,j+1,k+1) <= -1000.0)
        {
      	  origHist->SetBinContent(i+1,j+1,k+1,0.0);
        }
      }

      if(origHistCenter->GetBinContent(i+1,j+1) <= -1000.0)
      {
        origHistCenter->SetBinContent(i+1,j+1,0.0);
      }
    }
  }

  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      for(Int_t k = 0; k < numZbins; k++)
      {
        if(medianHist->GetBinContent(i+1,j+1,k+1) <= -1000.0)
        {
          medianHist->SetBinContent(i+1,j+1,k+1,0.0);
        }
      }

      if(medianHistCenter->GetBinContent(i+1,j+1) <= -1000.0)
      {
        medianHistCenter->SetBinContent(i+1,j+1,0.0);
      }
    }
  }

  for(Int_t i = 0; i < numXbins; i++)
  {
    for(Int_t j = 0; j < numYbins; j++)
    {
      for(Int_t k = 0; k < numZbins; k++)
      {
        if(paddedHist->GetBinContent(i+1,j+1,k+1) <= -1000.0)
        {
          paddedHist->SetBinContent(i+1,j+1,k+1,0.0);
        }
      }

      if(paddedHistCenter->GetBinContent(i+1,j+1) <= -1000.0)
      {
        paddedHistCenter->SetBinContent(i+1,j+1,0.0);
      }
    }
  }

  outputfile.cd();
  origHist->Write();
  origHistCenter->Write();
  medianHist->Write();
  medianHistCenter->Write();
  paddedHist->Write();
  paddedHistCenter->Write();
  filledHist->Write();
  filledHistCenter->Write();
}

vector<TH3F*> getCalibHists()
{
  const Double_t xScale = 2.56/Lx;
  const Double_t yScale = 2.33/Ly;
  const Double_t zScale = 10.37/Lz;

  TH3F* RecoBkwd_Displacement_X = new TH3F("RecoBkwd_Displacement_X","",numDispDivisions_x+1,-1.0*WF_cathode/(2.0*((Double_t) numDispDivisions_x)),WF_cathode*(1.0+1.0/(2.0*((Double_t) numDispDivisions_x))),numDispDivisions_y+1,-1.0*(WF_top-WF_bottom)/(2.0*((Double_t) numDispDivisions_y))+WF_bottom,(WF_top-WF_bottom)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_y)))+WF_bottom,numDispDivisions_z+1,-1.0*(WF_downstream-WF_upstream)/(2.0*((Double_t) numDispDivisions_z))+WF_upstream,(WF_downstream-WF_upstream)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_z)))+WF_upstream);
  TH3F* RecoBkwd_Displacement_Y = new TH3F("RecoBkwd_Displacement_Y","",numDispDivisions_x+1,-1.0*WF_cathode/(2.0*((Double_t) numDispDivisions_x)),WF_cathode*(1.0+1.0/(2.0*((Double_t) numDispDivisions_x))),numDispDivisions_y+1,-1.0*(WF_top-WF_bottom)/(2.0*((Double_t) numDispDivisions_y))+WF_bottom,(WF_top-WF_bottom)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_y)))+WF_bottom,numDispDivisions_z+1,-1.0*(WF_downstream-WF_upstream)/(2.0*((Double_t) numDispDivisions_z))+WF_upstream,(WF_downstream-WF_upstream)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_z)))+WF_upstream);
  TH3F* RecoBkwd_Displacement_Z = new TH3F("RecoBkwd_Displacement_Z","",numDispDivisions_x+1,-1.0*WF_cathode/(2.0*((Double_t) numDispDivisions_x)),WF_cathode*(1.0+1.0/(2.0*((Double_t) numDispDivisions_x))),numDispDivisions_y+1,-1.0*(WF_top-WF_bottom)/(2.0*((Double_t) numDispDivisions_y))+WF_bottom,(WF_top-WF_bottom)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_y)))+WF_bottom,numDispDivisions_z+1,-1.0*(WF_downstream-WF_upstream)/(2.0*((Double_t) numDispDivisions_z))+WF_upstream,(WF_downstream-WF_upstream)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_z)))+WF_upstream);

  TFile* fileCalib = new TFile(inputFileName);
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
  TTreeReaderValue<Int_t> numEntries_calib(readerCalib, "numEntries.data_calib");

  while (readerCalib.Next())
  {
    if(*elecFate_calib == 1) {
      RecoBkwd_Displacement_X->Fill(doInvCoordTransformX(*xReco_calib),doInvCoordTransformY(*yReco_calib),doInvCoordTransformZ(*zReco_calib),-100.0*xScale*(*Dx_calib));
      RecoBkwd_Displacement_Y->Fill(doInvCoordTransformX(*xReco_calib),doInvCoordTransformY(*yReco_calib),doInvCoordTransformZ(*zReco_calib),100.0*yScale*(*Dy_calib));
      RecoBkwd_Displacement_Z->Fill(doInvCoordTransformX(*xReco_calib),doInvCoordTransformY(*yReco_calib),doInvCoordTransformZ(*zReco_calib),100.0*zScale*(*Dz_calib));
    }
    else {
      RecoBkwd_Displacement_X->Fill(doInvCoordTransformX(*xReco_calib),doInvCoordTransformY(*yReco_calib),doInvCoordTransformZ(*zReco_calib),-1000.0);
      RecoBkwd_Displacement_Y->Fill(doInvCoordTransformX(*xReco_calib),doInvCoordTransformY(*yReco_calib),doInvCoordTransformZ(*zReco_calib),-1000.0);
      RecoBkwd_Displacement_Z->Fill(doInvCoordTransformX(*xReco_calib),doInvCoordTransformY(*yReco_calib),doInvCoordTransformZ(*zReco_calib),-1000.0);
    }
  }

  vector<TH3F*> histVec;
  histVec.push_back(RecoBkwd_Displacement_X);
  histVec.push_back(RecoBkwd_Displacement_Y);
  histVec.push_back(RecoBkwd_Displacement_Z);

  return histVec;
}

vector<TH3F*> getSimInterpHists()
{
  TFile* simInterpFile = new TFile(simInterpFileName,"READ");

  TH3F* bkwd_Dx;
  TH3F* bkwd_Dy;
  TH3F* bkwd_Dz;

  if(isMC == true) {
    bkwd_Dx = (TH3F*) simInterpFile->Get("TrueBkwd_Displacement_X");
    bkwd_Dy = (TH3F*) simInterpFile->Get("TrueBkwd_Displacement_Y");
    bkwd_Dz = (TH3F*) simInterpFile->Get("TrueBkwd_Displacement_Z");
  }
  else {
    bkwd_Dx = (TH3F*) simInterpFile->Get("RecoBkwd_Displacement_X");
    bkwd_Dy = (TH3F*) simInterpFile->Get("RecoBkwd_Displacement_Y");
    bkwd_Dz = (TH3F*) simInterpFile->Get("RecoBkwd_Displacement_Z");
  }

  vector<TH3F*> histVec;
  histVec.push_back(bkwd_Dx);
  histVec.push_back(bkwd_Dy);
  histVec.push_back(bkwd_Dz);

  return histVec;
}

void conditionHist(TH3F* &fixHist, Int_t driftSign, Bool_t doCornerFix, Int_t component)
{
  TH1F* tempX = (TH1F*) fixHist->ProjectionX();
  TH1F* tempY = (TH1F*) fixHist->ProjectionY();
  TH1F* tempZ = (TH1F*) fixHist->ProjectionZ();

  // Find First Problematic Bin in X
  Int_t badXbin[200][200];
  for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
    for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
      badXbin[j-1][k-1] = 10000000;
      for(Int_t i = 1; i <= fixHist->GetNbinsX(); i++) {
        if((fixHist->GetBinContent(i,j,k) == -1000.0) && (i < badXbin[j-1][k-1])) {
          badXbin[j-1][k-1] = i;
	}
      }
    }
  }
    
  // Interpolate/Extrapolate in Y Direction
  for(Int_t i = 1; i <= fixHist->GetNbinsX(); i++) {
    for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
      Int_t count = 0;
      for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
        if(fixHist->GetBinContent(i,j,k) != -1000.0) {
          count++;
	}
      }

      if(count > 0.75*fixHist->GetNbinsY()) {
        std::vector<Double_t> f, Y;
        tk::spline s;

        for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
          if(fixHist->GetBinContent(i,j,k) != -1000.0) {
            Y.push_back(tempY->GetBinCenter(j));
	    f.push_back(fixHist->GetBinContent(i,j,k));
	  }
        }

        s.set_points(Y,f);

	Int_t firstbin = 1;
	Int_t lastbin = fixHist->GetNbinsY();
	Int_t flag = 0;
	for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
          if((flag == 0) && (fixHist->GetBinContent(i,j,k) == -1000.0)) {
	    firstbin++;
	  }
	  else if(flag == 0) {
            flag = 1;
	  }
	}
	flag = 0;
	for(Int_t j = fixHist->GetNbinsY(); j >= 1; j--) {
          if((flag == 0) && (fixHist->GetBinContent(i,j,k) == -1000.0)) {
	    lastbin--;
	  }
	  else if(flag == 0) {
            flag = 1;
	  }
	}
	
        for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
          if(fixHist->GetBinContent(i,j,k) == -1000.0) {
	    //if ((j >= firstbin) && (j <= lastbin)) {
              fixHist->SetBinContent(i,j,k,s(tempY->GetBinCenter(j)));
	    //}
	    //else if (j < firstbin) {
            //  fixHist->SetBinContent(i,j,k,s(tempY->GetBinCenter(firstbin)));
	    //}
	    //else {
            //  fixHist->SetBinContent(i,j,k,s(tempY->GetBinCenter(lastbin)));
	    //}
	  }
        }
      }
    }
  }

  // Interpolate/Extrapolate in Z Direction
  for(Int_t i = 1; i <= fixHist->GetNbinsX(); i++) {
    for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
      Int_t count = 0;
      for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
        if(fixHist->GetBinContent(i,j,k) != -1000.0) {
          count++;
	}
      }

      if(count > 0.75*fixHist->GetNbinsZ()) {
        std::vector<Double_t> f, Z;
        tk::spline s;

        for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
          if(fixHist->GetBinContent(i,j,k) != -1000.0) {
            Z.push_back(tempZ->GetBinCenter(k));
	    f.push_back(fixHist->GetBinContent(i,j,k));
	  }
        }
	
        s.set_points(Z,f);

	Int_t firstbin = 1;
	Int_t lastbin = fixHist->GetNbinsZ();
	Int_t flag = 0;
        for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
          if((flag == 0) && (fixHist->GetBinContent(i,j,k) == -1000.0)) {
	    firstbin++;
	  }
	  else if(flag == 0) {
            flag = 1;
	  }
	}
	flag = 0;
        for(Int_t k = fixHist->GetNbinsZ(); k >= 1; k--) {
          if((flag == 0) && (fixHist->GetBinContent(i,j,k) == -1000.0)) {
	    lastbin--;
	  }
	  else if(flag == 0) {
            flag = 1;
	  }
	}
	
        for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
          if(fixHist->GetBinContent(i,j,k) == -1000.0) {
	    //if ((k >= firstbin) && (k <= lastbin)) {
              fixHist->SetBinContent(i,j,k,s(tempZ->GetBinCenter(k)));
	    //}
	    //else if (k < firstbin) {
            //  fixHist->SetBinContent(i,j,k,s(tempZ->GetBinCenter(firstbin)));
	    //}
	    //else {
            //  fixHist->SetBinContent(i,j,k,s(tempZ->GetBinCenter(lastbin)));
	    //}
	  }
        }
      }
    }
  }

  // Interpolate/Extrapolate in X Direction
  for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
    for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
      Int_t count = 0;
      for(Int_t i = 1; i <= fixHist->GetNbinsX(); i++) {
        if(fixHist->GetBinContent(i,j,k) != -1000.0) {
          count++;
	}
      }

      std::vector<Double_t> f, X;
      tk::spline s;

      for(Int_t i = 1; i <= fixHist->GetNbinsX(); i++) {
        if(fixHist->GetBinContent(i,j,k) != -1000.0) {
          X.push_back(tempX->GetBinCenter(i));
          f.push_back(fixHist->GetBinContent(i,j,k));
        }
      }

      s.set_points(X,f);

      Int_t firstbin = 1;
      Int_t lastbin = fixHist->GetNbinsX();
      Int_t flag = 0;
      for(Int_t i = 1; i <= fixHist->GetNbinsX(); i++) {
        if((flag == 0) && (fixHist->GetBinContent(i,j,k) == -1000.0)) {
          firstbin++;
        }
        else if(flag == 0) {
          flag = 1;
        }
      }
      flag = 0;
      for(Int_t i = fixHist->GetNbinsX(); i >= 1; i--) {
        if((flag == 0) && (fixHist->GetBinContent(i,j,k) == -1000.0)) {
          lastbin--;
        }
        else if(flag == 0) {
          flag = 1;
        }
      }
      
      for(Int_t i = 1; i <= fixHist->GetNbinsX(); i++) {
        if(fixHist->GetBinContent(i,j,k) == -1000.0) {
          //if ((i >= firstbin) && (i <= lastbin)) {
            fixHist->SetBinContent(i,j,k,s(tempX->GetBinCenter(i)));
          //}
          //else if (i < firstbin) {
          //  fixHist->SetBinContent(i,j,k,s(tempX->GetBinCenter(firstbin)));
          //}
          //else {
          //  fixHist->SetBinContent(i,j,k,s(tempX->GetBinCenter(lastbin)));
          //}
        }
      }
    }
  }

  // Fix Corners
  if(doCornerFix == kTRUE) {
    for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
      for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
    
        if((j > 3) && (j < fixHist->GetNbinsY()-2) && (k > 3) && (k < fixHist->GetNbinsZ()-2)) {
          continue;
        }
    
        std::vector<Double_t> f, X;
        tk::spline s;
    
        if(driftSign < 0) {
          for(Int_t i = 1; i <= fixHist->GetNbinsX()-1; i++) {
            X.push_back(tempX->GetBinCenter(i));
            f.push_back(fixHist->GetBinContent(i,j,k));
          }
    
          s.set_points(X,f);
      	  fixHist->SetBinContent(fixHist->GetNbinsX(),j,k,s(tempX->GetBinCenter(fixHist->GetNbinsX())));
        }
        else {
          for(Int_t i = 2; i <= fixHist->GetNbinsX(); i++) {
            X.push_back(tempX->GetBinCenter(i));
            f.push_back(fixHist->GetBinContent(i,j,k));
          }
    
          s.set_points(X,f);
    	  fixHist->SetBinContent(1,j,k,s(tempX->GetBinCenter(1)));
        }      
      }
    }
  }

//  // Extrapolate Offsets in X Direction at Y/Z Edges
//  for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
//    for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
//      if((component == 1) && (j > 3) && (j < fixHist->GetNbinsY()-2) && (k > 3) && (k < fixHist->GetNbinsZ()-2)) {
//        continue;
//      }
//      if((component == 2) && (k > 3) && (k < fixHist->GetNbinsZ()-2)) {
//        continue;
//      }
//      if((component == 3) && (j > 3) && (j < fixHist->GetNbinsY()-2)) {
//        continue;
//      }
//
//      for(Int_t i = badXbin[j-1][k-1]; i <= fixHist->GetNbinsX(); i++) {
//	if(badXbin[j-1][k-1] > 1) {
//          fixHist->SetBinContent(i,j,k,fixHist->GetBinContent(badXbin[j-1][k-1]-1,j,k));
//	}
//      }
//    }
//  }
//  for(Int_t i = 1; i <= fixHist->GetNbinsX(); i++) {
//    for(Int_t j = 1; j <= fixHist->GetNbinsY(); j++) {
//      fixHist->SetBinContent(i,j,1,fixHist->GetBinContent(i,j,2));
//      fixHist->SetBinContent(i,j,fixHist->GetNbinsZ(),fixHist->GetBinContent(i,j,fixHist->GetNbinsZ()-1));
//    }
//  }
//  for(Int_t i = 1; i <= fixHist->GetNbinsX(); i++) {
//    for(Int_t k = 1; k <= fixHist->GetNbinsZ(); k++) {
//      fixHist->SetBinContent(i,1,k,fixHist->GetBinContent(i,2,k));
//      fixHist->SetBinContent(i,fixHist->GetNbinsY(),k,fixHist->GetBinContent(i,fixHist->GetNbinsY()-1,k));
//    }
//  }

  return;
}

vector<TH3F*> makeEfieldHists(const vector<TH3F*> &spatialHists, const Int_t driftSign)
{
  TH3F* Reco_ElecField_X = new TH3F("Reco_ElecField_X","",numDispDivisions_x+1,-1.0*WF_cathode/(2.0*((Double_t) numDispDivisions_x)),WF_cathode*(1.0+1.0/(2.0*((Double_t) numDispDivisions_x))),numDispDivisions_y+1,-1.0*(WF_top-WF_bottom)/(2.0*((Double_t) numDispDivisions_y))+WF_bottom,(WF_top-WF_bottom)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_y)))+WF_bottom,numDispDivisions_z+1,-1.0*(WF_downstream-WF_upstream)/(2.0*((Double_t) numDispDivisions_z))+WF_upstream,(WF_downstream-WF_upstream)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_z)))+WF_upstream);
  TH3F* Reco_ElecField_Y = new TH3F("Reco_ElecField_Y","",numDispDivisions_x+1,-1.0*WF_cathode/(2.0*((Double_t) numDispDivisions_x)),WF_cathode*(1.0+1.0/(2.0*((Double_t) numDispDivisions_x))),numDispDivisions_y+1,-1.0*(WF_top-WF_bottom)/(2.0*((Double_t) numDispDivisions_y))+WF_bottom,(WF_top-WF_bottom)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_y)))+WF_bottom,numDispDivisions_z+1,-1.0*(WF_downstream-WF_upstream)/(2.0*((Double_t) numDispDivisions_z))+WF_upstream,(WF_downstream-WF_upstream)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_z)))+WF_upstream);
  TH3F* Reco_ElecField_Z = new TH3F("Reco_ElecField_Z","",numDispDivisions_x+1,-1.0*WF_cathode/(2.0*((Double_t) numDispDivisions_x)),WF_cathode*(1.0+1.0/(2.0*((Double_t) numDispDivisions_x))),numDispDivisions_y+1,-1.0*(WF_top-WF_bottom)/(2.0*((Double_t) numDispDivisions_y))+WF_bottom,(WF_top-WF_bottom)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_y)))+WF_bottom,numDispDivisions_z+1,-1.0*(WF_downstream-WF_upstream)/(2.0*((Double_t) numDispDivisions_z))+WF_upstream,(WF_downstream-WF_upstream)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_z)))+WF_upstream);
  TH3F* Reco_ElecField_Mag = new TH3F("Reco_ElecField_Mag","",numDispDivisions_x+1,-1.0*WF_cathode/(2.0*((Double_t) numDispDivisions_x)),WF_cathode*(1.0+1.0/(2.0*((Double_t) numDispDivisions_x))),numDispDivisions_y+1,-1.0*(WF_top-WF_bottom)/(2.0*((Double_t) numDispDivisions_y))+WF_bottom,(WF_top-WF_bottom)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_y)))+WF_bottom,numDispDivisions_z+1,-1.0*(WF_downstream-WF_upstream)/(2.0*((Double_t) numDispDivisions_z))+WF_upstream,(WF_downstream-WF_upstream)*(1.0+1.0/(2.0*((Double_t) numDispDivisions_z)))+WF_upstream);
  
  Fade_3D interpEngine;
  map<Point3,Point3> fieldMap;

  TH1F* tempX_bkwd = (TH1F*) spatialHists[0]->ProjectionX();
  TH1F* tempY_bkwd = (TH1F*) spatialHists[0]->ProjectionY();
  TH1F* tempZ_bkwd = (TH1F*) spatialHists[0]->ProjectionZ();
  Double_t L = TMath::Abs(tempX_bkwd->GetBinCenter(2)-tempX_bkwd->GetBinCenter(1));

  Int_t NbinsX = spatialHists[0]->GetNbinsX();
  Int_t NbinsY = spatialHists[0]->GetNbinsY();
  Int_t NbinsZ = spatialHists[0]->GetNbinsZ();

  Int_t signFactor;
  if(driftSign > 0)
  {
    signFactor = -1;
  }
  else
  {
    signFactor = 1;
  }

  for(Int_t i = 1; i <= NbinsX-1; i++) {
    for(Int_t j = 1; j <= NbinsY; j++) {
      for(Int_t k = 1; k <= NbinsZ; k++) {
        Double_t Lnew = TMath::Sqrt(TMath::Power(L+(spatialHists[0]->GetBinContent(i+1,j,k)-spatialHists[0]->GetBinContent(i,j,k)),2)+TMath::Power(spatialHists[1]->GetBinContent(i+1,j,k)-spatialHists[1]->GetBinContent(i,j,k),2)+TMath::Power(spatialHists[2]->GetBinContent(i+1,j,k)-spatialHists[2]->GetBinContent(i,j,k),2));

	Double_t EfieldCalc = findElecField(Lnew/L);

	Double_t projX = (L+(spatialHists[0]->GetBinContent(i+1,j,k)-spatialHists[0]->GetBinContent(i,j,k)))/Lnew;
	Double_t projY = (spatialHists[1]->GetBinContent(i+1,j,k)-spatialHists[1]->GetBinContent(i,j,k))/Lnew;
	Double_t projZ = (spatialHists[2]->GetBinContent(i+1,j,k)-spatialHists[2]->GetBinContent(i,j,k))/Lnew;

        Double_t EfieldX = projX*EfieldCalc;
        Double_t EfieldY = projY*EfieldCalc;
        Double_t EfieldZ = projZ*EfieldCalc;

        Point3 p(((tempX_bkwd->GetBinCenter(i)+spatialHists[0]->GetBinContent(i,j,k))+(tempX_bkwd->GetBinCenter(i+1)+spatialHists[0]->GetBinContent(i+1,j,k)))/2.0,((tempY_bkwd->GetBinCenter(j)+spatialHists[1]->GetBinContent(i,j,k))+(tempY_bkwd->GetBinCenter(j)+spatialHists[1]->GetBinContent(i+1,j,k)))/2.0,((tempZ_bkwd->GetBinCenter(k)+spatialHists[2]->GetBinContent(i,j,k))+(tempZ_bkwd->GetBinCenter(k)+spatialHists[2]->GetBinContent(i+1,j,k)))/2.0);
	Point3 q((EfieldX-Efield)/Efield,signFactor*EfieldY/Efield,signFactor*EfieldZ/Efield);
	fieldMap[p] = q;
	interpEngine.insert(p);
      }
    }
  }

  TH1F* tempX_field = (TH1F*) Reco_ElecField_X->ProjectionX();
  TH1F* tempY_field = (TH1F*) Reco_ElecField_X->ProjectionY();
  TH1F* tempZ_field = (TH1F*) Reco_ElecField_X->ProjectionZ();
  for(Int_t i = 1; i <= Reco_ElecField_X->GetNbinsX(); i++) {
    for(Int_t j = 1; j <= Reco_ElecField_X->GetNbinsY(); j++) {
      for(Int_t k = 1; k <= Reco_ElecField_X->GetNbinsZ(); k++) {
        Point3 p(tempX_field->GetBinCenter(i),tempY_field->GetBinCenter(j),tempZ_field->GetBinCenter(k));
        Tet3 *tet = interpEngine.locate(p);

        if(tet != NULL) {
          Point3 *pa;
          Point3 *pb;
          Point3 *pc;
          Point3 *pd;

          tet->getCorners(pa,pb,pc,pd);
          pair<Point3,Point3> v1(*pa,fieldMap[*pa]);
          pair<Point3,Point3> v2(*pb,fieldMap[*pb]);
          pair<Point3,Point3> v3(*pc,fieldMap[*pc]);
          pair<Point3,Point3> v4(*pd,fieldMap[*pd]);

          Point3 result = getInterpVal(p,v1,v2,v3,v4);
          Reco_ElecField_X->SetBinContent(i,j,k,result.x());
          Reco_ElecField_Y->SetBinContent(i,j,k,result.y());
  	  Reco_ElecField_Z->SetBinContent(i,j,k,result.z());
        }
        else {
          Reco_ElecField_X->SetBinContent(i,j,k,-1000.0);
          Reco_ElecField_Y->SetBinContent(i,j,k,-1000.0);
          Reco_ElecField_Z->SetBinContent(i,j,k,-1000.0);
        }
      }
    }
  }

  // Condition Edge Regions
  conditionHist(Reco_ElecField_X,driftSign,kFALSE,1);
  conditionHist(Reco_ElecField_Y,driftSign,kFALSE,2);
  conditionHist(Reco_ElecField_Z,driftSign,kFALSE,3);

  // Fix Region Near Anode and Cathode for DeltaX
  for(Int_t j = 1; j <= Reco_ElecField_X->GetNbinsY(); j++) {
    for(Int_t k = 1; k <= Reco_ElecField_X->GetNbinsZ(); k++) {
  
      std::vector<Double_t> dX, X;
      tk::spline s;
  
      for(Int_t i = 2; i <= Reco_ElecField_X->GetNbinsX()-1; i++) {
        X.push_back(tempX_field->GetBinCenter(i));
        dX.push_back(Reco_ElecField_X->GetBinContent(i,j,k));
      }
  
      s.set_points(X,dX);
      Reco_ElecField_X->SetBinContent(1,j,k,s(tempX_field->GetBinCenter(1)));
      Reco_ElecField_X->SetBinContent(Reco_ElecField_X->GetNbinsX(),j,k,s(tempX_field->GetBinCenter(Reco_ElecField_X->GetNbinsX())));
    }
  }

  // Fix Region Near Cathode for DeltaY, DeltaZ
  for(Int_t j = 1; j <= Reco_ElecField_X->GetNbinsY(); j++) {
    for(Int_t k = 1; k <= Reco_ElecField_X->GetNbinsZ(); k++) {
  
      std::vector<Double_t> dY, dZ, X;
      tk::spline s_dY, s_dZ;

      if(driftSign < 0) {
        for(Int_t i = 1; i <= Reco_ElecField_X->GetNbinsX()-2; i++) {
          X.push_back(tempX_field->GetBinCenter(i));
          dY.push_back(Reco_ElecField_Y->GetBinContent(i,j,k));
          dZ.push_back(Reco_ElecField_Z->GetBinContent(i,j,k));
        }
    
        s_dY.set_points(X,dY);
        s_dZ.set_points(X,dZ);
      
        Reco_ElecField_Y->SetBinContent(Reco_ElecField_X->GetNbinsX()-1,j,k,s_dY(tempX_field->GetBinCenter(Reco_ElecField_X->GetNbinsX()-1)));
        Reco_ElecField_Y->SetBinContent(Reco_ElecField_X->GetNbinsX(),j,k,s_dY(tempX_field->GetBinCenter(Reco_ElecField_X->GetNbinsX())));
        Reco_ElecField_Z->SetBinContent(Reco_ElecField_X->GetNbinsX()-1,j,k,s_dZ(tempX_field->GetBinCenter(Reco_ElecField_X->GetNbinsX()-1)));
        Reco_ElecField_Z->SetBinContent(Reco_ElecField_X->GetNbinsX(),j,k,s_dZ(tempX_field->GetBinCenter(Reco_ElecField_X->GetNbinsX())));
      }
      else {
        for(Int_t i = 3; i <= Reco_ElecField_X->GetNbinsX(); i++) {
          X.push_back(tempX_field->GetBinCenter(i));
          dY.push_back(Reco_ElecField_Y->GetBinContent(i,j,k));
          dZ.push_back(Reco_ElecField_Z->GetBinContent(i,j,k));
        }
    
        s_dY.set_points(X,dY);
        s_dZ.set_points(X,dZ);
      
        Reco_ElecField_Y->SetBinContent(1,j,k,s_dY(tempX_field->GetBinCenter(1)));
        Reco_ElecField_Y->SetBinContent(2,j,k,s_dY(tempX_field->GetBinCenter(2)));
        Reco_ElecField_Z->SetBinContent(1,j,k,s_dZ(tempX_field->GetBinCenter(1)));
        Reco_ElecField_Z->SetBinContent(2,j,k,s_dZ(tempX_field->GetBinCenter(2)));
      }      
    }
  }
  //
  // NOTE:  For above, can try simply averaging last and third to last bin to get second to last bin (keep values from going negative)
  //

  // Apply Median Filter
  doMedianFiltering(Reco_ElecField_X);
  doMedianFiltering(Reco_ElecField_Y);
  doMedianFiltering(Reco_ElecField_Z);

  // Fill E Field Magnitude Histogram
  for(Int_t i = 1; i <= Reco_ElecField_Mag->GetNbinsX(); i++) {
    for(Int_t j = 1; j <= Reco_ElecField_Mag->GetNbinsY(); j++) {
      for(Int_t k = 1; k <= Reco_ElecField_Mag->GetNbinsZ(); k++) {
        Reco_ElecField_Mag->SetBinContent(i,j,k,TMath::Sqrt(TMath::Power(1.0+Reco_ElecField_X->GetBinContent(i,j,k),2)+TMath::Power(Reco_ElecField_Y->GetBinContent(i,j,k),2)+TMath::Power(Reco_ElecField_Z->GetBinContent(i,j,k),2))-1.0);

	if(Reco_ElecField_X->GetBinContent(i,j,k) == 0.0) {
	  Reco_ElecField_X->SetBinContent(i,j,k,0.000001);
	}
        if(Reco_ElecField_Y->GetBinContent(i,j,k) == 0.0) {
	  Reco_ElecField_Y->SetBinContent(i,j,k,0.000001);
	}
        if(Reco_ElecField_Z->GetBinContent(i,j,k) == 0.0) {
	  Reco_ElecField_Z->SetBinContent(i,j,k,0.000001);
	}
        if(Reco_ElecField_Mag->GetBinContent(i,j,k) == 0.0) {
	  Reco_ElecField_Mag->SetBinContent(i,j,k,0.000001);
	}
      }
    }
  }

  vector<TH3F*> resultHists;
  resultHists.push_back(Reco_ElecField_X);
  resultHists.push_back(Reco_ElecField_Y);
  resultHists.push_back(Reco_ElecField_Z);
  resultHists.push_back(Reco_ElecField_Mag);
  
  return resultHists;
}

void doMedianFiltering(TH3F* &inputHist)
{
  TH3F *medianHist = (TH3F*) inputHist->Clone();

  for(Int_t i = 1; i <= inputHist->GetNbinsX(); i++) {
    for(Int_t j = 1; j <= inputHist->GetNbinsY(); j++) {
      for(Int_t k = 1; k <= inputHist->GetNbinsZ(); k++) {
        Double_t vals[27];
        Double_t weights[27];
      
        Int_t numCounts = 0;
        for(Int_t p = -1; p <= 1; p++) {
          for(Int_t q = -1; q <= 1; q++) {
            for(Int_t r = -1; r <= 1; r++)
            {
              if((i+p < 1) || (i+p > inputHist->GetNbinsX()) || (j+q < 1) || (j+q > inputHist->GetNbinsY()) || (k+r < 1) || (k+r > inputHist->GetNbinsZ())) {
                vals[numCounts] = 0.0;
                weights[numCounts] = 0.0;
	      }
	      else {
                vals[numCounts] = inputHist->GetBinContent(i+p,j+q,k+r);
                weights[numCounts] = 1.0;
	      }

      	      numCounts++;
	    }
          }
        }

        medianHist->SetBinContent(i,j,k,TMath::Median(numCounts,&vals[0],&weights[0]));
      }
    }
  }

  for(Int_t i = 1; i <= inputHist->GetNbinsX(); i++) {
    for(Int_t j = 1; j <= inputHist->GetNbinsY(); j++) {
      for(Int_t k = 1; k <= inputHist->GetNbinsZ(); k++) {
        inputHist->SetBinContent(i,j,k,medianHist->GetBinContent(i,j,k));
      }
    }
  }

  return;
}

Double_t findElecField(Double_t input)
{
  // MicroBooNE uses own fit (along with ICARUS data, taken at 89 K) for v(E) curve
  
  Double_t minE = 0.200;
  Double_t maxE = 0.400;
  Int_t numIter = 20;
  
  TF1 driftVelFit("driftVelFit","pol5",-0.05,1.1);
  driftVelFit.SetParameter(0,0.0);
  driftVelFit.SetParameter(1,5.52605);
  driftVelFit.SetParameter(2,-6.29745);
  driftVelFit.SetParameter(3,2.29095);
  driftVelFit.SetParameter(4,1.58989);
  driftVelFit.SetParameter(5,-1.06544);

  Double_t v0 = driftVelFit.Eval(Efield);
  Double_t minDiff = 9999999.0;
  Double_t Eval = Efield;
  
  for(Int_t h = 0; h <= numIter; h++)
  {
    Double_t Etest = minE+(((Double_t) h)/((Double_t) numIter))*(maxE-minE);
    Double_t newDiff = fabs((driftVelFit.Eval(Etest)/v0) - input);
    if(newDiff < minDiff)
    {
      minDiff = newDiff;
      Eval = Etest;
    }
  }

  Double_t minE2 = Eval-(maxE-minE)/(2.0*((Double_t) numIter));
  Double_t maxE2 = Eval+(maxE-minE)/(2.0*((Double_t) numIter));
  for(Int_t h = 0; h <= numIter; h++)
  {
    Double_t Etest = minE2+(((Double_t) h)/((Double_t) numIter))*(maxE2-minE2);
    Double_t newDiff = fabs((driftVelFit.Eval(Etest)/v0) - input);
    if(newDiff < minDiff)
    {
      minDiff = newDiff;
      Eval = Etest;
    }
  }

  Double_t minE3 = Eval-(maxE-minE)/(2.0*TMath::Power(((Double_t) numIter),2));
  Double_t maxE3 = Eval+(maxE-minE)/(2.0*TMath::Power(((Double_t) numIter),2));
  for(Int_t h = 0; h <= numIter; h++)
  {
    Double_t Etest = minE3+(((Double_t) h)/((Double_t) numIter))*(maxE3-minE3);
    Double_t newDiff = fabs((driftVelFit.Eval(Etest)/v0) - input);
    if(newDiff < minDiff)
    {
      minDiff = newDiff;
      Eval = Etest;
    }
  }

  return Eval;
}

Point3 getInterpVal(Point3 point, pair<Point3,Point3> vtx1, pair<Point3,Point3> vtx2, pair<Point3,Point3> vtx3, pair<Point3,Point3> vtx4)
{
  double a11 = vtx1.first.x()-vtx4.first.x();
  double a21 = vtx1.first.y()-vtx4.first.y();
  double a31 = vtx1.first.z()-vtx4.first.z();

  double a12 = vtx2.first.x()-vtx4.first.x();
  double a22 = vtx2.first.y()-vtx4.first.y();
  double a32 = vtx2.first.z()-vtx4.first.z();

  double a13 = vtx3.first.x()-vtx4.first.x();
  double a23 = vtx3.first.y()-vtx4.first.y();
  double a33 = vtx3.first.z()-vtx4.first.z();

  double det_a = a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31);

  double b11 = (a22*a33-a23*a32)/det_a;
  double b21 = (a23*a31-a21*a33)/det_a;
  double b31 = (a21*a32-a22*a31)/det_a;

  double b12 = (a13*a32-a12*a33)/det_a;
  double b22 = (a11*a33-a13*a31)/det_a;
  double b32 = (a12*a31-a11*a32)/det_a;

  double b13 = (a12*a23-a13*a22)/det_a;
  double b23 = (a13*a21-a11*a23)/det_a;
  double b33 = (a11*a22-a12*a21)/det_a;

  double lam1 = b11*(point.x()-vtx4.first.x())+b12*(point.y()-vtx4.first.y())+b13*(point.z()-vtx4.first.z());
  double lam2 = b21*(point.x()-vtx4.first.x())+b22*(point.y()-vtx4.first.y())+b23*(point.z()-vtx4.first.z());
  double lam3 = b31*(point.x()-vtx4.first.x())+b32*(point.y()-vtx4.first.y())+b33*(point.z()-vtx4.first.z());
  double lam4 = 1-lam1-lam2-lam3;

  Point3 f1 = vtx1.second;
  Point3 f2 = vtx2.second;
  Point3 f3 = vtx3.second;
  Point3 f4 = vtx4.second;

  Point3 result(lam1*f1.x()+lam2*f2.x()+lam3*f3.x()+lam4*f4.x(),lam1*f1.y()+lam2*f2.y()+lam3*f3.y()+lam4*f4.y(),lam1*f1.z()+lam2*f2.z()+lam3*f3.z()+lam4*f4.z());

  return result;
}
