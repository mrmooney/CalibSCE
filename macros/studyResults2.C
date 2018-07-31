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
#include "TTreeReaderArray.h"

#include "/usr/include/eigen3/Eigen/Dense"

const Int_t minTrackPoints = 50;
const Int_t numTrackSegPoints = 15;

const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;
//const Double_t Lx = 3.6
//const Double_t Ly = 6.0;
//const Double_t Lz = 7.2;

const Int_t numDivisions = 25;

Double_t corr_Dx[26][26][101];
Double_t corr_Dy[26][26][101];
Double_t corr_Dz[26][26][101];

struct Point {
  float x;
  float y;
  float z;
};

struct PCAResults {
  TVector3 centroid;
  pair<TVector3,TVector3> endPoints;
  float length;
  TVector3 eVals;
  vector<TVector3> eVecs;
};

typedef vector<Point> PointCloud;

PCAResults DoPCA(const PointCloud &points);
float LinInterp(float x, float x1, float x2, float q00, float q01);
float TrilinInterp(float x, float y, float z, float q000, float q001, float q010, float q011, float q100, float q101, float q110, float q111, float x1, float x2, float y1, float y2, float z1, float z2);

void studyResults2()
{
  const Double_t bufferLength = 0.05;
  //const Double_t bufferLength = 0.0;

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

  TH1F *origAngHist = new TH1F("origAngHist","",50,0.0,10.0);
  TH1F *corrAngHist = new TH1F("corrAngHist","",50,0.0,10.0);
  
  //TFile* inputFileInterp = new TFile("output_interp.root");
  //TFile* inputFileInterp = new TFile("output_MC.root");
  TFile* inputFileInterp = new TFile("output_MC_30k.root");
  //TFile* inputFileInterp = new TFile("output_data.root");
  //TFile* inputFileInterp = new TFile("output_data_new.root");
  //TFile* inputFileInterp = new TFile("output_data_30k.root");
  //TFile* inputFileInterp = new TFile("../results/Results100kTracks_MC_Dec22_2017_withInterp.root");
  //TFile* inputFileInterp = new TFile("../results/Results100kTracks_MC_Dec22_2017.root");
  //TFile* inputFileInterp = new TFile("../results/Results100kTracks_MC_Jan4_2018.root");
  TTreeReader readerCalib("SpaCEtree_calib", inputFileInterp);
  //TTreeReaderValue<Double_t> Dx_calib(readerCalib, "Dx");
  //TTreeReaderValue<Double_t> Dy_calib(readerCalib, "Dy");
  //TTreeReaderValue<Double_t> Dz_calib(readerCalib, "Dz");
  //TTreeReaderValue<Double_t> xReco_calib(readerCalib, "x_reco");
  //TTreeReaderValue<Double_t> yReco_calib(readerCalib, "y_reco");
  //TTreeReaderValue<Double_t> zReco_calib(readerCalib, "z_reco");
  //TTreeReaderValue<Double_t> xTrue_calib(readerCalib, "x_true");
  //TTreeReaderValue<Double_t> yTrue_calib(readerCalib, "y_true");
  //TTreeReaderValue<Double_t> zTrue_calib(readerCalib, "z_true");
  //TTreeReaderValue<Int_t> elecFate_calib(readerCalib, "elecFate");
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

  //TFile* inputFile = new TFile("output.root");
  TFile* inputFile = new TFile("output_MC_tracks.root");
  //TFile* inputFile = new TFile("output_MC_NoSCE_tracks.root");
  //TFile* inputFile = new TFile("output_data_tracks.root");
  //TFile* inputFile = new TFile("output_data_tracks_new.root");
  //TFile* inputFile = new TFile("output_data_tracks_noMCScut.root");
  TTreeReader readerTracks("SpaCEtree_tracks", inputFile);
  TTreeReaderValue<Int_t> nElec_tracks(readerTracks, "nElec_tracks");
  TTreeReaderArray<Double_t> elecX_tracks(readerTracks, "elecX_tracks");
  TTreeReaderArray<Double_t> elecY_tracks(readerTracks, "elecY_tracks");
  TTreeReaderArray<Double_t> elecZ_tracks(readerTracks, "elecZ_tracks");

  for (int x = 0; x <= numDivisions_x; x++ ) {
    for (int y = 0; y <= numDivisions_y; y++ ) {
      for (int z = 0; z <= numDivisions_z; z++ ) {
        corr_Dx[x][y][z] = 0.0;
        corr_Dy[x][y][z] = 0.0;
        corr_Dz[x][y][z] = 0.0;
      }
    }
  }

  while (readerCalib.Next())
  {
    corr_Dx[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += *Dx_calib;
    corr_Dy[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += *Dy_calib;
    corr_Dz[TMath::Nint(*xReco_calib*numDivisions_x/Lx)][TMath::Nint(*yReco_calib*numDivisions_y/Ly)][TMath::Nint(*zReco_calib*numDivisions_z/Lz)] += *Dz_calib;
  }

  int counter = 0;
  //while (readerTracks.Next())
  while ((readerTracks.Next()) && (counter < 100000))
  {
    if(counter % 1000 == 0) {cout << counter << endl;}
    counter++;
    
    if (*nElec_tracks >= minTrackPoints) {

      PointCloud startPoints;
      PointCloud endPoints;

      Int_t numBadPoints_start = 0;
      for (int i = 0; i < numTrackSegPoints; i++) {
	if ((elecX_tracks[i] <= 0.0-bufferLength) || (elecX_tracks[i] >= Lx+bufferLength) || (elecY_tracks[i] <= 0.0-bufferLength) || (elecY_tracks[i] >= Ly+bufferLength) || (elecZ_tracks[i] <= 0.0-bufferLength) || (elecZ_tracks[i] >= Lz+bufferLength)) {
          numBadPoints_start++;
	  continue;
	}
	
	Point tempPoint;
        tempPoint.x = elecX_tracks[i];
        tempPoint.y = elecY_tracks[i];
        tempPoint.z = elecZ_tracks[i];

        startPoints.push_back(tempPoint);
      }
      if (numBadPoints_start > 5) continue;

      Int_t numBadPoints_end = 0;
      for (int i = *nElec_tracks-1; i > *nElec_tracks-1-numTrackSegPoints; i--) {
        if ((elecX_tracks[i] <= 0.0-bufferLength) || (elecX_tracks[i] >= Lx+bufferLength) || (elecY_tracks[i] <= 0.0-bufferLength) || (elecY_tracks[i] >= Ly+bufferLength) || (elecZ_tracks[i] <= 0.0-bufferLength) || (elecZ_tracks[i] >= Lz+bufferLength)) {
          numBadPoints_end++;
	  continue;
	}

	Point tempPoint;
        tempPoint.x = elecX_tracks[i];
        tempPoint.y = elecY_tracks[i];
        tempPoint.z = elecZ_tracks[i];

        endPoints.push_back(tempPoint);
      }
      if (numBadPoints_end > 5) continue;

      PCAResults results_start = DoPCA(startPoints);
      PCAResults results_end = DoPCA(endPoints);

      Double_t dotProd = results_start.eVecs[0](0)*results_end.eVecs[0](0) + results_start.eVecs[0](1)*results_end.eVecs[0](1) + results_start.eVecs[0](2)*results_end.eVecs[0](2);
      Double_t startMag = sqrt(pow(results_start.eVecs[0](0),2) + pow(results_start.eVecs[0](1),2) + pow(results_start.eVecs[0](2),2));
      Double_t endMag = sqrt(pow(results_end.eVecs[0](0),2) + pow(results_end.eVecs[0](1),2) + pow(results_end.eVecs[0](2),2));

      Double_t dTheta = TMath::ACos(dotProd/(startMag*endMag));

      origAngHist->Fill(min(180.0*dTheta/3.14159265,180.0-(180.0*dTheta/3.14159265)));

      PointCloud startPointsCorr;
      PointCloud endPointsCorr;

      Int_t numBadPoints_start_corr = 0;
      for (int i = 0; i < numTrackSegPoints; i++) {
        if ((elecX_tracks[i] <= 0.0) || (elecX_tracks[i] >= Lx) || (elecY_tracks[i] <= 0.0) || (elecY_tracks[i] >= Ly) || (elecZ_tracks[i] <= 0.0) || (elecZ_tracks[i] >= Lz)) {
          numBadPoints_start_corr++;
	  continue;
	}

	/*
	cout << counter << " " << i << "   ";
	cout << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " " 
	     << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
             << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " " 
	     << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
	     << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " " 
	     << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
             << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " " 
	     << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << endl;
	*/
	  
	Point tempPoint;
        tempPoint.x = elecX_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
        tempPoint.y = elecY_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
	tempPoint.z = elecZ_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));

        startPointsCorr.push_back(tempPoint);
      }
      if (numBadPoints_start_corr > 5) continue;

      Int_t numBadPoints_end_corr = 0;
      for (int i = *nElec_tracks-1; i > *nElec_tracks-1-numTrackSegPoints; i--) {
        if ((elecX_tracks[i] <= 0.0) || (elecX_tracks[i] >= Lx) || (elecY_tracks[i] <= 0.0) || (elecY_tracks[i] >= Ly) || (elecZ_tracks[i] <= 0.0) || (elecZ_tracks[i] >= Lz)) {
	  numBadPoints_end_corr++;
	  continue;
	}

	/*
	cout << counter << " " << i << "   ";
	cout << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " " 
	     << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
             << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " " 
	     << corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
	     << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " " 
	     << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << " "
             << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)] << " " 
	     << corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)] << endl;
	*/
	
	Point tempPoint;
        tempPoint.x = elecX_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dx[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
        tempPoint.y = elecY_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dy[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));
	tempPoint.z = elecZ_tracks[i] + TrilinInterp(elecX_tracks[i],elecY_tracks[i],elecZ_tracks[i],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)],corr_Dz[(Int_t)TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)][(Int_t)TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)][(Int_t)TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)],TMath::Floor(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Ceil(elecX_tracks[i]*numDivisions_x/Lx)*(Lx/numDivisions_x),TMath::Floor(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Ceil(elecY_tracks[i]*numDivisions_y/Ly)*(Ly/numDivisions_y),TMath::Floor(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z),TMath::Ceil(elecZ_tracks[i]*numDivisions_z/Lz)*(Lz/numDivisions_z));

        endPointsCorr.push_back(tempPoint);
      }
      if (numBadPoints_end_corr > 5) continue;
      
      PCAResults results_start_corr = DoPCA(startPointsCorr);
      PCAResults results_end_corr = DoPCA(endPointsCorr);

      Double_t dotProdCorr = results_start_corr.eVecs[0](0)*results_end_corr.eVecs[0](0) + results_start_corr.eVecs[0](1)*results_end_corr.eVecs[0](1) + results_start_corr.eVecs[0](2)*results_end_corr.eVecs[0](2);
      Double_t startMagCorr = sqrt(pow(results_start_corr.eVecs[0](0),2) + pow(results_start_corr.eVecs[0](1),2) + pow(results_start_corr.eVecs[0](2),2));
      Double_t endMagCorr = sqrt(pow(results_end_corr.eVecs[0](0),2) + pow(results_end_corr.eVecs[0](1),2) + pow(results_end_corr.eVecs[0](2),2));

      Double_t dThetaCorr = TMath::ACos(dotProdCorr/(startMagCorr*endMagCorr));

      //cout << "    ANGLE:  " << min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265));

      corrAngHist->Fill(min(180.0*dThetaCorr/3.14159265,180.0-(180.0*dThetaCorr/3.14159265)));
    }
  }
  origAngHist->Scale(1.0/origAngHist->Integral());
  corrAngHist->Scale(1.0/corrAngHist->Integral());
  
  gStyle->SetTitleW(0.9);
  gStyle->SetOptStat(0);

  TCanvas *c_orig = new TCanvas();
  c_orig->cd();
  origAngHist->GetXaxis()->SetTitle("#Delta#theta [degrees]");
  origAngHist->GetXaxis()->SetTitleOffset(0.95);
  origAngHist->GetXaxis()->SetTitleSize(0.045);
  origAngHist->GetYaxis()->SetTitle("Arb. Units");
  origAngHist->GetYaxis()->SetTitleOffset(0.95);
  origAngHist->GetYaxis()->SetTitleSize(0.05);
  origAngHist->GetYaxis()->SetNoExponent(kTRUE);
  origAngHist->SetStats(0);
  origAngHist->SetLineWidth(2.0);
  origAngHist->SetLineColor(kRed);
  origAngHist->Draw("HIST");
  origAngHist->Draw("AXISsame");
  origAngHist->SetMinimum(0.0001);
  c_orig->SaveAs("origAngHist.png");

  TCanvas *c_corr = new TCanvas();
  c_corr->cd();
  corrAngHist->GetXaxis()->SetTitle("#Delta#theta [degrees]");
  corrAngHist->GetXaxis()->SetTitleOffset(0.95);
  corrAngHist->GetXaxis()->SetTitleSize(0.045);
  corrAngHist->GetYaxis()->SetTitle("Arb. Units");
  corrAngHist->GetYaxis()->SetTitleOffset(0.95);
  corrAngHist->GetYaxis()->SetTitleSize(0.05);
  corrAngHist->GetYaxis()->SetNoExponent(kTRUE);
  corrAngHist->SetStats(0);
  corrAngHist->SetLineWidth(2.0);
  corrAngHist->SetLineColor(kGreen+2);
  corrAngHist->Draw("HIST");
  corrAngHist->Draw("AXISsame");
  corrAngHist->SetMinimum(0.0001);
  c_corr->SaveAs("corrAngHist.png");
  
  TCanvas *c_combined = new TCanvas();
  c_combined->cd();
  TLegend *leg_combined = new TLegend(0.55,0.70,0.88,0.85);
  leg_combined->SetLineColor(kWhite);
  leg_combined->AddEntry(origAngHist,"Before Calibration","L");
  leg_combined->AddEntry(corrAngHist,"After Calibration","L");
  origAngHist->GetXaxis()->SetTitle("#Delta#theta [degrees]");
  origAngHist->Draw("HIST");
  corrAngHist->Draw("HISTsame");
  leg_combined->Draw("same");
  origAngHist->Draw("AXISsame");
  origAngHist->SetMaximum(1.1*max(origAngHist->GetMaximum(),corrAngHist->GetMaximum()));
  origAngHist->SetMinimum(0.0001);
  c_combined->SaveAs("combinedAngHist.png");
}

PCAResults DoPCA(const PointCloud &points) {

  TVector3 outputCentroid;
  pair<TVector3,TVector3> outputEndPoints;
  float outputLength;
  TVector3 outputEigenValues;
  vector<TVector3> outputEigenVecs;

  float meanPosition[3] = {0., 0., 0.};
  unsigned int nThreeDHits = 0;

  for (unsigned int i = 0; i < points.size(); i++) {
    meanPosition[0] += points[i].x;
    meanPosition[1] += points[i].y;
    meanPosition[2] += points[i].z;
    ++nThreeDHits;
  }

  if (nThreeDHits == 0) {
    PCAResults results;
    return results; // FAIL FROM NO INPUT POINTS
  }

  const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
  meanPosition[0] /= nThreeDHitsAsFloat;
  meanPosition[1] /= nThreeDHitsAsFloat;
  meanPosition[2] /= nThreeDHitsAsFloat;
  outputCentroid = TVector3(meanPosition[0], meanPosition[1], meanPosition[2]);

  // Define elements of our covariance matrix
  float xi2 = 0.0;
  float xiyi = 0.0;
  float xizi = 0.0;
  float yi2 = 0.0;
  float yizi = 0.0;
  float zi2 = 0.0;
  float weightSum = 0.0;

  for (unsigned int i = 0; i < points.size(); i++) {
      const float weight(1.);
      const float x((points[i].x - meanPosition[0]) * weight);
      const float y((points[i].y - meanPosition[1]) * weight);
      const float z((points[i].z - meanPosition[2]) * weight);

      xi2  += x * x;
      xiyi += x * y;
      xizi += x * z;
      yi2  += y * y;
      yizi += y * z;
      zi2  += z * z;
      weightSum += weight * weight;
  }

  // Using Eigen package
  Eigen::Matrix3f sig;

  sig << xi2, xiyi, xizi,
         xiyi, yi2, yizi,
         xizi, yizi, zi2;

  //if (std::fabs(weightSum) < std::numeric_limits<float>::epsilon())
  //{
  //    std::cout << "PCAShowerParticleBuildingAlgorithm::RunPCA - The total weight of three dimensional hits is 0!" << std::endl;
  //    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
  //}

  sig *= 1.0 / weightSum;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

  //if (eigenMat.info() != Eigen::ComputationInfo::Success)
  //{
  //    std::cout << "PCAShowerParticleBuildingAlgorithm::RunPCA - PCA decompose failure, number of three D hits = " << nThreeDHits << std::endl;
  //    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
  //}

  typedef std::pair<float,size_t> EigenValColPair;
  typedef std::vector<EigenValColPair> EigenValColVector;

  EigenValColVector eigenValColVector;
  const auto &resultEigenMat(eigenMat.eigenvalues());
  eigenValColVector.emplace_back(resultEigenMat(0), 0);
  eigenValColVector.emplace_back(resultEigenMat(1), 1);
  eigenValColVector.emplace_back(resultEigenMat(2), 2);

  std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;} );

  outputEigenValues = TVector3(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);

  const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

  for (const EigenValColPair &pair : eigenValColVector) {
     outputEigenVecs.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
  }

  PCAResults results;

  Eigen::ParametrizedLine<float,3> priAxis(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)),Eigen::Vector3f(outputEigenVecs[0](0),outputEigenVecs[0](1),outputEigenVecs[0](2)));

  Eigen::Vector3f endPoint1(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));
  Eigen::Vector3f endPoint2(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));

  Eigen::Vector3f testPoint;
  Eigen::Vector3f projTestPoint;
  float maxDist1 = -1.0;
  float maxDist2 = -1.0;
  float dist;
  float dotP;
  for (unsigned int i = 0; i < points.size(); i++) {
    testPoint = Eigen::Vector3f(points[i].x,points[i].y,points[i].z);
    projTestPoint = priAxis.projection(testPoint);
    dist = sqrt(pow(projTestPoint(0)-outputCentroid(0),2.0)+pow(projTestPoint(1)-outputCentroid(1),2.0)+pow(projTestPoint(2)-outputCentroid(2),2.0));
    dotP = (projTestPoint(0)-outputCentroid(0))*outputEigenVecs[0](0) + (projTestPoint(1)-outputCentroid(1))*outputEigenVecs[0](1) + (projTestPoint(2)-outputCentroid(2))*outputEigenVecs[0](2);

    if ((dotP < 0.0) && (dist > maxDist1)) {
      endPoint1 = projTestPoint;
      maxDist1 = dist;
    }
    else if ((dotP > 0.0) && (dist > maxDist2)) {
      endPoint2 = projTestPoint;
      maxDist2 = dist;
    }
  }

  outputEndPoints.first = TVector3(endPoint1(0),endPoint1(1),endPoint1(2));
  outputEndPoints.second = TVector3(endPoint2(0),endPoint2(1),endPoint2(2));

  outputLength = sqrt(pow(endPoint2(0)-endPoint1(0),2.0)+pow(endPoint2(1)-endPoint1(1),2.0)+pow(endPoint2(2)-endPoint1(2),2.0));

  results.centroid = outputCentroid;
  results.endPoints = outputEndPoints;
  results.length = outputLength;
  results.eVals = outputEigenValues;
  results.eVecs = outputEigenVecs;

  return results;
}

float LinInterp(float x, float x1, float x2, float q00, float q01) {
  return ((x2 - x) / (x2 - x1)) * q00 + ((x - x1) / (x2 - x1)) * q01;
}

float TrilinInterp(float x, float y, float z, float q000, float q001, float q010, float q011, float q100, float q101, float q110, float q111, float x1, float x2, float y1, float y2, float z1, float z2) {
  float x00 = LinInterp(x, x1, x2, q000, q100);
  float x10 = LinInterp(x, x1, x2, q010, q110);
  float x01 = LinInterp(x, x1, x2, q001, q101);
  float x11 = LinInterp(x, x1, x2, q011, q111);
  float r0 = LinInterp(y, y1, y2, x00, x01);
  float r1 = LinInterp(y, y1, y2, x10, x11);

  return LinInterp(z, z1, z2, r0, r1);
}
