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
#include <TPrincipal.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TVirtualFFT.h>
#include <TSystem.h>
#include <TGraph2D.h>

#include "normVec.hpp"
#include "SpaceChargeMicroBooNE.h"
#include "spline.h"

#include "/usr/include/eigen3/Eigen/Dense"

using namespace std;

//const Char_t *inputFileLaser = "data/laserDataSCE.root";
const Char_t *inputFileLaser = "data/laserDataSCE_NEW.root";
//const Char_t *inputFileLaser = "data/laserScan_data.root";
//const Char_t *inputFileLaser = "data/laserScan_MC.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_small.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_small_NEW.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_small_NEW_withMCS.root";
//const Char_t *inputFileCosmic = "data/mcsample_small_combined.root";
//const Char_t *inputFileCosmic = "data/oldMCsample_50000events.root";
//const Char_t *inputFileCosmic = "data/newMCsample_2Mevents.root";
const Char_t *inputFileCosmic = "data/MC_Cosmics.root";
//const Char_t *inputFileCosmic = "data/MC_Cosmics_NoSCE.root";
//const Char_t *inputFileCosmic = "data/Data_Run1_EXTBNB.root";
//const Char_t *inputFileCosmic = "data/first_set_of_EXTBNB_production.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_ProtoDUNESP.root";
//const Char_t *inputFileCosmic = "data/cosmicDataSCE_ProtoDUNESP_withMCS.root";

TFile* outputFile = new TFile("output.root","RECREATE");

const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;
//const Double_t Lx = 3.6;
//const Double_t Ly = 6.0;
//const Double_t Lz = 7.2;

const Double_t ScaleFactorX = Lx/2.56;
const Double_t ScaleFactorY = Ly/2.33;
const Double_t ScaleFactorZ = Lz/10.37;

const Double_t TrueAnode = Lx;
const Double_t TrueCathode = Lx*(2.56-2.548)/2.56;
//const Double_t TrueCathode = 0.0; // DATA
const Double_t TrueTop = Ly*(1.174+1.165)/(2.33);
const Double_t TrueBottom = Ly*(-1.154+1.165)/(2.33);
const Double_t TrueUpstream = Lz*(0.001)/(10.37);
const Double_t TrueDownstream = Lz*(10.369)/(10.37);
const Double_t ShiftedAnode = Lx*(2.56-0.0006)/2.56;
//const Double_t ShiftedAnode = Lx*(2.56-(-0.0056))/2.56; // DATA
const Double_t ShiftedCathode = Lx*(2.56-2.5524)/2.56;
//const Double_t ShiftedCathode = Lx*(2.56-2.5818)/2.56; // DATA
const Double_t OffsetCathode = 0.004;

const Bool_t isMC = true;
const Bool_t isSCEon = true;

const Double_t relAngleCut = 20.0;
const Double_t maxXdist = 0.05;
const Double_t maxYdist = 0.20;
//const Double_t maxYdist = 0.05;
const Double_t maxZdist = 0.20;
//const Double_t maxZdist = 0.05;

Int_t minInputTrackNum = 0;
//Int_t maxInputTrackNum = 1000000000;
Int_t maxInputTrackNum = 1150000;

Int_t nCalibDivisions = 25; // MICROBOONE
//Int_t nCalibDivisions = 18; // PROTODUNE-SP

const Double_t piVal = 3.14159265;

Int_t nCalibDivisions_x;
Int_t nCalibDivisions_y;
Int_t nCalibDivisions_z;

vector<Double_t> calibWeight[101][101][401];
vector<Double_t> calibDeltaX[101][101][401];
vector<Double_t> calibDeltaY[101][101][401];
vector<Double_t> calibDeltaZ[101][101][401];

Double_t trueDeltaX[101][101][401];
Double_t trueDeltaY[101][101][401];
Double_t trueDeltaZ[101][101][401];

Double_t trueFwdDeltaX[101][101][401];
Double_t trueFwdDeltaY[101][101][401];
Double_t trueFwdDeltaZ[101][101][401];

vector<vector<Double_t> > calibTopDeltaX;
vector<vector<Double_t> > calibTopDeltaY;
vector<vector<Double_t> > calibTopDeltaZ;

vector<vector<Double_t> > calibBottomDeltaX;
vector<vector<Double_t> > calibBottomDeltaY;
vector<vector<Double_t> > calibBottomDeltaZ;

vector<vector<Double_t> > calibDownstreamDeltaX;
vector<vector<Double_t> > calibDownstreamDeltaY;
vector<vector<Double_t> > calibDownstreamDeltaZ;

vector<vector<Double_t> > calibUpstreamDeltaX;
vector<vector<Double_t> > calibUpstreamDeltaY;
vector<vector<Double_t> > calibUpstreamDeltaZ;

vector<vector<Double_t> > calibCathodeDeltaX;
vector<vector<Double_t> > calibCathodeDeltaY;
vector<vector<Double_t> > calibCathodeDeltaZ;

TH3F *bulkCalibHistDeltaX;
TH3F *bulkCalibHistDeltaY;
TH3F *bulkCalibHistDeltaZ;

TH2F *faceCalibHistTopDeltaX;
TH2F *faceCalibHistTopDeltaY;
TH2F *faceCalibHistTopDeltaZ;

TH2F *faceCalibHistBottomDeltaX;
TH2F *faceCalibHistBottomDeltaY;
TH2F *faceCalibHistBottomDeltaZ;

TH2F *faceCalibHistUpstreamDeltaX;
TH2F *faceCalibHistUpstreamDeltaY;
TH2F *faceCalibHistUpstreamDeltaZ;

TH2F *faceCalibHistDownstreamDeltaX;
TH2F *faceCalibHistDownstreamDeltaY;
TH2F *faceCalibHistDownstreamDeltaZ;

TH2F *faceCalibHistCathodeDeltaX;
TH2F *faceCalibHistCathodeDeltaY;
TH2F *faceCalibHistCathodeDeltaZ;

TH2F *faceSimHistTopDeltaX;
TH2F *faceSimHistTopDeltaY;
TH2F *faceSimHistTopDeltaZ;

TH2F *faceSimHistBottomDeltaX;
TH2F *faceSimHistBottomDeltaY;
TH2F *faceSimHistBottomDeltaZ;

TH2F *faceSimHistUpstreamDeltaX;
TH2F *faceSimHistUpstreamDeltaY;
TH2F *faceSimHistUpstreamDeltaZ;

TH2F *faceSimHistDownstreamDeltaX;
TH2F *faceSimHistDownstreamDeltaY;
TH2F *faceSimHistDownstreamDeltaZ;

TH2F *faceSimHistCathodeDeltaX;
TH2F *faceSimHistCathodeDeltaY;
TH2F *faceSimHistCathodeDeltaZ;

struct elecInfo
{
  Double_t x;
  Double_t y;
  Double_t z;
  Double_t t;
  Double_t x_mod;
  Double_t y_mod;
  Double_t z_mod;
  Double_t t_mod;
  Int_t fate;
};

struct trackInfo
{
  Int_t pdgID;
  Double_t energy;
  Double_t x0;
  Double_t y0;
  Double_t z0;
  Double_t x1;
  Double_t y1;
  Double_t z1;
  Double_t theta;
  Double_t phi;
  vector<elecInfo> electrons;
};

struct calibTrackInfo
{
  Double_t x0_calib;
  Double_t y0_calib;
  Double_t z0_calib;
  Double_t x1_calib;
  Double_t y1_calib;
  Double_t z1_calib;
  Double_t theta_calib;
  Double_t phi_calib;
  vector<Double_t> DxVec;
  vector<Double_t> DyVec;
  vector<Double_t> DzVec;
  Bool_t calibFlag;
  trackInfo track;
};

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

Double_t doCoordTransformX(const Double_t inputX);
Double_t doCoordTransformY(const Double_t inputY);
Double_t doCoordTransformZ(const Double_t inputZ);
Double_t doInvCoordTransformX(const Double_t inputX);
Double_t doInvCoordTransformY(const Double_t inputY);
Double_t doInvCoordTransformZ(const Double_t inputZ);
vector<Double_t> getParabolaParameters(const vector<elecInfo> &parabola_points_track);
vector<Double_t> findClosestPOA(const calibTrackInfo &trackA, const calibTrackInfo &trackB);
vector<Double_t> findDistortedClosestPOA(const calibTrackInfo &trackA, const calibTrackInfo &trackB);
void getLArSoftTrackSet(vector<trackInfo> &tracks, Int_t inputType, Int_t calibMode, Int_t maxCosmicTracks, Double_t minTrackMCS_anode, Double_t minTrackMCS_cathode, Double_t minTrackMCS_crossing);
vector<trackInfo> getTrackSet(Int_t mapType);
vector<calibTrackInfo> makeCalibTracks(const vector<trackInfo> &tracks);
void doLaserLaserCalib(const vector<calibTrackInfo> &laserCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t saveInfo);
void doLaserCosmicCalib(const vector<calibTrackInfo> &laserCalibTracks, const vector<calibTrackInfo> &cosmicCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t saveInfo);
void doCosmicCosmicCalib(const vector<calibTrackInfo> &laserCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t saveInfo);
void doCalibration(const vector<trackInfo> &laserTracks, const vector<trackInfo> &cosmicTracks, Double_t distScale, Double_t maxDistFactor, Int_t numIterations, Int_t saveInfo);
void saveTrackInfo(const vector<trackInfo> &tracks);
void loadTruthMap();
Double_t getTruthOffset(Double_t xVal, Double_t yVal, Double_t zVal, int comp);
vector<Double_t> getTruthOffsets(Double_t xVal, Double_t yVal, Double_t zVal);
void loadTruthFwdMap();
Double_t getTruthFwdOffset(Double_t xVal, Double_t yVal, Double_t zVal, int comp);
vector<Double_t> getTruthFwdOffsets(Double_t xVal, Double_t yVal, Double_t zVal);
Double_t findCathodeOffset(Double_t yVal, Double_t zVal);
void doCalibFaces(const vector<trackInfo> &cosmicTracks, Int_t minTrackPoints, Int_t numTrackSegPoints, Double_t minWeight, Int_t minEntries, Int_t minMaxBin, Int_t saveInfo);
PCAResults DoPCA(const PointCloud &points);
Double_t findHistMedian(TH1F *hist);
void conditionFaceMap(vector<vector<Double_t> > &inputMap, Int_t faceNum);
void extrapolate(TH2F *interpHist, Int_t edgeXbins, Int_t edgeYbins, Int_t faceNum);
Double_t getBulkCorr(Double_t xVal, Double_t yVal, Double_t zVal, Int_t comp);
Double_t getFaceCorr(Double_t xVal, Double_t yVal, Double_t zVal, Int_t faceNum, Int_t comp);

Int_t main(Int_t argc, Char_t** argv)
{
  TStopwatch timer;
  timer.Start();

  gErrorIgnoreLevel = kError;
    
  if(argc > 1) {
    minInputTrackNum = atoi(argv[1]);
    maxInputTrackNum = atoi(argv[2]);
  }
  
  nCalibDivisions_x = nCalibDivisions;
  nCalibDivisions_y = TMath::Nint((Ly/Lx)*((Double_t)nCalibDivisions));
  nCalibDivisions_z = TMath::Nint((Lz/Lx)*((Double_t)nCalibDivisions));

  Double_t minTrackMCS_anode;
  Double_t minTrackMCS_cathode;
  Double_t minTrackMCS_crossing;
  if (isMC == true) {
    minTrackMCS_anode = 3.3;
    minTrackMCS_cathode = 1.7;
    minTrackMCS_crossing = 1.4;
  }
  else {
    minTrackMCS_anode = 1.5;
    minTrackMCS_cathode = 1.1;
    minTrackMCS_crossing = 0.0;
  }

  loadTruthMap();
  loadTruthFwdMap();

  calibTopDeltaX.resize(nCalibDivisions_x+1);
  calibTopDeltaY.resize(nCalibDivisions_x+1);
  calibTopDeltaZ.resize(nCalibDivisions_x+1);
  calibBottomDeltaX.resize(nCalibDivisions_x+1);
  calibBottomDeltaY.resize(nCalibDivisions_x+1);
  calibBottomDeltaZ.resize(nCalibDivisions_x+1);
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    calibTopDeltaX[x].resize(nCalibDivisions_z+1);
    calibTopDeltaY[x].resize(nCalibDivisions_z+1);
    calibTopDeltaZ[x].resize(nCalibDivisions_z+1);
    calibBottomDeltaX[x].resize(nCalibDivisions_z+1);
    calibBottomDeltaY[x].resize(nCalibDivisions_z+1);
    calibBottomDeltaZ[x].resize(nCalibDivisions_z+1);
  }
  calibUpstreamDeltaX.resize(nCalibDivisions_x+1);
  calibUpstreamDeltaY.resize(nCalibDivisions_x+1);
  calibUpstreamDeltaZ.resize(nCalibDivisions_x+1);
  calibDownstreamDeltaX.resize(nCalibDivisions_x+1);
  calibDownstreamDeltaY.resize(nCalibDivisions_x+1);
  calibDownstreamDeltaZ.resize(nCalibDivisions_x+1);
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    calibUpstreamDeltaX[x].resize(nCalibDivisions_y+1);
    calibUpstreamDeltaY[x].resize(nCalibDivisions_y+1);
    calibUpstreamDeltaZ[x].resize(nCalibDivisions_y+1);
    calibDownstreamDeltaX[x].resize(nCalibDivisions_y+1);
    calibDownstreamDeltaY[x].resize(nCalibDivisions_y+1);
    calibDownstreamDeltaZ[x].resize(nCalibDivisions_y+1);
  }
  calibCathodeDeltaX.resize(nCalibDivisions_y+1);
  calibCathodeDeltaY.resize(nCalibDivisions_y+1);
  calibCathodeDeltaZ.resize(nCalibDivisions_y+1);
  for(Int_t y = 0; y <= nCalibDivisions_y; y++)
  {
    calibCathodeDeltaX[y].resize(nCalibDivisions_z+1);
    calibCathodeDeltaY[y].resize(nCalibDivisions_z+1);
    calibCathodeDeltaZ[y].resize(nCalibDivisions_z+1);
  }

  bulkCalibHistDeltaX = new TH3F("bulkCalibHistDeltaX","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  bulkCalibHistDeltaY = new TH3F("bulkCalibHistDeltaY","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  bulkCalibHistDeltaZ = new TH3F("bulkCalibHistDeltaZ","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  
  faceCalibHistTopDeltaX = new TH2F("faceCalibHistTopDeltaX","",nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)),nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)));
  faceCalibHistTopDeltaY = new TH2F("faceCalibHistTopDeltaY","",nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)),nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)));
  faceCalibHistTopDeltaZ = new TH2F("faceCalibHistTopDeltaZ","",nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)),nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)));
  
  faceCalibHistBottomDeltaX = new TH2F("faceCalibHistBottomDeltaX","",nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)),nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)));
  faceCalibHistBottomDeltaY = new TH2F("faceCalibHistBottomDeltaY","",nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)),nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)));
  faceCalibHistBottomDeltaZ = new TH2F("faceCalibHistBottomDeltaZ","",nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)),nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)));

  faceCalibHistUpstreamDeltaX = new TH2F("faceCalibHistUpstreamDeltaX","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)));
  faceCalibHistUpstreamDeltaY = new TH2F("faceCalibHistUpstreamDeltaY","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)));
  faceCalibHistUpstreamDeltaZ = new TH2F("faceCalibHistUpstreamDeltaZ","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)));

  faceCalibHistDownstreamDeltaX = new TH2F("faceCalibHistDownstreamDeltaX","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)));
  faceCalibHistDownstreamDeltaY = new TH2F("faceCalibHistDownstreamDeltaY","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)));
  faceCalibHistDownstreamDeltaZ = new TH2F("faceCalibHistDownstreamDeltaZ","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)));

  faceCalibHistCathodeDeltaX = new TH2F("faceCalibHistCathodeDeltaX","",nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)));
  faceCalibHistCathodeDeltaY = new TH2F("faceCalibHistCathodeDeltaY","",nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)));
  faceCalibHistCathodeDeltaZ = new TH2F("faceCalibHistCathodeDeltaZ","",nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)));
  
  //vector<trackInfo> laserTracks = getTrackSet(1);
  //vector<trackInfo> cosmicTracks = getTrackSet(2);

  //////////////////////////////////////////////
  /// MAIN PART OF CODE (CHANGE THESE THINGS)
  //////////////////////////////////////////////
  
  vector<trackInfo> laserTracks;
  getLArSoftTrackSet(laserTracks,1,-1,-1,-1,-1,-1);

  vector<trackInfo> cosmicTracks;
  
  //////////getLArSoftTrackSet(cosmicTracks,2,1,10000,minTrackMCS_anode,minTrackMCS_cathode,minTrackMCS_crossing);
  //getLArSoftTrackSet(cosmicTracks,2,1,10000,minTrackMCS_anode,10000000.0,10000000.0);
  //doCalibration(laserTracks,cosmicTracks,0.01,3,1,0);
  //
  ////getLArSoftTrackSet(cosmicTracks,2,1,1000000,0.0,0.0,0.0);
  ////doCalibFaces(cosmicTracks,50,15,0.25,4,5,0);
  //
  ////getLArSoftTrackSet(cosmicTracks,2,1,10000,100000000.0,10000000.0,0.0);
  ////doCalibration(laserTracks,cosmicTracks,0.01,3,1,0);
  //
  //getLArSoftTrackSet(cosmicTracks,2,1,1000000,0.0,0.0,0.0);
  //doCalibFaces(cosmicTracks,50,15,0.25,4,5,1);
  //
  //getLArSoftTrackSet(cosmicTracks,2,3,10000,minTrackMCS_anode,minTrackMCS_cathode,minTrackMCS_crossing);
  //doCalibration(laserTracks,cosmicTracks,0.01,3,1,1);

  getLArSoftTrackSet(cosmicTracks,2,4,10000,minTrackMCS_anode,minTrackMCS_cathode,minTrackMCS_crossing);
  doCalibration(laserTracks,cosmicTracks,0.01,3,1,1);
  
  //saveTrackInfo(cosmicTracks);
  
  timer.Stop();
  cout << "Calibration Time:  " << timer.CpuTime() << " sec." << endl;
  
  outputFile->Write();
  outputFile->Close();
  
  return 0;
}

Double_t doCoordTransformX(const Double_t inputX)
{
  Double_t outputX;
  outputX = Lx - (Lx/2.56)*inputX/100.0;

  return outputX;
}

Double_t doCoordTransformY(const Double_t inputY)
{
  Double_t outputY;
  outputY = (Ly/2.33)*(inputY+116.5)/100.0;

  return outputY;
}

Double_t doCoordTransformZ(const Double_t inputZ)
{
  Double_t outputZ;
  outputZ = (Lz/10.37)*(inputZ)/100.0;

  return outputZ;
}

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

vector<Double_t> getParabolaParameters(const vector<elecInfo> &parabola_points_track)
{
  if (parabola_points_track.size() < 3) cout << "Less than three points provided for the parameters of parabola." << endl;
 
  Double_t x_middle = parabola_points_track.at(1).x_mod;
  Double_t y_middle = parabola_points_track.at(1).y_mod;
  Double_t z_middle = parabola_points_track.at(1).z_mod;
 
  //first/linear transformation
  Double_t x_0 = parabola_points_track.at(0).x_mod-x_middle;
  Double_t x_1 = parabola_points_track.at(1).x_mod-x_middle;
  Double_t x_2 = parabola_points_track.at(2).x_mod-x_middle;
 
  Double_t y_0 = parabola_points_track.at(0).y_mod-y_middle;
  Double_t y_1 = parabola_points_track.at(1).y_mod-y_middle;
  Double_t y_2 = parabola_points_track.at(2).y_mod-y_middle;
 
  Double_t z_0 = parabola_points_track.at(0).z_mod;
  Double_t z_1 = parabola_points_track.at(1).z_mod;
  Double_t z_2 = parabola_points_track.at(2).z_mod;
 
  //angle of rotation where y_0_2 == y_2_2
  Double_t phi = atan((y_2-y_0)/(x_0-x_2));
 
  //preform the second transformation
  Double_t x_0_2 = x_0*cos(phi)-y_0*sin(phi);
  Double_t x_1_2 = x_1*cos(phi)-y_1*sin(phi);
  Double_t x_2_2 = x_2*cos(phi)-y_2*sin(phi);
 
  Double_t y_0_2 = x_0*sin(phi)+y_0*cos(phi);
  Double_t y_1_2 = x_1*sin(phi)+y_1*cos(phi);
  Double_t y_2_2 = x_2*sin(phi)+y_2*cos(phi);
 
  //since x_1_2 = 0 and y_1_2 = 0, because that's the middle point, c = 0
  Double_t a = ((y_0_2*x_2_2)-(y_2_2*x_0_2))/(x_0_2*x_2_2*(x_0_2-x_2_2));
  Double_t b = (y_0_2-(a*pow(x_0_2,2)))/x_0_2;
 
  //find the plane defined by the three points without any transformations: z = d*x + e*y + f
  //any one of the three points dotted with the normal vector will give the plane: P_1*n = 0
 
  //comps. for a vector connecting point 2 to point 1
  Double_t p1p2_x = parabola_points_track.at(1).x_mod-parabola_points_track.at(0).x_mod;
  Double_t p1p2_y = parabola_points_track.at(1).y_mod-parabola_points_track.at(0).y_mod;
  Double_t p1p2_z = parabola_points_track.at(1).z_mod-parabola_points_track.at(0).z_mod;
 
  //comps. for a vector connecting point 3 to point 1
  Double_t p1p3_x = parabola_points_track.at(2).x_mod-parabola_points_track.at(0).x_mod;
  Double_t p1p3_y = parabola_points_track.at(2).y_mod-parabola_points_track.at(0).y_mod;
  Double_t p1p3_z = parabola_points_track.at(2).z_mod-parabola_points_track.at(0).z_mod;
 
  // normal vector, n, = P1P2xP1P3
  Double_t norm_x = (p1p2_y*p1p3_z)-(p1p3_y*p1p2_z);
  Double_t norm_y = (p1p3_x*p1p2_z)-(p1p2_x*p1p3_z);
  Double_t norm_z = (p1p2_x*p1p3_y)-(p1p3_x*p1p2_y);
 
  // z = d*x + e*y + f
  Double_t d = -norm_x/norm_z;
  Double_t e = -norm_y/norm_z;
  Double_t f = ((parabola_points_track.at(0).x_mod*norm_x)+(parabola_points_track.at(0).y_mod*norm_y)+(parabola_points_track.at(0).z_mod*norm_z))/norm_z;
 
  vector<Double_t> return_vector;
  return_vector.push_back(a);
  return_vector.push_back(b);
  return_vector.push_back(d);
  return_vector.push_back(e);
  return_vector.push_back(f);
  return_vector.push_back(phi);
  return_vector.push_back(x_middle);
  return_vector.push_back(y_middle);
  return return_vector;
}

vector<Double_t> findClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB)
{ 
  //   (xA,yA,zA)+t(xA_step,yA_step,zA_step)
  Double_t xA = calibTrackA.x0_calib, yA = calibTrackA.y0_calib, zA = calibTrackA.z0_calib;
  Double_t xB = calibTrackB.x0_calib, yB = calibTrackB.y0_calib, zB = calibTrackB.z0_calib;

  Double_t xA_step = -1.*sin(calibTrackA.theta_calib)*sin(calibTrackA.phi_calib);
  Double_t yA_step = cos(calibTrackA.theta_calib);
  Double_t zA_step = sin(calibTrackA.theta_calib)*cos(calibTrackA.phi_calib);

  Double_t xB_step = -1.*sin(calibTrackB.theta_calib)*sin(calibTrackB.phi_calib);
  Double_t yB_step = cos(calibTrackB.theta_calib);
  Double_t zB_step = sin(calibTrackB.theta_calib)*cos(calibTrackB.phi_calib);

  //perpendicular line between the two tracks
  Double_t x_prep = (yA_step*zB_step)-(yB_step*zA_step);
  Double_t y_prep = (xB_step*zA_step)-(xA_step*zB_step);
  Double_t z_prep = (xA_step*yB_step)-(xB_step*yA_step);
 
  // if cross product is zero then the lines are parallel so return distance = -2
  if (x_prep == 0 && y_prep == 0 && z_prep == 0) {
    vector<Double_t> return_vector;
    return_vector.push_back(-2.);
    return_vector.push_back((xA+xB)/2.);
    return_vector.push_back((yA+yB)/2.);
    return_vector.push_back((zA+zB)/2.);   
    return return_vector;
  }
  //normalize the perpendicular line
  Double_t mag_prep = sqrt(pow(x_prep,2)+pow(y_prep,2)+pow(z_prep,2));
  Double_t x_prep_norm = x_prep / mag_prep;
  Double_t y_prep_norm = y_prep / mag_prep;
  Double_t z_prep_norm = z_prep / mag_prep;
 
  //defined to make the math simplier
  Double_t a = y_prep_norm*(xA-xB);
  Double_t b = x_prep_norm*(yA-yB);
  Double_t c = z_prep_norm*(xA-xB);
  Double_t d = x_prep_norm*(zA-zB);
 
  Double_t g = y_prep_norm*xA_step;
  Double_t h = y_prep_norm*xB_step;
  Double_t i = x_prep_norm*yA_step;
  Double_t j = x_prep_norm*yB_step;
  Double_t k = z_prep_norm*xA_step;
  Double_t l = z_prep_norm*xB_step;
  Double_t m = x_prep_norm*zA_step;
  Double_t n = x_prep_norm*zB_step;
 
  Double_t chi = (l-n)/(h-j);
 
  //alpha: "t" for the first line //beta:  "t" for the second line
  Double_t alpha = (chi*(a-b)-c+d)/(k-m-(chi*(g-i)));
  Double_t beta = (c-d+alpha*(k-m))/(l-n);

  Double_t cpa_xA = xA+alpha*xA_step;
  Double_t cpa_yA = yA+alpha*yA_step;
  Double_t cpa_zA = zA+alpha*zA_step;

  Double_t cpa_xB = xB+beta*xB_step;
  Double_t cpa_yB = yB+beta*yB_step;
  Double_t cpa_zB = zB+beta*zB_step;
 
  //distance between the closest points on the lines\tracks
  //Double_t distance = sqrt(pow((xA+alpha*xA_step)-(xB+beta*xB_step),2)+pow((yA+alpha*yA_step)-(yB+beta*yB_step),2)+pow((zA+alpha*zA_step)-(zB+beta*zB_step),2));
  Double_t distance = sqrt(pow(cpa_xA-cpa_xB,2)+pow(cpa_yA-cpa_yB,2)+pow(cpa_zA-cpa_zB,2));

  //midpoint between points of closest approach on both lines
  Double_t x_mid = (cpa_xA+cpa_xB)/2.;
  Double_t y_mid = (cpa_yA+cpa_yB)/2.;
  Double_t z_mid = (cpa_zA+cpa_zB)/2.;
 
  //check to see if this point is outside the detector
  //  TODO: < or <=, Is 0 or Lx "outside"?
  //if (x_mid > Lx || x_mid < 0 || y_mid > Ly || y_mid < 0 || z_mid > Lz || z_mid < 0) distance = -1;
  if (cpa_xA > Lx || cpa_xA < 0 || cpa_yA > Ly || cpa_yA < 0 || cpa_zA > Lz || cpa_zA < 0 ||
      cpa_xB > Lx || cpa_xB < 0 || cpa_yB > Ly || cpa_yB < 0 || cpa_zB > Lz || cpa_zB < 0) distance = -1;
 
  vector<Double_t> return_vector;

  return_vector.push_back(distance);
  return_vector.push_back(x_mid);
  return_vector.push_back(y_mid);
  return_vector.push_back(z_mid);

  ////// NEW 12/5/2017 //////
  return_vector.push_back(cpa_xB-cpa_xA);
  return_vector.push_back(cpa_yB-cpa_yA);
  return_vector.push_back(cpa_zB-cpa_zA);
  ///////////////////////////
  
  return return_vector;
}

vector<Double_t> findDistortedClosestPOA(const calibTrackInfo &calibTrackA, const calibTrackInfo &calibTrackB)
{
  vector<Double_t> return_vector;
 
  if (calibTrackA.track.electrons.size() < 3) {
    cout << "Less than three track points are inside the detector for calibTrackA." << endl;
    return return_vector;
  }
  if (calibTrackB.track.electrons.size() < 3) {
    cout << "Less than three track points are inside the detector for calibTrackB." << endl;
    return return_vector;
  }

  Double_t min_distance = -1.;//-1 means no value has yet been set; for the first iteration
  Double_t min_a, min_b;
 
  for(int iter_a = 0; iter_a < calibTrackA.track.electrons.size(); iter_a++)
    for(int iter_b = 0; iter_b < calibTrackB.track.electrons.size(); iter_b++) {
      Double_t distance = sqrt(pow(calibTrackA.track.electrons.at(iter_a).x_mod-calibTrackB.track.electrons.at(iter_b).x_mod,2)+pow(calibTrackA.track.electrons.at(iter_a).y_mod-calibTrackB.track.electrons.at(iter_b).y_mod,2)+pow(calibTrackA.track.electrons.at(iter_a).z_mod-calibTrackB.track.electrons.at(iter_b).z_mod,2));
      if (min_distance == -1 || distance < min_distance) {
	min_distance = distance;
	min_a = iter_a; min_b = iter_b;
      }
    }
    
  //TODO if min_distance is zero then no need to fit a parabola

  vector<elecInfo> parabola_points_calibTrackA, parabola_points_calibTrackB;
 
  if (min_a != 0) parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a-1));
  else parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a+2));
  parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a));
  if (min_a != calibTrackA.track.electrons.size() - 1) parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a+1));
  else parabola_points_calibTrackA.push_back(calibTrackA.track.electrons.at(min_a-2));
 
  if (min_b != 0) parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b-1));
  else parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b+2));
  parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b));
  if (min_b != calibTrackB.track.electrons.size()-1) parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b+1));
  else parabola_points_calibTrackB.push_back(calibTrackB.track.electrons.at(min_b-2));
 
  vector<Double_t> parabolaParameters_calibTrackA = getParabolaParameters(parabola_points_calibTrackA);
  vector<Double_t> parabolaParameters_calibTrackB = getParabolaParameters(parabola_points_calibTrackB);
 
  Double_t aA = parabolaParameters_calibTrackA.at(0);
  Double_t bA = parabolaParameters_calibTrackA.at(1);
  Double_t dA = parabolaParameters_calibTrackA.at(2);
  Double_t eA = parabolaParameters_calibTrackA.at(3);
  Double_t fA = parabolaParameters_calibTrackA.at(4);
  Double_t phiA = parabolaParameters_calibTrackA.at(5);
  Double_t x_midA = parabolaParameters_calibTrackA.at(6);
  Double_t y_midA = parabolaParameters_calibTrackA.at(7);
 
  Double_t aB = parabolaParameters_calibTrackB.at(0);
  Double_t bB = parabolaParameters_calibTrackB.at(1);
  Double_t dB = parabolaParameters_calibTrackB.at(2);
  Double_t eB = parabolaParameters_calibTrackB.at(3);
  Double_t fB = parabolaParameters_calibTrackB.at(4);
  Double_t phiB = parabolaParameters_calibTrackB.at(5);
  Double_t x_midB = parabolaParameters_calibTrackB.at(6);
  Double_t y_midB = parabolaParameters_calibTrackB.at(7);
 
  Double_t stepSize_x2 = 2.5*pow(10,-4); //0.25mm
 
  Double_t min_distance_parabola = -1.;
 
  Double_t min_Ax = 0., min_Ay = 0., min_Az = 0., min_Bx = 0., min_By = 0., min_Bz = 0.;
  if (min_a != 0 && min_a != calibTrackA.track.electrons.size()-1 && min_b != 0 && min_b != calibTrackB.track.electrons.size()-1) {
    // x2 denotes x" //_0 is the first point and _2 is the third point //_1 is the middle point which is zero for x2 (x")
    Double_t A_x2_0 = (calibTrackA.track.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
    Double_t A_x2_2 = (calibTrackA.track.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
   
    Double_t B_x2_0 = (calibTrackB.track.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
    Double_t B_x2_2 = (calibTrackB.track.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);
   
    //starting and ending points for the loop //values determined below
    Double_t A_x2_start = 0., A_x2_end = 0., B_x2_start = 0., B_x2_end = 0.;
    if (A_x2_0 < A_x2_2) {
      A_x2_start = A_x2_0;
      A_x2_end = A_x2_2;
    } else {
      A_x2_start = A_x2_2;
      A_x2_end = A_x2_0;
    }
    if (B_x2_0 < B_x2_2) {
      B_x2_start = B_x2_0;
      B_x2_end = B_x2_2;
    } else {
      B_x2_start = B_x2_2;
      B_x2_end = B_x2_0;
    }
   
    for (Double_t A_x2 = A_x2_start; A_x2 < A_x2_end;  A_x2 += stepSize_x2) {
      Double_t A_y2 = aA*pow(A_x2,2)+bA*A_x2;
      //get x,y,z
      Double_t A_y = A_y2*cos(phiA) - A_x2*sin(phiA) + y_midA;
      Double_t A_x = (A_y2 - A_y*cos(phiA) + y_midA*cos(phiA) + x_midA*sin(phiA)) / sin(phiA);
      Double_t A_z = dA*A_x+eA*A_y+fA;
      for (Double_t B_x2 = B_x2_start; B_x2 < B_x2_end; B_x2 += stepSize_x2) {
	Double_t B_y2 = aB*pow(B_x2,2)+bB*B_x2;
	Double_t B_y = B_y2*cos(phiB) - B_x2*sin(phiB) + y_midB;
	Double_t B_x = (B_y2 - B_y*cos(phiB) + y_midB*cos(phiB) + x_midB*sin(phiB)) / sin(phiB);
	Double_t B_z = dB*B_x+eB*B_y+fB;
   
	Double_t distance = sqrt(pow(B_x-A_x,2)+pow(B_y-A_y,2)+pow(B_z-A_z,2));
	if (min_distance_parabola == -1 || distance < min_distance_parabola) {
	  min_distance_parabola = distance;
     
	  min_Ax = A_x; min_Ay = A_y; min_Az = A_z;
	  min_Bx = B_x; min_By = B_y; min_Bz = B_z;
	}
      }
    }
  }
  else {
    Double_t stepSize_x2_A = stepSize_x2;
    Double_t stepSize_x2_B = stepSize_x2;
   
    // x2 denotes x" //_0 is the first point and _2 is the third point //_1 is the middle point which is zero for x2 (x")
    Double_t A_x2_0 = min_a == 0 ? -1. : (calibTrackA.track.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
    Double_t A_x2_2 = min_a == calibTrackA.track.electrons.size()-1 ? -1. : (calibTrackA.track.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
   
    Double_t B_x2_0 = min_b == 0 ? -1. : (calibTrackB.track.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
    Double_t B_x2_2 = min_b == calibTrackB.track.electrons.size()-1 ? -1. : (calibTrackB.track.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);

    Int_t nearBoundary_A = 1, nearBoundary_B = 1;
   
    Double_t A_x2_start = 0., A_x2_end = 0., B_x2_start = 0., B_x2_end = 0.;
   
    //If the first point is set to be -1 then start the loop form the final point and move towards the first point until the detector boundary is reached.
    if (A_x2_0 == -1.) A_x2_start = A_x2_2;
    else if (A_x2_2 == -1.) A_x2_start = A_x2_0;
    else {
      nearBoundary_A = 0;
     
      A_x2_0 = (calibTrackA.track.electrons.at(min_a-1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a-1).y_mod-y_midA)*sin(phiA);
      A_x2_2 = (calibTrackA.track.electrons.at(min_a+1).x_mod-x_midA)*cos(phiA)-(calibTrackA.track.electrons.at(min_a+1).y_mod-y_midA)*sin(phiA);
      if (A_x2_0 < A_x2_2) {
	A_x2_start = A_x2_0;
	A_x2_end = A_x2_2;
      } else {
	A_x2_start = A_x2_2;
	A_x2_end = A_x2_0;
      }
    }
    if (B_x2_0 == -1.) B_x2_start = B_x2_2;
    else if (B_x2_2 == -1.) B_x2_start = B_x2_0;
    else {
      nearBoundary_B = 0;
     
      B_x2_0 = (calibTrackB.track.electrons.at(min_b-1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b-1).y_mod-y_midB)*sin(phiB);
      B_x2_2 = (calibTrackB.track.electrons.at(min_b+1).x_mod-x_midB)*cos(phiB)-(calibTrackB.track.electrons.at(min_b+1).y_mod-y_midB)*sin(phiB);
      if (B_x2_0 < B_x2_2) {
	B_x2_start = B_x2_0;
	B_x2_end = B_x2_2;
      } else {
	B_x2_start = B_x2_2;
	B_x2_end = B_x2_0;
      }
    }

    if (nearBoundary_A == 1) if (A_x2_start > 0) stepSize_x2_A *= -1.;
    if (nearBoundary_B == 1) if (B_x2_start > 0) stepSize_x2_B *= -1.;

    Double_t A_x = 0., A_y = 0., A_z = 0.;
    Int_t firstIter_A = 1; //1 means it's the first iteration
    for (Double_t A_x2 = A_x2_start; nearBoundary_A*((A_x > 0 && A_x < Lx && A_y > 0 && A_y < Ly && A_z > 0 && A_z < Lz)+firstIter_A)+!nearBoundary_A*(A_x2 < A_x2_end); A_x2 += stepSize_x2_A) {
      firstIter_A = 0;
      
      Double_t A_y2 = aA*pow(A_x2,2)+bA*A_x2;
      A_y = A_y2*cos(phiA) - A_x2*sin(phiA) + y_midA;
      A_x = (A_y2 - A_y*cos(phiA) + y_midA*cos(phiA) + x_midA*sin(phiA)) / sin(phiA);
      A_z = dA*A_x+eA*A_y+fA;
     
      Double_t B_y = 0., B_x = 0., B_z = 0.;
      Int_t firstIter_B = 1;
      for (Double_t B_x2 = B_x2_start; nearBoundary_B*((B_x > 0 && B_x < Lx && B_y > 0 && B_y < Ly && B_z > 0 && B_z < Lz)+firstIter_B)+!nearBoundary_B*(B_x2 < B_x2_end); B_x2 += stepSize_x2_B) {
	firstIter_B = 0;
   
	Double_t B_y2 = aB*pow(B_x2,2)+bB*B_x2;
	B_y = B_y2*cos(phiB) - B_x2*sin(phiB) + y_midB;
	B_x = (B_y2 - B_y*cos(phiB) + y_midB*cos(phiB) + x_midB*sin(phiB)) / sin(phiB);
	B_z = dB*B_x+eB*B_y+fB;

	Double_t distance = sqrt(pow(B_x-A_x,2)+pow(B_y-A_y,2)+pow(B_z-A_z,2));
	if (min_distance_parabola == -1 || distance < min_distance_parabola) {
	  min_distance_parabola = distance;     
	  min_Ax = A_x; min_Ay = A_y; min_Az = A_z;
	  min_Bx = B_x; min_By = B_y; min_Bz = B_z;
	}
      }
    }
  }

  //midpoint between points of closest approach on both lines
  Double_t x_mid = (min_Ax+min_Bx)/2.;
  Double_t y_mid = (min_Ay+min_By)/2.;
  Double_t z_mid = (min_Az+min_Bz)/2.;

  if(isnan(min_distance_parabola)) {
    return_vector.push_back(999999.0);
    return_vector.push_back(999999.0);
    return_vector.push_back(999999.0);
    return_vector.push_back(999999.0);    
  }
  else {
    return_vector.push_back(min_distance_parabola);
    return_vector.push_back(x_mid);
    return_vector.push_back(y_mid);
    return_vector.push_back(z_mid);
  }

  ////// NEW 12/5/2017 //////
  return_vector.push_back(min_Bx-min_Ax);
  return_vector.push_back(min_By-min_Ay);
  return_vector.push_back(min_Bz-min_Az);
  ///////////////////////////

  return return_vector;
}

void getLArSoftTrackSet(vector<trackInfo> &tracks, Int_t inputType, Int_t calibMode, Int_t maxCosmicTracks, Double_t minTrackMCS_anode, Double_t minTrackMCS_cathode, Double_t minTrackMCS_crossing)
{
  tracks.clear();

  TFile* mapFile;

  if(inputType == 1)
  {
    mapFile = new TFile(inputFileLaser,"READ");

    TTreeReader readerLasers("lasers", mapFile);
    TTreeReader readerTracks("tracks", mapFile);

    TTreeReaderValue<TVector3> track_start(readerLasers, "entry");
    TTreeReaderValue<TVector3> track_end(readerLasers, "exit");
    TTreeReaderArray<TVector3> track_points(readerTracks, "track");

    Int_t nTracks = 0;
    while (readerLasers.Next())
    {
      readerTracks.Next();
    
      nTracks++;
    
      Double_t x0, y0, z0;
      Double_t x1, y1, z1;
    
      if ((*track_start).y() > (*track_end).y()) {
        x0 = doCoordTransformX((*track_start).x());
        y0 = doCoordTransformY((*track_start).y());
        z0 = doCoordTransformZ((*track_start).z());
    
        x1 = doCoordTransformX((*track_end).x());
        y1 = doCoordTransformY((*track_end).y());
        z1 = doCoordTransformZ((*track_end).z());
      }
      else {
        x0 = doCoordTransformX((*track_end).x());
        y0 = doCoordTransformY((*track_end).y());
        z0 = doCoordTransformZ((*track_end).z());
    
        x1 = doCoordTransformX((*track_start).x());
        y1 = doCoordTransformY((*track_start).y());
        z1 = doCoordTransformZ((*track_start).z());
      }
      
      Double_t trackLength = sqrt(pow(x0-x1,2.0)+pow(y0-y1,2.0)+pow(z0-z1,2.0));
      
      Double_t theta = acos((y1-y0)/trackLength);
      Double_t phi = acos((z1-z0)/(trackLength*sin(theta)));
      if (x1 > x0) {
        phi = -1.0*fabs(phi);
      }
      else {
        phi = fabs(phi);
      }
    
      trackInfo track;
      elecInfo electron;
    
      track.energy = 1000000.0;
      track.pdgID = 22;
      track.x0 = x0;
      track.y0 = y0;
      track.z0 = z0;
      track.x1 = x1;
      track.y1 = y1;
      track.z1 = z1;
      track.theta = theta;
      track.phi = phi;
      track.electrons.clear();
    
      for(Int_t j = 0; j < track_points.GetSize(); j++)
      {
        if (((j % 10) != 0) && (j != track_points.GetSize()-1)) continue;
        
        electron.x = -1.0; // Dummy (currently don't save this info in file)
        electron.y = -1.0; // Dummy (currently don't save this info in file)
        electron.z = -1.0; // Dummy (currently don't save this info in file)
        electron.t = -1.0; // Dummy (currently don't save this info in file)
        electron.x_mod = min(max(0.0,doCoordTransformX(track_points[j].x())),Lx);
        electron.y_mod = min(max(0.0,doCoordTransformY(track_points[j].y())),Ly);
      	electron.z_mod = min(max(0.0,doCoordTransformZ(track_points[j].z())),Lz);
        electron.t_mod = -1.0; // Dummy (currently don't save this info in file)
        electron.fate = 0; // Dummy (currently don't save this info in file)

        track.electrons.push_back(electron);
      }
      
      tracks.push_back(track);
    }
  }
  else if(inputType == 2)
  {
    mapFile = new TFile(inputFileCosmic,"READ");
    
    char* filestring = (char*) "";
    char* filestring2 = (char*) "";
    if (isSCEon == false) {
      filestring = (char*) "_copy";
      filestring2 = (char*) "t0ana/";
    }
    
    TTreeReader reader(Form("%sSCEtree%s",filestring2,filestring), mapFile);
    TTreeReaderValue<Double_t> track_startX(reader, Form("track_startX%s",filestring));
    TTreeReaderValue<Double_t> track_startY(reader, Form("track_startY%s",filestring));
    TTreeReaderValue<Double_t> track_startZ(reader, Form("track_startZ%s",filestring));
    TTreeReaderValue<Double_t> track_endX(reader, Form("track_endX%s",filestring));
    TTreeReaderValue<Double_t> track_endY(reader, Form("track_endY%s",filestring));
    TTreeReaderValue<Double_t> track_endZ(reader, Form("track_endZ%s",filestring));
    TTreeReaderValue<Int_t> nPoints(reader, Form("track_nPoints%s",filestring));
    TTreeReaderArray<Double_t> pointX(reader, Form("track_pointX%s",filestring));
    TTreeReaderArray<Double_t> pointY(reader, Form("track_pointY%s",filestring));
    TTreeReaderArray<Double_t> pointZ(reader, Form("track_pointZ%s",filestring));
    TTreeReaderValue<Double_t> track_MCS(reader, Form("track_MCS_measurement%s",filestring));
    TTreeReaderValue<Double_t> track_t0(reader, Form("track_t0%s",filestring));
    
    trackInfo track;
    elecInfo electron;

    TH1F *Zhist_1s = new TH1F("Zhist_1s","",50,0.0,10.0);
    TH1F *Zhist_2s = new TH1F("Zhist_2s","",50,0.0,10.0);
    TH1F *Zhist_3s = new TH1F("Zhist_3s","",50,0.0,10.0);
    TH1F *Zhist_4s = new TH1F("Zhist_4s","",50,0.0,10.0);
    TH1F *Zhist_5s = new TH1F("Zhist_5s","",50,0.0,10.0);
    TH1F *Zhist_6s = new TH1F("Zhist_6s","",50,0.0,10.0);
    TH1F *Zhist_7s = new TH1F("Zhist_7s","",50,0.0,10.0);

    TH1F *Zhist_1e = new TH1F("Zhist_1e","",50,0.0,10.0);
    TH1F *Zhist_2e = new TH1F("Zhist_2e","",50,0.0,10.0);
    TH1F *Zhist_3e = new TH1F("Zhist_3e","",50,0.0,10.0);
    TH1F *Zhist_4e = new TH1F("Zhist_4e","",50,0.0,10.0);
    TH1F *Zhist_5e = new TH1F("Zhist_5e","",50,0.0,10.0);
    TH1F *Zhist_6e = new TH1F("Zhist_6e","",50,0.0,10.0);
    TH1F *Zhist_7e = new TH1F("Zhist_7e","",50,0.0,10.0);

    TF1 *weightFunc;
    if (isMC == true) {
      weightFunc = new TF1("weightFunc","(1.0/75.0)*(x-5.0)^2 + (2.0/3.0)",0.0,10.0);
    }
    else {
      weightFunc = new TF1("weightFunc","(2.0/75.0)*(x-5.0)^2 + (1.0/3.0)",0.0,10.0);
    }
    
    TRandom3 *rand = new TRandom3(0);

    Int_t inputTrackNum = -1;
    Int_t nTracks = 0;
    while (reader.Next())
    {
      inputTrackNum++;
      if ((inputTrackNum < minInputTrackNum) || (inputTrackNum > maxInputTrackNum)) continue;
      
      if ((maxCosmicTracks != -1) && (nTracks >= maxCosmicTracks)) continue;
      
      if (*nPoints < 3) continue;
      
      Double_t xS, yS, zS;
      Double_t xE, yE, zE;
      Double_t x0, y0, z0;
      Double_t x1, y1, z1;
      
      if (*track_startY > *track_endY) {
        xS = doCoordTransformX(*track_startX);
        yS = doCoordTransformY(*track_startY);
        zS = doCoordTransformZ(*track_startZ);
    
        xE = doCoordTransformX(*track_endX);
        yE = doCoordTransformY(*track_endY);
        zE = doCoordTransformZ(*track_endZ);
      }
      else {
        xS = doCoordTransformX(*track_endX);
        yS = doCoordTransformY(*track_endY);
        zS = doCoordTransformZ(*track_endZ);
    
        xE = doCoordTransformX(*track_startX);
        yE = doCoordTransformY(*track_startY);
        zE = doCoordTransformZ(*track_startZ);
      }

      Zhist_1s->Fill(zS);
      Zhist_1e->Fill(zE);
      
      if (((xS < maxXdist) && (xE < maxXdist)) || ((xS > (Lx - maxXdist)) && (xE > (Lx - maxXdist))) || ((xS < maxXdist) && (yS < maxYdist)) || ((xS < maxXdist) && (yS > (Ly - maxYdist))) || ((xS < maxXdist) && (zS < maxZdist)) || ((xS < maxXdist) && (zS > (Lz - maxZdist))) || ((xE < maxXdist) && (yE < maxYdist)) || ((xE < maxXdist) && (yE > (Ly - maxYdist))) || ((xE < maxXdist) && (zE < maxZdist)) || ((xE < maxXdist) && (zE > (Lz - maxZdist))) || ((yS < maxYdist) && (zS < maxZdist)) || ((yS > (Ly - maxYdist)) && (zS < maxZdist)) || ((yS < maxYdist) && (zS > (Lz - maxZdist))) || ((yS > (Ly - maxYdist)) && (zS > (Lz - maxZdist))) || ((yE < maxYdist) && (zE < maxZdist)) || ((yE > (Ly - maxYdist)) && (zE < maxZdist)) || ((yE < maxYdist) && (zE > (Lz - maxZdist))) || ((yE > (Ly - maxYdist)) && (zE > (Lz - maxZdist)))) continue;

      if (((xS > (Lx - maxXdist)) && ((yS < maxYdist)  || (yS > (Ly - maxYdist)) || (zS < maxZdist) || (zS > (Lz - maxZdist)))) || ((xE > (Lx - maxXdist)) && ((yE < maxYdist)  || (yE > (Ly - maxYdist)) || (zE < maxZdist) || (zE > (Lz - maxZdist))))) continue; // NEW
      
      Zhist_2s->Fill(zS);
      Zhist_2e->Fill(zE);
      
      if (((xS > maxXdist) && (xS < (Lx - maxXdist)) && (yS > maxYdist) && (yS < (Ly - maxYdist)) && (zS > maxZdist) && (zS < (Lz - maxZdist))) || ((xE > maxXdist) && (xE < (Lx - maxXdist)) && (yE > maxYdist) && (yE < (Ly - maxYdist)) && (zE > maxZdist) && (zE < (Lz - maxZdist)))) continue;

      Zhist_3s->Fill(zS);
      Zhist_3e->Fill(zE);
    
      if ((((xS > (Lx - maxXdist)) && (xE > maxXdist)) || ((xE > (Lx - maxXdist)) && (xS > maxXdist))) && (*track_MCS < 1000.0*minTrackMCS_anode)) continue;

      Zhist_4s->Fill(zS);
      Zhist_4e->Fill(zE);
      
      if ((((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) && (*track_MCS < 1000.0*minTrackMCS_cathode)) continue;

      Zhist_5s->Fill(zS);
      Zhist_5e->Fill(zE);

      if ((((xS > (Lx - maxXdist)) && (xE < maxXdist)) || ((xE > (Lx - maxXdist)) && (xS < maxXdist))) && (*track_MCS < 1000.0*minTrackMCS_crossing)) continue;

      Zhist_6s->Fill(zS);
      Zhist_6e->Fill(zE);

      double randNum = rand->Uniform(1.0);
      if (randNum > max(weightFunc->Eval(zS),weightFunc->Eval(zE))) continue;

      Zhist_7s->Fill(zS);
      Zhist_7e->Fill(zE);

      nTracks++;

      double SCEfactor = 1.0;
      if (isSCEon == false) {
        SCEfactor = 0.0;
      }
      
      // Correct track end point furthest from cathode/anode
      Double_t cathodeOffset = 0.0;
      if((calibMode == 2) || (calibMode == 4)) {
        if (((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) {
          if (xS < xE) {
            cathodeOffset = findCathodeOffset(yS,zS);
          }
          else {
            cathodeOffset = findCathodeOffset(yE,zE);
          }
        
          xS += SCEfactor*(TrueCathode-ShiftedCathode) + SCEfactor*cathodeOffset;
          xE += SCEfactor*(TrueCathode-ShiftedCathode) + SCEfactor*cathodeOffset;
        }
        else if (((xS > (Lx - maxXdist)) && (xE > maxXdist)) || ((xE > (Lx - maxXdist)) && (xS > maxXdist))) {
          xS += (TrueAnode-ShiftedAnode);
          xE += (TrueAnode-ShiftedAnode);
        }
      }
      
      vector<Double_t> StruthOffsets = getTruthOffsets(xS,yS,zS);
      vector<Double_t> EtruthOffsets = getTruthOffsets(xE,yE,zE);

      if(calibMode == 1) {
        if(xS < maxXdist) {
          x0 = TrueCathode;
          y0 = yS;
          z0 = zS;
        }
        else if (xS > (Lx - maxXdist)) {
          x0 = TrueAnode;
          y0 = yS;
          z0 = zS;
        }
        else if (min(fabs(yS),fabs(Ly-yS)) < min(fabs(zS),fabs(Lz-zS))) {
          if (fabs(yS) < fabs(Ly-yS)) {
            x0 = xS;// + 0.16*(yS-TrueBottom);//xS;
            y0 = TrueBottom;
            z0 = zS;
          }
          else {
            x0 = xS;// + 0.16*(TrueTop-yS);//xS;
            y0 = TrueTop;
            z0 = zS;
          }
        }
        else {
          if (fabs(zS) < fabs(Lz-zS)) {
            x0 = xS;// + 0.16*(zS-TrueUpstream);//xS;
            y0 = yS;
            z0 = TrueUpstream;
          }
          else {
            x0 = xS;// + 0.16*(TrueDownstream-zS);//xS;
            y0 = yS;
            z0 = TrueDownstream;
          }
        }
        
        if(xE < maxXdist) {
          x1 = TrueCathode;
          y1 = yE;
          z1 = zE;
        }
        else if (xE > (Lx - maxXdist)) {
          x1 = TrueAnode;
          y1 = yE;
          z1 = zE;
        }
        else if (min(fabs(yE),fabs(Ly-yE)) < min(fabs(zE),fabs(Lz-zE))) {
          if (fabs(yE) < fabs(Ly-yE)) {
            x1 = xE;// + 0.16*(yE-TrueBottom);//xE;
            y1 = TrueBottom;
            z1 = zE;
          }
          else {
            x1 = xE;// + 0.16*(TrueTop-yE);//xE;
            y1 = TrueTop;
            z1 = zE;
          }
        }
        else {
          if (fabs(zE) < fabs(Lz-zE)) {
            x1 = xE;// + 0.16*(zE-TrueUpstream);//xE;
            y1 = yE;
            z1 = TrueUpstream;
          }
          else {
            x1 = xE;// + 0.16*(TrueDownstream-zE);//xE;
            y1 = yE;
            z1 = TrueDownstream;
          }
        }
      }
      else if(calibMode == 2) {
        if(xS < maxXdist) {
          x0 = TrueCathode;
          y0 = yS + SCEfactor*StruthOffsets.at(1);
          z0 = zS + SCEfactor*StruthOffsets.at(2);
        }
        else if (xS > (Lx - maxXdist)) {
          x0 = TrueAnode;
          y0 = yS;
          z0 = zS;
        }
        else if (min(fabs(yS),fabs(Ly-yS)) < min(fabs(zS),fabs(Lz-zS))) {
          if (fabs(yS) < fabs(Ly-yS)) {
            x0 = xS + SCEfactor*StruthOffsets.at(0);
            y0 = TrueBottom;
            z0 = zS + SCEfactor*StruthOffsets.at(2);
          }
          else {
            x0 = xS + SCEfactor*StruthOffsets.at(0);
            y0 = TrueTop;
            z0 = zS + SCEfactor*StruthOffsets.at(2);
          }
        }
        else {
          if (fabs(zS) < fabs(Lz-zS)) {
            x0 = xS + SCEfactor*StruthOffsets.at(0);
            y0 = yS + SCEfactor*StruthOffsets.at(1);
            z0 = TrueUpstream;
          }
          else {
            x0 = xS + SCEfactor*StruthOffsets.at(0);
            y0 = yS + SCEfactor*StruthOffsets.at(1);
            z0 = TrueDownstream;
          }
        }
        
        if(xE < maxXdist) {
          x1 = TrueCathode;
          y1 = yE + SCEfactor*EtruthOffsets.at(1);
          z1 = zE + SCEfactor*EtruthOffsets.at(2);
        }
        else if (xE > (Lx - maxXdist)) {
          x1 = TrueAnode;
          y1 = yE;
          z1 = zE;
        }
        else if (min(fabs(yE),fabs(Ly-yE)) < min(fabs(zE),fabs(Lz-zE))) {
          if (fabs(yE) < fabs(Ly-yE)) {
            x1 = xE + SCEfactor*EtruthOffsets.at(0);
            y1 = TrueBottom;
            z1 = zE + SCEfactor*EtruthOffsets.at(2);
          }
          else {
            x1 = xE + SCEfactor*EtruthOffsets.at(0);
            y1 = TrueTop;
            z1 = zE + SCEfactor*EtruthOffsets.at(2);
          }
        }
        else {
          if (fabs(zE) < fabs(Lz-zE)) {
            x1 = xE + SCEfactor*EtruthOffsets.at(0);
            y1 = yE + SCEfactor*EtruthOffsets.at(1);
            z1 = TrueUpstream;
          }
          else {
            x1 = xE + SCEfactor*EtruthOffsets.at(0);
            y1 = yE + SCEfactor*EtruthOffsets.at(1);
            z1 = TrueDownstream;
          }
        }
      }
      else if(calibMode == 3) {
        if(xS < maxXdist) {
          x0 = TrueCathode;
          y0 = yS + getFaceCorr(xS,yS,zS,4,2);
          z0 = zS + getFaceCorr(xS,yS,zS,4,3);
        }
        else if (xS > (Lx - maxXdist)) {
          x0 = TrueAnode;
          y0 = yS;
          z0 = zS;
        }
        else if (min(fabs(yS),fabs(Ly-yS)) < min(fabs(zS),fabs(Lz-zS))) {
          if (fabs(yS) < fabs(Ly-yS)) {
            x0 = xS + getFaceCorr(xS,yS,zS,1,1);
            y0 = TrueBottom;
            z0 = zS + getFaceCorr(xS,yS,zS,1,3);
          }
          else {
            x0 = xS + getFaceCorr(xS,yS,zS,0,1);
            y0 = TrueTop;
            z0 = zS + getFaceCorr(xS,yS,zS,0,3);
          }
        }
        else {
          if (fabs(zS) < fabs(Lz-zS)) {
            x0 = xS + getFaceCorr(xS,yS,zS,2,1);
            y0 = yS + getFaceCorr(xS,yS,zS,2,2);
            z0 = TrueUpstream;
          }
          else {
            x0 = xS + getFaceCorr(xS,yS,zS,3,1);
            y0 = yS + getFaceCorr(xS,yS,zS,3,2);
            z0 = TrueDownstream;
          }
        }
        
        if(xE < maxXdist) {
          x1 = TrueCathode;
          y1 = yE + getFaceCorr(xE,yE,zE,4,2);
          z1 = zE + getFaceCorr(xE,yE,zE,4,3);
        }
        else if (xE > (Lx - maxXdist)) {
          x1 = TrueAnode;
          y1 = yE;
          z1 = zE;
        }
        else if (min(fabs(yE),fabs(Ly-yE)) < min(fabs(zE),fabs(Lz-zE))) {
          if (fabs(yE) < fabs(Ly-yE)) {
            x1 = xE + getFaceCorr(xE,yE,zE,1,1);
            y1 = TrueBottom;
            z1 = zE + getFaceCorr(xE,yE,zE,1,3);
          }
          else {
            x1 = xE + getFaceCorr(xE,yE,zE,0,1);
            y1 = TrueTop;
            z1 = zE + getFaceCorr(xE,yE,zE,0,3);
          }
        }
        else {
          if (fabs(zE) < fabs(Lz-zE)) {
            x1 = xE + getFaceCorr(xE,yE,zE,2,1);
            y1 = yE + getFaceCorr(xE,yE,zE,2,2);
            z1 = TrueUpstream;
          }
          else {
            x1 = xE + getFaceCorr(xE,yE,zE,3,1);
            y1 = yE + getFaceCorr(xE,yE,zE,3,2);
            z1 = TrueDownstream;
          }
        }
      }
      else if(calibMode == 4) {
        if(xS < maxXdist) {
          x0 = TrueCathode;
          y0 = yS + SCEfactor*StruthOffsets.at(1);
          z0 = zS + SCEfactor*StruthOffsets.at(2);
        }
        else if (xS > (Lx - maxXdist)) {
          x0 = TrueAnode;
          y0 = yS;
          z0 = zS;
        }
        else if (min(fabs(yS),fabs(Ly-yS)) < min(fabs(zS),fabs(Lz-zS))) {
          if (fabs(yS) < fabs(Ly-yS)) {
            x0 = xS + SCEfactor*(TrueBottom-yS)*StruthOffsets.at(0)/StruthOffsets.at(1);
            y0 = TrueBottom;
            z0 = zS + SCEfactor*(TrueBottom-yS)*StruthOffsets.at(2)/StruthOffsets.at(1);
          }
          else {
            x0 = xS + SCEfactor*(TrueTop-yS)*StruthOffsets.at(0)/StruthOffsets.at(1);
            y0 = TrueTop;
            z0 = zS + SCEfactor*(TrueTop-yS)*StruthOffsets.at(2)/StruthOffsets.at(1);
          }
        }
        else {
          if (fabs(zS) < fabs(Lz-zS)) {
            x0 = xS + SCEfactor*(TrueUpstream-zS)*StruthOffsets.at(0)/StruthOffsets.at(2);
            y0 = yS + SCEfactor*(TrueUpstream-zS)*StruthOffsets.at(1)/StruthOffsets.at(2);
            z0 = TrueUpstream;
          }
          else {
            x0 = xS + SCEfactor*(TrueDownstream-zS)*StruthOffsets.at(0)/StruthOffsets.at(2);
            y0 = yS + SCEfactor*(TrueDownstream-zS)*StruthOffsets.at(1)/StruthOffsets.at(2);
            z0 = TrueDownstream;
          }
        }
        
        if(xE < maxXdist) {
          x1 = TrueCathode;
          y1 = yE + SCEfactor*EtruthOffsets.at(1);
          z1 = zE + SCEfactor*EtruthOffsets.at(2);
        }
        else if (xE > (Lx - maxXdist)) {
          x1 = TrueAnode;
          y1 = yE;
          z1 = zE;
        }
        else if (min(fabs(yE),fabs(Ly-yE)) < min(fabs(zE),fabs(Lz-zE))) {
          if (fabs(yE) < fabs(Ly-yE)) {
            x1 = xE + SCEfactor*(TrueBottom-yE)*EtruthOffsets.at(0)/EtruthOffsets.at(1);
            y1 = TrueBottom;
            z1 = zE + SCEfactor*(TrueBottom-yE)*EtruthOffsets.at(2)/EtruthOffsets.at(1);
          }
          else {
            x1 = xE + SCEfactor*(TrueTop-yE)*EtruthOffsets.at(0)/EtruthOffsets.at(1);
            y1 = TrueTop;
            z1 = zE + SCEfactor*(TrueTop-yE)*EtruthOffsets.at(2)/EtruthOffsets.at(1);
          }
        }
        else {
          if (fabs(zE) < fabs(Lz-zE)) {
            x1 = xE + SCEfactor*(TrueUpstream-zE)*EtruthOffsets.at(0)/EtruthOffsets.at(2);
            y1 = yE + SCEfactor*(TrueUpstream-zE)*EtruthOffsets.at(1)/EtruthOffsets.at(2);
            z1 = TrueUpstream;
          }
          else {
            x1 = xE + SCEfactor*(TrueDownstream-zE)*EtruthOffsets.at(0)/EtruthOffsets.at(2);
            y1 = yE + SCEfactor*(TrueDownstream-zE)*EtruthOffsets.at(1)/EtruthOffsets.at(2);
            z1 = TrueDownstream;
          }
        }
      }
      
      Double_t trackLength = sqrt(pow(x0-x1,2.0)+pow(y0-y1,2.0)+pow(z0-z1,2.0));
      
      //Double_t theta = ArcCos((y1-y0)/trackLength);
      //Double_t phi = ArcSin((x0-x1)/(trackLength*Sin(theta)));
      Double_t theta = acos((y1-y0)/trackLength);
      Double_t phi = acos((z1-z0)/(trackLength*sin(theta)));
      if (x1 > x0) {
        phi = -1.0*fabs(phi);
      }
      else {
        phi = fabs(phi);
      }
      
      track.energy = *track_MCS/1000.0;
      track.pdgID = 13;
      track.x0 = x0;
      track.y0 = y0;
      track.z0 = z0;
      track.x1 = x1;
      track.y1 = y1;
      track.z1 = z1;
      track.theta = theta;
      track.phi = phi;
      track.electrons.clear();

      //if (((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) {
      //  if (x0 < x1) {
      //	  cathodeOffset = findCathodeOffset(y0,z0);
      //	}
      //	else {
      //	  cathodeOffset = findCathodeOffset(y1,z1);
      //  }
      //}

      for(Int_t j = 0; j < *nPoints; j++)
      {
        if (((j % 10) != 0) && (j != *nPoints-1)) continue;
	
        electron.x = -1.0; // Dummy (currently don't save this info in file)
        electron.y = -1.0; // Dummy (currently don't save this info in file)
        electron.z = -1.0; // Dummy (currently don't save this info in file)
        electron.t = -1.0; // Dummy (currently don't save this info in file)
        //electron.x_mod = doCoordTransformX(pointX[j]);
	if(calibMode == 1) {
          if (((xS > (Lx - maxXdist)) && (xE > maxXdist)) || ((xE > (Lx - maxXdist)) && (xS > maxXdist))) {
            electron.x_mod = doCoordTransformX(pointX[j])+(TrueAnode-ShiftedAnode);
          }
          else {
            electron.x_mod = doCoordTransformX(pointX[j]);
          }
	}
	else if(calibMode == 2) {
          if (((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) {
	    electron.x_mod = doCoordTransformX(pointX[j])+SCEfactor*(TrueCathode-ShiftedCathode)+SCEfactor*cathodeOffset;
	  }
          else if (((xS > (Lx - maxXdist)) && (xE > maxXdist)) || ((xE > (Lx - maxXdist)) && (xS > maxXdist))) {
            electron.x_mod = doCoordTransformX(pointX[j])+(TrueAnode-ShiftedAnode);
          }
          else {
            electron.x_mod = doCoordTransformX(pointX[j]);
          }
	}
	else if (calibMode == 3) {
          if (((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) {
	    if(xS < xE) {
	      electron.x_mod = doCoordTransformX(pointX[j])+SCEfactor*(TrueCathode-ShiftedCathode)-getFaceCorr(xS,yS,zS,4,1)-OffsetCathode; // not quite correct (improper inverse operation)
	    }
	    else {
	      electron.x_mod = doCoordTransformX(pointX[j])+SCEfactor*(TrueCathode-ShiftedCathode)-getFaceCorr(xE,yE,zE,4,1)-OffsetCathode; // not quite correct (improper inverse operation)
	    }
	  }
          else if (((xS > (Lx - maxXdist)) && (xE > maxXdist)) || ((xE > (Lx - maxXdist)) && (xS > maxXdist))) {
            electron.x_mod = doCoordTransformX(pointX[j])+(TrueAnode-ShiftedAnode);
          }
          else {
            electron.x_mod = doCoordTransformX(pointX[j]);
          }
	}
	else if(calibMode == 4) {
          if (((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) {
	    electron.x_mod = doCoordTransformX(pointX[j])+SCEfactor*(TrueCathode-ShiftedCathode)+SCEfactor*cathodeOffset;
	  }
          else if (((xS > (Lx - maxXdist)) && (xE > maxXdist)) || ((xE > (Lx - maxXdist)) && (xS > maxXdist))) {
            electron.x_mod = doCoordTransformX(pointX[j])+(TrueAnode-ShiftedAnode);
          }
          else {
            electron.x_mod = doCoordTransformX(pointX[j]);
          }
	}
        electron.y_mod = doCoordTransformY(pointY[j]);
        electron.z_mod = doCoordTransformZ(pointZ[j]);
        electron.t_mod = -1.0; // Dummy (currently don't save this info in file)
        electron.fate = 0; // Dummy (currently don't save this info in file)
    
        track.electrons.push_back(electron);
      }
      
      tracks.push_back(track);
    }

    outputFile->cd();
    
    Zhist_1s->Write();
    Zhist_2s->Write();
    Zhist_3s->Write();
    Zhist_4s->Write();
    Zhist_5s->Write();
    Zhist_6s->Write();
    Zhist_7s->Write();

    Zhist_1e->Write();
    Zhist_2e->Write();
    Zhist_3e->Write();
    Zhist_4e->Write();
    Zhist_5e->Write();
    Zhist_6e->Write();
    Zhist_7e->Write();
  }
}

vector<trackInfo> getTrackSet(Int_t inputType)
{
  vector<trackInfo> tracks;

  TFile* mapFile;
  Char_t* mapName = (Char_t*)"";

  if(inputType == 1)
  {
    mapFile = new TFile(inputFileLaser,"READ");
    mapName = (Char_t*)"laser";
  }
  else if(inputType == 2)
  {
    mapFile = new TFile(inputFileCosmic,"READ");
    mapName = (Char_t*)"cosmic";
  }
  else
  {
    return tracks;
  }

  TTreeReader reader(Form("SpaCEtree_%s",mapName), mapFile);
  TTreeReaderValue<Double_t> x0(reader, Form("x0_%s",mapName));
  TTreeReaderValue<Double_t> y0(reader, Form("y0_%s",mapName));
  TTreeReaderValue<Double_t> z0(reader, Form("z0_%s",mapName));
  TTreeReaderValue<Double_t> theta(reader, Form("theta_%s",mapName));
  TTreeReaderValue<Double_t> phi(reader, Form("phi_%s",mapName));
  TTreeReaderValue<Int_t> nElec(reader, Form("nElec_%s",mapName));
  TTreeReaderArray<Double_t> elecX(reader, Form("elecX_%s",mapName));
  TTreeReaderArray<Double_t> elecY(reader, Form("elecY_%s",mapName));
  TTreeReaderArray<Double_t> elecZ(reader, Form("elecZ_%s",mapName));

  trackInfo track;
  elecInfo electron;
  
  while (reader.Next())
  {
    track.energy = -1.0;
    track.pdgID = -1;
    track.x0 = *x0;
    track.y0 = *y0;
    track.z0 = *z0;
    track.x1 = -1.0;
    track.y1 = -1.0;
    track.z1 = -1.0;
    track.theta = *theta;
    track.phi = *phi;
    track.electrons.clear();
    
    for(Int_t j = 0; j < *nElec; j++)
    {
      electron.x = -1.0; // Dummy (currently don't save this info in file)
      electron.y = -1.0; // Dummy (currently don't save this info in file)
      electron.z = -1.0; // Dummy (currently don't save this info in file)
      electron.t = -1.0; // Dummy (currently don't save this info in file)
      electron.x_mod = elecX[j];
      electron.y_mod = elecY[j];
      electron.z_mod = elecZ[j];
      electron.t_mod = -1.0; // Dummy (currently don't save this info in file)
      electron.fate = 0; // Dummy (currently don't save this info in file)

      track.electrons.push_back(electron);
    }
    
    tracks.push_back(track);
  }

  return tracks;
}

vector<calibTrackInfo> makeCalibTracks(const vector<trackInfo> &tracks)
{
  vector<calibTrackInfo> calibTracks;
  calibTrackInfo calibTrack;
  trackInfo track;

  for(Int_t i = 0; i < tracks.size(); i++)
  {
    track = tracks.at(i);

    calibTrack.track = track;
    calibTrack.x0_calib = track.x0;
    calibTrack.y0_calib = track.y0;
    calibTrack.z0_calib = track.z0;
    calibTrack.x1_calib = track.x1;
    calibTrack.y1_calib = track.y1;
    calibTrack.z1_calib = track.z1;
    calibTrack.theta_calib = track.theta;
    calibTrack.phi_calib = track.phi;

    calibTrack.DxVec.clear();
    calibTrack.DyVec.clear();
    calibTrack.DzVec.clear();
    for(Int_t j = 0; j < track.electrons.size(); j++)
    {
      calibTrack.DxVec.push_back(0.0);
      calibTrack.DyVec.push_back(0.0);
      calibTrack.DzVec.push_back(0.0);
    }

    calibTrack.calibFlag = false;

    calibTracks.push_back(calibTrack);
  }

  return calibTracks;
}

void doLaserLaserCalib(const vector<calibTrackInfo> &laserCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t saveInfo)
{
  calibTrackInfo calibTrackA;
  calibTrackInfo calibTrackB;

  Int_t numLaserTracks = laserCalibTracks.size();

  vector<Double_t> POAparams;
  vector<Double_t> POAparamsDistorted;
  Double_t distVal;
  Double_t distValDistorted;
  Double_t distWeight;
  Double_t xVal;
  Double_t yVal;
  Double_t zVal;
  Double_t xValDistorted;
  Double_t yValDistorted;
  Double_t zValDistorted;
  Int_t xCalibLowIndex;
  Int_t xCalibHighIndex;
  Int_t yCalibLowIndex;
  Int_t yCalibHighIndex;
  Int_t zCalibLowIndex;
  Int_t zCalibHighIndex;
  Double_t xCalibFrac;
  Double_t yCalibFrac;
  Double_t zCalibFrac;
  Double_t tempFactor;

  Double_t crossDistX;
  Double_t crossDistY;
  Double_t crossDistZ;
  Double_t crossDistX_mod;
  Double_t crossDistY_mod;
  Double_t crossDistZ_mod;
  
  Int_t trackNum1;
  Int_t trackNum2;
  Int_t crossType = 1;
  
  TTree *T_crossings = new TTree("SpaCEtree_crossings","SpaCEtree_crossings");
  T_crossings->Branch("trackNum1",&trackNum1,"data_crossings/I");
  T_crossings->Branch("trackNum2",&trackNum2,"data_crossings/I");  
  T_crossings->Branch("crossX",&xVal,"data_crossings/D");
  T_crossings->Branch("crossY",&yVal,"data_crossings/D");
  T_crossings->Branch("crossZ",&zVal,"data_crossings/D");
  T_crossings->Branch("crossDist",&distVal,"data_crossings/D");
  T_crossings->Branch("crossDistX",&crossDistX,"data_crossings/D");
  T_crossings->Branch("crossDistY",&crossDistY,"data_crossings/D");
  T_crossings->Branch("crossDistZ",&crossDistZ,"data_crossings/D");
  T_crossings->Branch("crossX_mod",&xValDistorted,"data_crossings/D");
  T_crossings->Branch("crossY_mod",&yValDistorted,"data_crossings/D");
  T_crossings->Branch("crossZ_mod",&zValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDist_mod",&distValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDistX_mod",&crossDistX_mod,"data_crossings/D");
  T_crossings->Branch("crossDistY_mod",&crossDistY_mod,"data_crossings/D");
  T_crossings->Branch("crossDistZ_mod",&crossDistZ_mod,"data_crossings/D");
  T_crossings->Branch("distWeight",&distWeight,"data_crossings/D");
  T_crossings->Branch("crossType",&crossType,"data_crossings/I");
  if(saveInfo == 1) {
    T_crossings->SetDirectory(outputFile);
  }

  for(Int_t i = 0; i < numLaserTracks; i++)
  {
    cout << "LASER-LASER " << i << endl;

    calibTrackA = laserCalibTracks.at(i);
    if(calibTrackA.track.electrons.size() < 3) continue;

    for(Int_t j = i+1; j < numLaserTracks; j++)
    {
      calibTrackB = laserCalibTracks.at(j);
      if(calibTrackB.track.electrons.size() < 3) continue;

      double dTheta = fabs(ArcCos(((calibTrackA.track.x1-calibTrackA.track.x0)*(calibTrackB.track.x1-calibTrackB.track.x0) + (calibTrackA.track.y1-calibTrackA.track.y0)*(calibTrackB.track.y1-calibTrackB.track.y0) + (calibTrackA.track.z1-calibTrackA.track.z0)*(calibTrackB.track.z1-calibTrackB.track.z0))/(sqrt(pow((calibTrackA.track.x1-calibTrackA.track.x0),2.0)+pow((calibTrackA.track.y1-calibTrackA.track.y0),2.0)+pow((calibTrackA.track.z1-calibTrackA.track.z0),2.0))*sqrt(pow((calibTrackB.track.x1-calibTrackB.track.x0),2.0)+pow((calibTrackB.track.y1-calibTrackB.track.y0),2.0)+pow((calibTrackB.track.z1-calibTrackB.track.z0),2.0)))));
      if (dTheta > piVal/2.0) {
        dTheta = piVal - dTheta;
      }
      if(dTheta < relAngleCut*(piVal/180.0)) continue;

      POAparams = findClosestPOA(calibTrackA,calibTrackB);
      distVal = POAparams.at(0);
      if(isnan(distVal)) continue;
      
      if((distVal < 0.0) || (distVal > maxDistFactor*distScale)) continue;
      
      POAparamsDistorted = findDistortedClosestPOA(calibTrackA,calibTrackB);
      distValDistorted = POAparamsDistorted.at(0);
      if(isnan(distValDistorted)) continue;

      //if(distValDistorted > 3.0*distVal+0.005) continue;
      if(fabs(distValDistorted-distVal)>0.003) continue;

      distWeight = exp(-1.0*(distVal/distScale));
      
      xVal = POAparams.at(1);
      yVal = POAparams.at(2);
      zVal = POAparams.at(3);
      xValDistorted = POAparamsDistorted.at(1);
      yValDistorted = POAparamsDistorted.at(2);
      zValDistorted = POAparamsDistorted.at(3);

      crossDistX = POAparams.at(4);
      crossDistY = POAparams.at(5);
      crossDistZ = POAparams.at(6);      
      crossDistX_mod = POAparamsDistorted.at(4);
      crossDistY_mod = POAparamsDistorted.at(5);
      crossDistZ_mod = POAparamsDistorted.at(6);      
      
      trackNum1 = i;
      trackNum2 = j;
      if(saveInfo == 1) {
        T_crossings->Fill();
      }
      
      xCalibLowIndex = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
      xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
      yCalibLowIndex = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
      yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
      zCalibLowIndex = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
      zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);

      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xCalibLowIndex << " " << yCalibLowIndex << " " << zCalibLowIndex << " " << xCalibHighIndex << " " << yCalibHighIndex << " " << zCalibHighIndex << " " << xCalibFrac << " " << yCalibFrac << " " << zCalibFrac << " " << calibTrackA.track.electrons.size() << " " << calibTrackB.track.electrons.size() << endl;
      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << " Ax1 " << calibTrackA.track.x1 << " Ax0 " << calibTrackA.track.x0 << " Ay1 " << calibTrackA.track.y1 << " Ay0 " << calibTrackA.track.y0 << " Az1 " << calibTrackA.track.z1 << " Az0 " << calibTrackA.track.z0 <<  " Bx1 " << calibTrackB.track.x1 << " Bx0 " << calibTrackB.track.x0 << " By1 " << calibTrackB.track.y1 << " By0 " << calibTrackB.track.y0 << " Bz1 " << calibTrackB.track.z1 << " Bz0 " << calibTrackB.track.z0 << endl;
      cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << endl;

      xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((Double_t) xCalibLowIndex);
      yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((Double_t) yCalibLowIndex);
      zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((Double_t) zCalibLowIndex);

      if(xValDistorted < 0.0) {
        xCalibLowIndex = 0;
	xCalibHighIndex = 1;
	xCalibFrac = 0.0;
      }
      else if (xValDistorted > Lx) {
        xCalibLowIndex = nCalibDivisions_x - 1;
	xCalibHighIndex = nCalibDivisions_x;
	xCalibFrac = 1.0;
      }

      if(yValDistorted < 0.0) {
        yCalibLowIndex = 0;
	yCalibHighIndex = 1;
	yCalibFrac = 0.0;
      }
      else if (yValDistorted > Ly) {
        yCalibLowIndex = nCalibDivisions_y - 1;
	yCalibHighIndex = nCalibDivisions_y;
	yCalibFrac = 1.0;
      }

      if(zValDistorted < 0.0) {
        zCalibLowIndex = 0;
	zCalibHighIndex = 1;
	zCalibFrac = 0.0;
      }
      else if (zValDistorted > Lz) {
        zCalibLowIndex = nCalibDivisions_z - 1;
	zCalibHighIndex = nCalibDivisions_z;
	zCalibFrac = 1.0;
      }

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*yCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*yCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back((zVal-zValDistorted));
    }
  }

  return;
}

void doLaserCosmicCalib(const vector<calibTrackInfo> &laserCalibTracks, const vector<calibTrackInfo> &cosmicCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t saveInfo)
{
  calibTrackInfo calibTrackA;
  calibTrackInfo calibTrackB;

  Int_t numLaserTracks = laserCalibTracks.size();
  Int_t numCosmicTracks = cosmicCalibTracks.size();

  vector<Double_t> POAparams;
  vector<Double_t> POAparamsDistorted;
  Double_t distVal;
  Double_t distValDistorted;
  Double_t distWeight;
  Double_t xVal;
  Double_t yVal;
  Double_t zVal;
  Double_t xValDistorted;
  Double_t yValDistorted;
  Double_t zValDistorted;
  Int_t xCalibLowIndex;
  Int_t xCalibHighIndex;
  Int_t yCalibLowIndex;
  Int_t yCalibHighIndex;
  Int_t zCalibLowIndex;
  Int_t zCalibHighIndex;
  Double_t xCalibFrac;
  Double_t yCalibFrac;
  Double_t zCalibFrac;
  Double_t tempFactor;

  Double_t crossDistX;
  Double_t crossDistY;
  Double_t crossDistZ;
  Double_t crossDistX_mod;
  Double_t crossDistY_mod;
  Double_t crossDistZ_mod;
  
  Int_t trackNum1;
  Int_t trackNum2;
  Int_t crossType = 2;
  
  TTree *T_crossings = new TTree("SpaCEtree_crossings","SpaCEtree_crossings");
  T_crossings->Branch("trackNum1",&trackNum1,"data_crossings/I");
  T_crossings->Branch("trackNum2",&trackNum2,"data_crossings/I");  
  T_crossings->Branch("crossX",&xVal,"data_crossings/D");
  T_crossings->Branch("crossY",&yVal,"data_crossings/D");
  T_crossings->Branch("crossZ",&zVal,"data_crossings/D");
  T_crossings->Branch("crossDist",&distVal,"data_crossings/D");
  T_crossings->Branch("crossDistX",&crossDistX,"data_crossings/D");
  T_crossings->Branch("crossDistY",&crossDistY,"data_crossings/D");
  T_crossings->Branch("crossDistZ",&crossDistZ,"data_crossings/D");
  T_crossings->Branch("crossX_mod",&xValDistorted,"data_crossings/D");
  T_crossings->Branch("crossY_mod",&yValDistorted,"data_crossings/D");
  T_crossings->Branch("crossZ_mod",&zValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDist_mod",&distValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDistX_mod",&crossDistX_mod,"data_crossings/D");
  T_crossings->Branch("crossDistY_mod",&crossDistY_mod,"data_crossings/D");
  T_crossings->Branch("crossDistZ_mod",&crossDistZ_mod,"data_crossings/D");
  T_crossings->Branch("distWeight",&distWeight,"data_crossings/D");
  T_crossings->Branch("crossType",&crossType,"data_crossings/I");
  if(saveInfo == 1) {
    T_crossings->SetDirectory(outputFile);
  }

  for(Int_t i = 0; i < numLaserTracks; i++)
  {
    cout << "LASER-COSMIC " << i << endl;

    calibTrackA = laserCalibTracks.at(i);
    if(calibTrackA.track.electrons.size() < 3) continue;

    for(Int_t j = 0; j < numCosmicTracks; j++)
    {
      calibTrackB = cosmicCalibTracks.at(j);
      if(calibTrackB.track.electrons.size() < 3) continue;

      double dTheta = fabs(ArcCos(((calibTrackA.track.x1-calibTrackA.track.x0)*(calibTrackB.track.x1-calibTrackB.track.x0) + (calibTrackA.track.y1-calibTrackA.track.y0)*(calibTrackB.track.y1-calibTrackB.track.y0) + (calibTrackA.track.z1-calibTrackA.track.z0)*(calibTrackB.track.z1-calibTrackB.track.z0))/(sqrt(pow((calibTrackA.track.x1-calibTrackA.track.x0),2.0)+pow((calibTrackA.track.y1-calibTrackA.track.y0),2.0)+pow((calibTrackA.track.z1-calibTrackA.track.z0),2.0))*sqrt(pow((calibTrackB.track.x1-calibTrackB.track.x0),2.0)+pow((calibTrackB.track.y1-calibTrackB.track.y0),2.0)+pow((calibTrackB.track.z1-calibTrackB.track.z0),2.0)))));
      if (dTheta > piVal/2.0) {
        dTheta = piVal - dTheta;
      }
      if(dTheta < relAngleCut*(piVal/180.0)) continue;
      
      POAparams = findClosestPOA(calibTrackA,calibTrackB);
      distVal = POAparams.at(0);
      if(isnan(distVal)) continue;
      
      if((distVal < 0.0) || (distVal > maxDistFactor*distScale)) continue;

      POAparamsDistorted = findDistortedClosestPOA(calibTrackA,calibTrackB);
      distValDistorted = POAparamsDistorted.at(0);
      if(isnan(distValDistorted)) continue;

      //if(distValDistorted > 3.0*distVal+0.005) continue;
      if(fabs(distValDistorted-distVal)>0.003) continue;
      
      distWeight = exp(-1.0*(distVal/distScale));
      
      xVal = POAparams.at(1);
      yVal = POAparams.at(2);
      zVal = POAparams.at(3);
      xValDistorted = POAparamsDistorted.at(1);
      yValDistorted = POAparamsDistorted.at(2);
      zValDistorted = POAparamsDistorted.at(3);

      crossDistX = POAparams.at(4);
      crossDistY = POAparams.at(5);
      crossDistZ = POAparams.at(6);      
      crossDistX_mod = POAparamsDistorted.at(4);
      crossDistY_mod = POAparamsDistorted.at(5);
      crossDistZ_mod = POAparamsDistorted.at(6);      

      trackNum1 = i;
      trackNum2 = j;
      if(saveInfo == 1) {
        T_crossings->Fill();
      }

      xCalibLowIndex = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
      xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
      yCalibLowIndex = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
      yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
      zCalibLowIndex = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
      zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);

      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xCalibLowIndex << " " << yCalibLowIndex << " " << zCalibLowIndex << " " << xCalibHighIndex << " " << yCalibHighIndex << " " << zCalibHighIndex << " " << xCalibFrac << " " << yCalibFrac << " " << zCalibFrac << " " << calibTrackA.track.electrons.size() << " " << calibTrackB.track.electrons.size() << endl;
      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << " Ax1 " << calibTrackA.track.x1 << " Ax0 " << calibTrackA.track.x0 << " Ay1 " << calibTrackA.track.y1 << " Ay0 " << calibTrackA.track.y0 << " Az1 " << calibTrackA.track.z1 << " Az0 " << calibTrackA.track.z0 <<  " Bx1 " << calibTrackB.track.x1 << " Bx0 " << calibTrackB.track.x0 << " By1 " << calibTrackB.track.y1 << " By0 " << calibTrackB.track.y0 << " Bz1 " << calibTrackB.track.z1 << " Bz0 " << calibTrackB.track.z0 << endl;
      cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << endl;

      xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((Double_t) xCalibLowIndex);
      yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((Double_t) yCalibLowIndex);
      zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((Double_t) zCalibLowIndex);

      if(xValDistorted < 0.0) {
        xCalibLowIndex = 0;
	xCalibHighIndex = 1;
	xCalibFrac = 0.0;
      }
      else if (xValDistorted > Lx) {
        xCalibLowIndex = nCalibDivisions_x - 1;
	xCalibHighIndex = nCalibDivisions_x;
	xCalibFrac = 1.0;
      }

      if(yValDistorted < 0.0) {
        yCalibLowIndex = 0;
	yCalibHighIndex = 1;
	yCalibFrac = 0.0;
      }
      else if (yValDistorted > Ly) {
        yCalibLowIndex = nCalibDivisions_y - 1;
	yCalibHighIndex = nCalibDivisions_y;
	yCalibFrac = 1.0;
      }

      if(zValDistorted < 0.0) {
        zCalibLowIndex = 0;
	zCalibHighIndex = 1;
	zCalibFrac = 0.0;
      }
      else if (zValDistorted > Lz) {
        zCalibLowIndex = nCalibDivisions_z - 1;
	zCalibHighIndex = nCalibDivisions_z;
	zCalibFrac = 1.0;
      }

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*yCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*yCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back((zVal-zValDistorted));
    }
  }

  return;
}

void doCosmicCosmicCalib(const vector<calibTrackInfo> &cosmicCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t saveInfo)
{
  calibTrackInfo calibTrackA;
  calibTrackInfo calibTrackB;

  Int_t numCosmicTracks = cosmicCalibTracks.size();

  vector<Double_t> POAparams;
  vector<Double_t> POAparamsDistorted;
  Double_t distVal;
  Double_t distValDistorted;
  Double_t distWeight;
  Double_t xVal;
  Double_t yVal;
  Double_t zVal;
  Double_t xValDistorted;
  Double_t yValDistorted;
  Double_t zValDistorted;
  Int_t xCalibLowIndex;
  Int_t xCalibHighIndex;
  Int_t yCalibLowIndex;
  Int_t yCalibHighIndex;
  Int_t zCalibLowIndex;
  Int_t zCalibHighIndex;
  Double_t xCalibFrac;
  Double_t yCalibFrac;
  Double_t zCalibFrac;
  Double_t tempFactor;

  Double_t crossDistX;
  Double_t crossDistY;
  Double_t crossDistZ;
  Double_t crossDistX_mod;
  Double_t crossDistY_mod;
  Double_t crossDistZ_mod;
  
  Int_t trackNum1;
  Int_t trackNum2;
  Int_t crossType = 3;
  
  TTree *T_crossings = new TTree("SpaCEtree_crossings","SpaCEtree_crossings");
  T_crossings->Branch("trackNum1",&trackNum1,"data_crossings/I");
  T_crossings->Branch("trackNum2",&trackNum2,"data_crossings/I");  
  T_crossings->Branch("crossX",&xVal,"data_crossings/D");
  T_crossings->Branch("crossY",&yVal,"data_crossings/D");
  T_crossings->Branch("crossZ",&zVal,"data_crossings/D");
  T_crossings->Branch("crossDist",&distVal,"data_crossings/D");
  T_crossings->Branch("crossDistX",&crossDistX,"data_crossings/D");
  T_crossings->Branch("crossDistY",&crossDistY,"data_crossings/D");
  T_crossings->Branch("crossDistZ",&crossDistZ,"data_crossings/D");
  T_crossings->Branch("crossX_mod",&xValDistorted,"data_crossings/D");
  T_crossings->Branch("crossY_mod",&yValDistorted,"data_crossings/D");
  T_crossings->Branch("crossZ_mod",&zValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDist_mod",&distValDistorted,"data_crossings/D");
  T_crossings->Branch("crossDistX_mod",&crossDistX_mod,"data_crossings/D");
  T_crossings->Branch("crossDistY_mod",&crossDistY_mod,"data_crossings/D");
  T_crossings->Branch("crossDistZ_mod",&crossDistZ_mod,"data_crossings/D");
  T_crossings->Branch("distWeight",&distWeight,"data_crossings/D");
  T_crossings->Branch("crossType",&crossType,"data_crossings/I");
  if(saveInfo == 1) {
    T_crossings->SetDirectory(outputFile);
  }

  for(Int_t i = 0; i < numCosmicTracks; i++)
  {
    cout << "COSMIC-COSMIC " << i << endl;

    calibTrackA = cosmicCalibTracks.at(i);
    if(calibTrackA.track.electrons.size() < 3) continue;

    for(Int_t j = i+1; j < numCosmicTracks; j++)
    {
      calibTrackB = cosmicCalibTracks.at(j);
      if(calibTrackB.track.electrons.size() < 3) continue;

      double dTheta = fabs(ArcCos(((calibTrackA.track.x1-calibTrackA.track.x0)*(calibTrackB.track.x1-calibTrackB.track.x0) + (calibTrackA.track.y1-calibTrackA.track.y0)*(calibTrackB.track.y1-calibTrackB.track.y0) + (calibTrackA.track.z1-calibTrackA.track.z0)*(calibTrackB.track.z1-calibTrackB.track.z0))/(sqrt(pow((calibTrackA.track.x1-calibTrackA.track.x0),2.0)+pow((calibTrackA.track.y1-calibTrackA.track.y0),2.0)+pow((calibTrackA.track.z1-calibTrackA.track.z0),2.0))*sqrt(pow((calibTrackB.track.x1-calibTrackB.track.x0),2.0)+pow((calibTrackB.track.y1-calibTrackB.track.y0),2.0)+pow((calibTrackB.track.z1-calibTrackB.track.z0),2.0)))));
      if (dTheta > piVal/2.0) {
        dTheta = piVal - dTheta;
      }
      if(dTheta < relAngleCut*(piVal/180.0)) continue;
      
      POAparams = findClosestPOA(calibTrackA,calibTrackB);
      distVal = POAparams.at(0);
      if(isnan(distVal)) continue;

      if((distVal < 0.0) || (distVal > maxDistFactor*distScale)) continue;
      
      POAparamsDistorted = findDistortedClosestPOA(calibTrackA,calibTrackB);
      distValDistorted = POAparamsDistorted.at(0);
      if(isnan(distValDistorted)) continue;
      
      //if(distValDistorted > 3.0*distVal+0.005) continue;
      if(fabs(distValDistorted-distVal)>0.003) continue;
      
      distWeight = exp(-1.0*(distVal/distScale));
      
      xVal = POAparams.at(1);
      yVal = POAparams.at(2);
      zVal = POAparams.at(3);
      xValDistorted = POAparamsDistorted.at(1);
      yValDistorted = POAparamsDistorted.at(2);
      zValDistorted = POAparamsDistorted.at(3);

      crossDistX = POAparams.at(4);
      crossDistY = POAparams.at(5);
      crossDistZ = POAparams.at(6);      
      crossDistX_mod = POAparamsDistorted.at(4);
      crossDistY_mod = POAparamsDistorted.at(5);
      crossDistZ_mod = POAparamsDistorted.at(6);      

      trackNum1 = i;
      trackNum2 = j;
      if(saveInfo == 1) {
        T_crossings->Fill();
      }

      xCalibLowIndex = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
      xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
      yCalibLowIndex = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
      yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
      zCalibLowIndex = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
      zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);
      
      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xCalibLowIndex << " " << yCalibLowIndex << " " << zCalibLowIndex << " " << xCalibHighIndex << " " << yCalibHighIndex << " " << zCalibHighIndex << " " << xCalibFrac << " " << yCalibFrac << " " << zCalibFrac << " " << calibTrackA.track.electrons.size() << " " << calibTrackB.track.electrons.size() << endl;
      //cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << " Ax1 " << calibTrackA.track.x1 << " Ax0 " << calibTrackA.track.x0 << " Ay1 " << calibTrackA.track.y1 << " Ay0 " << calibTrackA.track.y0 << " Az1 " << calibTrackA.track.z1 << " Az0 " << calibTrackA.track.z0 <<  " Bx1 " << calibTrackB.track.x1 << " Bx0 " << calibTrackB.track.x0 << " By1 " << calibTrackB.track.y1 << " By0 " << calibTrackB.track.y0 << " Bz1 " << calibTrackB.track.z1 << " Bz0 " << calibTrackB.track.z0 << endl;
      cout << "  " << j << " " << distVal << " " << distValDistorted << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << dTheta*(180.0/piVal) << endl;

      xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((Double_t) xCalibLowIndex);
      yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((Double_t) yCalibLowIndex);
      zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((Double_t) zCalibLowIndex);

      if(xValDistorted < 0.0) {
        xCalibLowIndex = 0;
	xCalibHighIndex = 1;
	xCalibFrac = 0.0;
      }
      else if (xValDistorted > Lx) {
        xCalibLowIndex = nCalibDivisions_x - 1;
	xCalibHighIndex = nCalibDivisions_x;
	xCalibFrac = 1.0;
      }

      if(yValDistorted < 0.0) {
        yCalibLowIndex = 0;
	yCalibHighIndex = 1;
	yCalibFrac = 0.0;
      }
      else if (yValDistorted > Ly) {
        yCalibLowIndex = nCalibDivisions_y - 1;
	yCalibHighIndex = nCalibDivisions_y;
	yCalibFrac = 1.0;
      }

      if(zValDistorted < 0.0) {
        zCalibLowIndex = 0;
	zCalibHighIndex = 1;
	zCalibFrac = 0.0;
      }
      else if (zValDistorted > Lz) {
        zCalibLowIndex = nCalibDivisions_z - 1;
	zCalibHighIndex = nCalibDivisions_z;
	zCalibFrac = 1.0;
      }

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*(1.0-xCalibFrac)*(1.0-yCalibFrac)*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibLowIndex][zCalibHighIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*(1.0-zCalibFrac);
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*(1.0-xCalibFrac)*yCalibFrac*zCalibFrac;
      calibWeight[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibLowIndex][yCalibHighIndex][zCalibHighIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*(1.0-yCalibFrac)*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibLowIndex][zCalibHighIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*yCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibLowIndex].push_back((zVal-zValDistorted));

      tempFactor = distWeight*xCalibFrac*yCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      calibWeight[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back(tempFactor);
      calibDeltaX[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back((xVal-xValDistorted));
      calibDeltaY[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back((yVal-yValDistorted));
      calibDeltaZ[xCalibHighIndex][yCalibHighIndex][zCalibHighIndex].push_back((zVal-zValDistorted));
    }
  }

  return;
}

void doCalibration(const vector<trackInfo> &laserTracks, const vector<trackInfo> &cosmicTracks, Double_t distScale, Double_t maxDistFactor, Int_t numIterations, Int_t saveInfo)
{
  bulkCalibHistDeltaX->Reset();
  bulkCalibHistDeltaY->Reset();
  bulkCalibHistDeltaZ->Reset();

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
	calibDeltaX[x][y][z].clear();
	calibDeltaY[x][y][z].clear();
	calibDeltaZ[x][y][z].clear();
	calibWeight[x][y][z].clear();
      }
    }
  }
  
  vector<calibTrackInfo> laserCalibTracks = makeCalibTracks(laserTracks);
  vector<calibTrackInfo> cosmicCalibTracks = makeCalibTracks(cosmicTracks);

  Double_t x_true, y_true, z_true;
  Double_t x_reco, y_reco, z_reco;
  Double_t Dx, Dy, Dz;
  Int_t elecFate;

  TTree *T_calib = new TTree("SpaCEtree_calib","SpaCEtree_calib");
  T_calib->Branch("x_true",&x_true,"data_calib/D");
  T_calib->Branch("y_true",&y_true,"data_calib/D");
  T_calib->Branch("z_true",&z_true,"data_calib/D");
  T_calib->Branch("x_reco",&x_reco,"data_calib/D");
  T_calib->Branch("y_reco",&y_reco,"data_calib/D");
  T_calib->Branch("z_reco",&z_reco,"data_calib/D");
  T_calib->Branch("Dx",&Dx,"data_calib/D");
  T_calib->Branch("Dy",&Dy,"data_calib/D");
  T_calib->Branch("Dz",&Dz,"data_calib/D");
  T_calib->Branch("elecFate",&elecFate,"data_calib/I");
  if(saveInfo == 1) {
    T_calib->SetDirectory(outputFile);
  }

  //doLaserLaserCalib(laserCalibTracks,distScale,maxDistFactor,saveInfo);
  //doLaserCosmicCalib(laserCalibTracks,cosmicCalibTracks,distScale,maxDistFactor,saveInfo);
  doCosmicCosmicCalib(cosmicCalibTracks,distScale,maxDistFactor,saveInfo);

  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    y_reco = -1.0*Ly/nCalibDivisions_y;
  
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      y_reco += Ly/nCalibDivisions_y;
      z_reco = -1.0*Lz/nCalibDivisions_z;
  
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        z_reco += Lz/nCalibDivisions_z;

	Int_t length = calibWeight[x][y][z].size();
	Double_t sumWeights = 0.0;
	
	for(Int_t k = 0; k < length; k++)
	{
	  sumWeights += (calibWeight[x][y][z])[k];
	  
          if((calibWeight[x][y][z])[k] < 0.0)
	  {
	    (calibWeight[x][y][z])[k] = 0.0;
	  }
	}
	
        if((length > 0) && (sumWeights > 0.0))
	{
          elecFate = 1;

          Dx = TMath::Median(length,&(calibDeltaX[x][y][z])[0],&(calibWeight[x][y][z])[0]);
          Dy = TMath::Median(length,&(calibDeltaY[x][y][z])[0],&(calibWeight[x][y][z])[0]);
          Dz = TMath::Median(length,&(calibDeltaZ[x][y][z])[0],&(calibWeight[x][y][z])[0]);
         
          x_true = x_reco+Dx;
          y_true = y_reco+Dy;
          z_true = z_reco+Dz;
	}
	else if(x == nCalibDivisions_x)
	{
          elecFate = 1;

          Dx = 0.0;
          Dy = 0.0;
          Dz = 0.0;

          x_true = x_reco;
          y_true = y_reco;
          z_true = z_reco;
	}
	else
	{
          elecFate = 0;

          Dx = 0.0;
          Dy = 0.0;
          Dz = 0.0;

          x_true = x_reco;
          y_true = y_reco;
          z_true = z_reco;
	}
    
        bulkCalibHistDeltaX->SetBinContent(x+1,y+1,z+1,Dx);
        bulkCalibHistDeltaY->SetBinContent(x+1,y+1,z+1,Dy);
        bulkCalibHistDeltaZ->SetBinContent(x+1,y+1,z+1,Dz);

        if(saveInfo == 1) {
          T_calib->Fill();
	}
      }
    }
  }

  return;
}

void saveTrackInfo(const vector<trackInfo> &tracks)
{
  Int_t trackID;

  Double_t x0;
  Double_t y0;
  Double_t z0;
  Double_t x1;
  Double_t y1;
  Double_t z1;
  Double_t theta;
  Double_t phi;
  Double_t energy;

  Int_t nElec;
  Double_t elecX[10000];
  Double_t elecY[10000];
  Double_t elecZ[10000];

  TTree *T_tracks;
  T_tracks = new TTree("SpaCEtree_tracks","SpaCEtree_tracks");
  T_tracks->Branch("trackID_tracks",&trackID,"trackID_tracks/I");
  T_tracks->Branch("x0_tracks",&x0,"x0_tracks/D");  
  T_tracks->Branch("y0_tracks",&y0,"y0_tracks/D");
  T_tracks->Branch("z0_tracks",&z0,"z0_tracks/D");
  T_tracks->Branch("x1_tracks",&x1,"x1_tracks/D");
  T_tracks->Branch("y1_tracks",&y1,"y1_tracks/D");
  T_tracks->Branch("z1_tracks",&z1,"z1_tracks/D");
  T_tracks->Branch("theta_tracks",&theta,"theta_tracks/D");
  T_tracks->Branch("phi_tracks",&phi,"phi_tracks/D");
  T_tracks->Branch("energy_tracks",&energy,"energy_tracks/D");
  T_tracks->Branch("nElec_tracks",&nElec,"nElec_tracks/I");
  T_tracks->Branch("elecX_tracks",&elecX,"elecX_tracks[nElec_tracks]/D");
  T_tracks->Branch("elecY_tracks",&elecY,"elecY_tracks[nElec_tracks]/D");
  T_tracks->Branch("elecZ_tracks",&elecZ,"elecZ_tracks[nElec_tracks]/D");
  T_tracks->SetDirectory(outputFile);

  Int_t numTracks = tracks.size();
  for(Int_t j = 0; j < numTracks; j++)
  {
    trackID = j;
    
    x0 = tracks.at(j).x0;
    y0 = tracks.at(j).y0;
    z0 = tracks.at(j).z0;
    x1 = tracks.at(j).x1;
    y1 = tracks.at(j).y1;
    z1 = tracks.at(j).z1;
    theta = tracks.at(j).theta;
    phi = tracks.at(j).phi;
    energy = tracks.at(j).energy;
    
    //if() continue; // make check to see if track is of "high momentum" (which is?)

    nElec = 0;
    Int_t numElectrons = tracks.at(j).electrons.size();
    while(nElec < tracks.at(j).electrons.size())
    {        
      elecX[nElec] = tracks.at(j).electrons.at(nElec).x_mod;
      elecY[nElec] = tracks.at(j).electrons.at(nElec).y_mod;
      elecZ[nElec] = tracks.at(j).electrons.at(nElec).z_mod;
      nElec++;
    }

    T_tracks->Fill();
  }

  return;
}

void loadTruthMap()
{
  TFile* fileTruth = new TFile("data/dispOutput_MicroBooNE_E273.root");

  TTreeReader reader("SpaCEtree_bkwdDisp", fileTruth);
  TTreeReaderValue<Double_t> reco_x(reader, "x_reco.data_bkwdDisp");
  TTreeReaderValue<Double_t> reco_y(reader, "y_reco.data_bkwdDisp");
  TTreeReaderValue<Double_t> reco_z(reader, "z_reco.data_bkwdDisp");
  TTreeReaderValue<Double_t> Dx(reader, "Dx.data_bkwdDisp");
  TTreeReaderValue<Double_t> Dy(reader, "Dy.data_bkwdDisp");
  TTreeReaderValue<Double_t> Dz(reader, "Dz.data_bkwdDisp");
  TTreeReaderValue<Int_t> elecFate(reader, "elecFate.data_bkwdDisp");

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        trueDeltaX[x][y][z] = 0.0;
        trueDeltaY[x][y][z] = 0.0;
        trueDeltaZ[x][y][z] = 0.0;
      }
    }
  }

  while (reader.Next())
  {
    if (*elecFate == 1) {
      trueDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = *Dx;
      trueDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = *Dy;
      trueDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = *Dz;
    }
    else {
      trueDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = -999;
      trueDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = -999;
      trueDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(*reco_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*reco_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*reco_z/Lz))] = -999;
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueDeltaX[x][y][z] == -999) {
          if (y == nCalibDivisions_y-1) {
            trueDeltaX[x][y][z] = trueDeltaX[x][nCalibDivisions_y-2][z];
            trueDeltaY[x][y][z] = trueDeltaY[x][nCalibDivisions_y-2][z];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][nCalibDivisions_y-2][z];
          }
          else if (y == 1) {
            trueDeltaX[x][y][z] = trueDeltaX[x][2][z];
            trueDeltaY[x][y][z] = trueDeltaY[x][2][z];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][2][z];
          }
          else if (z == nCalibDivisions_z-1) {
            trueDeltaX[x][y][z] = trueDeltaX[x][y][nCalibDivisions_z-2];
            trueDeltaY[x][y][z] = trueDeltaY[x][y][nCalibDivisions_z-2];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][y][nCalibDivisions_z-2];
          }
          else if (z == 1) {
            trueDeltaX[x][y][z] = trueDeltaX[x][y][2];
            trueDeltaY[x][y][z] = trueDeltaY[x][y][2];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][y][2];
          }
	}
      }
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueDeltaX[x][y][z] == -999) {
          if ((y == nCalibDivisions_y) && (z == 0)) {
            trueDeltaX[x][y][z] = trueDeltaX[x][nCalibDivisions_y-1][1];
            trueDeltaY[x][y][z] = trueDeltaY[x][nCalibDivisions_y-1][1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][nCalibDivisions_y-1][1];
          }
          else if ((y == nCalibDivisions_y) && (z == nCalibDivisions_z)) {
            trueDeltaX[x][y][z] = trueDeltaX[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
            trueDeltaY[x][y][z] = trueDeltaY[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
          }
          else if ((y == 0) && (z == 0)) {
            trueDeltaX[x][y][z] = trueDeltaX[x][1][1];
            trueDeltaY[x][y][z] = trueDeltaY[x][1][1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][1][1];
          }
          else if ((y == 0) && (z == nCalibDivisions_z)) {
            trueDeltaX[x][y][z] = trueDeltaX[x][1][nCalibDivisions_z-1];
            trueDeltaY[x][y][z] = trueDeltaY[x][1][nCalibDivisions_z-1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][1][nCalibDivisions_z-1];
          }
	}
      }
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueDeltaX[x][y][z] == -999) {
          if (y == nCalibDivisions_y) {
            trueDeltaX[x][y][z] = trueDeltaX[x][nCalibDivisions_y-1][z];
            trueDeltaY[x][y][z] = trueDeltaY[x][nCalibDivisions_y-1][z];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][nCalibDivisions_y-1][z];
          }
          else if (y == 0) {
            trueDeltaX[x][y][z] = trueDeltaX[x][1][z];
            trueDeltaY[x][y][z] = trueDeltaY[x][1][z];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][1][z];
          }
          else if (z == nCalibDivisions_z) {
            trueDeltaX[x][y][z] = trueDeltaX[x][y][nCalibDivisions_z-1];
            trueDeltaY[x][y][z] = trueDeltaY[x][y][nCalibDivisions_z-1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][y][nCalibDivisions_z-1];
          }
          else if (z == 0) {
            trueDeltaX[x][y][z] = trueDeltaX[x][y][1];
            trueDeltaY[x][y][z] = trueDeltaY[x][y][1];
            trueDeltaZ[x][y][z] = trueDeltaZ[x][y][1];
          }
	}
      }
    }
  }

  return;
}

Double_t getTruthOffset(Double_t xVal, Double_t yVal, Double_t zVal, int comp)
{
  Double_t offset = 0.0;
  
  if (xVal < 0.0) {
    xVal = 0.0;
  }
  if (xVal > Lx) {
    xVal = Lx;
  }

  if (yVal < 0.0) {
    yVal = 0.0;
  }
  if (yVal > Ly) {
    yVal = Ly;
  }

  if (zVal < 0.0) {
    zVal = 0.0;
  }
  if (zVal > Lz) {
    zVal = Lz;
  }

  if (comp == 1) {
    offset = ScaleFactorX*trueDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  else if (comp == 2) {
    offset = ScaleFactorY*trueDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  else if (comp == 3) {
    offset = ScaleFactorZ*trueDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  
  if ((comp == 1) && (offset < 0.0)) {
    offset = 0.0;
  }

  return offset;
}

vector<Double_t> getTruthOffsets(Double_t xVal, Double_t yVal, Double_t zVal)
{
  vector<Double_t> offsets;
  offsets.push_back(getTruthOffset(xVal,yVal,zVal,1));
  offsets.push_back(getTruthOffset(xVal,yVal,zVal,2));
  offsets.push_back(getTruthOffset(xVal,yVal,zVal,3));

  return offsets;
}

void loadTruthFwdMap()
{
  TFile* fileTruth = new TFile("data/dispOutput_MicroBooNE_E273.root");

  TTreeReader reader("SpaCEtree_fwdDisp", fileTruth);
  TTreeReaderValue<Double_t> true_x(reader, "x_true.data_fwdDisp");
  TTreeReaderValue<Double_t> true_y(reader, "y_true.data_fwdDisp");
  TTreeReaderValue<Double_t> true_z(reader, "z_true.data_fwdDisp");
  TTreeReaderValue<Double_t> Dx(reader, "Dx.data_fwdDisp");
  TTreeReaderValue<Double_t> Dy(reader, "Dy.data_fwdDisp");
  TTreeReaderValue<Double_t> Dz(reader, "Dz.data_fwdDisp");
  TTreeReaderValue<Int_t> elecFate(reader, "elecFate.data_fwdDisp");

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        trueFwdDeltaX[x][y][z] = 0.0;
        trueFwdDeltaY[x][y][z] = 0.0;
        trueFwdDeltaZ[x][y][z] = 0.0;
      }
    }
  }

  while (reader.Next())
  {
    if (*elecFate == 1) {
      trueFwdDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = *Dx;
      trueFwdDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = *Dy;
      trueFwdDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = *Dz;
    }
    else {
      trueFwdDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = -999;
      trueFwdDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = -999;
      trueFwdDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(*true_x/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(*true_y/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(*true_z/Lz))] = -999;
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueFwdDeltaX[x][y][z] == -999) {
          if (y == nCalibDivisions_y-1) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][nCalibDivisions_y-2][z];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][nCalibDivisions_y-2][z];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][nCalibDivisions_y-2][z];
          }
          else if (y == 1) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][2][z];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][2][z];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][2][z];
          }
          else if (z == nCalibDivisions_z-1) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][y][nCalibDivisions_z-2];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][y][nCalibDivisions_z-2];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][y][nCalibDivisions_z-2];
          }
          else if (z == 1) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][y][2];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][y][2];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][y][2];
          }
	}
      }
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueFwdDeltaX[x][y][z] == -999) {
          if ((y == nCalibDivisions_y) && (z == 0)) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][nCalibDivisions_y-1][1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][nCalibDivisions_y-1][1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][nCalibDivisions_y-1][1];
          }
          else if ((y == nCalibDivisions_y) && (z == nCalibDivisions_z)) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][nCalibDivisions_y-1][nCalibDivisions_z-1];
          }
          else if ((y == 0) && (z == 0)) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][1][1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][1][1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][1][1];
          }
          else if ((y == 0) && (z == nCalibDivisions_z)) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][1][nCalibDivisions_z-1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][1][nCalibDivisions_z-1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][1][nCalibDivisions_z-1];
          }
	}
      }
    }
  }

  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        if (trueFwdDeltaX[x][y][z] == -999) {
          if (y == nCalibDivisions_y) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][nCalibDivisions_y-1][z];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][nCalibDivisions_y-1][z];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][nCalibDivisions_y-1][z];
          }
          else if (y == 0) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][1][z];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][1][z];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][1][z];
          }
          else if (z == nCalibDivisions_z) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][y][nCalibDivisions_z-1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][y][nCalibDivisions_z-1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][y][nCalibDivisions_z-1];
          }
          else if (z == 0) {
            trueFwdDeltaX[x][y][z] = trueFwdDeltaX[x][y][1];
            trueFwdDeltaY[x][y][z] = trueFwdDeltaY[x][y][1];
            trueFwdDeltaZ[x][y][z] = trueFwdDeltaZ[x][y][1];
          }
	}
      }
    }
  }

  return;
}

Double_t getTruthFwdOffset(Double_t xVal, Double_t yVal, Double_t zVal, int comp)
{
  Double_t offset = 0.0;
  
  if (xVal < 0.0) {
    xVal = 0.0;
  }
  if (xVal > Lx) {
    xVal = Lx;
  }

  if (yVal < 0.0) {
    yVal = 0.0;
  }
  if (yVal > Ly) {
    yVal = Ly;
  }

  if (zVal < 0.0) {
    zVal = 0.0;
  }
  if (zVal > Lz) {
    zVal = Lz;
  }

  if (comp == 1) {
    offset = ScaleFactorX*trueFwdDeltaX[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  else if (comp == 2) {
    offset = ScaleFactorY*trueFwdDeltaY[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  else if (comp == 3) {
    offset = ScaleFactorZ*trueFwdDeltaZ[(Int_t)TMath::Nint(nCalibDivisions_x*(xVal/Lx))][(Int_t)TMath::Nint(nCalibDivisions_y*(yVal/Ly))][(Int_t)TMath::Nint(nCalibDivisions_z*(zVal/Lz))];
  }
  
  if ((comp == 1) && (offset > 0.0)) {
    offset = 0.0;
  }

  return offset;
}

vector<Double_t> getTruthFwdOffsets(Double_t xVal, Double_t yVal, Double_t zVal)
{
  vector<Double_t> fwdOffsets;
  fwdOffsets.push_back(getTruthFwdOffset(xVal,yVal,zVal,1));
  fwdOffsets.push_back(getTruthFwdOffset(xVal,yVal,zVal,2));
  fwdOffsets.push_back(getTruthFwdOffset(xVal,yVal,zVal,3));

  return fwdOffsets;
}

Double_t findCathodeOffset(Double_t yVal, Double_t zVal)
{
  Double_t minMetric = 100000000.0;
  
  Double_t yValNewBest;
  Double_t zValNewBest;
  for(int j = 0; j <= 50; j++)
    for(int k = 0; k <= 50; k++)
    {
      Double_t yValNew = (j/50.0)*Ly;
      Double_t zValNew = (k/50.0)*Lz;
    	
      vector<Double_t> fwdOffsets = getTruthFwdOffsets(TrueCathode,yValNew,zValNew);
    	
      Double_t metric = sqrt(pow(yValNew+fwdOffsets.at(1)-yVal,2.0)+pow(zValNew+fwdOffsets.at(2)-zVal,2.0));

      if(metric < minMetric)
      {
        minMetric = metric;
        yValNewBest = yValNew;
        zValNewBest = zValNew;
      }
    }

  Double_t yValNewBest2;
  Double_t zValNewBest2;
  for(int j = -10; j <= 10; j++)
    for(int k = -10; k <= 10; k++)
    {
      Double_t yValNew = yValNewBest + (j/1000.0)*Ly;
      Double_t zValNew = zValNewBest + (k/1000.0)*Lz;
    	
      vector<Double_t> fwdOffsets = getTruthFwdOffsets(TrueCathode,yValNew,zValNew);
    	
      Double_t metric = sqrt(pow(yValNew+fwdOffsets.at(1)-yVal,2.0)+pow(zValNew+fwdOffsets.at(2)-zVal,2.0));

      if(metric < minMetric)
      {
        minMetric = metric;
        yValNewBest2 = yValNew;
        zValNewBest2 = zValNew;
      }
    }

  vector<Double_t> finalFwdOffsets = getTruthFwdOffsets(TrueCathode,yValNewBest2,zValNewBest2);

  return finalFwdOffsets.at(0);
}

void doCalibFaces(const vector<trackInfo> &cosmicTracks, Int_t minTrackPoints, Int_t numTrackSegPoints, Double_t minWeight, Int_t minEntries, Int_t minMaxBin, Int_t saveInfo)
{
  const double bufferDist = 0.05;

  faceCalibHistTopDeltaX->Reset();
  faceCalibHistTopDeltaY->Reset();
  faceCalibHistTopDeltaZ->Reset();
    
  faceCalibHistBottomDeltaX->Reset();
  faceCalibHistBottomDeltaY->Reset();
  faceCalibHistBottomDeltaZ->Reset();

  faceCalibHistUpstreamDeltaX->Reset();
  faceCalibHistUpstreamDeltaY->Reset();
  faceCalibHistUpstreamDeltaZ->Reset();

  faceCalibHistDownstreamDeltaX->Reset();
  faceCalibHistDownstreamDeltaY->Reset();
  faceCalibHistDownstreamDeltaZ->Reset();

  faceCalibHistCathodeDeltaX->Reset();
  faceCalibHistCathodeDeltaY->Reset();
  faceCalibHistCathodeDeltaZ->Reset();

  TH1F *hists_Top_dX[nCalibDivisions_x+1][nCalibDivisions_z+1];
  TH1F *hists_Top_dY[nCalibDivisions_x+1][nCalibDivisions_z+1];
  TH1F *hists_Top_dZ[nCalibDivisions_x+1][nCalibDivisions_z+1];
    
  TH1F *hists_Bottom_dX[nCalibDivisions_x+1][nCalibDivisions_z+1];
  TH1F *hists_Bottom_dY[nCalibDivisions_x+1][nCalibDivisions_z+1];
  TH1F *hists_Bottom_dZ[nCalibDivisions_x+1][nCalibDivisions_z+1];

  TH1F *hists_Upstream_dX[nCalibDivisions_x+1][nCalibDivisions_y+1];
  TH1F *hists_Upstream_dY[nCalibDivisions_x+1][nCalibDivisions_y+1];
  TH1F *hists_Upstream_dZ[nCalibDivisions_x+1][nCalibDivisions_y+1];
    
  TH1F *hists_Downstream_dX[nCalibDivisions_x+1][nCalibDivisions_y+1];
  TH1F *hists_Downstream_dY[nCalibDivisions_x+1][nCalibDivisions_y+1];
  TH1F *hists_Downstream_dZ[nCalibDivisions_x+1][nCalibDivisions_y+1];

  TH1F *hists_Cathode_dX[nCalibDivisions_y+1][nCalibDivisions_z+1];
  TH1F *hists_Cathode_dY[nCalibDivisions_y+1][nCalibDivisions_z+1];
  TH1F *hists_Cathode_dZ[nCalibDivisions_y+1][nCalibDivisions_z+1];
    
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      calibTopDeltaX[x][z] = 0.0;
      calibTopDeltaY[x][z] = 0.0;
      calibTopDeltaZ[x][z] = 0.0;
  
      TH1F *tempHist_Top_dX = new TH1F(Form("tempHist_Top_dX_%d_%d",x,z),"",251,-0.251,0.251);
      TH1F *tempHist_Top_dY = new TH1F(Form("tempHist_Top_dY_%d_%d",x,z),"",251,-0.251,0.251);
      TH1F *tempHist_Top_dZ = new TH1F(Form("tempHist_Top_dZ_%d_%d",x,z),"",251,-0.251,0.251);
      tempHist_Top_dX->SetDirectory(0);
      tempHist_Top_dY->SetDirectory(0);
      tempHist_Top_dZ->SetDirectory(0);
      hists_Top_dX[x][z] = tempHist_Top_dX;
      hists_Top_dY[x][z] = tempHist_Top_dY;
      hists_Top_dZ[x][z] = tempHist_Top_dZ;
      
      calibBottomDeltaX[x][z] = 0.0;
      calibBottomDeltaY[x][z] = 0.0;
      calibBottomDeltaZ[x][z] = 0.0;
  
      TH1F *tempHist_Bottom_dX = new TH1F(Form("tempHist_Bottom_dX_%d_%d",x,z),"",251,-0.251,0.251);
      TH1F *tempHist_Bottom_dY = new TH1F(Form("tempHist_Bottom_dY_%d_%d",x,z),"",251,-0.251,0.251);
      TH1F *tempHist_Bottom_dZ = new TH1F(Form("tempHist_Bottom_dZ_%d_%d",x,z),"",251,-0.251,0.251);
      tempHist_Bottom_dX->SetDirectory(0);
      tempHist_Bottom_dY->SetDirectory(0);
      tempHist_Bottom_dZ->SetDirectory(0);
      hists_Bottom_dX[x][z] = tempHist_Bottom_dX;
      hists_Bottom_dY[x][z] = tempHist_Bottom_dY;
      hists_Bottom_dZ[x][z] = tempHist_Bottom_dZ;
    }
  }
  
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      calibUpstreamDeltaX[x][y] = 0.0;
      calibUpstreamDeltaY[x][y] = 0.0;
      calibUpstreamDeltaZ[x][y] = 0.0;
  
      TH1F *tempHist_Upstream_dX = new TH1F(Form("tempHist_Upstream_dX_%d_%d",x,y),"",251,-0.251,0.251);
      TH1F *tempHist_Upstream_dY = new TH1F(Form("tempHist_Upstream_dY_%d_%d",x,y),"",251,-0.251,0.251);
      TH1F *tempHist_Upstream_dZ = new TH1F(Form("tempHist_Upstream_dZ_%d_%d",x,y),"",251,-0.251,0.251);
      tempHist_Upstream_dX->SetDirectory(0);
      tempHist_Upstream_dY->SetDirectory(0);
      tempHist_Upstream_dZ->SetDirectory(0);
      hists_Upstream_dX[x][y] = tempHist_Upstream_dX;
      hists_Upstream_dY[x][y] = tempHist_Upstream_dY;
      hists_Upstream_dZ[x][y] = tempHist_Upstream_dZ;
  
      calibDownstreamDeltaX[x][y] = 0.0;
      calibDownstreamDeltaY[x][y] = 0.0;
      calibDownstreamDeltaZ[x][y] = 0.0;
  
      TH1F *tempHist_Downstream_dX = new TH1F(Form("tempHist_Downstream_dX_%d_%d",x,y),"",251,-0.251,0.251);
      TH1F *tempHist_Downstream_dY = new TH1F(Form("tempHist_Downstream_dY_%d_%d",x,y),"",251,-0.251,0.251);
      TH1F *tempHist_Downstream_dZ = new TH1F(Form("tempHist_Downstream_dZ_%d_%d",x,y),"",251,-0.251,0.251);
      tempHist_Downstream_dX->SetDirectory(0);
      tempHist_Downstream_dY->SetDirectory(0);
      tempHist_Downstream_dZ->SetDirectory(0);
      hists_Downstream_dX[x][y] = tempHist_Downstream_dX;
      hists_Downstream_dY[x][y] = tempHist_Downstream_dY;
      hists_Downstream_dZ[x][y] = tempHist_Downstream_dZ;
    }
  }
  
  for(Int_t y = 0; y <= nCalibDivisions_y; y++)
  {
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      calibCathodeDeltaX[y][z] = 0.0;
      calibCathodeDeltaY[y][z] = 0.0;
      calibCathodeDeltaZ[y][z] = 0.0;
  
      TH1F *tempHist_Cathode_dX = new TH1F(Form("tempHist_Cathode_dX_%d_%d",y,z),"",251,-0.251,0.251);
      TH1F *tempHist_Cathode_dY = new TH1F(Form("tempHist_Cathode_dY_%d_%d",y,z),"",251,-0.251,0.251);
      TH1F *tempHist_Cathode_dZ = new TH1F(Form("tempHist_Cathode_dZ_%d_%d",y,z),"",251,-0.251,0.251);
      tempHist_Cathode_dX->SetDirectory(0);
      tempHist_Cathode_dY->SetDirectory(0);
      tempHist_Cathode_dZ->SetDirectory(0);
      hists_Cathode_dX[y][z] = tempHist_Cathode_dX;
      hists_Cathode_dY[y][z] = tempHist_Cathode_dY;
      hists_Cathode_dZ[y][z] = tempHist_Cathode_dZ;
    }
  }
  
  vector<calibTrackInfo> cosmicCalibTracks = makeCalibTracks(cosmicTracks);
  
  calibTrackInfo calibTrack;
  Int_t numCosmicTracks = cosmicCalibTracks.size();
  
  Double_t xVal;
  Double_t yVal;
  Double_t zVal;
  Double_t xValDistorted;
  Double_t yValDistorted;
  Double_t zValDistorted;
  Int_t xCalibLowIndex;
  Int_t xCalibHighIndex;
  Int_t yCalibLowIndex;
  Int_t yCalibHighIndex;
  Int_t zCalibLowIndex;
  Int_t zCalibHighIndex;
  Double_t xCalibFrac;
  Double_t yCalibFrac;
  Double_t zCalibFrac;
  Double_t tempFactor;
  Double_t startX;
  Double_t startY;
  Double_t startZ;
  Double_t theta;
  Double_t phi;
  
  Int_t trackNum;
  Int_t numElec;
  Int_t whichFace;
    
  TTree *T_faceOffsets = new TTree("SpaCEtree_faceOffsets","SpaCEtree_faceOffsets");
  T_faceOffsets->Branch("trackNum",&trackNum,"data_faceOffsets/I");
  T_faceOffsets->Branch("faceNum",&whichFace,"data_faceOffsets/I");
  T_faceOffsets->Branch("numElec",&numElec,"data_faceOffsets/I");
  T_faceOffsets->Branch("xVal",&xVal,"data_faceOffsets/D");
  T_faceOffsets->Branch("yVal",&yVal,"data_faceOffsets/D");
  T_faceOffsets->Branch("zVal",&zVal,"data_faceOffsets/D");
  T_faceOffsets->Branch("xValDistorted",&xValDistorted,"data_faceOffsets/D");
  T_faceOffsets->Branch("yValDistorted",&yValDistorted,"data_faceOffsets/D");
  T_faceOffsets->Branch("zValDistorted",&zValDistorted,"data_faceOffsets/D");
  T_faceOffsets->Branch("startX",&startX,"data_faceOffsets/D");
  T_faceOffsets->Branch("startY",&startY,"data_faceOffsets/D");
  T_faceOffsets->Branch("startZ",&startZ,"data_faceOffsets/D");
  T_faceOffsets->Branch("theta",&theta,"data_faceOffsets/D");
  T_faceOffsets->Branch("phi",&phi,"data_faceOffsets/D");
  T_faceOffsets->Branch("weight",&tempFactor,"data_faceOffsets/D");
  if(saveInfo == 1) {
    T_faceOffsets->SetDirectory(outputFile);
  }
  
  for(Int_t i = 0; i < numCosmicTracks; i++)
  {
    cout << "COSMIC-FACE " << i << endl;
  
    calibTrack = cosmicCalibTracks.at(i);
    if(calibTrack.track.electrons.size() < minTrackPoints) continue;
  
    trackNum = i;
    numElec = calibTrack.track.electrons.size();
    theta = calibTrack.track.theta;
    phi = calibTrack.track.phi;
    
    Int_t numNearAnode = 0;
    if (calibTrack.track.electrons.at(0).x_mod > (TrueAnode - maxXdist)) {
      numNearAnode++;
    }
    if (calibTrack.track.electrons.at(numElec-1).x_mod > (TrueAnode - maxXdist)) {
      numNearAnode++;
    }
    if((numNearAnode == 0) || (numNearAnode == 2)) continue;
  	
    PointCloud calibPoints;
    
    Double_t avgX1 = 0.0;
    Double_t avgX2 = 0.0;
  
    for(Int_t j = 0; j < 5; j++)
    {
      Double_t elecX = calibTrack.track.electrons.at(j).x_mod;
      avgX1 += elecX;
    }
    avgX1 /= 5.0;
  
    for(Int_t j = numElec-1; j > numElec-6; j--)
    {
      Double_t elecX = calibTrack.track.electrons.at(j).x_mod;
      avgX2 += elecX;
    }
    avgX2 /= 5.0;
  
    double SCEfactor = 1.0;
    if (isSCEon == false) {
      SCEfactor = 0.0;
    }
    
    Int_t numBadPoints = 0;
    if(avgX1 < avgX2) {
      xValDistorted = calibTrack.track.electrons.at(0).x_mod;
      yValDistorted = calibTrack.track.electrons.at(0).y_mod;
      zValDistorted = calibTrack.track.electrons.at(0).z_mod;
  
      for (int j = numElec-1; j > numElec-1-numTrackSegPoints; j--) {
        if ((calibTrack.track.electrons.at(j).x_mod < TrueCathode-bufferDist) || (calibTrack.track.electrons.at(j).x_mod > TrueAnode+bufferDist) || (calibTrack.track.electrons.at(j).y_mod < TrueBottom-bufferDist) || (calibTrack.track.electrons.at(j).y_mod > TrueTop+bufferDist) || (calibTrack.track.electrons.at(j).z_mod < TrueUpstream-bufferDist) || (calibTrack.track.electrons.at(j).z_mod > TrueDownstream+bufferDist)) {
          numBadPoints++;
  	  continue;
	}
  
  	Point tempPoint;
        //tempPoint.x = calibTrack.track.electrons.at(j).x_mod;
        //tempPoint.y = calibTrack.track.electrons.at(j).y_mod;
        //tempPoint.z = calibTrack.track.electrons.at(j).z_mod;
        //tempPoint.x = calibTrack.track.electrons.at(j).x_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,1); // TEMP USE OF TRUTH INFO
        //tempPoint.y = calibTrack.track.electrons.at(j).y_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,2); // TEMP USE OF TRUTH INFO
        //tempPoint.z = calibTrack.track.electrons.at(j).z_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,3); // TEMP USE OF TRUTH INFO
        tempPoint.x = calibTrack.track.electrons.at(j).x_mod + getBulkCorr(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,1);
        tempPoint.y = calibTrack.track.electrons.at(j).y_mod + getBulkCorr(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,2);
        tempPoint.z = calibTrack.track.electrons.at(j).z_mod + getBulkCorr(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,3);
  
        calibPoints.push_back(tempPoint);
      }
    }
    else {
      xValDistorted = calibTrack.track.electrons.at(numElec-1).x_mod;
      yValDistorted = calibTrack.track.electrons.at(numElec-1).y_mod;
      zValDistorted = calibTrack.track.electrons.at(numElec-1).z_mod;
  
      for (int j = 0; j < numTrackSegPoints; j++) {
  	if ((calibTrack.track.electrons.at(j).x_mod < TrueCathode-bufferDist) || (calibTrack.track.electrons.at(j).x_mod > TrueAnode+bufferDist) || (calibTrack.track.electrons.at(j).y_mod < TrueBottom-bufferDist) || (calibTrack.track.electrons.at(j).y_mod > TrueTop+bufferDist) || (calibTrack.track.electrons.at(j).z_mod < TrueUpstream-bufferDist) || (calibTrack.track.electrons.at(j).z_mod > TrueDownstream+bufferDist)) {
          numBadPoints++;
  	  continue;
  	}
  	
  	Point tempPoint;
        //tempPoint.x = calibTrack.track.electrons.at(j).x_mod;
        //tempPoint.y = calibTrack.track.electrons.at(j).y_mod;
        //tempPoint.z = calibTrack.track.electrons.at(j).z_mod;
        //tempPoint.x = calibTrack.track.electrons.at(j).x_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,1); // TEMP USE OF TRUTH INFO
        //tempPoint.y = calibTrack.track.electrons.at(j).y_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,2); // TEMP USE OF TRUTH INFO
        //tempPoint.z = calibTrack.track.electrons.at(j).z_mod + SCEfactor*getTruthOffset(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,3); // TEMP USE OF TRUTH INFO
        tempPoint.x = calibTrack.track.electrons.at(j).x_mod + getBulkCorr(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,1);
        tempPoint.y = calibTrack.track.electrons.at(j).y_mod + getBulkCorr(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,2);
        tempPoint.z = calibTrack.track.electrons.at(j).z_mod + getBulkCorr(calibTrack.track.electrons.at(j).x_mod,calibTrack.track.electrons.at(j).y_mod,calibTrack.track.electrons.at(j).z_mod,3);
  
        calibPoints.push_back(tempPoint);
      }
    }
    if (numBadPoints > 5) continue;
    
    PCAResults results = DoPCA(calibPoints);
  
    if (results.endPoints.first(0) > results.endPoints.second(0)) {
      startX = results.endPoints.first(0);
      startY = results.endPoints.first(1);
      startZ = results.endPoints.first(2);
    }
    else {
      startX = results.endPoints.second(0);
      startY = results.endPoints.second(1);
      startZ = results.endPoints.second(2);
    }
  
    Double_t unitX = results.eVecs[0](0);
    Double_t unitY = results.eVecs[0](1);
    Double_t unitZ = results.eVecs[0](2);
    if (unitX > 0.0) {
      unitX *= -1.0;
      unitY *= -1.0;
      unitZ *= -1.0;
    }
  
    whichFace = -1;
    
    Double_t t_Top = (TrueTop-startY)/unitY;
    Double_t x_Top = unitX*t_Top + startX;
    Double_t y_Top = TrueTop;
    Double_t z_Top = unitZ*t_Top + startZ;
    if((t_Top > 0.0) && (x_Top > TrueCathode) && (x_Top < TrueAnode) && (z_Top > TrueUpstream) && (z_Top < TrueDownstream)) {
      xVal = x_Top;
      yVal = y_Top;
      zVal = z_Top;
      whichFace = 0;
    }
  
    Double_t t_Bottom = (TrueBottom-startY)/unitY;
    Double_t x_Bottom = unitX*t_Bottom + startX;
    Double_t y_Bottom = TrueBottom;
    Double_t z_Bottom = unitZ*t_Bottom + startZ;
    if((t_Bottom > 0.0) && (x_Bottom > TrueCathode) && (x_Bottom < TrueAnode) && (z_Bottom > TrueUpstream) && (z_Bottom < TrueDownstream)) {
      xVal = x_Bottom;
      yVal = y_Bottom;
      zVal = z_Bottom;
      whichFace = 1;
    }
  
    Double_t t_Upstream = (TrueUpstream-startZ)/unitZ;
    Double_t x_Upstream = unitX*t_Upstream + startX;
    Double_t y_Upstream = unitY*t_Upstream + startY;
    Double_t z_Upstream = TrueUpstream;
    if((t_Upstream > 0.0) && (x_Upstream > TrueCathode) && (x_Upstream < TrueAnode) && (y_Upstream > TrueBottom) && (y_Upstream < TrueTop)) {
      xVal = x_Upstream;
      yVal = y_Upstream;
      zVal = z_Upstream;
      whichFace = 2;
    }
  
    Double_t t_Downstream = (TrueDownstream-startZ)/unitZ;
    Double_t x_Downstream = unitX*t_Downstream + startX;
    Double_t y_Downstream = unitY*t_Downstream + startY;
    Double_t z_Downstream = TrueDownstream;
    if((t_Downstream > 0.0) && (x_Downstream > TrueCathode) && (x_Downstream < TrueAnode) && (y_Downstream > TrueBottom) && (y_Downstream < TrueTop)) {
      xVal = x_Downstream;
      yVal = y_Downstream;
      zVal = z_Downstream;
      whichFace = 3;
    }
  
    Double_t t_Cathode = (TrueCathode-startX)/unitX;
    Double_t x_Cathode = TrueCathode;
    Double_t y_Cathode = unitY*t_Cathode + startY;
    Double_t z_Cathode = unitZ*t_Cathode + startZ;
    if((t_Cathode > 0.0) && (y_Cathode > TrueBottom) && (y_Cathode < TrueTop) && (z_Cathode > TrueUpstream) && (z_Cathode < TrueDownstream)) {
      xVal = x_Cathode;
      yVal = y_Cathode;
      zVal = z_Cathode;
      whichFace = 4;
    }
    
    xCalibLowIndex = TMath::Floor((xValDistorted/Lx)*nCalibDivisions_x);
    xCalibHighIndex = TMath::Ceil((xValDistorted/Lx)*nCalibDivisions_x);
    yCalibLowIndex = TMath::Floor((yValDistorted/Ly)*nCalibDivisions_y);
    yCalibHighIndex = TMath::Ceil((yValDistorted/Ly)*nCalibDivisions_y);
    zCalibLowIndex = TMath::Floor((zValDistorted/Lz)*nCalibDivisions_z);
    zCalibHighIndex = TMath::Ceil((zValDistorted/Lz)*nCalibDivisions_z);
    
    if(xValDistorted < 0.0) {
      xCalibLowIndex = 0;
      xCalibHighIndex = 1;
      xCalibFrac = 0.0;
    }
    else if (xValDistorted > Lx) {
      xCalibLowIndex = nCalibDivisions_x - 1;
      xCalibHighIndex = nCalibDivisions_x;
      xCalibFrac = 1.0;
    }
    
    if(yValDistorted < 0.0) {
      yCalibLowIndex = 0;
      yCalibHighIndex = 1;
      yCalibFrac = 0.0;
    }
    else if (yValDistorted > Ly) {
      yCalibLowIndex = nCalibDivisions_y - 1;
      yCalibHighIndex = nCalibDivisions_y;
      yCalibFrac = 1.0;
    }
    
    if(zValDistorted < 0.0) {
      zCalibLowIndex = 0;
      zCalibHighIndex = 1;
      zCalibFrac = 0.0;
    }
    else if (zValDistorted > Lz) {
      zCalibLowIndex = nCalibDivisions_z - 1;
      zCalibHighIndex = nCalibDivisions_z;
      zCalibFrac = 1.0;
    }
    
    xCalibFrac = ((xValDistorted/Lx)*nCalibDivisions_x)-((Double_t) xCalibLowIndex);
    yCalibFrac = ((yValDistorted/Ly)*nCalibDivisions_y)-((Double_t) yCalibLowIndex);
    zCalibFrac = ((zValDistorted/Lz)*nCalibDivisions_z)-((Double_t) zCalibLowIndex);
    
    if (whichFace == 0) {
      tempFactor = (1.0-xCalibFrac)*(1.0-zCalibFrac);
      if(tempFactor >= minWeight) {
        hists_Top_dX[xCalibLowIndex][zCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Top_dY[xCalibLowIndex][zCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Top_dZ[xCalibLowIndex][zCalibLowIndex]->Fill(zVal-zValDistorted);
      }
  
      tempFactor = (1.0-xCalibFrac)*zCalibFrac;
      if(tempFactor >= minWeight) {
        hists_Top_dX[xCalibLowIndex][zCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Top_dY[xCalibLowIndex][zCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Top_dZ[xCalibLowIndex][zCalibHighIndex]->Fill(zVal-zValDistorted);
      }
  
      tempFactor = xCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      if(tempFactor >= minWeight) {
        hists_Top_dX[xCalibHighIndex][zCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Top_dY[xCalibHighIndex][zCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Top_dZ[xCalibHighIndex][zCalibLowIndex]->Fill(zVal-zValDistorted);
      }
  
      tempFactor = xCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      if(tempFactor >= minWeight) {
        hists_Top_dX[xCalibHighIndex][zCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Top_dY[xCalibHighIndex][zCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Top_dZ[xCalibHighIndex][zCalibHighIndex]->Fill(zVal-zValDistorted);
      }
    }
    else if (whichFace == 1) {
      tempFactor = (1.0-xCalibFrac)*(1.0-zCalibFrac);
      if(tempFactor >= minWeight) {
        hists_Bottom_dX[xCalibLowIndex][zCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Bottom_dY[xCalibLowIndex][zCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Bottom_dZ[xCalibLowIndex][zCalibLowIndex]->Fill(zVal-zValDistorted);
      }

      tempFactor = (1.0-xCalibFrac)*zCalibFrac;
      if(tempFactor >= minWeight) {
        hists_Bottom_dX[xCalibLowIndex][zCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Bottom_dY[xCalibLowIndex][zCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Bottom_dZ[xCalibLowIndex][zCalibHighIndex]->Fill(zVal-zValDistorted);
      }

      tempFactor = xCalibFrac*(1.0-zCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      if(tempFactor >= minWeight) {
        hists_Bottom_dX[xCalibHighIndex][zCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Bottom_dY[xCalibHighIndex][zCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Bottom_dZ[xCalibHighIndex][zCalibLowIndex]->Fill(zVal-zValDistorted);
      }

      tempFactor = xCalibFrac*zCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      if(tempFactor >= minWeight) {
        hists_Bottom_dX[xCalibHighIndex][zCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Bottom_dY[xCalibHighIndex][zCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Bottom_dZ[xCalibHighIndex][zCalibHighIndex]->Fill(zVal-zValDistorted);
      }
    }
    else if (whichFace == 2) {
      tempFactor = (1.0-xCalibFrac)*(1.0-yCalibFrac);
      if(tempFactor >= minWeight) {
        hists_Upstream_dX[xCalibLowIndex][yCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Upstream_dY[xCalibLowIndex][yCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Upstream_dZ[xCalibLowIndex][yCalibLowIndex]->Fill(zVal-zValDistorted);
      }

      tempFactor = (1.0-xCalibFrac)*yCalibFrac;
      if(tempFactor >= minWeight) {
        hists_Upstream_dX[xCalibLowIndex][yCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Upstream_dY[xCalibLowIndex][yCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Upstream_dZ[xCalibLowIndex][yCalibHighIndex]->Fill(zVal-zValDistorted);
      }

      tempFactor = xCalibFrac*(1.0-yCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      if(tempFactor >= minWeight) {
        hists_Upstream_dX[xCalibHighIndex][yCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Upstream_dY[xCalibHighIndex][yCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Upstream_dZ[xCalibHighIndex][yCalibLowIndex]->Fill(zVal-zValDistorted);
      }

      tempFactor = xCalibFrac*yCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      if(tempFactor >= minWeight) {
        hists_Upstream_dX[xCalibHighIndex][yCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Upstream_dY[xCalibHighIndex][yCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Upstream_dZ[xCalibHighIndex][yCalibHighIndex]->Fill(zVal-zValDistorted);
      }
    }
    else if (whichFace == 3) {
      tempFactor = (1.0-xCalibFrac)*(1.0-yCalibFrac);
      if(tempFactor >= minWeight) {
        hists_Downstream_dX[xCalibLowIndex][yCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Downstream_dY[xCalibLowIndex][yCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Downstream_dZ[xCalibLowIndex][yCalibLowIndex]->Fill(zVal-zValDistorted);
      }

      tempFactor = (1.0-xCalibFrac)*yCalibFrac;
      if(tempFactor >= minWeight) {
        hists_Downstream_dX[xCalibLowIndex][yCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Downstream_dY[xCalibLowIndex][yCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Downstream_dZ[xCalibLowIndex][yCalibHighIndex]->Fill(zVal-zValDistorted);
      }

      tempFactor = xCalibFrac*(1.0-yCalibFrac);
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      if(tempFactor >= minWeight) {
        hists_Downstream_dX[xCalibHighIndex][yCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Downstream_dY[xCalibHighIndex][yCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Downstream_dZ[xCalibHighIndex][yCalibLowIndex]->Fill(zVal-zValDistorted);
      }
  
      tempFactor = xCalibFrac*yCalibFrac;
      if(xCalibHighIndex == nCalibDivisions_x)
        tempFactor = 0.0;
      if(tempFactor >= minWeight) {
        hists_Downstream_dX[xCalibHighIndex][yCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Downstream_dY[xCalibHighIndex][yCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Downstream_dZ[xCalibHighIndex][yCalibHighIndex]->Fill(zVal-zValDistorted);
      }
    }
    else if (whichFace == 4) {
      tempFactor = (1.0-yCalibFrac)*(1.0-zCalibFrac);
      if(tempFactor >= minWeight) {
        hists_Cathode_dX[yCalibLowIndex][zCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Cathode_dY[yCalibLowIndex][zCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Cathode_dZ[yCalibLowIndex][zCalibLowIndex]->Fill(zVal-zValDistorted);
      }
  
      tempFactor = (1.0-yCalibFrac)*zCalibFrac;
      if(tempFactor >= minWeight) {
        hists_Cathode_dX[yCalibLowIndex][zCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Cathode_dY[yCalibLowIndex][zCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Cathode_dZ[yCalibLowIndex][zCalibHighIndex]->Fill(zVal-zValDistorted);
      }
  
      tempFactor = yCalibFrac*(1.0-zCalibFrac);
      if(tempFactor >= minWeight) {
        hists_Cathode_dX[yCalibHighIndex][zCalibLowIndex]->Fill(xVal-xValDistorted);
        hists_Cathode_dY[yCalibHighIndex][zCalibLowIndex]->Fill(yVal-yValDistorted);
        hists_Cathode_dZ[yCalibHighIndex][zCalibLowIndex]->Fill(zVal-zValDistorted);
      }
  
      tempFactor = yCalibFrac*zCalibFrac;
      if(tempFactor >= minWeight) {
        hists_Cathode_dX[yCalibHighIndex][zCalibHighIndex]->Fill(xVal-xValDistorted);
        hists_Cathode_dY[yCalibHighIndex][zCalibHighIndex]->Fill(yVal-yValDistorted);
        hists_Cathode_dZ[yCalibHighIndex][zCalibHighIndex]->Fill(zVal-zValDistorted);
      }
    }
  
    if (whichFace >= 0) {
      cout << "  " << whichFace << " " << xVal << " " << yVal << " " << zVal << " " << xValDistorted << " " << yValDistorted << " " << zValDistorted << " " << endl;
      if(saveInfo == 1) {
        T_faceOffsets->Fill();
      }
    }
  }
  
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      if(hists_Top_dX[x][z]->GetEntries() >= minEntries)
      {
	if(hists_Top_dX[x][z]->GetMaximum() >= minMaxBin)
          calibTopDeltaX[x][z] = hists_Top_dX[x][z]->GetBinCenter(hists_Top_dX[x][z]->GetMaximumBin());
	else
          calibTopDeltaX[x][z] = findHistMedian(hists_Top_dX[x][z]);
      }
      else
      {
        calibTopDeltaX[x][z] = -10000.0;
      }

      if(hists_Top_dY[x][z]->GetEntries() >= minEntries)
      {
	if(hists_Top_dY[x][z]->GetMaximum() >= minMaxBin)
          calibTopDeltaY[x][z] = hists_Top_dY[x][z]->GetBinCenter(hists_Top_dY[x][z]->GetMaximumBin());
	else
          calibTopDeltaY[x][z] = findHistMedian(hists_Top_dY[x][z]);
      }
      else
      {
        calibTopDeltaY[x][z] = -10000.0;
      }
      
      if(hists_Top_dZ[x][z]->GetEntries() >= minEntries)
      {
	if(hists_Top_dZ[x][z]->GetMaximum() >= minMaxBin)
          calibTopDeltaZ[x][z] = hists_Top_dZ[x][z]->GetBinCenter(hists_Top_dZ[x][z]->GetMaximumBin());
	else
          calibTopDeltaZ[x][z] = findHistMedian(hists_Top_dZ[x][z]);
      }
      else
      {
        calibTopDeltaZ[x][z] = -10000.0;
      }

      if(hists_Bottom_dX[x][z]->GetEntries() >= minEntries)
      {
	if(hists_Bottom_dX[x][z]->GetMaximum() >= minMaxBin)
          calibBottomDeltaX[x][z] = hists_Bottom_dX[x][z]->GetBinCenter(hists_Bottom_dX[x][z]->GetMaximumBin());
	else
          calibBottomDeltaX[x][z] = findHistMedian(hists_Bottom_dX[x][z]);
      }
      else
      {
        calibBottomDeltaX[x][z] = -10000.0;
      }

      if(hists_Bottom_dY[x][z]->GetEntries() >= minEntries)
      {
	if(hists_Bottom_dY[x][z]->GetMaximum() >= minMaxBin)
          calibBottomDeltaY[x][z] = hists_Bottom_dY[x][z]->GetBinCenter(hists_Bottom_dY[x][z]->GetMaximumBin());
	else
          calibBottomDeltaY[x][z] = findHistMedian(hists_Bottom_dY[x][z]);
      }
      else
      {
        calibBottomDeltaY[x][z] = -10000.0;
      }

      if(hists_Bottom_dZ[x][z]->GetEntries() >= minEntries)
      {
	if(hists_Bottom_dZ[x][z]->GetMaximum() >= minMaxBin)
          calibBottomDeltaZ[x][z] = hists_Bottom_dZ[x][z]->GetBinCenter(hists_Bottom_dZ[x][z]->GetMaximumBin());
	else
          calibBottomDeltaZ[x][z] = findHistMedian(hists_Bottom_dZ[x][z]);
      }
      else
      {
        calibBottomDeltaZ[x][z] = -10000.0;
      }
    }
  }
  
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      if(hists_Upstream_dX[x][y]->GetEntries() >= minEntries)
      {
	if(hists_Upstream_dX[x][y]->GetMaximum() >= minMaxBin)
          calibUpstreamDeltaX[x][y] = hists_Upstream_dX[x][y]->GetBinCenter(hists_Upstream_dX[x][y]->GetMaximumBin());
	else
          calibUpstreamDeltaX[x][y] = findHistMedian(hists_Upstream_dX[x][y]);
      }
      else
      {
        calibUpstreamDeltaX[x][y] = -10000.0;
      }

      if(hists_Upstream_dY[x][y]->GetEntries() >= minEntries)
      {
	if(hists_Upstream_dY[x][y]->GetMaximum() >= minMaxBin)
          calibUpstreamDeltaY[x][y] = hists_Upstream_dY[x][y]->GetBinCenter(hists_Upstream_dY[x][y]->GetMaximumBin());
	else
          calibUpstreamDeltaY[x][y] = findHistMedian(hists_Upstream_dY[x][y]);
      }
      else
      {
        calibUpstreamDeltaY[x][y] = -10000.0;
      }

      if(hists_Upstream_dZ[x][y]->GetEntries() >= minEntries)
      {
	if(hists_Upstream_dZ[x][y]->GetMaximum() >= minMaxBin)
          calibUpstreamDeltaZ[x][y] = hists_Upstream_dZ[x][y]->GetBinCenter(hists_Upstream_dZ[x][y]->GetMaximumBin());
	else
          calibUpstreamDeltaZ[x][y] = findHistMedian(hists_Upstream_dZ[x][y]);
      }
      else
      {
        calibUpstreamDeltaZ[x][y] = -10000.0;
      }

      if(hists_Downstream_dX[x][y]->GetEntries() >= minEntries)
      {
	if(hists_Downstream_dX[x][y]->GetMaximum() >= minMaxBin)
          calibDownstreamDeltaX[x][y] = hists_Downstream_dX[x][y]->GetBinCenter(hists_Downstream_dX[x][y]->GetMaximumBin());
	else
          calibDownstreamDeltaX[x][y] = findHistMedian(hists_Downstream_dX[x][y]);
      }
      else
      {
        calibDownstreamDeltaX[x][y] = -10000.0;
      }

      if(hists_Downstream_dY[x][y]->GetEntries() >= minEntries)
      {
	if(hists_Downstream_dY[x][y]->GetMaximum() >= minMaxBin)
          calibDownstreamDeltaY[x][y] = hists_Downstream_dY[x][y]->GetBinCenter(hists_Downstream_dY[x][y]->GetMaximumBin());
	else
          calibDownstreamDeltaY[x][y] = findHistMedian(hists_Downstream_dY[x][y]);
      }
      else
      {
        calibDownstreamDeltaY[x][y] = -10000.0;
      }

      if(hists_Downstream_dZ[x][y]->GetEntries() >= minEntries)
      {
	if(hists_Downstream_dZ[x][y]->GetMaximum() >= minMaxBin)
          calibDownstreamDeltaZ[x][y] = hists_Downstream_dZ[x][y]->GetBinCenter(hists_Downstream_dZ[x][y]->GetMaximumBin());
	else
          calibDownstreamDeltaZ[x][y] = findHistMedian(hists_Downstream_dZ[x][y]);
      }
      else
      {
        calibDownstreamDeltaZ[x][y] = -10000.0;
      }
    }
  }
  
  for(Int_t y = 0; y <= nCalibDivisions_y; y++)
  {
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      if(hists_Cathode_dX[y][z]->GetEntries() >= minEntries)
      {
	if(hists_Cathode_dX[y][z]->GetMaximum() >= minMaxBin)
          calibCathodeDeltaX[y][z] = hists_Cathode_dX[y][z]->GetBinCenter(hists_Cathode_dX[y][z]->GetMaximumBin());
	else
          calibCathodeDeltaX[y][z] = findHistMedian(hists_Cathode_dX[y][z]);
      }
      else
      {
        calibCathodeDeltaX[y][z] = -10000.0;
      }

      if(hists_Cathode_dY[y][z]->GetEntries() >= minEntries)
      {
	if(hists_Cathode_dY[y][z]->GetMaximum() >= minMaxBin)
          calibCathodeDeltaY[y][z] = hists_Cathode_dY[y][z]->GetBinCenter(hists_Cathode_dY[y][z]->GetMaximumBin());
	else
          calibCathodeDeltaY[y][z] = findHistMedian(hists_Cathode_dY[y][z]);
      }
      else
      {
        calibCathodeDeltaY[y][z] = -10000.0;
      }

      if(hists_Cathode_dZ[y][z]->GetEntries() >= minEntries)
      {
	if(hists_Cathode_dZ[y][z]->GetMaximum() >= minMaxBin)
          calibCathodeDeltaZ[y][z] = hists_Cathode_dZ[y][z]->GetBinCenter(hists_Cathode_dZ[y][z]->GetMaximumBin());
	else
          calibCathodeDeltaZ[y][z] = findHistMedian(hists_Cathode_dZ[y][z]);
      }
      else
      {
        calibCathodeDeltaZ[y][z] = -10000.0;
      }
    }
  }

  conditionFaceMap(calibTopDeltaX,0);
  conditionFaceMap(calibTopDeltaY,0);
  conditionFaceMap(calibTopDeltaZ,0);
  conditionFaceMap(calibBottomDeltaX,1);
  conditionFaceMap(calibBottomDeltaY,1);
  conditionFaceMap(calibBottomDeltaZ,1);
  conditionFaceMap(calibUpstreamDeltaX,2);
  conditionFaceMap(calibUpstreamDeltaY,2);
  conditionFaceMap(calibUpstreamDeltaZ,2);
  conditionFaceMap(calibDownstreamDeltaX,3);
  conditionFaceMap(calibDownstreamDeltaY,3);
  conditionFaceMap(calibDownstreamDeltaZ,3);
  conditionFaceMap(calibCathodeDeltaX,4);
  conditionFaceMap(calibCathodeDeltaY,4);
  conditionFaceMap(calibCathodeDeltaZ,4);
    
  Double_t x_true, y_true, z_true;
  Double_t x_reco, y_reco, z_reco;
  Double_t Dx, Dy, Dz;
  Int_t numEntries;
  Int_t elecFate;
  Int_t faceNum;
  
  TTree *T_calibFaces = new TTree("SpaCEtree_calibFaces","SpaCEtree_calibFaces");
  T_calibFaces->Branch("x_true",&x_true,"data_calibFaces/D");
  T_calibFaces->Branch("y_true",&y_true,"data_calibFaces/D");
  T_calibFaces->Branch("z_true",&z_true,"data_calibFaces/D");
  T_calibFaces->Branch("x_reco",&x_reco,"data_calibFaces/D");
  T_calibFaces->Branch("y_reco",&y_reco,"data_calibFaces/D");
  T_calibFaces->Branch("z_reco",&z_reco,"data_calibFaces/D");
  T_calibFaces->Branch("Dx",&Dx,"data_calibFaces/D");
  T_calibFaces->Branch("Dy",&Dy,"data_calibFaces/D");
  T_calibFaces->Branch("Dz",&Dz,"data_calibFaces/D");
  T_calibFaces->Branch("numEntries",&numEntries,"data_calibFaces/I");
  T_calibFaces->Branch("elecFate",&elecFate,"data_calibFaces/I");
  T_calibFaces->Branch("faceNum",&faceNum,"data_calibFaces/I");
  if(saveInfo == 1) {
    T_calibFaces->SetDirectory(outputFile);
  }
  
  faceNum = 0;
  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    z_reco = -1.0*Lz/nCalibDivisions_z;
  
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      z_reco += Lz/nCalibDivisions_z;
  
      if(x == nCalibDivisions_x) {
        Dx = 0.0;
        Dy = 0.0;
        Dz = 0.0;
      }
      else {
        Dx = calibTopDeltaX[x][z];
        Dy = calibTopDeltaY[x][z];
        Dz = calibTopDeltaZ[x][z];
      }
  
      numEntries = hists_Top_dX[x][z]->GetEntries();
  
      x_true = x_reco+Dx;
      y_true = TrueTop;
      y_reco = y_true-Dy;
      z_true = z_reco+Dz;
  
      if((numEntries >= minEntries) || (x == nCalibDivisions_x))
        elecFate = 1;
      else
        elecFate = 0;

      faceCalibHistTopDeltaX->SetBinContent(z+1,x+1,Dx);
      faceCalibHistTopDeltaY->SetBinContent(z+1,x+1,Dy);
      faceCalibHistTopDeltaZ->SetBinContent(z+1,x+1,Dz);

      if(saveInfo == 1) {
        T_calibFaces->Fill();
      }
    }
  }

  faceNum = 1;
  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    z_reco = -1.0*Lz/nCalibDivisions_z;
  
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      z_reco += Lz/nCalibDivisions_z;
  
      if(x == nCalibDivisions_x) {
        Dx = 0.0;
        Dy = 0.0;
        Dz = 0.0;
      }
      else {
        Dx = calibBottomDeltaX[x][z];
        Dy = calibBottomDeltaY[x][z];
        Dz = calibBottomDeltaZ[x][z];
      }
  
      numEntries = hists_Bottom_dX[x][z]->GetEntries();
  
      x_true = x_reco+Dx;
      y_true = TrueBottom;
      y_reco = y_true-Dy;
      z_true = z_reco+Dz;
  
      if((numEntries >= minEntries) || (x == nCalibDivisions_x))
        elecFate = 1;
      else
        elecFate = 0;
    
      faceCalibHistBottomDeltaX->SetBinContent(z+1,x+1,Dx);
      faceCalibHistBottomDeltaY->SetBinContent(z+1,x+1,Dy);
      faceCalibHistBottomDeltaZ->SetBinContent(z+1,x+1,Dz);

      if(saveInfo == 1) {
        T_calibFaces->Fill();
      }
    }
  }
  
  faceNum = 2;
  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    y_reco = -1.0*Ly/nCalibDivisions_y;
  
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      y_reco += Ly/nCalibDivisions_y;
  
      if(x == nCalibDivisions_x) {
        Dx = 0.0;
        Dy = 0.0;
        Dz = 0.0;
      }
      else {
        Dx = calibUpstreamDeltaX[x][y];
        Dy = calibUpstreamDeltaY[x][y];
        Dz = calibUpstreamDeltaZ[x][y];
      }
  
      numEntries = hists_Upstream_dX[x][y]->GetEntries();
  
      x_true = x_reco+Dx;
      y_true = y_reco+Dy;
      z_true = TrueUpstream;
      z_reco = z_true-Dz;
  
      if((numEntries >= minEntries) || (x == nCalibDivisions_x))
        elecFate = 1;
      else
        elecFate = 0;
    
      faceCalibHistUpstreamDeltaX->SetBinContent(x+1,y+1,Dx);
      faceCalibHistUpstreamDeltaY->SetBinContent(x+1,y+1,Dy);
      faceCalibHistUpstreamDeltaZ->SetBinContent(x+1,y+1,Dz);

      if(saveInfo == 1) {
        T_calibFaces->Fill();
      }
    }
  }
  
  faceNum = 3;
  x_reco = -1.0*Lx/nCalibDivisions_x;
  for(Int_t x = 0; x <= nCalibDivisions_x; x++)
  {
    x_reco += Lx/nCalibDivisions_x;
    y_reco = -1.0*Ly/nCalibDivisions_y;
  
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      y_reco += Ly/nCalibDivisions_y;
  
      if(x == nCalibDivisions_x) {
        Dx = 0.0;
        Dy = 0.0;
        Dz = 0.0;
      }
      else {
        Dx = calibDownstreamDeltaX[x][y];
        Dy = calibDownstreamDeltaY[x][y];
        Dz = calibDownstreamDeltaZ[x][y];
      }
  
      numEntries = hists_Downstream_dX[x][y]->GetEntries();
  
      x_true = x_reco+Dx;
      y_true = y_reco+Dy;
      z_true = TrueDownstream;
      z_reco = z_true-Dz;
  
      if((numEntries >= minEntries) || (x == nCalibDivisions_x))
        elecFate = 1;
      else
        elecFate = 0;
    
      faceCalibHistDownstreamDeltaX->SetBinContent(x+1,y+1,Dx);
      faceCalibHistDownstreamDeltaY->SetBinContent(x+1,y+1,Dy);
      faceCalibHistDownstreamDeltaZ->SetBinContent(x+1,y+1,Dz);

      if(saveInfo == 1) {
        T_calibFaces->Fill();
      }
    }
  }
  
  faceNum = 4;
  y_reco = -1.0*Ly/nCalibDivisions_y;
  for(Int_t y = 0; y <= nCalibDivisions_y; y++)
  {
    y_reco += Ly/nCalibDivisions_y;
    z_reco = -1.0*Lz/nCalibDivisions_z;
  
    for(Int_t z = 0; z <= nCalibDivisions_z; z++)
    {
      z_reco += Lz/nCalibDivisions_z;
  
      Dx = calibCathodeDeltaX[y][z];
      Dy = calibCathodeDeltaY[y][z];
      Dz = calibCathodeDeltaZ[y][z];
  
      numEntries = hists_Cathode_dX[y][z]->GetEntries();
  
      x_true = TrueCathode;
      x_reco = x_true-Dx;
      y_true = y_reco+Dy;
      z_true = z_reco+Dz;
  
      if(numEntries >= minEntries)
        elecFate = 1;
      else
        elecFate = 0;
    
      faceCalibHistCathodeDeltaX->SetBinContent(z+1,y+1,Dx);
      faceCalibHistCathodeDeltaY->SetBinContent(z+1,y+1,Dy);
      faceCalibHistCathodeDeltaZ->SetBinContent(z+1,y+1,Dz);

      if(saveInfo == 1) {
        T_calibFaces->Fill();
      }
    }
  }
  
  return;
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

Double_t findHistMedian(TH1F *hist)
{
  Double_t q = 0.5;
  Double_t x = 0.0;
  hist->ComputeIntegral();
  hist->GetQuantiles(1, &x, &q);

  return x;
}

void conditionFaceMap(vector<vector<Double_t> > &inputMap, Int_t faceNum)
{
  const Int_t edgeXbins = 1;
  const Int_t edgeYbins = 1;
  const Int_t medianFiltXbins = 1;
  const Int_t medianFiltYbins = 1;
  const Double_t gaussFilterWidth = 0.10;

  TH2F *tempHist;
  if((faceNum == 0) || (faceNum == 1))
  {
    tempHist = new TH2F("tempHist","",nCalibDivisions_z+1,0,nCalibDivisions_z+1,nCalibDivisions_x+1,0,nCalibDivisions_x+1);
    for(Int_t x = 0; x <= nCalibDivisions_x; x++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        tempHist->SetBinContent(z+1,x+1,inputMap[x][z]);
      }
    }
  }
  else if((faceNum == 2) || (faceNum == 3))
  {
    tempHist = new TH2F("tempHist","",nCalibDivisions_x+1,0,nCalibDivisions_x+1,nCalibDivisions_y+1,0,nCalibDivisions_y+1);
    for(Int_t x = 0; x <= nCalibDivisions_x; x++)
    {
      for(Int_t y = 0; y <= nCalibDivisions_y; y++)
      {
        tempHist->SetBinContent(x+1,y+1,inputMap[x][y]);
      }
    }
  }
  else
  {
    tempHist = new TH2F("tempHist","",nCalibDivisions_z+1,0,nCalibDivisions_z+1,nCalibDivisions_y+1,0,nCalibDivisions_y+1);
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        tempHist->SetBinContent(z+1,y+1,inputMap[y][z]);
      }
    }
  }
  
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

  if((faceNum == 0) || (faceNum == 1))
  {
    for(Int_t x = 0; x <= nCalibDivisions_x; x++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        inputMap[x][z] = interpHist->GetBinContent(z+1,x+1);
      }
    }
  }
  else if((faceNum == 2) || (faceNum == 3))
  {
    for(Int_t x = 0; x <= nCalibDivisions_x; x++)
    {
      for(Int_t y = 0; y <= nCalibDivisions_y; y++)
      {
        inputMap[x][y] = interpHist->GetBinContent(x+1,y+1);
      }
    }
  }
  else
  {
    for(Int_t y = 0; y <= nCalibDivisions_y; y++)
    {
      for(Int_t z = 0; z <= nCalibDivisions_z; z++)
      {
        inputMap[y][z] = interpHist->GetBinContent(z+1,y+1);
      }
    }
  }
  
  return;
}

void extrapolate(TH2F *interpHist, Int_t edgeXbins, Int_t edgeYbins, Int_t faceNum)
{
  for(Int_t i = edgeXbins; i < interpHist->GetNbinsX()-edgeXbins; i++)
  {
    std::vector<Double_t> X, Y;
    tk::spline s;
    
    for(Int_t j = edgeYbins; j < interpHist->GetNbinsY()-edgeYbins; j++)
    {
      X.push_back(interpHist->GetYaxis()->GetBinCenter(j+1));
      Y.push_back(interpHist->GetBinContent(i+1,j+1));
    }

    if(((faceNum == 0) || (faceNum == 1)) && (edgeYbins != 0))
    {
      X.push_back(interpHist->GetYaxis()->GetBinCenter(interpHist->GetNbinsY()));
      Y.push_back(0.0);
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
    
    for(Int_t i = edgeXbins; i < interpHist->GetNbinsX()-edgeXbins; i++)
    {
      X.push_back(interpHist->GetXaxis()->GetBinCenter(i+1));
      Y.push_back(interpHist->GetBinContent(i+1,j+1));
    }

    if(((faceNum == 2) || (faceNum == 3)) && (edgeXbins != 0))
    {
      X.push_back(interpHist->GetXaxis()->GetBinCenter(interpHist->GetNbinsX()));
      Y.push_back(0.0);
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

      for(Int_t j = edgeYbins; j < interpHist->GetNbinsY()-edgeYbins; j++)
      {
        X1.push_back(interpHist->GetYaxis()->GetBinCenter(j+1));
        Y1.push_back(interpHist->GetBinContent(m+1,j+1));
      }

      if(((faceNum == 0) || (faceNum == 1)) && (edgeYbins != 0))
      {
        X1.push_back(interpHist->GetYaxis()->GetBinCenter(interpHist->GetNbinsY()));
        Y1.push_back(0.0);
      }

      s1.set_points(X1,Y1);

      for(Int_t i = edgeXbins; i < interpHist->GetNbinsX()-edgeXbins; i++)
      {
        X2.push_back(interpHist->GetXaxis()->GetBinCenter(i+1));
        Y2.push_back(interpHist->GetBinContent(i+1,n+1));	
      }

      if(((faceNum == 2) || (faceNum == 3)) && (edgeXbins != 0))
      {
        X2.push_back(interpHist->GetXaxis()->GetBinCenter(interpHist->GetNbinsX()));
        Y2.push_back(0.0);
      }

      s2.set_points(X2,Y2);

      for(Int_t j = edgeYbins; j < interpHist->GetNbinsY()-edgeYbins; j++)
      {
        X3.push_back(interpHist->GetYaxis()->GetBinCenter(j+1));
        Y3.push_back(interpHist->GetBinContent(interpHist->GetNbinsX()-m,j+1));
      }

      if(((faceNum == 0) || (faceNum == 1)) && (edgeYbins != 0))
      {
        X3.push_back(interpHist->GetYaxis()->GetBinCenter(interpHist->GetNbinsY()));
        Y3.push_back(0.0);
      }

      s3.set_points(X3,Y3);

      for(Int_t i = edgeXbins; i < interpHist->GetNbinsX()-edgeXbins; i++)
      {
        X4.push_back(interpHist->GetXaxis()->GetBinCenter(i+1));
        Y4.push_back(interpHist->GetBinContent(i+1,interpHist->GetNbinsY()-n));	
      }

      if(((faceNum == 2) || (faceNum == 3)) && (edgeXbins != 0))
      {
        X4.push_back(interpHist->GetXaxis()->GetBinCenter(interpHist->GetNbinsX()));
        Y4.push_back(0.0);
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

Double_t getFaceCorr(Double_t xVal, Double_t yVal, Double_t zVal, Int_t faceNum, Int_t comp)
{
  Double_t offset = 0.0;
  
  if (xVal < 0.00001) {
    xVal = 0.00001;
  }
  if (xVal > Lx-0.00001) {
    xVal = Lx-0.00001;
  }

  if (yVal < 0.00001) {
    yVal = 0.00001;
  }
  if (yVal > Ly-0.00001) {
    yVal = Ly-0.00001;
  }

  if (zVal < 0.00001) {
    zVal = 0.00001;
  }
  if (zVal > Lz-0.00001) {
    zVal = Lz-0.00001;
  }

  if ((faceNum == 0) && (comp == 1)) {
    offset = faceCalibHistTopDeltaX->Interpolate(zVal,xVal);
  }
  else if ((faceNum == 0) && (comp == 2)) {
    offset = faceCalibHistTopDeltaY->Interpolate(zVal,xVal);
  }
  else if ((faceNum == 0) && (comp == 3)) {
    offset = faceCalibHistTopDeltaZ->Interpolate(zVal,xVal);
  }
  if ((faceNum == 1) && (comp == 1)) {
    offset = faceCalibHistBottomDeltaX->Interpolate(zVal,xVal);
  }
  else if ((faceNum == 1) && (comp == 2)) {
    offset = faceCalibHistBottomDeltaY->Interpolate(zVal,xVal);
  }
  else if ((faceNum == 1) && (comp == 3)) {
    offset = faceCalibHistBottomDeltaZ->Interpolate(zVal,xVal);
  }
  if ((faceNum == 2) && (comp == 1)) {
    offset = faceCalibHistUpstreamDeltaX->Interpolate(xVal,yVal);
  }
  else if ((faceNum == 2) && (comp == 2)) {
    offset = faceCalibHistUpstreamDeltaY->Interpolate(xVal,yVal);
  }
  else if ((faceNum == 2) && (comp == 3)) {
    offset = faceCalibHistUpstreamDeltaZ->Interpolate(xVal,yVal);
  }
  if ((faceNum == 3) && (comp == 1)) {
    offset = faceCalibHistDownstreamDeltaX->Interpolate(xVal,yVal);
  }
  else if ((faceNum == 3) && (comp == 2)) {
    offset = faceCalibHistDownstreamDeltaY->Interpolate(xVal,yVal);
  }
  else if ((faceNum == 3) && (comp == 3)) {
    offset = faceCalibHistDownstreamDeltaZ->Interpolate(xVal,yVal);
  }
  if ((faceNum == 4) && (comp == 1)) {
    offset = faceCalibHistCathodeDeltaX->Interpolate(zVal,yVal);
  }
  else if ((faceNum == 4) && (comp == 2)) {
    offset = faceCalibHistCathodeDeltaY->Interpolate(zVal,yVal);
  }
  else if ((faceNum == 4) && (comp == 3)) {
    offset = faceCalibHistCathodeDeltaZ->Interpolate(zVal,yVal);
  }
  
  return offset;
}

Double_t getBulkCorr(Double_t xVal, Double_t yVal, Double_t zVal, Int_t comp)
{
  Double_t offset = 0.0;
  
  if (xVal < 0.00001) {
    xVal = 0.00001;
  }
  if (xVal > Lx-0.00001) {
    xVal = Lx-0.00001;
  }

  if (yVal < 0.00001) {
    yVal = 0.00001;
  }
  if (yVal > Ly-0.00001) {
    yVal = Ly-0.00001;
  }

  if (zVal < 0.00001) {
    zVal = 0.00001;
  }
  if (zVal > Lz-0.00001) {
    zVal = Lz-0.00001;
  }

  if (comp == 1) {
    offset = bulkCalibHistDeltaX->Interpolate(xVal,yVal,zVal);
  }
  else if (comp == 2) {
    offset = bulkCalibHistDeltaY->Interpolate(xVal,yVal,zVal);
  }
  else if (comp == 3) {
    offset = bulkCalibHistDeltaZ->Interpolate(xVal,yVal,zVal);
  }
  
  return offset;
}
