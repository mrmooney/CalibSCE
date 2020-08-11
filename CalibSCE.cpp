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

using namespace std;

const Char_t *inputFileName = "data/MC_Cosmics.root"; // NEED TO OBTAIN INPUT FILE FROM FNAL MACHINES FIRST
//const Char_t *inputFileName = "data/Data_Run1_EXTBNB.root"; // NEED TO OBTAIN INPUT FILE FROM FNAL MACHINES FIRST

const Char_t *simInterpFileName = "data/output_siminterp_MicroBooNE_4p5_gap.root";

const Bool_t isMC = true;
const Bool_t isSCEon = true;
const Int_t numCalibTracks = 200000;

Double_t VelRatio;
Double_t TrueAnode;
Double_t TrueCathode;
Double_t TrueTop;
Double_t TrueBottom;
Double_t TrueUpstream;
Double_t TrueDownstream;
Double_t ShiftedAnode;
Double_t ShiftedCathode;
Double_t OffsetCathodeReco;

TFile* outputFile = new TFile("output.root","RECREATE");

const Double_t Lx = 2.5;
const Double_t Ly = 2.5;
const Double_t Lz = 10.0;

const Double_t ScaleFactorX = Lx/2.56;
const Double_t ScaleFactorY = Ly/2.33;
const Double_t ScaleFactorZ = Lz/10.37;

const Double_t WF_top = 117.1;
const Double_t WF_bottom = -115.1;
const Double_t WF_upstream = 0.2;
const Double_t WF_downstream = 1036.6;
const Double_t WF_cathode = 254.4;

const Double_t relAngleCut = 20.0;
const Double_t maxXdist = 0.05;
const Double_t maxYdist = 0.20;
const Double_t maxZdist = 0.20;

Int_t minInputTrackNum = 0;
Int_t maxInputTrackNum = 1150000;

Int_t nCalibDivisions = 25;

const Double_t piVal = 3.14159265;

Int_t nCalibDivisions_x;
Int_t nCalibDivisions_y;
Int_t nCalibDivisions_z;

vector<Double_t> calibWeight[101][101][401];
vector<Double_t> calibDeltaX[101][101][401];
vector<Double_t> calibDeltaY[101][101][401];
vector<Double_t> calibDeltaZ[101][101][401];

TH3F *ResultDeltaX;
TH3F *ResultDeltaY;
TH3F *ResultDeltaZ;

TH3F *TrueFwdDeltaX;
TH3F *TrueFwdDeltaY;
TH3F *TrueFwdDeltaZ;

TH3F *TrueBkwdDeltaX;
TH3F *TrueBkwdDeltaY;
TH3F *TrueBkwdDeltaZ;

TH3F *RecoFwdDeltaX;
TH3F *RecoFwdDeltaY;
TH3F *RecoFwdDeltaZ;

TH3F *RecoBkwdDeltaX;
TH3F *RecoBkwdDeltaY;
TH3F *RecoBkwdDeltaZ;

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

Double_t doCoordTransformX(const Double_t inputX);
Double_t doCoordTransformY(const Double_t inputY);
Double_t doCoordTransformZ(const Double_t inputZ);
vector<Double_t> getParabolaParameters(const vector<elecInfo> &parabola_points_track);
vector<Double_t> findClosestPOA(const calibTrackInfo &trackA, const calibTrackInfo &trackB);
vector<Double_t> findDistortedClosestPOA(const calibTrackInfo &trackA, const calibTrackInfo &trackB);
void getLArSoftTrackSet(vector<trackInfo> &tracks, Int_t maxCosmicTracks, Double_t minTrackMCS_anode, Double_t minTrackMCS_cathode, Double_t minTrackMCS_crossing);
vector<calibTrackInfo> makeCalibTracks(const vector<trackInfo> &tracks);
void doCosmicCosmicCalib(const vector<calibTrackInfo> &laserCalibTracks, Double_t distScale, Double_t maxDistFactor, Int_t saveInfo);
void doCalibration(const vector<trackInfo> &cosmicTracks, Double_t distScale, Double_t maxDistFactor, Int_t numIterations, Int_t saveInfo);
void saveTrackInfo(const vector<trackInfo> &tracks);
void loadMaps();
Double_t getOffset(Double_t xVal, Double_t yVal, Double_t zVal, Int_t comp, Int_t calibMode);
  
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

  if(isMC == true) {
    VelRatio = 1.0; // MC
    TrueAnode = Lx; // MC
    TrueCathode = Lx*(2.56-2.548)/2.56; // MC
    TrueTop = Ly*(1.174+1.165)/(2.33); // MC
    TrueBottom = Ly*(-1.154+1.165)/(2.33); // MC
    TrueUpstream = Lz*(0.004)/(10.37); // MC
    TrueDownstream = Lz*(10.368)/(10.37); // MC
    ShiftedAnode = Lx*(2.56-0.0006)/2.56; // MC
    ShiftedCathode = Lx*(2.56-2.5524)/2.56; // MC
    OffsetCathodeReco = 0.004*Lx/2.56; // MC
  }
  else {
    VelRatio = 0.992; // DATA
    TrueAnode = Lx; // DATA
    TrueCathode = 0.0; // DATA
    TrueTop = Ly*(1.171+1.165)/(2.33); // DATA
    TrueBottom = Ly*(-1.151+1.165)/(2.33); // DATA
    TrueUpstream = Lz*(0.002)/(10.37); // DATA
    TrueDownstream = Lz*(10.366)/(10.37); // DATA
    ShiftedAnode = Lx*(2.56-(-0.0056))/2.56; // DATA
    ShiftedCathode = Lx*(2.56-2.5818*VelRatio)/2.56; // DATA
    OffsetCathodeReco = 0.004*Lx/2.56; // DATA
  }

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

  ResultDeltaX = new TH3F("ResultDeltaX","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  ResultDeltaY = new TH3F("ResultDeltaY","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  ResultDeltaZ = new TH3F("ResultDeltaZ","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));

  TrueFwdDeltaX = new TH3F("TrueFwdDeltaX","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  TrueFwdDeltaY = new TH3F("TrueFwdDeltaY","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  TrueFwdDeltaZ = new TH3F("TrueFwdDeltaZ","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));

  TrueBkwdDeltaX = new TH3F("TrueBkwdDeltaX","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  TrueBkwdDeltaY = new TH3F("TrueBkwdDeltaY","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  TrueBkwdDeltaZ = new TH3F("TrueBkwdDeltaZ","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));

  RecoFwdDeltaX = new TH3F("RecoFwdDeltaX","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  RecoFwdDeltaY = new TH3F("RecoFwdDeltaY","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  RecoFwdDeltaZ = new TH3F("RecoFwdDeltaZ","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));

  RecoBkwdDeltaX = new TH3F("RecoBkwdDeltaX","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  RecoBkwdDeltaY = new TH3F("RecoBkwdDeltaY","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));
  RecoBkwdDeltaZ = new TH3F("RecoBkwdDeltaZ","",nCalibDivisions_x+1,-Lx/(2.0*((Double_t) nCalibDivisions_x)),Lx+Lx/(2.0*((Double_t) nCalibDivisions_x)),nCalibDivisions_y+1,-Ly/(2.0*((Double_t) nCalibDivisions_y)),Ly+Ly/(2.0*((Double_t) nCalibDivisions_y)),nCalibDivisions_z+1,-Lz/(2.0*((Double_t) nCalibDivisions_z)),Lz+Lz/(2.0*((Double_t) nCalibDivisions_z)));

  loadMaps();
  
  outputFile->cd();

  //////////////////////////////////////////////
  /// MAIN PART OF CODE (CHANGE THESE THINGS)
  //////////////////////////////////////////////
  
  vector<trackInfo> cosmicTracks;
  
  getLArSoftTrackSet(cosmicTracks,numCalibTracks,minTrackMCS_anode,minTrackMCS_cathode,100000000.0);
  doCalibration(cosmicTracks,0.01,3,1,1);
  
  timer.Stop();
  cout << "Calibration Time:  " << timer.CpuTime() << " sec." << endl;
  
  outputFile->Write();
  outputFile->Close();
  
  return 0;
}

Double_t doCoordTransformX(const Double_t inputX)
{
  Double_t outputX;
  outputX = Lx - (Lx/2.56)*inputX*VelRatio/100.0;

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

void getLArSoftTrackSet(vector<trackInfo> &tracks, Int_t maxCosmicTracks, Double_t minTrackMCS_anode, Double_t minTrackMCS_cathode, Double_t minTrackMCS_crossing)
{
  tracks.clear();

  TFile* inputFile = new TFile(inputFileName,"READ");
  
  char* filestring = (char*) "";
  char* filestring2 = (char*) "";
  if (isSCEon == false) {
    filestring = (char*) "_copy";
    filestring2 = (char*) "t0ana/";
  }
  
  TTreeReader reader(Form("%sSCEtree%s",filestring2,filestring), inputFile);
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
    if (((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) {
      if (xS < xE) {
        cathodeOffset = getOffset(0.0,yS+getOffset(0.0,yS,zS,2,1),zS+getOffset(0.0,yS,zS,3,1),1,0);
      }
      else {
        cathodeOffset = getOffset(0.0,yE+getOffset(0.0,yE,zE,2,1),zE+getOffset(0.0,yE,zE,3,1),1,0);
      }
    
      xS += SCEfactor*(TrueCathode-ShiftedCathode) + SCEfactor*cathodeOffset;
      xE += SCEfactor*(TrueCathode-ShiftedCathode) + SCEfactor*cathodeOffset;
    }
    else if (((xS > (Lx - maxXdist)) && (xE > maxXdist)) || ((xE > (Lx - maxXdist)) && (xS > maxXdist))) {
      xS += (TrueAnode-ShiftedAnode);
      xE += (TrueAnode-ShiftedAnode);
    }

    x0 = xS + SCEfactor*getOffset(xS,yS,zS,1,1);
    y0 = yS + SCEfactor*getOffset(xS,yS,zS,2,1);
    z0 = zS + SCEfactor*getOffset(xS,yS,zS,3,1);

    x1 = xE + SCEfactor*getOffset(xE,yE,zE,1,1);
    y1 = yE + SCEfactor*getOffset(xE,yE,zE,2,1);
    z1 = zE + SCEfactor*getOffset(xE,yE,zE,3,1);
    
    Double_t trackLength = sqrt(pow(x0-x1,2.0)+pow(y0-y1,2.0)+pow(z0-z1,2.0));
    
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

    for(Int_t j = 0; j < *nPoints; j++)
    {
      if (((j % 10) != 0) && (j != *nPoints-1)) continue;
  	
      electron.x = -1.0; // Dummy (currently don't save this info in file)
      electron.y = -1.0; // Dummy (currently don't save this info in file)
      electron.z = -1.0; // Dummy (currently don't save this info in file)
      electron.t = -1.0; // Dummy (currently don't save this info in file)
      if (((xS < (Lx - maxXdist)) && (xE < maxXdist)) || ((xE < (Lx - maxXdist)) && (xS < maxXdist))) {
        electron.x_mod = doCoordTransformX(pointX[j])+SCEfactor*(TrueCathode-ShiftedCathode)+SCEfactor*cathodeOffset;
      }
      else if (((xS > (Lx - maxXdist)) && (xE > maxXdist)) || ((xE > (Lx - maxXdist)) && (xS > maxXdist))) {
        electron.x_mod = doCoordTransformX(pointX[j])+(TrueAnode-ShiftedAnode);
      }
      else {
        electron.x_mod = doCoordTransformX(pointX[j]);
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

void doCalibration(const vector<trackInfo> &cosmicTracks, Double_t distScale, Double_t maxDistFactor, Int_t numIterations, Int_t saveInfo)
{
  ResultDeltaX->Reset();
  ResultDeltaY->Reset();
  ResultDeltaZ->Reset();

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
  
  vector<calibTrackInfo> cosmicCalibTracks = makeCalibTracks(cosmicTracks);

  Double_t x_true, y_true, z_true;
  Double_t x_reco, y_reco, z_reco;
  Double_t Dx, Dy, Dz;
  Int_t elecFate;
  Int_t numEntries;

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
  T_calib->Branch("numEntries",&numEntries,"data_calib/I");
  if(saveInfo == 1) {
    T_calib->SetDirectory(outputFile);
  }

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

	numEntries = calibWeight[x][y][z].size();
	Double_t sumWeights = 0.0;
	
	for(Int_t k = 0; k < numEntries; k++)
	{
	  sumWeights += (calibWeight[x][y][z])[k];
	  
          if((calibWeight[x][y][z])[k] < 0.0)
	  {
	    (calibWeight[x][y][z])[k] = 0.0;
	  }
	}
	
        if((numEntries > 0) && (sumWeights > 0.0))
	{
          elecFate = 1;

          Dx = TMath::Median(numEntries,&(calibDeltaX[x][y][z])[0],&(calibWeight[x][y][z])[0]);
          Dy = TMath::Median(numEntries,&(calibDeltaY[x][y][z])[0],&(calibWeight[x][y][z])[0]);
          Dz = TMath::Median(numEntries,&(calibDeltaZ[x][y][z])[0],&(calibWeight[x][y][z])[0]);
         
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
    
        ResultDeltaX->SetBinContent(x+1,y+1,z+1,Dx);
        ResultDeltaY->SetBinContent(x+1,y+1,z+1,Dy);
        ResultDeltaZ->SetBinContent(x+1,y+1,z+1,Dz);

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

void loadMaps()
{
  TFile* simInterpFile = new TFile(simInterpFileName,"READ");

  TH3F* truefwd_Dx = (TH3F*) simInterpFile->Get("TrueFwd_Displacement_X");
  TH3F* truefwd_Dy = (TH3F*) simInterpFile->Get("TrueFwd_Displacement_Y");
  TH3F* truefwd_Dz = (TH3F*) simInterpFile->Get("TrueFwd_Displacement_Z");

  TH3F* truebkwd_Dx = (TH3F*) simInterpFile->Get("TrueBkwd_Displacement_X");
  TH3F* truebkwd_Dy = (TH3F*) simInterpFile->Get("TrueBkwd_Displacement_Y");
  TH3F* truebkwd_Dz = (TH3F*) simInterpFile->Get("TrueBkwd_Displacement_Z");

  TH3F* recofwd_Dx = (TH3F*) simInterpFile->Get("RecoFwd_Displacement_X");
  TH3F* recofwd_Dy = (TH3F*) simInterpFile->Get("RecoFwd_Displacement_Y");
  TH3F* recofwd_Dz = (TH3F*) simInterpFile->Get("RecoFwd_Displacement_Z");

  TH3F* recobkwd_Dx = (TH3F*) simInterpFile->Get("RecoBkwd_Displacement_X");
  TH3F* recobkwd_Dy = (TH3F*) simInterpFile->Get("RecoBkwd_Displacement_Y");
  TH3F* recobkwd_Dz = (TH3F*) simInterpFile->Get("RecoBkwd_Displacement_Z");


  for(Int_t x = 1; x <= nCalibDivisions_x+1; x++)
  {
    for(Int_t y = 1; y <= nCalibDivisions_y+1; y++)
    {
      for(Int_t z = 1; z <= nCalibDivisions_z+1; z++)
      {        
        TrueFwdDeltaX->SetBinContent(x,y,z,-0.01*truefwd_Dx->GetBinContent(nCalibDivisions_x+2-x,y,z));
        TrueFwdDeltaY->SetBinContent(x,y,z,0.01*truefwd_Dy->GetBinContent(nCalibDivisions_x+2-x,y,z));
        TrueFwdDeltaZ->SetBinContent(x,y,z,0.01*truefwd_Dz->GetBinContent(nCalibDivisions_x+2-x,y,z));

        TrueBkwdDeltaX->SetBinContent(x,y,z,-0.01*truebkwd_Dx->GetBinContent(nCalibDivisions_x+2-x,y,z));
        TrueBkwdDeltaY->SetBinContent(x,y,z,0.01*truebkwd_Dy->GetBinContent(nCalibDivisions_x+2-x,y,z));
        TrueBkwdDeltaZ->SetBinContent(x,y,z,0.01*truebkwd_Dz->GetBinContent(nCalibDivisions_x+2-x,y,z));

        RecoFwdDeltaX->SetBinContent(x,y,z,-0.01*recofwd_Dx->GetBinContent(nCalibDivisions_x+2-x,y,z));
        RecoFwdDeltaY->SetBinContent(x,y,z,0.01*recofwd_Dy->GetBinContent(nCalibDivisions_x+2-x,y,z));
        RecoFwdDeltaZ->SetBinContent(x,y,z,0.01*recofwd_Dz->GetBinContent(nCalibDivisions_x+2-x,y,z));

        RecoBkwdDeltaX->SetBinContent(x,y,z,-0.01*recobkwd_Dx->GetBinContent(nCalibDivisions_x+2-x,y,z));
        RecoBkwdDeltaY->SetBinContent(x,y,z,0.01*recobkwd_Dy->GetBinContent(nCalibDivisions_x+2-x,y,z));
        RecoBkwdDeltaZ->SetBinContent(x,y,z,0.01*recobkwd_Dz->GetBinContent(nCalibDivisions_x+2-x,y,z));
      }
    }
  }

  return;
}

Double_t getOffset(Double_t xVal, Double_t yVal, Double_t zVal, Int_t comp, Int_t calibMode)
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

  if ((isMC == true) && (isSCEon == true)) {
    if (calibMode == 0) {
      if (comp == 1) {
        offset = TrueFwdDeltaX->Interpolate(xVal,yVal,zVal);
      }
      else if (comp == 2) {
        offset = TrueFwdDeltaY->Interpolate(xVal,yVal,zVal);
      }
      else if (comp == 3) {
        offset = TrueFwdDeltaZ->Interpolate(xVal,yVal,zVal);
      }
    }
    else if (calibMode == 1) {
      if (comp == 1) {
        offset = TrueBkwdDeltaX->Interpolate(xVal,yVal,zVal);
      }
      else if (comp == 2) {
        offset = TrueBkwdDeltaY->Interpolate(xVal,yVal,zVal);
      }
      else if (comp == 3) {
        offset = TrueBkwdDeltaZ->Interpolate(xVal,yVal,zVal);
      }
    }
  }
  else if (isMC == false) {
    if (calibMode == 0) {
      if (comp == 1) {
        offset = RecoFwdDeltaX->Interpolate(xVal,yVal,zVal);
      }
      else if (comp == 2) {
        offset = RecoFwdDeltaY->Interpolate(xVal,yVal,zVal);
      }
      else if (comp == 3) {
        offset = RecoFwdDeltaZ->Interpolate(xVal,yVal,zVal);
      }
    }
    else if (calibMode == 1) {
      if (comp == 1) {
        offset = RecoBkwdDeltaX->Interpolate(xVal,yVal,zVal);
      }
      else if (comp == 2) {
        offset = RecoBkwdDeltaY->Interpolate(xVal,yVal,zVal);
      }
      else if (comp == 3) {
        offset = RecoBkwdDeltaZ->Interpolate(xVal,yVal,zVal);
      }
    }
  }
  
  return offset;
}
