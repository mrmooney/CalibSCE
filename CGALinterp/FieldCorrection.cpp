// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <cmath>
#include <cstring>
#include <thread>
#include <array>

// C headers
#include <pthread.h>
#include <unistd.h>

// ROOT headers
#include "TROOT.h"
#include "TObject.h"
#include "TApplication.h"
#include "TAttLine.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TProfile.h"
#include "TPolyLine3D.h"
#include "TStyle.h"
#include "TFrame.h"
#include "TFile.h"
#include "TVirtualPad.h"
#include "TView.h"
#include "TView3D.h"
#include "TTree.h"

// #include "TTreeReader.h"
// #include "TTreeReaderValue.h"

// Own Files
//#include "Geometry.hpp" // COMMENTED OUT BY MRM
#include "LaserTrack.hpp"
#include "ThreeVector.hpp"
#include "TPCVolumeHandler.hpp"
#include "Interpolation3D.hpp"
#include "Matrix3x3.hpp"
#include "Laser.hpp"

int main ();
void LaserInterpThread(Laser&, const Laser&, const Delaunay&);
void WriteRootFile(std::vector<ThreeVector<float>>&, TPCVolumeHandler&);
void WriteTextFile(std::vector<ThreeVector<float>>&);
std::vector<Laser> ReadMooneyTracks(const std::string&, TPCVolumeHandler&);
std::vector<Laser> ReadMooneyGrid(const std::string&, TPCVolumeHandler&);
void WriteMooneyGrid(const std::string&, Laser, LaserTrack);
// void DrawSpaceCharge();
// void DrawField();

const double kPi = 3.14159265358979323846;
const float kDetector[] = {2.5604,2.325,10.368}; // Size of the detectors in m (used from docdb 1821)
const int kResolution[] = {26,26,101}; // Field resolution in every coordinate
const int kTrackRes = 100; // Number of samples in every laser track
const int kSqrtNumberOfTracks = 101;
const float kLaserOffset[3] = {1.5,1.1625,-0.21}; // Laser head offset from field cage in [m]

int main ()
{
  
//   ThreeVector<float> DetectorOffset = {0.0,0.0,0.0};
  
//   ThreeVector<float> DetectorSize = {2.5604,2.325,10.368};
//   ThreeVector<float> DetectorOffset = {0.0,-DetectorSize[1]/(float)2.0,0.0};
  
  ThreeVector<float> DetectorSize = {2.5,2.5,10.0};
  ThreeVector<float> DetectorOffset = {0.0,0.0,0.0};
  ThreeVector<int> DetectorResolution = {26,26,101};
   
  TPCVolumeHandler Detector(DetectorSize, DetectorOffset, DetectorResolution);
  
  double HistRange[3][2];
  for (int coord = 0; coord < 3; coord++) for(int minmax = 0; minmax < 2; minmax++)
  {
    HistRange[coord][minmax] = (DetectorSize[coord]+pow(2,-25))*( minmax - pow(-1,minmax)/((double)DetectorResolution[coord]-1.0)/2.0 );
  }
  
//   auto start = std::chrono::high_resolution_clock::now();
//   auto elapsed = std::chrono::high_resolution_clock::now() - start;
//   std::cout << "Function Time [ns]: " << std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count() << std::endl;
    
//   time_t timer;
//   std::time(&timer);
  // Enter Functions here:
  unsigned BeamBins = 101;
  std::array<float,2> LaserAngles_1;
  std::array<float,2> LaserAngles_2;
  ThreeVector<float> LaserPosition_1 = {1.3,0.0,-0.21};
  ThreeVector<float> LaserPosition_2 = {1.3,0.0,Detector.GetDetectorSize().at(2)+ (float)0.21};
  
  
  std::vector<Laser> LaserTrackSets = ReadMooneyGrid("output.root", Detector);
  LaserTrack InterpolatedSet(LaserAngles_1,LaserPosition_1,Detector);
  
  Delaunay Mesh = TrackMesher(LaserTrackSets.at(0).GetTrackSet());
  
  ThreeVector<float> InterpCorr;
  ThreeVector<float> Pos;
  for(unsigned sample = 0; sample < LaserTrackSets.at(1).GetFirstTrack().GetNumberOfSamples(); sample++)
  {
    InterpCorr = InterpolateCGAL(LaserTrackSets.at(0).GetTrackSet(), Mesh, LaserTrackSets.at(1).GetFirstTrack().GetSamplePosition(sample));
    Pos = LaserTrackSets.at(1).GetFirstTrack().GetSamplePosition(sample);
    InterpolatedSet.AppendSample(Pos, InterpCorr);
  }
  
  WriteMooneyGrid("new_output.root", LaserTrackSets.at(0), InterpolatedSet);
  
  
//   std::vector<Laser> LaserTrackSets = ReadMooneyTracks("laserDataSCE.root",Detector);
  
//   std::vector<Laser> LaserTrackSets;
//   LaserTrackSets.resize(2);
  
//   std::vector<LaserTrack> TrackVector;
//   for(unsigned theta_entry = 0; theta_entry < BeamBins; theta_entry++)
//   {
//     LaserAngles_1[0] = -45.0 + 90.0/(float)(BeamBins-1)*theta_entry;
//     LaserAngles_2[0] = -45.0 + 90.0/(float)(BeamBins-1)*theta_entry;
//     for(unsigned phi_entry = 0; phi_entry < BeamBins; phi_entry++)
//     {
//       LaserAngles_1[1] = -45 + 90.0/(float)(BeamBins-1)*phi_entry;
//       LaserAngles_2[1] = -225 + 90.0/(float)(BeamBins-1)*phi_entry;
      
//       LaserTrackSets.at(0).AppendTrack(LaserTrack(50, LaserAngles_1, LaserPosition_1, Detector));
//       LaserTrackSets.at(1).AppendTrack(LaserTrack(50, LaserAngles_2, LaserPosition_2, Detector));
//     }
//   }
  
//   for(auto& laser : LaserTrackSets)
//   {
//     laser.DistortTrackSet("Field.root",Detector);
//   }
  
//   std::vector<Delaunay> MeshVector;
//   
//   for(unsigned set_no = 0; set_no < LaserTrackSets.size(); set_no++)
//   {
//     LaserTrackSets.at(set_no).CorrectTrackSet();
//     MeshVector.push_back( TrackMesher(LaserTrackSets.at(set_no).GetTrackSet()) );
//   }
  
//   LaserTrackSets.front().DrawTrack(1000);
  
//   Laser TempLaser_1 = LaserTrackSets.front();
//   Laser TempLaser_2 = LaserTrackSets.back();
  
//   std::cout << "Start 1st track interpolation" << std::endl;
//   LaserInterpThread(LaserTrackSets.front(),TempLaser_2,MeshVector.back());
//   std::thread Thread_0( LaserInterpThread, std::ref(LaserTrackSets.front()), std::ref(TempLaser_2), std::ref(MeshVector.back()) );
  
//   std::cout << "Start 2nd track interpolation" << std::endl;
//   LaserInterpThread(LaserTrackSets.back(),TempLaser_1,MeshVector.front());
//   std::thread Thread_1 ( LaserInterpThread, std::ref(LaserTrackSets.back()), std::ref(TempLaser_1), std::ref(MeshVector.front()) );
  
//   Thread_0.join();
//   Thread_1.join();
  
//   for(unsigned track_no = 0; track_no < LaserTrackSets.front().GetNumberOfTracks(); track_no++)
//   {
//     InterpolateTrack(LaserTrackSets.front().GetTrack(track_no),LaserTrackSets.back().GetTrackSet(),MeshVector.back());
//   }
  
  
  
//   for(unsigned track_no = 0; track_no < LaserTrackSets.back().GetNumberOfTracks(); track_no++)
//   {
//     InterpolateTrack(LaserTrackSets.back().GetTrack(track_no),TempLaser.GetTrackSet(),MeshVector.front());
//   }
  
//   std::cout << "Cleanup!" << std::endl;
  
//   MeshVector.clear();
  
//   std::vector<LaserTrack> TrackVector = LaserTrackSets.back().GetTrackSet();
//   LaserTrackSets.pop_back();
//   TrackVector.insert(TrackVector.begin(),LaserTrackSets.front().begin(), LaserTrackSets.front().end());
//   LaserTrackSets.clear();
  
  
//   std::cout << "Start field map interpolation" << std::endl;
  
//   Delaunay Mesh = TrackMesher(TrackVector);

  
//   for(unsigned iter = 0; iter < TrackVector.size(); iter++)
//   {
//     for(unsigned coord = 0; coord < 3; coord++)
//       std::cout << TrackVector[iter].LaserCorrection[20][coord] << " " << TrackVector[iter].LaserCorrection[20][coord] << " ";
//     std::cout << std::endl;
//   }

//   std::vector<ThreeVector<float>> Displacement;
  
//   std::cout << "Create Mesh..." << std::endl;
//   Delaunay Mesh_1 = TrackMesher(TrackVector);
//   std::cout << "Mesh done!" << std::endl;
  
//   Point bla(1.0,1.0,3.0);
//   ThreeVector<float> blabla = {1.0,1.0,3.0};
//   std::vector<ThreeVector<float>> shit;
//   shit.push_back(blabla);
//   Delaunay::Cell_handle Cell =  Mesh.locate(bla);
  
//   for(unsigned vertex_no = 0; vertex_no < 4; vertex_no++)
//   {
//     std::cout << Cell->vertex(vertex_no)->point()[0] << " " << Cell->vertex(vertex_no)->point()[1] << " " << Cell->vertex(vertex_no)->point()[2] << std::endl;
//   }
  
//   InterpolateCGAL(TrackVector,Mesh,Location);
  
//   ThreeVector<float> Location = {1.2,1.0,5.1};
  
//   ThreeVector<float> Location;
//   for(unsigned xbin = 0; xbin < DetectorResolution[0]; xbin++) 
//   {
//     std::cout << xbin << std::endl;
//     Location[0] = DetectorOffset[0]+(HistRange[0][1]-HistRange[0][0])/(float)DetectorResolution[0]*xbin;
//     for(unsigned ybin = 0; ybin < DetectorResolution[1]; ybin++) 
//     {
//       Location[1] = DetectorOffset[1]+(HistRange[1][1]-HistRange[1][0])/(float)DetectorResolution[1]*ybin;
//       for(unsigned zbin = 0; zbin < DetectorResolution[2]; zbin++)
//       {
// 	Location[2] = DetectorOffset[2]+(HistRange[2][1]-HistRange[2][0])/(float)DetectorResolution[2]*zbin;
// 	Displacement.push_back(InterpolateCGAL(TrackVector,Mesh,Location));
//       }
//     }
//   }
  
//   WriteRootFile(Displacement,Detector);
//   WriteTextFile(Displacement);
  
//   Displacement = InterpolateCGAL(TrackVector,Mesh,Location);
//   ThreeVector<float> Location = {1.2,1.0,5.1};
//   Displacement = InterpolateCGAL(TrackVector,Mesh,Location);
//   GetClosestTracksInfo(TrackVector,4);
  
//   std::cout << "End of program after "<< std::difftime(std::time(NULL),timer) << " s" << std::endl;
//   App->Run();
}

void LaserInterpThread(Laser& LaserTrackSet, const Laser& InterpolationLaser, const Delaunay& InterpolationMesh)
{
  LaserTrackSet.InterpolateTrackSet(InterpolationLaser, InterpolationMesh);
}

void WriteRootFile(std::vector<ThreeVector<float>>& InterpolationData, TPCVolumeHandler& TPCVolume)
{ 
  double HistRange[3][2];
  for (int coord = 0; coord < 3; coord++) for(int minmax = 0; minmax < 2; minmax++)
  {
    HistRange[coord][minmax] = (TPCVolume.GetDetectorSize()[coord]+pow(2,-25))*( minmax - pow(-1,minmax)/((double)TPCVolume.GetDetectorResolution()[coord]-1.0)/2.0 );
  }
  
  ThreeVector<float> DetectorOffset = {0.0,-TPCVolume.GetDetectorSize()[1]/(float)2.0,0.0};
  
  std::vector<TH3F*> RecoField;
  RecoField.push_back(new TH3F("Reco_Field_X","Reco Field X",TPCVolume.GetDetectorResolution()[0],HistRange[0][0],HistRange[0][1],TPCVolume.GetDetectorResolution()[1],HistRange[1][0],HistRange[1][1],TPCVolume.GetDetectorResolution()[2],HistRange[2][0],HistRange[2][1]));
  RecoField.push_back(new TH3F("Reco_Field_Y","Reco Field Y",TPCVolume.GetDetectorResolution()[0],HistRange[0][0],HistRange[0][1],TPCVolume.GetDetectorResolution()[1],HistRange[1][0],HistRange[1][1],TPCVolume.GetDetectorResolution()[2],HistRange[2][0],HistRange[2][1]));
  RecoField.push_back(new TH3F("Reco_Field_Z","Reco Field Z",TPCVolume.GetDetectorResolution()[0],HistRange[0][0],HistRange[0][1],TPCVolume.GetDetectorResolution()[1],HistRange[1][0],HistRange[1][1],TPCVolume.GetDetectorResolution()[2],HistRange[2][0],HistRange[2][1]));
  
  
  for(unsigned xbin = 0; xbin < TPCVolume.GetDetectorResolution()[0]; xbin++) for(unsigned ybin = 0; ybin < TPCVolume.GetDetectorResolution()[1]; ybin++) for(unsigned zbin = 0; zbin < TPCVolume.GetDetectorResolution()[2]; zbin++)
  {
    for(unsigned coord = 0; coord < 3; coord++)
    {
      RecoField[coord]->SetBinContent(xbin+1,ybin+1,zbin+1, InterpolationData[zbin+( TPCVolume.GetDetectorResolution()[2]*(ybin+TPCVolume.GetDetectorResolution()[1]*xbin) )][coord]);
      std::cout << InterpolationData[zbin+( TPCVolume.GetDetectorResolution()[2]*(ybin+TPCVolume.GetDetectorResolution()[1]*xbin) )][coord] << " ";
    }
    std::cout << std::endl;
  }
  
  TFile *OutputFile = new TFile("RecoField.root", "recreate");
  for(unsigned coord = 0; coord < RecoField.size(); coord++)
  {
    RecoField[coord]->Write();
  }
  OutputFile -> Close();
  delete OutputFile;
//   gDirectory->GetList()->Delete();
}

void WriteTextFile(std::vector<ThreeVector<float>>& InterpolationData)
{
  std::ofstream OutputFile;
  OutputFile.open("Reco.txt", std::ios::out);
  for(unsigned entry = 0; entry < InterpolationData.size(); entry++)
  {
    OutputFile << InterpolationData[entry][0] << InterpolationData[entry][1] << InterpolationData[entry][2];
  }
  OutputFile.close();
}

void WriteMooneyGrid(const std::string& OutputFileName, Laser ExistingSet, LaserTrack InterpolatedSet)
{
  double x_true;
  double y_true;
  double z_true;
  
  double x_reco;
  double y_reco;
  double z_reco;
  
  double Dx;
  double Dy;
  double Dz;
  
  int elecFate;
  
  TFile* OutputFile = new TFile(OutputFileName.c_str(),"RECREATE");
  
  TTree *T_laser = new TTree("SpaCEtree_calib","SpaCEtree_calib");
  T_laser->Branch("x_true",&x_true,"x_true/D");
  T_laser->Branch("y_true",&y_true,"y_true/D");
  T_laser->Branch("z_true",&z_true,"z_true/D");
  T_laser->Branch("x_reco",&x_reco,"x_reco/D");
  T_laser->Branch("y_reco",&y_reco,"y_reco/D");
  T_laser->Branch("z_reco",&z_reco,"z_reco/D");
  T_laser->Branch("Dx",&Dx,"Dx/D");
  T_laser->Branch("Dy",&Dy,"Dy/D");
  T_laser->Branch("Dz",&Dz,"Dz/D");
  T_laser->Branch("elecFate",&elecFate,"elecFate/I");
  
  unsigned long int TotNumberOfSamples = InterpolatedSet.GetNumberOfSamples() + ExistingSet.GetFirstTrack().GetNumberOfSamples();
  
  for(unsigned long int sample = 0; sample < ExistingSet.GetFirstTrack().GetNumberOfSamples(); sample++)
  {    
    x_reco = ExistingSet.GetFirstTrack().GetSamplePosition(sample).at(0);
    y_reco = ExistingSet.GetFirstTrack().GetSamplePosition(sample).at(1);
    z_reco = ExistingSet.GetFirstTrack().GetSamplePosition(sample).at(2);
    
    Dx = ExistingSet.GetFirstTrack().GetCorrection(sample).at(0);
    Dy = ExistingSet.GetFirstTrack().GetCorrection(sample).at(1);
    Dz = ExistingSet.GetFirstTrack().GetCorrection(sample).at(2);
    
    x_true = Dx + x_reco;
    y_true = Dy + y_reco;
    z_true = Dz + z_reco;
    
    elecFate = 1;
    if(Dx == 0.0 && Dy == 0.0 && Dz == 0.0)
    {
      elecFate = 0;
    }
    
    T_laser->Fill();
  }
  
  for(unsigned long int sample = 0; sample < InterpolatedSet.GetNumberOfSamples(); sample++)
  {
    elecFate = 1;
    x_reco = InterpolatedSet.GetSamplePosition(sample).at(0);
    y_reco = InterpolatedSet.GetSamplePosition(sample).at(1);
    z_reco = InterpolatedSet.GetSamplePosition(sample).at(2);
    
    Dx = InterpolatedSet.GetCorrection(sample).at(0);
    Dy = InterpolatedSet.GetCorrection(sample).at(1);
    Dz = InterpolatedSet.GetCorrection(sample).at(2);
    
    x_true = Dx + x_reco;
    y_true = Dy + y_reco;
    z_true = Dz + z_reco;
    
    T_laser->Fill();
  }
  OutputFile->Write();
  delete T_laser;
  OutputFile->Close();
  delete OutputFile;
  gDirectory->GetList()->Delete();
  
}

std::vector<Laser> ReadMooneyGrid(const std::string& InputFileName, TPCVolumeHandler& DetectorVolume)
{
  double x_reco;
  double y_reco;
  double z_reco;
  
  double Dx;
  double Dy;
  double Dz;
  
  int elecFate;
  
  // This is only a hoax but needed to construct LaserTrack
  std::array<float,2> Angles = {0.0, 0.0};
  ThreeVector<float> LaserPosition = {0.0, 0.0, 0.0};
  
  
  TFile* InFile = new TFile (InputFileName.c_str());
  if(InFile->IsZombie())
  {
    std::cerr << "ERROR: Input file " << InputFileName << " invalid!" << std::endl;
    exit(-1);
  }
    
  TTree* T_laser = new TTree();
  T_laser = (TTree*)InFile->Get("SpaCEtree_calib");
  
  T_laser->SetBranchAddress("x_reco",&x_reco);
  T_laser->SetBranchAddress("y_reco",&y_reco);
  T_laser->SetBranchAddress("z_reco",&z_reco);
  T_laser->SetBranchAddress("Dx",&Dx);
  T_laser->SetBranchAddress("Dy",&Dy);
  T_laser->SetBranchAddress("Dz",&Dz);
  T_laser->SetBranchAddress("elecFate",&elecFate);
  
  
  // This is just a hack to get the wanted structure (sorry)
  std::vector<std::vector<LaserTrack>> DataSet;
  DataSet.resize(2);
  
  DataSet.at(0).push_back(LaserTrack(Angles,LaserPosition,DetectorVolume));
  DataSet.at(1).push_back(LaserTrack(Angles,LaserPosition,DetectorVolume));
  
  for(unsigned iter = 0; iter < T_laser->GetEntries(); iter++)
  {
    T_laser->GetEntry(iter);
    
    if(elecFate == 1)
      DataSet.at(0).back().AppendSample(x_reco,y_reco,z_reco,Dx,Dy,Dz);
    else
      DataSet.at(1).back().AppendSample(x_reco,y_reco,z_reco,Dx,Dy,Dz);
  }
  delete T_laser;
  InFile->Close();
  delete InFile;
  gDirectory->GetList()->Delete();
  
  std::vector<Laser> ReturnVector;
  ReturnVector.push_back(DataSet.at(0));
  ReturnVector.push_back(DataSet.at(1));
  return ReturnVector;
}

std::vector<Laser> ReadMooneyTracks(const std::string& InputFileName, TPCVolumeHandler& DetectorVolume)
{
  double x0;
  double y0;
  double z0;
  double theta;
  double phi;

   
  int nElec;
  double elecX[1000];
  double elecY[1000];
  double elecZ[1000];

   
  float xNum;
  float yNum;
  float zNum;
  float thetaNum;
  float phiNum;

   
  TFile* InFile = new TFile (InputFileName.c_str());
  if(InFile->IsZombie())
  {
    std::cerr << "ERROR: Input file " << InputFileName << " invalid!" << std::endl;
    exit(-1);
  }
    
  TTree* T_laser = new TTree();
  T_laser = (TTree*)InFile->Get("SpaCEtree_laser");
  
  std::vector<std::vector<LaserTrack>> LaserTrackSets;  
  
//   T_laser->SetDirectory(InFile);
  
//   gDirectory->cd();
  
//   TTree* T_laser = (TTree*)InFile->Get("MCTree");
  
//   TTreeReader myReader("SpaCEtree_laser", InFile);
//   TTreeReaderValue<int> NumberOfElectrons(myReader, "nElec");
  
//   while (myReader.Next()) 
//   {
//     std::cout << *NumberOfElectrons << std::endl;
//   }
  
  
//   TTree* T_laser = (TTree*)InFile->Get("SpaCEtree_laser");
//   T_laser -> GetEntries();
  
//   if(T_laser->IsZombie())
//   {
//     std::cerr << "ERROR: Unable to get TTree!" << std::endl;
//     exit(-1);
//   }
  std::cout << InFile << " " << T_laser << std::endl;
  
  std::cout << "alive" << std::endl;
      
//   T_laser->Print();  
  
      
  T_laser->SetBranchAddress("x0_laser",&x0);
  T_laser->SetBranchAddress("y0_laser",&y0);
  T_laser->SetBranchAddress("z0_laser",&z0);
  T_laser->SetBranchAddress("theta_laser",&theta);
  T_laser->SetBranchAddress("phi_laser",&phi);
  T_laser->SetBranchAddress("nElec_laser",&nElec);
  T_laser->SetBranchAddress("elecX_laser",&elecX);
  T_laser->SetBranchAddress("elecY_laser",&elecY);
  T_laser->SetBranchAddress("elecZ_laser",&elecZ);
      
  std::array<float,2> Angles;
  ThreeVector<float> LaserPosition;
  ThreeVector<float> LaserPosition_check = {(float)0xDEADBEEF,(float)0xDEADBEEF,(float)0xDEADBEEF};
  
  ThreeVector<float> LaserSample;
//   LaserTrack Track;
  
  for(unsigned iter = 0; iter < T_laser->GetEntries(); iter++)
  {
    T_laser->GetEntry(iter);
    Angles = {-(float)theta + (float)M_PI/2, (float)phi};
    LaserPosition = {(float)(2.5-x0), (float)y0, (float)z0};
    if(LaserPosition != LaserPosition_check)
    {
      LaserPosition_check = LaserPosition;
      LaserTrackSets.push_back(std::vector<LaserTrack>());
    }
    if(nElec > 1)
    {
      LaserTrackSets.back().push_back( LaserTrack(Angles,LaserPosition,DetectorVolume) );
      for(unsigned sample_no = 0; sample_no < nElec; sample_no++)
      {
	if(sample_no == 0) std::cout << elecX[sample_no] << " " << elecY[sample_no] << " " << elecZ[sample_no] << " " << std::endl;
//       LaserSample = {(float)elecX[sample_no],(float)elecY[sample_no],(float)elecZ[sample_no]};
//       LaserTracks_1.back().AppendSample(LaserSample);
	LaserTrackSets.back().back().AppendSample((float)(2.5 - elecX[sample_no]),(float)elecY[sample_no],(float)elecZ[sample_no]);
      }
    }
    else std::cout << iter << " " << nElec << std::endl;
  }
  delete T_laser;
  InFile->Close();
  delete InFile;
  gDirectory->GetList()->Delete();
  
  std::vector<Laser> ReturnVector; 
  ReturnVector.push_back(LaserTrackSets.front());
  ReturnVector.push_back(LaserTrackSets.back());
  return ReturnVector;
}
// void DrawSpaceCharge()
// {
//   TH2F * hSpaceCharge = new TH2F("Space Charge","Space Charge",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]);
//   for (int ybin = 0; ybin < kResolution[1]; ybin++) 
//   {
//     for (int xbin = 0; xbin < kResolution[0]; xbin++) for (int zbin = 0; zbin < kResolution[2]; zbin++)
//     {
//       hSpaceCharge -> SetBinContent(zbin+1,xbin+1,-fChargeDistribution[xbin][ybin][zbin]);
//     }
//     hSpaceCharge -> SetStats(0);
//     hSpaceCharge -> SetMaximum(1e-8);
//     hSpaceCharge -> SetMinimum(0);
//     TCanvas * C0 = new TCanvas("Space Charge","Space Charge",1000,500);
//     hSpaceCharge -> Draw("colz");
//     C0 -> Print("SpaceCharge.gif+5","gif+5");
//   }
//   delete hSpaceCharge;
// }
// 
// void DrawField()
// {
//   TH2F * hFieldX = new TH2F("Field Map X","Field Map X",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]);
//   TH2F * hFieldY = new TH2F("Field Map Y","Field Map Y",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]);
//   TH2F * hFieldZ = new TH2F("Field Map Z","Field Map Z",kResolution[2],0,kDetector[2],kResolution[0],0,kDetector[0]);
//   for (int ybin = 0; ybin < kResolution[1]; ybin++) 
//   {
//     for (int xbin = 0; xbin < kResolution[0]; xbin++) for (int zbin = 0; zbin < kResolution[2]; zbin++)
//     {
//       hFieldX -> SetBinContent(zbin+1,xbin+1,fRecoField[0][xbin][ybin][zbin]);
//       hFieldY -> SetBinContent(zbin+1,xbin+1,fRecoField[1][xbin][ybin][zbin]);
//       hFieldZ -> SetBinContent(zbin+1,xbin+1,fRecoField[2][xbin][ybin][zbin]);
//     }
//     hFieldX -> SetStats(0);
//     hFieldY -> SetStats(0);
//     hFieldZ -> SetStats(0);
//     
//     hFieldX -> SetMaximum(0.02);
//     hFieldX -> SetMinimum(-0.02);
//     hFieldY -> SetMaximum(0.07);
//     hFieldY -> SetMinimum(-0.07);
//     hFieldZ-> SetMaximum(0.05);
//     hFieldZ -> SetMinimum(-0.05);
//     
//     TCanvas * C1 = new TCanvas("Field Map X","Field Map X",1000,500);
//     hFieldX -> Draw("colz");
//     C1 -> Print("FieldX.gif+5","gif+5");
//     TCanvas * C2 = new TCanvas("Field Map Y","Field Map Y",1000,500);
//     hFieldY -> Draw("colz");
//     C2 -> Print("FieldY.gif+5","gif+5");
//     TCanvas * C3 = new TCanvas("Field Map Z","Field Map Z",1000,500);
//     hFieldZ -> Draw("colz");
//     C3 -> Print("FieldZ.gif+5","gif+5");
//     
//     delete C1;
//     delete C2;
//     delete C3;
//   }
//   delete hFieldX;
//   delete hFieldY;
//   delete hFieldZ;
// }
// 
// void FillHisto ()
// { 
//   TH3F *DistMapX = new TH3F("Distortion_Field_X","Distortion Field X",kResolution[0],0,kDetector[0],kResolution[1],0,kDetector[1],kResolution[2],0,kDetector[2]);
//   TH3F *DistMapY = new TH3F("Distortion_Field_Y","Distortion Field Y",kResolution[0],0,kDetector[0],kResolution[1],0,kDetector[1],kResolution[2],0,kDetector[2]);
//   TH3F *DistMapZ = new TH3F("Distortion_Field_Z","Distortion Field Z",kResolution[0],0,kDetector[0],kResolution[1],0,kDetector[1],kResolution[2],0,kDetector[2]);
//   
//   for (int xbin = 0; xbin < kResolution[0]; xbin++) for (int ybin = 0; ybin < kResolution[1]; ybin++) for(int zbin = 0; zbin < kResolution[2]; zbin++)
//   {
//     DistMapX -> SetBinContent(xbin+1,ybin+1,zbin+1,fRecoField[0][xbin][ybin][zbin]);
//     DistMapY -> SetBinContent(xbin+1,ybin+1,zbin+1,fRecoField[1][xbin][ybin][zbin]);
//     DistMapZ -> SetBinContent(xbin+1,ybin+1,zbin+1,fRecoField[2][xbin][ybin][zbin]);
//   }
//   TFile *OutputFile = new TFile("Field.root", "recreate");
//   DistMapX -> Write();
//   DistMapY -> Write();
//   DistMapZ -> Write(); 
//   OutputFile -> Close();
//   
//   delete DistMapX;
//   delete DistMapY;
//   delete DistMapZ;
// }
