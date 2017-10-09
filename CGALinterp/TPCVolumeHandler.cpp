#include "TPCVolumeHandler.hpp"

TPCVolumeHandler::TPCVolumeHandler()
{
  
}

TPCVolumeHandler::TPCVolumeHandler(const ThreeVector<float>& TPCSize, const ThreeVector<float>& TPCOffset, const ThreeVector<int>& TPCResolution) : DetectorSize(TPCSize), DetectorOffset(TPCOffset), DetectorResolution(TPCResolution)
{
  NormalVector.resize(2*DetectorSize.size());
  NormalVectorOffSet.resize(2*DetectorOffset.size());
  
  CalcNormal(DetectorSize, DetectorOffset);
}

TPCVolumeHandler::TPCVolumeHandler(const std::array<float,3>& TPCSize, const std::array<float,3>& TPCOffset , const std::array<int,3>& TPCResolution) : DetectorSize(TPCSize), DetectorOffset(TPCOffset), DetectorResolution(TPCResolution)
{
  NormalVector.resize(2*DetectorSize.size());
  NormalVectorOffSet.resize(2*DetectorOffset.size());
  
  CalcNormal(DetectorSize, DetectorOffset);
}

TPCVolumeHandler::TPCVolumeHandler(float TPCSize[3], float TPCOffset[3], int TPCResolution[3]) : DetectorSize(TPCSize), DetectorOffset(TPCOffset), DetectorResolution(TPCResolution)
{
  NormalVector.resize(2*DetectorSize.size());
  NormalVectorOffSet.resize(2*DetectorOffset.size());
  
  CalcNormal(DetectorSize, DetectorOffset);
}

std::vector<ThreeVector<float>> TPCVolumeHandler::GetNormalVectors()
{
  return this->NormalVector;
}

std::vector< ThreeVector< float > > TPCVolumeHandler::GetNormVectorOffset()
{
  return this->NormalVectorOffSet;
}


void TPCVolumeHandler::CalcNormal(ThreeVector<float>& TPCSize, ThreeVector<float>& TPCOffset)
{
  std::vector<ThreeVector<float>> SingleCoordVec;
  SingleCoordVec.resize(TPCSize.size());
  
  // Split TPC size vector in single coordinate vectors
  SingleCoordVec[0] = {TPCSize[0],0.0,0.0};
  SingleCoordVec[1] = {0.0,TPCSize[1],0.0};
  SingleCoordVec[2] = {0.0,0.0,TPCSize[2]};
  
  // Fill normal vectors of all plains
  NormalVector[0] = ThreeVector<float>::VectorProduct(SingleCoordVec[0],SingleCoordVec[1]); // x-y Plane
  NormalVector[1] = NormalVector[0]; // far siede x-y Plane
  NormalVectorOffSet[0] = {0.0,0.0,0.0};
  NormalVectorOffSet[1] = SingleCoordVec[2];
  
  NormalVector[2] = ThreeVector<float>::VectorProduct(SingleCoordVec[0],SingleCoordVec[2]); // x-z Plane
  NormalVector[3] = NormalVector[2]; // far siede x-z Plane
  NormalVectorOffSet[2] = {0.0,0.0,0.0};
  NormalVectorOffSet[3] = SingleCoordVec[1];
  
  NormalVector[4] = ThreeVector<float>::VectorProduct(SingleCoordVec[1],SingleCoordVec[2]); // y-z Plane
  NormalVector[5] = NormalVector[4]; // far siede y-z Plane
  NormalVectorOffSet[4] = {0.0,0.0,0.0};
  NormalVectorOffSet[5] = SingleCoordVec[0];
  
  // coordinate system shift (zero to mid y-axis)
  for(unsigned norm_dim = 0; norm_dim < NormalVectorOffSet.size(); norm_dim++)
  {
    NormalVectorOffSet[norm_dim] += TPCOffset;
  }
  
//   for(auto NormIter : NormalVectorOffSet)
//   {
//     for(unsigned VectorIter = 0; VectorIter < NormIter.size(); VectorIter++)
//       std::cout << NormIter[VectorIter] << " ";
//     std::cout << std::endl;
//   }
}

ThreeVector<float> TPCVolumeHandler::GetDetectorSize()
{
  return this->DetectorSize;
}

ThreeVector<float> TPCVolumeHandler::GetDetectorOffset()
{
  return this->DetectorOffset;
}

ThreeVector<int> TPCVolumeHandler::GetDetectorResolution()
{
  return this->DetectorResolution;
}


