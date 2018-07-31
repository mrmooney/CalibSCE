// #define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <string>
#include <iomanip>
#include <chrono>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include "ThreeVector.hpp"
#include "LaserTrack.hpp"
#include "Matrix3x3.hpp"

#ifndef INTERPOLATION3D_H
#define INTERPOLATION3D_H

typedef CGAL::Exact_predicates_inexact_constructions_kernel InterpKernel;
typedef CGAL::Triangulation_3<InterpKernel> Triangulation;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, InterpKernel> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
// Use the Fast_location tag. Default or Compact_location works too.
// typedef CGAL::Delaunay_triangulation_3<InterpKernel, Tds, CGAL::Fast_location > Delaunay;
// Use the Lock for parallel computing
typedef CGAL::Delaunay_triangulation_3< InterpKernel, Tds > Delaunay;
typedef Delaunay::Point Point;

std::vector<std::pair<unsigned,float>> GetClosestTracksInfo(std::vector<LaserTrack>&,const unsigned);
std::vector<std::pair<unsigned int, unsigned int>> GetClosestLaserSample(std::vector<LaserTrack>&, const unsigned);
std::array<float,2> AnglesFromPoynting(ThreeVector<float>&);
ThreeVector<float> PoyntingFromAngles(const std::array<float,2>&);
bool PairSortFunction(std::pair<unsigned int,float>, std::pair<unsigned int,float>);
Delaunay TrackMesher(const std::vector<LaserTrack>&);
Point VectorToPoint(ThreeVector<float>&);
ThreeVector<float> PointToVector(Point&);

ThreeVector<float> InterpolateCGAL(const std::vector<LaserTrack>&, const Delaunay&, ThreeVector<float>);
void InterpolateTrack(LaserTrack& ,const std::vector<LaserTrack>& , const Delaunay&);


#endif