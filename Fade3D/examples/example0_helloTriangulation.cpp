#include <Fade_3D.h>
#include <stdio.h>

using namespace FADE3D;
using namespace std;

double getVal(Point3 point, pair<Point3,double> vtx1, pair<Point3,double> vtx2, pair<Point3,double> vtx3, pair<Point3,double> vtx4);

int example0_main()
{
	std::cout<<"\n";
	std::cout<<"example0: HelloTriangulation - 10 lines of code\n";
	std::cout<<"* Triangulate 4 points\n";
	std::cout<<"* Visualize the result\n\n";

	// Create a triangulation
	Fade_3D dt;

	// Create and insert 4 points
	Point3 p0(3.0,4.0,5.0);
	Point3 p1(4.0,4.0,5.0);
	Point3 p2(3.0,5.0,5.0);
	Point3 p3(3.0,4.0,6.0);
	dt.insert(p0);
	dt.insert(p1);
	dt.insert(p2);
	dt.insert(p3);

        map<Point3,double> themap;
	themap[p0] = 0;
	themap[p1] = 6;
	themap[p2] = 6;
	themap[p3] = 6;

	Point3 p(3.333,4.333,5.333); // 5.994
	//Point3 p(0.5,0.5,0.5);
	//Point3 p(0.0,0.0,0.0);	
        Tet3 *mytet = dt.locate(p);

	if(mytet != NULL) {
	  Point3 *pa;
	  Point3 *pb;
	  Point3 *pc;
	  Point3 *pd;
	  mytet->getCorners(pa,pb,pc,pd);
	  pair<Point3,double> v1(*pa,themap[*pa]);
	  pair<Point3,double> v2(*pb,themap[*pb]);
	  pair<Point3,double> v3(*pc,themap[*pc]);
	  pair<Point3,double> v4(*pd,themap[*pd]);
	  cout << pa->x() << " " << pa->y() << " " << pa->z() << endl;
	  
	  //pair<Point3,double> v1(p0,themap[p0]);
	  //pair<Point3,double> v2(p1,themap[p1]);
	  //pair<Point3,double> v3(p2,themap[p2]);
	  //pair<Point3,double> v4(p3,themap[p3]);
	  cout << endl << endl << "THEVAL:   " << getVal(p,v1,v2,v3,v4) << endl << endl;
	}
	else {
          cout << "HOUSTON WE HAVE A PROBLEM" << endl;
	}
	
	// Draw a file for the Geomview viewer
	dt.show("example0.list");
	return 0;
}

double getVal(Point3 point, pair<Point3,double> vtx1, pair<Point3,double> vtx2, pair<Point3,double> vtx3, pair<Point3,double> vtx4)
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

  double f1 = vtx1.second;
  double f2 = vtx2.second;
  double f3 = vtx3.second;
  double f4 = vtx4.second;

  return lam1*f1+lam2*f2+lam3*f3+lam4*f4;
}
