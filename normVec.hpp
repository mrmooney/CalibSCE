#include "mdefs.h"

using namespace std;

vector<Double_t> getNormVec(Double_t xA, Double_t yA, Double_t zA, Double_t xB, Double_t yB, Double_t zB, Double_t xC, Double_t yC, Double_t zC)
{
  vector<Double_t> normVec;

  if(xA == 0.0)
    xA = 0.00000001;
  if(yA == 0.0)
    yA = 0.00000001;
  if(zA == 0.0)
    zA = 0.00000001;
  if(xB == 0.0)
    xB = 0.00000001;
  if(yB == 0.0)
    yB = 0.00000001;
  if(zB == 0.0)
    zB = 0.00000001;
  if(xC == 0.0)
    xC = 0.00000001;
  if(yC == 0.0)
    yC = 0.00000001;
  if(zC == 0.0)
    zC = 0.00000001;

  Double_t a = -((-(yC*zB) + yA*zB + yB*zC - yA*zC - yB*zA + yC*zA)/
        (xA*yC*zB - xC*yA*zB - xA*yB*zC + xB*yA*zC + xC*yB*zA - xB*yC*zA));
  Double_t b = -((xC*zB - xA*zB - xB*zC + xA*zC + xB*zA - xC*zA)/
        (xA*yC*zB - xC*yA*zB - xA*yB*zC + xB*yA*zC + xC*yB*zA - xB*yC*zA));
  Double_t c = -((xC*yB - xA*yB - xB*yC + xA*yC + xB*yA - xC*yA)/
        (-(xA*yC*zB) + xC*yA*zB + xA*yB*zC - xB*yA*zC - xC*yB*zA + xB*yC*zA));

  Double_t x1;
  Double_t y1;
  Double_t x2;
  Double_t y2;

  if(xB < xA)
  {
    x1 = xB-xA;
    y1 = yB-yA;
    x2 = xC-xA;
    y2 = yC-yA;
  }
  else
  {
    x1 = xC-xA;
    y1 = yC-yA;
    x2 = xB-xA;
    y2 = yB-yA;
  }

  Double_t slope = -(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2)/(3.*(Power(x2,2)*y1 - Power(x1,2)*y2)) - 
      (Power(2,0.3333333333333333)*(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
           3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2))))/
       (3.*(Power(x2,2)*y1 - Power(x1,2)*y2)*Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 
           6*Power(x1,4)*Power(x2,5) - 2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
           18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
           18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
           9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
           18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
           6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
           27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
           18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
           45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
           2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
           30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 27*Power(x1,5)*Power(y2,4) - 
           18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
           Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
               2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 
               27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 18*Power(x1,4)*Power(x2,3)*y1*y2 - 
               18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 
               45*x1*Power(x2,4)*Power(y1,3)*y2 - 18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 
               18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 
               27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 
               6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 
               18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 45*Power(x1,4)*x2*y1*Power(y2,3) + 
               9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 
               30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 
               2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
               18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
             4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
          0.3333333333333333)) + Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
         2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 
         27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 18*Power(x1,4)*Power(x2,3)*y1*y2 - 
         18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 
         45*x1*Power(x2,4)*Power(y1,3)*y2 - 18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
         6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
         27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
         18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
         45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 
         30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 
         2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
         18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
         Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
             2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 
             27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 18*Power(x1,4)*Power(x2,3)*y1*y2 - 
             18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 
             45*x1*Power(x2,4)*Power(y1,3)*y2 - 18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 
             18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 
             27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 
             6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 
             18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 45*Power(x1,4)*x2*y1*Power(y2,3) + 
             9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 
             30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 
             2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
             18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
           4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
              3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
				       0.3333333333333333)/(3.*Power(2,0.3333333333333333)*(Power(x2,2)*y1 - Power(x1,2)*y2));

  Double_t compactness = (-(Power(y1,2)*y2) + y1*Power(y2,2) + 2*x1*y1*y2*
         (-(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2)/(3.*(Power(x2,2)*y1 - Power(x1,2)*y2)) - 
           (Power(2,0.3333333333333333)*(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2))))/
            (3.*(Power(x2,2)*y1 - Power(x1,2)*y2)*Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 
                6*Power(x1,4)*Power(x2,5) - 2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
                Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
                    2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                    18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                    18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                    9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                    18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                    6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                    27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                    18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                    45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                    2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                    30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                    27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                    18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
                  4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                     3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
               0.3333333333333333)) + Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 
              6*Power(x1,4)*Power(x2,5) - 2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
              18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
              18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
              9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
              18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
              6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
              27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
              18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
              45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
              2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
              30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 27*Power(x1,5)*Power(y2,4) - 
              18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
              Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
                  2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                  18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                  18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                  9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                  18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                  6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                  27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                  18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                  45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                  2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                  30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                  27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                  18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
                4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                   3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
             0.3333333333333333)/(3.*Power(2,0.3333333333333333)*(Power(x2,2)*y1 - Power(x1,2)*y2))) - 
        2*x2*y1*y2*(-(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2)/(3.*(Power(x2,2)*y1 - Power(x1,2)*y2)) - 
           (Power(2,0.3333333333333333)*(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2))))/
            (3.*(Power(x2,2)*y1 - Power(x1,2)*y2)*Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 
                6*Power(x1,4)*Power(x2,5) - 2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
                Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
                    2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                    18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                    18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                    9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                    18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                    6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                    27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                    18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                    45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                    2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                    30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                    27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                    18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
                  4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                     3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
               0.3333333333333333)) + Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 
              6*Power(x1,4)*Power(x2,5) - 2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
              18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
              18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
              9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
              18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
              6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
              27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
              18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
              45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
              2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
              30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 27*Power(x1,5)*Power(y2,4) - 
              18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
              Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
                  2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                  18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                  18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                  9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                  18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                  6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                  27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                  18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                  45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                  2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                  30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                  27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                  18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
                4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                   3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
             0.3333333333333333)/(3.*Power(2,0.3333333333333333)*(Power(x2,2)*y1 - Power(x1,2)*y2))) + 
        Power(x2,2)*y1*Power(-(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2)/
            (3.*(Power(x2,2)*y1 - Power(x1,2)*y2)) - 
           (Power(2,0.3333333333333333)*(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2))))/
            (3.*(Power(x2,2)*y1 - Power(x1,2)*y2)*Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 
                6*Power(x1,4)*Power(x2,5) - 2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
                Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
                    2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                    18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                    18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                    9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                    18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                    6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                    27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                    18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                    45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                    2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                    30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                    27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                    18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
                  4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                     3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
               0.3333333333333333)) + Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 
              6*Power(x1,4)*Power(x2,5) - 2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
              18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
              18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
              9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
              18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
              6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
              27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
              18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
              45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
              2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
              30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 27*Power(x1,5)*Power(y2,4) - 
              18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
              Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
                  2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                  18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                  18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                  9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                  18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                  6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                  27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                  18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                  45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                  2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                  30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                  27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                  18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
                4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                   3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
             0.3333333333333333)/(3.*Power(2,0.3333333333333333)*(Power(x2,2)*y1 - Power(x1,2)*y2)),2) - 
        Power(x1,2)*y2*Power(-(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2)/
            (3.*(Power(x2,2)*y1 - Power(x1,2)*y2)) - 
           (Power(2,0.3333333333333333)*(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2))))/
            (3.*(Power(x2,2)*y1 - Power(x1,2)*y2)*Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 
                6*Power(x1,4)*Power(x2,5) - 2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
                Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
                    2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                    18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                    18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                    9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                    18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                    6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                    27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                    18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                    45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                    2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                    30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                    27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                    18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
                  4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                     3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
               0.3333333333333333)) + Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 
              6*Power(x1,4)*Power(x2,5) - 2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
              18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
              18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
              9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
              18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
              6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
              27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
              18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
              45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
              2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
              30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 27*Power(x1,5)*Power(y2,4) - 
              18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4) + 
              Sqrt(Power(2*Power(x1,6)*Power(x2,3) - 6*Power(x1,5)*Power(x2,4) + 6*Power(x1,4)*Power(x2,5) - 
                  2*Power(x1,3)*Power(x2,6) - 18*Power(x1,3)*Power(x2,4)*Power(y1,2) + 
                  18*Power(x1,2)*Power(x2,5)*Power(y1,2) + 27*Power(x2,5)*Power(y1,4) + 6*Power(x1,5)*Power(x2,2)*y1*y2 + 
                  18*Power(x1,4)*Power(x2,3)*y1*y2 - 18*Power(x1,3)*Power(x2,4)*y1*y2 - 6*Power(x1,2)*Power(x2,5)*y1*y2 - 
                  9*Power(x1,2)*Power(x2,3)*Power(y1,3)*y2 - 45*x1*Power(x2,4)*Power(y1,3)*y2 - 
                  18*Power(x1,5)*Power(x2,2)*Power(y2,2) + 18*Power(x1,4)*Power(x2,3)*Power(y2,2) + 
                  6*Power(x1,4)*x2*Power(y1,2)*Power(y2,2) - 27*Power(x1,3)*Power(x2,2)*Power(y1,2)*Power(y2,2) + 
                  27*Power(x1,2)*Power(x2,3)*Power(y1,2)*Power(y2,2) - 6*x1*Power(x2,4)*Power(y1,2)*Power(y2,2) - 
                  18*x1*Power(x2,2)*Power(y1,4)*Power(y2,2) + 18*Power(x2,3)*Power(y1,4)*Power(y2,2) + 
                  45*Power(x1,4)*x2*y1*Power(y2,3) + 9*Power(x1,3)*Power(x2,2)*y1*Power(y2,3) + 
                  2*Power(x1,3)*Power(y1,3)*Power(y2,3) + 30*Power(x1,2)*x2*Power(y1,3)*Power(y2,3) - 
                  30*x1*Power(x2,2)*Power(y1,3)*Power(y2,3) - 2*Power(x2,3)*Power(y1,3)*Power(y2,3) - 
                  27*Power(x1,5)*Power(y2,4) - 18*Power(x1,3)*Power(y1,2)*Power(y2,4) + 
                  18*Power(x1,2)*x2*Power(y1,2)*Power(y2,4),2) + 
                4*Power(-Power(-(Power(x1,2)*x2) + x1*Power(x2,2) + 2*x1*y1*y2 - 2*x2*y1*y2,2) + 
                   3*(Power(x2,2)*y1 - Power(x1,2)*y2)*(2*x1*x2*y1 - 2*x1*x2*y2 - Power(y1,2)*y2 + y1*Power(y2,2)),3)),
					    0.3333333333333333)/(3.*Power(2,0.3333333333333333)*(Power(x2,2)*y1 - Power(x1,2)*y2)),2))/(-(x2*y1) + x1*y2);

  if(slope == TMath::Infinity())
    slope = 10000000.0;
  if(slope == -1.0*TMath::Infinity())
    slope = -10000000.0;
  if(compactness == TMath::Infinity())
    compactness = 10000000.0;
  if(compactness == -1.0*TMath::Infinity())
    compactness = -10000000.0;

  if((TMath::IsNaN(slope) == 1) && (fabs(x1) < 0.0000001) && (fabs(x2) < 0.0000001))
  {
    slope = 0.0000001;
    compactness = 10000000.0;
  }
  else if((TMath::IsNaN(slope) == 1) && (fabs(y1) < 0.0000001) && (fabs(y2) < 0.0000001))
  {
    slope = 1000000.0;
    compactness = 100000000000.0;
  }

  Double_t aFactor;
  Double_t bFactor;
  if((slope > 0.0) && (y1 > (-1.0/slope)*x1))
  {
    aFactor = 1.0;
    bFactor = 1.0;
  }
  else if((slope > 0.0) && (y1 < (-1.0/slope)*x1))
  {
    aFactor = -1.0;
    bFactor = -1.0;
  }
  else if((slope < 0.0) && (y1 > (-1.0/slope)*x1))
  {
    aFactor = -1.0;
    bFactor = -1.0;
  }
  else if((slope < 0.0) && (y1 < (-1.0/slope)*x1))
  {
    aFactor = 1.0;
    bFactor = 1.0;
  }
  else
  {
    aFactor = 1.0;
    bFactor = 1.0;
  }

  Double_t zSlope = (1.0 - aFactor*a*1.0 - bFactor*b*slope)/c;
  Double_t norm = sqrt(1.0 + pow(slope,2) + pow(zSlope-zA,2));

  normVec.push_back(aFactor/norm);
  normVec.push_back((bFactor*slope)/norm);
  normVec.push_back((zSlope-zA)/norm);

  return normVec;
}
