/** Example of a script to link the package and compute LD from bare root
 * 
 * Usage: to run as an interpreted script:
 *
 * root -q -b testKD.C
 *
 * to run as a compiled script, you also need to load loadMELA.C:
 *
 * root -q -b loadMELA.C testKD.C+
 *
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#include "ZZMatrixElement/MELA/interface/Mela.h"
#endif

#include <iostream>
using namespace std;


void testKD() 
{

#ifdef __CINT__
  gROOT->Macro("loadMELA.C");
#endif

  Mela mela(true);

  float kd, ps, pb;

  mela.computeKD(210.437,91.872,51.459,0.981,0.175,0.934,-1.980,0.269,kd,ps,pb);
 
  cout << endl << "KD=" << kd << endl;
  
}
