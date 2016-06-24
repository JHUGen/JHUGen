#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "TROOT.h"
#include "TF1.h"
#include <interface/MELAHXSWidth.h>

using namespace std;


MELAHXSWidth::MELAHXSWidth(std::string fileLoc, std::string strAppend) :
xmhW(0),
sigW(0),
graphW(0),
gsW(0)
{
  ifstream file;
  fileName = fileLoc + "/HiggsTotalWidth_" + strAppend + ".txt";
  file.open(fileName.c_str());
  while (!file.eof()){
    double mass=0, br=0;
    file >> mass >> br;
    if (mass>0. && br>0.){
      mass_BR.push_back(mass);
      BR.push_back(br);
      //cout << "mass, width = " << mass << '\t' << br << endl;
    }
  }
  file.close();
  int indexW = (int)mass_BR.size();
  if (indexW>0){
    double* xmhW = new double[indexW];
    double* sigW = new double[indexW];
    for (int ix=0; ix<indexW; ix++){
      xmhW[ix] = mass_BR.at(ix);
      sigW[ix] = BR.at(ix);
    }
    double dbegin = (sigW[1]-sigW[0])/(xmhW[1]-xmhW[0]);
    double dend = (sigW[indexW-1]-sigW[indexW-2])/(xmhW[indexW-1]-xmhW[indexW-2]);
    graphW = new TGraph(indexW, xmhW, sigW);
    gsW = new TSpline3("gsW", graphW, "b1e1", dbegin, dend);
  }
}


MELAHXSWidth::~MELAHXSWidth(){
  if (gsW!=0) delete gsW; gsW=0;
  if (graphW!=0) delete graphW; graphW=0;
  if (xmhW!=0) delete[] xmhW; xmhW=0;
  if (sigW!=0) delete[] sigW; sigW=0;
}

double MELAHXSWidth::HiggsWidth(double mH){
  double result = (double)gsW->Eval(mH);
  //cout << "MELAHXSWidth::HiggsWidth: mH requested = " << mH << endl;
  //cout << "MELAHXSWidth::HiggsWidth: GaH requested = " << result << endl;
  return result;
}

