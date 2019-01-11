#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "MELAHXSWidth.h"
#include "TROOT.h"
#include "TF1.h"

using namespace std;


MELAHXSWidth::MELAHXSWidth(std::string fileLoc, std::string strAppend) :
  xmhW(nullptr),
  sigW(nullptr),
  graphW(nullptr),
  gsW(nullptr)
{
  fileName = fileLoc + "/HiggsTotalWidth_" + strAppend + ".txt";
  build();
}
MELAHXSWidth::MELAHXSWidth(const MELAHXSWidth& other) :
  fileName(other.fileName),
  xmhW(nullptr),
  sigW(nullptr),
  graphW(nullptr),
  gsW(nullptr)
{
  build();
}
void MELAHXSWidth::build(){
  ifstream file;
  file.open(fileName.c_str());
  while (!file.eof()){
    double mass=0, br=0;
    file >> mass >> br;
    if (mass>0. && br>0.){
      mass_BR.push_back(mass);
      BR.push_back(br);
    }
  }
  file.close();
  const unsigned int indexW = mass_BR.size();
  if (indexW>1){
    xmhW = new double[indexW];
    sigW = new double[indexW];
    for (unsigned int ix=0; ix<indexW; ix++){
      xmhW[ix] = mass_BR.at(ix);
      sigW[ix] = BR.at(ix);
    }
    double dbegin = (sigW[1]-sigW[0])/(xmhW[1]-xmhW[0]);
    double cB = (sigW[indexW-1]-sigW[indexW-2])/(pow(xmhW[indexW-1], 3)-pow(xmhW[indexW-2], 3));
    double dend = 3.*cB*pow(xmhW[indexW-1], 2);
    graphW = new TGraph(indexW, xmhW, sigW);
    gsW = new TSpline3("gsW", graphW, "b1e1", dbegin, dend);
  }
}


MELAHXSWidth::~MELAHXSWidth(){
  delete gsW; gsW=nullptr;
  delete graphW; graphW=nullptr;
  delete[] xmhW; xmhW=nullptr;
  delete[] sigW; sigW=nullptr;
}

double MELAHXSWidth::HiggsWidth(double mH) const{
  double result = 0;
  const unsigned int indexW = mass_BR.size();
  if (gsW){
    if (mH<xmhW[indexW-1]) result = (double)gsW->Eval(mH);
    else{
      double cB = (sigW[indexW-1]-sigW[indexW-2])/(pow(xmhW[indexW-1], 3)-pow(xmhW[indexW-2], 3));
      double cA = sigW[indexW-1] - cB*pow(xmhW[indexW-1], 3);
      result = cA + cB*pow(mH, 3);
    }
  }
  return result;
}

