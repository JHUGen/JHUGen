#include "MELANCSplinePdfCore.h" 
#include <cmath>
#include "Riostream.h" 
#include "TMath.h"

using namespace TMath;
using namespace RooFit;
using namespace std;
using namespace TNumericUtil;


MELANCSplinePdfCore::MELANCSplinePdfCore() :
RooAbsPdf(),
verbosity(MELANCSplinePdfCore::kSilent),
theXVar("theXVar", "theXVar", this)
{}

MELANCSplinePdfCore::MELANCSplinePdfCore(
  const char* name,
  const char* title
  ) :
  RooAbsPdf(name, title),
  verbosity(MELANCSplinePdfCore::kSilent),
  useFloor(true),
  floorEval(0),
  floorInt(0),
  theXVar("theXVar", "theXVar", this)
{}

MELANCSplinePdfCore::MELANCSplinePdfCore(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  const std::vector<T>& inXList,
  Bool_t inUseFloor,
  T inFloorEval,
  T inFloorInt
  ) :
  RooAbsPdf(name, title),
  verbosity(MELANCSplinePdfCore::kSilent),
  useFloor(inUseFloor),
  floorEval(inFloorEval),
  floorInt(inFloorInt),
  theXVar("theXVar", "theXVar", this, inXVar),
  XList(inXList)
{}

MELANCSplinePdfCore::MELANCSplinePdfCore(
  const MELANCSplinePdfCore& other,
  const char* name
  ) :
  RooAbsPdf(other, name),
  verbosity(other.verbosity),
  useFloor(other.useFloor),
  floorEval(other.floorEval),
  floorInt(other.floorInt),
  theXVar("theXVar", this, other.theXVar),
  XList(other.XList)
{}

void MELANCSplinePdfCore::setVerbosity(VerbosityLevel flag){ verbosity = flag; }
void MELANCSplinePdfCore::setEvalFloor(MELANCSplinePdfCore::T val){ floorEval = val; }
void MELANCSplinePdfCore::setIntFloor(MELANCSplinePdfCore::T val){ floorInt = val; }
void MELANCSplinePdfCore::doFloor(Bool_t flag){ useFloor = flag; }

void MELANCSplinePdfCore::getBArray(const std::vector<MELANCSplinePdfCore::T>& kappas, const vector<MELANCSplinePdfCore::T>& fcnList, std::vector<MELANCSplinePdfCore::T>& BArray)const{
  BArray.clear();
  int npoints=kappas.size();
  if (npoints!=(int)fcnList.size()){
    coutE(InputArguments) << "MELANCSplinePdfCore::getBArray: Dim(kappas)=" << npoints << " != Dim(fcnList)=" << fcnList.size() << endl;
    assert(0);
  }
  if (npoints>1){
    BArray.push_back(3.*(fcnList.at(1)-fcnList.at(0)));
    for (int j=1; j<npoints-1; j++){
      MELANCSplinePdfCore::T val_j = fcnList.at(j);
      MELANCSplinePdfCore::T val_jpo = fcnList.at(j+1);
      MELANCSplinePdfCore::T val_jmo = fcnList.at(j-1);
      MELANCSplinePdfCore::T kappa_j = kappas.at(j);
      MELANCSplinePdfCore::T kappa_jmo = kappas.at(j-1);
      MELANCSplinePdfCore::T rsq = pow(kappa_j/kappa_jmo, 2);
      MELANCSplinePdfCore::T Bval = val_jpo*rsq + val_j*(1.-rsq) - val_jmo;
      Bval *= 3.;
      BArray.push_back(Bval);
    }
    BArray.push_back(3.*(fcnList.at(npoints-1)-fcnList.at(npoints-2)));
  }
  else if (npoints==1) BArray.push_back(0);
}
void MELANCSplinePdfCore::getAArray(const vector<MELANCSplinePdfCore::T>& kappas, vector<vector<MELANCSplinePdfCore::T>>& AArray)const{
  AArray.clear();
  Int_t npoints = kappas.size();
  for (int i=0; i<npoints; i++){
    vector<MELANCSplinePdfCore::T> Ai(npoints, 0.);
    if (npoints==1) Ai[0]=1;
    else if (i==0){ Ai[0]=2; Ai[1]=kappas.at(1)/kappas.at(0); }
    else if (i==npoints-1){ Ai[npoints-2]=1; Ai[npoints-1]=2.*kappas.at(npoints-1)/kappas.at(npoints-2); }
    else{
      MELANCSplinePdfCore::T kappa_j = kappas.at(i);
      MELANCSplinePdfCore::T kappa_jmo = kappas.at(i-1);
      MELANCSplinePdfCore::T kappa_jpo = kappas.at(i+1);

      Ai[i-1]=1;
      Ai[i]=2.*kappa_j/kappa_jmo*(1.+kappa_j/kappa_jmo);
      Ai[i+1]=kappa_j*kappa_jpo/pow(kappa_jmo, 2);
    }
    AArray.push_back(Ai);
  }
}

vector<MELANCSplinePdfCore::T> MELANCSplinePdfCore::getCoefficients(const TVector_t& S, const vector<MELANCSplinePdfCore::T>& kappas, const vector<MELANCSplinePdfCore::T>& fcnList, const Int_t& bin)const{
  DefaultAccumulator<MELANCSplinePdfCore::T> A, B, C, D;
  vector<MELANCSplinePdfCore::T> res;

  const int fcnsize = fcnList.size();
  if (fcnsize>bin){
    A=fcnList.at(bin);
    B=S[bin];
    if (fcnsize>(bin+1)){
      DefaultAccumulator<MELANCSplinePdfCore::T> dFcn = fcnList.at(bin+1);
      dFcn -= A;

      C += MELANCSplinePdfCore::T(3.)*dFcn;
      C -= 2.*B;
      C -= S[bin+1]*kappas.at(bin+1)/kappas.at(bin);

      D += MELANCSplinePdfCore::T(-2.)*dFcn;
      D += B;
      D += S[bin+1]*kappas.at(bin+1)/kappas.at(bin);
    }
  }

  res.push_back(A);
  res.push_back(B);
  res.push_back(C);
  res.push_back(D);
  return res;
}
vector<vector<MELANCSplinePdfCore::T>> MELANCSplinePdfCore::getCoefficientsAlongDirection(const std::vector<MELANCSplinePdfCore::T>& kappas, const TMatrix_t& Ainv, const vector<MELANCSplinePdfCore::T>& fcnList, const Int_t pickBin)const{
  vector<MELANCSplinePdfCore::T> BArray;
  getBArray(kappas, fcnList, BArray);

  Int_t npoints = BArray.size();
  TVector_t Btrans(npoints);
  for (int i=0; i<npoints; i++) Btrans[i]=BArray.at(i);
  TVector_t Strans = Ainv*Btrans;

  vector<vector<MELANCSplinePdfCore::T>> coefs;
  for (Int_t bin=0; bin<(npoints>1 ? npoints-1 : 1); bin++){
    if (pickBin>=0 && bin!=pickBin) continue;
    vector<MELANCSplinePdfCore::T> coef = getCoefficients(Strans, kappas, fcnList, bin);
    coefs.push_back(coef);
  }
  return coefs;
}

MELANCSplinePdfCore::T MELANCSplinePdfCore::evalSplineSegment(const std::vector<MELANCSplinePdfCore::T>& coefs, const MELANCSplinePdfCore::T& kappa, const MELANCSplinePdfCore::T& tup, const MELANCSplinePdfCore::T& tdn, Bool_t doIntegrate)const{
  DefaultAccumulator<MELANCSplinePdfCore::T> res;
  for (unsigned int ic=0; ic<coefs.size(); ic++){
    if (doIntegrate) res += coefs.at(ic)*(pow(tup, (int)(ic+1))-pow(tdn, (int)(ic+1)))/((MELANCSplinePdfCore::T)(ic+1));
    else res += coefs.at(ic)*pow(tup, (int)ic);
  }
  if (doIntegrate) res /= kappa;
  return res;
}
