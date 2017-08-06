#include "MELANCSplineCore.h" 
#include <cmath>
#include "TMath.h"
#include "TIterator.h"
#include "Riostream.h"

using namespace TMath;
using namespace RooFit;
using namespace std;
using namespace TNumericUtil;


ClassImp(MELANCSplineCore)

MELANCSplineCore::MELANCSplineCore() :
RooAbsReal(),
verbosity(MELANCSplineCore::kSilent),
useFloor(true), floorEval(0), floorInt(0),
rangeXmin(1), rangeXmax(-1),
theXVar("theXVar", "theXVar", this),
leafDepsList("leafDepsList", "leafDepsList", this)
{}

MELANCSplineCore::MELANCSplineCore(
  const char* name,
  const char* title
  ) :
  RooAbsReal(name, title),
  verbosity(MELANCSplineCore::kSilent),
  useFloor(true), floorEval(0), floorInt(0),
  rangeXmin(1), rangeXmax(-1),
  theXVar("theXVar", "theXVar", this),
  leafDepsList("leafDepsList", "leafDepsList", this)
{}

MELANCSplineCore::MELANCSplineCore(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  const std::vector<T>& inXList,
  Bool_t inUseFloor,
  T inFloorEval,
  T inFloorInt
  ) :
  RooAbsReal(name, title),
  verbosity(MELANCSplineCore::kSilent),
  useFloor(inUseFloor), floorEval(inFloorEval), floorInt(inFloorInt),
  rangeXmin(1), rangeXmax(-1),
  theXVar("theXVar", "theXVar", this, inXVar),
  leafDepsList("leafDepsList", "leafDepsList", this),
  XList(inXList)
{}

MELANCSplineCore::MELANCSplineCore(
  const MELANCSplineCore& other,
  const char* name
  ) :
  RooAbsReal(other, name),
  verbosity(other.verbosity),
  useFloor(other.useFloor), floorEval(other.floorEval), floorInt(other.floorInt),
  rangeXmin(other.rangeXmin), rangeXmax(other.rangeXmax),
  theXVar("theXVar", this, other.theXVar),
  leafDepsList("leafDepsList", this, other.leafDepsList),
  XList(other.XList)
{}

void MELANCSplineCore::setVerbosity(VerbosityLevel flag){ verbosity = flag; }
void MELANCSplineCore::setEvalFloor(MELANCSplineCore::T val){ floorEval = val; }
void MELANCSplineCore::setIntFloor(MELANCSplineCore::T val){ floorInt = val; }
void MELANCSplineCore::doFloor(Bool_t flag){ useFloor = flag; }

void MELANCSplineCore::getBArray(const std::vector<MELANCSplineCore::T>& kappas, const vector<MELANCSplineCore::T>& fcnList, std::vector<MELANCSplineCore::T>& BArray)const{
  BArray.clear();
  int npoints=kappas.size();
  if (npoints!=(int)fcnList.size()){
    coutE(InputArguments) << "MELANCSplineCore::getBArray: Dim(kappas)=" << npoints << " != Dim(fcnList)=" << fcnList.size() << endl;
    assert(0);
  }
  if (npoints>1){
    BArray.push_back(3.*(fcnList.at(1)-fcnList.at(0)));
    for (int j=1; j<npoints-1; j++){
      MELANCSplineCore::T val_j = fcnList.at(j);
      MELANCSplineCore::T val_jpo = fcnList.at(j+1);
      MELANCSplineCore::T val_jmo = fcnList.at(j-1);
      MELANCSplineCore::T kappa_j = kappas.at(j);
      MELANCSplineCore::T kappa_jmo = kappas.at(j-1);
      MELANCSplineCore::T rsq = pow(kappa_j/kappa_jmo, 2);
      MELANCSplineCore::T Bval = val_jpo*rsq + val_j*(1.-rsq) - val_jmo;
      Bval *= 3.;
      BArray.push_back(Bval);
    }
    BArray.push_back(3.*(fcnList.at(npoints-1)-fcnList.at(npoints-2)));
  }
  else if (npoints==1) BArray.push_back(0);
}
void MELANCSplineCore::getAArray(const vector<MELANCSplineCore::T>& kappas, vector<vector<MELANCSplineCore::T>>& AArray)const{
  AArray.clear();
  Int_t npoints = kappas.size();
  for (int i=0; i<npoints; i++){
    vector<MELANCSplineCore::T> Ai(npoints, 0.);
    if (npoints==1) Ai[0]=1;
    else if (i==0){ Ai[0]=2; Ai[1]=kappas.at(1)/kappas.at(0); }
    else if (i==npoints-1){ Ai[npoints-2]=1; Ai[npoints-1]=2.*kappas.at(npoints-1)/kappas.at(npoints-2); }
    else{
      MELANCSplineCore::T kappa_j = kappas.at(i);
      MELANCSplineCore::T kappa_jmo = kappas.at(i-1);
      MELANCSplineCore::T kappa_jpo = kappas.at(i+1);

      Ai[i-1]=1;
      Ai[i]=2.*kappa_j/kappa_jmo*(1.+kappa_j/kappa_jmo);
      Ai[i+1]=kappa_j*kappa_jpo/pow(kappa_jmo, 2);
    }
    AArray.push_back(Ai);
  }
}

vector<MELANCSplineCore::T> MELANCSplineCore::getCoefficients(const TVector_t& S, const vector<MELANCSplineCore::T>& kappas, const vector<MELANCSplineCore::T>& fcnList, const Int_t& bin)const{
  DefaultAccumulator<MELANCSplineCore::T> A, B, C, D;
  vector<MELANCSplineCore::T> res;

  const int fcnsize = fcnList.size();
  if (fcnsize>bin){
    A=fcnList.at(bin);
    B=S[bin];
    if (fcnsize>(bin+1)){
      DefaultAccumulator<MELANCSplineCore::T> dFcn = fcnList.at(bin+1);
      dFcn -= A;

      C += MELANCSplineCore::T(3.)*dFcn;
      C -= 2.*B;
      C -= S[bin+1]*kappas.at(bin+1)/kappas.at(bin);

      D += MELANCSplineCore::T(-2.)*dFcn;
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
vector<vector<MELANCSplineCore::T>> MELANCSplineCore::getCoefficientsAlongDirection(const std::vector<MELANCSplineCore::T>& kappas, const TMatrix_t& Ainv, const vector<MELANCSplineCore::T>& fcnList, const Int_t pickBin)const{
  vector<MELANCSplineCore::T> BArray;
  getBArray(kappas, fcnList, BArray);

  Int_t npoints = BArray.size();
  TVector_t Btrans(npoints);
  for (int i=0; i<npoints; i++) Btrans[i]=BArray.at(i);
  TVector_t Strans = Ainv*Btrans;

  vector<vector<MELANCSplineCore::T>> coefs;
  for (Int_t bin=0; bin<(npoints>1 ? npoints-1 : 1); bin++){
    if (pickBin>=0 && bin!=pickBin) continue;
    vector<MELANCSplineCore::T> coef = getCoefficients(Strans, kappas, fcnList, bin);
    coefs.push_back(coef);
  }
  return coefs;
}

MELANCSplineCore::T MELANCSplineCore::evalSplineSegment(const std::vector<MELANCSplineCore::T>& coefs, const MELANCSplineCore::T& kappa, const MELANCSplineCore::T& tup, const MELANCSplineCore::T& tdn, Bool_t doIntegrate)const{
  DefaultAccumulator<MELANCSplineCore::T> res;
  for (unsigned int ic=0; ic<coefs.size(); ic++){
    if (doIntegrate) res += coefs.at(ic)*(pow(tup, (int)(ic+1))-pow(tdn, (int)(ic+1)))/((MELANCSplineCore::T)(ic+1));
    else res += coefs.at(ic)*pow(tup, (int)ic);
  }
  if (doIntegrate) res /= kappa;
  return res;
}

void MELANCSplineCore::getLeafDependents(RooRealProxy& proxy, RooArgSet& set){
  RooArgSet deps;
  proxy.absArg()->leafNodeServerList(&deps, 0, true);
  set.add(deps);
}
void MELANCSplineCore::addLeafDependents(RooArgSet& set){
  TIterator* iter = set.createIterator();
  RooAbsArg* absarg;
  while ((absarg = (RooAbsArg*)iter->Next())){ if (dynamic_cast<RooRealVar*>(absarg)) leafDepsList.add(*absarg); }
  delete iter;
}
