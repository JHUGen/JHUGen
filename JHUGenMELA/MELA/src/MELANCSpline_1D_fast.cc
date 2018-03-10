#include "MELANCSpline_1D_fast.h" 
#include <cmath>
#include "TMath.h"
#include "Riostream.h" 
#include "RooAbsReal.h" 

using namespace TMath;
using namespace RooFit;
using namespace std;
using namespace TNumericUtil;


ClassImp(MELANCSpline_1D_fast)

MELANCSpline_1D_fast::MELANCSpline_1D_fast() :
  MELANCSplineCore(),
  bcBeginX(MELANCSplineCore::bcNaturalSpline), bcEndX(MELANCSplineCore::bcNaturalSpline)
{}

MELANCSpline_1D_fast::MELANCSpline_1D_fast(
  const char* name,
  const char* title
  ) :
  MELANCSplineCore(name, title),
  bcBeginX(MELANCSplineCore::bcNaturalSpline), bcEndX(MELANCSplineCore::bcNaturalSpline)
{}

MELANCSpline_1D_fast::MELANCSpline_1D_fast(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  const std::vector<T>& inXList,
  const std::vector<T>& inFcnList,
  MELANCSplineCore::BoundaryCondition const bcBeginX_,
  MELANCSplineCore::BoundaryCondition const bcEndX_,
  Bool_t inUseFloor,
  T inFloorEval,
  T inFloorInt
  ) :
  MELANCSplineCore(name, title, inXVar, inXList, inUseFloor, inFloorEval, inFloorInt),
  bcBeginX(bcBeginX_), bcEndX(bcEndX_),
  FcnList(inFcnList)
{
  if (npointsX()>1){
    int npoints;

    vector<vector<MELANCSplineCore::T>> xA; getKappas(kappaX, 0); getAArray(kappaX, xA, bcBeginX, bcEndX);
    npoints=kappaX.size();
    TMatrix_t xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    Double_t det=0;
    TMatrix_t xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "MELANCSpline_1D::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    coefficients = getCoefficientsAlongDirection(kappaX, xAinv, FcnList, bcBeginX, bcEndX, -1);
  }
  else assert(0);

  RooArgSet leafset;
  getLeafDependents(theXVar, leafset);
  addLeafDependents(leafset);

  emptyFcnList();
}

MELANCSpline_1D_fast::MELANCSpline_1D_fast(
  const MELANCSpline_1D_fast& other,
  const char* name
  ) :
  MELANCSplineCore(other, name),
  bcBeginX(other.bcBeginX), bcEndX(other.bcEndX),
  FcnList(other.FcnList),
  kappaX(other.kappaX),
  coefficients(other.coefficients)
{}

MELANCSplineCore::T MELANCSpline_1D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  DefaultAccumulator<MELANCSplineCore::T> res;

  if (verbosity==MELANCSplineCore::kVerbose) cout << "MELANCSpline_1D_fast(" << GetName() << ")::interpolateFcn begin with code: " << code << endl;

  // Get bins
  Int_t xbin=-1, xbinmin=-1, xbinmax=-1;
  MELANCSplineCore::T tx=0, txmin=0, txmax=0;
  if (code==0 || code%2!=0){ // Case to just compute the value at x
    if (!testRangeValidity(theXVar)) return 0;
    xbin = getWhichBin(theXVar, 0);
    tx = getTVar(kappaX, theXVar, xbin, 0);
  }
  else{ // Case to integrate along x
    MELANCSplineCore::T coordmin = theXVar.min(rangeName); cropValueForRange(coordmin);
    MELANCSplineCore::T coordmax = theXVar.max(rangeName); cropValueForRange(coordmax);
    xbinmin = getWhichBin(coordmin, 0);
    txmin = getTVar(kappaX, coordmin, xbinmin, 0);
    xbinmax = getWhichBin(coordmax, 0);
    txmax = getTVar(kappaX, coordmax, xbinmax, 0);
  }

  int nxbins = (int)coefficients.size();
  for (int ix=0; ix<nxbins; ix++){
    if (
      (xbin>=0 && ix!=xbin)
      ||
      (xbinmin>=0 && xbinmax>=xbinmin && !(xbinmin<=ix && ix<=xbinmax))
      ) continue;

    MELANCSplineCore::T txlow=0, txhigh=1;
    if (code>0 && code%2==0){
      if (ix==xbinmin) txlow=txmin;
      if (ix==xbinmax) txhigh=txmax;
    }
    else txhigh=tx;

    // Get the x coefficients at bin ix and evaluate value of spline at x
    res += evalSplineSegment(coefficients.at(ix), kappaX.at(ix), txhigh, txlow, (code>0 && code%2==0));
  }

  return res;
}

void MELANCSpline_1D_fast::getKappas(vector<MELANCSplineCore::T>& kappas, const Int_t /*whichDirection*/){
  kappas.clear();
  MELANCSplineCore::T kappa=1;

  Int_t npoints;
  vector<MELANCSplineCore::T> const* coord;
  npoints=npointsX();
  coord=&XList;

  for (Int_t j=0; j<npoints-1; j++){
    MELANCSplineCore::T val_j = coord->at(j);
    MELANCSplineCore::T val_jpo = coord->at(j+1);
    MELANCSplineCore::T val_diff = (val_jpo-val_j);
    if (fabs(val_diff)>MELANCSplineCore::T(0)) kappa = 1./val_diff;
    else kappa = 0;
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
Int_t MELANCSpline_1D_fast::getWhichBin(const MELANCSplineCore::T& val, const Int_t /*whichDirection*/)const{
  Int_t bin=-1;
  MELANCSplineCore::T valj, valjpo;
  Int_t npoints;
  vector<MELANCSplineCore::T> const* coord;
  coord=&XList;
  npoints=npointsX();

  if (npoints<=1) bin=0;
  else{
    valjpo = coord->at(0);
    for (Int_t j=0; j<npoints-1; j++){
      valj = coord->at(j);
      valjpo = coord->at(j+1);
      if (val<valjpo && val>=valj){ bin=j; break; }
    }
    if (bin==-1 && val>=valjpo) bin=npoints-2;
    else if (bin==-1) bin=0;
  }

  return bin;
}
MELANCSplineCore::T MELANCSpline_1D_fast::getTVar(const vector<MELANCSplineCore::T>& kappas, const MELANCSplineCore::T& val, const Int_t& bin, const Int_t /*whichDirection*/)const{
  const MELANCSplineCore::T& K=kappas.at(bin);
  return (val-XList.at(bin))*K;
}

Double_t MELANCSpline_1D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (useFloor && value<floorEval){
    if (verbosity>=MELANCSplineCore::kError) coutE(Eval) << "MELANCSpline_1D_fast ERROR::MELANCSpline_1D_fast(" << GetName() << ") evaluation returned " << value << " at x = " << theXVar << endl;
    value = floorEval;
  }
  if (verbosity==MELANCSplineCore::kVerbose){
    cout << "MELANCSpline_1D_fast(" << GetName() << ")::evaluate = " << value << " at x = " << theXVar << endl;
    RooArgSet Xdeps; theXVar.absArg()->leafNodeServerList(&Xdeps, 0, true);
    TIterator* iter = Xdeps.createIterator();
    RooAbsArg* var;
    while ((var = (RooAbsArg*)iter->Next())){
      cout << var->GetName() << " value = " << dynamic_cast<RooAbsReal*>(var)->getVal() << endl;
    }
    delete iter;
    cout << endl;
  }
  return value;
}
Int_t MELANCSpline_1D_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  if (_forceNumInt) return 0;
  Int_t code=0;
  if (dynamic_cast<RooRealVar*>(theXVar.absArg())!=0){
    if (matchArgs(allVars, analVars, theXVar)) code=2;
  }
  return code;
}
Double_t MELANCSpline_1D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (useFloor && value<floorInt){
    if (verbosity>=MELANCSplineCore::kError) coutE(Integration) << "MELANCSpline_1D_fast ERROR::MELANCSpline_1D_fast(" << GetName() << ") integration returned " << value << " for code = " << code << endl;
    value = floorInt;
  }
  if (verbosity==MELANCSplineCore::kVerbose){ cout << "MELANCSpline_1D_fast(" << GetName() << ")::analyticalIntegral = " << value << " for code = " << code << endl; }
  return value;
}

Bool_t MELANCSpline_1D_fast::testRangeValidity(const T& val, const Int_t /*whichDirection*/) const{
  const T* range[2];
  range[0] = &rangeXmin;
  range[1] = &rangeXmax;
  return (*(range[0])>*(range[1]) || (val>=*(range[0]) && val<=*(range[1])));
}
void MELANCSpline_1D_fast::setRangeValidity(const T valmin, const T valmax, const Int_t /*whichDirection*/){
  T* range[2];
  range[0] = &rangeXmin;
  range[1] = &rangeXmax;
  *(range[0])=valmin;
  *(range[1])=valmax;
}
void MELANCSpline_1D_fast::cropValueForRange(T& val, const Int_t /*whichDirection*/)const{
  if (testRangeValidity(val)) return;
  const T* range[2];
  range[0] = &rangeXmin;
  range[1] = &rangeXmax;
  if (val<*(range[0])) val = *(range[0]);
  if (val>*(range[1])) val = *(range[1]);
}
