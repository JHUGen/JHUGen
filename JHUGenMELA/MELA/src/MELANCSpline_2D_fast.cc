#include "MELANCSpline_2D_fast.h" 
#include <cmath>
#include "TMath.h"
#include "Riostream.h" 
#include "RooAbsReal.h" 

using namespace TMath;
using namespace RooFit;
using namespace std;
using namespace TNumericUtil;


ClassImp(MELANCSpline_2D_fast)

MELANCSpline_2D_fast::MELANCSpline_2D_fast() :
MELANCSplineCore(),
rangeYmin(1), rangeYmax(-1),
theYVar("theYVar", "theYVar", this)
{}

MELANCSpline_2D_fast::MELANCSpline_2D_fast(
  const char* name,
  const char* title
  ) :
  MELANCSplineCore(name, title),
  rangeYmin(1), rangeYmax(-1),
  theYVar("theYVar", "theYVar", this)
{}

MELANCSpline_2D_fast::MELANCSpline_2D_fast(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  RooAbsReal& inYVar,
  const std::vector<T>& inXList,
  const std::vector<T>& inYList,
  const std::vector<std::vector<T>>& inFcnList,
  Bool_t inUseFloor,
  T inFloorEval,
  T inFloorInt
  ) :
  MELANCSplineCore(name, title, inXVar, inXList, inUseFloor, inFloorEval, inFloorInt),
  rangeYmin(1), rangeYmax(-1),
  theYVar("theYVar", "theYVar", this, inYVar),
  YList(inYList),
  FcnList(inFcnList)
{
  if (npointsX()>1 && npointsY()>1){
    // Prepare A and kappa arrays for x and y coordinates
    int npoints;
    Double_t det;

    vector<vector<MELANCSplineCore::T>> xA; getKappas(kappaX, 0); getAArray(kappaX, xA);
    npoints=kappaX.size();
    TMatrix_t xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    det=0;
    TMatrix_t xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "MELANCSpline_2D_fast::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    vector<vector<MELANCSplineCore::T>> yA; getKappas(kappaY, 1); getAArray(kappaY, yA);
    npoints=kappaY.size();
    TMatrix_t yAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ yAtrans[i][j]=yA.at(i).at(j); } }
    det=0;
    TMatrix_t yAinv = yAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "MELANCSpline_2D_fast::interpolateFcn: Matrix yA could not be inverted. Something is wrong with the y coordinates of points!" << endl;
      assert(0);
    }

    // Get the grid of coefficients
    vector<vector<vector<MELANCSplineCore::T>>> coefsAlongY; // [Ax(y),Bx(y),Cx(y),Dx(y)][xbin][ybin]
    int npoldim=0;
    int nxbins=0;
    for (unsigned int j=0; j<npointsY(); j++){
      vector<vector<MELANCSplineCore::T>> xcoefsAtYj = getCoefficientsPerY(kappaX, xAinv, j, -1); // [ix][Ax,Bx,Cx,Dx] at each y_j
      if (j==0){
        nxbins=xcoefsAtYj.size();
        npoldim=xcoefsAtYj.at(0).size();
        for (int ipow=0; ipow<npoldim; ipow++){
          vector<vector<MELANCSplineCore::T>> dum_xycoefarray;
          for (int ix=0; ix<nxbins; ix++){
            vector<MELANCSplineCore::T> dum_ycoefarray;
            dum_xycoefarray.push_back(dum_ycoefarray);
          }
          coefsAlongY.push_back(dum_xycoefarray);
        }
      }
      if (nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()){
        coutE(InputArguments) << "MELANCSpline_2D_fast::interpolateFcn: nxbins!=(int)xcoefsAtYj.size() || npoldim!=(int)xcoefsAtYj.at(0).size()!" << endl;
        assert(0);
      }
      for (int ix=0; ix<nxbins; ix++){
        for (int ipow=0; ipow<npoldim; ipow++) coefsAlongY.at(ipow).at(ix).push_back(xcoefsAtYj.at(ix).at(ipow));
      }
    }

    for (int ix=0; ix<nxbins; ix++){
      // Get the x coefficients interpolated across y
      vector<vector<vector<MELANCSplineCore::T>>> xCoefs;
      for (int ic=0; ic<npoldim; ic++){
        vector<vector<MELANCSplineCore::T>> yCoefs = getCoefficientsAlongDirection(kappaY, yAinv, coefsAlongY.at(ic).at(ix), -1); // [iy][A,B,C,D]
        xCoefs.push_back(yCoefs);
      }
      coefficients.push_back(xCoefs);
    }
  }
  else assert(0);

  RooArgSet leafset;
  getLeafDependents(theXVar, leafset);
  getLeafDependents(theYVar, leafset);
  addLeafDependents(leafset);

  emptyFcnList();
}

MELANCSpline_2D_fast::MELANCSpline_2D_fast(
  const MELANCSpline_2D_fast& other,
  const char* name
  ) :
  MELANCSplineCore(other, name),
  rangeYmin(other.rangeYmin), rangeYmax(other.rangeYmax),
  theYVar("theYVar", this, other.theYVar),
  YList(other.YList),
  FcnList(other.FcnList),
  kappaX(other.kappaX),
  kappaY(other.kappaY),
  coefficients(other.coefficients)
{}


MELANCSplineCore::T MELANCSpline_2D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  DefaultAccumulator<MELANCSplineCore::T> res;

  if (verbosity==MELANCSplineCore::kVerbose){ cout << "MELANCSpline_2D_fast(" << GetName() << ")::interpolateFcn begin with code: " << code << endl; }

  // Get bins
  Int_t xbin=-1, xbinmin=-1, xbinmax=-1, ybin=-1, ybinmin=-1, ybinmax=-1;
  MELANCSplineCore::T tx=0, txmin=0, txmax=0, ty=0, tymin=0, tymax=0;
  if (code==0 || code%2!=0){ // Case to just compute the value at x
    if (!testRangeValidity(theXVar, 0)) return 0;
    xbin = getWhichBin(theXVar, 0);
    tx = getTVar(kappaX, theXVar, xbin, 0);
  }
  else{ // Case to integrate along x
    MELANCSplineCore::T coordmin = theXVar.min(rangeName); cropValueForRange(coordmin, 0);
    MELANCSplineCore::T coordmax = theXVar.max(rangeName); cropValueForRange(coordmax, 0);
    xbinmin = getWhichBin(coordmin, 0);
    txmin = getTVar(kappaX, coordmin, xbinmin, 0);
    xbinmax = getWhichBin(coordmax, 0);
    txmax = getTVar(kappaX, coordmax, xbinmax, 0);
  }
  if (code==0 || code%3!=0){ // Case to just compute the value at y
    if (!testRangeValidity(theYVar, 1)) return 0;
    ybin = getWhichBin(theYVar, 1);
    ty = getTVar(kappaY, theYVar, ybin, 1);
  }
  else{ // Case to integrate along y
    MELANCSplineCore::T coordmin = theYVar.min(rangeName); cropValueForRange(coordmin, 1);
    MELANCSplineCore::T coordmax = theYVar.max(rangeName); cropValueForRange(coordmax, 1);
    ybinmin = getWhichBin(coordmin, 1);
    tymin = getTVar(kappaY, coordmin, ybinmin, 1);
    ybinmax = getWhichBin(coordmax, 1);
    tymax = getTVar(kappaY, coordmax, ybinmax, 1);
  }

  for (int ix=0; ix<(int)coefficients.size(); ix++){
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

    if (verbosity==MELANCSplineCore::kVerbose){
      if (code==0 || code%2!=0) cout << "Evaluating tx=" << txhigh << " in bin " << ix << endl;
      else cout << "Evaluating tx[" << txlow << ", " << txhigh << "] in bin " << ix << endl;
    }

    // Get the x coefficients interpolated across y
    vector<MELANCSplineCore::T> xCoefs;
    for (int ic=0; ic<(int)coefficients.at(ix).size(); ic++){
      const vector<vector<MELANCSplineCore::T>>& yCoefs = coefficients.at(ix).at(ic);

      if (verbosity==MELANCSplineCore::kVerbose) cout << "\tCoefficient " << ic << ":\n";

      DefaultAccumulator<MELANCSplineCore::T> theCoef;
      for (int iy=0; iy<(int)yCoefs.size(); iy++){
        if (
          (ybin>=0 && iy!=ybin)
          ||
          (ybinmin>=0 && ybinmax>=ybinmin && !(ybinmin<=iy && iy<=ybinmax))
          ) continue;

        MELANCSplineCore::T tylow=0, tyhigh=1;
        if (code>0 && code%3==0){
          if (iy==ybinmin) tylow=tymin;
          if (iy==ybinmax) tyhigh=tymax;
        }
        else tyhigh=ty;

        if (verbosity==MELANCSplineCore::kVerbose){
          if (code==0 || code%3!=0) cout << "\tEvaluating ty=" << tyhigh << " in bin " << iy << endl;
          else cout << "\tEvaluating ty[" << tylow << ", " << tyhigh << "] in bin " << iy << endl;
        }

        theCoef += evalSplineSegment(yCoefs.at(iy), kappaY.at(iy), tyhigh, tylow, (code>0 && code%3==0));
      }

      //if (code==0) cout << "\tCoefficient is " << theCoef << endl;

      xCoefs.push_back(theCoef);
    }

    // Evaluate value of spline at x with coefficients evaluated at y
    res += evalSplineSegment(xCoefs, kappaX.at(ix), txhigh, txlow, (code>0 && code%2==0));
  }

  return res;
}

void MELANCSpline_2D_fast::getKappas(vector<MELANCSplineCore::T>& kappas, const Int_t whichDirection){
  kappas.clear();
  MELANCSplineCore::T kappa=1;

  Int_t npoints;
  vector<MELANCSplineCore::T> const* coord;
  if (whichDirection==0){
    npoints=npointsX();
    coord=&XList;
  }
  else{
    npoints=npointsY();
    coord=&YList;
  }

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
Int_t MELANCSpline_2D_fast::getWhichBin(const MELANCSplineCore::T& val, const Int_t whichDirection)const{
  Int_t bin=-1;
  MELANCSplineCore::T valj, valjpo;
  Int_t npoints;
  vector<MELANCSplineCore::T> const* coord;
  if (whichDirection==0){
    coord=&XList;
    npoints=npointsX();
  }
  else{
    coord=&YList;
    npoints=npointsY();
  }

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
MELANCSplineCore::T MELANCSpline_2D_fast::getTVar(const vector<MELANCSplineCore::T>& kappas, const MELANCSplineCore::T& val, const Int_t& bin, const Int_t whichDirection)const{
  const MELANCSplineCore::T& K=kappas.at(bin);
  vector<MELANCSplineCore::T> const* coord;
  if (whichDirection==0) coord=&XList;
  else coord=&YList;
  return (val-coord->at(bin))*K;
}

vector<vector<MELANCSplineCore::T>> MELANCSpline_2D_fast::getCoefficientsPerY(const std::vector<MELANCSplineCore::T>& kappaX, const TMatrix_t& xAinv, const Int_t& ybin, const Int_t xbin)const{
  vector<MELANCSplineCore::T> fcnList;
  for (unsigned int bin=0; bin<npointsX(); bin++){ fcnList.push_back(FcnList.at(ybin).at(bin)); }
  vector<vector<MELANCSplineCore::T>> coefs = getCoefficientsAlongDirection(kappaX, xAinv, fcnList, xbin);
  return coefs;
}

Double_t MELANCSpline_2D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (useFloor && value<floorEval){
    if (verbosity>=MELANCSplineCore::kError) coutE(Eval) << "MELANCSpline_2D_fast ERROR::MELANCSpline_2D_fast(" << GetName() << ") evaluation returned " << value << " at (x, y) = (" << theXVar << ", " << theYVar << ")" << endl;
    value = floorEval;
  }
  if (verbosity==MELANCSplineCore::kVerbose){ cout << "MELANCSpline_2D_fast(" << GetName() << ")::evaluate = " << value << " at (x, y) = (" << theXVar << ", " << theYVar << ")" << endl; }
  return value;
}
Int_t MELANCSpline_2D_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  if (_forceNumInt) return 0;

  Int_t code=1;

  RooArgSet Xdeps, Ydeps;
  RooRealVar* rrv_x = dynamic_cast<RooRealVar*>(theXVar.absArg());
  RooRealVar* rrv_y = dynamic_cast<RooRealVar*>(theYVar.absArg());
  if (rrv_x==0) theXVar.absArg()->leafNodeServerList(&Xdeps, 0, true);
  if (rrv_y==0) theYVar.absArg()->leafNodeServerList(&Ydeps, 0, true);

  if (rrv_x!=0){
    if (Ydeps.find(*rrv_x)==0 || rrv_y!=0){
      if (matchArgs(allVars, analVars, theXVar)) code*=2;
    }
  }
  if (rrv_y!=0){
    if (Xdeps.find(*rrv_y)==0 || rrv_x!=0){
      if (matchArgs(allVars, analVars, theYVar)) code*=3;
    }
  }

  if (code==1) code=0;
  return code;
}
Double_t MELANCSpline_2D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (useFloor && value<floorInt){
    if (verbosity>=MELANCSplineCore::kError) coutE(Integration) << "MELANCSpline_2D_fast ERROR::MELANCSpline_2D_fast(" << GetName() << ") integration returned " << value << " for code = " << code << endl;
    value = floorInt;
  }
  if (verbosity==MELANCSplineCore::kVerbose){ cout << "MELANCSpline_2D_fast(" << GetName() << ")::analyticalIntegral = " << value << " for code = " << code << endl; }
  return value;
}

Bool_t MELANCSpline_2D_fast::testRangeValidity(const T& val, const Int_t whichDirection) const{
  const T* range[2];
  if (whichDirection==0){
    range[0] = &rangeXmin;
    range[1] = &rangeXmax;
  }
  else{
    range[0] = &rangeYmin;
    range[1] = &rangeYmax;
  }
  return (*(range[0])>*(range[1]) || (val>=*(range[0]) && val<=*(range[1])));
}
void MELANCSpline_2D_fast::setRangeValidity(const T valmin, const T valmax, const Int_t whichDirection){
  T* range[2];
  if (whichDirection==0){
    range[0] = &rangeXmin;
    range[1] = &rangeXmax;
  }
  else{
    range[0] = &rangeYmin;
    range[1] = &rangeYmax;
  }
  *(range[0])=valmin;
  *(range[1])=valmax;
}
void MELANCSpline_2D_fast::cropValueForRange(T& val, const Int_t whichDirection)const{
  if (testRangeValidity(val, whichDirection)) return;
  const T* range[2];
  if (whichDirection==0){
    range[0] = &rangeXmin;
    range[1] = &rangeXmax;
  }
  else{
    range[0] = &rangeYmin;
    range[1] = &rangeYmax;
  }
  if (val<*(range[0])) val = *(range[0]);
  if (val>*(range[1])) val = *(range[1]);
}
