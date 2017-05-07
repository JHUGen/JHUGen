#include "MELANCSplinePdf_3D_fast.h" 
#include <cmath>
#include "TMath.h"
#include "Riostream.h" 
#include "RooAbsReal.h" 

using namespace TMath;
using namespace RooFit;
using namespace std;
using namespace TNumericUtil;


MELANCSplinePdf_3D_fast::MELANCSplinePdf_3D_fast() :
MELANCSplinePdfCore(),
theYVar("theYVar", "theYVar", this),
theZVar("theZVar", "theZVar", this)
{}

MELANCSplinePdf_3D_fast::MELANCSplinePdf_3D_fast(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  RooAbsReal& inYVar,
  RooAbsReal& inZVar,
  const std::vector<T>& inXList,
  const std::vector<T>& inYList,
  const std::vector<T>& inZList,
  const std::vector<std::vector<std::vector<T>>>& inFcnList,
  Bool_t inUseFloor,
  T inFloorEval,
  T inFloorInt
  ) :
  MELANCSplinePdfCore(name, title, inXVar, inXList, inUseFloor, inFloorEval, inFloorInt),
  theYVar("theYVar", "theYVar", this, inYVar),
  theZVar("theZVar", "theZVar", this, inZVar),
  YList(inYList),
  ZList(inZList),
  FcnList(inFcnList)
{
  if (npointsX()>1 && npointsY()>1 && npointsZ()>1){
    // Prepare A and kappa arrays for x, y and z coordinates
    int npoints;
    Double_t det;

    vector<vector<MELANCSplinePdfCore::T>> xA; getKappas(kappaX, 0); getAArray(kappaX, xA);
    npoints=kappaX.size();
    TMatrix_t xAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ xAtrans[i][j]=xA.at(i).at(j); } }
    det=0;
    TMatrix_t xAinv = xAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "MELANCSplinePdf_3D_fast::interpolateFcn: Matrix xA could not be inverted. Something is wrong with the x coordinates of points!" << endl;
      assert(0);
    }

    vector<vector<MELANCSplinePdfCore::T>> yA; getKappas(kappaY, 1); getAArray(kappaY, yA);
    npoints=kappaY.size();
    TMatrix_t yAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ yAtrans[i][j]=yA.at(i).at(j); } }
    det=0;
    TMatrix_t yAinv = yAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "MELANCSplinePdf_3D_fast::interpolateFcn: Matrix yA could not be inverted. Something is wrong with the y coordinates of points!" << endl;
      assert(0);
    }

    vector<vector<MELANCSplinePdfCore::T>> zA; getKappas(kappaZ, 2); getAArray(kappaZ, zA);
    npoints=kappaZ.size();
    TMatrix_t zAtrans(npoints, npoints);
    for (int i=0; i<npoints; i++){ for (int j=0; j<npoints; j++){ zAtrans[i][j]=zA.at(i).at(j); } }
    det=0;
    TMatrix_t zAinv = zAtrans.Invert(&det);
    if (det==0.){
      coutE(InputArguments) << "MELANCSplinePdf_3D_fast::interpolateFcn: Matrix zA could not be inverted. Something is wrong with the z coordinates of points!" << endl;
      assert(0);
    }

    //cout << "MELANCSplinePdf_3D_fast::MELANCSplinePdf_3D_fast: Initial kappa arrays are set up" << endl;

    // Get the grid of coefficients
    int nxpoldim=0;
    int nxbins=0;
    int nybins=0;
    int nypoldim=0;
    std::vector<
      std::vector<std::vector<
      std::vector<std::vector<MELANCSplinePdfCore::T>>
      >>
      > coefsAlongZ; // [ix][A_x,B_x,C_x,D_x][iy][A_x_y,B_x_y,C_x_y,D_x_y][z_k] at each z_k
    for (unsigned int k=0; k<npointsZ(); k++){
      //cout << "Finding coefficients along z line " << k << endl;

      std::vector<std::vector<
        std::vector<std::vector<MELANCSplinePdfCore::T>>
        >> coefficients_perZ; // [ix][A_x,B_x,C_x,D_x][iy][A_x_y,B_x_y,C_x_y,D_x_y] at each z_k

      vector<vector<vector<MELANCSplinePdfCore::T>>> coefsAlongY; // [xbin][Ax(y),Bx(y),Cx(y),Dx(y)][ybin] in each z
      for (unsigned int j=0; j<npointsY(); j++){
        vector<vector<MELANCSplinePdfCore::T>> xcoefsAtYjZk = getCoefficientsPerYPerZ(kappaX, xAinv, j, k, -1); // [ix][Ax,Bx,Cx,Dx] at each y_j z_k
        //cout << "\tCoefficients in y line " << j << " are found" << endl;
        if (j==0){
          if (k==0){
            nxbins=xcoefsAtYjZk.size();
            nxpoldim=xcoefsAtYjZk.at(0).size();
          }
          for (int ix=0; ix<nxbins; ix++){
            vector<vector<MELANCSplinePdfCore::T>> dum_xycoefarray;
            for (int icx=0; icx<nxpoldim; icx++){
              vector<MELANCSplinePdfCore::T> dum_ycoefarray;
              dum_xycoefarray.push_back(dum_ycoefarray);
            }
            coefsAlongY.push_back(dum_xycoefarray);
          }
        }
        if (nxbins!=(int)xcoefsAtYjZk.size() || nxpoldim!=(int)xcoefsAtYjZk.at(0).size()){
          coutE(InputArguments) << "MELANCSplinePdf_3D_fast::interpolateFcn: nxbins!=(int)xcoefsAtYjZk.size() || nxpoldim!=(int)xcoefsAtYjZk.at(0).size()!" << endl;
          assert(0);
        }
        for (int ix=0; ix<nxbins; ix++){
          for (int icx=0; icx<nxpoldim; icx++) coefsAlongY.at(ix).at(icx).push_back(xcoefsAtYjZk.at(ix).at(icx));
        }
      }

      //cout << "\tCoefficients in each y line are found" << endl;

      for (int ix=0; ix<nxbins; ix++){
        // Get the x coefficients interpolated across y
        vector<vector<vector<MELANCSplinePdfCore::T>>> xCoefs;
        for (int icx=0; icx<nxpoldim; icx++){
          vector<vector<MELANCSplinePdfCore::T>> yCoefs = getCoefficientsAlongDirection(kappaY, yAinv, coefsAlongY.at(ix).at(icx), -1); // [iy][A,B,C,D]
          xCoefs.push_back(yCoefs);
        }
        coefficients_perZ.push_back(xCoefs);
      }

      if (k==0){
        nybins = coefficients_perZ.at(0).at(0).size();
        nypoldim = coefficients_perZ.at(0).at(0).at(0).size();
        for (int ix=0; ix<nxbins; ix++){
          vector<vector<vector<vector<MELANCSplinePdfCore::T>>>> xCoefs;
          for (int icx=0; icx<nxpoldim; icx++){
            vector<vector<vector<MELANCSplinePdfCore::T>>> xCoefsAlongY;
            for (int iy=0; iy<nybins; iy++){
              vector<vector<MELANCSplinePdfCore::T>> yCoefs;
              for (int icy=0; icy<nypoldim; icy++){
                vector<MELANCSplinePdfCore::T> yCoefAtZj;
                yCoefs.push_back(yCoefAtZj);
              }
              xCoefsAlongY.push_back(yCoefs);
            }
            xCoefs.push_back(xCoefsAlongY);
          }
          coefsAlongZ.push_back(xCoefs);
        }
      }
      for (int ix=0; ix<nxbins; ix++){
        for (int icx=0; icx<nxpoldim; icx++){
          for (int iy=0; iy<nybins; iy++){
            for (int icy=0; icy<nypoldim; icy++){
              coefsAlongZ.at(ix).at(icx).at(iy).at(icy).push_back(coefficients_perZ.at(ix).at(icx).at(iy).at(icy));
            }
          }
        }
      }
    }

    //cout << "Finding the final coefficients" << endl;
    for (int ix=0; ix<nxbins; ix++){
      vector<vector<vector<vector<vector<MELANCSplinePdfCore::T>>>>> xCoefs;
      for (int icx=0; icx<nxpoldim; icx++){
        vector<vector<vector<vector<MELANCSplinePdfCore::T>>>> xCoefsAlongY;
        for (int iy=0; iy<nybins; iy++){
          vector<vector<vector<MELANCSplinePdfCore::T>>> yCoefs;
          for (int icy=0; icy<nypoldim; icy++){
            vector<vector<MELANCSplinePdfCore::T>> yCoefsAlongZ = getCoefficientsAlongDirection(kappaZ, zAinv, coefsAlongZ.at(ix).at(icx).at(iy).at(icy), -1); // [iz][A,B,C,D]
            yCoefs.push_back(yCoefsAlongZ);
          }
          xCoefsAlongY.push_back(yCoefs);
        }
        xCoefs.push_back(xCoefsAlongY);
      }
      coefficients.push_back(xCoefs);
    }

  }
  else assert(0);

  emptyFcnList();
}

MELANCSplinePdf_3D_fast::MELANCSplinePdf_3D_fast(
  const MELANCSplinePdf_3D_fast& other,
  const char* name
  ) :
  MELANCSplinePdfCore(other, name),
  theYVar("theYVar", this, other.theYVar),
  theZVar("theZVar", this, other.theZVar),
  YList(other.YList),
  ZList(other.ZList),
  FcnList(other.FcnList),
  kappaX(other.kappaX),
  kappaY(other.kappaY),
  kappaZ(other.kappaZ),
  coefficients(other.coefficients)
{}

MELANCSplinePdfCore::T MELANCSplinePdf_3D_fast::interpolateFcn(Int_t code, const char* rangeName)const{
  DefaultAccumulator<MELANCSplinePdfCore::T> res;

  if (verbosity==MELANCSplinePdfCore::kVerbose) cout << "MELANCSplinePdf_3D_fast(" << GetName() << ")::interpolateFcn begin with code: " << code << endl;

  // Get bins
  vector<Int_t> varprime; varprime.push_back(2); varprime.push_back(3); varprime.push_back(5);
  const Int_t ndims=(const Int_t)varprime.size();
  vector<Int_t> varbin;
  vector<Int_t> varbinmin;
  vector<Int_t> varbinmax;
  vector<Float_t> tvar;
  vector<Float_t> tvarmin;
  vector<Float_t> tvarmax;
  vector<const vector<MELANCSplinePdfCore::T>*> varkappa; varkappa.push_back(&kappaX); varkappa.push_back(&kappaY); varkappa.push_back(&kappaZ);
  vector<const RooRealProxy*> varcoord; varcoord.push_back(&theXVar); varcoord.push_back(&theYVar); varcoord.push_back(&theZVar);
  for (Int_t idim=0; idim<ndims; idim++){
    Int_t binval=-1;
    Int_t binmin=-1;
    Int_t binmax=-1;
    Float_t tval=0;
    Float_t tmin=0;
    Float_t tmax=0;
    if (code==0 || code%varprime.at(idim)!=0){
      binval = getWhichBin(*(varcoord.at(idim)), idim);
      tval = getTVar(*(varkappa.at(idim)), *(varcoord.at(idim)), binval, idim);
    }
    else{
      binmin = getWhichBin(varcoord.at(idim)->min(rangeName), idim);
      tmin = getTVar(*(varkappa.at(idim)), varcoord.at(idim)->min(rangeName), binmin, idim);
      binmax = getWhichBin(varcoord.at(idim)->max(rangeName), idim);
      tmax = getTVar(*(varkappa.at(idim)), varcoord.at(idim)->max(rangeName), binmax, idim);
    }
    varbin.push_back(binval);
    varbinmin.push_back(binmin);
    varbinmax.push_back(binmax);
    tvar.push_back(tval);
    tvarmin.push_back(tmin);
    tvarmax.push_back(tmax);
  }

  for (int ix=0; ix<(int)coefficients.size(); ix++){
    if (
      (varbin[0]>=0 && ix!=varbin[0])
      ||
      (varbinmin[0]>=0 && varbinmax[0]>=varbinmin[0] && !(varbinmin[0]<=ix && ix<=varbinmax[0]))
      ) continue;
    MELANCSplinePdfCore::T txlow=0, txhigh=1;
    if (code>0 && code%varprime[0]==0){
      if (ix==varbinmin[0]) txlow=tvarmin[0];
      if (ix==varbinmax[0]) txhigh=tvarmax[0];
    }
    else txhigh=tvar[0];
    // Get the x coefficients interpolated across y
    vector<MELANCSplinePdfCore::T> xCoefs;
    for (int icx=0; icx<(int)coefficients.at(ix).size(); icx++){
      DefaultAccumulator<MELANCSplinePdfCore::T> theXCoef;
      for (int iy=0; iy<(int)coefficients.at(ix).at(icx).size(); iy++){
        if (
          (varbin[1]>=0 && iy!=varbin[1])
          ||
          (varbinmin[1]>=0 && varbinmax[1]>=varbinmin[1] && !(varbinmin[1]<=iy && iy<=varbinmax[1]))
          ) continue;
        MELANCSplinePdfCore::T tylow=0, tyhigh=1;
        if (code>0 && code%varprime[1]==0){
          if (iy==varbinmin[1]) tylow=tvarmin[1];
          if (iy==varbinmax[1]) tyhigh=tvarmax[1];
        }
        else tyhigh=tvar[1];
        // Get the y coefficients interpolated across z
        vector<MELANCSplinePdfCore::T> yCoefs;
        for (int icy=0; icy<(int)coefficients.at(ix).at(icx).at(iy).size(); icy++){
          DefaultAccumulator<MELANCSplinePdfCore::T> theYCoef;
          for (int iz=0; iz<(int)coefficients.at(ix).at(icx).at(iy).at(icy).size(); iz++){
            if (
              (varbin[2]>=0 && iz!=varbin[2])
              ||
              (varbinmin[2]>=0 && varbinmax[2]>=varbinmin[2] && !(varbinmin[2]<=iz && iz<=varbinmax[2]))
              ) continue;

            MELANCSplinePdfCore::T tzlow=0, tzhigh=1;
            if (code>0 && code%varprime[2]==0){
              if (iz==varbinmin[2]) tzlow=tvarmin[2];
              if (iz==varbinmax[2]) tzhigh=tvarmax[2];
            }
            else tzhigh=tvar[2];

            theYCoef += evalSplineSegment(coefficients.at(ix).at(icx).at(iy).at(icy).at(iz), varkappa[2]->at(iz), tzhigh, tzlow, (code>0 && code%varprime[2]==0));
          }
          yCoefs.push_back(theYCoef);
        }
        theXCoef += evalSplineSegment(yCoefs, varkappa[1]->at(iy), tyhigh, tylow, (code>0 && code%varprime[1]==0));
      }
      xCoefs.push_back(theXCoef);
    }
    // Evaluate value of spline at x with coefficients evaluated at y
    res += evalSplineSegment(xCoefs, varkappa[0]->at(ix), txhigh, txlow, (code>0 && code%varprime[0]==0));
  }

  return res;
}

void MELANCSplinePdf_3D_fast::getKappas(vector<MELANCSplinePdfCore::T>& kappas, const Int_t whichDirection){
  kappas.clear();
  MELANCSplinePdfCore::T kappa=1;

  Int_t npoints;
  vector<MELANCSplinePdfCore::T> const* coord;
  if (whichDirection==0){
    npoints=npointsX();
    coord=&XList;
  }
  else if (whichDirection==1){
    npoints=npointsY();
    coord=&YList;
  }
  else{
    npoints=npointsZ();
    coord=&ZList;
  }

  for (Int_t j=0; j<npoints-1; j++){
    MELANCSplinePdfCore::T val_j = coord->at(j);
    MELANCSplinePdfCore::T val_jpo = coord->at(j+1);
    kappa = 1./(val_jpo-val_j);
    kappas.push_back(kappa);
  }
  kappas.push_back(kappa); // Push the same kappa_(N-1)=kappa_(N-2) at the end point
}
Int_t MELANCSplinePdf_3D_fast::getWhichBin(const MELANCSplinePdfCore::T& val, const Int_t whichDirection)const{
  Int_t bin=-1;
  MELANCSplinePdfCore::T valj, valjpo;
  Int_t npoints;
  vector<MELANCSplinePdfCore::T> const* coord;
  if (whichDirection==0){
    npoints=npointsX();
    coord=&XList;
  }
  else if (whichDirection==1){
    npoints=npointsY();
    coord=&YList;
  }
  else{
    npoints=npointsZ();
    coord=&ZList;
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
MELANCSplinePdfCore::T MELANCSplinePdf_3D_fast::getTVar(const vector<MELANCSplinePdfCore::T>& kappas, const MELANCSplinePdfCore::T& val, const Int_t& bin, const Int_t whichDirection)const{
  const MELANCSplinePdfCore::T& K=kappas.at(bin);
  vector<MELANCSplinePdfCore::T> const* coord;
  if (whichDirection==0) coord=&XList;
  else if (whichDirection==1) coord=&YList;
  else coord=&ZList;
  return (val-coord->at(bin))*K;
}

vector<vector<MELANCSplinePdfCore::T>> MELANCSplinePdf_3D_fast::getCoefficientsPerYPerZ(
  const std::vector<MELANCSplinePdfCore::T>& kappaX, const TMatrix_t& xAinv,
  const Int_t& ybin, const Int_t& zbin,
  const Int_t xbin
  )const{
  vector<MELANCSplinePdfCore::T> fcnList;
  for (unsigned int bin=0; bin<npointsX(); bin++){ fcnList.push_back(FcnList.at(zbin).at(ybin).at(bin)); }
  vector<vector<MELANCSplinePdfCore::T>> coefs = getCoefficientsAlongDirection(kappaX, xAinv, fcnList, xbin);
  return coefs;
}

Double_t MELANCSplinePdf_3D_fast::evaluate() const{
  Double_t value = interpolateFcn(0);
  if (useFloor && value<floorEval){
    if (verbosity>=MELANCSplinePdfCore::kError) coutE(Eval) << "MELANCSplinePdf_3D_fast ERROR::MELANCSplinePdf_3D_fast(" << GetName() << ") evaluation returned " << value << " at (x, y, z) = (" << theXVar << ", " << theYVar << ", " << theZVar << ")" << endl;
    value = floorEval;
  }
  if (verbosity==MELANCSplinePdfCore::kVerbose){ cout << "MELANCSplinePdf_3D_fast(" << GetName() << ")::evaluate = " << value << " at (x, y, z) = (" << theXVar << ", " << theYVar << ", " << theZVar << ")" << endl; }
  return value;
}
Int_t MELANCSplinePdf_3D_fast::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  if (_forceNumInt) return 0;

  Int_t code=1;

  RooArgSet Xdeps, Ydeps, Zdeps;
  RooRealVar* rrv_x = dynamic_cast<RooRealVar*>(theXVar.absArg());
  RooRealVar* rrv_y = dynamic_cast<RooRealVar*>(theYVar.absArg());
  RooRealVar* rrv_z = dynamic_cast<RooRealVar*>(theZVar.absArg());
  if (rrv_x==0) theXVar.absArg()->leafNodeServerList(&Xdeps, 0, true);
  if (rrv_y==0) theYVar.absArg()->leafNodeServerList(&Ydeps, 0, true);
  if (rrv_z==0) theZVar.absArg()->leafNodeServerList(&Zdeps, 0, true);

  if (rrv_x!=0){
    if (
      (Ydeps.find(*rrv_x)==0 || rrv_y!=0)
      &&
      (Zdeps.find(*rrv_x)==0 || rrv_z!=0)
      ){
      if (matchArgs(allVars, analVars, theXVar)) code*=2;
    }
  }
  if (rrv_y!=0){
    if (
      (Xdeps.find(*rrv_y)==0 || rrv_x!=0)
      &&
      (Zdeps.find(*rrv_y)==0 || rrv_z!=0)
      ){
      if (matchArgs(allVars, analVars, theYVar)) code*=3;
    }
  }
  if (rrv_z!=0){
    if (
      (Xdeps.find(*rrv_z)==0 || rrv_x!=0)
      &&
      (Ydeps.find(*rrv_z)==0 || rrv_y!=0)
      ){
      if (matchArgs(allVars, analVars, theZVar)) code*=5;
    }
  }

  if (code==1) code=0;
  return code;
}
Double_t MELANCSplinePdf_3D_fast::analyticalIntegral(Int_t code, const char* rangeName) const{
  Double_t value = interpolateFcn(code, rangeName);
  if (useFloor && value<floorInt){
    if (verbosity>=MELANCSplinePdfCore::kError) coutE(Integration) << "MELANCSplinePdf_3D_fast ERROR::MELANCSplinePdf_3D_fast(" << GetName() << ") integration returned " << value << " for code = " << code << endl;
    value = floorInt;
  }
  if (verbosity==MELANCSplinePdfCore::kVerbose){ cout << "MELANCSplinePdf_3D_fast(" << GetName() << ")::analyticalIntegral = " << value << " for code = " << code << endl; }
  return value;
}
