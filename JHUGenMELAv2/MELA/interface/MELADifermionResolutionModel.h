#ifndef MELADIFERMIONRESOLUTIONMODEL
#define MELADIFERMIONRESOLUTIONMODEL

#include "TUtil.hh"
#include "MELANCSplineFactory_1D.h"


class MELADifermionResolutionModel{
protected:
  bool valid;
  TVar::Production prod;
  float sqrts;
  RooRealVar* varReco;
  MELANCSplineFactory_1D* splineFactory;
  MELAFuncPdf* recoBW;

public:
  MELADifermionResolutionModel(TVar::Production prod_, float sqrts_, TString filename, TString extraName);
  ~MELADifermionResolutionModel();
  float getVal(float val);
  bool isProduction(TVar::Production prod_);
  bool isValid(){ return valid; }

};


#endif
