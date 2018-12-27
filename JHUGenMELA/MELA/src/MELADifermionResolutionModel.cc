#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include "TGraph.h"
#include "MELADifermionResolutionModel.h"


MELADifermionResolutionModel::MELADifermionResolutionModel(TVar::Production prod_, float sqrts_, TString filename, TString extraName) :
valid(false),
prod(prod_),
sqrts(sqrts_),
varReco(0),
splineFactory(0),
recoBW(0)
{
  TFile* finput = TFile::Open(filename, "read");
  if (finput!=0){
    TGraph* tg = (TGraph*)finput->Get("tg_recoBW");
    if (tg!=0){
      TString appendName = extraName + TVar::ProductionName(prod);
      varReco = new RooRealVar(Form("%s_varReco", appendName.Data()), "", 0.1, sqrts*1000.);
      splineFactory = new MELANCSplineFactory_1D(*varReco, appendName);
      splineFactory->setPoints(tg);
      recoBW = splineFactory->getPDF();
      valid=(recoBW!=0);
    }
    finput->Close();
  }
}
MELADifermionResolutionModel::~MELADifermionResolutionModel(){
  // Do not delete recoBW. It is owned by splineFactory.
  delete splineFactory;
  delete varReco;
}
float MELADifermionResolutionModel::getVal(float val){
  float result = -1;
  if (recoBW!=0){
    if (val<0.){ /*result=-1;*/ }
    else if (val<0.1) result=0;
    else{
      varReco->setVal(val);
      result = recoBW->getVal();
      result /= 2.*val; // Result is BW(s), not BW(m)
    }
  }
  return result;
}
bool MELADifermionResolutionModel::isProduction(TVar::Production prod_){ return (prod_==prod); }

