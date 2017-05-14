#include "MELALinearInterpPdf.h" 
#include <cmath>
#include "Riostream.h" 
#include "TMath.h"

using namespace TMath;
using namespace RooFit;
using namespace std;
using namespace TNumericUtil;


MELALinearInterpPdf::MELALinearInterpPdf() :
RooAbsPdf(),
verbosity(MELALinearInterpPdf::kSilent),
useFloor(true), floorEval(0), floorInt(0),
rangeXmin(1), rangeXmax(-1),
theXVar("theXVar", "theXVar", this),
FcnList("FcnList", "FcnList", this),
leafDepsList("leafDepsList", "leafDepsList", this)
{}

MELALinearInterpPdf::MELALinearInterpPdf(
  const char* name,
  const char* title
  ) :
  RooAbsPdf(name, title),
  verbosity(MELALinearInterpPdf::kSilent),
  useFloor(true), floorEval(0), floorInt(0),
  rangeXmin(1), rangeXmax(-1),
  theXVar("theXVar", "theXVar", this),
  FcnList("FcnList", "FcnList", this),
  leafDepsList("leafDepsList", "leafDepsList", this)
{}

MELALinearInterpPdf::MELALinearInterpPdf(
  const char* name,
  const char* title,
  RooAbsReal& inXVar,
  const std::vector<MELALinearInterpPdf::T>& inXList,
  const RooArgList& inFcnList,
  Bool_t inUseFloor,
  T inFloorEval,
  T inFloorInt
  ) :
  RooAbsPdf(name, title),
  verbosity(MELALinearInterpPdf::kSilent),
  useFloor(inUseFloor), floorEval(inFloorEval), floorInt(inFloorInt),
  rangeXmin(1), rangeXmax(-1),
  XList(inXList),
  theXVar("theXVar", "theXVar", this, inXVar),
  FcnList("FcnList", "FcnList", this),
  leafDepsList("leafDepsList", "leafDepsList", this)
{
  RooArgSet leafset;
  if (!dynamic_cast<RooRealVar*>(theXVar.absArg())){
    RooArgSet tmpdeps;
    theXVar.absArg()->leafNodeServerList(&tmpdeps, 0, true);
    leafset.add(tmpdeps);
  }

  TIterator* iter = inFcnList.createIterator();
  RooAbsArg* absarg;
  while ((absarg = (RooAbsArg*)iter->Next())){
    if (!dynamic_cast<RooAbsReal*>(absarg)){
      coutE(InputArguments) << "ERROR::MELALinearInterpPdf(" << GetName() << ") function " << absarg->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    FcnList.add(*absarg);

    if (!dynamic_cast<RooRealVar*>(absarg)){
      RooArgSet tmpdeps;
      absarg->leafNodeServerList(&tmpdeps, 0, true);
      leafset.add(tmpdeps);
    }
  }
  delete iter;

  iter = leafset.createIterator();
  while ((absarg = (RooAbsArg*)iter->Next())){ if (dynamic_cast<RooRealVar*>(absarg)) leafDepsList.add(*absarg); }
  delete iter;

  if (FcnList.getSize()!=(int)XList.size()){
    coutE(InputArguments) << "MELALinearInterpPdf ERROR::MELALinearInterpPdf(" << GetName() << ") input XList size not the same as FcnList size!" << endl;
    assert(0);
  }
}

MELALinearInterpPdf::MELALinearInterpPdf(
  const MELALinearInterpPdf& other,
  const char* name
  ) :
  RooAbsPdf(other, name),
  verbosity(other.verbosity),
  useFloor(other.useFloor), floorEval(other.floorEval), floorInt(other.floorInt),
  rangeXmin(other.rangeXmin), rangeXmax(other.rangeXmax),
  XList(other.XList),
  theXVar("theXVar", this, other.theXVar),
  FcnList("FcnList", this, other.FcnList),
  leafDepsList("leafDepsList", this, other.leafDepsList)
{}

void MELALinearInterpPdf::setVerbosity(VerbosityLevel flag){ verbosity = flag; }
void MELALinearInterpPdf::setEvalFloor(MELALinearInterpPdf::T val){ floorEval = val; }
void MELALinearInterpPdf::setIntFloor(MELALinearInterpPdf::T val){ floorInt = val; }
void MELALinearInterpPdf::doFloor(Bool_t flag){ useFloor = flag; }

Int_t MELALinearInterpPdf::getWhichBin(const MELALinearInterpPdf::T& val)const{
  Int_t bin=-1;
  MELALinearInterpPdf::T valj, valjpo;
  Int_t np;
  vector<MELALinearInterpPdf::T> const* coord=&XList;
  np=npoints();

  if (np<=1) bin=0;
  else{
    valjpo = coord->at(0);
    for (Int_t j=0; j<np-1; j++){
      valj = coord->at(j);
      valjpo = coord->at(j+1);
      if (val<valjpo && val>=valj){ bin=j; break; }
    }
    if (bin==-1 && val>=valjpo) bin=np-2;
    else if (bin==-1) bin=0;
  }

  return bin;
}
MELALinearInterpPdf::T MELALinearInterpPdf::getKappa(const Int_t& bin)const{
  if (bin>=0){
    vector<MELALinearInterpPdf::T> const* coord=&XList;
    MELALinearInterpPdf::T diff = 1;
    if (Int_t(coord->size())>(bin+1)) diff = coord->at(bin+1)-coord->at(bin);
    if (fabs(diff)>MELALinearInterpPdf::T(0)) return MELALinearInterpPdf::T(1)/diff;
    else return 0;
  }
  else return 1;
}
MELALinearInterpPdf::T MELALinearInterpPdf::getTVar(const MELALinearInterpPdf::T& val)const{
  const Int_t bin = getWhichBin(val);
  vector<MELALinearInterpPdf::T> const* coord=&XList;
  const MELALinearInterpPdf::T kappa = getKappa(bin);
  return (val - coord->at(bin))*kappa;
}

MELALinearInterpPdf::T MELALinearInterpPdf::interpolateFcn(Int_t code, const char* rangeName)const{
  const Int_t xprime = 101;
  const bool intAlongX=(code>0 && code%xprime==0);

  DefaultAccumulator<MELALinearInterpPdf::T> res;
  Int_t coderem = (intAlongX ? code/xprime : code);
  if (coderem==1) coderem=0;

  if (verbosity==MELALinearInterpPdf::kVerbose){
    if (intAlongX) cout << "MELALinearInterpPdf(" << GetName() << ")::interpolateFcn integrates using code = " << code << ", coderem = " << coderem << endl;
    else cout << "MELALinearInterpPdf(" << GetName() << ")::interpolateFcn evaluates using code = " << code << ", coderem = " << coderem << endl;
  }

  // Get bins
  Int_t xbin=-1, xbinmin=-1, xbinmax=-1;
  MELALinearInterpPdf::T tx=0, txmin=0, txmax=0;
  if (!intAlongX){ // Case to just compute the value at x
    if (!testRangeValidity(theXVar)) return 0;
    xbin = getWhichBin(theXVar);
    tx = getTVar(theXVar);
  }
  else{ // Case to integrate along x
    MELALinearInterpPdf::T coordmin = theXVar.min(rangeName); cropValueForRange(coordmin);
    MELALinearInterpPdf::T coordmax = theXVar.max(rangeName); cropValueForRange(coordmax);
    xbinmin = getWhichBin(coordmin);
    txmin = getTVar(coordmin);
    xbinmax = getWhichBin(coordmax);
    txmax = getTVar(coordmax);
  }

  int nxbins = npoints()-1;
  if (nxbins==0){
    MELALinearInterpPdf::T fcnval;
    if (coderem==0) fcnval = dynamic_cast<const RooAbsReal*>(FcnList.at(0))->getVal();
    else fcnval = dynamic_cast<const RooAbsReal*>(FcnList.at(0))->analyticalIntegral(coderem, rangeName);
    res = fcnval;
    if (intAlongX) res *= (txmax-txmin);
  }
  else{
    for (int ix=0; ix<nxbins; ix++){
      if (
        (xbin>=0 && ix!=xbin)
        ||
        (xbinmin>=0 && xbinmax>=xbinmin && !(xbinmin<=ix && ix<=xbinmax))
        ) continue;

      MELALinearInterpPdf::T txlow=0, txhigh=1;
      if (intAlongX){
        if (ix==xbinmin) txlow=txmin;
        if (ix==xbinmax) txhigh=txmax;
      }
      else txhigh=tx;

      MELALinearInterpPdf::T fcnval[2]={ 0 };
      if (coderem==0){
        for (unsigned int j=0; j<2; j++) fcnval[j] = dynamic_cast<const RooAbsReal*>(FcnList.at(ix+j))->getVal();
      }
      else{
        for (unsigned int j=0; j<2; j++) fcnval[j] = dynamic_cast<const RooAbsReal*>(FcnList.at(ix+j))->analyticalIntegral(coderem, rangeName);
      }

      DefaultAccumulator<MELALinearInterpPdf::T> segval = fcnval[0];
      if (intAlongX) segval += (fcnval[1] - fcnval[0])*(pow(txhigh, 2)-pow(txlow, 2))/2.;
      else segval += (fcnval[1] - fcnval[0])*txhigh;
      res += segval;

      if (verbosity==MELALinearInterpPdf::kVerbose) cout
        << "MELALinearInterpPdf(" << GetName() << ")::interpolateFcn evaluated bin " << ix
        << " with txlow = " << txlow << ", txhigh = " << txhigh << ", fcnval[0] = " << fcnval[0] << ", fcnval[1] = " << fcnval[1]
        << endl;
    }
  }

  return res;
}
Double_t MELALinearInterpPdf::evaluate()const{
  Double_t value = interpolateFcn(0);
  if (useFloor && value<floorEval){
    if (verbosity>=MELALinearInterpPdf::kError) coutE(Eval) << "MELALinearInterpPdf ERROR::MELALinearInterpPdf(" << GetName() << ") evaluation returned " << value << " at x = " << theXVar << endl;
    value = floorEval;
  }
  if (verbosity==MELALinearInterpPdf::kVerbose){ cout << "MELALinearInterpPdf(" << GetName() << ")::evaluate = " << value << " at x = " << theXVar << endl; }
  return value;
}
Int_t MELALinearInterpPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName)const{
  if (_forceNumInt) return 0;

  const Int_t xprime = 101;
  Int_t code=1;

  RooArgSet Xdeps, Fcndeps;
  RooRealVar* rrv_x = dynamic_cast<RooRealVar*>(theXVar.absArg());
  if (rrv_x==0) theXVar.absArg()->leafNodeServerList(&Xdeps, 0, true);

  for (int ix=0; ix<FcnList.getSize(); ix++){
    RooArgSet tmpdeps;
    RooAbsReal const* fcn = dynamic_cast<const RooAbsReal*>(FcnList.at(ix));
    fcn->leafNodeServerList(&tmpdeps, 0, true);
    Fcndeps.add(tmpdeps);
    if (ix==0){
      RooArgSet mirrorvars;
      TIterator* iter = allVars.createIterator();
      RooAbsArg* var;
      while ((var = (RooAbsArg*)iter->Next())){
        RooRealVar* rrv_tmp = dynamic_cast<RooRealVar*>(var);
        if (rrv_tmp==0) continue;
        bool isFound=false;
        if (rrv_x!=0){
          isFound = (rrv_tmp==rrv_x);
        }
        else{
          TIterator* iter_deps = Xdeps.createIterator();
          RooAbsArg* var_deps;
          while ((var_deps = (RooAbsArg*)iter_deps->Next())){
            RooRealVar* depvar = dynamic_cast<RooRealVar*>(var_deps);
            if (depvar==0) continue;
            if (depvar==rrv_tmp){ isFound=true; break; }
          }
          delete iter_deps;
        }
        if (!isFound) mirrorvars.add(*var);
      }
      delete iter;

      Int_t codetmp = fcn->getAnalyticalIntegral(mirrorvars, analVars, rangeName);
      if (codetmp!=0) code *= codetmp;
    }
  }

  if (rrv_x!=0){
    if (Fcndeps.find(*rrv_x)==0){
      if (matchArgs(allVars, analVars, theXVar)) code *= xprime;
    }
  }

  if (code==1) code=0;
  if (verbosity==MELALinearInterpPdf::kVerbose){
    cout << "MELALinearInterpPdf(" << GetName() << ")::getAnalyticalIntegral code = " << code << endl;
    allVars.Print("v");
    analVars.Print("v");
  }
  return code;
}
Double_t MELALinearInterpPdf::analyticalIntegral(Int_t code, const char* rangeName)const{
  Double_t value = interpolateFcn(code, rangeName);
  if (useFloor && value<floorInt){
    if (verbosity>=MELALinearInterpPdf::kError) coutE(Integration) << "MELALinearInterpPdf ERROR::MELALinearInterpPdf(" << GetName() << ") integration returned " << value << " for code = " << code << endl;
    value = floorInt;
  }
  if (verbosity==MELALinearInterpPdf::kVerbose){ cout << "MELALinearInterpPdf(" << GetName() << ")::analyticalIntegral = " << value << " for code = " << code << endl; }
  return value;
}

Bool_t MELALinearInterpPdf::testRangeValidity(const T& val) const{
  const T* range[2];
  range[0] = &rangeXmin;
  range[1] = &rangeXmax;
  return (*(range[0])>*(range[1]) || (val>=*(range[0]) && val<=*(range[1])));
}
void MELALinearInterpPdf::setRangeValidity(const T valmin, const T valmax){
  T* range[2];
  range[0] = &rangeXmin;
  range[1] = &rangeXmax;
  *(range[0])=valmin;
  *(range[1])=valmax;
}
void MELALinearInterpPdf::cropValueForRange(T& val)const{
  if (testRangeValidity(val)) return;
  const T* range[2];
  range[0] = &rangeXmin;
  range[1] = &rangeXmax;
  if (val<*(range[0])) val = *(range[0]);
  if (val>*(range[1])) val = *(range[1]);
}
