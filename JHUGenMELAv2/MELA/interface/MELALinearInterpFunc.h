#ifndef MELALINEARINTERPFUNC
#define MELALINEARINTERPFUNC

#include <vector>
#include "MELAAccumulators.h"
#include "RooRealVar.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "RooArgList.h"
#include "RooMsgService.h"


class MELALinearInterpFunc : public RooAbsReal{
public:
  typedef Float_t T;

  enum VerbosityLevel{
    kSilent,
    kError,
    kVerbose
  };

public:
  MELALinearInterpFunc();
  MELALinearInterpFunc(
    const char* name,
    const char* title
    );
  MELALinearInterpFunc(
    const char* name,
    const char* title,
    RooAbsReal& inXVar,
    const std::vector<T>& inXList,
    const RooArgList& inFcnList,
    Bool_t inUseFloor=true,
    T inFloorEval=1e-15,
    T inFloorInt=1e-10
    );
  MELALinearInterpFunc(const MELALinearInterpFunc& other, const char* name=0);
  virtual TObject* clone(const char* newname)const { return new MELALinearInterpFunc(*this, newname); }
  inline virtual ~MELALinearInterpFunc(){}

  void setVerbosity(VerbosityLevel flag);
  void setEvalFloor(T val);
  void setIntFloor(T val);
  void doFloor(Bool_t flag);

  void setRangeValidity(const T valmin, const T valmax);

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0)const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

protected:
  VerbosityLevel verbosity;
  Bool_t useFloor;
  T floorEval;
  T floorInt;

  T rangeXmin;
  T rangeXmax;

  std::vector<T> XList;
  RooRealProxy theXVar;
  RooListProxy FcnList; // List of function values
  RooListProxy leafDepsList; // List of all leaf dependents

  unsigned int npoints()const{ return XList.size(); }

  Int_t getWhichBin(const T& val)const;
  T getKappa(const Int_t& bin)const;
  T getTVar(const T& val)const;

  Bool_t testRangeValidity(const T& val)const;
  void cropValueForRange(T& val)const;

  T interpolateFcn(Int_t code, const char* rangeName=0)const;
  Double_t evaluate()const;


  ClassDef(MELALinearInterpFunc, 0)

};
 
#endif
