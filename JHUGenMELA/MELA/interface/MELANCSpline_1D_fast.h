#ifndef MELANCSPLINE_1D_FAST
#define MELANCSPLINE_1D_FAST

#include <vector>
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "MELANCSplineCore.h"


class MELANCSpline_1D_fast : public MELANCSplineCore{
protected:
  std::vector<T> FcnList; // List of function values

  std::vector<T> kappaX;
  std::vector<std::vector<T>> coefficients;

public:
  MELANCSpline_1D_fast();
  MELANCSpline_1D_fast(
    const char* name,
    const char* title
    );
  MELANCSpline_1D_fast(
    const char* name,
    const char* title,
    RooAbsReal& inXVar,
    const std::vector<T>& inXList,
    const std::vector<T>& inFcnList,
    Bool_t inUseFloor=true,
    T inFloorEval=0,
    T inFloorInt=0
    );
  MELANCSpline_1D_fast(const MELANCSpline_1D_fast& other, const char* name=0);
	virtual TObject* clone(const char* newname)const { return new MELANCSpline_1D_fast(*this, newname); }
	inline virtual ~MELANCSpline_1D_fast(){}

  void setRangeValidity(const T valmin, const T valmax, const Int_t whichDirection=0);

  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0)const;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const;

protected:
  virtual void emptyFcnList(){ std::vector<T> tmp; FcnList.swap(tmp); }

  virtual Int_t getWhichBin(const T& val, const Int_t whichDirection)const;
  virtual T getTVar(const std::vector<T>& kappas, const T& val, const Int_t& bin, const Int_t whichDirection)const;
  virtual void getKappas(std::vector<T>& kappas, const Int_t whichDirection);

  Bool_t testRangeValidity(const T& val, const Int_t whichDirection=0)const;
  void cropValueForRange(T& val, const Int_t whichDirection=0)const;

  virtual T interpolateFcn(Int_t code, const char* rangeName=0)const;

  virtual Double_t evaluate()const;


  ClassDef(MELANCSpline_1D_fast, 1)

};
 
#endif
