#ifndef MELANCSPLINEPDFCORE
#define MELANCSPLINEPDFCORE  

#include <vector>
#include "MELAAccumulators.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooConstVar.h"
#include "RooArgList.h"

class MELANCSplinePdfCore : public RooAbsPdf{
public:
  typedef Float_t T;
  typedef TMatrixT<T> TMatrix_t;
  typedef TVectorT<T> TVector_t;

  enum VerbosityLevel{
    kSilent,
    kVerbose
  };

  MELANCSplinePdfCore();
  MELANCSplinePdfCore(
    const char* name,
    const char* title
    );
  MELANCSplinePdfCore(
    const char* name,
    const char* title,
    RooAbsReal& inXVar,
    std::vector<T>& inXList
    );
  MELANCSplinePdfCore(const MELANCSplinePdfCore& other, const char* name=0);
  virtual TObject* clone(const char* newname)const=0;
  inline virtual ~MELANCSplinePdfCore(){}

  virtual void setVerbosity(VerbosityLevel flag);

protected:
  VerbosityLevel verbosity;

  RooRealProxy theXVar;
  std::vector<T> XList;
  
  virtual void emptyFcnList() = 0;

  unsigned int npointsX()const{ return XList.size(); }

  virtual Int_t getWhichBin(const T& val, const Int_t whichDirection)const = 0;
  virtual void getKappas(std::vector<T>& kappas, const Int_t whichDirection)const = 0;
  virtual T getTVar(const std::vector<T>& kappas, const T& val, const Int_t& bin, const Int_t whichDirection)const = 0;

  virtual T interpolateFcn(Int_t code, const char* rangeName=0)const = 0;
  virtual Double_t evaluate()const = 0;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0)const = 0;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const = 0;

  virtual void getBArray(const std::vector<T>& kappas, const std::vector<T>& fcnList, std::vector<T>& BArray)const;
  virtual void getAArray(const std::vector<T>& kappas, std::vector<std::vector<T>>& AArray)const;
  virtual std::vector<std::vector<T>> getCoefficientsAlongDirection(const std::vector<T>& kappas, const TMatrix_t& Ainv, const std::vector<T>& fcnList, const Int_t pickBin)const;
  virtual std::vector<T> getCoefficients(const TVector_t& S, const std::vector<T>& kappas, const std::vector<T>& fcnList, const Int_t& bin)const;

  virtual T evalSplineSegment(const std::vector<T>& coefs, const T kappa, T tup, T tdn=0, Bool_t doIntegrate=false)const;

};
 
#endif
