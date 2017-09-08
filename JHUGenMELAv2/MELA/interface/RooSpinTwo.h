#ifndef ROOSPINTWO
#define ROOSPINTWO

#include "RooSpin.h"


class RooSpinTwo : public RooSpin {
public:

  struct modelCouplings{
    RooAbsReal* bList[SIZE_GVV][2];
    RooRealVar* Lambda;
    RooRealVar* f_spinz1; // Set to 1 for qqb production
    RooRealVar* f_spinz2; // Set to 0 for qqb production
    // There is no equivalent to graviton_qq_left/right yet!
  };

  RooSpinTwo(){};
  RooSpinTwo(
    const char* name, const char* title,
    modelMeasurables _measurables,
    modelParameters _parameters,
    modelCouplings _couplings,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll
    );

  RooSpinTwo(const RooSpinTwo& other, const char* name=0);
  virtual TObject* clone(const char* newname) const = 0;
  inline virtual ~RooSpinTwo(){}

  virtual Double_t evaluate() const = 0;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const = 0;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const = 0;

protected:

  RooRealProxy  b1Val;
  RooRealProxy  b2Val;
  RooRealProxy  b3Val;
  RooRealProxy  b4Val;
  RooRealProxy  b5Val;
  RooRealProxy  b6Val;
  RooRealProxy  b7Val;
  RooRealProxy  b8Val;
  RooRealProxy  b9Val;
  RooRealProxy  b10Val;

  RooRealProxy  b1ValIm;
  RooRealProxy  b2ValIm;
  RooRealProxy  b3ValIm;
  RooRealProxy  b4ValIm;
  RooRealProxy  b5ValIm;
  RooRealProxy  b6ValIm;
  RooRealProxy  b7ValIm;
  RooRealProxy  b8ValIm;
  RooRealProxy  b9ValIm;
  RooRealProxy  b10ValIm;

  RooRealProxy Lambda;

  RooRealProxy f_spinz1;
  RooRealProxy f_spinz2;

  virtual void evaluatePolarizationTerms(std::vector<Double_t>& Axxyyterm, const Int_t code, bool isGammaV1=false, bool isGammaV2=false) const = 0;
  
  virtual void calculateCi(std::vector<Double_t>& ciRe, std::vector<Double_t>& ciIm, bool isGammaV1=false, bool isGammaV2=false) const;
  virtual void calculateAmplitudes(
    Double_t& A00Re, Double_t& A00Im,
    Double_t& AppRe, Double_t& AppIm, Double_t& A0pRe, Double_t& A0pIm, Double_t& Ap0Re, Double_t& Ap0Im,
    Double_t& AmmRe, Double_t& AmmIm, Double_t& A0mRe, Double_t& A0mIm, Double_t& Am0Re, Double_t& Am0Im,
    Double_t& ApmRe, Double_t& ApmIm, Double_t& AmpRe, Double_t& AmpIm,
    bool isGammaV1=false, bool isGammaV2=false
    ) const;

};

#endif
