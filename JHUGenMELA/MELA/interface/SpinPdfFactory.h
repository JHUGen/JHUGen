#ifndef SPIN_PDF_FACTORY
#define SPIN_PDF_FACTORY

#include "RooSpin.h"
#include "RooFormulaVar.h"
#include "TString.h"


class SpinPdfFactory{
public:
  RooSpin::modelMeasurables measurables;
  RooSpin::modelParameters parameters;

  SpinPdfFactory(RooSpin::modelMeasurables measurables_, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true);
  virtual ~SpinPdfFactory();

  virtual void getMVGamV(Double_t* mV=0, Double_t* gamV=0) const;
  virtual void getMVprimeGamVprime(Double_t* mV=0, Double_t* gamV=0) const;

  virtual void makeParamsConst(bool yesNo);
  virtual void makeCouplingsConst(bool yesNo)=0;
  virtual void resetHypotheses()=0;
  virtual void resetVdecay(RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_);
  virtual void resetVdecayParams(Double_t mWval, Double_t gamWval, Double_t mZval, Double_t gamZval, Double_t Sin2ThetaWval, Double_t vevval);
  virtual void resetVprimeffCouplings(Double_t gVprimeff_decay1[2], Double_t gVprimeff_decay2[2]); // 0 = left, 1 = right
  virtual void resetVprimeMasses(Double_t mWval, Double_t gamWval, Double_t mZval, Double_t gamZval);
  virtual void resetHiggsMassWidth(Double_t mXval, Double_t gamXval=0.);

  virtual void defaultIntegration(){ PDF_base->defaultIntegration(); }
  virtual void alwaysIntegrate(Int_t code=1){ PDF_base->alwaysIntegrate(code); }

  static void setVariableValue(RooRealVar* var, Double_t val);

protected:

  RooSpin::VdecayType V1decay;
  RooSpin::VdecayType V2decay;

  Bool_t OnshellH;

  RooSpin* PDF_base;

  virtual void initMeasurables(RooSpin::modelMeasurables measurables_);

  virtual void initMassPole();
  virtual void initVdecayParams();

  virtual void destroyMassPole();
  virtual void destroyVdecayParams();

  virtual void initGVals()=0;
  virtual void destroyGVals()=0;
  virtual void initPDF()=0;
  virtual void destroyPDF()=0;

};




#endif



