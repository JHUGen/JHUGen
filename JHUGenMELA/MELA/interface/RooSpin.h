#ifndef ROOSPIN
#define ROOSPIN

#include <cmath>
#include <vector>
#include "TVar.hh"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsCategory.h"
#include "Riostream.h" 
#include "TMath.h"
#include "TCouplingsBase.hh"

using namespace TMath;
using namespace std;


namespace AnaMelaHelpers{
  void multiplyComplexNumbers(std::vector<Double_t> reals, std::vector<Double_t> imags, Double_t& resRe, Double_t& resIm);
}


class RooSpin : public RooAbsPdf {
public:

  enum VdecayType{
    kVdecayType_Wany=-1,
    kVdecayType_GammaOnshell=0,
    kVdecayType_Zll=1,
    kVdecayType_Znn=2,
    kVdecayType_Zuu=3,
    kVdecayType_Zdd=4,
    kVdecayType_Zud=5
  };

  enum{
    prime_h1=2,
    prime_h2=3,
    prime_hs=5,
    prime_Phi=7,
    prime_Phi1=11,
    prime_m1=13,
    prime_m2=17,
    prime_m12=19,
    prime_Y=23
  };

  struct modelMeasurables{
    RooAbsReal* h1;
    RooAbsReal* h2;
    RooAbsReal* hs;
    RooAbsReal* Phi;
    RooAbsReal* Phi1;
    RooAbsReal* m1;
    RooAbsReal* m2;
    RooAbsReal* m12;
    RooAbsReal* Y;
  };
  struct modelParameters{
    RooAbsReal* mX;
    RooAbsReal* gamX;
    RooAbsReal* mW;
    RooAbsReal* gamW;
    RooAbsReal* mZ;
    RooAbsReal* gamZ;
    RooAbsReal* mWprime;
    RooAbsReal* gamWprime;
    RooAbsReal* mZprime;
    RooAbsReal* gamZprime;
    RooAbsReal* Sin2ThetaW;
    RooAbsReal* vev;
    RooAbsReal* gVprimeff_decay1_left;
    RooAbsReal* gVprimeff_decay1_right;
    RooAbsReal* gVprimeff_decay2_left;
    RooAbsReal* gVprimeff_decay2_right;
  };

  RooSpin();
  RooSpin(
    const char* name, const char* title,
    modelMeasurables _measurables,
    modelParameters _parameters,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll
    );
  RooSpin(const RooSpin& other, const char* name=0);
  inline virtual ~RooSpin(){}

  virtual TObject* clone(const char* newname) const = 0;

  virtual Double_t evaluate() const = 0;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const = 0;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const = 0;

  virtual void setDecayModes(RooSpin::VdecayType Vdecay1_, RooSpin::VdecayType Vdecay2_){ Vdecay1=Vdecay1_; Vdecay2=Vdecay2_; }
  virtual void getMVGamV(Double_t* mV=0, Double_t* gamV=0) const;
  virtual void getMVprimeGamVprime(Double_t* mV=0, Double_t* gamV=0) const;

  virtual void defaultIntegration(){ intCodeStart=1; }
  virtual void alwaysIntegrate(Int_t code=1);


protected:

  RooRealProxy h1;
  RooRealProxy h2;
  RooRealProxy Phi;
  RooRealProxy m1;
  RooRealProxy m2;
  RooRealProxy m12;
  RooRealProxy hs;
  RooRealProxy Phi1;
  RooRealProxy Y;

  RooRealProxy mX;
  RooRealProxy gamX;
  RooRealProxy mW;
  RooRealProxy gamW;
  RooRealProxy mZ;
  RooRealProxy gamZ;
  RooRealProxy mWprime;
  RooRealProxy gamWprime;
  RooRealProxy mZprime;
  RooRealProxy gamZprime;
  RooRealProxy Sin2ThetaW;
  RooRealProxy vev;
  RooRealProxy gVprimeff_decay1_left;
  RooRealProxy gVprimeff_decay1_right;
  RooRealProxy gVprimeff_decay2_left;
  RooRealProxy gVprimeff_decay2_right;

  RooSpin::VdecayType Vdecay1;
  RooSpin::VdecayType Vdecay2;

  Int_t intCodeStart;
  const Double_t GeVunit;

  virtual void calculatePropagator(Double_t& propRe, Double_t& propIm, Double_t mass, Int_t propType=1) const;
  virtual void calculateVffGVGA(Double_t& gV, Double_t& gA, RooSpin::VdecayType Vdecay, bool isGamma=false) const;
  virtual void calculateVffR1R2(Double_t& R1Val, Double_t& R2Val, bool isGammaV1=false, bool isGammaV2=false) const;
  virtual Double_t calculateAmplitudeScale(int VGammaVpmode1=0, int VGammaVpmode2=0) const;

  virtual void calculateVprimeffGVGA(Double_t& gV, Double_t& gA, int whichVprime/*1 or 2*/) const;
  virtual void calculateVprimeffR1R2(Double_t& R1Val, Double_t& R2Val) const;

  virtual void setProxies(modelMeasurables _measurables);
  virtual void setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr);
  virtual Bool_t checkFundamentalType(const RooRealProxy& proxy) const;

  // Check if some amplitudes are needed, otherwise don't even compute them
  virtual Bool_t computeNeededAmplitude(int VGammaVpmode1, int VGammaVpmode2) const { return true; }

};

#endif
