#ifndef SCALAR_PDF_FACTORY
#define SCALAR_PDF_FACTORY

#include "RooSpinZero.h"
#include "SpinPdfFactory.h"
#include "RooFormulaVar.h"
#include "TString.h"


class ScalarPdfFactory : public SpinPdfFactory{
public:
  RooSpinZero::modelCouplings couplings;

  RooRealVar* g1Frac[7]; // f_a1 = 1.-sum_{i, i>1}{fabs(f_ai)}
  RooRealVar* g2Frac[8];
  RooRealVar* g3Frac[8];
  RooRealVar* g4Frac[8];
  RooRealVar* g1Phase[7]; // phi_a1=0
  RooRealVar* g2Phase[8];
  RooRealVar* g3Phase[8];
  RooRealVar* g4Phase[8];

  RooRealVar* gzgs1Frac[1];
  RooRealVar* gzgs2Frac[1];
  RooRealVar* gzgs3Frac[1];
  RooRealVar* gzgs4Frac[1];
  RooRealVar* gzgs1Phase[1];
  RooRealVar* gzgs2Phase[1];
  RooRealVar* gzgs3Phase[1];
  RooRealVar* gzgs4Phase[1];

  RooRealVar* ggsgs2Frac[1];
  RooRealVar* ggsgs3Frac[1];
  RooRealVar* ggsgs4Frac[1];
  RooRealVar* ggsgs2Phase[1];
  RooRealVar* ggsgs3Phase[1];
  RooRealVar* ggsgs4Phase[1];

  RooRealVar* gvvp1Frac[1];
  RooRealVar* gvvp1Phase[1];

  RooRealVar* gvpvp1Frac[1];
  RooRealVar* gvpvp1Phase[1];

  RooFormulaVar* gFracSum; // sum_{i, i>1}{fabs(f_ai)}
  RooFormulaVar* g1FracInterp[8]; // f_a1^interp = (f_a1<0 ? 0 : f_ai)
  RooFormulaVar* g2FracInterp[8];
  RooFormulaVar* g3FracInterp[8];
  RooFormulaVar* g4FracInterp[8];

  RooFormulaVar* gzgs1FracInterp[1];
  RooFormulaVar* gzgs2FracInterp[1];
  RooFormulaVar* gzgs3FracInterp[1];
  RooFormulaVar* gzgs4FracInterp[1];
  RooFormulaVar* ggsgs2FracInterp[1];
  RooFormulaVar* ggsgs3FracInterp[1];
  RooFormulaVar* ggsgs4FracInterp[1];

  RooFormulaVar* gvvp1FracInterp[1];

  RooFormulaVar* gvpvp1FracInterp[1];

  RooRealVar* gRatioVal[4][8];
  RooRealVar* gZGsRatioVal[4][1];
  RooRealVar* gGsGsRatioVal[3][1];
  RooRealVar* gVVpRatioVal[1][1];
  RooRealVar* gVpVpRatioVal[1][1];

  ScalarPdfFactory(RooSpin::modelMeasurables measurables_, bool acceptance_=false, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true);
  ScalarPdfFactory(
    RooSpin::modelMeasurables measurables_,
    double gRatio_[4][8],
    double gZGsRatio_[4][1],
    double gGsGsRatio_[3][1],
    double gVVpRatio_[1][1],
    double gVpVpRatio_[1][1],
    bool pmf_applied_=false, bool acceptance_=false,
    RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true
  );
  virtual ~ScalarPdfFactory();

  virtual void makeCouplingsConst(bool yesNo=true);
  virtual void addHypothesis(int ig, int ilam, double iphase=0, double altparam_fracval=0);
  virtual void resetHypotheses();
  virtual RooSpinZero* getPDF()=0;

protected:
  int parameterization;
  bool pmf_applied;
  bool acceptance;

  double gRatio[4][8];
  double gZGsRatio[4][1];
  double gGsGsRatio[3][1];
  double gVVpRatio[1][1];
  double gVpVpRatio[1][1];

  virtual void initFractionsPhases();
  virtual void initGVals();

  virtual void destroyFractionsPhases();
  virtual void destroyGVals();

  virtual void initPDF()=0;
  virtual void destroyPDF()=0;

};




#endif



