#ifndef ROOSPINZERO
#define ROOSPINZERO

#include "RooSpin.h"


class RooSpinZero : public RooSpin {
public:

  struct modelCouplings{
    RooAbsReal* g1List[8][2];
    RooAbsReal* g2List[8][2];
    RooAbsReal* g3List[8][2];
    RooAbsReal* g4List[8][2];

    RooAbsReal* gzgs1List[1][2]; // ghzgs_prime2
    RooAbsReal* gzgs2List[1][2];
    RooAbsReal* gzgs3List[1][2];
    RooAbsReal* gzgs4List[1][2];
    RooAbsReal* ggsgs2List[1][2];
    RooAbsReal* ggsgs3List[1][2];
    RooAbsReal* ggsgs4List[1][2];

    RooRealVar* Lambda;
    RooRealVar* Lambda_zgs1;
    RooRealVar* Lambda_z1;
    RooRealVar* Lambda_z2;
    RooRealVar* Lambda_z3;
    RooRealVar* Lambda_z4;
    RooRealVar* Lambda_Q;

    RooAbsReal* Lambda_z1qsq[SIZE_HVV_CQSQ];
    RooAbsReal* Lambda_z2qsq[SIZE_HVV_CQSQ];
    RooAbsReal* Lambda_z3qsq[SIZE_HVV_CQSQ];
    RooAbsReal* Lambda_z4qsq[SIZE_HVV_CQSQ];
    RooAbsReal* cLambda_qsq[SIZE_HVV_CQSQ];

    RooAbsReal* gvvp1List[1][2];
    RooAbsReal* gvpvp1List[1][2];
  };

  RooSpinZero();
  RooSpinZero(
    const char* name, const char* title,
    modelMeasurables const& _measurables,
    modelParameters const& _parameters,
    modelCouplings const& _couplings,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll,
    TVar::VerbosityLevel verbosity_=TVar::ERROR
  );

  RooSpinZero(const RooSpinZero& other, const char* name=0);
  virtual TObject* clone(const char* newname) const = 0;
  inline virtual ~RooSpinZero(){}

  virtual Double_t evaluate() const = 0;
  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const = 0;
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const = 0;
  virtual void printParameters() const;

protected:

  RooRealProxy g1Val;
  RooRealProxy g2Val;
  RooRealProxy g3Val;
  RooRealProxy g4Val;
  RooRealProxy g1_primeVal;
  RooRealProxy g2_primeVal;
  RooRealProxy g3_primeVal;
  RooRealProxy g4_primeVal;
  RooRealProxy g1_prime2Val;
  RooRealProxy g2_prime2Val;
  RooRealProxy g3_prime2Val;
  RooRealProxy g4_prime2Val;
  RooRealProxy g1_prime3Val;
  RooRealProxy g2_prime3Val;
  RooRealProxy g3_prime3Val;
  RooRealProxy g4_prime3Val;
  RooRealProxy g1_prime4Val;
  RooRealProxy g2_prime4Val;
  RooRealProxy g3_prime4Val;
  RooRealProxy g4_prime4Val;
  RooRealProxy g1_prime5Val;
  RooRealProxy g2_prime5Val;
  RooRealProxy g3_prime5Val;
  RooRealProxy g4_prime5Val;
  RooRealProxy g1_prime6Val;
  RooRealProxy g2_prime6Val;
  RooRealProxy g3_prime6Val;
  RooRealProxy g4_prime6Val;
  RooRealProxy g1_prime7Val;
  RooRealProxy g2_prime7Val;
  RooRealProxy g3_prime7Val;
  RooRealProxy g4_prime7Val;

  RooRealProxy gzgs1_prime2Val;
  RooRealProxy gzgs2Val;
  RooRealProxy gzgs3Val;
  RooRealProxy gzgs4Val;
  RooRealProxy ggsgs2Val;
  RooRealProxy ggsgs3Val;
  RooRealProxy ggsgs4Val;


  RooRealProxy g1ValIm;
  RooRealProxy g2ValIm;
  RooRealProxy g3ValIm;
  RooRealProxy g4ValIm;
  RooRealProxy g1_primeValIm;
  RooRealProxy g2_primeValIm;
  RooRealProxy g3_primeValIm;
  RooRealProxy g4_primeValIm;
  RooRealProxy g1_prime2ValIm;
  RooRealProxy g2_prime2ValIm;
  RooRealProxy g3_prime2ValIm;
  RooRealProxy g4_prime2ValIm;
  RooRealProxy g1_prime3ValIm;
  RooRealProxy g2_prime3ValIm;
  RooRealProxy g3_prime3ValIm;
  RooRealProxy g4_prime3ValIm;
  RooRealProxy g1_prime4ValIm;
  RooRealProxy g2_prime4ValIm;
  RooRealProxy g3_prime4ValIm;
  RooRealProxy g4_prime4ValIm;
  RooRealProxy g1_prime5ValIm;
  RooRealProxy g2_prime5ValIm;
  RooRealProxy g3_prime5ValIm;
  RooRealProxy g4_prime5ValIm;
  RooRealProxy g1_prime6ValIm;
  RooRealProxy g2_prime6ValIm;
  RooRealProxy g3_prime6ValIm;
  RooRealProxy g4_prime6ValIm;
  RooRealProxy g1_prime7ValIm;
  RooRealProxy g2_prime7ValIm;
  RooRealProxy g3_prime7ValIm;
  RooRealProxy g4_prime7ValIm;

  RooRealProxy gzgs1_prime2ValIm;
  RooRealProxy gzgs2ValIm;
  RooRealProxy gzgs3ValIm;
  RooRealProxy gzgs4ValIm;
  RooRealProxy ggsgs2ValIm;
  RooRealProxy ggsgs3ValIm;
  RooRealProxy ggsgs4ValIm;


  RooRealProxy Lambda;
  RooRealProxy Lambda_zgs1;
  RooRealProxy Lambda_z1;
  RooRealProxy Lambda_z2;
  RooRealProxy Lambda_z3;
  RooRealProxy Lambda_z4;
  RooRealProxy Lambda_Q;

  RooRealProxy Lambda_z11;
  RooRealProxy Lambda_z21;
  RooRealProxy Lambda_z31;
  RooRealProxy Lambda_z41;
  RooRealProxy Lambda_z12;
  RooRealProxy Lambda_z22;
  RooRealProxy Lambda_z32;
  RooRealProxy Lambda_z42;
  RooRealProxy Lambda_z10;
  RooRealProxy Lambda_z20;
  RooRealProxy Lambda_z30;
  RooRealProxy Lambda_z40;
  RooRealProxy cz_q1sq;
  RooRealProxy cz_q2sq;
  RooRealProxy cz_q12sq;


  RooRealProxy gvvp1Val;
  RooRealProxy gvpvp1Val;

  RooRealProxy gvvp1ValIm;
  RooRealProxy gvpvp1ValIm;

  virtual void evaluatePolarizationTerms(
    Double_t& A00term, Double_t& Appterm, Double_t& Ammterm,
    Double_t& A00ppterm, Double_t& A00mmterm, Double_t& Appmmterm,
    const Int_t code,
    int VGammaVpmode1=0, int VGammaVpmode2=0
  ) const = 0;

  virtual void calculateAi(
    Double_t& a1Re, Double_t& a1Im, Double_t& a2Re, Double_t& a2Im, Double_t& a3Re, Double_t& a3Im,
    int VGammaVpmode1=0, int VGammaVpmode2=0
  ) const;
  virtual void calculateAmplitudes(
    Double_t& A00Re, Double_t& A00Im, Double_t& AppRe, Double_t& AppIm, Double_t& AmmRe, Double_t& AmmIm,
    int VGammaVpmode1=0, int VGammaVpmode2=0
  ) const;


  virtual Bool_t computeNeededAmplitude(int VGammaVpmode1, int VGammaVpmode2) const final;

};

#endif
