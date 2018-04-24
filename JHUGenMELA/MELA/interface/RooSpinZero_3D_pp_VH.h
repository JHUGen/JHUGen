#ifndef ROOSPINZERO_3D_PP_VH
#define ROOSPINZERO_3D_PP_VH

#include "RooSpinZero.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"


class RooSpinZero_3D_pp_VH : public RooSpinZero {
public:

  Double_t sqrts;

  RooSpinZero_3D_pp_VH(){}
  RooSpinZero_3D_pp_VH(
    const char *name, const char *title,
    modelMeasurables _measurables,
    modelParameters _parameters,
    modelCouplings _couplings,
    Double_t _sqrts,
    RooSpin::VdecayType _Vdecay1=RooSpin::kVdecayType_Zll, RooSpin::VdecayType _Vdecay2=RooSpin::kVdecayType_Zll
    );

  RooSpinZero_3D_pp_VH(const RooSpinZero_3D_pp_VH& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooSpinZero_3D_pp_VH(*this, newname); }
  inline virtual ~RooSpinZero_3D_pp_VH(){}

  Double_t evaluate() const;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:


  void evaluatePolarizationTerms(
    Double_t& A00term, Double_t& Appterm, Double_t& Ammterm,
    Double_t& A00ppterm, Double_t& A00mmterm, Double_t& Appmmterm,
    const Int_t code,
    int VGammaVpmode1=0, int VGammaVpmode2=0
  ) const;

  Double_t partonicLuminosity(Double_t mVal, Double_t YVal, Double_t sqrts) const;

};

#endif
