#include "RooSpinZero_5D_VH.h"


RooSpinZero_5D_VH::RooSpinZero_5D_VH(
  const char *name, const char *title,
  modelMeasurables _measurables,
  modelParameters _parameters,
  modelCouplings _couplings,
  RooSpin::VdecayType _Vdecay1, RooSpin::VdecayType _Vdecay2
  ) : RooSpinZero(
  name, title,
  _measurables,
  _parameters,
  _couplings,
  _Vdecay1, _Vdecay2
  )
{}


RooSpinZero_5D_VH::RooSpinZero_5D_VH(
  const RooSpinZero_5D_VH& other, const char* name
  ) : RooSpinZero(other, name)
{}

void RooSpinZero_5D_VH::evaluatePolarizationTerms(
  Double_t& A00term, Double_t& Appterm, Double_t& Ammterm,
  Double_t& A00ppterm, Double_t& A00mmterm, Double_t& Appmmterm,
  const Int_t code,
  int VGammaVpmode1, int VGammaVpmode2
) const{
  const Double_t Pi = TMath::Pi();

  Double_t R1Val, R2Val;
  calculateVffR1R2(R1Val, R2Val, VGammaVpmode1==1, VGammaVpmode2==1);
  if (VGammaVpmode1==2 || VGammaVpmode2==2){
    Double_t RVp1Val=0, RVp2Val=0;
    calculateVprimeffR1R2(RVp1Val, RVp2Val);
    if (VGammaVpmode1==2) R1Val=RVp1Val;
    if (VGammaVpmode2==2) R1Val=RVp2Val;
  }

  Double_t A00Re, A00Im, AppRe, AppIm, AmmRe, AmmIm;
  calculateAmplitudes(A00Re, A00Im, AppRe, AppIm, AmmRe, AmmIm, VGammaVpmode1, VGammaVpmode2);

  Double_t A00 = A00Im*A00Im + A00Re*A00Re;
  Double_t App = AppIm*AppIm + AppRe*AppRe;
  Double_t Amm = AmmIm*AmmIm + AmmRe*AmmRe;

  Double_t phi00=atan2(A00Im, A00Re);
  Double_t phipp=atan2(AppIm, AppRe)-phi00;
  Double_t phimm=atan2(AmmIm, AmmRe)-phi00;

  Double_t A00_prefactor = 1.;
  Double_t Amm_pp_prefactor = 1.;
  Double_t A00mm_prefactor = 2.;
  Double_t A00pp_prefactor = 2.;
  Double_t Ammpp_prefactor = 2.;

  A00term = A00*A00_prefactor;
  Appterm = App*Amm_pp_prefactor;
  Ammterm = Amm*Amm_pp_prefactor;
  A00ppterm = sqrt(A00*App)*A00pp_prefactor;
  A00mmterm = sqrt(A00*Amm)*A00mm_prefactor;
  Appmmterm = sqrt(Amm*App)*Ammpp_prefactor;

  if ((code % prime_h1)==0){
    Double_t A00_h1int = 4./3.;
    Double_t Ammpp_h1int = 8./3.;
    Double_t A0m_h1int = Pi/2.*R1Val;
    Double_t A0p_h1int = Pi/2.*R1Val;
    Double_t Amp_h1int = 4./3.;

    if (A00term!=0) A00term *= A00_h1int;
    if (Appterm!=0) Appterm *= Ammpp_h1int;
    if (Ammterm!=0) Ammterm *= Ammpp_h1int;
    if (A00ppterm!=0) A00ppterm *= A0p_h1int;
    if (A00mmterm!=0) A00mmterm *= A0m_h1int;
    if (Appmmterm!=0) Appmmterm *= Amp_h1int;
  }
  else{
    if (A00term!=0) A00term *= (1. - pow(h1, 2));
    if (Appterm!=0) Appterm *= (1. + pow(h1, 2) - 2.*h1*R1Val);
    if (Ammterm!=0) Ammterm *= (1. + pow(h1, 2) + 2.*h1*R1Val);
    if (A00ppterm!=0) A00ppterm *= sqrt(1. - pow(h1, 2))*(R1Val - h1);
    if (A00mmterm!=0) A00mmterm *= sqrt(1. - pow(h1, 2))*(R1Val + h1);
    if (Appmmterm!=0) Appmmterm *= (1. - pow(h1, 2));
  }

  if ((code % prime_h2)==0){
    Double_t A00_h2int = 4./3.;
    Double_t Ammpp_h2int = 8./3.;
    Double_t A0m_h2int = Pi/2.*R2Val;
    Double_t A0p_h2int = Pi/2.*R2Val;
    Double_t Amp_h2int = 4./3.;

    if (A00term!=0) A00term *= A00_h2int;
    if (Appterm!=0) Appterm *= Ammpp_h2int;
    if (Ammterm!=0) Ammterm *= Ammpp_h2int;
    if (A00ppterm!=0) A00ppterm *= A0p_h2int;
    if (A00mmterm!=0) A00mmterm *= A0m_h2int;
    if (Appmmterm!=0) Appmmterm *= Amp_h2int;
  }
  else{
    if (A00term!=0) A00term *= (1. - pow(h2, 2));
    if (Appterm!=0) Appterm *= (1. + pow(h2, 2) + 2.*h2*R2Val);
    if (Ammterm!=0) Ammterm *= (1. + pow(h2, 2) - 2.*h2*R2Val);
    if (A00ppterm!=0) A00ppterm *= sqrt(1. - pow(h2, 2))*(R2Val + h2);
    if (A00mmterm!=0) A00mmterm *= sqrt(1. - pow(h2, 2))*(R2Val - h2);
    if (Appmmterm!=0) Appmmterm *= (1. - pow(h2, 2));
  }

  if ((code % prime_hs)==0){
    Double_t A_hsint = 2.;

    if (A00term!=0) A00term *= A_hsint;
    if (Appterm!=0) Appterm *= A_hsint;
    if (Ammterm!=0) Ammterm *= A_hsint;
    if (A00ppterm!=0) A00ppterm *= A_hsint;
    if (A00mmterm!=0) A00mmterm *= A_hsint;
    if (Appmmterm!=0) Appmmterm *= A_hsint;
  } // else *= 1

  if ((code % prime_Phi)==0){
    Double_t A00mmpp_phiint = 2.*Pi;
    Double_t A0p_phiint = 0;
    Double_t A0m_phiint = 0;
    Double_t Amp_phiint = 0;

    if (A00term!=0) A00term *= A00mmpp_phiint;
    if (Appterm!=0) Appterm *= A00mmpp_phiint;
    if (Ammterm!=0) Ammterm *= A00mmpp_phiint;
    if (A00ppterm!=0) A00ppterm *= A0p_phiint;
    if (A00mmterm!=0) A00mmterm *= A0m_phiint;
    if (Appmmterm!=0) Appmmterm *= Amp_phiint;
  }
  else{
    //if (A00term!=0) A00term *= 1.;
    //if (Appterm!=0) Appterm *= 1.;
    //if (Ammterm!=0) Ammterm *= 1.;
    if (A00ppterm!=0) A00ppterm *= cos(Phi + phipp);
    if (A00mmterm!=0) A00mmterm *= cos(Phi - phimm);
    if (Appmmterm!=0) Appmmterm *= cos(2*Phi - phimm + phipp);
  }

  if ((code % prime_Phi1)==0){
    Double_t A_phi1int = 2.*Pi;

    if (A00term!=0) A00term *= A_phi1int;
    if (Appterm!=0) Appterm *= A_phi1int;
    if (Ammterm!=0) Ammterm *= A_phi1int;
    if (A00ppterm!=0) A00ppterm *= A_phi1int;
    if (A00mmterm!=0) A00mmterm *= A_phi1int;
    if (Appmmterm!=0) Appmmterm *= A_phi1int;
  }
  // else *= 1
}

Double_t RooSpinZero_5D_VH::evaluate() const{
  Double_t epsilon=1e-15;
  Double_t m1_=m1; if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) return epsilon;
  Double_t m2_=m2; if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) m2_=0;
  if ((m12+m2_) > m1_ || (m2_ <= 0. && Vdecay2!=RooSpin::kVdecayType_GammaOnshell) || m1_ <= 0.) return epsilon;

  Int_t code = intCodeStart;
  if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_h2*prime_Phi;

  Double_t betaValSq = (1.-(pow(m12-m2_, 2)/pow(m1_, 2)))*(1.-(pow(m12+m2_, 2)/pow(m1_, 2)));
  if (betaValSq<=0.) return epsilon;
  Double_t betaVal = sqrt(betaValSq);

  Double_t term1Coeff = 1;
  Double_t term2Coeff = 1;
  term1Coeff = pow(m1_*GeVunit, -2);
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) term2Coeff = 2.*m2_*GeVunit;

  Double_t value = 0;
  for (int VGammaVpmode1=0; VGammaVpmode1<=2; VGammaVpmode1++){
    for (int VGammaVpmode2=0; VGammaVpmode2<=2; VGammaVpmode2++){
      if (!(
        Vdecay1!=RooSpin::kVdecayType_GammaOnshell
        &&
        (VGammaVpmode2==1 || Vdecay2!=RooSpin::kVdecayType_GammaOnshell)
        )
        ||
        VGammaVpmode1==1
        ||
        (VGammaVpmode1==2 && VGammaVpmode2==1)
        ||
        !computeNeededAmplitude(VGammaVpmode1, VGammaVpmode2)
        ) continue;
      Double_t val_A00=0, val_App=0, val_Amm=0, val_A0p=0, val_A0m=0, val_Amp=0;
      evaluatePolarizationTerms(val_A00, val_App, val_Amm, val_A0p, val_A0m, val_Amp, code, VGammaVpmode1, VGammaVpmode2);
      value += val_A00 + val_App + val_Amm + val_A0p + val_A0m + val_Amp;
    }
  }
  value *= betaVal*term1Coeff*term2Coeff;
  return value;
}

Int_t RooSpinZero_5D_VH::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code = intCodeStart;
  if (checkFundamentalType(h1)){ if (matchArgs(allVars, analVars, h1)) code *= prime_h1; }
  if (checkFundamentalType(h2)){ if (matchArgs(allVars, analVars, h2) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_h2; }
  if (checkFundamentalType(hs)){ if (matchArgs(allVars, analVars, hs)) code *= prime_hs; }
  if (checkFundamentalType(Phi)){ if (matchArgs(allVars, analVars, Phi) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_Phi; }
  if (checkFundamentalType(Phi1)){ if (matchArgs(allVars, analVars, Phi1) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_Phi1; }
  if (code==1) code=0;
  return code;
}
Double_t RooSpinZero_5D_VH::analyticalIntegral(Int_t code, const char* /*rangeName*/) const{
  Double_t epsilon=1e-10;
  Double_t m1_=m1; if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) return epsilon;
  Double_t m2_=m2; if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) m2_=0;
  if ((m12+m2_) > m1_ || (m2_ <= 0. && Vdecay2!=RooSpin::kVdecayType_GammaOnshell) || m1_ <= 0.) return epsilon;

  Double_t betaValSq = (1.-(pow(m12-m2_, 2)/pow(m1_, 2)))*(1.-(pow(m12+m2_, 2)/pow(m1_, 2)));
  if (betaValSq<=0.) return epsilon;
  Double_t betaVal = sqrt(betaValSq);

  Double_t term1Coeff = 1;
  Double_t term2Coeff = 1;
  term1Coeff = pow(m1_*GeVunit, -2);
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) term2Coeff = 2.*m2_*GeVunit;

  Double_t value = 0;
  for (int VGammaVpmode1=0; VGammaVpmode1<=2; VGammaVpmode1++){
    for (int VGammaVpmode2=0; VGammaVpmode2<=2; VGammaVpmode2++){
      if (!(
        (VGammaVpmode1==0 && VGammaVpmode2==0 && Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell)
        ||
        (VGammaVpmode1==0 && VGammaVpmode2==1 && Vdecay1!=RooSpin::kVdecayType_GammaOnshell)
        ||
        (VGammaVpmode1==1 && VGammaVpmode2==0 && Vdecay2!=RooSpin::kVdecayType_GammaOnshell)
        ||
        (VGammaVpmode1==1 && VGammaVpmode2==1)
        ||
        (VGammaVpmode1==0 && VGammaVpmode2==2 && Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell)
        ||
        (VGammaVpmode1==2 && VGammaVpmode2==0 && Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell)
        ||
        (VGammaVpmode1==2 && VGammaVpmode2==2 && Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell)
        )) continue;
      Double_t val_A00=0, val_App=0, val_Amm=0, val_A0p=0, val_A0m=0, val_Amp=0;
      evaluatePolarizationTerms(val_A00, val_App, val_Amm, val_A0p, val_A0m, val_Amp, code, VGammaVpmode1, VGammaVpmode2);
      value += val_A00 + val_App + val_Amm + val_A0p + val_A0m + val_Amp;
    }
  }
  value *= betaVal*term1Coeff*term2Coeff;
  return value;
}




