#include "RooSpinZero_7DComplex_withAccep_HVV.h"


using namespace std;
using namespace MELAStreamHelpers;


RooSpinZero_7DComplex_withAccep_HVV::RooSpinZero_7DComplex_withAccep_HVV() : RooSpinZero(){}
RooSpinZero_7DComplex_withAccep_HVV::RooSpinZero_7DComplex_withAccep_HVV(
  const char *name, const char *title,
  modelMeasurables const& _measurables,
  modelParameters const& _parameters,
  modelCouplings const& _couplings,
  accepParameters const& _accepParams,
  RooSpin::VdecayType _Vdecay1, RooSpin::VdecayType _Vdecay2,
  TVar::VerbosityLevel verbosity_
) : RooSpinZero(
  name, title,
  _measurables,
  _parameters,
  _couplings,
  _Vdecay1, _Vdecay2,
  verbosity_
),
aPhi("aPhi", "aPhi", this, (RooAbsReal&)*(_accepParams.aPhi)),
bPhi("bPhi", "bPhi", this, (RooAbsReal&)*(_accepParams.bPhi)),
cPhi("cPhi", "cPhi", this, (RooAbsReal&)*(_accepParams.cPhi)),
dPhi("dPhi", "dPhi", this, (RooAbsReal&)*(_accepParams.dPhi)),
ePhi("ePhi", "ePhi", this, (RooAbsReal&)*(_accepParams.ePhi)),
aPhi1("aPhi1", "aPhi1", this, (RooAbsReal&)*(_accepParams.aPhi1)),
bPhi1("bPhi1", "bPhi1", this, (RooAbsReal&)*(_accepParams.bPhi1)),
cPhi1("cPhi1", "cPhi1", this, (RooAbsReal&)*(_accepParams.cPhi1)),
dPhi1("dPhi1", "dPhi1", this, (RooAbsReal&)*(_accepParams.dPhi1)),
ePhi1("ePhi1", "ePhi1", this, (RooAbsReal&)*(_accepParams.ePhi1)),
aH1("aH1", "aH1", this, (RooAbsReal&)*(_accepParams.aH1)),
bH1("bH1", "bH1", this, (RooAbsReal&)*(_accepParams.bH1)),
cH1("cH1", "cH1", this, (RooAbsReal&)*(_accepParams.cH1)),
dH1("dH1", "dH1", this, (RooAbsReal&)*(_accepParams.dH1)),
eH1("eH1", "eH1", this, (RooAbsReal&)*(_accepParams.eH1)),
aH2("aH2", "aH2", this, (RooAbsReal&)*(_accepParams.aH2)),
bH2("bH2", "bH2", this, (RooAbsReal&)*(_accepParams.bH2)),
cH2("cH2", "cH2", this, (RooAbsReal&)*(_accepParams.cH2)),
dH2("dH2", "dH2", this, (RooAbsReal&)*(_accepParams.dH2)),
eH2("eH2", "eH2", this, (RooAbsReal&)*(_accepParams.eH2)),
aHs("aHs", "aHs", this, (RooAbsReal&)*(_accepParams.aHs)),
bHs("bHs", "bHs", this, (RooAbsReal&)*(_accepParams.bHs)),
cHs("cHs", "cHs", this, (RooAbsReal&)*(_accepParams.cHs)),
dHs("dHs", "dHs", this, (RooAbsReal&)*(_accepParams.dHs)),
eHs("eHs", "eHs", this, (RooAbsReal&)*(_accepParams.eHs)),
aM1("aM1", "aM1", this, (RooAbsReal&)*(_accepParams.aM1)),
bM1("bM1", "bM1", this, (RooAbsReal&)*(_accepParams.bM1)),
cM1("cM1", "cM1", this, (RooAbsReal&)*(_accepParams.cM1)),
dM1("dM1", "dM1", this, (RooAbsReal&)*(_accepParams.dM1)),
aM2("aM2", "aM2", this, (RooAbsReal&)*(_accepParams.aM2)),
bM2("bM2", "bM2", this, (RooAbsReal&)*(_accepParams.bM2)),
cM2("cM2", "cM2", this, (RooAbsReal&)*(_accepParams.cM2)),
dM2("dM2", "dM2", this, (RooAbsReal&)*(_accepParams.dM2)),
ZZ4fOrdering(true)
{}


RooSpinZero_7DComplex_withAccep_HVV::RooSpinZero_7DComplex_withAccep_HVV(
  const RooSpinZero_7DComplex_withAccep_HVV& other, const char* name
) : RooSpinZero(other, name),
aPhi("aPhi", this, other.aPhi),
bPhi("bPhi", this, other.bPhi),
cPhi("cPhi", this, other.cPhi),
dPhi("dPhi", this, other.dPhi),
ePhi("ePhi", this, other.ePhi),
aPhi1("aPhi1", this, other.aPhi1),
bPhi1("bPhi1", this, other.bPhi1),
cPhi1("cPhi1", this, other.cPhi1),
dPhi1("dPhi1", this, other.dPhi1),
ePhi1("ePhi1", this, other.ePhi1),
aH1("aH1", this, other.aH1),
bH1("bH1", this, other.bH1),
cH1("cH1", this, other.cH1),
dH1("dH1", this, other.dH1),
eH1("eH1", this, other.eH1),
aH2("aH2", this, other.aH2),
bH2("bH2", this, other.bH2),
cH2("cH2", this, other.cH2),
dH2("dH2", this, other.dH2),
eH2("eH2", this, other.eH2),
aHs("aHs", this, other.aHs),
bHs("bHs", this, other.bHs),
cHs("cHs", this, other.cHs),
dHs("dHs", this, other.dHs),
eHs("eHs", this, other.eHs),
aM1("aM1", this, other.aM1),
bM1("bM1", this, other.bM1),
cM1("cM1", this, other.cM1),
dM1("dM1", this, other.dM1),
aM2("aM2", this, other.aM2),
bM2("bM2", this, other.bM2),
cM2("cM2", this, other.cM2),
dM2("dM2", this, other.dM2),
ZZ4fOrdering(other.ZZ4fOrdering)
{}

void RooSpinZero_7DComplex_withAccep_HVV::evaluatePolarizationTerms(
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
    Double_t A00_h1int = 4.*aH1/3. + 4.*bH1/15. + 4.*cH1/35. + 4.*dH1/63. + 4.*eH1/99.;
    Double_t Ammpp_h1int = 8.*aH1/3. + 16.*bH1/15. + 24.*cH1/35. + 32.*dH1/63. + 40.*eH1/99.;
    Double_t A0m_h1int = (128.*aH1 + 32.*bH1 + 16.*cH1 + 10.*dH1 + 7.*eH1)*Pi/256.*R1Val;
    Double_t A0p_h1int = (128.*aH1 + 32.*bH1 + 16.*cH1 + 10.*dH1 + 7.*eH1)*Pi/256.*R1Val;
    Double_t Amp_h1int = 4.*aH1/3. + 4.*bH1/15. + 4.*cH1/35. + 4.*dH1/63. + 4.*eH1/99.;

    if (A00term!=0) A00term *= A00_h1int;
    if (Appterm!=0) Appterm *= Ammpp_h1int;
    if (Ammterm!=0) Ammterm *= Ammpp_h1int;
    if (A00ppterm!=0) A00ppterm *= A0p_h1int;
    if (A00mmterm!=0) A00mmterm *= A0m_h1int;
    if (Appmmterm!=0) Appmmterm *= Amp_h1int;
  }
  else{
    Double_t common_fac = (aH1 + bH1*pow(h1, 2) + cH1*pow(h1, 4) +dH1*pow(h1, 6) + eH1*pow(h1, 8));

    if (A00term!=0) A00term *= (1. - pow(h1, 2))*common_fac;
    if (Appterm!=0) Appterm *= (1. + pow(h1, 2) - 2.*h1*R1Val)*common_fac;
    if (Ammterm!=0) Ammterm *= (1. + pow(h1, 2) + 2.*h1*R1Val)*common_fac;
    if (A00ppterm!=0) A00ppterm *= sqrt(1. - pow(h1, 2))*(R1Val - h1)*common_fac;
    if (A00mmterm!=0) A00mmterm *= sqrt(1. - pow(h1, 2))*(R1Val + h1)*common_fac;
    if (Appmmterm!=0) Appmmterm *= (1. - pow(h1, 2))*common_fac;
  }

  if ((code % prime_h2)==0){
    Double_t A00_h2int = 4.*aH2/3. + 4.*bH2/15. + 4.*cH2/35. + 4.*dH2/63. + 4.*eH2/99.;
    Double_t Ammpp_h2int = 8.*aH2/3. + 16.*bH2/15. + 24.*cH2/35. + 32.*dH2/63. + 40.*eH2/99.;
    Double_t A0m_h2int = (128.*aH2 + 32.*bH2 + 16.*cH2 + 10.*dH2 + 7.*eH2)*Pi/256.*R2Val;
    Double_t A0p_h2int = (128.*aH2 + 32.*bH2 + 16.*cH2 + 10.*dH2 + 7.*eH2)*Pi/256.*R2Val;
    Double_t Amp_h2int = 4.*aH2/3. + 4.*bH2/15. + 4.*cH2/35. + 4.*dH2/63. + 4.*eH2/99.;

    if (A00term!=0) A00term *= A00_h2int;
    if (Appterm!=0) Appterm *= Ammpp_h2int;
    if (Ammterm!=0) Ammterm *= Ammpp_h2int;
    if (A00ppterm!=0) A00ppterm *= A0p_h2int;
    if (A00mmterm!=0) A00mmterm *= A0m_h2int;
    if (Appmmterm!=0) Appmmterm *= Amp_h2int;
  }
  else{
    Double_t common_fac = (aH2 + bH2*pow(h2, 2) + cH2*pow(h2, 4) +dH2*pow(h2, 6) + eH2*pow(h2, 8));

    if (A00term!=0) A00term *= (1. - pow(h2, 2))*common_fac;
    if (Appterm!=0) Appterm *= (1. + pow(h2, 2) - 2.*h2*R2Val)*common_fac;
    if (Ammterm!=0) Ammterm *= (1. + pow(h2, 2) + 2.*h2*R2Val)*common_fac;
    if (A00ppterm!=0) A00ppterm *= sqrt(1. - pow(h2, 2))*(R2Val - h2)*common_fac;
    if (A00mmterm!=0) A00mmterm *= sqrt(1. - pow(h2, 2))*(R2Val + h2)*common_fac;
    if (Appmmterm!=0) Appmmterm *= (1. - pow(h2, 2))*common_fac;
  }

  if ((code % prime_hs)==0){
    Double_t A_hsint = 2.*aHs + 2.*bHs/3. + 2.*cHs/5. + 2.*dHs/7. + 2.*eHs/9.; // All of them

    if (A00term!=0) A00term *= A_hsint;
    if (Appterm!=0) Appterm *= A_hsint;
    if (Ammterm!=0) Ammterm *= A_hsint;
    if (A00ppterm!=0) A00ppterm *= A_hsint;
    if (A00mmterm!=0) A00mmterm *= A_hsint;
    if (Appmmterm!=0) Appmmterm *= A_hsint;
  }
  else{
    Double_t common_fac = (aHs + bHs*pow(hs, 2) + cHs*pow(hs, 4) + dHs*pow(hs, 6) +eHs*pow(hs, 8));

    if (A00term!=0) A00term *= common_fac;
    if (Appterm!=0) Appterm *= common_fac;
    if (Ammterm!=0) Ammterm *= common_fac;
    if (A00ppterm!=0) A00ppterm *= common_fac;
    if (A00mmterm!=0) A00mmterm *= common_fac;
    if (Appmmterm!=0) Appmmterm *= common_fac;
  }

  if ((code % prime_Phi)==0){
    Double_t A00mmpp_phiint = 2.*aPhi*Pi;
    Double_t A0p_phiint = bPhi*Pi*cos(phipp);
    Double_t A0m_phiint = bPhi*Pi*cos(phimm);
    Double_t Amp_phiint = cPhi*Pi*cos(phimm - phipp);

    if (A00term!=0) A00term *= A00mmpp_phiint;
    if (Appterm!=0) Appterm *= A00mmpp_phiint;
    if (Ammterm!=0) Ammterm *= A00mmpp_phiint;
    if (A00ppterm!=0) A00ppterm *= A0p_phiint;
    if (A00mmterm!=0) A00mmterm *= A0m_phiint;
    if (Appmmterm!=0) Appmmterm *= Amp_phiint;
  }
  else{
    Double_t common_fac = (aPhi + bPhi*cos(Phi) + cPhi*cos(2*Phi) + dPhi*cos(3*Phi) +ePhi*cos(4*Phi));

    if (A00term!=0) A00term *= common_fac;
    if (Appterm!=0) Appterm *= common_fac;
    if (Ammterm!=0) Ammterm *= common_fac;
    if (A00ppterm!=0) A00ppterm *= cos(Phi + phipp)*common_fac;
    if (A00mmterm!=0) A00mmterm *= cos(Phi - phimm)*common_fac;
    if (Appmmterm!=0) Appmmterm *= cos(2*Phi - phimm + phipp)*common_fac;
  }

  if ((code % prime_Phi1)==0){
    Double_t A_phi1int = 2.*aPhi1*Pi;

    if (A00term!=0) A00term *= A_phi1int;
    if (Appterm!=0) Appterm *= A_phi1int;
    if (Ammterm!=0) Ammterm *= A_phi1int;
    if (A00ppterm!=0) A00ppterm *= A_phi1int;
    if (A00mmterm!=0) A00mmterm *= A_phi1int;
    if (Appmmterm!=0) Appmmterm *= A_phi1int;
  }
  else{
    Double_t common_fac = (aPhi1 + bPhi1*cos(Phi1) + cPhi1*cos(2*Phi1) + dPhi1*cos(3*Phi1) +ePhi1*cos(4*Phi1));

    if (A00term!=0) A00term *= common_fac;
    if (Appterm!=0) Appterm *= common_fac;
    if (Ammterm!=0) Ammterm *= common_fac;
    if (A00ppterm!=0) A00ppterm *= common_fac;
    if (A00mmterm!=0) A00mmterm *= common_fac;
    if (Appmmterm!=0) Appmmterm *= common_fac;
  }

  if (verbosity>=TVar::DEBUG){
    MELAout << "RooSpinZero_7DComplex_withAccep_HVV::evaluatePolarizationTerms( " << code << " , " << VGammaVpmode1 << " , " << VGammaVpmode2 << " ):" << endl;
    MELAout << "\t- |A00|**2 = " << A00term << endl;
    MELAout << "\t- |A++|**2 = " << Appterm << endl;
    MELAout << "\t- |A--|**2 = " << Ammterm << endl;
    MELAout << "\t- |A00||A++| = " << A00ppterm << endl;
    MELAout << "\t- |A00||A--| = " << A00mmterm << endl;
    MELAout << "\t- |A++||A--| = " << Appmmterm << endl;
  }
}

Double_t RooSpinZero_7DComplex_withAccep_HVV::evaluate() const{
  Double_t mV;
  getMVGamV(&mV);
  bool isZZ = (mV >= 90.);
  Double_t epsilon=1e-15;
  Double_t m1_=m1; if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) m1_=0;
  Double_t m2_=m2; if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) m2_=0;
  if (
    (m1_+m2_)>m12 ||
    (isZZ && Vdecay1==Vdecay2 && ZZ4fOrdering && fabs(m2_-mV)<fabs(m1_-mV) && Vdecay2!=RooSpin::kVdecayType_GammaOnshell) ||
    (m1_<=0. && Vdecay1!=RooSpin::kVdecayType_GammaOnshell) ||
    (m2_<=0. && Vdecay2!=RooSpin::kVdecayType_GammaOnshell)
    ) return epsilon;

  Int_t code = intCodeStart;
  if (Vdecay1==RooSpin::kVdecayType_GammaOnshell || Vdecay2==RooSpin::kVdecayType_GammaOnshell){
    code *= prime_Phi;
    if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) code *= prime_h1;
    if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_h2;
    if (Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_Phi1;
  }

  Double_t betaValSq = (1.-(pow(m1_-m2_, 2)/pow(m12, 2)))*(1.-(pow(m1_+m2_, 2)/pow(m12, 2)));
  if (betaValSq<=0.) return epsilon;
  Double_t betaVal = sqrt(betaValSq);

  Double_t term1Coeff = 1;
  Double_t term2Coeff = 1;
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell) term1Coeff = 2.*m1_*GeVunit; // dm**2 = 2m dm
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) term2Coeff = 2.*m2_*GeVunit;

  Double_t value = 0;
  for (int VGammaVpmode1=0; VGammaVpmode1<=2; VGammaVpmode1++){
    for (int VGammaVpmode2=0; VGammaVpmode2<=2; VGammaVpmode2++){
      if (!(
        (VGammaVpmode1==1 || Vdecay1!=RooSpin::kVdecayType_GammaOnshell)
        &&
        (VGammaVpmode2==1 || Vdecay2!=RooSpin::kVdecayType_GammaOnshell)
        )
        ||
        (VGammaVpmode1==1 && VGammaVpmode2==2) || (VGammaVpmode1==2 && VGammaVpmode2==1)
        ||
        !computeNeededAmplitude(VGammaVpmode1, VGammaVpmode2)
        ) continue;
      Double_t val_A00=0, val_App=0, val_Amm=0, val_A0p=0, val_A0m=0, val_Amp=0;
      evaluatePolarizationTerms(val_A00, val_App, val_Amm, val_A0p, val_A0m, val_Amp, code, VGammaVpmode1, VGammaVpmode2);
      value += val_A00 + val_App + val_Amm + val_A0p + val_A0m + val_Amp;
    }
  }
  value *= betaVal*term1Coeff*term2Coeff
    *(1+aM1*m1_+bM1*m1_*m1_+cM1*m1_*m1_*m1_+dM1*m1_*m1_*m1_*m1_)
    *(1+aM2*m2_+bM2*m2_*m2_+cM2*m2_*m2_*m2_+dM2*m2_*m2_*m2_*m2_);

  if (!(value==value)){
    MELAout << "Evaluate NaN=" << value << " at "
      << "h1=" << h1 << '\t'
      << "h2=" << h2 << '\t'
      << "hs=" << hs << '\t'
      << "Phi1=" << Phi1 << '\t'
      << "Phi=" << Phi << '\t'
      << "m1=" << m1 << '\t'
      << "m2=" << m2 << '\t'
      << "m12=" << m12 << '\t'
      << endl;
    MELAout << "Possible sources:\n"
      << "betaVal=" << betaVal << '\t'
      << "term1Coeff=" << term1Coeff << '\t'
      << "term2Coeff=" << term2Coeff
      << endl;
  }
  return value;
}

Int_t RooSpinZero_7DComplex_withAccep_HVV::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code = intCodeStart;
  if (checkFundamentalType(h1)){ if (matchArgs(allVars, analVars, h1) || Vdecay1==RooSpin::kVdecayType_GammaOnshell) code *= prime_h1; }
  if (checkFundamentalType(h2)){ if (matchArgs(allVars, analVars, h2) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_h2; }
  if (checkFundamentalType(hs)){ if (matchArgs(allVars, analVars, hs)) code *= prime_hs; }
  if (checkFundamentalType(Phi)){ if (matchArgs(allVars, analVars, Phi) || Vdecay1==RooSpin::kVdecayType_GammaOnshell || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_Phi; }
  if (checkFundamentalType(Phi1)){ if (matchArgs(allVars, analVars, Phi1) || (Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell)) code *= prime_Phi1; }
  if (code==1) code=0;
  return code;
}
Double_t RooSpinZero_7DComplex_withAccep_HVV::analyticalIntegral(Int_t code, const char* /*rangeName*/) const{
  Double_t mV;
  getMVGamV(&mV);
  bool isZZ = (mV >= 90.);
  Double_t epsilon=1e-10;
  Double_t m1_=m1; if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) m1_=0;
  Double_t m2_=m2; if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) m2_=0;
  if (
    (m1_+m2_)>m12 ||
    (isZZ && Vdecay1==Vdecay2 && ZZ4fOrdering && fabs(m2_-mV)<fabs(m1_-mV) && Vdecay2!=RooSpin::kVdecayType_GammaOnshell) ||
    (m1_<=0. && Vdecay1!=RooSpin::kVdecayType_GammaOnshell) ||
    (m2_<=0. && Vdecay2!=RooSpin::kVdecayType_GammaOnshell)
    ) return epsilon;

  Double_t betaValSq = (1.-(pow(m1_-m2_, 2)/pow(m12, 2)))*(1.-(pow(m1_+m2_, 2)/pow(m12, 2)));
  if (betaValSq<=0.) return epsilon;
  Double_t betaVal = sqrt(betaValSq);

  Double_t term1Coeff = 1;
  Double_t term2Coeff = 1;
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell) term1Coeff = 2.*m1_*GeVunit; // dm**2 = 2m dm
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
  value *= betaVal*term1Coeff*term2Coeff
    *(1+aM1*m1_+bM1*m1_*m1_+cM1*m1_*m1_*m1_+dM1*m1_*m1_*m1_*m1_)
    *(1+aM2*m2_+bM2*m2_*m2_+cM2*m2_*m2_*m2_+dM2*m2_*m2_*m2_*m2_);

  if (!(value==value)){
    MELAout << "Integral NaN=" << value << " at "
      << "h1=" << h1 << '\t'
      << "h2=" << h2 << '\t'
      << "hs=" << hs << '\t'
      << "Phi1=" << Phi1 << '\t'
      << "Phi=" << Phi << '\t'
      << "m1=" << m1 << '\t'
      << "m2=" << m2 << '\t'
      << "m12=" << m12 << '\t'
      << endl;
    MELAout << "Possible sources:\n"
      << "betaVal=" << betaVal << '\t'
      << "term1Coeff=" << term1Coeff << '\t'
      << "term2Coeff=" << term2Coeff
      << endl;
  }
  return value;
}

void RooSpinZero_7DComplex_withAccep_HVV::setZZ4fOrdering(Bool_t flag){ ZZ4fOrdering=flag; }
