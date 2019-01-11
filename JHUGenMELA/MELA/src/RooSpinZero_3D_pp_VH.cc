#include "RooSpinZero_3D_pp_VH.h"


RooSpinZero_3D_pp_VH::RooSpinZero_3D_pp_VH() : RooSpinZero(){}
RooSpinZero_3D_pp_VH::RooSpinZero_3D_pp_VH(
  const char *name, const char *title,
  modelMeasurables const& _measurables,
  modelParameters const& _parameters,
  modelCouplings const& _couplings,
  Double_t _sqrts,
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
sqrts(_sqrts)
{}


RooSpinZero_3D_pp_VH::RooSpinZero_3D_pp_VH(
  const RooSpinZero_3D_pp_VH& other, const char* name
) : RooSpinZero(other, name),
sqrts(other.sqrts)
{}


void RooSpinZero_3D_pp_VH::evaluatePolarizationTerms(
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


Double_t RooSpinZero_3D_pp_VH::evaluate() const{
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
  Double_t plumi = partonicLuminosity(m1_, Y, sqrts);

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
  value *= betaVal*term1Coeff*term2Coeff*plumi;
  return value;
}

Int_t RooSpinZero_3D_pp_VH::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code = intCodeStart;
  if (checkFundamentalType(h1)){ if (matchArgs(allVars, analVars, h1)) code *= prime_h1; }
  if (checkFundamentalType(h2)){ if (matchArgs(allVars, analVars, h2) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_h2; }
  if (checkFundamentalType(hs)){ if (matchArgs(allVars, analVars, hs)) code *= prime_hs; }
  if (checkFundamentalType(Phi)){ if (matchArgs(allVars, analVars, Phi) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_Phi; }
  if (checkFundamentalType(Phi1)){ if (matchArgs(allVars, analVars, Phi1) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_Phi1; }
  if (code==1) code=0;
  return code;
}
Double_t RooSpinZero_3D_pp_VH::analyticalIntegral(Int_t code, const char* /*rangeName*/) const{
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
  Double_t plumi = partonicLuminosity(m1_, Y, sqrts);

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
  value *= betaVal*term1Coeff*term2Coeff*plumi;
  return value;
}

Double_t RooSpinZero_3D_pp_VH::partonicLuminosity(Double_t mVal, Double_t YVal, Double_t sqrtsVal) const{
  Double_t Q = mVal;
  Double_t xa0 = mVal/sqrtsVal; // -> E E <- => xa0=xb0=p/Ebeam (Ebeam = sqrts/2)
  Double_t xa = exp(YVal)*xa0; // -> Ea Eb <- => xa=Ea/Ebeam
  Double_t xb = exp(-YVal)*xa0; // -> Ea Eb <- => xb=Eb/Ebeam
  Double_t prefactor = 1.;
  if (
    (
    (mVal <= 600. && fabs(YVal) > 20.*pow(mVal, -0.32))
    ||
    (mVal > 600. && fabs(YVal) > 21.*pow(mVal, -0.34))
    ) && Vdecay1==RooSpin::kVdecayType_Wany // !WH
    ){
    xa = xa0;
    xb = xa0;
    prefactor = 1e-5;
  }

  Double_t weightu = 0.5;
  Double_t weightd = 0.5;
  Double_t weightc = 1.0;
  Double_t weights = 1.0;
  Double_t weightb = 1.0;


  // PDF parameters
  // up params
  Double_t u0par0 = 0.03134; Double_t u0par1 =-2.068e-05; Double_t u0par2 = 1.283e-08;
  Double_t u1par0 = 0.9; Double_t u1par1 =-0.0004307; Double_t u1par2 = 2.458e-07;
  Double_t u2par0 =-0.1369; Double_t u2par1 = 0.003423; Double_t u2par2 =-2.155e-06;
  Double_t u3par0 =-0.4013; Double_t u3par1 =-0.0002574; Double_t u3par2 = 1.561e-07;
  Double_t u4par0 = 0.5782; Double_t u4par1 =-0.004728; Double_t u4par2 = 2.906e-06;
  Double_t ubar0par0 = 0.02856; Double_t ubar0par1 =-2.112e-05; Double_t ubar0par2 = 1.272e-08;
  Double_t ubar1par0 =-0.06822; Double_t ubar1par1 = 3.172e-05; Double_t ubar1par2 =-2.008e-08;
  Double_t ubar2par0 = 0.1967; Double_t ubar2par1 =-0.000118; Double_t ubar2par2 = 6.871e-08;
  Double_t ubar3par0 =-0.2251; Double_t ubar3par1 = 0.0001295; Double_t ubar3par2 =-7.181e-08;
  Double_t ubar4par0 =-0.4068; Double_t ubar4par1 =-0.0002956; Double_t ubar4par2 = 1.783e-07;
  Double_t ubar5par0 =-2.251; Double_t ubar5par1 =-0.0001699; Double_t ubar5par2 = 1.492e-07;
  // down params
  Double_t d0par0 = 0.03278; Double_t d0par1 =-2.915e-05; Double_t d0par2 = 1.809e-08;
  Double_t d1par0 = 0.479; Double_t d1par1 =-0.0002559; Double_t d1par2 = 1.557e-07;
  Double_t d2par0 =-0.5972; Double_t d2par1 = 0.0003118; Double_t d2par2 =-1.905e-07;
  Double_t d3par0 =-0.3892; Double_t d3par1 =-0.000317; Double_t d3par2 = 1.944e-07;
  Double_t d4par0 = 0.5007; Double_t d4par1 =-0.001665; Double_t d4par2 = 9.895e-07;
  Double_t dbar0par0 = 0.02328; Double_t dbar0par1 =-1.367e-05; Double_t dbar0par2 = 8.246e-09;
  Double_t dbar1par0 = 0.09422; Double_t dbar1par1 =-0.0001019; Double_t dbar1par2 = 6.375e-08;
  Double_t dbar2par0 =-0.5296; Double_t dbar2par1 = 0.000466; Double_t dbar2par2 =-2.896e-07;
  Double_t dbar3par0 = 0.5354; Double_t dbar3par1 =-0.0004404; Double_t dbar3par2 = 2.728e-07;
  Double_t dbar4par0 =-0.4386; Double_t dbar4par1 =-0.0002605; Double_t dbar4par2 = 1.582e-07;
  Double_t dbar5par0 =-1.289; Double_t dbar5par1 =-0.001618; Double_t dbar5par2 = 9.601e-07;
  // charm, strange, bottom params
  Double_t c0par0 = 0.01829; Double_t c0par1 =-6.93e-06; Double_t c0par2 = 3.796e-09;
  Double_t c1par0 = 0.03081; Double_t c1par1 = 4.325e-05; Double_t c1par2 =-3.95e-08;
  Double_t c2par0 = 0.5398; Double_t c2par1 =-4.284e-05; Double_t c2par2 =-1.362e-08;
  Double_t c3par0 =-0.5986; Double_t c3par1 = 0.002565; Double_t c3par2 =-1.937e-06;
  Double_t c4par0 =-0.4534; Double_t c4par1 =-0.0002329; Double_t c4par2 = 1.343e-07;
  Double_t c5par0 =-8.657; Double_t c5par1 =-0.005157; Double_t c5par2 = 3.68e-06;
  Double_t s0par0 = 0.01312; Double_t s0par1 =-3.743e-06; Double_t s0par2 = 2.076e-09;
  Double_t s1par0 =-0.001416; Double_t s1par1 =-7.649e-06; Double_t s1par2 = 4.757e-09;
  Double_t s2par0 = 0.2864; Double_t s2par1 =-6.693e-05; Double_t s2par2 = 3.566e-08;
  Double_t s3par0 =-0.4857; Double_t s3par1 =-0.000253; Double_t s3par2 = 1.541e-07;
  Double_t s4par0 =-10.33; Double_t s4par1 =-0.001601; Double_t s4par2 = 9.718e-07;
  Double_t b0par0 = 0.005934; Double_t b0par1 = 2.516e-06; Double_t b0par2 =-1.828e-09;
  Double_t b1par0 =-0.003063; Double_t b1par1 =-6.761e-06; Double_t b1par2 = 4.298e-09;
  Double_t b2par0 = 0.1174; Double_t b2par1 = 3.752e-05; Double_t b2par2 =-2.863e-08;
  Double_t b3par0 =-0.5549; Double_t b3par1 =-0.0002205; Double_t b3par2 = 1.334e-07;
  Double_t b4par0 =-10.18; Double_t b4par1 =-0.001136; Double_t b4par2 = 6.931e-07;

  // PDF definition
  Double_t up0 = u0par0 + u0par1*Q + u0par2*Q*Q;
  Double_t up1 = u1par0 + u1par1*Q + u1par2*Q*Q;
  Double_t up2 = u2par0 + u2par1*Q + u2par2*Q*Q;
  Double_t up3 = u3par0 + u3par1*Q + u3par2*Q*Q;
  Double_t up4 = u4par0 + u4par1*Q + u4par2*Q*Q;
  Double_t antiup0 = ubar0par0 + ubar0par1*Q + ubar0par2*Q*Q;
  Double_t antiup1 = ubar1par0 + ubar1par1*Q + ubar1par2*Q*Q;
  Double_t antiup2 = ubar2par0 + ubar2par1*Q + ubar2par2*Q*Q;
  Double_t antiup3 = ubar3par0 + ubar3par1*Q + ubar3par2*Q*Q;
  Double_t antiup4 = ubar4par0 + ubar4par1*Q + ubar4par2*Q*Q;
  Double_t antiup5 = ubar5par0 + ubar5par1*Q + ubar5par2*Q*Q;
  Double_t down0 = d0par0 + d0par1*Q + d0par2*Q*Q;
  Double_t down1 = d1par0 + d1par1*Q + d1par2*Q*Q;
  Double_t down2 = d2par0 + d2par1*Q + d2par2*Q*Q;
  Double_t down3 = d3par0 + d3par1*Q + d3par2*Q*Q;
  Double_t down4 = d4par0 + d4par1*Q + d4par2*Q*Q;
  Double_t antidown0 = dbar0par0 + dbar0par1*Q + dbar0par2*Q*Q;
  Double_t antidown1 = dbar1par0 + dbar1par1*Q + dbar1par2*Q*Q;
  Double_t antidown2 = dbar2par0 + dbar2par1*Q + dbar2par2*Q*Q;
  Double_t antidown3 = dbar3par0 + dbar3par1*Q + dbar3par2*Q*Q;
  Double_t antidown4 = dbar4par0 + dbar4par1*Q + dbar4par2*Q*Q;
  Double_t antidown5 = dbar5par0 + dbar5par1*Q + dbar5par2*Q*Q;
  Double_t charm0 = c0par0 + c0par1*Q + c0par2*Q*Q;
  Double_t charm1 = c1par0 + c1par1*Q + c1par2*Q*Q;
  Double_t charm2 = c2par0 + c2par1*Q + c2par2*Q*Q;
  Double_t charm3 = c3par0 + c3par1*Q + c3par2*Q*Q;
  Double_t charm4 = c4par0 + c4par1*Q + c4par2*Q*Q;
  Double_t charm5 = c5par0 + c5par1*Q + c5par2*Q*Q;
  Double_t strange0 = s0par0 + s0par1*Q + s0par2*Q*Q;
  Double_t strange1 = s1par0 + s1par1*Q + s1par2*Q*Q;
  Double_t strange2 = s2par0 + s2par1*Q + s2par2*Q*Q;
  Double_t strange3 = s3par0 + s3par1*Q + s3par2*Q*Q;
  Double_t strange4 = s4par0 + s4par1*Q + s4par2*Q*Q;
  Double_t bottom0 = b0par0 + b0par1*Q + b0par2*Q*Q;
  Double_t bottom1 = b1par0 + b1par1*Q + b1par2*Q*Q;
  Double_t bottom2 = b2par0 + b2par1*Q + b2par2*Q*Q;
  Double_t bottom3 = b3par0 + b3par1*Q + b3par2*Q*Q;
  Double_t bottom4 = b4par0 + b4par1*Q + b4par2*Q*Q;

  Double_t FuncAu1 = (up0+up1*xa+up2*pow(xa, 2))*pow((1-xa), 4)*pow(xa, up3)*exp(1.0+up4*xa);
  Double_t FuncBu1 = (antiup0+antiup1*xb+antiup2*pow(xb, 2)+antiup3*pow(xb, 3))*pow((1-xb), 4)*pow(xb, antiup4)*exp(1.0+antiup5*xb);
  Double_t FuncAu2 = (up0+up1*xb+up2*pow(xb, 2))*pow((1-xb), 4)*pow(xb, up3)*exp(1.0+up4*xb);
  Double_t FuncBu2 = (antiup0+antiup1*xa+antiup2*pow(xa, 2)+antiup3*pow(xa, 3))*pow((1-xa), 4)*pow(xa, antiup4)*exp(1.0+antiup5*xa);

  Double_t FuncAd1 = (down0+down1*xa+down2*pow(xa, 2))*pow((1-xa), 4)*pow(xa, down3)*exp(1.0+down4*xa);
  Double_t FuncBd1 = (antidown0+antidown1*xb+antidown2*pow(xb, 2)+antidown3*pow(xb, 3))*pow((1-xb), 4)*pow(xb, antidown4)*exp(1.0+antidown5*xb);
  Double_t FuncAd2 = (down0+down1*xb+down2*pow(xb, 2))*pow((1-xb), 4)*pow(xb, down3)*exp(1.0+down4*xb);
  Double_t FuncBd2 = (antidown0+antidown1*xa+antidown2*pow(xa, 2)+antidown3*pow(xa, 3))*pow((1-xa), 4)*pow(xa, antidown4)*exp(1.0+antidown5*xa);

  Double_t Funcca = (charm0+charm1*xa+charm2*pow(xa, 2)+charm3*pow(xa, 3))*pow((1-xa), 4)*pow(xa, charm4)*exp(1.0+charm5*xa);
  Double_t Funccb = (charm0+charm1*xb+charm2*pow(xb, 2)+charm3*pow(xb, 3))*pow((1-xb), 4)*pow(xb, charm4)*exp(1.0+charm5*xb);
  Double_t Funcsa = (strange0+strange1*xa+strange2*pow(xa, 2))*pow((1-xa), 4)*pow(xa, strange3)*exp(1.0+strange4*xa);
  Double_t Funcsb = (strange0+strange1*xb+strange2*pow(xb, 2))*pow((1-xb), 4)*pow(xb, strange3)*exp(1.0+strange4*xb);
  Double_t Funcba = (bottom0+bottom1*xa+bottom2*pow(xa, 2))*pow((1-xa), 4)*pow(xa, bottom3)*exp(1.0+bottom4*xa);
  Double_t Funcbb = (bottom0+bottom1*xb+bottom2*pow(xb, 2))*pow((1-xb), 4)*pow(xb, bottom3)*exp(1.0+bottom4*xb);

  Double_t totSec = 0;
  
  if (Vdecay1!=RooSpin::kVdecayType_Wany){ // ZH or gammaH
    Double_t FuncABu = FuncAu1/xa*FuncBu1/xb+FuncAu2/xa*FuncBu2/xb;
    Double_t FuncABd = FuncAd1/xa*FuncBd1/xb+FuncAd2/xa*FuncBd2/xb;
    Double_t FuncABc = Funcsa*Funcsb/xa/xb;
    Double_t FuncABs = Funcca*Funccb/xa/xb;
    Double_t FuncABb = Funcba*Funcbb/xa/xb;

    totSec = 2.*prefactor*mVal*(
      (FuncABu)*weightu
      +(FuncABd)*weightd
      +(FuncABc)*weightc
      +(FuncABs)*weights
      +(FuncABb)*weightb
      );
  }
  else{ // WH
    Double_t FuncAB_udbar = FuncAu1/xa*FuncBd1/xb + FuncAu2/xa*FuncBd2/xb;
    Double_t FuncAB_csbar = Funcca*Funcsb/xa/xb;
    Double_t FuncAB_dubar = FuncAd1/xa*FuncBu1/xb + FuncAd2/xa*FuncBu2/xb;
    Double_t FuncAB_scbar = Funccb*Funcsa/xa/xb;

    totSec = 2.*prefactor*mVal*(
      (FuncAB_udbar)*weightu
      +(FuncAB_dubar)*weightd
      +(FuncAB_csbar)*weightc
      +(FuncAB_scbar)*weights
      );
  }

  if (totSec<=0) totSec = 1e-5;
  return totSec;
}



