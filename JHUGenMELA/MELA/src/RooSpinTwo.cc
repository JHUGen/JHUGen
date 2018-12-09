#include "RooSpinTwo.h"

using namespace std;
using namespace MELAStreamHelpers;


RooSpinTwo::RooSpinTwo() : RooSpin(){}
RooSpinTwo::RooSpinTwo(
  const char* name, const char* title,
  modelMeasurables const& _measurables,
  modelParameters const& _parameters,
  modelCouplings const& _couplings,
  RooSpin::VdecayType _Vdecay1, RooSpin::VdecayType _Vdecay2,
  TVar::VerbosityLevel verbosity_
) : RooSpin(
  name, title,
  _measurables, _parameters,
  _Vdecay1, _Vdecay2,
  verbosity_
),

b1Val("b1Val", "b1Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_1][0])),
b2Val("b2Val", "b2Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_2][0])),
b3Val("b3Val", "b3Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_3][0])),
b4Val("b4Val", "b4Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_4][0])),
b5Val("b5Val", "b5Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_5][0])),
b6Val("b6Val", "b6Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_6][0])),
b7Val("b7Val", "b7Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_7][0])),
b8Val("b8Val", "b8Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_8][0])),
b9Val("b9Val", "b9Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_9][0])),
b10Val("b10Val", "b10Val", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_10][0])),

b1ValIm("b1ValIm", "b1ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_1][1])),
b2ValIm("b2ValIm", "b2ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_2][1])),
b3ValIm("b3ValIm", "b3ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_3][1])),
b4ValIm("b4ValIm", "b4ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_4][1])),
b5ValIm("b5ValIm", "b5ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_5][1])),
b6ValIm("b6ValIm", "b6ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_6][1])),
b7ValIm("b7ValIm", "b7ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_7][1])),
b8ValIm("b8ValIm", "b8ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_8][1])),
b9ValIm("b9ValIm", "b9ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_9][1])),
b10ValIm("b10ValIm", "b10ValIm", this, (RooAbsReal&)*(_couplings.bList[gGRAVITON_VV_10][1])),

Lambda("Lambda", "Lambda", this, (RooAbsReal&)*(_couplings.Lambda)),

f_spinz1("f_spinz1", "f_spinz1", this, (RooAbsReal&)*(_couplings.f_spinz1)),
f_spinz2("f_spinz2", "f_spinz2", this, (RooAbsReal&)*(_couplings.f_spinz2))
{}

RooSpinTwo::RooSpinTwo(const RooSpinTwo& other, const char* name) :
RooSpin(other, name),

b1Val("b1Val", this, other.b1Val),
b2Val("a2Val", this, other.b2Val),
b3Val("b3Val", this, other.b3Val),
b4Val("b4Val", this, other.b4Val),
b5Val("b5Val", this, other.b5Val),
b6Val("b6Val", this, other.b6Val),
b7Val("b7Val", this, other.b7Val),
b8Val("b8Val", this, other.b8Val),
b9Val("b9Val", this, other.b9Val),
b10Val("b10Val", this, other.b10Val),

b1ValIm("b1ValIm", this, other.b1ValIm),
b2ValIm("a2ValIm", this, other.b2ValIm),
b3ValIm("b3ValIm", this, other.b3ValIm),
b4ValIm("b4ValIm", this, other.b4ValIm),
b5ValIm("b5ValIm", this, other.b5ValIm),
b6ValIm("b6ValIm", this, other.b6ValIm),
b7ValIm("b7ValIm", this, other.b7ValIm),
b8ValIm("b8ValIm", this, other.b8ValIm),
b9ValIm("b9ValIm", this, other.b9ValIm),
b10ValIm("b10ValIm", this, other.b10ValIm),

Lambda("Lambda", this, other.Lambda),

f_spinz1("f_spinz1", this, other.f_spinz1),
f_spinz2("f_spinz2", this, other.f_spinz2)
{}

void RooSpinTwo::calculateCi(std::vector<Double_t>& ciRe, std::vector<Double_t>& ciIm, bool isGammaV1, bool isGammaV2) const{
  Double_t mV;
  getMVGamV(&mV);

  Double_t m1_=m1; if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) m1_=0;
  Double_t m2_=m2; if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) m2_=0;
  Double_t m1sq = pow(m1_, 2);
  Double_t m2sq = pow(m2_, 2);
  Double_t mVsq = pow(mV, 2);
  Double_t m12sq = pow(m12, 2);

  Double_t s = (m12sq-m1sq-m2sq)/2.;
  if (m1_>m2_+m12 || m2_>m1_+m12) s = -s;
  Double_t kappa = s/pow(Lambda, 2);

  if (!isGammaV1 && !isGammaV2 && !(Vdecay1==RooSpin::kVdecayType_GammaOnshell || Vdecay2==RooSpin::kVdecayType_GammaOnshell)){ // ZZ/WW
    Double_t c1Re = 2.*(b1Val + b2Val*kappa*(1.+m1sq/s)*(1.+m2sq/s) + b5Val*mVsq/s); ciRe.push_back(c1Re);
    Double_t c2Re = -0.5*b1Val + b3Val*kappa*(1.-(m1sq+m2sq)/(2*s)) + 2.*b4Val*kappa + b7Val*kappa*mVsq/s; ciRe.push_back(c2Re);
    Double_t c3Re = -(b2Val/2.+b3Val+2.*b4Val)*kappa*m12sq/s; ciRe.push_back(c3Re);
    Double_t c41Re = -b1Val - b2Val*kappa - (b2Val*m1sq+b3Val*m2sq+2.*b6Val*mVsq)*kappa/s; ciRe.push_back(c41Re);
    Double_t c42Re = -b1Val - b2Val*kappa - (b2Val*m2sq+b3Val*m1sq+2.*b6Val*mVsq)*kappa/s; ciRe.push_back(c42Re);
    Double_t c5Re = 2.*b8Val*kappa*(m12sq)/s; ciRe.push_back(c5Re);
    Double_t c6Re = b9Val*kappa*mVsq/s; ciRe.push_back(c6Re);
    Double_t c7Re = b10Val*m12sq*mVsq*pow(kappa/s, 2); ciRe.push_back(c7Re);

    Double_t c1Im = 2.*(b1ValIm + b2ValIm*kappa*(1.+m1sq/s)*(1.+m2sq/s) + b5ValIm*mVsq/s); ciIm.push_back(c1Im);
    Double_t c2Im = -0.5*b1ValIm + b3ValIm*kappa*(1.-(m1sq+m2sq)/(2*s)) + 2.*b4ValIm*kappa + b7ValIm*kappa*mVsq/s; ciIm.push_back(c2Im);
    Double_t c3Im = -(b2ValIm/2.+b3ValIm+2.*b4ValIm)*kappa*m12sq/s; ciIm.push_back(c3Im);
    Double_t c41Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m1sq+b3ValIm*m2sq+2.*b6ValIm*mVsq)*kappa/s; ciIm.push_back(c41Im);
    Double_t c42Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m2sq+b3ValIm*m1sq+2.*b6ValIm*mVsq)*kappa/s; ciIm.push_back(c42Im);
    Double_t c5Im = 2.*b8ValIm*kappa*(m12sq)/s; ciIm.push_back(c5Im);
    Double_t c6Im = b9ValIm*kappa*mVsq/s; ciIm.push_back(c6Im);
    Double_t c7Im = b10ValIm*m12sq*mVsq*pow(kappa/s, 2); ciIm.push_back(c7Im);
  }
  //else if ((!isGammaV1 || !isGammaV2) && !(Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell)){ // ZGs/ZG
  //???
  //}
  //else{ // GG/GGs
  else{ // Z/G/Gs - G/Gs
    Double_t c1Re = 2.*(b1Val + b2Val*kappa*(1.+m1sq/s)*(1.+m2sq/s)); ciRe.push_back(c1Re);
    Double_t c2Re = -0.5*b1Val + b3Val*kappa*(1.-(m1sq+m2sq)/(2*s)) + 2.*b4Val*kappa; ciRe.push_back(c2Re);
    Double_t c3Re = -(b2Val/2.+b3Val+2.*b4Val)*kappa*m12sq/s; ciRe.push_back(c3Re);
    Double_t c41Re = -b1Val - b2Val*kappa - (b2Val*m1sq+b3Val*m2sq)*kappa/s; ciRe.push_back(c41Re);
    Double_t c42Re = -b1Val - b2Val*kappa - (b2Val*m2sq+b3Val*m1sq)*kappa/s; ciRe.push_back(c42Re);
    Double_t c5Re = 2.*b8Val*kappa*(m12sq)/s; ciRe.push_back(c5Re);
    Double_t c6Re = 0; ciRe.push_back(c6Re);
    Double_t c7Re = 0; ciRe.push_back(c7Re);

    Double_t c1Im = 2.*(b1ValIm + b2ValIm*kappa*(1.+m1sq/s)*(1.+m2sq/s)); ciIm.push_back(c1Im);
    Double_t c2Im = -0.5*b1ValIm + b3ValIm*kappa*(1.-(m1sq+m2sq)/(2*s)) + 2.*b4ValIm*kappa; ciIm.push_back(c2Im);
    Double_t c3Im = -(b2ValIm/2.+b3ValIm+2.*b4ValIm)*kappa*m12sq/s; ciIm.push_back(c3Im);
    Double_t c41Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m1sq+b3ValIm*m2sq)*kappa/s; ciIm.push_back(c41Im);
    Double_t c42Im = -b1ValIm - b2ValIm*kappa - (b2ValIm*m2sq+b3ValIm*m1sq)*kappa/s; ciIm.push_back(c42Im);
    Double_t c5Im = 2.*b8ValIm*kappa*(m12sq)/s; ciIm.push_back(c5Im);
    Double_t c6Im = 0; ciIm.push_back(c6Im);
    Double_t c7Im = 0; ciIm.push_back(c7Im);
  }
}
void RooSpinTwo::calculateAmplitudes(
  Double_t& A00Re, Double_t& A00Im,
  Double_t& AppRe, Double_t& AppIm, Double_t& A0pRe, Double_t& A0pIm, Double_t& Ap0Re, Double_t& Ap0Im,
  Double_t& AmmRe, Double_t& AmmIm, Double_t& A0mRe, Double_t& A0mIm, Double_t& Am0Re, Double_t& Am0Im,
  Double_t& ApmRe, Double_t& ApmIm, Double_t& AmpRe, Double_t& AmpIm,
  bool isGammaV1, bool isGammaV2
  )const{
  Double_t m1_=m1; if (Vdecay1==RooSpin::kVdecayType_GammaOnshell) m1_=0;
  Double_t m2_=m2; if (Vdecay2==RooSpin::kVdecayType_GammaOnshell) m2_=0;

  std::vector<Double_t> ciRe;
  std::vector<Double_t> ciIm;
  calculateCi(ciRe, ciIm, isGammaV1, isGammaV2);

  Double_t propV1Re=1, propV2Re=1, propHRe=1;
  Double_t propV1Im=0, propV2Im=0, propHIm=0;
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell) calculatePropagator(propV1Re, propV1Im, m1_, (isGammaV1 ? 0 : 1));
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) calculatePropagator(propV2Re, propV2Im, m2_, (isGammaV2 ? 0 : 1));
  calculatePropagator(propHRe, propHIm, m12, 2);

  Double_t ampScale = calculateAmplitudeScale((int) isGammaV1, (int) isGammaV2)*pow(GeVunit/Lambda, 2);

  Double_t c1Re = ciRe.at(0);
  Double_t c2Re = ciRe.at(1);
  Double_t c3Re = ciRe.at(2);
  Double_t c41Re = ciRe.at(3);
  Double_t c42Re = ciRe.at(4);
  Double_t c5Re = ciRe.at(5);
  Double_t c6Re = ciRe.at(6);
  Double_t c7Re = ciRe.at(7);
  Double_t c1Im = ciIm.at(0);
  Double_t c2Im = ciIm.at(1);
  Double_t c3Im = ciIm.at(2);
  Double_t c41Im = ciIm.at(3);
  Double_t c42Im = ciIm.at(4);
  Double_t c5Im = ciIm.at(5);
  Double_t c6Im = ciIm.at(6);
  Double_t c7Im = ciIm.at(7);

  Double_t eta1 = m1_ / m12;
  Double_t eta2 = m2_ / m12;
  Double_t eta1p2 = eta1*eta2;

  Double_t eta1sq = eta1*eta1;
  Double_t eta2sq = eta2*eta2;
  Double_t eta1p2sq = pow(eta1p2, 2);

  Double_t etas = (1. - eta1sq - eta2sq)/2.;
  if (pow(eta1+eta2, 2)>1.) etas = -etas;
  Double_t x = etas;
  Double_t xsq = x*x;
  Double_t xxp = (pow(etas, 2)-eta1p2sq);
  if (xxp<0) xxp=0;

  Double_t A00Re_tmp=0, A00Im_tmp=0,
    AppRe_tmp=0, AppIm_tmp=0, A0pRe_tmp=0, A0pIm_tmp=0, Ap0Re_tmp=0, Ap0Im_tmp=0,
    AmmRe_tmp=0, AmmIm_tmp=0, A0mRe_tmp=0, A0mIm_tmp=0, Am0Re_tmp=0, Am0Im_tmp=0,
    ApmRe_tmp=0, ApmIm_tmp=0, AmpRe_tmp=0, AmpIm_tmp=0;

  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell){
    A00Re_tmp =
      pow(m12, 4)*sqrt(2./3.)*
      (
      c1Re*(
      eta1p2sq * (xsq - eta1p2sq/4.)
      - (pow(eta1, 4)+pow(eta2, 4)) * xsq/2.
      + (pow(eta1, 8)+pow(eta2, 8))/8.
      + xsq/2.
      - (pow(eta1, 4)+pow(eta2, 4))/4.
      + 1.0/8.
      )
      + c2Re*2.*xxp*(
      (pow(eta1, 4) + pow(eta2, 4))
      - 2.*(eta1p2sq + 2.*xxp)
      - 1.
      )
      -c3Re*8.*pow(xxp, 2)
      + c41Re*2.*xxp*(1.+eta1sq-eta2sq)
      + c42Re*2.*xxp*(1.-eta1sq+eta2sq)
      );

    A00Im_tmp =
      pow(m12, 4)*sqrt(2./3.)*
      (
      c1Im*(
      eta1p2sq * (xsq - eta1p2sq/4.)
      - (pow(eta1, 4)+pow(eta2, 4)) * xsq/2.
      + (pow(eta1, 8)+pow(eta2, 8))/8.
      + xsq/2.
      - (pow(eta1, 4)+pow(eta2, 4))/4.
      + 1.0/8.
      )
      + c2Im*2.*xxp*(
      (pow(eta1, 4) + pow(eta2, 4))
      - 2.*(eta1p2sq + 2.*xxp)
      - 1.
      )
      -c3Im*8.*pow(xxp, 2)
      + c41Im*2.*xxp*(1.+eta1sq-eta2sq)
      + c42Im*2.*xxp*(1.-eta1sq+eta2sq)
      );
  }

  //-----------------------------------------------------------------------
  // No m1_/m2_ singularities in A--

  AmmRe_tmp =
    pow(m12, 4)/sqrt(6.)*
    (
    c1Re/4.*(
    1. + 4.*xxp - pow(eta1sq-eta2sq, 2)
    )
    + c2Re*8.*xxp
    + c5Im*8.*pow(xxp, 1.5)
    - c6Im*4.*sqrt(xxp)
    );

  AmmIm_tmp =
    pow(m12, 4)/sqrt(6.)*
    (
    c1Im/4.*(
    1. + 4.*xxp - pow(eta1sq-eta2sq, 2)
    )
    + c2Im*8.*xxp
    - c5Re*8.*pow(xxp, 1.5)
    + c6Re*4.*sqrt(xxp)
    );

  //-----------------------------------------------------------------------
  // No m1_/m2_ singularities in A++

  AppRe_tmp =
    pow(m12, 4)/sqrt(6.)*
    (
    c1Re/4.*(
    1. + 4.*xxp - pow(eta1sq-eta2sq, 2)
    )
    + c2Re*8.*xxp
    - c5Im*8.*pow(xxp, 1.5)
    + c6Im*4.*sqrt(xxp)
    );

  AppIm_tmp =
    pow(m12, 4)/sqrt(6.)*
    (
    c1Im/4.*(
    1. + 4.*xxp - pow(eta1sq-eta2sq, 2)
    )
    + c2Im*8.*xxp
    + c5Re*8.*pow(xxp, 1.5)
    - c6Re*4.*sqrt(xxp)
    );

  //-----------------------------------------------------------------------
  // No m1_/m2_ singularities in A+- or A-+

  AmpRe_tmp = pow(m12, 4)*(c1Re/4.*(1.+4.*xxp-pow(eta1sq-eta2sq, 2)));
  AmpIm_tmp = pow(m12, 4)*(c1Im/4.*(1.+4.*xxp-pow(eta1sq-eta2sq, 2)));
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell) { AmmIm_tmp *= eta1; AmmRe_tmp *= eta1; AppIm_tmp *= eta1; AppRe_tmp *= eta1; AmpIm_tmp *= eta1; AmpRe_tmp *= eta1; } // Do these multiplications here...
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) { AmmIm_tmp *= eta2; AmmRe_tmp *= eta2; AppIm_tmp *= eta2; AppRe_tmp *= eta2; AmpIm_tmp *= eta2; AmpRe_tmp *= eta2; }
  ApmRe_tmp = AmpRe_tmp; // ...for this reason!
  ApmIm_tmp = AmpIm_tmp;

  //-----------------------------------------------------------------------

  Double_t A0m_0p_m0_p0_c4factor = 2.*xxp;
  Double_t A0m_0p_m0_p0_c7factor = 4.*pow(xxp, 1.5); // x+-i

  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell){
    Double_t A0m_0p_c1factor = (
      -(pow(eta1, 6)-pow(eta2, 6))/8.
      + (eta1sq-eta2sq)*(3.*eta1p2sq + 4.*xxp)/8.
      - pow(eta1sq-eta2sq, 2)/8.
      + xxp/2.
      + (1. + (eta1sq-eta2sq))/8.
      );
    Double_t A0m_0p_c6factor = sqrt(xxp)*(1.+eta1sq-eta2sq); // x+-i

    A0mRe_tmp =
      pow(m12, 4)*
      (
      c1Re*A0m_0p_c1factor
      + c42Re*A0m_0p_m0_p0_c4factor
      - c6Im*A0m_0p_c6factor
      - c7Im*A0m_0p_m0_p0_c7factor
      );

    A0mIm_tmp =
      pow(m12, 4)*
      (
      c1Im*A0m_0p_c1factor
      + c42Im*A0m_0p_m0_p0_c4factor
      + c6Re*A0m_0p_c6factor
      + c7Re*A0m_0p_m0_p0_c7factor
      );

    //-----------------------------------------------------------------------

    A0pRe_tmp =
      pow(m12, 4)*
      (
      c1Re*A0m_0p_c1factor
      + c42Re*A0m_0p_m0_p0_c4factor
      + c6Im*A0m_0p_c6factor
      + c7Im*A0m_0p_m0_p0_c7factor
      );

    A0pIm_tmp =
      pow(m12, 4)*
      (
      c1Im*A0m_0p_c1factor
      + c42Im*A0m_0p_m0_p0_c4factor
      - c6Re*A0m_0p_c6factor
      - c7Re*A0m_0p_m0_p0_c7factor
      );

    if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell){
      A0mIm_tmp *= eta2; A0mRe_tmp *= eta2;
      A0pIm_tmp *= eta2; A0pRe_tmp *= eta2;
    }
  }

  //-----------------------------------------------------------------------

  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell){
    Double_t Am0_p0_c1factor = (
      (pow(eta1, 6)-pow(eta2, 6))/8.
      - (eta1sq-eta2sq)*(3.*eta1p2sq + 4.*xxp)/8.
      - pow(eta1sq-eta2sq, 2)/8.
      + xxp/2.
      + (1. - (eta1sq-eta2sq))/8.
      );
    Double_t Am0_p0_c6factor = sqrt(xxp)*(1.-eta1sq+eta2sq); // x+-i

    Am0Re_tmp =
      pow(m12, 4)*
      (
      c1Re*Am0_p0_c1factor
      + c41Re*A0m_0p_m0_p0_c4factor
      - c6Im*Am0_p0_c6factor
      - c7Im*A0m_0p_m0_p0_c7factor
      );

    Am0Im_tmp =
      pow(m12, 4)*
      (
      c1Im*Am0_p0_c1factor
      + c41Im*A0m_0p_m0_p0_c4factor
      + c6Re*Am0_p0_c6factor
      + c7Re*A0m_0p_m0_p0_c7factor
      );

    //-----------------------------------------------------------------------

    Ap0Re_tmp =
      pow(m12, 4)*
      (
      c1Re*Am0_p0_c1factor
      + c41Re*A0m_0p_m0_p0_c4factor
      + c6Im*Am0_p0_c6factor
      + c7Im*A0m_0p_m0_p0_c7factor
      );

    Ap0Im_tmp =
      pow(m12, 4)*
      (
      c1Im*Am0_p0_c1factor
      + c41Im*A0m_0p_m0_p0_c4factor
      - c6Re*Am0_p0_c6factor
      - c7Re*A0m_0p_m0_p0_c7factor
      );

    if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell){
      Am0Im_tmp *= eta1; Am0Re_tmp *= eta1;
      Ap0Im_tmp *= eta1; Ap0Re_tmp *= eta1;
    }
  }

  //-----------------------------------------------------------------------

  A00Re_tmp *= ampScale;
  AmmRe_tmp *= ampScale;
  AppRe_tmp *= ampScale;
  A0mRe_tmp *= ampScale;
  A0pRe_tmp *= ampScale;
  Am0Re_tmp *= ampScale;
  Ap0Re_tmp *= ampScale;
  AmpRe_tmp *= ampScale;
  ApmRe_tmp *= ampScale;
  A00Im_tmp *= ampScale;
  AmmIm_tmp *= ampScale;
  AppIm_tmp *= ampScale;
  A0mIm_tmp *= ampScale;
  A0pIm_tmp *= ampScale;
  Am0Im_tmp *= ampScale;
  Ap0Im_tmp *= ampScale;
  AmpIm_tmp *= ampScale;
  ApmIm_tmp *= ampScale;

  std::vector<Double_t> A00_reals, A00_imags; A00_reals.push_back(A00Re_tmp); A00_imags.push_back(A00Im_tmp); A00_reals.push_back(propV1Re); A00_imags.push_back(propV1Im); A00_reals.push_back(propV2Re); A00_imags.push_back(propV2Im); A00_reals.push_back(propHRe); A00_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(A00_reals, A00_imags, A00Re, A00Im);
  std::vector<Double_t> Amm_reals, Amm_imags; Amm_reals.push_back(AmmRe_tmp); Amm_imags.push_back(AmmIm_tmp); Amm_reals.push_back(propV1Re); Amm_imags.push_back(propV1Im); Amm_reals.push_back(propV2Re); Amm_imags.push_back(propV2Im); Amm_reals.push_back(propHRe); Amm_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(Amm_reals, Amm_imags, AmmRe, AmmIm);
  std::vector<Double_t> App_reals, App_imags; App_reals.push_back(AppRe_tmp); App_imags.push_back(AppIm_tmp); App_reals.push_back(propV1Re); App_imags.push_back(propV1Im); App_reals.push_back(propV2Re); App_imags.push_back(propV2Im); App_reals.push_back(propHRe); App_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(App_reals, App_imags, AppRe, AppIm);
  std::vector<Double_t> A0m_reals, A0m_imags; A0m_reals.push_back(A0mRe_tmp); A0m_imags.push_back(A0mIm_tmp); A0m_reals.push_back(propV1Re); A0m_imags.push_back(propV1Im); A0m_reals.push_back(propV2Re); A0m_imags.push_back(propV2Im); A0m_reals.push_back(propHRe); A0m_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(A0m_reals, A0m_imags, A0mRe, A0mIm);
  std::vector<Double_t> A0p_reals, A0p_imags; A0p_reals.push_back(A0pRe_tmp); A0p_imags.push_back(A0pIm_tmp); A0p_reals.push_back(propV1Re); A0p_imags.push_back(propV1Im); A0p_reals.push_back(propV2Re); A0p_imags.push_back(propV2Im); A0p_reals.push_back(propHRe); A0p_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(A0p_reals, A0p_imags, A0pRe, A0pIm);
  std::vector<Double_t> Am0_reals, Am0_imags; Am0_reals.push_back(Am0Re_tmp); Am0_imags.push_back(Am0Im_tmp); Am0_reals.push_back(propV1Re); Am0_imags.push_back(propV1Im); Am0_reals.push_back(propV2Re); Am0_imags.push_back(propV2Im); Am0_reals.push_back(propHRe); Am0_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(Am0_reals, Am0_imags, Am0Re, Am0Im);
  std::vector<Double_t> Ap0_reals, Ap0_imags; Ap0_reals.push_back(Ap0Re_tmp); Ap0_imags.push_back(Ap0Im_tmp); Ap0_reals.push_back(propV1Re); Ap0_imags.push_back(propV1Im); Ap0_reals.push_back(propV2Re); Ap0_imags.push_back(propV2Im); Ap0_reals.push_back(propHRe); Ap0_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(Ap0_reals, Ap0_imags, Ap0Re, Ap0Im);
  std::vector<Double_t> Amp_reals, Amp_imags; Amp_reals.push_back(AmpRe_tmp); Amp_imags.push_back(AmpIm_tmp); Amp_reals.push_back(propV1Re); Amp_imags.push_back(propV1Im); Amp_reals.push_back(propV2Re); Amp_imags.push_back(propV2Im); Amp_reals.push_back(propHRe); Amp_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(Amp_reals, Amp_imags, AmpRe, AmpIm);
  std::vector<Double_t> Apm_reals, Apm_imags; Apm_reals.push_back(ApmRe_tmp); Apm_imags.push_back(ApmIm_tmp); Apm_reals.push_back(propV1Re); Apm_imags.push_back(propV1Im); Apm_reals.push_back(propV2Re); Apm_imags.push_back(propV2Im); Apm_reals.push_back(propHRe); Apm_imags.push_back(propHIm); AnaMelaHelpers::multiplyComplexNumbers(Apm_reals, Apm_imags, ApmRe, ApmIm);


  if (
    A00Re!=A00Re || A00Im!=A00Im ||
    AppRe!=AppRe || AppIm!=AppIm ||
    AmmRe!=AmmRe || AmmIm!=AmmIm ||
    A0pRe!=A0pRe || A0pIm!=A0pIm ||
    A0mRe!=A0mRe || A0mIm!=A0mIm ||
    Ap0Re!=Ap0Re || Ap0Im!=Ap0Im ||
    Am0Re!=Am0Re || Am0Im!=Am0Im ||
    ApmRe!=ApmRe || ApmIm!=ApmIm ||
    AmpRe!=AmpRe || AmpIm!=AmpIm ||
    (
    A00Re==0 && A00Im==0 &&
    AppRe==0 && AppIm==0 &&
    AmmRe==0 && AmmIm==0 &&
    A0pRe==0 && A0pIm==0 &&
    A0mRe==0 && A0mIm==0 &&
    Ap0Re==0 && Ap0Im==0 &&
    Am0Re==0 && Am0Im==0 &&
    ApmRe==0 && ApmIm==0 &&
    AmpRe==0 && AmpIm==0
    )
    ){
    MELAerr << "Some of the amplitudes are NaN or all are 0:" << endl;
    MELAerr << "A00Re=" << A00Re << ", A00Im=" << A00Im << endl;
    MELAerr << "AppRe=" << AppRe << ", AppIm=" << AppIm << endl;
    MELAerr << "AmmRe=" << AmmRe << ", AmmIm=" << AmmIm << endl;
    MELAerr << "A0pRe=" << A0pRe << ", A0pIm=" << A0pIm << endl;
    MELAerr << "A0mRe=" << A0mRe << ", A0mIm=" << A0mIm << endl;
    MELAerr << "Ap0Re=" << Ap0Re << ", Ap0Im=" << Ap0Im << endl;
    MELAerr << "Am0Re=" << Am0Re << ", Am0Im=" << Am0Im << endl;
    MELAerr << "ApmRe=" << ApmRe << ", ApmIm=" << ApmIm << endl;
    MELAerr << "AmpRe=" << AmpRe << ", AmpIm=" << AmpIm << endl;

    MELAerr << "Possible causes are" << endl;
    MELAerr << "m12=" << m12 << endl;
    MELAerr << "m1=" << m1_ << endl;
    MELAerr << "m2=" << m2_ << endl;
    MELAerr << "x=" << x << endl;
    MELAerr << "xxp=" << xxp << endl;

    MELAerr << "c1 = (" << c1Re << ", " << c1Im << ")" << endl;
    MELAerr << "c2 = (" << c2Re << ", " << c2Im << ")" << endl;
    MELAerr << "c3 = (" << c3Re << ", " << c3Im << ")" << endl;
    MELAerr << "c41 = (" << c41Re << ", " << c41Im << ")" << endl;
    MELAerr << "c42 = (" << c42Re << ", " << c42Im << ")" << endl;
    MELAerr << "c5 = (" << c5Re << ", " << c5Im << ")" << endl;
    MELAerr << "c6 = (" << c6Re << ", " << c6Im << ")" << endl;
    MELAerr << "c7 = (" << c7Re << ", " << c7Im << ")" << endl;
  }

  return;
}

void RooSpinTwo::printParameters() const{
  MELAout << "b1: (" << b1Val << ", " << b1ValIm << ")" << endl;
  MELAout << "b2: (" << b2Val << ", " << b2ValIm << ")" << endl;
  MELAout << "b3: (" << b3Val << ", " << b3ValIm << ")" << endl;
  MELAout << "b4: (" << b4Val << ", " << b4ValIm << ")" << endl;
  MELAout << "b5: (" << b5Val << ", " << b5ValIm << ")" << endl;
  MELAout << "b6: (" << b6Val << ", " << b6ValIm << ")" << endl;
  MELAout << "b7: (" << b7Val << ", " << b7ValIm << ")" << endl;
  MELAout << "b8: (" << b8Val << ", " << b8ValIm << ")" << endl;
  MELAout << "b9: (" << b9Val << ", " << b9ValIm << ")" << endl;
  MELAout << "b10: (" << b10Val << ", " << b10ValIm << ")" << endl;
  MELAout << "Lambda: " << Lambda << endl;
  MELAout << "f_spinz1: " << f_spinz1 << endl;
  MELAout << "f_spinz2: " << f_spinz2 << endl;

  RooSpin::printParameters();
}
