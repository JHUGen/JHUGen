#include "RooSpinTwo_7DComplex_ppHVV.h"


using namespace std;
using namespace MELAStreamHelpers;


RooSpinTwo_7DComplex_ppHVV::RooSpinTwo_7DComplex_ppHVV() : RooSpinTwo(){}
RooSpinTwo_7DComplex_ppHVV::RooSpinTwo_7DComplex_ppHVV(
  const char *name, const char *title,
  modelMeasurables const& _measurables,
  modelParameters const& _parameters,
  modelCouplings const& _couplings,
  RooSpin::VdecayType _Vdecay1, RooSpin::VdecayType _Vdecay2,
  TVar::VerbosityLevel verbosity_
) : RooSpinTwo(
  name, title,
  _measurables,
  _parameters,
  _couplings,
  _Vdecay1, _Vdecay2,
  verbosity_
),
ZZ4fOrdering(true)
{}


RooSpinTwo_7DComplex_ppHVV::RooSpinTwo_7DComplex_ppHVV(
  const RooSpinTwo_7DComplex_ppHVV& other, const char* name
) : RooSpinTwo(other, name),
ZZ4fOrdering(other.ZZ4fOrdering)
{}


Double_t RooSpinTwo_7DComplex_ppHVV::evaluateH1Factor(Int_t i1, Int_t j1, Int_t helicity, Int_t code) const{
  const Double_t Pi = TMath::Pi();
  Double_t dHel = (Double_t)helicity;
  Double_t result = 0;
  if ((code % prime_h1)==0){
    if ((i1==0 && j1==0) || (i1==-1 && j1==1) || (i1==1 && j1==-1)) result = 4./3.; // 15 amps
    else if (i1==1 && j1==1) result = 8./3.; // 6 amps
    else if (i1==-1 && j1==-1) result = 8./3.; // 6 amps
    else if ((i1==0 && j1==1) || (i1==1 && j1==0)) result = Pi*dHel/2.; // 9 amps
    else if ((i1==0 && j1==-1) || (i1==-1 && j1==0)) result = Pi*dHel/2.; // 9 amps
  }
  else{
    if ((i1==0 && j1==0) || (i1==-1 && j1==1) || (i1==1 && j1==-1)) result = 1.-pow(h1, 2); // 15 amps
    else if (i1==1 && j1==1) result = 1.+pow(h1, 2)-2.*h1*dHel; // 6 amps
    else if (i1==-1 && j1==-1) result = 1.+pow(h1, 2)+2.*h1*dHel; // 6 amps
    else if ((i1==0 && j1==1) || (i1==1 && j1==0)) result = sqrt(fabs(1.-pow(h1, 2)))*(dHel-h1); // 9 amps
    else if ((i1==0 && j1==-1) || (i1==-1 && j1==0)) result = sqrt(fabs(1.-pow(h1, 2)))*(dHel+h1); // 9 amps
  }
  return result;
}
Double_t RooSpinTwo_7DComplex_ppHVV::evaluateH2Factor(Int_t i2, Int_t j2, Int_t helicity, Int_t code) const{
  const Double_t Pi = TMath::Pi();
  Double_t dHel = (Double_t)helicity;
  Double_t result = 0;
  if ((code % prime_h2)==0){
    if ((i2==0 && j2==0) || (i2==-1 && j2==1) || (i2==1 && j2==-1)) result = 4./3.; // 15 amps
    else if (i2==1 && j2==1) result = 8./3.; // 6 amps
    else if (i2==-1 && j2==-1) result = 8./3.; // 6 amps
    else if ((i2==0 && j2==1) || (i2==1 && j2==0)) result = Pi*dHel/2.; // 9 amps
    else if ((i2==0 && j2==-1) || (i2==-1 && j2==0)) result = Pi*dHel/2.; // 9 amps
  }
  else{
    if ((i2==0 && j2==0) || (i2==-1 && j2==1) || (i2==1 && j2==-1)) result = 1.-pow(h2, 2); // 15 amps
    else if (i2==1 && j2==1) result = 1.+pow(h2, 2)-2.*h2*dHel; // 6 amps
    else if (i2==-1 && j2==-1) result = 1.+pow(h2, 2)+2.*h2*dHel; // 6 amps
    else if ((i2==0 && j2==1) || (i2==1 && j2==0)) result = sqrt(fabs(1.-pow(h2, 2)))*(dHel-h2); // 9 amps
    else if ((i2==0 && j2==-1) || (i2==-1 && j2==0)) result = sqrt(fabs(1.-pow(h2, 2)))*(dHel+h2); // 9 amps
  }
  return result;
}
Double_t RooSpinTwo_7DComplex_ppHVV::evaluateHSFactor(Int_t di, Int_t dj, Int_t code) const{
  Double_t f_spinz0 = 1. - f_spinz1 - f_spinz2;
  if (f_spinz0<0.) f_spinz0=0;
  Double_t hsneg = -hs; if (fabs(hsneg)>1.) hsneg *= 1./fabs(hsneg);

  Double_t AF200 = 0;
  Double_t AF201 = 0;
  Double_t AF202 = 0;
  Double_t AF211 = 0;
  Double_t AF2m11 = 0;
  Double_t AF212 = 0;
  Double_t AF2m12 = 0;
  Double_t AF222 = 0;
  Double_t AF2m22 = 0;
  if ((code % prime_hs)==0){
    AF200 = 2.*((2.*f_spinz0 + 3.*f_spinz2) - 6.*(2.*f_spinz0 - 2.*f_spinz1 + f_spinz2)/3. +3.*(6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)/5.)/8.; // 6
    //AF201 = 0; // 12
    AF202 = ((2.*f_spinz0 - f_spinz2)*(4./3.) - (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)*(4./15.))*(-sqrt(3./2.)/8.); // 6

    AF211 = 2.*((f_spinz1 + f_spinz2) + 3.*(2.*f_spinz0 - f_spinz1)/3. - (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)/5.)/4.; // 6
    AF2m11 = ((f_spinz1 - f_spinz2)*(4./3.) + (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)*(4./15.))/(-4.); // 4

    //AF212 = 0; // 4
    //AF2m12 = 0; // 4

    AF222 = 2.*((6.*f_spinz0 + 4.*f_spinz1 + f_spinz2) - 6.*(2.*f_spinz0 - f_spinz2)/3. + (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)/5.)/16.; // 2
    AF2m22 = 16./15.*(6*f_spinz0 - 4*f_spinz1 + f_spinz2)/16.; // 1
  }
  else{
    AF200 = ((2.*f_spinz0 + 3.*f_spinz2) - 6.*(2.*f_spinz0 - 2.*f_spinz1 + f_spinz2)*pow(hsneg, 2) +3.*(6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)*pow(hsneg, 4))/8.; // 6 // F(2)00
    AF201 = (hsneg*sqrt(1 - pow(hsneg, 2))*((2.*f_spinz0 - 2.*f_spinz1 + f_spinz2) - (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)*pow(hsneg, 2)))*(-sqrt(6.)/8.); // 12 // F(2)01
    AF202 = ((1. - pow(hsneg, 2))*((2.*f_spinz0 - f_spinz2) - (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)*pow(hsneg, 2)))*(-sqrt(3./2.)/8.); // 6 // F(2)02

    AF211 = ((f_spinz1 + f_spinz2) + 3.*(2.*f_spinz0 - f_spinz1)*pow(hsneg, 2) - (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)*pow(hsneg, 4))/4.; // 6 // F(2)11
    AF2m11 = ((1. - pow(hsneg, 2))*((f_spinz1 - f_spinz2) + (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)*pow(hsneg, 2)))/(-4.); // 4 // F(2)-11

    AF212 = (hsneg*sqrt(1 - pow(hsneg, 2))*((6.*f_spinz0 - 3.*f_spinz2) - (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)*pow(hsneg, 2)))/(-8.); // 4 // F(2)12
    AF2m12 = ((6*f_spinz0 - 4*f_spinz1 + f_spinz2)*hsneg*pow(1. - pow(hsneg, 2), 1.5))/(-8.); // 4 // F(2)-12

    AF222 = ((6.*f_spinz0 + 4.*f_spinz1 + f_spinz2) - 6.*(2.*f_spinz0 - f_spinz2)*pow(hsneg, 2) + (6.*f_spinz0 - 4.*f_spinz1 + f_spinz2)*pow(hsneg, 4))/16.; // 2 // F(2)22
    AF2m22 = ((6*f_spinz0 - 4*f_spinz1 + f_spinz2)*pow(1. - pow(hsneg, 2), 2))/16.; // 1 // F(2)-22
  }

  Double_t result = 0;
  if (di==0 && dj==0) result = AF200;
  else if ((di==1 && dj==1) || (di==-1 && dj==-1)) result = AF211;
  else if ((di==2 && dj==2) || (di==-2 && dj==-2)) result = AF222;

  else if ((di==1 && dj==0) || (di==0 && dj==1)) result = AF201;
  else if ((di==-1 && dj==0) || (di==0 && dj==-1)) result = -AF201;
  else if ((di==1 && dj==2) || (di==2 && dj==1)) result = AF212;
  else if ((di==-1 && dj==-2) || (di==-2 && dj==-1)) result = -AF212;

  else if ((di==2 && dj==0) || (di==0 && dj==2) || (di==-2 && dj==0) || (di==0 && dj==-2)) result = AF202;
  else if ((di==1 && dj==-1) || (di==-1 && dj==1)) result = AF2m11;

  else if ((di==-1 && dj==2) || (di==2 && dj==-1)) result = AF2m12;
  else if ((di==-2 && dj==1) || (di==1 && dj==-2)) result = -AF2m12;

  else if ((di==-2 && dj==2) || (di==2 && dj==-2)) result = AF2m22;
  return result;
}
Double_t RooSpinTwo_7DComplex_ppHVV::evaluatePhi1PhiFactor(Int_t i1, Int_t i2, Int_t j1, Int_t j2, Int_t code, Double_t extraPhase1, Double_t extraPhase2) const{
  const Double_t Pi = TMath::Pi();

  Double_t result = 0;
  Double_t phase = 0;
  Double_t phasePhi = 0;
  Double_t phasePhi1 = 0;

  phasePhi1 += -i1;
  phasePhi1 += i2;
  phasePhi += i2;

  phasePhi1 += j1;
  phasePhi1 += -j2;
  phasePhi += -j2;

  if ((code % prime_Phi)==0 && (code % prime_Phi1)==0){
    if (i1==j1 && i2==j2) result = 4.*pow(Pi, 2);
    // Everything else is 0!
  }
  else if ((code % prime_Phi)==0){
    if (i1==j1 && i2==j2) result = 2.*Pi;
    else if (phasePhi==0.){
      phase = Phi1*phasePhi1+extraPhase1-extraPhase2;
      result = cos(phase)*2.*Pi;
    }
  }
  else if ((code % prime_Phi1)==0){
    if (i1==j1 && i2==j2) result = 2.*Pi;
    else if (phasePhi1==0.){
      phase = Phi*phasePhi+extraPhase1-extraPhase2;
      result = cos(phase)*2.*Pi;
    }
  }
  else{
    phase = Phi1*phasePhi1+Phi*phasePhi+extraPhase1-extraPhase2;
    result = cos(phase);
  }
  return result;
}

void RooSpinTwo_7DComplex_ppHVV::evaluatePolarizationTerms(std::vector<Double_t>& Axxyyterm, const Int_t code, bool isGammaV1, bool isGammaV2) const{
  Double_t R1Val, R2Val;
  calculateVffR1R2(R1Val, R2Val, isGammaV1, isGammaV2);

  Double_t
    A00Re, A00Im,
    AppRe, AppIm,
    A0pRe, A0pIm, Ap0Re, Ap0Im,
    AmmRe, AmmIm, A0mRe, A0mIm, Am0Re, Am0Im,
    ApmRe, ApmIm, AmpRe, AmpIm;

  calculateAmplitudes(
    A00Re, A00Im,
    AppRe, AppIm,
    A0pRe, A0pIm, Ap0Re, Ap0Im,
    AmmRe, AmmIm, A0mRe, A0mIm, Am0Re, Am0Im,
    ApmRe, ApmIm, AmpRe, AmpIm,
    isGammaV1, isGammaV2
    );

  std::vector<Double_t> ARexy, AImxy;
  std::vector<Int_t> index_x,index_y;
  ARexy.push_back(AmmRe);
  ARexy.push_back(Am0Re);
  ARexy.push_back(AmpRe);
  ARexy.push_back(A0mRe);
  ARexy.push_back(A00Re);
  ARexy.push_back(A0pRe);
  ARexy.push_back(ApmRe);
  ARexy.push_back(Ap0Re);
  ARexy.push_back(AppRe);
  AImxy.push_back(AmmIm);
  AImxy.push_back(Am0Im);
  AImxy.push_back(AmpIm);
  AImxy.push_back(A0mIm);
  AImxy.push_back(A00Im);
  AImxy.push_back(A0pIm);
  AImxy.push_back(ApmIm);
  AImxy.push_back(Ap0Im);
  AImxy.push_back(AppIm);
  index_x.push_back(-1);
  index_x.push_back(-1);
  index_x.push_back(-1);
  index_x.push_back(0);
  index_x.push_back(0);
  index_x.push_back(0);
  index_x.push_back(1);
  index_x.push_back(1);
  index_x.push_back(1);
  index_y.push_back(-1);
  index_y.push_back(0);
  index_y.push_back(1);
  index_y.push_back(-1);
  index_y.push_back(0);
  index_y.push_back(1);
  index_y.push_back(-1);
  index_y.push_back(0);
  index_y.push_back(1);

  for (unsigned int ii=0; ii<index_x.size(); ii++){
    Int_t i1 = index_x.at(ii);
    Int_t i2 = index_y.at(ii);
    Int_t i12 = i1-i2;
    Double_t ARexyi = ARexy.at(ii);
    Double_t AImxyi = AImxy.at(ii);

    Double_t Axyi = sqrt(ARexyi*ARexyi + AImxyi*AImxyi);
    Double_t phixyi = atan2(AImxyi, ARexyi);

    for (unsigned int jj=ii; jj<index_x.size(); jj++){
      Int_t j1 = index_x.at(jj);
      Int_t j2 = index_y.at(jj);
      Int_t j12 = j1-j2;
      Double_t ARexyj = ARexy.at(jj);
      Double_t AImxyj = AImxy.at(jj);

      Double_t Axyj = sqrt(ARexyj*ARexyj + AImxyj*AImxyj);
      Double_t phixyj = atan2(AImxyj, ARexyj);

      Double_t globalFactor = 1;
      if (ii!=jj) globalFactor = 2;
      Double_t phifactor = evaluatePhi1PhiFactor(i1, i2, j1, j2, code, phixyi, phixyj);
      if (phifactor==0) continue;
      Double_t h1h2factor[4]={ 0 };
      for (int ih1=0; ih1<2; ih1++){ // (LL, LR), (RL, RR)
        Int_t hh1 = 1-2*ih1; // L/R
        Double_t h1factor = evaluateH1Factor(i1, j1, hh1, code);
        for (int ih2=0; ih2<2; ih2++){
          Int_t hh2 = 1-2*ih2; // L/R
          Double_t h2factor = evaluateH2Factor(i2, j2, hh2, code);

          Double_t fraction = (1.+((Double_t)hh1)*R1Val+((Double_t)hh2)*R2Val+((Double_t)hh1)*R1Val*((Double_t)hh2)*R2Val)/4.;
          h1h2factor[2*ih1+ih2] = fraction*h1factor*h2factor;
        }
      }
      Double_t hsfactor=evaluateHSFactor(i12, j12, code);

      globalFactor *= Axyi*Axyj*phifactor*hsfactor;
      Double_t result = 0;
      for (int ih1=0; ih1<2; ih1++){ // (LL, LR), (RL, RR)
        for (int ih2=0; ih2<2; ih2++) result += globalFactor*h1h2factor[2*ih1+ih2];
      }
      if (result!=0.) Axxyyterm.push_back(result);
    }
  }
  return;
}

Double_t RooSpinTwo_7DComplex_ppHVV::evaluate() const{
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
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell) term1Coeff = 2.*m1_*GeVunit;
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) term2Coeff = 2.*m2_*GeVunit;

  std::vector<Double_t> Axxyyterm;
  evaluatePolarizationTerms(Axxyyterm, code);
  Double_t value = 0;
  for (unsigned int s=0; s<Axxyyterm.size(); s++) value += Axxyyterm.at(s);
  value *= term1Coeff*term2Coeff*betaVal;

  if (!(value==value)) MELAout << "Evaluate NaN=" << value << endl;
  if (value<=0.){
    MELAout << "Evaluated value<=0: " << value << endl;
    value=epsilon;
  }
  return value;
}

Int_t RooSpinTwo_7DComplex_ppHVV::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const{
  Int_t code = intCodeStart;
  if (checkFundamentalType(h1)){ if (matchArgs(allVars, analVars, h1) || Vdecay1==RooSpin::kVdecayType_GammaOnshell) code *= prime_h1; }
  if (checkFundamentalType(h2)){ if (matchArgs(allVars, analVars, h2) || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_h2; }
  if (checkFundamentalType(hs)){ if (matchArgs(allVars, analVars, hs)) code *= prime_hs; }
  if (checkFundamentalType(Phi)){ if (matchArgs(allVars, analVars, Phi) || Vdecay1==RooSpin::kVdecayType_GammaOnshell || Vdecay2==RooSpin::kVdecayType_GammaOnshell) code *= prime_Phi; }
  if (checkFundamentalType(Phi1)){ if (matchArgs(allVars, analVars, Phi1) || (Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell)) code *= prime_Phi1; }
  if (code==1) code=0;
  return code;
}
Double_t RooSpinTwo_7DComplex_ppHVV::analyticalIntegral(Int_t code, const char* /*rangeName*/) const{
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
  if (Vdecay1!=RooSpin::kVdecayType_GammaOnshell) term1Coeff = 2.*m1_*GeVunit;
  if (Vdecay2!=RooSpin::kVdecayType_GammaOnshell) term2Coeff = 2.*m2_*GeVunit;

  std::vector<Double_t> Axxyyterm;
  evaluatePolarizationTerms(Axxyyterm, code);
  Double_t value = 0;
  for (unsigned int s=0; s<Axxyyterm.size(); s++) value += Axxyyterm.at(s);
  value *= term1Coeff*term2Coeff*betaVal;

  if (!(value==value)){
    MELAout << "Integral NaN=" << value << " at "
      << "h1=" << h1 << '\t'
      << "h2=" << h2 << '\t'
      << "hs=" << hs << '\t'
      << "Phi1=" << Phi1 << '\t'
      << "Phi=" << Phi << '\t'
      << "m1=" << m1_ << '\t'
      << "m2=" << m2_ << '\t'
      << "m12=" << m12 << '\t'
      << endl;
    MELAout << "Possible sources:\n"
      << "betaVal=" << betaVal << '\t'
      << "term1Coeff=" << term1Coeff << '\t'
      << "term2Coeff=" << term2Coeff << '\t'
      << endl;
  }
  if (value<=0.){
    MELAout << "Evaluated integral<=0: " << value << endl;
    value=epsilon;
  }
  return value;
}

void RooSpinTwo_7DComplex_ppHVV::setZZ4fOrdering(Bool_t flag){ ZZ4fOrdering=flag; }
