#include "RooSpin.h"

using namespace std;


void AnaMelaHelpers::multiplyComplexNumbers(std::vector<Double_t> reals, std::vector<Double_t> imags, Double_t& resRe, Double_t& resIm){
  resRe=0; resIm=0;
  const unsigned int nreals = reals.size();
  const unsigned int nimags = imags.size();
  const unsigned int nloops = pow(2, min(nreals, nimags));
  const unsigned int nterms = max(nreals, nimags);

  for (unsigned int ic=0; ic<nloops; ic++){
    Double_t termval=1;
    unsigned int nimagsused=0;
    bool doSkip=false;
    for (unsigned int it=0; it<nterms; it++){
      unsigned int code = (nreals>=nimags ? 0 : 1);
      if (it<nloops) code = (ic >> it);
      if (code%2==1){
        nimagsused++;
        termval*=imags.at(it);
      }
      else termval*=reals.at(it);
      if (termval==0.){ doSkip=true; break; }
    }
    if (doSkip) continue;
    if (nimagsused%4==0) resRe += termval;
    else if (nimagsused%4==1) resIm += termval;
    else if (nimagsused%4==2) resRe -= termval;
    else resIm -= termval;
  }
}

RooSpin::RooSpin() : RooAbsPdf(),
GeVunit(1e-2)
{}

RooSpin::RooSpin(
  const char* name, const char* title,
  modelMeasurables _measurables,
  modelParameters _parameters,
  RooSpin::VdecayType _Vdecay1, RooSpin::VdecayType _Vdecay2
) : RooAbsPdf(name, title),

h1("h1", "h1", this),
h2("h2", "h2", this),
Phi("Phi", "Phi", this),
m1("m1", "m1", this),
m2("m2", "m2", this),
m12("m12", "m12", this),
hs("hs", "hs", this),
Phi1("Phi1", "Phi1", this),
Y("Y", "Y", this),

mX("mX", "mX", this, (RooAbsReal&)*(_parameters.mX)),
gamX("gamX", "gamX", this, (RooAbsReal&)*(_parameters.gamX)),
mW("mW", "mW", this, (RooAbsReal&)*(_parameters.mW)),
gamW("gamW", "gamW", this, (RooAbsReal&)*(_parameters.gamW)),
mZ("mZ", "mZ", this, (RooAbsReal&)*(_parameters.mZ)),
gamZ("gamZ", "gamZ", this, (RooAbsReal&)*(_parameters.gamZ)),
mWprime("mWprime", "mWprime", this, (RooAbsReal&)*(_parameters.mWprime)),
gamWprime("gamWprime", "gamWprime", this, (RooAbsReal&)*(_parameters.gamWprime)),
mZprime("mZprime", "mZprime", this, (RooAbsReal&)*(_parameters.mZprime)),
gamZprime("gamZprime", "gamZprime", this, (RooAbsReal&)*(_parameters.gamZprime)),
Sin2ThetaW("Sin2ThetaW", "Sin2ThetaW", this, (RooAbsReal&)*(_parameters.Sin2ThetaW)),
vev("vev", "vev", this, (RooAbsReal&)*(_parameters.vev)),
gVprimeff_decay1_left("gVprimeff_decay1_left", "gVprimeff_decay1_left", this, (RooAbsReal&)*(_parameters.gVprimeff_decay1_left)),
gVprimeff_decay1_right("gVprimeff_decay1_right", "gVprimeff_decay1_right", this, (RooAbsReal&)*(_parameters.gVprimeff_decay1_right)),
gVprimeff_decay2_left("gVprimeff_decay2_left", "gVprimeff_decay2_left", this, (RooAbsReal&)*(_parameters.gVprimeff_decay2_left)),
gVprimeff_decay2_right("gVprimeff_decay2_right", "gVprimeff_decay2_right", this, (RooAbsReal&)*(_parameters.gVprimeff_decay2_right)),


Vdecay1(_Vdecay1), Vdecay2(_Vdecay2),
intCodeStart(1),

GeVunit(1e-2)
{
  setProxies(_measurables);
}

RooSpin::RooSpin(const RooSpin& other, const char* name) :
  RooAbsPdf(other, name),

  h1("h1", this, other.h1),
  h2("h2", this, other.h2),
  Phi("Phi", this, other.Phi),
  m1("m1", this, other.m1),
  m2("m2", this, other.m2),
  m12("m12", this, other.m12),
  hs("hs", this, other.hs),
  Phi1("Phi1", this, other.Phi1),
  Y("Y", this, other.Y),

  mX("mX", this, other.mX),
  gamX("gamX", this, other.gamX),
  mW("mW", this, other.mW),
  gamW("gamW", this, other.gamW),
  mZ("mZ", this, other.mZ),
  gamZ("gamZ", this, other.gamZ),
  mWprime("mWprime", this, other.mWprime),
  gamWprime("gamWprime", this, other.gamWprime),
  mZprime("mZprime", this, other.mZprime),
  gamZprime("gamZprime", this, other.gamZprime),
  Sin2ThetaW("Sin2ThetaW", this, other.Sin2ThetaW),
  vev("vev", this, other.vev),
  gVprimeff_decay1_left("gVprimeff_decay1_left", this, other.gVprimeff_decay1_left),
  gVprimeff_decay1_right("gVprimeff_decay1_right", this, other.gVprimeff_decay1_right),
  gVprimeff_decay2_left("gVprimeff_decay2_left", this, other.gVprimeff_decay2_left),
  gVprimeff_decay2_right("gVprimeff_decay2_right", this, other.gVprimeff_decay2_right),

  Vdecay1(other.Vdecay1), Vdecay2(other.Vdecay2),
  intCodeStart(other.intCodeStart),

  GeVunit(other.GeVunit)
{}

void RooSpin::alwaysIntegrate(Int_t code){
  intCodeStart=1;
  if (code%prime_h1==0)intCodeStart *= prime_h1;
  if (code%prime_h2==0)intCodeStart *= prime_h2;
  if (code%prime_hs==0)intCodeStart *= prime_hs;
  if (code%prime_Phi==0)intCodeStart *= prime_Phi;
  if (code%prime_Phi1==0)intCodeStart *= prime_Phi1;
  if (code%prime_m1==0)intCodeStart *= prime_m1;
  if (code%prime_m2==0)intCodeStart *= prime_m2;
  if (code%prime_m12==0)intCodeStart *= prime_m12;
  if (code%prime_Y==0)intCodeStart *= prime_Y;
}

void RooSpin::calculatePropagator(Double_t& propRe, Double_t& propIm, Double_t mass, Int_t propType)const{
  // prop = -i / ((m**2-mV**2) + i*mV*GaV) = - ( mV*GaV + i*(m**2-mV**2) ) / ((m**2-mV**2)**2 + (mV*GaV)**2)
  if (propType==0){ // Photon
    propRe = 0.;
    propIm = (mass!=0. ? -1./pow(mass, 2) : 0.);
    propIm *= pow(GeVunit, -2);
  }
  else if (propType==1){ // Massive vector boson
    Double_t mV, gamV;
    getMVGamV(&mV, &gamV);
    if (gamV>0){
      Double_t denominator = pow(mV*gamV, 2)+pow(pow(mass, 2)-pow(mV, 2), 2);
      propRe = -mV*gamV/denominator;
      propIm = -(pow(mass, 2)-pow(mV, 2))/denominator;
      propRe *= pow(GeVunit, -2);
      propIm *= pow(GeVunit, -2);
    }
    else{
      propRe = (mass==mV ? 1. : 0.);
      propIm = 0.;
    }
  }
  else if (propType==2){ // Higgs prop = i / ((m**2-mX**2) + i*mX*GaX) = - ( mX*GaX + i*(m**2-mX**2) ) / ((m**2-mX**2)**2 + (mX*GaX)**2)
    if (gamX>0.){
      Double_t denominator = pow(mX*gamX, 2)+pow(pow(mass, 2)-pow(mX, 2), 2);
      propRe = mX*gamX/denominator;
      propIm = (pow(mass, 2)-pow(mX, 2))/denominator;
      propRe *= pow(GeVunit, -2);
      propIm *= pow(GeVunit, -2);
    }
    else{
      propRe = (mass==mX ? 1. : 0.);
      propIm = 0.;
    }
  }
  else if (propType==3){ // Massive vector boson Vprime
    Double_t mV, gamV;
    getMVprimeGamVprime(&mV, &gamV);
    if (gamV>=0. && mV>=0.){
      Double_t denominator = pow(mV*gamV, 2)+pow(pow(mass, 2)-pow(mV, 2), 2);
      propRe = -mV*gamV/denominator;
      propIm = -(pow(mass, 2)-pow(mV, 2))/denominator;
      propRe *= pow(GeVunit, -2);
      propIm *= pow(GeVunit, -2);
    }
    else if (mV<0.){
      getMVGamV(&mV, &gamV);
      calculatePropagator(propRe, propIm, mV, 0);
    }
    else{
      propRe = 0;
      propIm = 0;
    }
  }
  else{
    propRe = 1.;
    propIm = 0.;
  }
}
void RooSpin::calculateVffGVGA(Double_t& gV, Double_t& gA, RooSpin::VdecayType Vdecay, bool isGamma)const{
  const Double_t atomicT3 = 0.5;
  const Double_t atomicCharge = 1.;

  const Double_t gW = 2.*mW/vev;
  const Double_t overallFactorZ = gW*0.5/sqrt(1.-Sin2ThetaW); // i*g/(2*cos(thetaW))
  const Double_t overallFactorGamma = -gW*sqrt(Sin2ThetaW); // -i*e*Qf
  const Double_t overallFactorW = gW*0.5/sqrt(2.); // i*g/(2*sqrt(2))

  const Double_t Q_up = 2.*atomicCharge/3.;
  const Double_t Q_dn = -atomicCharge/3.;
  const Double_t Q_l = -atomicCharge;
  const Double_t Q_nu = 0;

  // gV = T3 - 2*Qf*sintW**2
  const Double_t gV_up = atomicT3 - 2.*Q_up*Sin2ThetaW; // ~0.19
  const Double_t gV_dn = -atomicT3 - 2.*Q_dn*Sin2ThetaW; // ~-0.35
  const Double_t gV_l = -atomicT3 - 2.*Q_l*Sin2ThetaW;
  const Double_t gV_nu = atomicT3 - 2.*Q_nu*Sin2ThetaW;

  // gA = T3
  const Double_t gA_up = atomicT3;
  const Double_t gA_dn = -atomicT3;
  const Double_t gA_l = -atomicT3;
  const Double_t gA_nu = atomicT3;

  if (Vdecay==RooSpin::kVdecayType_Zud){
    if (!isGamma){
      double yy_up = pow(gV_up, 2) + pow(gA_up, 2);
      double yy_dn = pow(gV_dn, 2) + pow(gA_dn, 2);
      double xx_up = gV_up*gA_up;
      double xx_dn = gV_dn*gA_dn;
      double yy = (2.*yy_up+3.*yy_dn)/5.;
      double xx = (2.*xx_up+3.*xx_dn)/5.;
      double discriminant = pow(yy, 2)-4.*pow(xx, 2);
      gV = (yy+sqrt(fabs(discriminant)))/2.;
      gA = pow(xx, 2)/gV;
      gV = -overallFactorZ*sqrt(gV); // ~-0.28
      gA = -overallFactorZ*sqrt(gA); // ~-0.5
    }
    else{
      gV = overallFactorGamma*sqrt((2.*pow(Q_up, 2) + 3.*pow(Q_dn, 2))/5.);
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Zdd){
    if (!isGamma){
      gV = overallFactorZ*gV_dn;
      gA = overallFactorZ*gA_dn;
    }
    else{
      gV = overallFactorGamma*Q_dn;
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Zuu){
    if (!isGamma){
      gV = overallFactorZ*gV_up;
      gA = overallFactorZ*gA_up;
    }
    else{
      gV = overallFactorGamma*Q_up;
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Znn){
    if (!isGamma){
      gV = overallFactorZ*gV_nu;
      gA = overallFactorZ*gA_nu;
    }
    else{
      gV = overallFactorGamma*Q_nu;
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Zll){
    if (!isGamma){
      gV = overallFactorZ*gV_l;
      gA = overallFactorZ*gA_l;
    }
    else{
      gV = overallFactorGamma*Q_l;
      gA = 0;
    }
  }
  else if (Vdecay==RooSpin::kVdecayType_Wany){
    if (!isGamma){
      gV = overallFactorW;
      gA = overallFactorW;
    }
    else{
      gV = 0;
      gA = 0;
    }
  }
  else{
    gV = 1;
    gA = 0;
  }
}
void RooSpin::calculateVffR1R2(Double_t& R1Val, Double_t& R2Val, bool isGammaV1, bool isGammaV2)const{
  R1Val=0; R2Val=0;
  Double_t gV1=0, gV2=0, gA1=0, gA2=0;
  calculateVffGVGA(gV1, gA1, Vdecay1, isGammaV1);
  if (gV1!=0. || gA1!=0.) R1Val = 2.*gV1*gA1/(pow(gV1, 2) + pow(gA1, 2));
  calculateVffGVGA(gV2, gA2, Vdecay2, isGammaV2);
  if (gV2!=0. || gA2!=0.) R2Val = 2.*gV2*gA2/(pow(gV2, 2) + pow(gA2, 2));
}
void RooSpin::getMVGamV(Double_t* mV, Double_t* gamV)const{
  if (Vdecay1==RooSpin::kVdecayType_Wany){
    if (mV!=0) (*mV)=mW;
    if (gamV!=0) (*gamV)=gamW;
  }
  else if (!(Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell)){
    if (mV!=0) (*mV)=mZ;
    if (gamV!=0) (*gamV)=gamZ;
  }
  else{
    if (mV!=0) (*mV)=0;
    if (gamV!=0) (*gamV)=0;
  }
}

void RooSpin::calculateVprimeffGVGA(Double_t& gV, Double_t& gA, int whichVprime/*1 or 2*/) const{
  const Double_t gW = 2.*mW/vev;
  const Double_t overallFactor = gW*0.5/sqrt(1.-Sin2ThetaW);

  Double_t gL=0, gR=0;
  switch (whichVprime){
  case 1:
    gL=gVprimeff_decay1_left;
    gR=gVprimeff_decay1_right;
    break;
  case 2:
    gL=gVprimeff_decay2_left;
    gR=gVprimeff_decay2_right;
    break;
  default:
    break;
  };

  gV=(gL+gR)/2.*overallFactor;
  gA=(gL-gR)/2.*overallFactor;
}
void RooSpin::calculateVprimeffR1R2(Double_t& R1Val, Double_t& R2Val) const{
  R1Val=0; R2Val=0;
  Double_t gV1=0, gV2=0, gA1=0, gA2=0;
  calculateVprimeffGVGA(gV1, gA1, 1);
  if (gV1!=0. || gA1!=0.) R1Val = 2.*gV1*gA1/(pow(gV1, 2) + pow(gA1, 2));
  calculateVprimeffGVGA(gV2, gA2, 2);
  if (gV2!=0. || gA2!=0.) R2Val = 2.*gV2*gA2/(pow(gV2, 2) + pow(gA2, 2));
}
void RooSpin::getMVprimeGamVprime(Double_t* mV, Double_t* gamV)const{
  if (Vdecay1==RooSpin::kVdecayType_Wany){
    if (mV!=0) (*mV)=mWprime;
    if (gamV!=0) (*gamV)=gamWprime;
  }
  else if (!(Vdecay1==RooSpin::kVdecayType_GammaOnshell && Vdecay2==RooSpin::kVdecayType_GammaOnshell)){
    if (mV!=0) (*mV)=mZprime;
    if (gamV!=0) (*gamV)=gamZprime;
  }
  else{
    if (mV!=0) (*mV)=0;
    if (gamV!=0) (*gamV)=0;
  }
}

Double_t RooSpin::calculateAmplitudeScale(int VGammaVpmode1, int VGammaVpmode2)const{
  Double_t gV1=0, gV2=0, gA1=0, gA2=0;
  if (VGammaVpmode1<2) calculateVffGVGA(gV1, gA1, Vdecay1, (VGammaVpmode1==1));
  else calculateVprimeffR1R2(gV1, gA1);
  if (VGammaVpmode2<2) calculateVffGVGA(gV2, gA2, Vdecay2, (VGammaVpmode2==1));
  else calculateVprimeffR1R2(gV2, gA2);
  Double_t ampScale = sqrt((pow(gV1, 2) + pow(gA1, 2))*(pow(gV2, 2) + pow(gA2, 2)));
  return ampScale;
}

void RooSpin::setProxies(modelMeasurables _measurables){
  setProxy(h1, (RooAbsReal*) _measurables.h1);
  setProxy(h2, (RooAbsReal*) _measurables.h2);
  setProxy(Phi, (RooAbsReal*) _measurables.Phi);
  setProxy(m1, (RooAbsReal*) _measurables.m1);
  setProxy(m2, (RooAbsReal*) _measurables.m2);
  setProxy(m12, (RooAbsReal*) _measurables.m12);
  setProxy(hs, (RooAbsReal*) _measurables.hs);
  setProxy(Phi1, (RooAbsReal*) _measurables.Phi1);
  setProxy(Y, (RooAbsReal*) _measurables.Y);
}
void RooSpin::setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr){
  if (objectPtr!=0) proxy.setArg((RooAbsReal&) *objectPtr);
}
Bool_t RooSpin::checkFundamentalType(const RooRealProxy& proxy)const{
  RooAbsArg* arg = proxy.absArg();
  return (dynamic_cast<RooRealVar*>(arg)!=0);
}
