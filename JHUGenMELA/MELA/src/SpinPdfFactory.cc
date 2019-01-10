#include "SpinPdfFactory.h"


using namespace std;
using namespace MELAStreamHelpers;


SpinPdfFactory::SpinPdfFactory(RooSpin::modelMeasurables const& measurables_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_, Bool_t OnshellH_) :
V1decay(V1decay_),
V2decay(V2decay_),
OnshellH(OnshellH_)
{
  initMeasurables(measurables_);
  initMassPole();
  initVdecayParams();
}

SpinPdfFactory::~SpinPdfFactory(){
  destroyVdecayParams();
  destroyMassPole();
}
void SpinPdfFactory::initMeasurables(RooSpin::modelMeasurables const& measurables_){
  measurables.h1 = (RooAbsReal*)measurables_.h1;
  measurables.h2 = (RooAbsReal*)measurables_.h2;
  measurables.Phi = (RooAbsReal*)measurables_.Phi;
  measurables.m1 = (RooAbsReal*)measurables_.m1;
  measurables.m2 = (RooAbsReal*)measurables_.m2;
  measurables.m12 = (RooAbsReal*)measurables_.m12;
  measurables.hs = (RooAbsReal*)measurables_.hs;
  measurables.Phi1 = (RooAbsReal*)measurables_.Phi1;
  measurables.Y = (RooAbsReal*)measurables_.Y;
}
void SpinPdfFactory::initMassPole(){
  if (!OnshellH) parameters.mX = new RooRealVar("mX", "mX", (measurables.m12)->getVal());
  else parameters.mX = measurables.m12;
  parameters.gamX = new RooRealVar("gamX", "gamX", 0);
}
void SpinPdfFactory::resetHiggsMassWidth(Double_t mXval, Double_t gamXval){
  if (!OnshellH){
    if (dynamic_cast<RooRealVar*>(parameters.mX)!=0){ ((RooRealVar*)parameters.mX)->removeMin(); ((RooRealVar*)parameters.mX)->removeMax(); ((RooRealVar*)parameters.mX)->setConstant(false); ((RooRealVar*)parameters.mX)->setVal(mXval); ((RooRealVar*)parameters.mX)->setRange(mXval, mXval); ((RooRealVar*)parameters.mX)->setConstant(true); }
    if (dynamic_cast<RooRealVar*>(parameters.gamX)!=0){ ((RooRealVar*)parameters.gamX)->removeMin(); ((RooRealVar*)parameters.gamX)->removeMax(); ((RooRealVar*)parameters.gamX)->setConstant(false); ((RooRealVar*)parameters.gamX)->setVal(gamXval); ((RooRealVar*)parameters.gamX)->setRange(gamXval, gamXval); ((RooRealVar*)parameters.gamX)->setConstant(true); }
  }
  else MELAerr << "SpinPdfFactory::resetHiggsMassWidth: Higgs mass is already determined by the virtuality with zero width. Cannot set Higgs mass or width" << endl;
}
void SpinPdfFactory::initVdecayParams(){
  if ((((int)V1decay)>0 && ((int)V2decay)<0) || (((int)V1decay)<0 && ((int)V2decay)>0)) MELAerr << "SpinPdfFactory::initVdecayParams: V1 and V2 decays are inconsistent!" << endl;

  const Double_t GfVal = 1.16639e-5;
  const Double_t vevVal = 1./sqrt(GfVal*sqrt(2.));

  parameters.mW = new RooRealVar("mW", "mW", 80.399);
  parameters.gamW = new RooRealVar("gamW", "gamW", 2.085);
  parameters.mZ = new RooRealVar("mZ", "mZ", 91.1876);
  parameters.gamZ = new RooRealVar("gamZ", "gamZ", 2.4952);
  parameters.mWprime = new RooRealVar("mWprime", "mWprime", 0);
  parameters.gamWprime = new RooRealVar("gamWprime", "gamWprime", -1); // Turns off Wprime
  parameters.mZprime = new RooRealVar("mZprime", "mZprime", 0);
  parameters.gamZprime = new RooRealVar("gamZprime", "gamZprime", -1); // Turns off Zprime
  parameters.Sin2ThetaW = new RooRealVar("Sin2ThetaW", "Sin2ThetaW", 0.23119);
  parameters.vev = new RooRealVar("vev", "vev", vevVal);
  parameters.gVprimeff_decay1_left = new RooRealVar("gVprimeff_decay1_left", "gVprimeff_decay1_left", 0);
  parameters.gVprimeff_decay1_right = new RooRealVar("gVprimeff_decay1_right", "gVprimeff_decay1_right", 0);
  parameters.gVprimeff_decay2_left = new RooRealVar("gVprimeff_decay2_left", "gVprimeff_decay2_left", 0);
  parameters.gVprimeff_decay2_right = new RooRealVar("gVprimeff_decay2_right", "gVprimeff_decay2_right", 0);
}
void SpinPdfFactory::resetVdecayParams(Double_t mWval, Double_t gamWval, Double_t mZval, Double_t gamZval, Double_t Sin2ThetaWval, Double_t vevval){
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.mW), mWval);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.gamW), gamWval);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.mZ), mZval);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.gamZ), gamZval);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.Sin2ThetaW), Sin2ThetaWval);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.vev), vevval);
}
void SpinPdfFactory::resetVprimeffCouplings(Double_t gVprimeff_decay1[2], Double_t gVprimeff_decay2[2]){
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.gVprimeff_decay1_left), gVprimeff_decay1[0]);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.gVprimeff_decay1_right), gVprimeff_decay1[1]);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.gVprimeff_decay2_left), gVprimeff_decay2[0]);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.gVprimeff_decay2_right), gVprimeff_decay2[1]);
}
void SpinPdfFactory::resetVprimeMasses(Double_t mWval, Double_t gamWval, Double_t mZval, Double_t gamZval){
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.mWprime), mWval);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.gamWprime), gamWval);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.mZprime), mZval);
  SpinPdfFactory::setVariableValue(dynamic_cast<RooRealVar*>(parameters.gamZprime), gamZval);
}
void SpinPdfFactory::getMVGamV(Double_t* mV, Double_t* gamV)const{
  if (V1decay==RooSpin::kVdecayType_Wany){
    if (mV!=0) (*mV)=(parameters.mW)->getVal();
    if (gamV!=0) (*gamV)=(parameters.gamW)->getVal();
  }
  else if (!(V1decay==RooSpin::kVdecayType_GammaOnshell && V2decay==RooSpin::kVdecayType_GammaOnshell)){
    if (mV!=0) (*mV)=(parameters.mZ)->getVal();
    if (gamV!=0) (*gamV)=(parameters.gamZ)->getVal();
  }
  else{
    if (mV!=0) (*mV)=0;
    if (gamV!=0) (*gamV)=0;
  }
}
void SpinPdfFactory::getMVprimeGamVprime(Double_t* mV, Double_t* gamV)const{
  if (V1decay==RooSpin::kVdecayType_Wany){
    if (mV!=0) (*mV)=(parameters.mWprime)->getVal();
    if (gamV!=0) (*gamV)=(parameters.gamWprime)->getVal();
  }
  else if (!(V1decay==RooSpin::kVdecayType_GammaOnshell && V2decay==RooSpin::kVdecayType_GammaOnshell)){
    if (mV!=0) (*mV)=(parameters.mZprime)->getVal();
    if (gamV!=0) (*gamV)=(parameters.gamZprime)->getVal();
  }
  else{
    if (mV!=0) (*mV)=0;
    if (gamV!=0) (*gamV)=0;
  }
}
void SpinPdfFactory::resetVdecay(RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_){
  if ((((int)V1decay)>0 && ((int)V2decay)<0) || (((int)V1decay)<0 && ((int)V2decay)>0)) MELAerr << "SpinPdfFactory::resetVdecay: V1 and V2 decays are inconsistent!" << endl;

  V1decay=V1decay_;
  V2decay=V2decay_;
  PDF_base->setDecayModes(V1decay, V2decay);
}

void SpinPdfFactory::destroyMassPole(){
  delete parameters.gamX;
  if (!OnshellH) delete parameters.mX;
}
void SpinPdfFactory::destroyVdecayParams(){
  delete parameters.gVprimeff_decay2_right;
  delete parameters.gVprimeff_decay2_left;
  delete parameters.gVprimeff_decay1_right;
  delete parameters.gVprimeff_decay1_left;
  delete parameters.vev;
  delete parameters.Sin2ThetaW;
  delete parameters.gamZ;
  delete parameters.mZ;
  delete parameters.gamW;
  delete parameters.mW;
  delete parameters.gamZprime;
  delete parameters.mZprime;
  delete parameters.gamWprime;
  delete parameters.mWprime;
}

void SpinPdfFactory::setVariableValue(RooRealVar* var, Double_t val){
  if (var){
    var->removeMin();
    var->removeMax();
    var->setConstant(false);
    var->setVal(val);
    var->setRange(val, val);
    var->setConstant(true);
  }
}

void SpinPdfFactory::makeParamsConst(bool yesNo){
  ((RooRealVar*) parameters.mX)->setConstant(yesNo);
  ((RooRealVar*) parameters.gamX)->setConstant(yesNo);
  ((RooRealVar*) parameters.mW)->setConstant(yesNo);
  ((RooRealVar*) parameters.gamW)->setConstant(yesNo);
  ((RooRealVar*) parameters.mZ)->setConstant(yesNo);
  ((RooRealVar*) parameters.gamZ)->setConstant(yesNo);
  ((RooRealVar*) parameters.mWprime)->setConstant(yesNo);
  ((RooRealVar*) parameters.gamWprime)->setConstant(yesNo);
  ((RooRealVar*) parameters.mZprime)->setConstant(yesNo);
  ((RooRealVar*) parameters.gamZprime)->setConstant(yesNo);
  ((RooRealVar*) parameters.Sin2ThetaW)->setConstant(yesNo);
  ((RooRealVar*) parameters.vev)->setConstant(yesNo);
  ((RooRealVar*) parameters.gVprimeff_decay1_left)->setConstant(yesNo);
  ((RooRealVar*) parameters.gVprimeff_decay1_right)->setConstant(yesNo);
  ((RooRealVar*) parameters.gVprimeff_decay2_left)->setConstant(yesNo);
  ((RooRealVar*) parameters.gVprimeff_decay2_right)->setConstant(yesNo);
}

void SpinPdfFactory::setVerbosity(TVar::VerbosityLevel verbosity){ PDF_base->setVerbosity(verbosity); }
