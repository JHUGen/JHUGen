#include "SpinPdfFactory.h"


SpinPdfFactory::SpinPdfFactory(RooSpin::modelMeasurables measurables_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_, Bool_t OnshellH_) :
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
void SpinPdfFactory::initMeasurables(RooSpin::modelMeasurables measurables_){
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
  else cerr << "SpinPdfFactory::resetHiggsMassWidth: Higgs mass is already determined by the virtuality with zero width. Cannot set Higgs mass or width" << endl;
}
void SpinPdfFactory::initVdecayParams(){
  if ((((int)V1decay)>0 && ((int)V2decay)<0) || (((int)V1decay)<0 && ((int)V2decay)>0)) cerr << "SpinPdfFactory::initVdecayParams: V1 and V2 decays are inconsistent!" << endl;

  const Double_t GfVal = 1.16639e-5;
  const Double_t vevVal = 1./sqrt(GfVal*sqrt(2.));

  parameters.mW = new RooRealVar("mW", "mW", 80.399);
  parameters.gamW = new RooRealVar("gamW", "gamW", 2.085);
  parameters.mZ = new RooRealVar("mZ", "mZ", 91.1876);
  parameters.gamZ = new RooRealVar("gamZ", "gamZ", 2.4952);
  parameters.Sin2ThetaW = new RooRealVar("Sin2ThetaW", "Sin2ThetaW", 0.23119);
  parameters.vev = new RooRealVar("vev", "vev", vevVal);
}
void SpinPdfFactory::resetVdecayParams(Double_t mWval, Double_t gamWval, Double_t mZval, Double_t gamZval, Double_t Sin2ThetaWval, Double_t vevval){
  if (dynamic_cast<RooRealVar*>(parameters.mW)!=0){ ((RooRealVar*)parameters.mW)->removeMin(); ((RooRealVar*)parameters.mW)->removeMax(); ((RooRealVar*)parameters.mW)->setConstant(false); ((RooRealVar*)parameters.mW)->setVal(mWval); ((RooRealVar*)parameters.mW)->setRange(mWval, mWval); ((RooRealVar*)parameters.mW)->setConstant(true); }
  if (dynamic_cast<RooRealVar*>(parameters.gamW)!=0){ ((RooRealVar*)parameters.gamW)->removeMin(); ((RooRealVar*)parameters.gamW)->removeMax(); ((RooRealVar*)parameters.gamW)->setConstant(false); ((RooRealVar*)parameters.gamW)->setVal(gamWval); ((RooRealVar*)parameters.gamW)->setRange(gamWval, gamWval); ((RooRealVar*)parameters.gamW)->setConstant(true); }
  if (dynamic_cast<RooRealVar*>(parameters.mZ)!=0){ ((RooRealVar*)parameters.mZ)->removeMin(); ((RooRealVar*)parameters.mZ)->removeMax(); ((RooRealVar*)parameters.mZ)->setConstant(false); ((RooRealVar*)parameters.mZ)->setVal(mZval); ((RooRealVar*)parameters.mZ)->setRange(mZval, mZval); ((RooRealVar*)parameters.mZ)->setConstant(true); }
  if (dynamic_cast<RooRealVar*>(parameters.gamZ)!=0){ ((RooRealVar*)parameters.gamZ)->removeMin(); ((RooRealVar*)parameters.gamZ)->removeMax(); ((RooRealVar*)parameters.gamZ)->setConstant(false); ((RooRealVar*)parameters.gamZ)->setVal(gamZval); ((RooRealVar*)parameters.gamZ)->setRange(gamZval, gamZval); ((RooRealVar*)parameters.gamZ)->setConstant(true); }
  if (dynamic_cast<RooRealVar*>(parameters.Sin2ThetaW)!=0){ ((RooRealVar*)parameters.Sin2ThetaW)->removeMin(); ((RooRealVar*)parameters.Sin2ThetaW)->removeMax(); ((RooRealVar*)parameters.Sin2ThetaW)->setConstant(false); ((RooRealVar*)parameters.Sin2ThetaW)->setVal(Sin2ThetaWval); ((RooRealVar*)parameters.Sin2ThetaW)->setRange(Sin2ThetaWval, Sin2ThetaWval); ((RooRealVar*)parameters.Sin2ThetaW)->setConstant(true); }
  if (dynamic_cast<RooRealVar*>(parameters.vev)!=0){ ((RooRealVar*)parameters.vev)->removeMin(); ((RooRealVar*)parameters.vev)->removeMax(); ((RooRealVar*)parameters.vev)->setConstant(false); ((RooRealVar*)parameters.vev)->setVal(vevval); ((RooRealVar*)parameters.vev)->setRange(vevval, vevval); ((RooRealVar*)parameters.vev)->setConstant(true); }
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
void SpinPdfFactory::resetVdecay(RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_){
  if ((((int)V1decay)>0 && ((int)V2decay)<0) || (((int)V1decay)<0 && ((int)V2decay)>0)) cerr << "SpinPdfFactory::resetVdecay: V1 and V2 decays are inconsistent!" << endl;

  V1decay=V1decay_;
  V2decay=V2decay_;
  PDF_base->setDecayModes(V1decay, V2decay);
}

void SpinPdfFactory::destroyMassPole(){
  delete parameters.gamX;
  if (!OnshellH) delete parameters.mX;
}
void SpinPdfFactory::destroyVdecayParams(){
  delete parameters.vev;
  delete parameters.Sin2ThetaW;
  delete parameters.gamZ;
  delete parameters.mZ;
  delete parameters.gamW;
  delete parameters.mW;
}


