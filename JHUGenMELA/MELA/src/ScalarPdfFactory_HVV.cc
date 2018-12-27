#include "ScalarPdfFactory_HVV.h"


ScalarPdfFactory_HVV::ScalarPdfFactory_HVV(RooSpinZero::modelMeasurables const& measurables_, bool acceptance_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_, Bool_t OnshellH_) :
ScalarPdfFactory(measurables_, acceptance_, V1decay_, V2decay_, OnshellH_)
{
  measurables.Y=0;
  initAcceptanceParams();
  makeParamsConst(true);
  initPDF();
}
ScalarPdfFactory_HVV::ScalarPdfFactory_HVV(
  RooSpinZero::modelMeasurables const& measurables_,
  double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], double gVVpRatio_[1][1], double gVpVpRatio_[1][1],
  bool pmf_applied_, bool acceptance_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_, Bool_t OnshellH_
) : ScalarPdfFactory(measurables_, gRatio_, gZGsRatio_, gGsGsRatio_, gVVpRatio_, gVpVpRatio_, pmf_applied_, acceptance_, V1decay_, V2decay_, OnshellH_)
{
  measurables.Y=0;
  initAcceptanceParams();
  makeParamsConst(true);
  initPDF();
}

ScalarPdfFactory_HVV::~ScalarPdfFactory_HVV(){
  destroyPDF();
  destroyAcceptanceParams();
}

void ScalarPdfFactory_HVV::initAcceptanceParams(){
  if (acceptance){
    accepParams.aPhi = new RooRealVar("aPhi", "aPhi", 1.);
    accepParams.bPhi = new RooRealVar("bPhi", "bPhi", 4.88199e-03);
    accepParams.cPhi = new RooRealVar("cPhi", "cPhi", 3.69579e-02);
    accepParams.dPhi = new RooRealVar("dPhi", "dPhi", 0.);
    accepParams.ePhi = new RooRealVar("ePhi", "ePhi", 0.);

    accepParams.aPhi1 = new RooRealVar("aPhi1", "aPhi1", 1.);
    accepParams.bPhi1 = new RooRealVar("bPhi1", "bPhi1", -1.27958e-02);
    accepParams.cPhi1 = new RooRealVar("cPhi1", "cPhi1", -1.64892e-01);
    accepParams.dPhi1 = new RooRealVar("dPhi1", "dPhi1", 0.);
    accepParams.ePhi1 = new RooRealVar("ePhi1", "ePhi1", 0.);

    accepParams.aH1 = new RooRealVar("aH1", "aH1", 1.);
    accepParams.bH1 = new RooRealVar("bH1", "bH1", 2.64540e-02);
    accepParams.cH1 = new RooRealVar("cH1", "cH1", 0.);
    accepParams.dH1 = new RooRealVar("dH1", "dH1", 0.);
    accepParams.eH1 = new RooRealVar("eH1", "eH1", 0.);

    accepParams.aH2 = new RooRealVar("aH2", "aH2", 1.);
    accepParams.bH2 = new RooRealVar("bH2", "bH2", -3.73167e-01);
    accepParams.cH2 = new RooRealVar("cH2", "cH2", 0.);
    accepParams.dH2 = new RooRealVar("dH2", "dH2", 0.);
    accepParams.eH2 = new RooRealVar("eH2", "eH2", 0.);

    accepParams.aHs = new RooRealVar("aHs", "aHs", 1.);
    accepParams.bHs = new RooRealVar("bHs", "bHs", -1.55528e-01);
    accepParams.cHs = new RooRealVar("cHs", "cHs", 0.);
    accepParams.dHs = new RooRealVar("dHs", "dHs", 0.);
    accepParams.eHs = new RooRealVar("eHs", "eHs", 0.);

    accepParams.aM1 = new RooRealVar("aM1", "aM1", 1.);
    accepParams.bM1 = new RooRealVar("bM1", "bM1", -1.26554e-02);
    accepParams.cM1 = new RooRealVar("cM1", "cM1", 3.13526e-05);
    accepParams.dM1 = new RooRealVar("dM1", "dM1", 0.);

    accepParams.aM2 = new RooRealVar("aM2", "aM2", 1.);
    accepParams.bM2 = new RooRealVar("bM2", "bM2", 5.75519e-04);
    accepParams.cM2 = new RooRealVar("cM2", "cM2", -7.74696e-05);
    accepParams.dM2 = new RooRealVar("dM2", "dM2", 0.);
  }
  else{
    accepParams.aPhi = new RooRealVar("aPhi", "aPhi", 1.);
    accepParams.bPhi = new RooRealVar("bPhi", "bPhi", 0.);
    accepParams.cPhi = new RooRealVar("cPhi", "cPhi", 0.);
    accepParams.dPhi = new RooRealVar("dPhi", "dPhi", 0.);
    accepParams.ePhi = new RooRealVar("ePhi", "ePhi", 0.);

    accepParams.aPhi1 = new RooRealVar("aPhi1", "aPhi1", 1.);
    accepParams.bPhi1 = new RooRealVar("bPhi1", "bPhi1", 0.);
    accepParams.cPhi1 = new RooRealVar("cPhi1", "cPhi1", 0.);
    accepParams.dPhi1 = new RooRealVar("dPhi1", "dPhi1", 0.);
    accepParams.ePhi1 = new RooRealVar("ePhi1", "ePhi1", 0.);

    accepParams.aH1 = new RooRealVar("aH1", "aH1", 1.);
    accepParams.bH1 = new RooRealVar("bH1", "bH1", 0.);
    accepParams.cH1 = new RooRealVar("cH1", "cH1", 0.);
    accepParams.dH1 = new RooRealVar("dH1", "dH1", 0.);
    accepParams.eH1 = new RooRealVar("eH1", "eH1", 0.);

    accepParams.aH2 = new RooRealVar("aH2", "aH2", 1.);
    accepParams.bH2 = new RooRealVar("bH2", "bH2", 0.);
    accepParams.cH2 = new RooRealVar("cH2", "cH2", 0.);
    accepParams.dH2 = new RooRealVar("dH2", "dH2", 0.);
    accepParams.eH2 = new RooRealVar("eH2", "eH2", 0.);

    accepParams.aHs = new RooRealVar("aHs", "aHs", 1.);
    accepParams.bHs = new RooRealVar("bHs", "bHs", 0.);
    accepParams.cHs = new RooRealVar("cHs", "cHs", 0.);
    accepParams.dHs = new RooRealVar("dHs", "dHs", 0.);
    accepParams.eHs = new RooRealVar("eHs", "eHs", 0.);

    accepParams.aM1 = new RooRealVar("aM1", "aM1", 0.);
    accepParams.bM1 = new RooRealVar("bM1", "bM1", 0.);
    accepParams.cM1 = new RooRealVar("cM1", "cM1", 0.);
    accepParams.dM1 = new RooRealVar("dM1", "dM1", 0.);

    accepParams.aM2 = new RooRealVar("aM2", "aM2", 0.);
    accepParams.bM2 = new RooRealVar("bM2", "bM2", 0.);
    accepParams.cM2 = new RooRealVar("cM2", "cM2", 0.);
    accepParams.dM2 = new RooRealVar("dM2", "dM2", 0.);
  }
}
void ScalarPdfFactory_HVV::destroyAcceptanceParams(){
  delete accepParams.aM1;
  delete accepParams.bM1;
  delete accepParams.cM1;
  delete accepParams.dM1;
  delete accepParams.aM2;
  delete accepParams.bM2;
  delete accepParams.cM2;
  delete accepParams.dM2;
  delete accepParams.aPhi;
  delete accepParams.bPhi;
  delete accepParams.cPhi;
  delete accepParams.dPhi;
  delete accepParams.ePhi;
  delete accepParams.aPhi1;
  delete accepParams.bPhi1;
  delete accepParams.cPhi1;
  delete accepParams.dPhi1;
  delete accepParams.ePhi1;
  delete accepParams.aH1;
  delete accepParams.bH1;
  delete accepParams.cH1;
  delete accepParams.dH1;
  delete accepParams.eH1;
  delete accepParams.aH2;
  delete accepParams.bH2;
  delete accepParams.cH2;
  delete accepParams.dH2;
  delete accepParams.eH2;
  delete accepParams.aHs;
  delete accepParams.bHs;
  delete accepParams.cHs;
  delete accepParams.dHs;
  delete accepParams.eHs;
}

void ScalarPdfFactory_HVV::makeParamsConst(bool yesNo){
  SpinPdfFactory::makeParamsConst(yesNo);
  if (acceptance && !yesNo){
    accepParams.aPhi->setConstant(kFALSE);
    accepParams.bPhi->setConstant(kFALSE);
    accepParams.cPhi->setConstant(kFALSE);
    accepParams.dPhi->setConstant(kFALSE);
    accepParams.ePhi->setConstant(kFALSE);

    accepParams.aPhi1->setConstant(kFALSE);
    accepParams.bPhi1->setConstant(kFALSE);
    accepParams.cPhi1->setConstant(kFALSE);
    accepParams.dPhi1->setConstant(kFALSE);
    accepParams.ePhi1->setConstant(kFALSE);

    accepParams.aH1->setConstant(kFALSE);
    accepParams.bH1->setConstant(kFALSE);
    accepParams.cH1->setConstant(kFALSE);
    accepParams.dH1->setConstant(kFALSE);
    accepParams.eH1->setConstant(kFALSE);

    accepParams.aH2->setConstant(kFALSE);
    accepParams.bH2->setConstant(kFALSE);
    accepParams.cH2->setConstant(kFALSE);
    accepParams.dH2->setConstant(kFALSE);
    accepParams.eH2->setConstant(kFALSE);

    accepParams.aHs->setConstant(kFALSE);
    accepParams.bHs->setConstant(kFALSE);
    accepParams.cHs->setConstant(kFALSE);
    accepParams.dHs->setConstant(kFALSE);
    accepParams.eHs->setConstant(kFALSE);

    accepParams.aM1->setConstant(kFALSE);
    accepParams.bM1->setConstant(kFALSE);
    accepParams.cM1->setConstant(kFALSE);
    accepParams.dM1->setConstant(kFALSE);

    accepParams.aM2->setConstant(kFALSE);
    accepParams.bM2->setConstant(kFALSE);
    accepParams.cM2->setConstant(kFALSE);
    accepParams.dM2->setConstant(kFALSE);
  }
  else{
    accepParams.aPhi->setConstant(kTRUE);
    accepParams.bPhi->setConstant(kTRUE);
    accepParams.cPhi->setConstant(kTRUE);
    accepParams.dPhi->setConstant(kTRUE);
    accepParams.ePhi->setConstant(kTRUE);

    accepParams.aPhi1->setConstant(kTRUE);
    accepParams.bPhi1->setConstant(kTRUE);
    accepParams.cPhi1->setConstant(kTRUE);
    accepParams.dPhi1->setConstant(kTRUE);
    accepParams.ePhi1->setConstant(kTRUE);

    accepParams.aH1->setConstant(kTRUE);
    accepParams.bH1->setConstant(kTRUE);
    accepParams.cH1->setConstant(kTRUE);
    accepParams.dH1->setConstant(kTRUE);
    accepParams.eH1->setConstant(kTRUE);

    accepParams.aH2->setConstant(kTRUE);
    accepParams.bH2->setConstant(kTRUE);
    accepParams.cH2->setConstant(kTRUE);
    accepParams.dH2->setConstant(kTRUE);
    accepParams.eH2->setConstant(kTRUE);

    accepParams.aHs->setConstant(kTRUE);
    accepParams.bHs->setConstant(kTRUE);
    accepParams.cHs->setConstant(kTRUE);
    accepParams.dHs->setConstant(kTRUE);
    accepParams.eHs->setConstant(kTRUE);

    accepParams.aM1->setConstant(kTRUE);
    accepParams.bM1->setConstant(kTRUE);
    accepParams.cM1->setConstant(kTRUE);
    accepParams.dM1->setConstant(kTRUE);

    accepParams.aM2->setConstant(kTRUE);
    accepParams.bM2->setConstant(kTRUE);
    accepParams.cM2->setConstant(kTRUE);
    accepParams.dM2->setConstant(kTRUE);
  }
}
void ScalarPdfFactory_HVV::setZZ4fOrdering(bool flag){ PDF->setZZ4fOrdering(flag); }

void ScalarPdfFactory_HVV::initPDF(){
  PDF = new RooSpinZero_7DComplex_withAccep_HVV(
    "PDF", "PDF",
    measurables,
    parameters,
    couplings,
    accepParams,
    V1decay,V2decay
    );
  PDF_base = (RooSpin*)PDF;
}




