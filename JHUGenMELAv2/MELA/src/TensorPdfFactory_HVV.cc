#include "TensorPdfFactory_HVV.h"


TensorPdfFactory_HVV::TensorPdfFactory_HVV(RooSpin::modelMeasurables measurables_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_, Bool_t OnshellH_) :
TensorPdfFactory(measurables_, V1decay_, V2decay_, OnshellH_)
{
  measurables.Y=0;
  makeParamsConst(true);
  initPDF();
}

TensorPdfFactory_HVV::~TensorPdfFactory_HVV(){
  destroyPDF();
}

void TensorPdfFactory_HVV::makeParamsConst(bool yesNo){
  couplings.Lambda->setConstant(true);

  ((RooRealVar*)parameters.mX)->setConstant(yesNo);
  ((RooRealVar*)parameters.gamX)->setConstant(yesNo);
  ((RooRealVar*)parameters.mW)->setConstant(yesNo);
  ((RooRealVar*)parameters.gamW)->setConstant(yesNo);
  ((RooRealVar*)parameters.mZ)->setConstant(yesNo);
  ((RooRealVar*)parameters.gamZ)->setConstant(yesNo);
  ((RooRealVar*)parameters.Sin2ThetaW)->setConstant(yesNo);
  ((RooRealVar*)parameters.vev)->setConstant(yesNo);
}

void TensorPdfFactory_HVV::initPDF(){
  PDF = new RooSpinTwo_7DComplex_HVV(
    "PDF", "PDF",
    measurables,
    parameters,
    couplings,
    V1decay, V2decay
    );
  PDF_base = (RooSpin*)PDF;
}




