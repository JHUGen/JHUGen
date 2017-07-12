#include "ScalarPdfFactory_VH.h"


ScalarPdfFactory_VH::ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double sqrts_, RooSpin::VdecayType VHmode1_, RooSpin::VdecayType VHmode2_, Bool_t OnshellH_) :
ScalarPdfFactory(measurables_, false, VHmode1_, VHmode2_, OnshellH_),
sqrts(sqrts_)
{
  if (VHmode1_==RooSpin::kVdecayType_Wany || VHmode1_==RooSpin::kVdecayType_Zuu || VHmode1_==RooSpin::kVdecayType_Zdd || VHmode1_==RooSpin::kVdecayType_Zud) PDFType = 1;
  else PDFType = 2;
  if (PDFType==2) measurables.Y=0;

  makeParamsConst(true);
  initPDF();
}
ScalarPdfFactory_VH::ScalarPdfFactory_VH(RooSpinZero::modelMeasurables measurables_, double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], double sqrts_, bool pmf_applied_, RooSpin::VdecayType VHmode1_, RooSpin::VdecayType VHmode2_, Bool_t OnshellH_) :
ScalarPdfFactory(measurables_, gRatio_, gZGsRatio_, gGsGsRatio_, pmf_applied_, false, VHmode1_, VHmode2_, OnshellH_),
sqrts(sqrts_)
{
  if (VHmode1_==RooSpin::kVdecayType_Wany || VHmode1_==RooSpin::kVdecayType_Zuu || VHmode1_==RooSpin::kVdecayType_Zdd || VHmode1_==RooSpin::kVdecayType_Zud) PDFType = 1;
  else PDFType = 2;
  if (PDFType==2) measurables.Y=0;

  makeParamsConst(true);
  initPDF();
}
ScalarPdfFactory_VH::~ScalarPdfFactory_VH(){
  destroyPDF();
}

void ScalarPdfFactory_VH::makeParamsConst(bool yesNo){
  couplings.Lambda->setConstant(true);
  couplings.Lambda_zgs1->setConstant(true);
  couplings.Lambda_z1->setConstant(true);
  couplings.Lambda_z2->setConstant(true);
  couplings.Lambda_z3->setConstant(true);
  couplings.Lambda_z4->setConstant(true);
  couplings.Lambda_Q->setConstant(true);

  ((RooRealVar*)parameters.mX)->setConstant(yesNo);
  ((RooRealVar*)parameters.gamX)->setConstant(yesNo);
  ((RooRealVar*)parameters.mW)->setConstant(yesNo);
  ((RooRealVar*)parameters.gamW)->setConstant(yesNo);
  ((RooRealVar*)parameters.mZ)->setConstant(yesNo);
  ((RooRealVar*)parameters.gamZ)->setConstant(yesNo);
  ((RooRealVar*)parameters.Sin2ThetaW)->setConstant(yesNo);
  ((RooRealVar*)parameters.vev)->setConstant(yesNo);
}

void ScalarPdfFactory_VH::initPDF(){
  PDF_ILC_5D=0;
  PDF_LHC_3D=0;
  if (PDFType==2){
    PDF_ILC_5D = new RooSpinZero_5D_VH(
      "PDF", "PDF",
      measurables,
      parameters,
      couplings,
      V1decay, V2decay
      );
    PDF_base = (RooSpinZero*)PDF_ILC_5D;
  }
  else if (PDFType==1){
    PDF_LHC_3D = new RooSpinZero_3D_pp_VH(
      "PDF", "PDF",
      measurables,
      parameters,
      couplings,
      sqrts,
      V1decay, V2decay
      );
    PDF_base = (RooSpin*)PDF_LHC_3D;
  }
}

RooSpinZero* ScalarPdfFactory_VH::getPDF(){
  if (PDFType==2) return (RooSpinZero*)PDF_ILC_5D;
  else if (PDFType==1) return (RooSpinZero*)PDF_LHC_3D;
  else return 0;
}

void ScalarPdfFactory_VH::destroyPDF(){
  if (PDF_ILC_5D!=0) delete PDF_ILC_5D;
  if (PDF_LHC_3D!=0) delete PDF_LHC_3D;
}


