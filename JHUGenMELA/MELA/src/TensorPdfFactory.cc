#include "TensorPdfFactory.h"


TensorPdfFactory::TensorPdfFactory(RooSpin::modelMeasurables measurables_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_, Bool_t OnshellH_) :
SpinPdfFactory(measurables_, V1decay_, V2decay_, OnshellH_)
{
  initGVals();
}

TensorPdfFactory::~TensorPdfFactory(){
  destroyGVals();
}
void TensorPdfFactory::initGVals(){
  couplings.Lambda = new RooRealVar("Lambda", "Lambda", 1000.);

  for (int v=0; v<(int)SIZE_GVV; v++){
    for (int im=0; im<2; im++){
      TString strcore;
      double initval = 0;
      TString strapp = "Val";
      if (im==1) strapp.Append("Im");
      strapp.Prepend(Form("%i", v+1));

      strcore = "g";
      strcore.Append(strapp);
      RooRealVar* gVal = new RooRealVar(strcore, strcore, initval, -1e15, 1e15);
      gVal->removeMin();
      gVal->removeMax();
      couplings.bList[v][im] = (RooAbsReal*)gVal;
    }
  }

  for (int f=0; f<2; f++){
    TString strcore = Form("f_spinz%i", f+1);
    RooRealVar* fVal = new RooRealVar(strcore, strcore, 0., 0., 1.);
    if (f==0) couplings.f_spinz1 = (RooRealVar*)fVal;
    else couplings.f_spinz2 = (RooRealVar*)fVal;
  }
}
void TensorPdfFactory::destroyGVals(){
  for (int v=0; v<(int)SIZE_GVV; v++){
    for (int im=0; im<2; im++){
      delete couplings.bList[v][im];
    }
  }
  delete couplings.Lambda;

  delete couplings.f_spinz1;
  delete couplings.f_spinz2;
}

void TensorPdfFactory::addHypothesis(int ig, double initval, double iphase){
  if (ig>=(int)SIZE_GVV || ig<0) cerr << "Invalid g" << ig << endl;
  else{
    ((RooRealVar*)couplings.bList[ig][0])->setVal(initval*cos(iphase));
    ((RooRealVar*)couplings.bList[ig][1])->setVal(initval*sin(iphase));
  }
}
void TensorPdfFactory::setTensorPolarization(int ig, double initval){
  if (ig>2 || ig<=0) cerr << "Cannot set f_spinz" << ig << ". Please st f_spinz1 or f_spinz2 only." << endl;
  else{
    if (ig==1) ((RooRealVar*)couplings.f_spinz1)->setVal(initval);
    else ((RooRealVar*)couplings.f_spinz2)->setVal(initval);
  }
}
void TensorPdfFactory::resetHypotheses(){
  for (int ig=0; ig<(int)SIZE_GVV; ig++){
    for (int im=0; im<2; im++) ((RooRealVar*)couplings.bList[ig][im])->setVal(0.);
  }
  ((RooRealVar*)couplings.f_spinz1)->setVal(0.);
  ((RooRealVar*)couplings.f_spinz2)->setVal(0.);
}

void TensorPdfFactory::makeCouplingsConst(bool yesNo){
  couplings.Lambda->setConstant(true); // The user is not allowed to change this value!
  
  // Set fqq, fz2
  couplings.f_spinz1->setConstant(yesNo);
  couplings.f_spinz2->setConstant(yesNo);

  // Set the b decay couplings
  for (int ig=0; ig<(int)SIZE_GVV; ig++){
    for (int im=0; im<2; im++){
      if (dynamic_cast<RooRealVar*>(couplings.bList[ig][im])!=0) ((RooRealVar*)couplings.bList[ig][im])->setConstant(yesNo);
    }
  }
}


