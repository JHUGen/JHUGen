#ifndef TENSOR_PDF_FACTORY
#define TENSOR_PDF_FACTORY

#ifdef _def_melatools_
#include <RooSpinTwo.h>
#include <SpinPdfFactory.h>
#else
#include "RooSpinTwo.h"
#include "SpinPdfFactory.h"
#endif
#include "TString.h"
#include "RooFormulaVar.h"


class TensorPdfFactory : public SpinPdfFactory {
public:
  RooSpinTwo::modelCouplings couplings;

  TensorPdfFactory(RooSpin::modelMeasurables measurables_, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true);
  virtual ~TensorPdfFactory();

  virtual void makeParamsConst(bool yesNo)=0;
  virtual void makeCouplingsConst(bool yesNo=true);
  virtual void addHypothesis(int ig, double initval, double iphase=0);
  virtual void setTensorPolarization(int ig, double initval);
  virtual void resetHypotheses();
  virtual RooSpinTwo* getPDF()=0;

protected:
  virtual void initGVals();
  virtual void destroyGVals();

  virtual void initPDF()=0;
  virtual void destroyPDF()=0;
};




#endif



