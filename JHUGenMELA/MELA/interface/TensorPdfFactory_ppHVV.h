#ifndef TENSOR_PDF_FACTORY_PPHVV
#define TENSOR_PDF_FACTORY_PPHVV

#include "RooSpinTwo_7DComplex_ppHVV.h"
#include "TensorPdfFactory.h"


class TensorPdfFactory_ppHVV : public TensorPdfFactory{
public:

  TensorPdfFactory_ppHVV(RooSpin::modelMeasurables measurables_, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true);
  ~TensorPdfFactory_ppHVV();

  void setZZ4fOrdering(bool flag=true);
  RooSpinTwo* getPDF(){ return (RooSpinTwo*)PDF; }

protected:
  RooSpinTwo_7DComplex_ppHVV* PDF;

  virtual void initPDF();
  virtual void destroyPDF(){ delete PDF; PDF=0; PDF_base=0; }

};


#endif



