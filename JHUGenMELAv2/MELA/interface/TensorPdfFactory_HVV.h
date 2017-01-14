#ifndef TENSOR_PDF_FACTORY_HZZ
#define TENSOR_PDF_FACTORY_HZZ

#include "RooSpinTwo_7DComplex_HVV.h"
#include "TensorPdfFactory.h"


class TensorPdfFactory_HVV : public TensorPdfFactory {
public:

  TensorPdfFactory_HVV(RooSpin::modelMeasurables measurables_, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true);
  ~TensorPdfFactory_HVV();

  void makeParamsConst(bool yesNo=true);
  RooSpinTwo* getPDF(){ return (RooSpinTwo*)PDF; }

protected:
  RooSpinTwo_7DComplex_HVV* PDF;

  void initPDF();
  void destroyPDF(){ delete PDF; }
};


#endif



