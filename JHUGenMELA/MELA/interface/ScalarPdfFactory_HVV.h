#ifndef SCALAR_PDF_FACTORY_HVV
#define SCALAR_PDF_FACTORY_HVV

#include "ScalarPdfFactory.h"
#include "RooSpinZero_7DComplex_withAccep_HVV.h"


class ScalarPdfFactory_HVV : public ScalarPdfFactory{
public:
  RooSpinZero_7DComplex_withAccep_HVV::accepParameters accepParams;

  ScalarPdfFactory_HVV(RooSpin::modelMeasurables measurables_, bool acceptance_=false, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true);
  ScalarPdfFactory_HVV(
    RooSpin::modelMeasurables measurables_,
    double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], double gVVpRatio_[1][1], double gVpVpRatio_[1][1],
    bool pmf_applied_=false, bool acceptance_=false, RooSpin::VdecayType V1decay_=RooSpin::kVdecayType_Zll, RooSpin::VdecayType V2decay_=RooSpin::kVdecayType_Zll, Bool_t OnshellH_=true
  );
  ~ScalarPdfFactory_HVV();

  void makeParamsConst(bool yesNo=true);
  void setZZ4fOrdering(bool flag=true);
  RooSpinZero* getPDF(){ return (RooSpinZero*)PDF; }

protected:
  RooSpinZero_7DComplex_withAccep_HVV* PDF;

  virtual void initAcceptanceParams();
  virtual void destroyAcceptanceParams();

  virtual void initPDF();
  virtual void destroyPDF(){ delete PDF; PDF=0; PDF_base=0; }

};


#endif



