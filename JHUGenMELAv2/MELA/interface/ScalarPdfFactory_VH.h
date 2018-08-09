#ifndef SCALAR_PDF_FACTORY_VH
#define SCALAR_PDF_FACTORY_VH

#ifdef _def_melatools_
#include <RooSpinZero_5D_VH.h>
#include <RooSpinZero_3D_pp_VH.h>
#include <ScalarPdfFactory.h>
#else
#include "RooSpinZero_5D_VH.h"
#include "RooSpinZero_3D_pp_VH.h"
#include "ScalarPdfFactory.h"
#endif


class ScalarPdfFactory_VH : public ScalarPdfFactory {
public:

  ScalarPdfFactory_VH(RooSpin::modelMeasurables measurables_, double sqrts_, RooSpin::VdecayType VHmode1_=RooSpin::kVdecayType_Zud, RooSpin::VdecayType VHmode2_=RooSpin::kVdecayType_Zud, Bool_t OnshellH_=true);
  ScalarPdfFactory_VH(RooSpin::modelMeasurables measurables_, double gRatio_[4][8], double gZGsRatio_[4][1], double gGsGsRatio_[3][1], double sqrts_, bool pmf_applied_=false, RooSpin::VdecayType VHmode1_=RooSpin::kVdecayType_Zud, RooSpin::VdecayType VHmode2_=RooSpin::kVdecayType_Zud, Bool_t OnshellH_=true);
  ~ScalarPdfFactory_VH();

  void makeParamsConst(bool yesNo=true);
  RooSpinZero* getPDF();

protected:
  RooSpinZero_5D_VH* PDF_ILC_5D;
  RooSpinZero_3D_pp_VH* PDF_LHC_3D;
  double sqrts;
  int PDFType;

  void initPDF();
  void destroyPDF();
};

#endif



