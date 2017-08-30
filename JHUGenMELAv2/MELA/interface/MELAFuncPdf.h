#ifndef MELAFUNCPDF
#define MELAFUNCPDF

#include "RooRealProxy.h"
#include "RooAbsPdf.h"

class MELAFuncPdf : public RooAbsPdf{
protected:
  RooRealProxy theFunc;

public:
  MELAFuncPdf() : RooAbsPdf(){}
  MELAFuncPdf(
    const char* name,
    const char* title
    ) : RooAbsPdf(name, title), theFunc("theFunc","theFunc",this){}
  MELAFuncPdf(
    const char* name,
    const char* title,
    RooAbsReal& inFunc
    ) : RooAbsPdf(name, title), theFunc("theFunc", "theFunc", this, inFunc){}
  MELAFuncPdf(const MELAFuncPdf& other, const char* name=0) : RooAbsPdf(other, name), theFunc("theFunc", this, other.theFunc){}
  TObject* clone(const char* newname)const{ return new MELAFuncPdf(*this, newname); }
  inline virtual ~MELAFuncPdf(){}

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0)const{ return dynamic_cast<RooAbsReal*>(theFunc.absArg())->getAnalyticalIntegral(allVars, analVars, rangeName); }
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0)const{ return dynamic_cast<RooAbsReal*>(theFunc.absArg())->analyticalIntegral(code, rangeName); }

  Int_t getAnalyticalIntegralWN(RooArgSet& allVars, RooArgSet& analVars, const RooArgSet *normSet, const char* rangeName=0)const{ return dynamic_cast<RooAbsReal*>(theFunc.absArg())->getAnalyticalIntegralWN(allVars, analVars, normSet, rangeName); }
  Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet, const char* rangeName=0)const{
    RooAbsPdf* pdfcast = dynamic_cast<RooAbsPdf*>(theFunc.absArg());
    if (pdfcast!=0) return pdfcast->analyticalIntegralWN(code, normSet, rangeName);
    else return RooAbsPdf::analyticalIntegralWN(code, normSet, rangeName);
  }

protected:
  Double_t evaluate()const{ return theFunc; }


  ClassDef(MELAFuncPdf, 0)

};

ClassImp(MELAFuncPdf)

#endif
