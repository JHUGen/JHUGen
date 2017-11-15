#ifndef MELACOMBINEPDFS_H
#define MELACOMBINEPDFS_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include "TDirectory.h"
#include "TH2F.h"

#include <iostream>
#include <fstream>
#include <string>

/*** MELA ***/
class MELAqqZZPdf_v2 : public RooAbsPdf {
public:
  MELAqqZZPdf_v2(){};
  MELAqqZZPdf_v2(
    const char *name, const char *title,
    RooAbsReal& _m4l,
    RooAbsReal& _a0,
    RooAbsReal& _a1,
    RooAbsReal& _a2,
    RooAbsReal& _a3,
    RooAbsReal& _a4,
    RooAbsReal& _a5,
    RooAbsReal& _a6,
    RooAbsReal& _a7,
    RooAbsReal& _a8,
    RooAbsReal& _a9,
    RooAbsReal& _a10,
    RooAbsReal& _a11,
    RooAbsReal& _a12,
    RooAbsReal& _a13
    );
  MELAqqZZPdf_v2(const MELAqqZZPdf_v2& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new MELAqqZZPdf_v2(*this, newname); }
  inline virtual ~MELAqqZZPdf_v2(){}

protected:

  RooRealProxy m4l;
  RooRealProxy a0;
  RooRealProxy a1;
  RooRealProxy a2;
  RooRealProxy a3;
  RooRealProxy a4;
  RooRealProxy a5;
  RooRealProxy a6;
  RooRealProxy a7;
  RooRealProxy a8;
  RooRealProxy a9;
  RooRealProxy a10;
  RooRealProxy a11;
  RooRealProxy a12;
  RooRealProxy a13;


  Double_t evaluate() const;
};

class MELAggZZPdf_v2 : public RooAbsPdf {
public:
  MELAggZZPdf_v2(){};
  MELAggZZPdf_v2(
    const char *name, const char *title,
    RooAbsReal& _m4l,
    RooAbsReal& _a0,
    RooAbsReal& _a1,
    RooAbsReal& _a2,
    RooAbsReal& _a3,
    RooAbsReal& _a4,
    RooAbsReal& _a5,
    RooAbsReal& _a6,
    RooAbsReal& _a7,
    RooAbsReal& _a8,
    RooAbsReal& _a9
    );
  MELAggZZPdf_v2(const MELAggZZPdf_v2& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new MELAggZZPdf_v2(*this, newname); }
  inline virtual ~MELAggZZPdf_v2(){}

protected:

  RooRealProxy m4l;
  RooRealProxy a0;
  RooRealProxy a1;
  RooRealProxy a2;
  RooRealProxy a3;
  RooRealProxy a4;
  RooRealProxy a5;
  RooRealProxy a6;
  RooRealProxy a7;
  RooRealProxy a8;
  RooRealProxy a9;

  Double_t evaluate() const;
};
class MELADoubleCB : public RooAbsPdf {
public:
  MELADoubleCB();
  MELADoubleCB(
    const char *name, const char *title,
    RooAbsReal& _x,
    RooAbsReal& _mean,
    RooAbsReal& _width,
    RooAbsReal& _alpha1,
    RooAbsReal& _n1,
    RooAbsReal& _alpha2,
    RooAbsReal& _n2
    );
  MELADoubleCB(const MELADoubleCB& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new MELADoubleCB(*this, newname); }
  inline virtual ~MELADoubleCB(){}
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:

  RooRealProxy x;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  RooRealProxy n2;

  Double_t evaluate() const;
};

class MELARelBWUFParam : public RooAbsPdf {
public:
  MELARelBWUFParam(){};
  MELARelBWUFParam(
    const char *name, const char *title,
    RooAbsReal& _m4l,
    RooAbsReal& _mH,
    RooAbsReal& _scaleParam,
    Double_t _widthSF=1.
    );
  MELARelBWUFParam(const MELARelBWUFParam& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new MELARelBWUFParam(*this, newname); }
  inline virtual ~MELARelBWUFParam(){}

protected:

  RooRealProxy m4l;
  RooRealProxy mH;
  RooRealProxy scaleParam;
  Double_t widthSF;

  Double_t evaluate() const;
};



#endif
