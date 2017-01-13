/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_TSALLISEXP
#define ROO_TSALLISEXP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;
 
class RooTsallisExp : public RooAbsPdf {
public:
  RooTsallisExp(const char *name, const char *title,
	          RooAbsReal& _x, 
                  RooAbsReal& _mzz,     
        	  RooAbsReal& _m,    
              	  RooAbsReal& _n0,   
                  RooAbsReal& _n1,   
	          RooAbsReal& _n2,   
                  RooAbsReal& _ndue, 
                  RooAbsReal& _bb0,  
                  RooAbsReal& _bb1,  
	          RooAbsReal& _bb2,  
                  RooAbsReal& _T0,   
                  RooAbsReal& _T1,   
		  RooAbsReal& _T2,
                  RooAbsReal& _bbdue0,
		  RooAbsReal& _bbdue1,
		  RooAbsReal& _bbdue2,
		  RooAbsReal& _fexp0, 
		  RooAbsReal& _fexp1, 
		  RooAbsReal& _fexp2); 

  RooTsallisExp(const RooTsallisExp& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooTsallisExp(*this,newname); }
  inline virtual ~RooTsallisExp() { }
  /* Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
     Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;*/

protected:

  RooRealProxy x;
  RooRealProxy mzz;
  RooRealProxy m;    
  RooRealProxy n0;   
  RooRealProxy n1;   
  RooRealProxy n2;   
  RooRealProxy ndue; 
  RooRealProxy bb0;  
  RooRealProxy bb1;  
  RooRealProxy bb2;  
  RooRealProxy T0;   
  RooRealProxy T1;   
  RooRealProxy T2;  
  RooRealProxy bbdue0; 
  RooRealProxy bbdue1; 
  RooRealProxy bbdue2; 
  RooRealProxy fexp0;  
  RooRealProxy fexp1;  
  RooRealProxy fexp2;  

  Double_t evaluate() const ;

private:

};
 
#endif
