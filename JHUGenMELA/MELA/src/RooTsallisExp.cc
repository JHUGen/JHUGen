#include <iostream>
#include <cmath>
#include "RooTsallisExp.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "TH1.h"

//ClassImp(RooTsallisExp) 

 RooTsallisExp::RooTsallisExp(const char *name, const char *title, 
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
			      RooAbsReal& _fexp2): 

   
   RooAbsPdf(name,title), 
   x("x","x",this,_x),		    
   mzz("mzz","mzz",this,_mzz),	    
   m("m","m",this,_m),		    
   n0("n0","n0",this,_n0),	    
   n1("n1","n1",this,_n1),	    
   n2("n2","n2",this,_n2),	    
   ndue("ndue","ndue",this,_ndue),  
   bb0("bb0","bb0",this,_bb0),	    
   bb1("bb1","bb1",this,_bb1),	    
   bb2("bb2","bb2",this,_bb2),	    
   T0("T0","T0",this,_T0),	    
   T1("T1","T1",this,_T1),	    
   T2("T2","T2",this,_T2),         
   bbdue0("bbdue0","bbdue0",this,_bbdue0),	  
   bbdue1("bbdue1","bbdue1",this,_bbdue1),	  
   bbdue2("bbdue2","bbdue2",this,_bbdue2),	  
   fexp0("fexp0","fexp0",this,_fexp0),	  
   fexp1("fexp1","fexp1",this,_fexp1),	  
   fexp2("fexp2","fexp2",this,_fexp2)           
 { 
 } 


 RooTsallisExp::RooTsallisExp(const RooTsallisExp& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),		      
   mzz("mzz",this,other.mzz),	    
   m("m",this,other.m),		    
   n0("n0",this,other.n0),	    
   n1("n1",this,other.n1),	    
   n2("n2",this,other.n2),	    
   ndue("ndue",this,other.ndue),  
   bb0("bb0",this,other.bb0),	    
   bb1("bb1",this,other.bb1),	    
   bb2("bb2",this,other.bb2),	    
   T0("T0",this,other.T0),	    
   T1("T1",this,other.T1),	    
   T2("T2",this,other.T2), 
   bbdue0("bbdue0",this,other.bbdue0),	    
   bbdue1("bbdue1",this,other.bbdue1),	    
   bbdue2("bbdue2",this,other.bbdue2),	    
   fexp0("fexp0",this,other.fexp0),	    
   fexp1("fexp1",this,other.fexp1),	    
   fexp2("fexp2",this,other.fexp2) 
          
 {
 } 



 double RooTsallisExp::evaluate() const 
 { 
   // cout<<"In rooTsallis::evaluate()"<<endl;
   double T = T0 + mzz*T1 + mzz*mzz*T2;
   double n = n0 + mzz*n1 + mzz*mzz*n2;
   double bb = bb0 + mzz*bb1 + mzz*mzz*bb2;
   double fexp = fexp0 + mzz*fexp1 + mzz*mzz*fexp2;
   double bbdue = bbdue0 + mzz*bbdue1 + mzz*mzz*bbdue2;
   double result = pow(x,ndue)*exp(-bb*x)*pow(1 + (sqrt(x*x + m*m) - m)/(fabs(n*T)),-n) + fexp*exp(-bbdue*x);;
   return result;
 } 


// LET ROOFIT COMPUTE IT

/* int RooTsallisExp::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  return 0 ;
}

double RooTsallisExp::analyticalIntegral(int code, const char* rangeName) const
{
  switch(code)
    {
    case 1:
      {
	// mathematica dixit
	float term1 = x*pow(1 + (sqrt(x*x + m*m) - m)/(n*T),-n);
        float term2 = m*n*(sqrt(x*x + m*m) + 2*T) - n*n*T*(sqrt(x*x + m*m) + T) -n*m*m - n*x*x + x*x;
	return -term1*term2/((n-2)*(n-1));
      }
    }

assert(0) ;
return 0 ;
} */





