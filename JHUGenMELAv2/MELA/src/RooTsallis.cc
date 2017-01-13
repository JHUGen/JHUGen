#include <iostream>
#include <cmath>
#include "RooTsallis.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "TH1.h"

//ClassImp(RooTsallis) 

RooTsallisSM::RooTsallisSM(const char *name, const char *title, 
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
			RooAbsReal& _T2):
   
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
   T2("T2","T2",this,_T2)           
 { 
 } 


 RooTsallisSM::RooTsallisSM(const RooTsallisSM& other, const char* name) :  
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
   T2("T2",this,other.T2)           
 {
 } 



 double RooTsallisSM::evaluate() const 
 { 
   // cout<<"In rooTsallisSM::evaluate()"<<endl;
   double T = T0 + mzz*T1 + mzz*mzz*T2;
   double n = n0 + mzz*n1 + mzz*mzz*n2;
   double bb = bb0 + mzz*bb1 + mzz*mzz*bb2;
   double result =  pow(x,ndue)*exp(-bb*x)*pow(1 + (sqrt(x*x + m*m) - m)/(fabs(n*T)),-n);
   return result;
 } 


// LET ROOFIT COMPUTE IT

/* int RooTsallisSM::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  return 0 ;
}

double RooTsallisSM::analyticalIntegral(int code, const char* rangeName) const
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





