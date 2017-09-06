#ifndef VECTOR_PDF_FACTORY
#define VECTOR_PDF_FACTORY

#include <cmath>
#include "TVar.hh"
#include "RooSpinOne_7D.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "TF1.h"


class VectorPdfFactory{

public:
    
  RooRealVar* mZ;
  RooRealVar* gamZ;

  RooRealVar* R1Val;
  RooRealVar* R2Val;

  RooSpinOne_7D* PDF;

  RooRealVar* g1Val;
  RooRealVar* g2Val;

  RooRealVar* aParam;
  
  VectorPdfFactory(){};

  VectorPdfFactory(RooRealVar* m1,RooRealVar* m2,RooRealVar* hs,RooRealVar* h1,RooRealVar* h2,RooRealVar* Phi,RooRealVar* Phi1,RooRealVar* mZZ){

    // Parameters
    mZ     = new RooRealVar("mZ","mZ",91.188);
    gamZ   = new RooRealVar("gamZ","gamZ",2.5);

    // related to tensor structure of V decays
    R1Val  = new RooRealVar("R1Val","R1Val",0.15);
    R2Val  = new RooRealVar("R2Val","R2Val",0.15);

    // dimensionless couplings
    g1Val = new RooRealVar("g1Val", "g1Val", 0.0);        
    g2Val = new RooRealVar("g2Val", "g2Val", 0.0);

    // random paramter (?)
    aParam = new RooRealVar("aParam","aParam",0.0);

    PDF = new RooSpinOne_7D("PDF","PDF", *mZZ, *m1, *m2, *h1, *h2,*hs, *Phi, *Phi1, 
			    *g1Val, *g2Val, *R1Val, *R2Val, *aParam, *mZ, *gamZ);

  };

  ~VectorPdfFactory(){

    delete g1Val;
    delete g2Val; 

    delete aParam;

    delete mZ;
    delete gamZ;

    delete R1Val;
    delete R2Val;

    delete PDF;

  };

  int configure(TVar::Process model_){

    switch (model_){
    case TVar::H1plus: makePseudoZprime(); return 0; break;
    case TVar::H1minus: makeZprime(); return 0; break;
    case TVar::SelfDefine_spin1 : return 0; break;
    default: makeZprime(); return 1; break;
    }

  };


  void makePseudoZprime(){  // NEED TO CALCULATE NORMALIZATIONS

    g1Val->setVal(0.0);
    g2Val->setVal(1.0); 

  };

  void makeZprime(){  // NEED TO CALCULATE NORMALIZATIONS

    g1Val->setVal(1.0); 
    g2Val->setVal(0.0); 

  };

  void makeParamsConst(bool yesNo=true){
    if(yesNo){

      g1Val->setConstant(kTRUE);
      g2Val->setConstant(kTRUE);

      gamZ->setConstant(kTRUE);
      mZ->setConstant(kTRUE);
      R1Val->setConstant(kTRUE);
      R2Val->setConstant(kTRUE);

    }else{

      g1Val->setConstant(kFALSE);
      g2Val->setConstant(kFALSE);

      gamZ->setConstant(kFALSE);
      mZ->setConstant(kFALSE);
      R1Val->setConstant(kFALSE);
      R2Val->setConstant(kFALSE);
    }
  };

  void setZZ4fOrdering(bool flag=true){ PDF->setZZ4fOrdering(flag); }

};

#endif


