#ifndef _TEVTPROB_HH_
#define _TEVTPROB_HH_
//-----------------------------------------------------------------------------
// Description: Class TEvtProb: EvtProb base class
// ------------
//
//      Event Probability Density Calculation
//
// Feb 21 2011
// Sergo Jindariani
// Yanyan Gao
//-----------------------------------------------------------------------------
#include <sstream>
#include <string>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

#include "TObject.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "assert.h"
#include "TROOT.h"
// ME related
#include "TMCFM.hh"
#include "TVar.hh"
#include "TUtil.hh"


//----------------------------------------
// Class TEvtProb
//----------------------------------------
class TEvtProb : public TObject {
  
public:
  //--------------------
  // Variables
  //--------------------
  TVar::Process _process;
  TVar::MatrixElement _matrixElement;
  TVar::Production _production;
  double _hmass;
  double _hwidth;
  
  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TEvtProb();
  ~TEvtProb();
  
  //----------------------
  // Function
  //----------------------
  void SetProcess(TVar::Process tmp) { _process = tmp; }
  void SetMatrixElement(TVar::MatrixElement tmp){ _matrixElement = tmp; }
  void SetProduction(TVar::Production tmp){ _production = tmp; }
  double XsecCalc(TVar::Process proc,
		  TVar::Production production,
		  const hzz4l_event_type &hzz4l_event,
		  TVar::VerbosityLevel verbosity);
  double XsecCalcXJJ(TVar::Process proc, TLorentzVector p4[3],
		  TVar::VerbosityLevel verbosity);

  // this appears to be some kind of 
  // way of setting MCFM parameters through
  // an interface defined in TMCFM.hh
  void SetHiggsMass(double mass);
  ClassDef(TEvtProb,0);
};

#endif

