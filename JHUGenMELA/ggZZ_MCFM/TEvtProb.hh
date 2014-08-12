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
  TVar::LeptonInterference _leptonInterf;
  double _hmass;
  double _hwidth;
  double EBEAM;
  
  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TEvtProb(double ebeam = 4000); // sqrt(s)/2*1000
  ~TEvtProb();
  
  //----------------------
  // Function
  //----------------------
  void SetProcess(TVar::Process tmp) { _process = tmp; }
  void SetMatrixElement(TVar::MatrixElement tmp){ _matrixElement = tmp; }
  void SetProduction(TVar::Production tmp){ _production = tmp; }
  void SetLeptonInterf(TVar::LeptonInterference tmp){ _leptonInterf = tmp; }
  void ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ);

  double XsecCalc(TVar::Process proc,
		 TVar::Production production,
		 const hzz4l_event_type &hzz4l_event,
		 TVar::VerbosityLevel verbosity,
		 double couplingvals[SIZE_HVV_FREENORM],
		 double selfDHvvcoupl[SIZE_HVV][2],
		 double selfDZqqcoupl[SIZE_ZQQ][2],
		 double selfDZvvcoupl[SIZE_ZVV][2],
		 double selfDGqqcoupl[SIZE_GQQ][2],
		 double selfDGggcoupl[SIZE_GGG][2],
		 double selfDGvvcoupl[SIZE_GVV][2]);

  double XsecCalcXJJ(TVar::Process proc, TVar::Production production, TLorentzVector p4[3],
		     TVar::VerbosityLevel verbosity,
			 double selfDHggcoupl[SIZE_HGG][2],
			 double selfDHvvcoupl[SIZE_HVV_VBF][2],
			 double selfDHwwcoupl[SIZE_HWW_VBF][2]);

  double XsecCalcXJ(TVar::Process proc, TVar::Production production, TLorentzVector p4[2],
		     TVar::VerbosityLevel verbosity);

  double XsecCalc_VX(TVar::Process proc, TVar::Production production, vh_event_type &vh_event,
			TVar::VerbosityLevel verbosity,
			double selfDHvvcoupl[SIZE_HVV_VBF][2]);

  // this appears to be some kind of 
  // way of setting MCFM parameters through
  // an interface defined in TMCFM.hh
  void SetHiggsMass(double mass=125.6, float wHiggs=-1);
  ClassDef(TEvtProb,0);
};

#endif

