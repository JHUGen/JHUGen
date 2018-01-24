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
// 
// June 2016
// Ulascan Sarica
//
//-----------------------------------------------------------------------------
//
// STD includes
#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <string>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <cassert>
// ME related includes
#include "TMCFM.hh"
#include "TCouplings.hh"
#include "TVar.hh"
#include "TUtil.hh"
#include "MELAHXSWidth.h"
// ROOT includes
#include "TObject.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"


//----------------------------------------
// Class TEvtProb
//----------------------------------------
class TEvtProb : public TObject{
public:
  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TEvtProb() {};
  TEvtProb(const char* pathtoXSW, double ebeam, const char* pathtoPDFSet, int PDFMember=0, TVar::VerbosityLevel verbosity_=TVar::ERROR);
  TEvtProb(const TEvtProb& other);
  ~TEvtProb();

  //----------------------
  // Functions
  //----------------------
  void Set_LHAgrid(const char* path, int pdfmember=0);
  void SetProcess(TVar::Process proc, TVar::MatrixElement me, TVar::Production prod);
  void SetVerbosity(TVar::VerbosityLevel tmp);
  void SetLeptonInterf(TVar::LeptonInterference tmp);

  void SetCandidateDecayMode(TVar::CandidateDecayMode mode);
  void SetCurrentCandidateFromIndex(unsigned int icand);
  void SetCurrentCandidate(MELACandidate* cand);

  void AllowSeparateWWCouplings(bool doAllow=false);
  void ResetMass(double inmass, int ipart);
  void ResetWidth(double inwidth, int ipart);
  void SetZprimeMassWidth(double inmass, double inwidth);
  void SetWprimeMassWidth(double inmass, double inwidth);
  void ResetQuarkMasses();
  void ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  void ResetCouplings();

  void SetPrimaryHiggsMass(double mass);
  void SetHiggsMass(double mass, double wHiggs=-1., int whichResonance=-1);

  void SetRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);
  void ResetRenFacScaleMode();

  // Convert std::vectors to MELAPArticle* and MELACandidate* objects, stored in particleList and candList, respectively.
  // Also set melaCand to this candidsate if it is valid.
  void SetInputEvent(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    );
  void AppendTopCandidate(SimpleParticleCollection_t* TopDaughters);
  void ResetInputEvent();

  // Reset the IO record, called at te beginning of each comoputation
  void ResetIORecord();

  double XsecCalc_XVV();

  double XsecCalc_VVXVV();

  double XsecCalcXJJ();

  double XsecCalcXJ();

  double XsecCalc_VX(
    bool includeHiggsDecay
    );

  double XsecCalc_TTX(
    int topProcess, int topDecay
    );

  double GetXPropagator(TVar::ResonancePropagatorScheme scheme);

  // Get-functions
  MELAHXSWidth const* GetHXSWidthEstimator() const;
  SpinZeroCouplings* GetSelfDSpinZeroCouplings();
  SpinOneCouplings* GetSelfDSpinOneCouplings();
  SpinTwoCouplings* GetSelfDSpinTwoCouplings();
  VprimeCouplings* GetSelfDVprimeCouplings();
  double GetPrimaryHiggsMass();
  double GetPrimaryMass(int ipart);
  double GetPrimaryWidth(int ipart);
  double GetHiggsWidthAtPoleMass(double mass);
  MelaIO* GetIORecord();
  MELACandidate* GetCurrentCandidate();
  int GetCurrentCandidateIndex(); // Return the index of current melaCand in the candList array, or -1 if it does not exist
  int GetNCandidates();
  std::vector<MELATopCandidate*>* GetTopCandidates();

protected:
  //--------------------
  // Variables
  //--------------------
  const char* pathtoPDFSet_;
  int PDFMember_;
  TVar::Process process;
  TVar::MatrixElement matrixElement;
  TVar::Production production;
  TVar::VerbosityLevel verbosity;
  TVar::LeptonInterference leptonInterf;
  double PrimaryHMass;
  double _hmass;
  double _hwidth;
  double _h2mass;
  double _h2width;
  double EBEAM;
  MELAHXSWidth myCSW_;
  TVar::event_scales_type event_scales;

  SpinZeroCouplings selfDSpinZeroCoupl;
  SpinOneCouplings selfDSpinOneCoupl;
  SpinTwoCouplings selfDSpinTwoCoupl;
  VprimeCouplings selfDVprimeCoupl;
  MelaIO RcdME;

  MELACandidate* melaCand; // Only a pointer to the top-level (input) candList object
  std::vector<MELAParticle*> particleList; // Container of intermediate objects, for bookkeeping to delete later
  std::vector<MELACandidate*> candList; // Container of candidate objects, for bookkeeping to delete later
  std::vector<MELATopCandidate*> topCandList; // Container of candidate objects, for bookkeeping to delete later

  // Initialization functions
  void InitializeMCFM();
  void InitializeJHUGen(const char* pathtoPDFSet, int PDFMember);
  void CrossInitialize();

  // Check if at least one input candidate is present
  bool CheckInputPresent();
  void SetRcdCandPtr();

  // Check if self-defined couplings are specified
  bool CheckSelfDCouplings_Hgg();
  bool CheckSelfDCouplings_Hqq();
  bool CheckSelfDCouplings_Htt();
  bool CheckSelfDCouplings_Hbb();
  bool CheckSelfDCouplings_HVV();

  // Constructor wrapper
  void Build();

  ClassDef(TEvtProb, 0);
};

#endif

