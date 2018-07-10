// March 28 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)


#ifndef ZZ_COMMON
#define ZZ_COMMON
#include <string>
#include <vector>
// MelaIO class
#include "MelaIO.h"
// Couplings classes
#include "TCouplings.hh"
// MCFM utilities
#include "TMCFMUtils.hh"
// Mod_Parameters
#include "TModParameters.hh"
// NNPDF Driver for JHUGen
#include "TNNPDFDriver.hh"
// Mod Kinematics
#include "TModKinematics.hh"
// JHUGenMELA
#include "TModJHUGen.hh"
#include "TModJHUGenMELA.hh"
// Higgs + 0 jet
#include "TModHiggsMatEl.hh"
#include "TModGravitonMatEl.hh"
#include "TModZprimeMatEl.hh"
// Higgs + 1/2 jets
#include "TModHiggsJJMatEl.hh"
#include "TModHiggsJMatEl.hh"
// VH
#include "TModVHiggsMatEl.hh"
// ttH
#include "TModTTBHMatEl.hh"
// ROOT includes
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TF1.h"


namespace TUtil{
  /// Remove fermion mass if the flag is set to true
  extern bool forbidMassiveLeptons;
  extern bool forbidMassiveJets;
  extern TVar::FermionMassRemoval LeptonMassScheme;
  extern TVar::FermionMassRemoval JetMassScheme;

  void applyLeptonMassCorrection(bool flag=true);
  void applyJetMassCorrection(bool flag=true);
  void setLeptonMassScheme(TVar::FermionMassRemoval scheme=TVar::ConserveDifermionMass);
  void setJetMassScheme(TVar::FermionMassRemoval scheme=TVar::ConserveDifermionMass);
  // This version makes the masses of p1 and p2 to be m1 and m2, leaving p1+p2 unchanged.
  void constrainedRemovePairMass(TLorentzVector& p1, TLorentzVector& p2, double m1=0, double m2=0);
  // This version simply scales momentum to match energy for the desired mass.
  void scaleMomentumToEnergy(const TLorentzVector& massiveJet, TLorentzVector& masslessJet, double mass=0);
  // Function that has generic removal features
  std::pair<TLorentzVector, TLorentzVector> removeMassFromPair(
    TLorentzVector jet1, int jet1Id,
    TLorentzVector jet2, int jet2Id,
    double m1=0, double m2=0
    );
  // Function that adjusts top daughter kinematics
  void adjustTopDaughters(SimpleParticleCollection_t& daughters); // Daughters are arranged as b, Wf, Wfb
  // Compute a fake jet from the massless jets
  void computeFakeJet(TLorentzVector realJet, TLorentzVector others, TLorentzVector& fakeJet); // Input massive + higgs -> output massless fake jet

  // TLorentzVector::Boost in complex plane
  std::pair<TLorentzVector, TLorentzVector> ComplexBoost(TVector3 beta, TLorentzVector p4);

  /// Compute decay angles from the lepton four-vectors and pdgIds.  
  /// Theta1 is the angle corresponding to Z1.
  /// Z1_lept1 and  Z1_lept2 are supposed to come from the same Z.
  /// Leptons are re-ordered internally according to a standard convention:
  /// lept1 = negative-charged lepton (for OS pairs).
  void computeAngles(
    float& costhetastar,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& Phi1,
    TLorentzVector Z1_lept1, int Z1_lept1Id,
    TLorentzVector Z1_lept2, int Z1_lept2Id,
    TLorentzVector Z2_lept1, int Z2_lept1Id,
    TLorentzVector Z2_lept2, int Z2_lept2Id
    );
  void computeAnglesCS(
    float pbeam,
    float& costhetastar,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& Phi1,
    TLorentzVector Z1_lept1, int Z1_lept1Id,
    TLorentzVector Z1_lept2, int Z1_lept2Id,
    TLorentzVector Z2_lept1, int Z2_lept1Id,
    TLorentzVector Z2_lept2, int Z2_lept2Id
    );
  // Angles of associated production
  void computeVBFAngles(
    float& costhetastar,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& Phi1,
    float& Q2V1,
    float& Q2V2,
    TLorentzVector p4M11, int Z1_lept1Id,
    TLorentzVector p4M12, int Z1_lept2Id,
    TLorentzVector p4M21, int Z2_lept1Id,
    TLorentzVector p4M22, int Z2_lept2Id,
    TLorentzVector jet1, int jet1Id,
    TLorentzVector jet2, int jet2Id,
    TLorentzVector* injet1=0, int injet1Id=0, // Gen. partons in lab frame
    TLorentzVector* injet2=0, int injet2Id=0
    );
  void computeVBFAngles_ComplexBoost(
    float& costhetastar,
    float& costheta1_real, float& costheta1_imag,
    float& costheta2_real, float& costheta2_imag,
    float& Phi,
    float& Phi1,
    float& Q2V1,
    float& Q2V2,
    TLorentzVector p4M11, int Z1_lept1Id,
    TLorentzVector p4M12, int Z1_lept2Id,
    TLorentzVector p4M21, int Z2_lept1Id,
    TLorentzVector p4M22, int Z2_lept2Id,
    TLorentzVector jet1, int jet1Id,
    TLorentzVector jet2, int jet2Id,
    TLorentzVector* injet1=0, int injet1Id=0, // Gen. partons in lab frame
    TLorentzVector* injet2=0, int injet2Id=0
    );
  void computeVHAngles(
    float& costhetastar,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& Phi1,
    TLorentzVector p4M11, int Z1_lept1Id,
    TLorentzVector p4M12, int Z1_lept2Id,
    TLorentzVector p4M21, int Z2_lept1Id,
    TLorentzVector p4M22, int Z2_lept2Id,
    TLorentzVector jet1, int jet1Id,
    TLorentzVector jet2, int jet2Id,
    TLorentzVector* injet1=0, int injet1Id=0, // Gen. partons in lab frame
    TLorentzVector* injet2=0, int injet2Id=0
    );

  // Parameter settings
  void SetEwkCouplingParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme);
  void SetMass(double inmass, int ipart);
  void SetDecayWidth(double inwidth, int ipart);
  double GetMass(int ipart);
  double GetDecayWidth(int ipart);
  double GetMass(const MELAParticle* part);
  double GetDecayWidth(const MELAParticle* part);
  void GetMassWidth(int ipart, double& m, double& ga);
  void GetMassWidth(const MELAParticle* part, double& m, double& ga);
  void SetCKMElements(double* invckm_ud, double* invckm_us, double* invckm_cd, double* invckm_cs, double* invckm_ts, double* invckm_tb, double* invckm_ub=0, double* invckm_cb=0, double* invckm_td=0);
  double GetCKMElement(int iquark, int jquark);
  double InterpretScaleScheme(const TVar::Production& production, const TVar::MatrixElement& matrixElement, const TVar::EventScaleScheme& scheme, TLorentzVector p[mxpart]);
  void SetAlphaS(double& Q_ren, double& Q_fac, double multiplier_ren, double multiplier_fac, int mynloop, int mynflav, std::string mypartons); // Q_ren/fac -> Q_ren/fac * multiplier_ren/fac
  void GetAlphaS(double* alphas_, double* alphasmz_); // Get last alpha_s value set
 
  // chooser.f split into 3 different functions
  bool MCFM_chooser(
    const TVar::Process& process, const TVar::Production& production, const TVar::LeptonInterference& leptonInterf,
    const TVar::VerbosityLevel& verbosity,
    const TVar::simple_event_record& mela_event
    );
  bool MCFM_SetupParticleCouplings(
    const TVar::Process& process, const TVar::Production& production,
    const TVar::VerbosityLevel& verbosity,
    const TVar::simple_event_record& mela_event,
    std::vector<int>* partOrder, std::vector<int>* apartOrder
    );
  TString GetMCFMParticleLabel(const int& pid, bool useQJ, bool useExtendedConventions);

  // JHUGen-specific wrappers
  void InitJHUGenMELA(const char* pathtoPDFSet, int PDFMember);
  void SetJHUGenHiggsMassWidth(double MReso, double GaReso);
  void SetJHUGenDistinguishWWCouplings(bool doAllow);
  void ResetAmplitudeIncludes();

  // Spin-0 couplings
  void SetMCFMSpinZeroCouplings(bool useBSM, SpinZeroCouplings* Hcouplings, bool forceZZ);
  void SetJHUGenSpinZeroVVCouplings(double Hvvcoupl[SIZE_HVV][2], double Hvvpcoupl[SIZE_HVV][2], double Hvpvpcoupl[SIZE_HVV][2], int Hvvcoupl_cqsq[SIZE_HVV_CQSQ], double HvvLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ], bool useWWcoupl);
  void SetJHUGenSpinZeroGGCouplings(double Hggcoupl[SIZE_HGG][2]);
  void SetJHUGenSpinZeroQQCouplings(double Hqqcoupl[SIZE_HQQ][2]);
  // Spin-1 couplings
  void SetJHUGenSpinOneCouplings(double Zqqcoupl[SIZE_ZQQ][2], double Zvvcoupl[SIZE_ZVV][2]);
  // Spin-2 couplings
  void SetJHUGenSpinTwoCouplings(double Gacoupl[SIZE_GGG][2], double Gvvcoupl[SIZE_GVV][2], double Gvvpcoupl[SIZE_GVV][2], double Gvpvpcoupl[SIZE_GVV][2], double qLeftRightcoupl[SIZE_GQQ][2]);
  // Vprime / contact couplings
  void SetJHUGenVprimeContactCouplings(double Zpffcoupl[SIZE_Vpff][2], double Wpffcoupl[SIZE_Vpff][2]);

  // PS cuts, unused
  bool MCFM_masscuts(double s[][mxpart], const TVar::Process& process);
  bool MCFM_smalls(double s[][mxpart], int npart);

  // ME computations
  double SumMatrixElementPDF(
    const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement, const TVar::LeptonInterference& leptonInterf,
    TVar::event_scales_type* event_scales, MelaIO* RcdME,
    const double& EBEAM,
    TVar::VerbosityLevel verbosity
    );
  double JHUGenMatEl(
    const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
    TVar::event_scales_type* event_scales, MelaIO* RcdME,
    const double& EBEAM,
    TVar::VerbosityLevel verbosity
    );
  double HJJMatEl(
    const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
    TVar::event_scales_type* event_scales, MelaIO* RcdME,
    const double& EBEAM,
    TVar::VerbosityLevel verbosity
    );
  double VHiggsMatEl(
    const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
    TVar::event_scales_type* event_scales, MelaIO* RcdME,
    const double& EBEAM,
    bool includeHiggsDecay,
    TVar::VerbosityLevel verbosity
    );
  double TTHiggsMatEl(
    const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
    TVar::event_scales_type* event_scales, MelaIO* RcdME,
    const double& EBEAM,
    int topDecay, int topProcess,
    TVar::VerbosityLevel verbosity
    );
  double BBHiggsMatEl(
    const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
    TVar::event_scales_type* event_scales, MelaIO* RcdME,
    const double& EBEAM,
    int botProcess,
    TVar::VerbosityLevel verbosity
    );

  int WipeMEArray(const TVar::Process& process, const TVar::Production& production, const int id[mxpart], double msq[nmsq][nmsq], const TVar::VerbosityLevel& verbosity);
  bool CheckPartonMomFraction(const TLorentzVector& p0, const TLorentzVector& p1, double xx[2], const double& EBEAM, const TVar::VerbosityLevel& verbosity);
  void ComputePDF(const TLorentzVector& p0, const TLorentzVector& p1, double fx1[nmsq], double fx2[nmsq], const double& EBEAM, const TVar::VerbosityLevel& verbosity);
  double SumMEPDF(const TLorentzVector& p0, const TLorentzVector& p1, double msq[nmsq][nmsq], MelaIO* RcdME, const double& EBEAM, const TVar::VerbosityLevel& verbosity);

  // Propagator reweighting
  double ResonancePropagator(double sqrts, TVar::ResonancePropagatorScheme scheme);

  // Boost the particles with or without associated ones to pT=0 frame and return std::vectors filled with (id, momentum) pairs
  void GetBoostedParticleVectors(
    MELACandidate* melaCand,
    TVar::simple_event_record& mela_event,
    TVar::VerbosityLevel verbosity=TVar::DEBUG
    );

  // Convert vectors of simple particles to MELAParticles and create a MELACandidate
  // The output lists could be members of TEvtProb directly.
  MELACandidate* ConvertVectorFormat(
    // Inputs
    SimpleParticleCollection_t* pDaughters, // Cannot be 0
    SimpleParticleCollection_t* pAssociated, // Allowed to be 0
    SimpleParticleCollection_t* pMothers, // Allowed to be 0
    bool isGen,
    // Outputs
    std::vector<MELAParticle*>* particleList,
    std::vector<MELACandidate*>* candList
    );
  // Convert the vector of top daughters (as simple particles) to MELAParticles and create a MELATopCandidate
  // The output lists could be members of TEvtProb directly.
  MELATopCandidate* ConvertTopCandidate(
    // Input
    SimpleParticleCollection_t* TopDaughters,
    // Outputs
    std::vector<MELAParticle*>* particleList,
    std::vector<MELATopCandidate*>* topCandList
    );
  void PrintCandidateSummary(MELACandidate* cand);
  void PrintCandidateSummary(TVar::simple_event_record* cand);

}

#endif

