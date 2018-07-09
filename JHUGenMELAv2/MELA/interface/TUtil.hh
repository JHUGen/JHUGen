// March 28 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)


#ifndef ZZ_COMMON
#define ZZ_COMMON
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TF1.h"
#include "TVar.hh"
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
  void scaleMomentumToEnergy(TLorentzVector massiveJet, TLorentzVector& masslessJet, double mass=0);
  // Function that has generic removal features
  pair<TLorentzVector, TLorentzVector> removeMassFromPair(
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
    TLorentzVector Z1_lept1, int Z1_lept1Id,
    TLorentzVector Z1_lept2, int Z1_lept2Id,
    TLorentzVector Z2_lept1, int Z2_lept1Id,
    TLorentzVector Z2_lept2, int Z2_lept2Id,
    float& costhetastar,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& Phi1
    );
  void computeAnglesCS(
    TLorentzVector Z1_lept1, int Z1_lept1Id,
    TLorentzVector Z1_lept2, int Z1_lept2Id,
    TLorentzVector Z2_lept1, int Z2_lept1Id,
    TLorentzVector Z2_lept2, int Z2_lept2Id,
    float pbeam,
    float& costhetastar,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& Phi1
    );
  // Angles of associated production
  void computeVBFangles(
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
  void computeVBFangles_ComplexBoost(
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
  void computeVHangles(
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
  void SetEwkCouplingParameters();
  double InterpretScaleScheme(TVar::Production production, TVar::MatrixElement matrixElement, TVar::EventScaleScheme scheme, TLorentzVector p[mxpart]);
  void SetAlphaS(double& Q_ren, double& Q_fac, double multiplier_ren, double multiplier_fac, int mynloop, int mynflav, string mypartons); // Q_ren/fac -> Q_ren/fac * multiplier_ren/fac
  void GetAlphaS(double* alphas_, double* alphasmz_); // Get last alpha_s value set
  bool MCFM_chooser(TVar::Process process, TVar::Production production, TVar::LeptonInterference leptonInterf, MELACandidate* cand);

  // JHUGen-specific wrappers
  void InitJHUGenMELA(const char* pathtoPDFSet, int PDFMember);
  void SetJHUGenHiggsMassWidth(double MReso, double GaReso);
  void SetJHUGenDistinguishWWCouplings(bool doAllow);

  // Spin-0 couplings
  void SetMCFMSpinZeroVVCouplings(bool useBSM, SpinZeroCouplings* Hcouplings, bool forceZZ);
  void SetJHUGenSpinZeroVVCouplings(double Hvvcoupl[SIZE_HVV][2], int Hvvcoupl_cqsq[3], double HvvLambda_qsq[4][3], bool useWWcoupl);
  void SetJHUGenSpinZeroGGCouplings(double Hggcoupl[SIZE_HGG][2]);
  void SetJHUGenSpinZeroQQCouplings(double Hqqcoupl[SIZE_HQQ][2]);
  // Spin-1 couplings
  void SetJHUGenSpinOneCouplings(double Zqqcoupl[SIZE_ZQQ][2], double Zvvcoupl[SIZE_ZVV][2]);
  // Spin-2 couplings
  void SetJHUGenSpinTwoCouplings(double Gacoupl[SIZE_GGG][2], double Gbcoupl[SIZE_GVV][2], double qLeftRightcoupl[SIZE_GQQ][2]);

  // PS cuts, unused
  bool MCFM_masscuts(double s[][mxpart], TVar::Process process);
  bool MCFM_smalls(double s[][mxpart], int npart);

  // ME computations
  double SumMatrixElementPDF(
    TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
    event_scales_type* event_scales, MelaIO* RcdME,
    double EBEAM,
    double coupling[SIZE_HVV_FREENORM],
    TVar::VerbosityLevel verbosity
    );
  double JHUGenMatEl(
    TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
    event_scales_type* event_scales, MelaIO* RcdME,
    double EBEAM,
    TVar::VerbosityLevel verbosity
    );
  double HJJMatEl(
    TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
    event_scales_type* event_scales, MelaIO* RcdME,
    double EBEAM,
    TVar::VerbosityLevel verbosity
    );
  double VHiggsMatEl(
    TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
    event_scales_type* event_scales, MelaIO* RcdME,
    double EBEAM,
    bool includeHiggsDecay,
    TVar::VerbosityLevel verbosity
    );
  double TTHiggsMatEl(
    TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
    event_scales_type* event_scales, MelaIO* RcdME,
    double EBEAM,
    int topDecay, int topProcess,
    TVar::VerbosityLevel verbosity
    );
  double BBHiggsMatEl(
    TVar::Process process, TVar::Production production, TVar::MatrixElement matrixElement,
    event_scales_type* event_scales, MelaIO* RcdME,
    double EBEAM,
    int botProcess,
    TVar::VerbosityLevel verbosity
    );

  bool CheckPartonMomFraction(const TLorentzVector p0, const TLorentzVector p1, double xx[2], double EBEAM, TVar::VerbosityLevel verbosity);
  void ComputePDF(const TLorentzVector p0, const TLorentzVector p1, double fx1[nmsq], double fx2[nmsq], double EBEAM, TVar::VerbosityLevel verbosity);
  double SumMEPDF(const TLorentzVector p0, const TLorentzVector p1, double msq[nmsq][nmsq], MelaIO* RcdME, double EBEAM, TVar::VerbosityLevel verbosity);

  // Propagator reweighting
  double ResonancePropagator(double shat, TVar::ResonancePropagatorScheme scheme);

  // Boost the particles with or without associated ones to pT=0 frame and return std::vectors filled with (id, momentum) pairs
  void GetBoostedParticleVectors(
    MELACandidate* melaCand,
    simple_event_record& mela_event,
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

}

#endif

