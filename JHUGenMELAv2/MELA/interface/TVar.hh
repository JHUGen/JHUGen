#ifndef EvtProb_VAR
#define EvtProb_VAR

#define xstr_lit(s) str_lit(s)
#define str_lit(s) #s
#ifndef _melapkgpathstr_
#ifndef _melapkgpath_
#define _melapkgpath_ ./
#endif
#define _melapkgpathstr_ xstr_lit(_melapkgpath_)
#endif

#include <cstring>
#include <string>
#include <vector>
#include <utility>
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

#define fbGeV2 0.389379E12
#define smallnumber 1e-15
#define sixteen_2Pi_to_8 3.88650230418250561e+07
#define eight_2Pi_to_5 7.83410393050320417e+04
#define four_2Pi_to_2 39.478417604357432

// typedefs for use in simple_event_record
typedef std::pair<int, TLorentzVector> SimpleParticle_t;
typedef std::vector<SimpleParticle_t> SimpleParticleCollection_t;


namespace TVar{
  enum{
    kNoAssociated=1,
    kUseAssociated_Leptons=2, // l or nu
    kUseAssociated_Photons=3,
    kUseAssociated_Jets=5,
    kUseAssociated_UnstableTops=7,
    kUseAssociated_StableTops=11
  };
  enum CandidateDecayMode{
    CandidateDecay_Stable,
    CandidateDecay_ff,
    CandidateDecay_WW,
    CandidateDecay_ZZ,
    CandidateDecay_ZW, // Untested
    CandidateDecay_ZG,
    CandidateDecay_WG,
    CandidateDecay_GG
  };
  enum VerbosityLevel {
    SILENT = 0,
    ERROR = 1,
    INFO = 2,
    DEBUG = 3,
    DEBUG_VERBOSE = 4,
    DEBUG_MECHECK = 5
  };
  enum MatrixElement{
    MCFM = 0,
    JHUGen = 1,
    ANALYTICAL = 2
  };
  enum Production{
    ZZGG,
    ZZQQB,
    ZZQQB_STU, // Should be the same as ZZQQB, just for crosscheck
    ZZINDEPENDENT,

    ttH, // ttH
    bbH, // bbH

    JQCD, // ? + 1 jet

    JJQCD, // SBF
    JJVBF, // VBF
    JJEW, // VBF+VH (had.)
    JJEWQCD, // VBF+VH+QCD, all hadronic
    Had_ZH, // ZH, Z->uu/dd
    Had_WH, // W(+/-)H, W->ud
    Lep_ZH, // ZH, Z->ll/nunu
    Lep_WH, // W(+/-)H, W->lnu

    // s-channel contributions
    ZZQQB_S,
    JJQCD_S,
    JJVBF_S,
    JJEW_S,
    JJEWQCD_S,
    Had_ZH_S,
    Had_WH_S,
    Lep_ZH_S,
    Lep_WH_S,

    // t+u-channel contributions
    ZZQQB_TU,
    JJQCD_TU,
    JJVBF_TU,
    JJEW_TU,
    JJEWQCD_TU,
    Had_ZH_TU,
    Had_WH_TU,
    Lep_ZH_TU,
    Lep_WH_TU,

    GammaH, // gammaH, stable A (could implement S and TU in the future
    //
    nProductions
  };
  enum LeptonInterference{
    DefaultLeptonInterf,
    InterfOn,
    InterfOff
  };
  enum FermionMassRemoval{
    NoRemoval,
    ConserveDifermionMass,
    MomentumToEnergy,
    nFermionMassRemovalSchemes
  };
  enum ResonancePropagatorScheme{ // Assigned specific integer value on purpose, translated directly to the JHUGen propagator indices
    NoPropagator=0,
    RunningWidth=1,
    FixedWidth=2,
    CPS=3
  };

  enum Process{
    HSMHiggs, // Call this for any MCFM |H|**2-only ME.
    H0_g1prime2,
    H0hplus,
    H0minus,
    H0_Zgsg1prime2,
    H0_Zgs,
    H0_Zgs_PS,
    H0_gsgs,
    H0_gsgs_PS,

    D_g1g1prime2,
    D_g1g2,
    D_g1g2_pi_2,
    D_g1g4,
    D_g1g4_pi_2,
    D_zzzg,
    D_zzgg,
    D_zzzg_PS,
    D_zzgg_PS,
    D_zzzg_g1prime2,
    D_zzzg_g1prime2_pi_2,

    H1minus, // 1-
    H1plus, // 1+

    H2_g1, // 2m+, Zg, gg
    H2_g2, // 2h2+
    H2_g3, // 2h3+
    H2_g4, // 2h+
    H2_g5, // 2b+
    H2_g1g5, // 2m+
    H2_g6, // 2h6+
    H2_g7, // 2h7+
    H2_g8, // 2h-
    H2_g9, // 2h9-
    H2_g10, // 2h10-

    bkgZGamma, // Z+gamma cont.
    bkgZJets, // Z + 0/1/2 jets (ZZGG, JQCD, JJQCD)
    bkgZZ, // qq/gg->ZZ cont.
    bkgWW, // qq/gg->WW cont.
    bkgWWZZ, // gg->ZZ+WW cont.

    bkgZZ_SMHiggs, // ggZZ cont. + SMHigg
    bkgWW_SMHiggs, // ggWW cont. + SMHiggs
    bkgWWZZ_SMHiggs, // ggZZ+WW cont. + SMHiggs

    HSMHiggs_WWZZ, // MCFM |H|**2 ZZ+WW with ZZ-WW interference

    /**** For width ***/
    D_gg10,

    /***** Self Defined******/
    SelfDefine_spin0,
    SelfDefine_spin1,
    SelfDefine_spin2,

    nProcesses
  };
  enum SuperMelaSyst{
    // Nominal value
    SMSyst_None      = 0,
    // Scale uncertainties
    SMSyst_ScaleUp   = 1,
    SMSyst_ScaleDown = 2,
    // Resolution uncertainties
    SMSyst_ResUp     = 3,
    SMSyst_ResDown   = 4
  };
  enum EventScaleScheme{
    DefaultScaleScheme,
    Fixed_mH,
    Fixed_mW,
    Fixed_mZ,
    Fixed_mWPlusmH,
    Fixed_mZPlusmH,
    Fixed_TwomtPlusmH,
    Fixed_mtPlusmH,
    Dynamic_qH,
    Dynamic_qJJH,
    Dynamic_qJJ_qH,
    Dynamic_qJ_qJ_qH,
    Dynamic_HT,

    nEventScaleSchemes
  };

  //---------------------------------
  // Functions
  //---------------------------------
  TString ProcessName(TVar::Process temp);
  TString ProductionName(TVar::Production temp);
  TString MatrixElementName(TVar::MatrixElement temp);

  //---------------------------------
  // Structs
  //---------------------------------
  struct simple_event_record{ // Somewhat not-so-simple particles
    int AssociationCode;
    int AssociationVCompatibility; // Z=23, W+-=|+-24| or none=0
    int nRequested_AssociatedJets;
    int nRequested_AssociatedLeptons;
    int nRequested_AssociatedPhotons;
    int nRequested_Tops;
    int nRequested_Antitops;

    // Output 4-vectors
    std::vector<int> intermediateVid; // Origin of daughters, not associated particles
    SimpleParticleCollection_t pDaughters;
    SimpleParticleCollection_t pAssociated;
    SimpleParticleCollection_t pMothers;

    std::vector<SimpleParticleCollection_t> pTopDaughters;
    std::vector<SimpleParticleCollection_t> pAntitopDaughters;
    SimpleParticleCollection_t pStableTops;
    SimpleParticleCollection_t pStableAntitops;

    // Constructor
    simple_event_record() :
      AssociationCode(TVar::kNoAssociated),
      AssociationVCompatibility(0),
      nRequested_AssociatedJets(0),
      nRequested_AssociatedLeptons(0),
      nRequested_AssociatedPhotons(0),
      nRequested_Tops(0),
      nRequested_Antitops(0)
    {}

  };

  struct event_scales_type{
    TVar::EventScaleScheme renomalizationScheme;
    TVar::EventScaleScheme factorizationScheme;
    double ren_scale_factor;
    double fac_scale_factor;
  };

}


#endif
