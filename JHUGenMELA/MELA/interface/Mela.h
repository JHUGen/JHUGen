/*
************* HEADER: CMS MELA interface to MCFM/JHUGen-MELA *************
Please see the ../src/Mela.cc file for the instructions.
*/

#ifndef MELA_Mela_h
#define MELA_Mela_h

#include <vector>
#include "TLorentzVector.h"
#include "TRandom3.h"


class TFile; 
class TGraph;
class TH1F;
class TH2F;
class TH3F;
class RooRealVar;
class RooAbsPdf;
class RooArgSet;
class ScalarPdfFactory_HVV;
class VectorPdfFactory;
class TensorPdfFactory;
class RooqqZZ_JHU_ZgammaZZ_fast;
class ZZMatrixElement;
class SuperMELA;

#include "TVar.hh"
#include "TEvtProb.hh"
#include "MelaPConstant.h"
#include "SuperDijetMela.h"
#include "ScalarPdfFactory_HVV.h"
#include "VectorPdfFactory.h"
#include "TensorPdfFactory_ppHVV.h"
#include "RooqqZZ_JHU_ZgammaZZ_fast.h"

class Mela{

public:

  Mela(double LHCsqrts_=13., double mh_=125., TVar::VerbosityLevel verbosity_=TVar::ERROR); // Higgs mass for supermela
  Mela(const Mela& other);
  ~Mela();

  // Constructor wrapper
  void build(double mh_);

  void setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction);
  void setVerbosity(TVar::VerbosityLevel verbosity_=TVar::ERROR);
  TVar::VerbosityLevel getVerbosity();
  void setMelaLeptonInterference(TVar::LeptonInterference myLepInterf=TVar::DefaultLeptonInterf);
  void setRemoveLeptonMasses(bool MasslessLeptonSwitch=true);
  void setRemoveJetMasses(bool MasslessLeptonSwitch=true);
  void setMelaPrimaryHiggsMass(double myHiggsMass);
  void setMelaHiggsMass(double myHiggsMass, int index=0);
  void setMelaHiggsWidth(double myHiggsWidth=-1, int index=0);
  void setMelaHiggsMassWidth(double myHiggsMass, double myHiggsWidth, int index);
  void setRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);
  void setCandidateDecayMode(TVar::CandidateDecayMode mode);
  void setCurrentCandidateFromIndex(unsigned int icand); // Switches to another candidate
  void setCurrentCandidate(MELACandidate* cand); // Switches to another candidate
  void setInputEvent(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Adds another candidate
  void resetInputEvent(); // Reset the input candidates. Important to call in order to clean up TEvtProb!
  void setTempCandidate(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Adds a temp. candidate
  void appendTopCandidate(SimpleParticleCollection_t* TopDaughters); // Adds a top

  // Function to set EW parameters in MCFM/JHUGen
  void resetMass(double inmass, int ipart);
  void resetWidth(double inwidth, int ipart);
  void resetQuarkMasses();
  void resetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  // Function to get current primary EW/QCD parameters from MCFM/JHUGen (notice Higgs mass/width used in the ME could be different)
  double getPrimaryMass(int ipart);
  double getPrimaryWidth(int ipart);
  double getHiggsWidthAtPoleMass(double mass);


  MelaIO* getIORecord(); // Full parton-by-parton ME record
  MELACandidate* getCurrentCandidate();
  int getCurrentCandidateIndex();
  int getNCandidates();
  std::vector<MELATopCandidate*>* getTopCandidateCollection();


  void getConstant(float& prob); // <ME> constants
  void getPAux(float& prob); // SuperProb


  RooSpin::modelMeasurables getMeasurablesRRV();


  void computeDecayAngles(
    float& qH,
    float& m1,
    float& m2,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& costhetastar,
    float& Phi1
    );
  void computeVBFAngles(
    float& Q2V1,
    float& Q2V2,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& costhetastar,
    float& Phi1
  );
  void computeVBFAngles_ComplexBoost(
    float& Q2V1,
    float& Q2V2,
    float& costheta1_real, float& costheta1_imag,
    float& costheta2_real, float& costheta2_imag,
    float& Phi,
    float& costhetastar,
    float& Phi1
  );
  void computeVHAngles(
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& costhetastar,
    float& Phi1
  );

  void computeP_selfDspin0(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin1(
    double selfDZqqcoupl_input[SIZE_ZQQ][2],
    double selfDZvvcoupl_input[SIZE_ZVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin1(
    double selfDZvvcoupl_input[SIZE_ZVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin2(
    double selfDGggcoupl_input[SIZE_GGG][2],
    double selfDGqqcoupl_input[SIZE_GQQ][2],
    double selfDGvvcoupl_input[SIZE_GVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin2(
    double selfDGggcoupl_input[SIZE_GGG][2],
    double selfDGvvcoupl_input[SIZE_GVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP(
    float& prob,
    bool useConstant=true
    );

  void computeD_CP(
    TVar::MatrixElement myME,
    TVar::Process myType,
    float& prob
    );

  //****VVH Spin-0****//
  void computeProdDecP(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeProdDecP(
    float& prob,
    bool useConstant=true
    );

  //****HJ/HJJ/VBF Spin-0****//
  void computeProdP(
    double selfDHggcoupl_input[SIZE_HGG][2],
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeProdP(
    float& prob,
    bool useConstant=true
    );

  //****VH Spin-0****//
  void computeProdP_VH(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool includeHiggsDecay=false,
    bool useConstant=true
    );
  void computeProdP_VH(
    float& prob,
    bool includeHiggsDecay=false,
    bool useConstant=true
    );

  //***ttH Spin-0****//
  void computeProdP_ttH(
    float& prob,
    int topProcess=2,
    int topDecay=0,
    bool useConstant=true
    );

  // Calculation weight to correct for fermion interference
  void compute4FermionWeight(float& w);

  // Calculation of X propagator
  void getXPropagator(TVar::ResonancePropagatorScheme scheme, float& prop);

  //*** SuperMela ***//
  void computePM4l(
    TVar::SuperMelaSyst syst,
    float& prob
    );

  //*** SuperJJMela ***//
  void computeDijetConvBW(float& prob, bool useTrueBW=false);

  //*** Dgg10 ***//
  void computeD_gg(
    TVar::MatrixElement myME,
    TVar::Process myType,
    float& prob
    );

  // Access ZZMEs Calculate4Momentum
  std::vector<TLorentzVector> calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi);

  /********************/
  /*** Data members ***/
  /********************/
  TRandom3 melaRandomNumber; // Used in SuperMELA smearing
  RooRealVar* mzz_rrv;
  RooRealVar* z1mass_rrv;
  RooRealVar* z2mass_rrv;
  RooRealVar* costhetastar_rrv;
  RooRealVar* costheta1_rrv;
  RooRealVar* costheta2_rrv;
  RooRealVar* phi_rrv;
  RooRealVar* phi1_rrv;
  RooRealVar* Y_rrv;
  RooRealVar* upFrac_rrv;

  RooAbsPdf* pdf;
  ScalarPdfFactory_HVV* ggSpin0Model;
  VectorPdfFactory* spin1Model;
  TensorPdfFactory_ppHVV* spin2Model;
  RooqqZZ_JHU_ZgammaZZ_fast* qqZZmodel;

  SuperMELA* super;

  // Self-define arrays are now members of MELA.
  // There are a lot of them!
  //****Spin-0****//
  // The first dimension (of size [nSupportedHiggses=2]) supports a second resonance present in MCFM
  double selfDHggcoupl[nSupportedHiggses][SIZE_HGG][2];
  double selfDHg4g4coupl[nSupportedHiggses][SIZE_HGG][2];
  double selfDHqqcoupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHbbcoupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHttcoupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHb4b4coupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHt4t4coupl[nSupportedHiggses][SIZE_HQQ][2];
  double selfDHzzcoupl[nSupportedHiggses][SIZE_HVV][2];
  double selfDHwwcoupl[nSupportedHiggses][SIZE_HVV][2];
  double selfDHzzLambda_qsq[nSupportedHiggses][SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ];
  double selfDHwwLambda_qsq[nSupportedHiggses][SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ];
  int selfDHzzCLambda_qsq[nSupportedHiggses][SIZE_HVV_CQSQ];
  int selfDHwwCLambda_qsq[nSupportedHiggses][SIZE_HVV_CQSQ];
  bool differentiate_HWW_HZZ;
  double selfDHzzpcoupl[SIZE_HVV][2];
  double selfDHzpzpcoupl[SIZE_HVV][2];
  double selfDZpffcoupl[SIZE_Vpff][2];
  double selfDHwwpcoupl[SIZE_HVV][2];
  double selfDHwpwpcoupl[SIZE_HVV][2];
  double selfDWpffcoupl[SIZE_Vpff][2];
  double selfDM_Zprime;
  double selfDGa_Zprime;
  double selfDM_Wprime;
  double selfDGa_Wprime;
  //****Spin-1****//
  double selfDZqqcoupl[SIZE_ZQQ][2];
  double selfDZvvcoupl[SIZE_ZVV][2];
  //****Spin-2****//
  double selfDGqqcoupl[SIZE_GQQ][2];
  double selfDGggcoupl[SIZE_GGG][2];
  double selfDGvvcoupl[SIZE_GVV][2];
  double selfDGvvpcoupl[SIZE_GVV][2];
  double selfDGvpvpcoupl[SIZE_GVV][2];
  // That is a lot of them!

protected:
  /********************/
  /*** Data members ***/
  /********************/
  double LHCsqrts;
  TVar::Process myModel_;
  TVar::MatrixElement myME_;
  TVar::Production myProduction_;
  TVar::LeptonInterference myLepInterf_;
  TVar::VerbosityLevel myVerbosity_;

  ZZMatrixElement* ZZME;
  SuperDijetMela* superDijet;


  float auxiliaryProb;

  MELACandidate* melaCand; // Pointer to persistent TEvtProb object

  /***** ME CONSTANT HANDLES *****/
  // Constants that vary with sqrts due to application of PDFs
  //
  MelaPConstant* pAvgSmooth_JHUGen_JQCD_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  //
  MelaPConstant* pAvgSmooth_JHUGen_JJQCD_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  //
  MelaPConstant* pAvgSmooth_JHUGen_JJVBF_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  //
  MelaPConstant* pAvgSmooth_JHUGen_Had_ZH_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  //
  MelaPConstant* pAvgSmooth_JHUGen_Had_WH_HSMHiggs[TVar::nFermionMassRemovalSchemes-1];
  // Decay ME constants that do not use PDFs
  //
  MelaPConstant* pAvgSmooth_JHUGen_ZZGG_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_JHUGen_ZZGG_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_JHUGen_ZZGG_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_S_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_S_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_S_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_S_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_S_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_S_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_S_HSMHiggs_4mu;
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_S_HSMHiggs_4e;
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_S_HSMHiggs_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_ZZGG_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_ZZQQB_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_ZZQQB_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_ZZQQB_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_JJVBF_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_Had_ZH_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_Had_WH_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_JJQCD_bkgZZ_4mu;
  MelaPConstant* pAvgSmooth_MCFM_JJQCD_bkgZZ_4e;
  MelaPConstant* pAvgSmooth_MCFM_JJQCD_bkgZZ_2mu2e;
  //
  MelaPConstant* pAvgSmooth_MCFM_JJQCD_bkgZJets_2l2q;


  /*****************/
  /*** Functions ***/
  /*****************/
  void printLogo() const;

  void setSpinZeroCouplings();
  void setSpinOneCouplings();
  void setSpinTwoCouplings();

  bool configureAnalyticalPDFs();
  void reset_SelfDCouplings();
  void reset_PAux(); // SuperProb reset
  void reset_CandRef();

  void constructDggr(
    float bkg_VAMCFM_noscale,
    float ggzz_VAMCFM_noscale,
    float ggHZZ_prob_pure_noscale,
    float ggHZZ_prob_int_noscale,
    float widthScale,
    float& myDggr
    );

  void getPConstantHandles();
  void deletePConstantHandles();
  MelaPConstant* getPConstantHandle(
    TVar::MatrixElement me_,
    TVar::Production prod_,
    TVar::Process proc_,
    TString relpath,
    TString spname,
    const bool useSqrts=false
    );
  void deletePConstantHandle(MelaPConstant*& handle);
  void computeConstant(float& prob);
  void setConstant();
  float getConstant_JHUGenUndecayed();
  float getConstant_4l();
  float getConstant_2l2q();
  float getConstant_4q();
  float getConstant_FourFermionDecay(const int& decid);

};

#endif

