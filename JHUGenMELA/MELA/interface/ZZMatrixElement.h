#ifndef ZZMatrixElement_ZZMatrixElement_h
#define ZZMatrixElement_ZZMatrixElement_h

#include <vector>
#include "TVar.hh"
#include "TEvtProb.hh"


class  ZZMatrixElement{
public:

  ZZMatrixElement(
    const char* pathtoPDFSet,
    int PDFMember,
    const char* pathtoHiggsCSandWidth, // path to the textfiles of the HiggsCSandWidth package
    double ebeam,
    TVar::VerbosityLevel verbosity
    );
  ZZMatrixElement(const ZZMatrixElement& other);
  ~ZZMatrixElement();

  void computeXS(
    float &mevalue
    );
  void computeProdXS_VVHVV(
    float& mevalue
    );
  void computeProdXS_JJH(
    float &mevalue
    );
  void computeProdXS_JH(
    float &mevalue
    );
  void computeProdXS_VH(
    float &mevalue,
    bool includeHiggsDecay=false
    );
  void computeProdXS_ttH(
    float &mevalue,
    int topProcess,
    int topDecay=0
    );

  void get_XPropagator(TVar::ResonancePropagatorScheme scheme, float& prop);

  // Set-functions
  void set_Process(TVar::Process process_, TVar::MatrixElement me_, TVar::Production production_); // Sets variables in Xcal2 as well
  void set_Verbosity(TVar::VerbosityLevel verbosity_); // Sets variables in Xcal2 as well
  void set_LeptonInterference(TVar::LeptonInterference myLepInterf); // Sets variables in Xcal2 as well
  //
  void set_TempCandidate(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Sets melaCand in Xcal2 to a temporary candidate, without pushing this candidate to candList of Xcal2 for storage and deletion at a later stage
  //
  void set_RenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf); // Sets variables exclusive to Xcal2
  void set_LHAgrid(const char* path, int pdfmember=0); // Sets variable exclusive to Xcal2
  void set_PrimaryHiggsMass(double mh);
  void set_CandidateDecayMode(TVar::CandidateDecayMode mode); // Sets variables exclusive to Xcal2
  void set_CurrentCandidateFromIndex(unsigned int icand); // Sets variables exclusive to Xcal2
  void set_CurrentCandidate(MELACandidate* cand); // Sets variables exclusive to Xcal2
  void set_InputEvent(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Sets variables exclusive to Xcal2
  void append_TopCandidate(SimpleParticleCollection_t* TopDaughters); // Sets variable exclusive to Xcal2
  //
  void set_mHiggs(double mh_, int index); // Does not set any variables in Xcal2!
  void set_wHiggs(double gah_, int index); // Does not set any variables in Xcal2!
  void set_mHiggs_wHiggs(double mh_, double gah_, int index); // Does not set any variables in Xcal2!
  //

  // Reset-functions
  void reset_Mass(double inmass, int ipart);
  void reset_Width(double inmass, int ipart);
  void reset_QuarkMasses();
  void reset_MCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  void resetPerEvent(); // Resets variables and owned objects that are per-event
  void reset_InputEvent(); // Resets all candidates in Xcal2, to be called at the end of each event after all computations are done

  void set_SpinZeroCouplings(
    double selfDHggcoupl[nSupportedHiggses][SIZE_HGG][2],
    double selfDHg4g4coupl[nSupportedHiggses][SIZE_HGG][2],
    double selfDHqqcoupl[nSupportedHiggses][SIZE_HQQ][2],
    double selfDHbbcoupl[nSupportedHiggses][SIZE_HQQ][2],
    double selfDHttcoupl[nSupportedHiggses][SIZE_HQQ][2],
    double selfDHb4b4coupl[nSupportedHiggses][SIZE_HQQ][2],
    double selfDHt4t4coupl[nSupportedHiggses][SIZE_HQQ][2],
    double selfDHzzcoupl[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl[nSupportedHiggses][SIZE_HVV][2],
    double selfDHzzLambda_qsq[nSupportedHiggses][SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ],
    double selfDHwwLambda_qsq[nSupportedHiggses][SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ],
    int selfDHzzCLambda_qsq[nSupportedHiggses][SIZE_HVV_CQSQ],
    int selfDHwwCLambda_qsq[nSupportedHiggses][SIZE_HVV_CQSQ],
    bool diffHWW = false
    );
  void set_SpinZeroContact(
    double selfDHzzpcoupl[SIZE_HVV][2],
    double selfDHzpzpcoupl[SIZE_HVV][2],
    double selfDHwwpcoupl[SIZE_HVV][2],
    double selfDHwpwpcoupl[SIZE_HVV][2]
    );
  void set_SpinOneCouplings(
    double selfDZqqcoupl[SIZE_ZQQ][2],
    double selfDZvvcoupl[SIZE_ZVV][2]
    );
  void set_SpinTwoCouplings(
    double selfDGqqcoupl[SIZE_GQQ][2],
    double selfDGggcoupl[SIZE_GGG][2],
    double selfDGvvcoupl[SIZE_GVV][2]
    );
  void set_SpinTwoContact(
    double selfDGvvpcoupl[SIZE_GVV][2],
    double selfDGvpvpcoupl[SIZE_GVV][2]
  );
  void set_VprimeContactCouplings(
    double selfDZpffcoupl[SIZE_Vpff][2],
    double selfDWpffcoupl[SIZE_Vpff][2],
    double M_Zprime,
    double Ga_Zprime,
    double M_Wprime,
    double Ga_Wprime
  );

  // Compute four-momenta from angles - not cos(theta) - only 
  std::vector<TLorentzVector> Calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi);

  // Get-functions
  MelaIO* get_IORecord();
  double get_PrimaryMass(int ipart);
  double get_PrimaryHiggsMass(){ return get_PrimaryMass(25); }
  double get_PrimaryWidth(int ipart);
  double get_HiggsWidthAtPoleMass(double mass);
  MELACandidate* get_CurrentCandidate();
  int get_CurrentCandidateIndex();
  int get_NCandidates();
  std::vector<MELATopCandidate*>* get_TopCandidateCollection(); // Just so that the user can set MELATopCandidate::passSelection=true or false to omit some tops, in case tere are more than two

protected:
  
  TVar::Process processModel;
  TVar::MatrixElement processME;
  TVar::Production processProduction;
  TVar::VerbosityLevel processVerbosity;
  TVar::LeptonInterference processLeptonInterference;

  double EBEAM;
  double mHiggs[nSupportedHiggses];
  double wHiggs[nSupportedHiggses];
  TEvtProb Xcal2;

  SpinZeroCouplings* selfD_SpinZeroCouplings;
  SpinOneCouplings* selfD_SpinOneCouplings;
  SpinTwoCouplings* selfD_SpinTwoCouplings;
  VprimeCouplings* selfD_VprimeCouplings;

  MELACandidate* melaCand; // Pointer to current candidate object of Xcal2
  std::vector<MELAParticle*> tmpPartList; // Vector of pointers to the owned, temporary MELAParticles
  // Having a temporary top candidate list does not make much sense at the moment
  //std::vector<MELATopCandidate*> tmpTopCandList; // Vector of pointers to the owned, temporary MELATopCandidates
  std::vector<MELACandidate*> tmpCandList; // Vector of pointers to the owned, temporary MELACandidates

  // Constructor wrapper
  void build();

};
#endif
