//-----------------------------------------------------------------------------
//
// Class EventProb Module
//
//   EventProb Module
//
// March 21 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)
//-----------------------------------------------------------------------------

#include "TEvtProb.hh"
#include "MELAStreamHelpers.hh"


ClassImp(TEvtProb)


using namespace std;
using namespace TUtil;
using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;


//-----------------------------------------------------------------------------
// Constructors and Destructor
//-----------------------------------------------------------------------------
TEvtProb::TEvtProb(
  const char* pathtoXSW, double ebeam, const char* pathtoPDFSet, int PDFMember, TVar::VerbosityLevel verbosity_
  ) :
  pathtoPDFSet_(pathtoPDFSet),
  PDFMember_(PDFMember),
  verbosity(verbosity_),
  EBEAM(ebeam),
  myCSW_(pathtoXSW)
{
  if (verbosity>=TVar::DEBUG) MELAout << "Begin TEvtProb constructor" << endl;

  /***** Build everything except MELAHXSWidth *****/
  Build();

  if (verbosity>=TVar::DEBUG) MELAout << "End TEvtProb constructor" << endl;
}
TEvtProb::TEvtProb(const TEvtProb& other) :
pathtoPDFSet_(other.pathtoPDFSet_),
PDFMember_(other.PDFMember_),
verbosity(other.verbosity),
EBEAM(other.EBEAM),
myCSW_(other.myCSW_)
{

}
void TEvtProb::Build(){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin TEvtProb::Build" << endl;

  /***** Initialize lepton interference scheme *****/
  SetLeptonInterf(TVar::DefaultLeptonInterf);

  /***** Initialize MCFM *****/
  InitializeMCFM();

  /***** Initialize JHUGen *****/
  InitializeJHUGen(pathtoPDFSet_, PDFMember_);

  /***** Cross-initializations *****/
  CrossInitialize();

  /***** Initialize input event properties *****/
  ResetInputEvent();
  SetCandidateDecayMode(TVar::CandidateDecay_ZZ);

  if (verbosity>=TVar::DEBUG) MELAout << "End TEvtProb::Build" << endl;
}


TEvtProb::~TEvtProb(){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin TEvtProb destructor" << endl;

  ResetInputEvent();

  if (verbosity>=TVar::DEBUG) MELAout << "End TEvtProb destructor" << endl;
}

void TEvtProb::InitializeMCFM(){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin TEvtProb::InitializeMCFM" << endl;

  mcfm_init_((char *)"input.DAT", (char *)"./");
  if (verbosity>=TVar::DEBUG) MELAout << "TEvtProb::TEvtProb: mcfm_init successful" << endl;
  energy_.sqrts = 2.*EBEAM;
  ResetQuarkMasses();
  ResetMCFM_EWKParameters(
    1.16639E-05, // Gf
    1./128., // alpha_EW
    80.399, // mW
    91.1876, // mZ
    0.23119, // xw=sin**2(thetaW)
    3 // MCFM EW scheme
    );
  qlinit_();
  // First resonance constant parameters
  spinzerohiggs_anomcoupl_.LambdaBSM=1000;
  spinzerohiggs_anomcoupl_.Lambda_Q=10000;
  spinzerohiggs_anomcoupl_.Lambda_z1=10000;
  spinzerohiggs_anomcoupl_.Lambda_z2=10000;
  spinzerohiggs_anomcoupl_.Lambda_z3=10000;
  spinzerohiggs_anomcoupl_.Lambda_z4=10000;
  spinzerohiggs_anomcoupl_.Lambda_z11=100; spinzerohiggs_anomcoupl_.Lambda_z21=100; spinzerohiggs_anomcoupl_.Lambda_z31=100; spinzerohiggs_anomcoupl_.Lambda_z41=100;
  spinzerohiggs_anomcoupl_.Lambda_z12=100; spinzerohiggs_anomcoupl_.Lambda_z22=100; spinzerohiggs_anomcoupl_.Lambda_z32=100; spinzerohiggs_anomcoupl_.Lambda_z42=100;
  spinzerohiggs_anomcoupl_.Lambda_z10=100; spinzerohiggs_anomcoupl_.Lambda_z20=100; spinzerohiggs_anomcoupl_.Lambda_z30=100; spinzerohiggs_anomcoupl_.Lambda_z40=100;
  spinzerohiggs_anomcoupl_.Lambda_zgs1=10000;
  spinzerohiggs_anomcoupl_.Lambda_w1=10000;
  spinzerohiggs_anomcoupl_.Lambda_w2=10000;
  spinzerohiggs_anomcoupl_.Lambda_w3=10000;
  spinzerohiggs_anomcoupl_.Lambda_w4=10000;
  spinzerohiggs_anomcoupl_.Lambda_w11=100; spinzerohiggs_anomcoupl_.Lambda_w21=100; spinzerohiggs_anomcoupl_.Lambda_w31=100; spinzerohiggs_anomcoupl_.Lambda_w41=100;
  spinzerohiggs_anomcoupl_.Lambda_w12=100; spinzerohiggs_anomcoupl_.Lambda_w22=100; spinzerohiggs_anomcoupl_.Lambda_w32=100; spinzerohiggs_anomcoupl_.Lambda_w42=100;
  spinzerohiggs_anomcoupl_.Lambda_w10=100; spinzerohiggs_anomcoupl_.Lambda_w20=100; spinzerohiggs_anomcoupl_.Lambda_w30=100; spinzerohiggs_anomcoupl_.Lambda_w40=100;
  // Second resonance constant parameters
  spinzerohiggs_anomcoupl_.Lambda2BSM=1000;
  spinzerohiggs_anomcoupl_.Lambda2_Q=10000;
  spinzerohiggs_anomcoupl_.Lambda2_z1=10000;
  spinzerohiggs_anomcoupl_.Lambda2_z2=10000;
  spinzerohiggs_anomcoupl_.Lambda2_z3=10000;
  spinzerohiggs_anomcoupl_.Lambda2_z4=10000;
  spinzerohiggs_anomcoupl_.Lambda2_z11=100; spinzerohiggs_anomcoupl_.Lambda2_z21=100; spinzerohiggs_anomcoupl_.Lambda2_z31=100; spinzerohiggs_anomcoupl_.Lambda2_z41=100;
  spinzerohiggs_anomcoupl_.Lambda2_z12=100; spinzerohiggs_anomcoupl_.Lambda2_z22=100; spinzerohiggs_anomcoupl_.Lambda2_z32=100; spinzerohiggs_anomcoupl_.Lambda2_z42=100;
  spinzerohiggs_anomcoupl_.Lambda2_z10=100; spinzerohiggs_anomcoupl_.Lambda2_z20=100; spinzerohiggs_anomcoupl_.Lambda2_z30=100; spinzerohiggs_anomcoupl_.Lambda2_z40=100;
  spinzerohiggs_anomcoupl_.Lambda2_zgs1=10000;
  spinzerohiggs_anomcoupl_.Lambda2_w1=10000;
  spinzerohiggs_anomcoupl_.Lambda2_w2=10000;
  spinzerohiggs_anomcoupl_.Lambda2_w3=10000;
  spinzerohiggs_anomcoupl_.Lambda2_w4=10000;
  spinzerohiggs_anomcoupl_.Lambda2_w11=100; spinzerohiggs_anomcoupl_.Lambda2_w21=100; spinzerohiggs_anomcoupl_.Lambda2_w31=100; spinzerohiggs_anomcoupl_.Lambda2_w41=100;
  spinzerohiggs_anomcoupl_.Lambda2_w12=100; spinzerohiggs_anomcoupl_.Lambda2_w22=100; spinzerohiggs_anomcoupl_.Lambda2_w32=100; spinzerohiggs_anomcoupl_.Lambda2_w42=100;
  spinzerohiggs_anomcoupl_.Lambda2_w10=100; spinzerohiggs_anomcoupl_.Lambda2_w20=100; spinzerohiggs_anomcoupl_.Lambda2_w30=100; spinzerohiggs_anomcoupl_.Lambda2_w40=100;
  // Switches for spin-0 Higgs couplings
  spinzerohiggs_anomcoupl_.channeltoggle_stu=2; spinzerohiggs_anomcoupl_.vvhvvtoggle_vbfvh=2;
  spinzerohiggs_anomcoupl_.AnomalCouplPR=1; spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
  // Constant parameters for all processes
  qlfirst_.qlfirst=false;

  if (verbosity>=TVar::DEBUG) MELAout << "End TEvtProb::InitializeMCFM" << endl;
}
void TEvtProb::InitializeJHUGen(const char* pathtoPDFSet, int PDFMember){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin TEvtProb::InitializeJHUGen" << endl;

  InitJHUGenMELA(pathtoPDFSet, PDFMember);

  if (verbosity>=TVar::DEBUG) MELAout << "End TEvtProb::InitializeJHUGen" << endl;
}
void TEvtProb::CrossInitialize(){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin TEvtProb::CrossInitialize" << endl;

  /**** Initialize alphaS *****/
  const double GeV=1./100.;
  double init_Q;
  int init_nl, init_nf;
  __modjhugenmela_MOD_getpdfconstants(&init_Q, &init_nl, &init_nf);
  init_Q /= GeV;
  if (verbosity>=TVar::DEBUG) MELAout << "TEvtProb::TEvtProb: Initializing the PDF with initial Q=" << init_Q << ", nloops=" << init_nl << ", nf=" << init_nf << endl;
  SetAlphaS(init_Q, init_Q, 1, 1, init_nl, init_nf, "cteq6_l"); // "cteq6_l" is just some dummy variable

  /**** Initialize MCFM CKM initializers in cabib.f using JHUGen defaults *****/
  int i, j;
  i=2; j=1;
  cabib_.Vud = __modparameters_MOD_ckmbare(&i, &j);
  i=2; j=3;
  cabib_.Vus = __modparameters_MOD_ckmbare(&i, &j);
  i=2; j=5;
  cabib_.Vub = __modparameters_MOD_ckmbare(&i, &j);
  i=4; j=1;
  cabib_.Vcd = __modparameters_MOD_ckmbare(&i, &j);
  i=4; j=3;
  cabib_.Vcs = __modparameters_MOD_ckmbare(&i, &j);
  i=4; j=5;
  cabib_.Vcb = __modparameters_MOD_ckmbare(&i, &j);

  /***** Initialize schemes *****/
  ResetCouplings();
  ResetRenFacScaleMode();
  SetPrimaryHiggsMass(125.); // Should come after InitializeMCFM and InitializeJHUGen

  if (verbosity>=TVar::DEBUG) MELAout << "End TEvtProb::CrossInitialize" << endl;
}

// Set NNPDF driver path
void TEvtProb::Set_LHAgrid(const char* path, int pdfmember){
  char path_nnpdf_c[200];
  sprintf(path_nnpdf_c, "%s", path);
  int pathLength = strlen(path_nnpdf_c);
  nnpdfdriver_(path_nnpdf_c, &pathLength);
  nninitpdf_(&pdfmember);
}
void TEvtProb::SetProcess(TVar::Process proc, TVar::MatrixElement me, TVar::Production prod){
  matrixElement = me;
  production = prod;
  // In case s-channel processes are passed for JHUGen ME, flip them back to JHUGen-specific productions.
  if (matrixElement==TVar::JHUGen){
    if (production==TVar::Had_ZH_S) production=TVar::Had_ZH;
    else if (production==TVar::Had_WH_S) production=TVar::Had_WH;
    else if (production==TVar::Lep_ZH_S) production=TVar::Lep_ZH;
    else if (production==TVar::Lep_WH_S) production=TVar::Lep_WH;
    else if (production==TVar::JJVBF_S) production=TVar::JJVBF;
    else if (production==TVar::JJQCD_S) production=TVar::JJQCD;
  }
  process = proc;
}
void TEvtProb::SetVerbosity(TVar::VerbosityLevel tmp){ verbosity = tmp; }
void TEvtProb::SetLeptonInterf(TVar::LeptonInterference tmp){ leptonInterf = tmp; }
void TEvtProb::SetCandidateDecayMode(TVar::CandidateDecayMode mode){ PDGHelpers::setCandidateDecayMode(mode); }
void TEvtProb::SetRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf){
  event_scales.renomalizationScheme = renormalizationSch;
  event_scales.factorizationScheme = factorizationSch;
  event_scales.ren_scale_factor = ren_sf;
  event_scales.fac_scale_factor = fac_sf;
}
void TEvtProb::AllowSeparateWWCouplings(bool doAllow){ SetJHUGenDistinguishWWCouplings(doAllow); selfDSpinZeroCoupl.allow_WWZZSeparation(doAllow); }
void TEvtProb::SetPrimaryHiggsMass(double mass){ PrimaryHMass=mass; SetHiggsMass(PrimaryHMass, -1., -1); }
void TEvtProb::SetHiggsMass(double mass, double wHiggs, int whichResonance){
  // Regular, first resonance
  if (whichResonance==1 || whichResonance==-1){
    if (mass<0.){
      _hmass = -1;
      _hwidth = 0;
    }
    else if (wHiggs<0.){
      _hmass = mass;
      _hwidth = myCSW_.HiggsWidth(_hmass);
    }
    else{
      _hmass = mass;
      _hwidth = wHiggs;
    }
    masses_mcfm_.hmass = _hmass;
    masses_mcfm_.hwidth = _hwidth;

    if (_hmass<0.) SetJHUGenHiggsMassWidth(0., _hwidth);
    else SetJHUGenHiggsMassWidth(_hmass, _hwidth);
  }

  // Second resonance
  if (whichResonance==2){
    if (mass<0.){
      _h2mass = -1;
      _h2width = 0;
    }
    else if (wHiggs<0.){
      _h2mass = mass;
      _h2width = myCSW_.HiggsWidth(_h2mass);
    }
    else{
      _h2mass = mass;
      _h2width = wHiggs;
    }
    spinzerohiggs_anomcoupl_.h2mass = _h2mass;
    spinzerohiggs_anomcoupl_.h2width = _h2width;
    //SetJHUGenHiggsMassWidth(_h2mass, _h2width); // Second resonance is not implemented in JHUGen yet.
  }
  else if (whichResonance==-1){
    _h2mass = -1;
    _h2width = 0;
    spinzerohiggs_anomcoupl_.h2mass = _h2mass;
    spinzerohiggs_anomcoupl_.h2width = _h2width;
    //SetJHUGenHiggsMassWidth(_h2mass, _h2width); // Second resonance is not implemented in JHUGen yet.
  }

  if (verbosity>=TVar::DEBUG) MELAout
    << "TEvtProb::SetHiggsMass(" << mass << ", " << wHiggs << ", " << whichResonance << "):\n"
    << '\t' << "hmass: " << _hmass << ", " << _hwidth << '\n'
    << '\t' << "h2mass: " << _h2mass << ", " << _h2width << '\n'
    << '\t' << "MCFM hmass: " << masses_mcfm_.hmass << ", " << masses_mcfm_.hwidth << '\n'
    << '\t' << "MCFM h2mass: " << spinzerohiggs_anomcoupl_.h2mass << ", " << spinzerohiggs_anomcoupl_.h2width << endl;
}
void TEvtProb::SetInputEvent(
  SimpleParticleCollection_t* pDaughters,
  SimpleParticleCollection_t* pAssociated,
  SimpleParticleCollection_t* pMothers,
  bool isGen
  ){
  MELACandidate* cand = ConvertVectorFormat(
    pDaughters,
    pAssociated,
    pMothers,
    isGen,
    &particleList, &candList // push_back is done automatically
    );
  if (cand!=0) melaCand=cand;
}
void TEvtProb::AppendTopCandidate(SimpleParticleCollection_t* TopDaughters){
  if (!CheckInputPresent()){
    MELAerr << "TEvtProb::AppendTopCandidate: No MELACandidates are present to append this top!" << endl;
    return;
  }
  MELATopCandidate* cand = ConvertTopCandidate(
    TopDaughters,
    &particleList, &topCandList // push_back is done automatically
    );
  if (cand!=0) melaCand->addAssociatedTops(cand);
}
void TEvtProb::SetRcdCandPtr(){ RcdME.melaCand = melaCand; }
void TEvtProb::SetCurrentCandidateFromIndex(unsigned int icand){
  if (candList.size()>icand) melaCand = candList.at(icand);
  else MELAerr << "TEvtProb::SetCurrentCandidateFromIndex: icand=" << icand << ">=candList.size()=" << candList.size() << endl;
}
void TEvtProb::SetCurrentCandidate(MELACandidate* cand){
  melaCand = cand;
  if (verbosity>=TVar::INFO && melaCand==0) MELAout << "TEvtProb::SetCurrentCandidate: BE CAREFUL! melaCand==0!" << endl;
  if (verbosity>=TVar::INFO && GetCurrentCandidateIndex()<0) MELAout << "TEvtProb::SetCurrentCandidate: The current candidate is not in the list of candidates. It is the users' responsibility to delete this candidate and all of its associated particles." << endl;
}

// Reset functions
void TEvtProb::ResetIORecord(){ RcdME.reset(); }
void TEvtProb::ResetRenFacScaleMode(){ SetRenFacScaleMode(TVar::DefaultScaleScheme, TVar::DefaultScaleScheme, 0.5, 0.5); }
void TEvtProb::ResetMass(double inmass, int ipart){ TUtil::SetMass(inmass, ipart); }
void TEvtProb::SetZprimeMassWidth(double inmass, double inwidth){ this->ResetMass(inmass, 32); this->ResetWidth(inwidth, 32); }
void TEvtProb::SetWprimeMassWidth(double inmass, double inwidth){ this->ResetMass(inmass, 34); this->ResetWidth(inwidth, 34); }
void TEvtProb::ResetWidth(double inwidth, int ipart){ TUtil::SetDecayWidth(inwidth, ipart); }
void TEvtProb::ResetQuarkMasses(){
  ResetMass(1e-3, 1); // d
  ResetMass(5e-3, 2); // u
  ResetMass(1e-1, 3); // s
  ResetMass(1.275, 4); // c
  ResetMass(4.75, 5); // b
  ResetMass(173.2, 6); // t
  ResetMass(1e5, 7); // bprime
  ResetMass(1e5, 8); // tprime
  //coupling_(); // Already called at b and t twice!
}
void TEvtProb::ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme){ SetEwkCouplingParameters(ext_Gf, ext_aemmz, ext_mW, ext_mZ, ext_xW, ext_ewscheme); }
void TEvtProb::ResetCouplings(){
  selfDSpinZeroCoupl.reset();
  selfDSpinOneCoupl.reset();
  selfDSpinTwoCoupl.reset();
  AllowSeparateWWCouplings(false);
  ResetAmplitudeIncludes();
}
void TEvtProb::ResetInputEvent(){
  RcdME.melaCand = 0;
  melaCand = 0;

  // Clear bookkeeping objects
  // Clean MELACandidates first since they contain all other objects
  for (MELACandidate* tmpPart:candList) delete tmpPart;
  candList.clear();
  // Clean MELATopCandidates next since they contain MELAParticles
  for (MELATopCandidate* tmpPart:topCandList) delete tmpPart;
  topCandList.clear();
  // Clean all remaining MELAPArticles
  for (MELAParticle* tmpPart:particleList) delete tmpPart;
  particleList.clear();
}

// Get-functions
MELAHXSWidth const* TEvtProb::GetHXSWidthEstimator() const{ return &myCSW_; }
SpinZeroCouplings* TEvtProb::GetSelfDSpinZeroCouplings(){ return selfDSpinZeroCoupl.getRef(); }
SpinOneCouplings* TEvtProb::GetSelfDSpinOneCouplings(){ return selfDSpinOneCoupl.getRef(); }
SpinTwoCouplings* TEvtProb::GetSelfDSpinTwoCouplings(){ return selfDSpinTwoCoupl.getRef(); }
VprimeCouplings* TEvtProb::GetSelfDVprimeCouplings(){ return selfDVprimeCoupl.getRef(); }
double TEvtProb::GetPrimaryHiggsMass(){ return PrimaryHMass; }
double TEvtProb::GetPrimaryMass(int ipart){
  if (PDGHelpers::isAHiggs(ipart)) return GetPrimaryHiggsMass();
  else return TUtil::GetMass(ipart);
}
double TEvtProb::GetPrimaryWidth(int ipart){
  if (PDGHelpers::isAHiggs(ipart)) return myCSW_.HiggsWidth(GetPrimaryHiggsMass());
  else return TUtil::GetDecayWidth(ipart);
}
double TEvtProb::GetHiggsWidthAtPoleMass(double mass){
  if (mass>0.) return myCSW_.HiggsWidth(mass);
  else return -1.;
}
MelaIO* TEvtProb::GetIORecord(){ return RcdME.getRef(); }
MELACandidate* TEvtProb::GetCurrentCandidate(){ return melaCand; }
int TEvtProb::GetCurrentCandidateIndex(){
  if (melaCand==0) return -1;
  for (unsigned int icand=0; icand<candList.size(); icand++){
    if (candList.at(icand)==melaCand) return (int)icand;
  }
  return -1;
}
int TEvtProb::GetNCandidates(){ return (static_cast<int>(candList.size())); }
std::vector<MELATopCandidate*>* TEvtProb::GetTopCandidates(){ return &topCandList; }


// Check/test functions
bool TEvtProb::CheckInputPresent(){
  if (melaCand==0){
    MELAerr
      << "TEvtProb::CheckInputPresent: melaCand==" << melaCand << " is nullPtr!"
      << endl;
    if (candList.size()==0) return false;
    else{
      SetCurrentCandidateFromIndex(candList.size()-1);
      MELAerr << "TEvtProb::CheckInputPresent: melaCand now points to the latest candidate (cand" << (candList.size()-1) << ")" << endl;
    }
  }
  SetRcdCandPtr(); // If input event is present, set the RcdME pointer to melaCand
  return true;
}


// ME functions
// Cross-section calculations for H + 0 jet
double TEvtProb::XsecCalc_XVV(){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin XsecCalc_XVV" << endl;
  double dXsec=0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  bool useMCFM = matrixElement == TVar::MCFM;
  bool calculateME=false;
  bool needBSMHiggs=false;
  if (useMCFM){
    if (verbosity>=TVar::DEBUG) MELAout << "TEvtProb::XsecCalc_XVV: Try MCFM" << endl;
    needBSMHiggs = CheckSelfDCouplings_Hgg() || CheckSelfDCouplings_Htt() || CheckSelfDCouplings_Hbb() || CheckSelfDCouplings_HVV();
    if (needBSMHiggs) SetLeptonInterf(TVar::InterfOn); // All anomalous coupling computations have to use lepton interference

    calculateME = (
      production==TVar::ZZGG ||
      production==TVar::ZZQQB ||
      production==TVar::ZZQQB_STU ||
      production==TVar::ZZQQB_S ||
      production==TVar::ZZQQB_TU ||
      production==TVar::ZZINDEPENDENT ||
      (production == TVar::JJQCD && process == TVar::bkgZJets) // Allow Z+jets to be computed here.
      );
    if (calculateME){
      SetMCFMSpinZeroCouplings(needBSMHiggs, &selfDSpinZeroCoupl, false);
      dXsec = SumMatrixElementPDF(process, production, matrixElement, leptonInterf, &event_scales, &RcdME, EBEAM, verbosity);
    }
    else if (verbosity>=TVar::INFO) MELAout << "TEvtProb::XsecCalc_XVV: MCFM_chooser failed to determine the process configuration." << endl;
  }
  else if (matrixElement == TVar::JHUGen){
    if (verbosity>=TVar::DEBUG) MELAout << "TEvtProb::XsecCalc_XVV::Try JHUGen" << endl;
    AllowSeparateWWCouplings(false); // HZZ couplings are used for both in spin-0
    // all the possible couplings
    double Hggcoupl[SIZE_HGG][2] ={ { 0 } };
    double Hvvcoupl[SIZE_HVV][2] ={ { 0 } };
    double HvvLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ] ={ { 0 } };
    int HvvCLambda_qsq[SIZE_HVV_CQSQ] ={ 0 };

    double Hvvpcoupl[SIZE_HVV][2] = { { 0 } };
    double Hvpvpcoupl[SIZE_HVV][2] = { { 0 } };

    double Zqqcoupl[SIZE_ZQQ][2] ={ { 0 } };
    double Zvvcoupl[SIZE_ZVV][2] ={ { 0 } };

    double Gqqcoupl[SIZE_GQQ][2] ={ { 0 } };
    double Gggcoupl[SIZE_GGG][2] ={ { 0 } };
    double Gvvcoupl[SIZE_GVV][2] ={ { 0 } };
    double Gvvpcoupl[SIZE_GVV][2] ={ { 0 } };
    double Gvpvpcoupl[SIZE_GVV][2] ={ { 0 } };

    double Zpffcoupl[SIZE_Vpff][2] ={ { 0 } };
    double Wpffcoupl[SIZE_Vpff][2] ={ { 0 } };
    double M_Zprime = -1;
    double Ga_Zprime = 0;
    double M_Wprime = -1;
    double Ga_Wprime = 0;

    //
    // set spin 0 default numbers
    //
    // By default set the Spin 0 couplings for SM case (0+m)
    Hggcoupl[gHIGGS_GG_2][0]=1.0;
    for (int ic=0; ic<SIZE_HVV_LAMBDAQSQ; ic++){ for (int ik=0; ik<SIZE_HVV_CQSQ; ik++) HvvLambda_qsq[ic][ik]=100.; }
    //
    // set spin 1 default numbers
    //
    Zqqcoupl[gZPRIME_QQ_LEFT][0]=1.0;
    Zqqcoupl[gZPRIME_QQ_RIGHT][0]=1.0;
    //
    // set spin 2 default numbers
    //
    Gqqcoupl[gGRAVITON_QQ_LEFT][0]=1.0;
    Gqqcoupl[gGRAVITON_QQ_RIGHT][0]=1.0;

    bool isSpinZero = false;
    bool isSpinOne = false;
    bool isSpinTwo = false;

    // 0m+
    if (process == TVar::HSMHiggs){
      Hvvcoupl[gHIGGS_VV_1][0]=1.;
      isSpinZero = true;
    }
    // 0+L1
    else if (process == TVar::H0_g1prime2){
      Hvvcoupl[gHIGGS_VV_1_PRIME2][0] = 1.;
      isSpinZero = true;
    }
    // 0h+
    else if (process == TVar::H0hplus){
      Hvvcoupl[gHIGGS_VV_2][0] = 1.;
      isSpinZero = true;
    }
    // 0-
    else if (process == TVar::H0minus){
      Hvvcoupl[gHIGGS_VV_4][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_Zgsg1prime2){
      Hvvcoupl[gHIGGS_ZA_1_PRIME2][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_Zgs){
      Hvvcoupl[gHIGGS_ZA_2][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_Zgs_PS){
      Hvvcoupl[gHIGGS_ZA_4][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_gsgs){
      Hvvcoupl[gHIGGS_AA_2][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_gsgs_PS){
      Hvvcoupl[gHIGGS_AA_4][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::SelfDefine_spin0){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_HGG; i++) Hggcoupl[i][j] = (selfDSpinZeroCoupl.Hggcoupl)[i][j];
        for (int i=0; i<SIZE_HVV; i++){
          Hvvcoupl[i][j] = (selfDSpinZeroCoupl.Hzzcoupl)[i][j];

          Hvvpcoupl[i][j] = (selfDSpinZeroCoupl.Hzzpcoupl)[i][j];
          Hvpvpcoupl[i][j] = (selfDSpinZeroCoupl.Hzpzpcoupl)[i][j];
        }
      }
      for (int j=0; j<SIZE_HVV_CQSQ; j++){
        for (int i=0; i<SIZE_HVV_LAMBDAQSQ; i++) HvvLambda_qsq[i][j] = (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j];
        HvvCLambda_qsq[j] = (selfDSpinZeroCoupl.HzzCLambda_qsq)[j];
      }

      if (verbosity>=TVar::DEBUG_VERBOSE){
        for (int j=0; j<2; j++){
          for (int i=0; i<SIZE_HGG; i++){ if ((selfDSpinZeroCoupl.Hggcoupl)[i][j]!=0.) MELAout << "Hggcoupl[" << i << "][" << j << "] = " << (selfDSpinZeroCoupl.Hggcoupl)[i][j] << endl; }
          for (int i=0; i<SIZE_HVV; i++){
            if ((selfDSpinZeroCoupl.Hzzcoupl)[i][j]!=0.) MELAout << "Hvvcoupl[" << i << "][" << j << "] = " << (selfDSpinZeroCoupl.Hzzcoupl)[i][j] << endl;

            if ((selfDSpinZeroCoupl.Hzzpcoupl)[i][j]!=0.) MELAout << "Hzzpcoupl[" << i << "][" << j << "] = " << (selfDSpinZeroCoupl.Hzzpcoupl)[i][j] << endl;
            if ((selfDSpinZeroCoupl.Hzpzpcoupl)[i][j]!=0.) MELAout << "Hzpzpcoupl[" << i << "][" << j << "] = " << (selfDSpinZeroCoupl.Hzpzpcoupl)[i][j] << endl;
          }
        }
        for (int j=0; j<SIZE_HVV_CQSQ; j++){
          for (int i=0; i<SIZE_HVV_LAMBDAQSQ; i++){ if ((selfDSpinZeroCoupl.HzzLambda_qsq)[i][j]!=0.) MELAout << "HvvLambda_qsq[" << i << "][" << j << "] = " << (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j] << endl; }
          if ((selfDSpinZeroCoupl.HzzCLambda_qsq)[j]!=0.) MELAout << "HvvCLambda_qsq[" << j << "] = " << (selfDSpinZeroCoupl.HzzCLambda_qsq)[j] << endl;
        }
      }

      isSpinZero = true;
    }

    // 1-
    else if (process == TVar::H1minus){
      Zvvcoupl[gZPRIME_VV_1][0]=1.;
      isSpinOne = true;
    }
    // 1+
    else if (process == TVar::H1plus){
      Zvvcoupl[gZPRIME_VV_2][0]=1.;
      isSpinOne = true;
    }
    else if (process == TVar::SelfDefine_spin1){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_ZQQ; i++) Zqqcoupl[i][j] = (selfDSpinOneCoupl.Zqqcoupl)[i][j];
        for (int i=0; i<SIZE_ZVV; i++) Zvvcoupl[i][j] = (selfDSpinOneCoupl.Zvvcoupl)[i][j];
      }
      isSpinOne = true;
    }

    // 2m+
    else if (process == TVar::H2_g1g5){
      Gggcoupl[gGRAVITON_GG_1][0]=1.;
      Gvvcoupl[gGRAVITON_VV_1][0]=1.;
      Gvvcoupl[gGRAVITON_VV_5][0]=1.;
      isSpinTwo = true;
    }
    // 2b+
    else if (process == TVar::H2_g1){
      Gggcoupl[gGRAVITON_GG_1][0]=1.;
      Gvvcoupl[gGRAVITON_VV_1][0]=1.;
      isSpinTwo = true;
    }
    else if (process == TVar::H2_g2){
      Gggcoupl[gGRAVITON_GG_2][0]=1.;
      Gvvcoupl[gGRAVITON_VV_2][0]=1.;
      isSpinTwo = true;
    }
    // 2h3plus
    else if (process == TVar::H2_g3){
      Gggcoupl[gGRAVITON_GG_3][0]=1.;
      Gvvcoupl[gGRAVITON_VV_3][0]=1.;
      isSpinTwo = true;
    }
    // 2h+
    else if (process == TVar::H2_g4){
      Gggcoupl[gGRAVITON_GG_4][0]=1.;
      Gvvcoupl[gGRAVITON_VV_4][0]=1.;
      isSpinTwo = true;
    }
    // 2b+
    else if (process == TVar::H2_g5){
      Gggcoupl[gGRAVITON_GG_1][0]=1.;
      Gvvcoupl[gGRAVITON_VV_5][0]=1.;
      isSpinTwo = true;
    }
    // 2h6+
    else if (process == TVar::H2_g6){
      Gggcoupl[gGRAVITON_GG_1][0]=1.;
      Gvvcoupl[gGRAVITON_VV_6][0]=1.;
      isSpinTwo = true;
    }
    // 2h7plus
    else if (process == TVar::H2_g7){
      Gggcoupl[gGRAVITON_GG_1][0]=1.;
      Gvvcoupl[gGRAVITON_VV_7][0]=1.;
      isSpinTwo = true;
    }
    // 2h-
    else if (process == TVar::H2_g8){
      Gggcoupl[gGRAVITON_GG_5][0]=1.;
      Gvvcoupl[gGRAVITON_VV_8][0]=1.;
      isSpinTwo = true;
    }
    // 2h9minus
    else if (process == TVar::H2_g9){
      Gggcoupl[gGRAVITON_GG_5][0]=1.;
      Gvvcoupl[gGRAVITON_VV_9][0]=1.;
      isSpinTwo = true;
    }
    // 2h10minus
    else if (process == TVar::H2_g10){
      Gggcoupl[gGRAVITON_GG_5][0]=1.;
      Gvvcoupl[gGRAVITON_VV_10][0]=1.;
      isSpinTwo = true;
    }
    else if (process == TVar::SelfDefine_spin2){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_GGG; i++) Gggcoupl[i][j] = (selfDSpinTwoCoupl.Gggcoupl)[i][j];
        for (int i=0; i<SIZE_GQQ; i++) Gqqcoupl[i][j] = (selfDSpinTwoCoupl.Gqqcoupl)[i][j];
        for (int i=0; i<SIZE_GVV; i++){
          Gvvcoupl[i][j] = (selfDSpinTwoCoupl.Gvvcoupl)[i][j];
          Gvvpcoupl[i][j] = (selfDSpinTwoCoupl.Gvvpcoupl)[i][j];
          Gvpvpcoupl[i][j] = (selfDSpinTwoCoupl.Gvpvpcoupl)[i][j];
        }
      }
      isSpinTwo = true;
    }

    // Vprime / contact couplings
    if (process == TVar::SelfDefine_spin0 || process == TVar::SelfDefine_spin2){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_Vpff; i++){
          Zpffcoupl[i][j] = (selfDVprimeCoupl.Zpffcoupl)[i][j];
          Wpffcoupl[i][j] = (selfDVprimeCoupl.Wpffcoupl)[i][j];
        }
      }
      M_Zprime = selfDVprimeCoupl.M_Zprime;
      Ga_Zprime = selfDVprimeCoupl.Ga_Zprime;
      M_Wprime = selfDVprimeCoupl.M_Wprime;
      Ga_Wprime = selfDVprimeCoupl.Ga_Wprime;

      if (verbosity>=TVar::DEBUG_VERBOSE){
        for (int j=0; j<2; j++){
          for (int i=0; i<SIZE_Vpff; i++){
            if ((selfDVprimeCoupl.Zpffcoupl)[i][j]!=0.) MELAout << "Zpffcoupl[" << i << "][" << j << "] = " << (selfDVprimeCoupl.Zpffcoupl)[i][j] << endl;
            if ((selfDVprimeCoupl.Wpffcoupl)[i][j]!=0.) MELAout << "Wpffcoupl[" << i << "][" << j << "] = " << (selfDVprimeCoupl.Wpffcoupl)[i][j] << endl;
          }
        }
        MELAout << "M_Zprime = " << selfDVprimeCoupl.M_Zprime << endl;
        MELAout << "Ga_Zprime = " << selfDVprimeCoupl.Ga_Zprime << endl;
        MELAout << "M_Wprime = " << selfDVprimeCoupl.M_Wprime << endl;
        MELAout << "Ga_Wprime = " << selfDVprimeCoupl.Ga_Wprime << endl;
      }
    }

    if (isSpinZero){
      SetJHUGenSpinZeroGGCouplings(Hggcoupl);
      SetJHUGenSpinZeroVVCouplings(Hvvcoupl, Hvvpcoupl, Hvpvpcoupl, HvvCLambda_qsq, HvvLambda_qsq, false);
      SetJHUGenVprimeContactCouplings(Zpffcoupl, Wpffcoupl);
      SetZprimeMassWidth(M_Zprime, Ga_Zprime);
      SetWprimeMassWidth(M_Wprime, Ga_Wprime);
    }
    else if (isSpinOne) SetJHUGenSpinOneCouplings(Zqqcoupl, Zvvcoupl);
    else if (isSpinTwo){
      SetJHUGenSpinTwoCouplings(Gggcoupl, Gvvcoupl, Gvvpcoupl, Gvpvpcoupl, Gqqcoupl);
      SetJHUGenVprimeContactCouplings(Zpffcoupl, Wpffcoupl);
      SetZprimeMassWidth(M_Zprime, Ga_Zprime);
      SetWprimeMassWidth(M_Wprime, Ga_Wprime);
    }

    if (isSpinZero || isSpinOne || isSpinTwo) dXsec = JHUGenMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, verbosity);
    else MELAerr
      << "TEvtProb::XsecCalc_XVV: JHUGen ME is not spin zero, one or two! The process is described by "
      << "Process: " << TVar::ProcessName(process) << ", Production: " << TVar::ProductionName(production) << ", and ME: " << TVar::MatrixElementName(matrixElement)
      << endl;
  } // end of JHUGen matrix element calculations

  if (verbosity >= TVar::DEBUG) MELAout << "TEvtProb::XsecCalc_XVV: Process " << TVar::ProcessName(process) << " dXsec=" << dXsec << endl;

  if (verbosity>=TVar::DEBUG) MELAout << "TEvtProb::XsecCalc_XVV::Reset couplings" << endl;
  ResetCouplings(); // Should come first
  if (useMCFM){ // Set defaults. Should come next...
    if (needBSMHiggs) SetLeptonInterf(TVar::DefaultLeptonInterf);
    if (calculateME) SetMCFMSpinZeroCouplings(false, &selfDSpinZeroCoupl, true); // ... because of this!
  }
  ResetRenFacScaleMode();
  if (verbosity>=TVar::DEBUG) MELAout << "End XsecCalc_XVV" << endl;
  return dXsec;
}

// Cross-section calculations for H(->VV) + 2 jets
double TEvtProb::XsecCalc_VVXVV(){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin XsecCalc_VVXVV" << endl;
  double dXsec=0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  bool useMCFM = matrixElement == TVar::MCFM;
  bool needBSMHiggs=false;
  bool calculateME=false;
  if (useMCFM){
    if (verbosity>=TVar::DEBUG) MELAout << "TEvtProb::XsecCalc_VVXVV: Try MCFM" << endl;
    needBSMHiggs = CheckSelfDCouplings_HVV();
    if (needBSMHiggs) SetLeptonInterf(TVar::InterfOn); // All anomalous coupling computations have to use lepton interference

    calculateME = (
      production==TVar::Had_WH || production==TVar::Had_ZH || production==TVar::Lep_WH || production==TVar::Lep_ZH || production==TVar::JJVBF || production==TVar::JJEW || production==TVar::JJQCD || production==TVar::JJEWQCD
      || production==TVar::Had_WH_S || production==TVar::Had_ZH_S || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S || production==TVar::JJVBF_S || production==TVar::JJEW_S/* || production==TVar::JJQCD_S*/ || production==TVar::JJEWQCD_S
      || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU || production==TVar::JJVBF_TU || production==TVar::JJEW_TU/* || production==TVar::JJQCD_TU*/ || production==TVar::JJEWQCD_TU
      );
    if (calculateME){
      SetMCFMSpinZeroCouplings(needBSMHiggs, &selfDSpinZeroCoupl, false);
      dXsec = SumMatrixElementPDF(process, production, matrixElement, leptonInterf, &event_scales, &RcdME, EBEAM, verbosity);
    }
    else if (verbosity>=TVar::INFO) MELAout << "TEvtProb::XsecCalc_VVXVV: MCFM_chooser failed to determine the process configuration." << endl;
  }
  else MELAerr << "TEvtProb::XsecCalc_VVXVV: Non-MCFM Mes are not supported." << endl;

  if (verbosity >= TVar::DEBUG) MELAout
    << "Process " << TVar::ProcessName(process)
    << " TEvtProb::XsecCalc_VVXVV: dXsec=" << dXsec
    << endl;

  ResetCouplings(); // Should come first
  if (useMCFM){ // Set defaults. Should come next...
    if (needBSMHiggs) SetLeptonInterf(TVar::DefaultLeptonInterf);
    if (calculateME) SetMCFMSpinZeroCouplings(false, &selfDSpinZeroCoupl, true); // ... because of this!
  }
  ResetRenFacScaleMode();
  if (verbosity>=TVar::DEBUG) MELAout << "End XsecCalc_VVXVV" << endl;
  return dXsec;
}

// Cross-section calculations for H + 2 jets
double TEvtProb::XsecCalcXJJ(){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin XsecCalc_XJJ" << endl;
  if (matrixElement == TVar::MCFM) return XsecCalc_VVXVV();
  double dXsec = 0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  if (matrixElement == TVar::JHUGen){
    if (production == TVar::JJQCD){
      double Hggcoupl[SIZE_HGG][2] ={ { 0 } };
      if (process == TVar::HSMHiggs) Hggcoupl[gHIGGS_GG_2][0] = 1.;
      else if (process == TVar::H0minus) Hggcoupl[gHIGGS_GG_4][0] = 1.;
      else if (process == TVar::SelfDefine_spin0){
        for (int j=0; j<2; j++){
          for (int i=0; i<SIZE_HGG; i++) Hggcoupl[i][j] = (selfDSpinZeroCoupl.Hggcoupl)[i][j];
        }
      }
      SetJHUGenSpinZeroGGCouplings(Hggcoupl);
    }
    else if (production == TVar::JJVBF){
      double Hzzcoupl[SIZE_HVV][2] ={ { 0 } };
      double Hwwcoupl[SIZE_HVV][2] ={ { 0 } };
      double HzzLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ] ={ { 0 } };
      double HwwLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ] ={ { 0 } };
      int HzzCLambda_qsq[SIZE_HVV_CQSQ] ={ 0 };
      int HwwCLambda_qsq[SIZE_HVV_CQSQ] ={ 0 };
      double Hzzpcoupl[SIZE_HVV][2] = { { 0 } };
      double Hzpzpcoupl[SIZE_HVV][2] = { { 0 } };
      double Zpffcoupl[SIZE_Vpff][2] = { { 0 } };
      double Hwwpcoupl[SIZE_HVV][2] = { { 0 } };
      double Hwpwpcoupl[SIZE_HVV][2] = { { 0 } };
      double Wpffcoupl[SIZE_Vpff][2] = { { 0 } };
      double M_Zprime = -1;
      double Ga_Zprime = 0;
      double M_Wprime = -1;
      double Ga_Wprime = 0;

      for (int ic=0; ic<SIZE_HVV_LAMBDAQSQ; ic++){ for (int ik=0; ik<SIZE_HVV_CQSQ; ik++){ HzzLambda_qsq[ic][ik]=100.; HwwLambda_qsq[ic][ik]=100.; } }
      SetJHUGenDistinguishWWCouplings(false);

      if (process == TVar::HSMHiggs){ Hzzcoupl[gHIGGS_VV_1][0] = 1.; Hwwcoupl[gHIGGS_VV_1][0] = 1.; }
      else if (process == TVar::H0_g1prime2){ Hzzcoupl[gHIGGS_VV_1_PRIME2][0] = 1.; Hwwcoupl[gHIGGS_VV_1_PRIME2][0] = 1.; }
      else if (process == TVar::H0hplus){ Hzzcoupl[gHIGGS_VV_2][0] = 1.; Hwwcoupl[gHIGGS_VV_2][0] = 1.; }
      else if (process == TVar::H0minus){ Hzzcoupl[gHIGGS_VV_4][0] = 1.; Hwwcoupl[gHIGGS_VV_4][0] = 1.; }
      else if (process == TVar::H0_Zgsg1prime2){ Hzzcoupl[gHIGGS_ZA_1_PRIME2][0] = 1.; }
      else if (process == TVar::H0_Zgs){ Hzzcoupl[gHIGGS_ZA_2][0] = 1.; }
      else if (process == TVar::H0_Zgs_PS){ Hzzcoupl[gHIGGS_ZA_4][0] = 1.; }
      else if (process == TVar::H0_gsgs){ Hzzcoupl[gHIGGS_AA_2][0] = 1.; }
      else if (process == TVar::H0_gsgs_PS){ Hzzcoupl[gHIGGS_AA_4][0] = 1.; }
      else if (process == TVar::SelfDefine_spin0){
        for (int j=0; j<2; j++){
          for (int i=0; i<SIZE_HVV; i++){
            Hzzcoupl[i][j] = (selfDSpinZeroCoupl.Hzzcoupl)[i][j]; Hwwcoupl[i][j] = (selfDSpinZeroCoupl.Hwwcoupl)[i][j];
            Hzzpcoupl[i][j] = (selfDSpinZeroCoupl.Hzzpcoupl)[i][j]; Hwwpcoupl[i][j] = (selfDSpinZeroCoupl.Hwwpcoupl)[i][j];
            Hzpzpcoupl[i][j] = (selfDSpinZeroCoupl.Hzpzpcoupl)[i][j]; Hwpwpcoupl[i][j] = (selfDSpinZeroCoupl.Hwpwpcoupl)[i][j];
          }
          for (int i=0; i<SIZE_Vpff; i++){
            Zpffcoupl[i][j] = (selfDVprimeCoupl.Zpffcoupl)[i][j]; Wpffcoupl[i][j] = (selfDVprimeCoupl.Wpffcoupl)[i][j];
          }
        }
        M_Zprime = selfDVprimeCoupl.M_Zprime;
        Ga_Zprime = selfDVprimeCoupl.Ga_Zprime;
        M_Wprime = selfDVprimeCoupl.M_Wprime;
        Ga_Wprime = selfDVprimeCoupl.Ga_Wprime;
        for (int j=0; j<SIZE_HVV_CQSQ; j++){
          for (int i=0; i<SIZE_HVV_LAMBDAQSQ; i++){ HzzLambda_qsq[i][j] = (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j]; HwwLambda_qsq[i][j] = (selfDSpinZeroCoupl.HwwLambda_qsq)[i][j]; }
          HzzCLambda_qsq[j] = (selfDSpinZeroCoupl.HzzCLambda_qsq)[j]; HwwCLambda_qsq[j] = (selfDSpinZeroCoupl.HwwCLambda_qsq)[j];
        }
        SetJHUGenDistinguishWWCouplings(selfDSpinZeroCoupl.separateWWZZcouplings);
      }
      SetJHUGenSpinZeroVVCouplings(Hzzcoupl, Hzzpcoupl, Hzpzpcoupl, HzzCLambda_qsq, HzzLambda_qsq, false);
      SetJHUGenSpinZeroVVCouplings(Hwwcoupl, Hwwpcoupl, Hwpwpcoupl, HwwCLambda_qsq, HwwLambda_qsq, true); // Set the WW couplings regardless of SetJHUGenDistinguishWWCouplings(false/true) because of how JHUGen handles this true flag.
      SetJHUGenVprimeContactCouplings(Zpffcoupl, Wpffcoupl);
      SetZprimeMassWidth(M_Zprime, Ga_Zprime);
      SetWprimeMassWidth(M_Wprime, Ga_Wprime);
    }

    dXsec = HJJMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, verbosity);
    if (verbosity >= TVar::DEBUG) MELAout << "TEvtProb::XsecCalc_XJJ: Process " << TVar::ProcessName(process) << " dXsec=" << dXsec << endl;
  }
  else MELAerr << "TEvtProb::XsecCalc_XJJ: Non-JHUGen vbfMELA MEs are not supported!" << endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  if (verbosity>=TVar::DEBUG) MELAout << "End XsecCalc_XJJ" << endl;
  return dXsec;
}

// Cross-section calculations for H + 1 jet (only SM for the moment)
double TEvtProb::XsecCalcXJ(){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin XsecCalc_XJ" << endl;
  double dXsec = 0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  // Calculate the ME
  if (matrixElement == TVar::JHUGen){
    // Not currently supported, but implement anyway
    double Hggcoupl[SIZE_HGG][2] ={ { 0 } };
    if (process == TVar::HSMHiggs) Hggcoupl[gHIGGS_GG_2][0] = 1.;
    else if (process == TVar::H0minus) Hggcoupl[gHIGGS_GG_4][0] = 1.;
    else if (process == TVar::SelfDefine_spin0){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_HGG; i++) Hggcoupl[i][j] = (selfDSpinZeroCoupl.Hggcoupl)[i][j];
      }
    }
    SetJHUGenSpinZeroGGCouplings(Hggcoupl);

    dXsec = HJJMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, verbosity);
    if (verbosity >= TVar::DEBUG) MELAout << "TEvtProb::XsecCalc_XJ: Process " << TVar::ProcessName(process) << " dXsec=" << dXsec << endl;
  }
  else MELAerr << "TEvtProb::XsecCalc_XJ: Non-JHUGen HJ MEs are not supported!" << endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  if (verbosity>=TVar::DEBUG) MELAout << "End XsecCalc_XJ" << endl;
  return dXsec;
}

// Cross-section calculations for VH
double TEvtProb::XsecCalc_VX(
  bool includeHiggsDecay
  ){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin XsecCalc_VX" << endl;
  if (matrixElement == TVar::MCFM) return XsecCalc_VVXVV();
  double dXsec = 0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  AllowSeparateWWCouplings(false); // HZZ couplings are used for both in spin-0

  if (matrixElement == TVar::JHUGen){
    double Hvvcoupl[SIZE_HVV][2] ={ { 0 } };
    double HvvLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ] ={ { 0 } };
    int HvvCLambda_qsq[SIZE_HVV_CQSQ] ={ 0 };
    double Hvvpcoupl[SIZE_HVV][2] = { { 0 } };
    double Hvpvpcoupl[SIZE_HVV][2] = { { 0 } };
    double Zpffcoupl[SIZE_Vpff][2] = { { 0 } };
    double Wpffcoupl[SIZE_Vpff][2] = { { 0 } };
    double M_Zprime = -1;
    double Ga_Zprime = 0;
    double M_Wprime = -1;
    double Ga_Wprime = 0;

    for (int ic=0; ic<SIZE_HVV_LAMBDAQSQ; ic++){ for (int ik=0; ik<SIZE_HVV_CQSQ; ik++) HvvLambda_qsq[ic][ik]=100.; }

    if (process == TVar::HSMHiggs) Hvvcoupl[gHIGGS_VV_1][0] = 1.;
    else if (process == TVar::H0_g1prime2) Hvvcoupl[gHIGGS_VV_1_PRIME2][0] = 1.;
    else if (process == TVar::H0hplus) Hvvcoupl[gHIGGS_VV_2][0] = 1.;
    else if (process == TVar::H0minus) Hvvcoupl[gHIGGS_VV_4][0] = 1.;
    else if (process == TVar::H0_Zgsg1prime2) Hvvcoupl[gHIGGS_ZA_1_PRIME2][0] = 1.;
    else if (process == TVar::H0_Zgs) Hvvcoupl[gHIGGS_ZA_2][0] = 1.;
    else if (process == TVar::H0_Zgs_PS) Hvvcoupl[gHIGGS_ZA_4][0] = 1.;
    else if (process == TVar::H0_gsgs) Hvvcoupl[gHIGGS_AA_2][0] = 1.;
    else if (process == TVar::H0_gsgs_PS) Hvvcoupl[gHIGGS_AA_4][0] = 1.;
    else if (process == TVar::SelfDefine_spin0){
      for (int i=0; i<SIZE_HVV; i++){
        for (int j=0; j<2; j++){
          Hvvcoupl[i][j] = (selfDSpinZeroCoupl.Hzzcoupl)[i][j];

          Hvvpcoupl[i][j] = (selfDSpinZeroCoupl.Hzzpcoupl)[i][j];
          Hvpvpcoupl[i][j] = (selfDSpinZeroCoupl.Hzpzpcoupl)[i][j];
        }
      }
      for (int i=0; i<SIZE_Vpff; i++){
        for (int j=0; j<2; j++){
          Zpffcoupl[i][j] = (selfDVprimeCoupl.Zpffcoupl)[i][j];
          Wpffcoupl[i][j] = (selfDVprimeCoupl.Wpffcoupl)[i][j];
        }
      }
      M_Zprime = selfDVprimeCoupl.M_Zprime;
      Ga_Zprime = selfDVprimeCoupl.Ga_Zprime;
      M_Wprime = selfDVprimeCoupl.M_Wprime;
      Ga_Wprime = selfDVprimeCoupl.Ga_Wprime;
      for (int j=0; j<SIZE_HVV_CQSQ; j++){
        for (int i=0; i<SIZE_HVV_LAMBDAQSQ; i++) HvvLambda_qsq[i][j] = (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j];
        HvvCLambda_qsq[j] = (selfDSpinZeroCoupl.HzzCLambda_qsq)[j];
      }
    }
    SetJHUGenSpinZeroVVCouplings(Hvvcoupl, Hvvpcoupl, Hvpvpcoupl, HvvCLambda_qsq, HvvLambda_qsq, false);
    SetJHUGenVprimeContactCouplings(Zpffcoupl, Wpffcoupl);
    SetZprimeMassWidth(M_Zprime, Ga_Zprime);
    SetWprimeMassWidth(M_Wprime, Ga_Wprime);

    dXsec = VHiggsMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, includeHiggsDecay, verbosity);
    if (verbosity >= TVar::DEBUG) MELAout << "TEVtProb::XsecCalc_VX: Process " << TVar::ProcessName(process) << " dXsec=" << dXsec << endl;
  } // end of JHUGen matrix element calculations
  else MELAerr << "TEVtProb::XsecCalc_VX: Non-JHUGen VH MEs are not supported!" << endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  if (verbosity>=TVar::DEBUG) MELAout << "End XsecCalc_VX" << endl;
  return dXsec;
}

// Cross-section calculations for ttH or bbH
double TEvtProb::XsecCalc_TTX(
  int topProcess, int topDecay
  ){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin XsecCalc_TTX" << endl;
  double dXsec = 0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  if (matrixElement == TVar::JHUGen){
    double Hqqcoupl[SIZE_HQQ][2]={ { 0 } };
    if (process == TVar::HSMHiggs) Hqqcoupl[gHIGGS_KAPPA][0] = 1.;
    else if (process == TVar::H0minus) Hqqcoupl[gHIGGS_KAPPA_TILDE][0] = 1.;
    else if (process == TVar::SelfDefine_spin0){
      for (int i=0; i<SIZE_HQQ; i++){
        for (int j=0; j<2; j++){
          if ((selfDSpinZeroCoupl.Httcoupl)[i][j]!=0. && production==TVar::ttH) Hqqcoupl[i][j] = (selfDSpinZeroCoupl.Httcoupl)[i][j];
          else if ((selfDSpinZeroCoupl.Hbbcoupl)[i][j]!=0. && production==TVar::bbH) Hqqcoupl[i][j] = (selfDSpinZeroCoupl.Hbbcoupl)[i][j];
          else Hqqcoupl[i][j] = (selfDSpinZeroCoupl.Hqqcoupl)[i][j];
        }
      }
    }
    SetJHUGenSpinZeroQQCouplings(Hqqcoupl);

    if (production==TVar::ttH) dXsec = TTHiggsMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, topDecay, topProcess, verbosity);
    else if (production==TVar::bbH) dXsec = BBHiggsMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, topProcess, verbosity);
    else if (verbosity >= TVar::ERROR) MELAout << "TEvtProb::XsecCalc_TTX only supports ttH and bbH productions for the moment." << endl;
    if (verbosity >= TVar::DEBUG) MELAout << "TEvtProb::XsecCalc_TTX: Process " << TVar::ProcessName(process) << " dXsec=" << dXsec << endl;
  }
  else MELAerr << "TEvtProb::XsecCalc_TTX: Non-JHUGen ttH MEs are not supported!" << endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  if (verbosity>=TVar::DEBUG) MELAout << "End XsecCalc_TTX" << endl;
  return dXsec;
}

double TEvtProb::GetXPropagator(
  TVar::ResonancePropagatorScheme scheme
  ){
  if (!CheckInputPresent()) return 0.;
  return TUtil::ResonancePropagator(melaCand->m(), scheme);
}


bool TEvtProb::CheckSelfDCouplings_Hgg(){
  if (_hmass>=0. && _hwidth>0.){
    for (int vv = 0; vv < SIZE_HGG; vv++){
      if (
        (selfDSpinZeroCoupl.Hggcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hggcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.Hg4g4coupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hg4g4coupl)[vv][0] != 0
        ){
        return true;
      }
    }
  }
  if (_h2mass>=0. && _h2width>0.){
    for (int vv = 0; vv < SIZE_HGG; vv++){
      if (
        (selfDSpinZeroCoupl.H2ggcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2ggcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.H2g4g4coupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2g4g4coupl)[vv][0] != 0
        ){
        return true;
      }
    }
  }
  return false;
}
bool TEvtProb::CheckSelfDCouplings_Hqq(){
  if (_hmass>=0. && _hwidth>0.){
    for (int vv = 0; vv < SIZE_HQQ; vv++){
      if (
        (selfDSpinZeroCoupl.Hqqcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hqqcoupl)[vv][0] != 0
        ){
        return true;
      }
    }
  }
  if (_h2mass>=0. && _h2width>0.){
    for (int vv = 0; vv < SIZE_HQQ; vv++){
      if (
        (selfDSpinZeroCoupl.H2qqcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2qqcoupl)[vv][0] != 0
        ){
        return true;
      }
    }
  }
  return false;
}
bool TEvtProb::CheckSelfDCouplings_Htt(){
  if (_hmass>=0. && _hwidth>0.){
    for (int vv = 0; vv < SIZE_HQQ; vv++){
      if (
        (selfDSpinZeroCoupl.Httcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Httcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.Ht4t4coupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Ht4t4coupl)[vv][0] != 0
        ){
        return true;
      }
    }
  }
  if (_h2mass>=0. && _h2width>0.){
    for (int vv = 0; vv < SIZE_HQQ; vv++){
      if (
        (selfDSpinZeroCoupl.H2ttcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2ttcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.H2t4t4coupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2t4t4coupl)[vv][0] != 0
        ){
        return true;
      }
    }
  }
  return false;
}
bool TEvtProb::CheckSelfDCouplings_Hbb(){
  if (_hmass>=0. && _hwidth>0.){
    for (int vv = 0; vv < SIZE_HQQ; vv++){
      if (
        (selfDSpinZeroCoupl.Hbbcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hbbcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.Hb4b4coupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hb4b4coupl)[vv][0] != 0
        ){
        return true;
      }
    }
  }
  if (_h2mass>=0. && _h2width>0.){
    for (int vv = 0; vv < SIZE_HQQ; vv++){
      if (
        (selfDSpinZeroCoupl.H2bbcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2bbcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.H2b4b4coupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2b4b4coupl)[vv][0] != 0
        ){
        return true;
      }
    }
  }
  return false;
}
bool TEvtProb::CheckSelfDCouplings_HVV(){
  if (_hmass>=0. && _hwidth>0.){
    for (int vv = 0; vv < SIZE_HVV; vv++){
      if (
        (selfDSpinZeroCoupl.Hzzcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hzzcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.Hwwcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hwwcoupl)[vv][0] != 0
        ){
        return true;
      }
    } // No need to check c_q**2. If these are 0, z_q**2 do not have any effect.
  }
  if (_h2mass>=0. && _h2width>0.){
    for (int vv = 0; vv < SIZE_HVV; vv++){
      if (
        (selfDSpinZeroCoupl.H2zzcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2zzcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.H2wwcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2wwcoupl)[vv][0] != 0
        ){
        return true;
      }
    } // No need to check c_q**2. If these are 0, z_q**2 do not have any effect.
  }
  return false;
}


