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

#include <interface/TEvtProb.hh>


ClassImp(TEvtProb)


using namespace std;
using namespace TUtil;


//-----------------------------------------------------------------------------
// Constructors and Destructor
//-----------------------------------------------------------------------------
TEvtProb::TEvtProb(
  const char* path, double ebeam, const char* pathtoPDFSet, int PDFMember, TVar::VerbosityLevel verbosity_
  ) :
  verbosity(verbosity_),
  EBEAM(ebeam)
{
  if (verbosity>=TVar::DEBUG) cout << "Begin TEvtProb constructor" << endl;

  SetLeptonInterf(TVar::DefaultLeptonInterf);

  /***** Initialize Higgs width reader *****/
  string path_string = path;
  myCSW_ = new MELAHXSWidth(path_string);
  if (verbosity>=TVar::DEBUG) cout << "TEvtProb::TEvtProb: HXS successful" << endl;

  /***** Initialize MCFM *****/
  InitializeMCFM();

  /***** Initialize JHUGen *****/
  InitializeJHUGen(pathtoPDFSet, PDFMember);

  /***** Initialize schemes *****/
  ResetCouplings();
  ResetRenFacScaleMode();
  ResetInputEvent();
  SetPrimaryHiggsMass(125.); // Should come after InitializeMCFM and InitializeJHUGen
  SetCandidateDecayMode(TVar::CandidateDecay_ZZ);

  if (verbosity>=TVar::DEBUG) cout << "End TEvtProb constructor" << endl;
}

TEvtProb::~TEvtProb(){
  if (verbosity>=TVar::DEBUG) cout << "Begin TEvtProb destructor" << endl;

  ResetInputEvent();
  if (myCSW_!=0) delete myCSW_;

  if (verbosity>=TVar::DEBUG) cout << "End TEvtProb destructor" << endl;
}


void TEvtProb::InitializeMCFM(){
  if (verbosity>=TVar::DEBUG) cout << "Begin TEvtProb::InitializeMCFM" << endl;

  mcfm_init_((char *)"input.DAT", (char *)"./");
  if (verbosity>=TVar::DEBUG) cout << "TEvtProb::TEvtProb: mcfm_init successful" << endl;
  SetEwkCouplingParameters();
  energy_.sqrts = 2.*EBEAM;
  coupling_();
  qlinit_();
  // First resonance constant parameters
  spinzerohiggs_anomcoupl_.LambdaBSM=1000;
  spinzerohiggs_anomcoupl_.Lambda_z1=10000;
  spinzerohiggs_anomcoupl_.Lambda_z2=10000;
  spinzerohiggs_anomcoupl_.Lambda_z3=10000;
  spinzerohiggs_anomcoupl_.Lambda_z4=10000;
  spinzerohiggs_anomcoupl_.Lambda_zgs1=10000;
  spinzerohiggs_anomcoupl_.Lambda_Q=10000;
  // Second resonance constant parameters
  spinzerohiggs_anomcoupl_.Lambda2BSM=1000;
  spinzerohiggs_anomcoupl_.Lambda2_z1=10000;
  spinzerohiggs_anomcoupl_.Lambda2_z2=10000;
  spinzerohiggs_anomcoupl_.Lambda2_z3=10000;
  spinzerohiggs_anomcoupl_.Lambda2_z4=10000;
  spinzerohiggs_anomcoupl_.Lambda2_zgs1=10000;
  spinzerohiggs_anomcoupl_.Lambda2_Q=10000;
  // Constant parameters for all processes
  qlfirst_.qlfirst=false;

  if (verbosity>=TVar::DEBUG) cout << "End TEvtProb::InitializeMCFM" << endl;
}
void TEvtProb::InitializeJHUGen(const char* pathtoPDFSet, int PDFMember){
  if (verbosity>=TVar::DEBUG) cout << "Begin TEvtProb::InitializeJHUGen" << endl;

  InitJHUGenMELA(pathtoPDFSet, PDFMember);

  if (verbosity>=TVar::DEBUG) cout << "End TEvtProb::InitializeJHUGen" << endl;
}



// Set NNPDF driver path
void TEvtProb::Set_LHAgrid(const char* path, int pdfmember){
  char path_nnpdf_c[200];
  sprintf(path_nnpdf_c, "%s", path);
  int pathLength = strlen(path_nnpdf_c);
  nnpdfdriver_(path_nnpdf_c, &pathLength);
  nninitpdf_(&pdfmember);
}
void TEvtProb::SetProcess(TVar::Process tmp) { process = tmp; }
void TEvtProb::SetMatrixElement(TVar::MatrixElement tmp){ matrixElement = tmp; }
void TEvtProb::SetProduction(TVar::Production tmp){ production = tmp; }
void TEvtProb::SetVerbosity(TVar::VerbosityLevel tmp){ verbosity = tmp; }
void TEvtProb::SetLeptonInterf(TVar::LeptonInterference tmp){ leptonInterf = tmp; }
void TEvtProb::SetCandidateDecayMode(TVar::CandidateDecayMode mode){
  if (mode==TVar::CandidateDecay_WW) PDGHelpers::setHVVmass(PDGHelpers::Wmass);
  else if (mode==TVar::CandidateDecay_ff || mode==TVar::CandidateDecay_Stable || mode==TVar::CandidateDecay_GG) PDGHelpers::setHVVmass(PDGHelpers::Zeromass);
  else PDGHelpers::setHVVmass(PDGHelpers::Zmass); // Anything that contains a Z
}
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
      _hwidth = myCSW_->HiggsWidth(_hmass);
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
      _h2width = myCSW_->HiggsWidth(_h2mass);
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

  if (verbosity>=TVar::DEBUG) cout
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
    cerr << "TEvtProb::AppendTopCandidate: No MELACandidates are present to append this top!" << endl;
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
  else cerr << "TEvtProb::SetCurrentCandidateFromIndex: icand=" << icand << ">=candList.size()=" << candList.size() << endl;
}
void TEvtProb::SetCurrentCandidate(MELACandidate* cand){
  melaCand = cand;
  if (verbosity>=TVar::INFO && melaCand==0) cout << "TEvtProb::SetCurrentCandidate: BE CAREFUL! melaCand==0!" << endl;
  if (verbosity>=TVar::INFO && GetCurrentCandidateIndex()<0) cout << "TEvtProb::SetCurrentCandidate: The current candidate is not in the list of candidates. It is the users' responsibility to delete this candidate and all of its associated particles." << endl;
}

// Reset functions
void TEvtProb::ResetIORecord(){ RcdME.reset(); }
void TEvtProb::ResetRenFacScaleMode(){ SetRenFacScaleMode(TVar::DefaultScaleScheme, TVar::DefaultScaleScheme, 0.5, 0.5); }
void TEvtProb::ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme){
  if (ext_ewscheme<-1 || ext_ewscheme>3) ext_ewscheme=3;
  ewinput_.Gf_inp = ext_Gf;
  ewinput_.aemmz_inp = ext_aemmz;
  ewinput_.wmass_inp = ext_mW;
  ewinput_.zmass_inp = ext_mZ;
  ewinput_.xw_inp = ext_xW;
  ewscheme_.ewscheme = ext_ewscheme;
  coupling_();
}
void TEvtProb::ResetCouplings(){
  selfDSpinZeroCoupl.reset();
  selfDSpinOneCoupl.reset();
  selfDSpinTwoCoupl.reset();
  AllowSeparateWWCouplings(false);
}
void TEvtProb::ResetInputEvent(){
  RcdME.melaCand = 0;
  melaCand = 0;

  // Clear bookkeeping objects
  // Clean MELACandidates first since they contain all other objects
  for (unsigned int p=0; p<candList.size(); p++){
    MELACandidate* tmpCand = (MELACandidate*)candList.at(p);
    if (tmpCand!=0) delete tmpCand;
  }
  candList.clear();
  // Clean MELATopCandidates next since they contain MELAParticles
  for (unsigned int p=0; p<topCandList.size(); p++){
    MELATopCandidate* tmpCand = (MELATopCandidate*)topCandList.at(p);
    if (tmpCand!=0) delete tmpCand;
  }
  topCandList.clear();
  // Clean all remaining MELAPArticles
  for (unsigned int p=0; p<particleList.size(); p++){
    MELAParticle* tmpPart = (MELAParticle*)particleList.at(p);
    if (tmpPart!=0) delete tmpPart;
  }
  particleList.clear();
}

// Get-functions
SpinZeroCouplings* TEvtProb::GetSelfDSpinZeroCouplings(){ return selfDSpinZeroCoupl.getRef(); }
SpinOneCouplings* TEvtProb::GetSelfDSpinOneCouplings(){ return selfDSpinOneCoupl.getRef(); }
SpinTwoCouplings* TEvtProb::GetSelfDSpinTwoCouplings(){ return selfDSpinTwoCoupl.getRef(); }
double TEvtProb::GetPrimaryHiggsMass(){ return PrimaryHMass; }
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
  if (melaCand==0 || candList.size()==0){
    cerr
      << "TEvtProb::XsecCalc_XVV: melaCand==" << melaCand << " is nullPtr"
      << " or candList.size()==" << candList.size() << " is problematic!"
      << endl;
    if (candList.size()==0) return false;
    else{
      SetCurrentCandidateFromIndex(candList.size()-1);
      cerr << "TEvtProb::XsecCalc_XVV: melaCand now points to the latest candidate (cand" << (candList.size()-1) << ")" << endl;
    }
  }
  SetRcdCandPtr(); // If input event is present, set the RcdME pointer to melaCand
  return true;
}


// ME functions

//
// Directly calculate the VV->4f differential cross-section
//
double TEvtProb::XsecCalc_XVV(
  TVar::Process process_, TVar::Production production_,
  TVar::VerbosityLevel verbosity_
  ){
  if (verbosity_>=TVar::DEBUG) cout << "Begin XsecCalc_XVV"<< endl;
  double dXsec=0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  //Initialize Process
  SetProcess(process_);
  SetProduction(production_);
  SetVerbosity(verbosity_);

  bool forceUseMCFM = (matrixElement == TVar::MCFM || process == TVar::bkgZZ_SMHiggs);
  bool calculateME=true;

  bool needBSMHiggs=false;
  if (forceUseMCFM){
    if (verbosity_>=TVar::DEBUG) cout << "TEvtProb::XsecCalc_XVV: Try MCFM"<< endl;
    // Check self-defined couplings are specified.
    for (int vv = 0; vv < SIZE_HVV; vv++){
      if (
        (selfDSpinZeroCoupl.Hzzcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hzzcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.Hwwcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hwwcoupl)[vv][0] != 0
        ){
        needBSMHiggs = true; break;
      } // Only possible if selfDefine is called.
    }
    if (_h2mass>=0. && _h2width>0. && !needBSMHiggs){
      for (int vv = 0; vv < SIZE_HVV; vv++){
        if (
          (selfDSpinZeroCoupl.H2zzcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2zzcoupl)[vv][0] != 0
          ||
          (selfDSpinZeroCoupl.H2wwcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2wwcoupl)[vv][0] != 0
          ){
          needBSMHiggs = true; break;
        } // Only possible if selfDefine is called.
      }
    }
    if (needBSMHiggs) SetLeptonInterf(TVar::InterfOn); // All anomalous coupling computations have to use lepton interference
    SetMCFMSpinZeroVVCouplings(needBSMHiggs, &selfDSpinZeroCoupl, false);
    calculateME = MCFM_chooser(process, production, leptonInterf, melaCand);
    if (verbosity_>=TVar::INFO && !calculateME) cout << "TEvtProb::XsecCalc_XVV: MCFM_chooser failed to determine the process configuration."<< endl;
  }

  // Final check before moving on to ME calculations
  if (!calculateME){
    if (verbosity_>=TVar::DEBUG) cout << "XsecCalc_XVV::calculateME is false"<< endl;
    return dXsec;
  }
  // ME calculations
  if (forceUseMCFM) dXsec = SumMatrixElementPDF(process, production, matrixElement, &event_scales, &RcdME, EBEAM, (selfDSpinZeroCoupl.Hvvcoupl_freenorm), verbosity);
  else if (matrixElement == TVar::JHUGen){
    if (verbosity_>=TVar::DEBUG) cout << "XsecCalc_XVV::Try JHUGen"<< endl;
    AllowSeparateWWCouplings(false); // HZZ couplings are used for both in spin-0
    // all the possible couplings
    double Hggcoupl[SIZE_HGG][2] ={ { 0 } };
    double Hvvcoupl[SIZE_HVV][2] ={ { 0 } };
    double HvvLambda_qsq[4][3] ={ { 0 } };
    int HvvCLambda_qsq[3] ={ 0 };

    double Zqqcoupl[SIZE_ZQQ][2] ={ { 0 } };
    double Zvvcoupl[SIZE_ZVV][2] ={ { 0 } };

    double Gqqcoupl[SIZE_GQQ][2] ={ { 0 } };
    double Gggcoupl[SIZE_GGG][2] ={ { 0 } };
    double Gvvcoupl[SIZE_GVV][2] ={ { 0 } };

    //
    // set spin 0 default numbers
    //
    // By default set the Spin 0 couplings for SM case (0+m)
    Hggcoupl[0][0]=1.0;  Hggcoupl[0][1]=0.0;
    for (int ic=0; ic<4; ic++){ for (int ik=0; ik<3; ik++) HvvLambda_qsq[ic][ik]=100.; }
    //
    // set spin 2 default numbers
    //
    Gqqcoupl[0][0]=1.0;  Gqqcoupl[0][1]=0.0;
    Gqqcoupl[1][0]=1.0;  Gqqcoupl[1][1]=0.0;
    //
    // set spin 1 default numbers
    //
    Zqqcoupl[0][0]=1.0;  Zqqcoupl[0][1]=0.0;
    Zqqcoupl[1][0]=1.0;  Zqqcoupl[1][1]=0.0;

    bool isSpinZero = false;
    bool isSpinOne = false;
    bool isSpinTwo = false;

    // 0m+
    if (process == TVar::HSMHiggs){
      Hvvcoupl[0][0]=1.;
      isSpinZero = true;
    }
    // 0-
    else if (process == TVar::H0minus){
      Hvvcoupl[3][0] = 1.;
      isSpinZero = true;
    }
    // 0h+
    else if (process == TVar::H0hplus){
      Hvvcoupl[1][0] = 1.;
      isSpinZero = true;
    }
    // 0+L1
    else if (process == TVar::H0_g1prime2){
      Hvvcoupl[11][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_Zgs){
      Hvvcoupl[4][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_gsgs){
      Hvvcoupl[7][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_Zgs_PS){
      Hvvcoupl[6][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_gsgs_PS){
      Hvvcoupl[9][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::H0_Zgsg1prime2){
      Hvvcoupl[30][0] = 1.;
      isSpinZero = true;
    }
    else if (process == TVar::SelfDefine_spin0){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_HGG; i++) Hggcoupl[i][j] = (selfDSpinZeroCoupl.Hggcoupl)[i][j];
        for (int i=0; i<SIZE_HVV; i++) Hvvcoupl[i][j] = (selfDSpinZeroCoupl.Hzzcoupl)[i][j];
      }
      for (int j=0; j<3; j++){
        for (int i=0; i<4; i++) HvvLambda_qsq[i][j] = (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j];
        HvvCLambda_qsq[j] = (selfDSpinZeroCoupl.HzzCLambda_qsq)[j];
      }
      if (verbosity_>=TVar::DEBUG_VERBOSE){
        for (int j=0; j<2; j++){
          for (int i=0; i<SIZE_HGG; i++) cout << "Hggcoupl[" << i << "][" << j << "] = " << (selfDSpinZeroCoupl.Hggcoupl)[i][j] << endl;
          for (int i=0; i<SIZE_HVV; i++) cout << "Hvvcoupl[" << i << "][" << j << "] = " << (selfDSpinZeroCoupl.Hzzcoupl)[i][j] << endl;
        }
        for (int j=0; j<3; j++){
          for (int i=0; i<4; i++) cout << "HvvLambda_qsq[" << i << "][" << j << "] = " << (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j] << endl;
          cout << "HvvCLambda_qsq[" << j << "] = " << (selfDSpinZeroCoupl.HzzCLambda_qsq)[j] << endl;
        }
      }
      isSpinZero = true;
    }

    // 2m+
    else if (process == TVar::H2_g1g5){
      Gggcoupl[0][0]=1.;
      Gvvcoupl[0][0]=1.;
      Gvvcoupl[4][0]=1.;
      isSpinTwo = true;
    }
    // 2b+
    else if (process == TVar::H2_g1){
      Gggcoupl[0][0]=1.;
      Gvvcoupl[0][0]=1.;
      isSpinTwo = true;
    }
    else if (process == TVar::H2_g2){
      Gggcoupl[1][0]=1.;
      Gvvcoupl[1][0]=1.;
      isSpinTwo = true;
    }
    // 2h3plus
    else if (process == TVar::H2_g3){
      Gggcoupl[2][0]=1.;
      Gvvcoupl[2][0]=1.;
      isSpinTwo = true;
    }
    // 2h+
    else if (process == TVar::H2_g4){
      Gggcoupl[3][0]=1.;
      Gvvcoupl[3][0]=1.;
      isSpinTwo = true;
    }
    // 2b+
    else if (process == TVar::H2_g5){
      Gggcoupl[0][0]=1.;
      Gvvcoupl[4][0]=1.;
      isSpinTwo = true;
    }
    // 2h6+
    else if (process == TVar::H2_g6){
      Gggcoupl[0][0]=1.;
      Gvvcoupl[5][0]=1.;
      isSpinTwo = true;
    }
    // 2h7plus
    else if (process == TVar::H2_g7){
      Gggcoupl[0][0]=1.;
      Gvvcoupl[6][0]=1.;
      isSpinTwo = true;
    }
    // 2h-
    else if (process == TVar::H2_g8){
      Gggcoupl[4][0]=1.;
      Gvvcoupl[7][0]=1.;
      isSpinTwo = true;
    }
    // 2h9minus
    else if (process == TVar::H2_g9){
      Gggcoupl[4][0]=1.;
      Gvvcoupl[8][0]=1.;
      isSpinTwo = true;
    }
    // 2h10minus
    else if (process == TVar::H2_g10){
      Gggcoupl[4][0]=1.;
      Gvvcoupl[9][0]=1.;
      isSpinTwo = true;
    }
    else if (process == TVar::SelfDefine_spin2){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_GGG; i++) Gggcoupl[i][j] = (selfDSpinTwoCoupl.Gggcoupl)[i][j];
        for (int i=0; i<SIZE_GQQ; i++) Gqqcoupl[i][j] = (selfDSpinTwoCoupl.Gqqcoupl)[i][j];
        for (int i=0; i<SIZE_GVV; i++) Gvvcoupl[i][j] = (selfDSpinTwoCoupl.Gvvcoupl)[i][j];
      }
      isSpinTwo = true;
    }

    // 1-
    else if (process == TVar::H1minus){
      Zvvcoupl[0][0]=1.;
      isSpinOne = true;
    }
    // 1+
    else if (process == TVar::H1plus){
      Zvvcoupl[1][0]=1.;
      isSpinOne = true;
    }
    else if (process == TVar::SelfDefine_spin1){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_ZQQ; i++) Zqqcoupl[i][j] = (selfDSpinOneCoupl.Zqqcoupl)[i][j];
        for (int i=0; i<SIZE_ZVV; i++) Zvvcoupl[i][j] = (selfDSpinOneCoupl.Zvvcoupl)[i][j];
      }
      isSpinOne = true;
    }

    if (isSpinZero){
      SetJHUGenSpinZeroGGCouplings(Hggcoupl);
      SetJHUGenSpinZeroVVCouplings(Hvvcoupl, HvvCLambda_qsq, HvvLambda_qsq, false);
    }
    else if (isSpinOne) SetJHUGenSpinOneCouplings(Zqqcoupl, Zvvcoupl);
    else if (isSpinTwo) SetJHUGenSpinTwoCouplings(Gggcoupl, Gvvcoupl, Gqqcoupl);

    if (isSpinZero || isSpinOne || isSpinTwo) dXsec = JHUGenMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, verbosity);
    else cerr
      << "TEvtProb::XsecCalc_XVV: JHUGen ME is not spin zero, one or two! The process is described by "
      << "Process: " << process << ", Production: " << production << ", and ME: " << matrixElement
      << endl;
  } // end of JHUGen matrix element calculations

  if (verbosity >= TVar::DEBUG) cout
    << "Process " << TVar::ProcessName(process)
    << " TEvtProb::XsecCalc(): dXsec=" << dXsec
    << endl;

  if (verbosity_>=TVar::DEBUG) cout << "XsecCalc_XVV::Reset couplings"<< endl;
  ResetCouplings(); // Should come first
  if (forceUseMCFM){ // Set defaults. Should come next...
    if (needBSMHiggs) SetLeptonInterf(TVar::DefaultLeptonInterf);
    SetMCFMSpinZeroVVCouplings(false, &selfDSpinZeroCoupl, true); // ... because of this!
  }
  ResetRenFacScaleMode();
  if (verbosity_>=TVar::DEBUG) cout << "End XsecCalc_XVV"<< endl;
  return dXsec;
}

double TEvtProb::XsecCalc_VVXVV(
  TVar::Process process_, TVar::Production production_,
  TVar::VerbosityLevel verbosity_
  ){
  if (verbosity_>=TVar::DEBUG) cout << "Begin XsecCalc_VVXVV"<< endl;

  double dXsec=0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  //Initialize Process
  SetProcess(process_);
  SetProduction(production_);
  SetVerbosity(verbosity_);

  bool forceUseMCFM = (matrixElement == TVar::MCFM || process == TVar::bkgZZ_SMHiggs);
  bool needBSMHiggs=false;
  bool calculateME=true;
  // process == TVar::bkgZZ_SMHiggs && matrixElement == TVar::JHUGen is still MCFM
  if (forceUseMCFM){ // Always uses MCFM
    // Check self-defined couplings are specified.
    if (verbosity_>=TVar::DEBUG) cout << "TEvtProb::XsecCalc_VVXVV: Try MCFM"<< endl;
    for (int vv = 0; vv < SIZE_HVV; vv++){
      if (
        (selfDSpinZeroCoupl.Hzzcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hzzcoupl)[vv][0] != 0
        ||
        (selfDSpinZeroCoupl.Hwwcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hwwcoupl)[vv][0] != 0
        ){
        needBSMHiggs = true; break;
      } // Only possible if selfDefine is called.
    }
    if (_h2mass>=0. && _h2width>0. && !needBSMHiggs){
      for (int vv = 0; vv < SIZE_HVV; vv++){
        if (
          (selfDSpinZeroCoupl.H2zzcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2zzcoupl)[vv][0] != 0
          ||
          (selfDSpinZeroCoupl.H2wwcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.H2wwcoupl)[vv][0] != 0
          ){
          needBSMHiggs = true; break;
        } // Only possible if selfDefine is called.
      }
    }
    if (needBSMHiggs) SetLeptonInterf(TVar::InterfOn); // All anomalous coupling computations have to use lepton interference
    SetMCFMSpinZeroVVCouplings(needBSMHiggs, &selfDSpinZeroCoupl, false);
    calculateME = MCFM_chooser(process, production, leptonInterf, melaCand);
    if (verbosity_>=TVar::INFO && !calculateME) cout << "TEvtProb::XsecCalc_VVXVV: MCFM_chooser failed to determine the process configuration."<< endl;
  }

  // Last check before ME calculations
  if (!calculateME) return dXsec;

  // ME calculation
  if (forceUseMCFM) dXsec = SumMatrixElementPDF(process, production, matrixElement, &event_scales, &RcdME, EBEAM, (selfDSpinZeroCoupl.Hvvcoupl_freenorm), verbosity);
  else cerr << "Non-MCFM Mes are not supported." << endl;

  if (verbosity >= TVar::DEBUG) cout
    << "Process " << TVar::ProcessName(process)
    << " TEvtProb::XsecCalc_VVXVV(): dXsec=" << dXsec
    << endl;

  ResetCouplings(); // Should come first
  if (forceUseMCFM){ // Set defaults. Should come next...
    if (needBSMHiggs) SetLeptonInterf(TVar::DefaultLeptonInterf);
    SetMCFMSpinZeroVVCouplings(false, &selfDSpinZeroCoupl, true); // ... because of this!
  }
  ResetRenFacScaleMode();
  if (verbosity_>=TVar::DEBUG) cout << "End XsecCalc_VVXVV"<< endl;
  return dXsec;
}

// Cross-section calculations for H + 2 jets
double TEvtProb::XsecCalcXJJ(
  TVar::Process process_, TVar::Production production_,
  TVar::VerbosityLevel verbosity_
  ){
  if (verbosity_>=TVar::DEBUG) cout << "Begin XsecCalc_XJJ"<< endl;

  double dXsec = 0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  // Initialize Process
  SetProcess(process_);
  SetProduction(production_);
  SetVerbosity(verbosity_);

  // first/second number is the real/imaginary part
  double Hggcoupl[SIZE_HGG][2] ={ { 0 } };
  double Hzzcoupl[SIZE_HVV][2] ={ { 0 } };
  double Hwwcoupl[SIZE_HVV][2] ={ { 0 } };
  double HzzLambda_qsq[4][3] ={ { 0 } };
  double HwwLambda_qsq[4][3] ={ { 0 } };
  int HzzCLambda_qsq[3] ={ 0 };
  int HwwCLambda_qsq[3] ={ 0 };

  Hggcoupl[0][0]=1.0;  Hggcoupl[0][1]=0.0; // g2
  Hzzcoupl[0][0]=1.0;  Hzzcoupl[0][1]=0.0; // g1
  Hwwcoupl[0][0]=1.0;  Hwwcoupl[0][1]=0.0; // g1
  for (int ic=0; ic<4; ic++){
    for (int ik=0; ik<3; ik++){
      HzzLambda_qsq[ic][ik]=100.;
      HwwLambda_qsq[ic][ik]=100.;
    }
  }

  // 0-
  if (process == TVar::H0minus){
    Hggcoupl[0][0] = 0.;
    Hggcoupl[2][0] = 1.;

    Hzzcoupl[0][0] = 0.;
    Hzzcoupl[3][0] = 1.;
    Hwwcoupl[0][0] = 0.;
    Hwwcoupl[3][0] = 1.;
  }
  // 0+h
  else if (process == TVar::H0hplus) { // No need to re-set ggcoupl
    Hzzcoupl[0][0] = 0.;
    Hzzcoupl[1][0] = 1.;
    Hwwcoupl[0][0] = 0.;
    Hwwcoupl[1][0] = 1.;
  }
  // 0+L1
  else if (process == TVar::H0_g1prime2){ // No need to re-set ggcoupl
    Hzzcoupl[0][0] = 0.;
    Hzzcoupl[11][0] = 1.;
    Hwwcoupl[0][0] = 0.;
    Hwwcoupl[11][0] = 1.;
  }
  else if (process == TVar::SelfDefine_spin0){
    for (int j=0; j<2; j++){
      for (int i=0; i<SIZE_HGG; i++) Hggcoupl[i][j] = (selfDSpinZeroCoupl.Hggcoupl)[i][j];
      for (int i=0; i<SIZE_HVV; i++){
        Hzzcoupl[i][j] = (selfDSpinZeroCoupl.Hzzcoupl)[i][j];
        Hwwcoupl[i][j] = (selfDSpinZeroCoupl.Hwwcoupl)[i][j];
      }
    }
    for (int j=0; j<3; j++){
      for (int i=0; i<4; i++){
        HzzLambda_qsq[i][j] = (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j];
        HwwLambda_qsq[i][j] = (selfDSpinZeroCoupl.HwwLambda_qsq)[i][j];
      }
      HzzCLambda_qsq[j] = (selfDSpinZeroCoupl.HzzCLambda_qsq)[j];
      HwwCLambda_qsq[j] = (selfDSpinZeroCoupl.HwwCLambda_qsq)[j];
    }
  }
  SetJHUGenSpinZeroGGCouplings(Hggcoupl);
  SetJHUGenSpinZeroVVCouplings(Hzzcoupl, HzzCLambda_qsq, HzzLambda_qsq, false);
  SetJHUGenSpinZeroVVCouplings(Hwwcoupl, HwwCLambda_qsq, HwwLambda_qsq, true);

  // ME calculations
  if (matrixElement == TVar::JHUGen){
    dXsec = HJJMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, verbosity);
    if (verbosity >= TVar::DEBUG) cout <<"Process " << TVar::ProcessName(process) << " TEvtProb::XsecCalc_XJJ(): dXsec=" << dXsec << endl;
  }
  else cerr << "Non-JHUGen vbfMELA MEs are not supported!" << endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  if (verbosity_>=TVar::DEBUG) cout << "End XsecCalc_XJJ"<< endl;
  return dXsec;
}

// Cross-section calculations for H (SM) + 1 jet
double TEvtProb::XsecCalcXJ(TVar::Process process_, TVar::Production production_, TVar::VerbosityLevel verbosity_){
  double dXsec = 0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  // Initialize Process
  SetProcess(process_);
  SetProduction(production_);
  SetVerbosity(verbosity_);

  // Calculate the ME
  if (matrixElement == TVar::JHUGen){
    dXsec = HJJMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, verbosity);
    if (verbosity >= TVar::DEBUG) std::cout << "Process " << TVar::ProcessName(process) << " TEvtProb::XsecCalc_XJ(): dXsec=" << dXsec << endl;
  }
  else cerr << "Non-JHUGen HJ MEs are not supported!" << endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  return dXsec;
}


double TEvtProb::XsecCalc_VX(
  TVar::Process process_, TVar::Production production_,
  TVar::VerbosityLevel verbosity_,
  bool includeHiggsDecay
  ){
  if (verbosity_>=TVar::DEBUG) cout << "Begin XsecCalc_VX"<< endl;

  double dXsec = 0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  //Initialize Process
  SetProcess(process_);
  SetProduction(production_);
  SetVerbosity(verbosity_);
  AllowSeparateWWCouplings(false); // HZZ couplings are used for both in spin-0

  if (matrixElement == TVar::JHUGen){
    // Set Couplings at the HVV* vertex
    double Hvvcoupl[SIZE_HVV][2] ={ { 0 } };
    double HvvLambda_qsq[4][3] ={ { 0 } };
    int HvvCLambda_qsq[3] ={ 0 };

    // By default set the Spin 0 couplings for SM case
    Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
    for (int ic=0; ic<4; ic++){ for (int ik=0; ik<3; ik++) HvvLambda_qsq[ic][ik]=100.; }

    // 0-
    if (process == TVar::H0minus) {
      Hvvcoupl[0][0] = 0.;
      Hvvcoupl[1][0] = 0.;
      Hvvcoupl[2][0] = 0.;
      Hvvcoupl[3][0] = 1.;
    }
    // 0h+
    else if (process == TVar::H0hplus) {
      Hvvcoupl[0][0] = 0.;
      Hvvcoupl[1][0] = 1.;
      Hvvcoupl[2][0] = 0.;
      Hvvcoupl[3][0] = 0.;
    }
    // 0+L1
    else if (process == TVar::H0_g1prime2){
      Hvvcoupl[0][0] = 0.;
      Hvvcoupl[11][0] = 1.;
    }
    else if (process == TVar::SelfDefine_spin0){
      for (int i=0; i<SIZE_HVV; i++){ for (int j=0; j<2; j++) Hvvcoupl[i][j] = (selfDSpinZeroCoupl.Hzzcoupl)[i][j]; }
      for (int j=0; j<3; j++){
        for (int i=0; i<4; i++) HvvLambda_qsq[i][j] = (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j];
        HvvCLambda_qsq[j] = (selfDSpinZeroCoupl.HzzCLambda_qsq)[j];
      }
    }
    SetJHUGenSpinZeroVVCouplings(Hvvcoupl, HvvCLambda_qsq, HvvLambda_qsq, false);

    dXsec = VHiggsMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, includeHiggsDecay, verbosity);
    if (verbosity >= TVar::DEBUG) std::cout << "Process " << TVar::ProcessName(process) << " TEvtProb::XsecCalc_VX(): dXsec=" << dXsec << endl;
  } // end of JHUGen matrix element calculations
  else cerr << "Non-JHUGen VH MEs are not supported!" << endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  if (verbosity_>=TVar::DEBUG) cout << "End XsecCalc_VX"<< endl;
  return dXsec;
}

// Cross-section calculations for ttbar -> H
double TEvtProb::XsecCalc_TTX(
  TVar::Process process_, TVar::Production production_,
  TVar::VerbosityLevel verbosity_,
  int topProcess, int topDecay
  ){
  if (verbosity_>=TVar::DEBUG) cout << "Begin XsecCalc_TTX"<< endl;

  double dXsec = 0;
  ResetIORecord();
  if (!CheckInputPresent()) return dXsec;

  //Initialize Process
  SetProcess(process_);
  SetProduction(production_);
  SetVerbosity(verbosity_);

  // Set couplings at the qqH vertex
  double Hqqcoupl[SIZE_HQQ][2]={ { 0 } };

  if (matrixElement == TVar::JHUGen){
    // By default set the spin-0 couplings for SM case
    Hqqcoupl[0][0]=1.;  Hqqcoupl[0][1]=0.;   // first/second number is the real/imaginary part
    for (int i = 1; i<SIZE_HQQ; i++){ for (int com=0; com<2; com++) Hqqcoupl[i][com] = 0; }

    // 0-
    if (process == TVar::H0minus) {
      Hqqcoupl[0][0] = 0.;
      Hqqcoupl[1][0] = 1.;
    }
    else if (process == TVar::SelfDefine_spin0){
      for (int i=0; i<SIZE_HQQ; i++){ for (int j=0; j<2; j++) Hqqcoupl[i][j] = (selfDSpinZeroCoupl.Hqqcoupl)[i][j]; }
    }
    SetJHUGenSpinZeroQQCouplings(Hqqcoupl);

    if (production==TVar::ttH) dXsec = TTHiggsMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, topDecay, topProcess, verbosity);
    else if (production==TVar::bbH) dXsec = BBHiggsMatEl(process, production, matrixElement, &event_scales, &RcdME, EBEAM, topProcess, verbosity);
    if (verbosity >= TVar::DEBUG) std::cout << "Process " << TVar::ProcessName(process) << " TEvtProb::XsecCalc_TTX(): dXsec=" << dXsec << endl;
  }
  else cerr << "Non-JHUGen ttH MEs are not supported!" << endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  if (verbosity_>=TVar::DEBUG) cout << "End XsecCalc_TTX"<< endl;
  return dXsec;
}

double TEvtProb::GetXPropagator(
  TVar::ResonancePropagatorScheme scheme
  ){
  if (!CheckInputPresent()) return 0.;
  return TUtil::ResonancePropagator(melaCand->m(), scheme);
}




