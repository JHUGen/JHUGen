#include "MELAStreamHelpers.hh"
#include "ZZMatrixElement.h"
#include "TLorentzRotation.h"


using namespace std;
using namespace TUtil;
using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;


ZZMatrixElement::ZZMatrixElement(
  const char* pathtoPDFSet,
  int PDFMember,
  const char* pathtoHiggsCSandWidth,
  double ebeam,
  TVar::VerbosityLevel verbosity
  ) :
  processVerbosity(verbosity),
  processLeptonInterference(TVar::DefaultLeptonInterf),
  EBEAM(ebeam),
  Xcal2(pathtoHiggsCSandWidth, EBEAM, pathtoPDFSet, PDFMember, verbosity),
  melaCand(0)
{
  if (processVerbosity>=TVar::DEBUG) MELAout << "Begin ZZMatrixElement constructor" << endl;
  build();
  if (processVerbosity>=TVar::DEBUG) MELAout << "End ZZMatrixElement constructor" << endl;
}
ZZMatrixElement::ZZMatrixElement(const ZZMatrixElement& other) :
processVerbosity(other.processVerbosity),
processLeptonInterference(other.processLeptonInterference),
EBEAM(other.EBEAM),
Xcal2(other.Xcal2),
melaCand(0) // 0 is correct in the copy constructor
{
  if (processVerbosity>=TVar::DEBUG) MELAout << "Begin ZZMatrixElement copy constructor" << endl;
  build();
  if (processVerbosity>=TVar::DEBUG) MELAout << "End ZZMatrixElement copy constructor" << endl;
}

void ZZMatrixElement::build(){
  if (processVerbosity>=TVar::DEBUG) MELAout << "Begin ZZMatrixElement::build" << endl;

  // Set default parameters explicitly
  set_PrimaryHiggsMass(125.);
  set_mHiggs(125., 0); set_wHiggs(-1., 0);
  set_mHiggs(-1., 1); set_wHiggs(-1, 1);

  selfD_SpinZeroCouplings = Xcal2.GetSelfDSpinZeroCouplings();
  selfD_SpinOneCouplings = Xcal2.GetSelfDSpinOneCouplings();
  selfD_SpinTwoCouplings = Xcal2.GetSelfDSpinTwoCouplings();
  selfD_VprimeCouplings = Xcal2.GetSelfDVprimeCouplings();
  selfD_aTQGCCouplings = Xcal2.GetSelfDaTQGCCouplings();

  if (processVerbosity>=TVar::DEBUG) MELAout << "End ZZMatrixElement::build" << endl;
}

ZZMatrixElement::~ZZMatrixElement(){
  if (processVerbosity>=TVar::DEBUG) MELAout << "Begin ZZMatrixElement destructor" << endl;
  resetPerEvent();
  reset_InputEvent();
  if (processVerbosity>=TVar::DEBUG) MELAout << "End ZZMatrixElement destructor" << endl;
}


std::vector<TLorentzVector> ZZMatrixElement::Calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi){
  double phi1, phi2;
  phi1=TMath::Pi()-Phi1;
  phi2=Phi1+Phi;

  double gamma1=1, gamma2=1, beta1=0, beta2=0;

  if (M1>0. && Mx>0.){
    gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
    beta1=sqrt(1.-1./(gamma1*gamma1));
  }
  if (M2>0. && Mx>0.){
    gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
    beta2=sqrt(1.-1./(gamma2*gamma2));
  }

  //gluon 4 vectors
  TLorentzVector p1CM(0, 0, Mx/2, Mx/2);
  TLorentzVector p2CM(0, 0, -Mx/2, Mx/2);

  //vector boson 4 vectors
  TLorentzVector kZ1(gamma1*M1*sin(theta)*beta1, 0, gamma1*M1*cos(theta)*beta1, gamma1*M1);
  TLorentzVector kZ2(-gamma2*M2*sin(theta)*beta2, 0, -gamma2*M2*cos(theta)*beta2, gamma2*M2);

  //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.

  TLorentzRotation Z1ToZ, Z2ToZ;

  Z1ToZ.Boost(0, 0, beta1);
  Z2ToZ.Boost(0, 0, beta2);
  Z1ToZ.RotateY(theta);
  Z2ToZ.RotateY(TMath::Pi()+theta);


  //fermion 4 vectors in vector boson rest frame

  TLorentzVector p3Z1((M1/2)*sin(theta1)*cos(phi1), (M1/2)*sin(theta1)*sin(phi1), (M1/2)*cos(theta1), (M1/2)*1);
  TLorentzVector p4Z1(-(M1/2)*sin(theta1)*cos(phi1), -(M1/2)*sin(theta1)*sin(phi1), -(M1/2)*cos(theta1), (M1/2)*1);
  TLorentzVector p5Z2((M2/2)*sin(theta2)*cos(phi2), (M2/2)*sin(theta2)*sin(phi2), (M2/2)*cos(theta2), (M2/2)*1);
  TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2), -(M2/2)*sin(theta2)*sin(phi2), -(M2/2)*cos(theta2), (M2/2)*1);

  // fermions 4 vectors in CM frame

  TLorentzVector p3CM, p4CM, p5CM, p6CM;

  p3CM=Z1ToZ*p3Z1;
  p4CM=Z1ToZ*p4Z1;
  p5CM=Z2ToZ*p5Z2;
  p6CM=Z2ToZ*p6Z2;

  vector<TLorentzVector> p;

  p.push_back(p3CM);
  p.push_back(p4CM);
  p.push_back(p5CM);
  p.push_back(p6CM);

  return p;
}

// Set-functions that set variables that belong to Xcal2 in addition to setting variables that belong to ZZMatrixElement
void ZZMatrixElement::set_Process(TVar::Process process_, TVar::MatrixElement me_, TVar::Production production_){
  processModel = process_;
  processME = me_;
  processProduction = production_;
  Xcal2.SetProcess(processModel, processME, processProduction);
}
void ZZMatrixElement::set_Verbosity(TVar::VerbosityLevel verbosity_){ processVerbosity = verbosity_; Xcal2.SetVerbosity(verbosity_); }
void ZZMatrixElement::set_LeptonInterference(TVar::LeptonInterference lepInterf_){ processLeptonInterference = lepInterf_; Xcal2.SetLeptonInterf(processLeptonInterference); }
// Set-functions that set variables that exclusively belong to Xcal2
void ZZMatrixElement::set_LHAgrid(const char* path, int pdfmember){ Xcal2.Set_LHAgrid(path, pdfmember); }
void ZZMatrixElement::set_RenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf){
  Xcal2.SetRenFacScaleMode(renormalizationSch, factorizationSch, ren_sf, fac_sf);
}
void ZZMatrixElement::set_CandidateDecayMode(TVar::CandidateDecayMode mode){ Xcal2.SetCandidateDecayMode(mode); }
void ZZMatrixElement::set_PrimaryHiggsMass(double mh){ Xcal2.SetPrimaryHiggsMass(mh); }
void ZZMatrixElement::set_CurrentCandidateFromIndex(unsigned int icand){ Xcal2.SetCurrentCandidateFromIndex(icand); }
void ZZMatrixElement::set_CurrentCandidate(MELACandidate* cand){ Xcal2.SetCurrentCandidate(cand); }
void ZZMatrixElement::set_InputEvent(
  SimpleParticleCollection_t* pDaughters,
  SimpleParticleCollection_t* pAssociated,
  SimpleParticleCollection_t* pMothers,
  bool isGen
  ){ // Adds a new candidate to Xcal2
  Xcal2.SetInputEvent(
    pDaughters,
    pAssociated,
    pMothers,
    isGen
    );
}
// Sets melaCand in Xcal2 to a temporary candidate, without pushing this candidate to candList of Xcal2
void ZZMatrixElement::set_TempCandidate(
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
    &tmpPartList, &tmpCandList // push_back is done automatically
    );
  if (cand){
    melaCand=cand;
    set_CurrentCandidate(melaCand);
  }
}
// Adds a top candidate
void ZZMatrixElement::append_TopCandidate(SimpleParticleCollection_t* TopDaughters){ Xcal2.AppendTopCandidate(TopDaughters); }
// Set-functions that do not set anything that belongs to Xcal2
void ZZMatrixElement::set_mHiggs(double mh_, int index){
  if (index<nSupportedHiggses && index>=0) mHiggs[index] = mh_;
  else MELAerr << "ZZMatrixElement::set_mHiggs: Only resonances 0 (regular) and 1 (additional, possibly high-mass) are supported" << endl;
}
void ZZMatrixElement::set_wHiggs(double gah_, int index){
  if (index<nSupportedHiggses && index>=0) wHiggs[index] = (double)gah_;
  else MELAerr << "ZZMatrixElement::set_wHiggs: Only resonances 0 (regular) and 1 (additional, possibly high-mass) are supported" << endl;
}
void ZZMatrixElement::set_mHiggs_wHiggs(double mh_, double gah_, int index){
  if (index<nSupportedHiggses && index>=0){
    mHiggs[index] = mh_;
    wHiggs[index] = gah_;
  }
  else MELAerr << "ZZMatrixElement::set_mHiggs_wHiggs: Only resonances 0 (regular) and 1 (additional, possibly high-mass) are supported" << endl;
}

// reset_MCFM_EWKParameters resets the MCFM EW parameters to those specified. This function is a wrapper around the TEvtProb version.
void ZZMatrixElement::reset_Mass(double inmass, int ipart){ Xcal2.ResetMass(inmass, ipart); }
void ZZMatrixElement::reset_Width(double inwidth, int ipart){ Xcal2.ResetWidth(inwidth, ipart); }
void ZZMatrixElement::reset_QuarkMasses(){ Xcal2.ResetQuarkMasses(); }
void ZZMatrixElement::reset_MCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme){
  Xcal2.ResetMCFM_EWKParameters(ext_Gf, ext_aemmz, ext_mW, ext_mZ, ext_xW, ext_ewscheme);
}
//
// resetPerEvent resets the mass, width and lepton interference settings and deletes the temporary input objects ZZMatrixElement owns.
void ZZMatrixElement::resetPerEvent(){
  // Protection against forgetfulness; custom width has to be set per-computation
  set_mHiggs(Xcal2.GetPrimaryHiggsMass(), 0); // Sets mHiggs[0]
  if (wHiggs[0]>=0.) set_wHiggs(-1., 0);

  if (mHiggs[1]>=0.) set_mHiggs(-1., 1);
  if (wHiggs[1]>=0.) set_wHiggs(-1., 1);

  Xcal2.SetHiggsMass(mHiggs[0], -1, -1);

  // Return back to default lepton interference settings after each calculation
  if (processLeptonInterference != TVar::DefaultLeptonInterf) set_LeptonInterference(TVar::DefaultLeptonInterf);

  // Delete the temporary input objects owned
  for (unsigned int ic=0; ic<tmpCandList.size(); ic++){ if (tmpCandList.at(ic)) delete tmpCandList.at(ic); }
  //for (unsigned int itc=0; itc<tmpTopCandList.size(); itc++){ if (tmpTopCandList.at(itc)) delete tmpTopCandList.at(itc); }
  for (unsigned int ip=0; ip<tmpPartList.size(); ip++){ if (tmpPartList.at(ip)) delete tmpPartList.at(ip); }
  melaCand=nullptr;
}
// Resets all candidates in Xcal2, to be called at the end of each event after all computations are done
void ZZMatrixElement::reset_InputEvent(){ Xcal2.ResetInputEvent(); }


MelaIO* ZZMatrixElement::get_IORecord(){ return Xcal2.GetIORecord(); }
double ZZMatrixElement::get_PrimaryMass(int ipart){ return Xcal2.GetPrimaryMass(ipart); }
double ZZMatrixElement::get_PrimaryWidth(int ipart){ return Xcal2.GetPrimaryWidth(ipart); }
double ZZMatrixElement::get_HiggsWidthAtPoleMass(double mass){ return Xcal2.GetHiggsWidthAtPoleMass(mass); }
MELACandidate* ZZMatrixElement::get_CurrentCandidate(){ return Xcal2.GetCurrentCandidate(); }
int ZZMatrixElement::get_CurrentCandidateIndex(){ return Xcal2.GetCurrentCandidateIndex(); }
int ZZMatrixElement::get_NCandidates(){ return Xcal2.GetNCandidates(); }
std::vector<MELATopCandidate_t*>* ZZMatrixElement::get_TopCandidateCollection(){ return Xcal2.GetTopCandidates(); }


void ZZMatrixElement::set_SpinZeroCouplings(
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
  bool diffHWW
  ){
  Xcal2.AllowSeparateWWCouplings(diffHWW);
  for (int jh=1; jh<=(int)nSupportedHiggses; jh++){
    for (int ic=0; ic<SIZE_HGG; ic++) selfD_SpinZeroCouplings->SetHGGCouplings(ic, selfDHggcoupl[jh-1][ic][0], selfDHggcoupl[jh-1][ic][1], 1, jh);
    for (int ic=0; ic<SIZE_HGG; ic++) selfD_SpinZeroCouplings->SetHGGCouplings(ic, selfDHg4g4coupl[jh-1][ic][0], selfDHg4g4coupl[jh-1][ic][1], 2, jh);

    for (int ic=0; ic<SIZE_HQQ; ic++) selfD_SpinZeroCouplings->SetHQQCouplings(ic, selfDHqqcoupl[jh-1][ic][0], selfDHqqcoupl[jh-1][ic][1], 0, jh);
    for (int ic=0; ic<SIZE_HQQ; ic++) selfD_SpinZeroCouplings->SetHQQCouplings(ic, selfDHbbcoupl[jh-1][ic][0], selfDHbbcoupl[jh-1][ic][1], 5, jh);
    for (int ic=0; ic<SIZE_HQQ; ic++) selfD_SpinZeroCouplings->SetHQQCouplings(ic, selfDHttcoupl[jh-1][ic][0], selfDHttcoupl[jh-1][ic][1], 6, jh);
    for (int ic=0; ic<SIZE_HQQ; ic++) selfD_SpinZeroCouplings->SetHQQCouplings(ic, selfDHb4b4coupl[jh-1][ic][0], selfDHb4b4coupl[jh-1][ic][1], 7, jh);
    for (int ic=0; ic<SIZE_HQQ; ic++) selfD_SpinZeroCouplings->SetHQQCouplings(ic, selfDHt4t4coupl[jh-1][ic][0], selfDHt4t4coupl[jh-1][ic][1], 8, jh);

    for (int ic=0; ic<SIZE_HVV; ic++) selfD_SpinZeroCouplings->SetHVVCouplings(ic, selfDHzzcoupl[jh-1][ic][0], selfDHzzcoupl[jh-1][ic][1], false, jh);
    for (int ic=0; ic<SIZE_HVV; ic++) selfD_SpinZeroCouplings->SetHVVCouplings(ic, selfDHwwcoupl[jh-1][ic][0], selfDHwwcoupl[jh-1][ic][1], true, jh);
    for (int ik=0; ik<SIZE_HVV_CQSQ; ik++){
      for (int ig=0; ig<SIZE_HVV_LAMBDAQSQ; ig++){
        selfD_SpinZeroCouplings->SetHVVLambdaQ2(ig, ik, selfDHzzLambda_qsq[jh-1][ig][ik], false, jh);
        selfD_SpinZeroCouplings->SetHVVLambdaQ2(ig, ik, selfDHwwLambda_qsq[jh-1][ig][ik], true, jh);
      }
      selfD_SpinZeroCouplings->SetHVVSignCQ2(ik, selfDHzzCLambda_qsq[jh-1][ik], false, jh);
      selfD_SpinZeroCouplings->SetHVVSignCQ2(ik, selfDHwwCLambda_qsq[jh-1][ik], true, jh);
    }
  }
}
void ZZMatrixElement::set_SpinZeroContact(
  double selfDHzzpcoupl[SIZE_HVV][2],
  double selfDHzpzpcoupl[SIZE_HVV][2],
  double selfDHwwpcoupl[SIZE_HVV][2],
  double selfDHwpwpcoupl[SIZE_HVV][2]
  ){
  for (int ic=0; ic<SIZE_HVV; ic++){
    selfD_SpinZeroCouplings->SetHVVpCouplings(ic, selfDHzzpcoupl[ic][0], selfDHzzpcoupl[ic][1]);
    selfD_SpinZeroCouplings->SetHVpVpCouplings(ic, selfDHzpzpcoupl[ic][0], selfDHzpzpcoupl[ic][1]);
    selfD_SpinZeroCouplings->SetHVVpCouplings(ic, selfDHwwpcoupl[ic][0], selfDHwwpcoupl[ic][1], true);
    selfD_SpinZeroCouplings->SetHVpVpCouplings(ic, selfDHwpwpcoupl[ic][0], selfDHwpwpcoupl[ic][1], true);
  }
}
void ZZMatrixElement::set_SpinOneCouplings(
  double selfDZqqcoupl[SIZE_ZQQ][2],
  double selfDZvvcoupl[SIZE_ZVV][2]
  ){
  for (int ic=0; ic<SIZE_ZQQ; ic++) selfD_SpinOneCouplings->SetZQQCouplings(ic, selfDZqqcoupl[ic][0], selfDZqqcoupl[ic][1]);
  for (int ic=0; ic<SIZE_ZVV; ic++) selfD_SpinOneCouplings->SetZVVCouplings(ic, selfDZvvcoupl[ic][0], selfDZvvcoupl[ic][1]);
}
void ZZMatrixElement::set_SpinTwoCouplings(
  double selfDGqqcoupl[SIZE_GQQ][2],
  double selfDGggcoupl[SIZE_GGG][2],
  double selfDGvvcoupl[SIZE_GVV][2]
  ){
  for (int ic=0; ic<SIZE_GQQ; ic++) selfD_SpinTwoCouplings->SetGQQCouplings(ic, selfDGqqcoupl[ic][0], selfDGqqcoupl[ic][1]);
  for (int ic=0; ic<SIZE_GGG; ic++) selfD_SpinTwoCouplings->SetGGGCouplings(ic, selfDGggcoupl[ic][0], selfDGggcoupl[ic][1]);
  for (int ic=0; ic<SIZE_GVV; ic++) selfD_SpinTwoCouplings->SetGVVCouplings(ic, selfDGvvcoupl[ic][0], selfDGvvcoupl[ic][1]);
}
void ZZMatrixElement::set_SpinTwoContact(
  double selfDGvvpcoupl[SIZE_GVV][2],
  double selfDGvpvpcoupl[SIZE_GVV][2]
){
  for (int ic=0; ic<SIZE_GVV; ic++){
    selfD_SpinTwoCouplings->SetGVVpCouplings(ic, selfDGvvpcoupl[ic][0], selfDGvvpcoupl[ic][1]);
    selfD_SpinTwoCouplings->SetGVpVpCouplings(ic, selfDGvpvpcoupl[ic][0], selfDGvpvpcoupl[ic][1]);
  }
}
void ZZMatrixElement::set_VprimeContactCouplings(
  double selfDZpffcoupl[SIZE_Vpff][2],
  double selfDWpffcoupl[SIZE_Vpff][2],
  double M_Zprime,
  double Ga_Zprime,
  double M_Wprime,
  double Ga_Wprime
){
  for (int ic=0; ic<SIZE_Vpff; ic++){
    selfD_VprimeCouplings->SetVpffCouplings(ic, selfDZpffcoupl[ic][0], selfDZpffcoupl[ic][1]);
    selfD_VprimeCouplings->SetVpffCouplings(ic, selfDWpffcoupl[ic][0], selfDWpffcoupl[ic][1], true);
  }
  selfD_VprimeCouplings->SetZPrimeMassWidth(M_Zprime, Ga_Zprime);
  selfD_VprimeCouplings->SetWPrimeMassWidth(M_Wprime, Ga_Wprime);
}
void ZZMatrixElement::set_aTQGCCouplings(
  double selfDaTQGCcoupl[SIZE_ATQGC][2]
){
  for (int ic=0; ic<SIZE_ATQGC; ic++) selfD_aTQGCCouplings->SetATQGCCouplings(ic, selfDaTQGCcoupl[ic][0], selfDaTQGCcoupl[ic][1]);
}



// Higgs + 0 jets dedicated function (with no Higgs decay)
void ZZMatrixElement::computeXS(
  float &mevalue
  ){
  melaCand = get_CurrentCandidate();

  if (melaCand){
    double zzmass = melaCand->m();
    if (processME == TVar::MCFM){
      for (int jh=0; jh<(int)nSupportedHiggses; jh++) Xcal2.SetHiggsMass(mHiggs[jh], wHiggs[jh], jh+1);
    }
    else Xcal2.SetHiggsMass(zzmass, wHiggs[0], -1);

    mevalue = Xcal2.XsecCalc_XVV();
  }

  resetPerEvent();
  return;
}

// VBF+VH dedicated function (production(+)decay)
void ZZMatrixElement::computeProdXS_VVHVV(
  float& mevalue
  ){
  melaCand = get_CurrentCandidate();

  if (melaCand){
    double zzmass = melaCand->m();
    if (processME == TVar::MCFM){
      for (int jh=0; jh<(int)nSupportedHiggses; jh++) Xcal2.SetHiggsMass(mHiggs[jh], wHiggs[jh], jh+1);
    }
    else Xcal2.SetHiggsMass(zzmass, wHiggs[0], -1);

    mevalue = Xcal2.XsecCalc_VVXVV();
  }

  resetPerEvent();
  return;
}

// Higgs + 2 jet dedicated function (with no Higgs decay)
void ZZMatrixElement::computeProdXS_JJH(
  float &mevalue
  ){
  melaCand = get_CurrentCandidate();

  if (melaCand){
    mevalue  = Xcal2.XsecCalcXJJ();
  }

  resetPerEvent();
  return;
}

// Higgs + 1 jet: Only SM is supported for now.
void ZZMatrixElement::computeProdXS_JH(
  float &mevalue
  ){
  melaCand = get_CurrentCandidate();

  if (melaCand){
    mevalue  = Xcal2.XsecCalcXJ();
  }

  resetPerEvent();
  return;
}

// Dedicated VH function (with no Higgs decay)
void ZZMatrixElement::computeProdXS_VH(
  float &mevalue,
  bool includeHiggsDecay
  ){
  melaCand = get_CurrentCandidate();

  if (melaCand){
    double zzmass = melaCand->m();
    if (processME == TVar::MCFM){ for (int jh=0; jh<(int)nSupportedHiggses; jh++) Xcal2.SetHiggsMass(mHiggs[jh], wHiggs[jh], jh+1); }
    else Xcal2.SetHiggsMass(zzmass, wHiggs[0], -1);

    mevalue  = Xcal2.XsecCalc_VX(
      includeHiggsDecay
      );
  }

  resetPerEvent();
  return;
}

// Dedicated ttH/bbH function (with no Higgs decay)
// topProcess=0 for gg, =1 for qqb initial states
// topDecay=0 (default) for no top decay, =1 to include t/tb->b/bb + W+/W-(->ffb). =1 not relevant for the bbH process.
void ZZMatrixElement::computeProdXS_ttH(
  float &mevalue,
  int topProcess,
  int topDecay
  ){
  melaCand = get_CurrentCandidate();

  if (melaCand){
    double zzmass = melaCand->m();
    if (processME == TVar::MCFM){ for (int jh=0; jh<(int)nSupportedHiggses; jh++) Xcal2.SetHiggsMass(mHiggs[jh], wHiggs[jh], jh+1); }
    else Xcal2.SetHiggsMass(zzmass, wHiggs[0], -1);

    if (processProduction == TVar::ttH || processProduction == TVar::bbH){
      mevalue  = Xcal2.XsecCalc_TTX(
        topProcess, topDecay
        );
    }
  }

  resetPerEvent();
  return;
}

// Higgs propagator
void ZZMatrixElement::get_XPropagator(TVar::ResonancePropagatorScheme scheme, float& prop){
  prop=0.;
  melaCand = get_CurrentCandidate();

  if (melaCand){
    Xcal2.SetHiggsMass(mHiggs[0], wHiggs[0], -1);
    prop=Xcal2.GetXPropagator(scheme);
  }

  resetPerEvent();
  return;
}


