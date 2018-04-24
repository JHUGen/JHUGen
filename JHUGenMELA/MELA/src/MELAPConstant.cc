#include <iostream>
#include <vector>
#include <cmath>
#include "MELAStreamHelpers.hh"
#include "MelaPConstant.h"
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"


using namespace std;
using TVar::simple_event_record;
using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;


MelaPConstant::MelaPConstant(
  TVar::MatrixElement me_,
  TVar::Production prod_,
  TVar::Process proc_,
  const char* path,
  const char* spname
  ) :
  processME(me_),
  processProd(prod_),
  processProc(proc_),
  fcnLow(0),
  fcnHigh(0),
  fcnMid(0),
  fname(path)
{
  GetFcnFromFile(path, spname);
}

MelaPConstant::~MelaPConstant(){}

void MelaPConstant::GetFcnFromFile(const char* path, const char* spname){
  TString spname_core = spname;
  spname_core.Append("_Smooth");
  spname_core.Prepend("sp_tg_");

  TFile* fin = TFile::Open(path, "read");
  gROOT->cd();

  if (fin!=0 && !fin->IsZombie() && fin->IsOpen()){
    TString spname_low = spname_core; spname_low.Prepend("lowFcn_");
    fcnLow = (TF1*)fin->Get(spname_low);

    TString spname_high = spname_core; spname_high.Prepend("highFcn_");
    fcnHigh = (TF1*)fin->Get(spname_high);

    TString spname_mid = spname_core;
    fcnMid = (TSpline3*)fin->Get(spname_mid);
  }
  else MELAerr << "MelaPConstant::GetFcnFromFile: Failed to open file in path " << path << endl;
  if (fin!=0 && fin->IsOpen()) fin->Close();
}

double MelaPConstant::Eval(const MelaIO* RcdME, TVar::VerbosityLevel verbosity)const{
  if (verbosity>=TVar::DEBUG) MELAout << "Begin MelaPConstant::Eval" << endl;

  double result=1;
  if (RcdME->melaCand->genStatus==-1) return result;

  double candMass = RcdME->melaCand->m();
  if (verbosity>=TVar::DEBUG) MELAout << "MelaPConstant::Eval: Candidate mass is " << candMass << endl;
  if (candMass<=0.) return 0.;
  else if (fcnLow!=0 && candMass<fcnLow->GetXmax()) result = fcnLow->Eval(candMass);
  else if (fcnHigh!=0 && candMass>fcnHigh->GetXmin()) result = fcnHigh->Eval(candMass);
  else if (fcnMid!=0){
    double var = candMass;
    if (var<fcnMid->GetXmin()) var = fcnMid->GetXmin();
    else if (var>fcnMid->GetXmax()) var = fcnMid->GetXmax();
    result = fcnMid->Eval(candMass);
  }
  else return result;

  if (verbosity>=TVar::DEBUG) MELAout << "MelaPConstant::Eval: Spline evaluated to " << result << endl;

  bool multiplyALR=false;
  bool multiplyHprop=false;
  bool multiplyAVjjprop=false;
  int sigmaZ4lprop=0;
  unsigned int powAlphaSatMZ=0;
  result = pow(10., result);
  if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::ZZINDEPENDENT
    &&
    processProc==TVar::HSMHiggs
    ){
    multiplyALR=true;
  }
  else if (
    (processME==TVar::JHUGen || processME==TVar::MCFM)
    &&
    processProd==TVar::ZZGG
    &&
    processProc==TVar::HSMHiggs
    ){
    multiplyALR=true;
    powAlphaSatMZ=2;
    multiplyHprop=true;
  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::ZZGG
    &&
    processProc==TVar::bkgZZ
    ){
    multiplyALR=true;
    powAlphaSatMZ=2;
  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::ZZQQB
    &&
    processProc==TVar::bkgZZ
    ){
    sigmaZ4lprop=4;
    multiplyALR=true;
  }
  else if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::JQCD
    &&
    processProc==TVar::HSMHiggs
    ){
    powAlphaSatMZ=3;
  }
  else if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::JJQCD
    &&
    processProc==TVar::HSMHiggs
    ){
    powAlphaSatMZ=4;
  }
  else if (
    processME==TVar::JHUGen
    &&
    (processProd==TVar::JJVBF || processProd==TVar::Had_WH || processProd==TVar::Had_ZH)
    &&
    processProc==TVar::HSMHiggs
    ){
    multiplyHprop=(processProd==TVar::Had_WH || processProd==TVar::Had_ZH);
    //multiplyAVjjprop=(processProd==TVar::Had_WH || processProd==TVar::Had_ZH); // Disabled to preserve discrimination power
  }
  else if (
    processME==TVar::MCFM
    &&
    (
    (processProc==TVar::HSMHiggs && (processProd==TVar::JJVBF_S || processProd==TVar::Had_WH_S || processProd==TVar::Had_ZH_S))
    ||
    (processProc==TVar::bkgZZ && (processProd==TVar::JJVBF || processProd==TVar::Had_WH || processProd==TVar::Had_ZH))
    )
    ){
    multiplyALR=true;
    multiplyHprop=(processProc==TVar::HSMHiggs);
    sigmaZ4lprop=4*(processProc==TVar::bkgZZ);
    //multiplyAVjjprop=(processProd==TVar::Had_WH || processProd==TVar::Had_ZH || processProd==TVar::Had_WH_S || processProd==TVar::Had_ZH_S); // Disabled to preserve discrimination power
  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::JJQCD
    &&
    processProc==TVar::bkgZZ
    ){
    multiplyALR=true;
    sigmaZ4lprop=4;
    powAlphaSatMZ=2;
  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::JJQCD
    &&
    processProc==TVar::bkgZJets
    ){
    multiplyALR=true;
  }
  else{
    if (verbosity>=TVar::DEBUG) MELAout << "MelaPConstant::Eval: Multiplying " << result << " by nothing." << endl;
  }

  if (multiplyALR) result *= this->GetVDaughterCouplings(RcdME, verbosity);
  if (multiplyHprop) result *= this->GetHPropagator(RcdME, verbosity);
  if (sigmaZ4lprop!=0) result *= this->GetZPropagator(RcdME, sigmaZ4lprop, verbosity);
  if (powAlphaSatMZ!=0) result *= this->GetAlphaSatMZ(RcdME, powAlphaSatMZ, verbosity);
  if (multiplyAVjjprop!=0) result *= this->GetAssociatedVjjPropagator(RcdME, verbosity);

  if (verbosity>=TVar::DEBUG) MELAout << "End MelaPConstant::Eval with " << result << endl;
  return result;
}

double MelaPConstant::GetHPropagator(const MelaIO* RcdME, const TVar::VerbosityLevel& verbosity)const{
  double propagator=1.;
  double sh = pow(RcdME->melaCand->m(), 2);
  double mh, gah;
  RcdME->getHiggsMassWidth(mh, gah, 0);
  propagator = 1./(pow(sh-pow(mh, 2), 2) + pow(mh*gah, 2));
  if (verbosity>=TVar::DEBUG) MELAout
    << "MelaPConstant::GetHPropagator: "
    << "propagatorH=" << propagator << " "
    << endl;
  return propagator;
}
double MelaPConstant::GetZPropagator(const MelaIO* RcdME, const int& sigmaZ4lprop, const TVar::VerbosityLevel& verbosity)const{
  double propagator=1.;
  double candMass = RcdME->melaCand->m();
  double sh = pow(candMass, 2);
  double mz = TUtil::GetMass(23);
  double gaz = TUtil::GetMass(23);
  double prop_sh = 1./(pow(sh-pow(mz, 2), 2) + pow(mz*gaz, 2));
  if (sigmaZ4lprop>0.){
    const double dsigmaZ4lprop=sigmaZ4lprop;
    if (fabs(candMass-mz)<=dsigmaZ4lprop*gaz){
      double shdn = pow(mz-dsigmaZ4lprop*gaz, 2);
      double shup = pow(mz+dsigmaZ4lprop*gaz, 2);
      double prop_shdn = 1./(pow(shdn-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double prop_shup = 1./(pow(shup-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double fsh = (sh-shdn)/(shup-shdn);
      propagator = prop_sh / (prop_shdn*(1.-fsh) + prop_shup*fsh);
    }
  }
  else if (sigmaZ4lprop<0.) propagator = prop_sh;
  if (verbosity>=TVar::DEBUG) MELAout
    << "MelaPConstant::GetZPropagator: "
    << "propagatorZ=" << propagator << " "
    << endl;
  return propagator;
}
double MelaPConstant::GetAssociatedVjjPropagator(const MelaIO* RcdME, const TVar::VerbosityLevel& verbosity)const{
  double propagator=1.;
  if (
    processProd==TVar::Had_WH || processProd==TVar::Had_WH_S
    ||
    processProd==TVar::Had_ZH || processProd==TVar::Had_ZH_S
    ){
    // Get a simple event record, safest way to handle jet mass corrections
    int nRequested_AssociatedJets=2;
    int AssociationVCompatibility=0;
    int partIncCode=TVar::kUseAssociated_Jets;
    double mv=0, gav=0;
    if (processProd==TVar::Had_WH || processProd==TVar::Had_WH_S){
      AssociationVCompatibility=24;
      mv = TUtil::GetMass(24);
      gav = TUtil::GetMass(24);
    }
    else{
      AssociationVCompatibility=23;
      mv = TUtil::GetMass(23);
      gav = TUtil::GetMass(23);
    }
    simple_event_record mela_event;
    mela_event.AssociationCode=partIncCode;
    mela_event.AssociationVCompatibility=AssociationVCompatibility;
    mela_event.nRequested_AssociatedJets=nRequested_AssociatedJets;
    TUtil::GetBoostedParticleVectors(
      RcdME->melaCand,
      mela_event,
      verbosity
      );

    float mJJval=-1;
    vector<TLorentzVector> pJets;
    const SimpleParticleCollection_t& pAssociated = mela_event.pAssociated;
    for (auto& part : pAssociated){
      if (PDGHelpers::isAJet(part.first)) pJets.push_back(part.second);
      if (pJets.size()==2) break;
    }
    if (pJets.size()==2) mJJval = (pJets[0] + pJets[1]).M();
    if (mJJval>=0.) propagator = 1./(pow(pow(mJJval, 2)-pow(mv, 2), 2) + pow(mv*gav, 2));
  }
  if (verbosity>=TVar::DEBUG) MELAout
    << "MelaPConstant::GetAssociatedVjjPropagator: "
    << "propagatorV=" << propagator << " "
    << endl;
  return propagator;
}
double MelaPConstant::GetVDaughterCouplings(const MelaIO* RcdME, const TVar::VerbosityLevel& verbosity)const{
  double result=1;
  double aL1, aR1, aL2, aR2;
  RcdME->getVDaughterCouplings(aL1, aR1, 0);
  RcdME->getVDaughterCouplings(aL2, aR2, 1);
  if (verbosity>=TVar::DEBUG) MELAout
    << "MelaPConstant::GetVDaughterCouplings: "
    << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
    << endl;
  if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
  if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);
  return result;
}
double MelaPConstant::GetAlphaSatMZ(const MelaIO* RcdME, const unsigned int& powAlphaSatMZ, const TVar::VerbosityLevel& verbosity)const{
  double result;
  double alphasVal = RcdME->getAlphaSatMZ();
  result = pow(alphasVal, powAlphaSatMZ);
  if (verbosity>=TVar::DEBUG) MELAout
    << "MelaPConstant::GetAlphaSatMZ: "
    << "alphas(MZ)=" << alphasVal << " (**" << powAlphaSatMZ << ") "
    << endl;
  return result;
}


