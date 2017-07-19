#include <iostream>
#include <vector>
#include <cmath>
#include "MelaPConstant.h"
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"


using namespace std;
using TVar::simple_event_record;


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
  else cerr << "MelaPConstant::GetFcnFromFile: Failed to open file in path " << path << endl;
  if (fin!=0 && fin->IsOpen()) fin->Close();
}

double MelaPConstant::Eval(const MelaIO* RcdME, TVar::VerbosityLevel verbosity)const{
  if (verbosity>=TVar::DEBUG) cout << "Begin MelaPConstant::Eval" << endl;

  double result=1;
  if (RcdME->melaCand->genStatus==-1) return result;

  double candMass = RcdME->melaCand->m();
  if (verbosity>=TVar::DEBUG) cout << "MelaPConstant::Eval: Candidate mass is " << candMass << endl;
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

  if (verbosity>=TVar::DEBUG) cout << "MelaPConstant::Eval: Spline evaluated to " << result << endl;

  if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::ZZINDEPENDENT
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
      << endl;

    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);

  }
  else if (
    (processME==TVar::JHUGen || processME==TVar::MCFM)
    &&
    processProd==TVar::ZZGG
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();
    double propagator, mh, gah;
    RcdME->getHiggsMassWidth(mh, gah, 0);
    propagator = 1./(pow(pow(candMass, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "alphas(MZ)=" << alphasVal << "**2 "
      << "propagator=" << propagator << " "
      << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
      << endl;

    result *= pow(alphasVal, 2);
    result *= propagator;
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);

  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::ZZGG
    &&
    processProc==TVar::bkgZZ
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();
    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "alphas(MZ)=" << alphasVal << "**2 "
      << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
      << endl;

    result *= pow(alphasVal, 2);
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);

  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::ZZQQB
    &&
    processProc==TVar::bkgZZ
    ){
    result = pow(10., result);

    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);

    double propagator;
    double mz = TUtil::GetMass(23);
    double gaz = TUtil::GetMass(23);
    if (fabs(candMass-mz)<=4.*gaz){
      double sh = pow(candMass, 2);
      double shdn = pow(mz-4.*gaz, 2);
      double shup = pow(mz+4.*gaz, 2);
      double prop_sh = 1./(pow(sh-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double prop_shdn = 1./(pow(shdn-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double prop_shup = 1./(pow(shup-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double fsh = (sh-shdn)/(shup-shdn);
      propagator = prop_sh / (prop_shdn*(1.-fsh) + prop_shup*fsh);
    }
    else propagator=1.;

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "propagator=" << propagator << " "
      << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
      << endl;

    result *= propagator;
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);

  }
  else if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::JQCD
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "alphas(MZ)=" << alphasVal << "**3 "
      << endl;

    result *= pow(alphasVal, 3);
  }
  else if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::JJQCD
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "alphas(MZ)=" << alphasVal << "**4 "
      << endl;

    result *= pow(alphasVal, 4);

  }
  else if (
    processME==TVar::JHUGen
    &&
    (processProd==TVar::JJVBF || processProd==TVar::Had_WH || processProd==TVar::Had_ZH)
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double propagator=1.;
    double propagatorV=1;
    if (processProd==TVar::Had_WH || processProd==TVar::Had_ZH){
      double mh, gah;
      RcdME->getHiggsMassWidth(mh, gah, 0);
      propagator = 1./(pow(pow(candMass, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));

      // Get a simple event record, safest way to handle jet mass corrections
      int nRequested_AssociatedJets=2;
      int AssociationVCompatibility=0;
      int partIncCode=TVar::kUseAssociated_Jets;
      double mv=0, gav=0;
      if (processProd==TVar::Had_WH){
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
      if (mJJval>=0.) propagatorV = 1./(pow(pow(mJJval, 2)-pow(mv, 2), 2) + pow(mv*gav, 2));
    }

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "propagator=" << propagator << " "
      << "propagatorV=" << propagatorV << " "
      << endl;

    result *= propagator*propagatorV;

  }
  else if (
    processME==TVar::MCFM
    &&
    (processProd==TVar::JJVBF || processProd==TVar::Had_WH || processProd==TVar::Had_ZH)
    &&
    (processProc==TVar::HSMHiggs || processProc==TVar::bkgZZ)
    ){
    result = pow(10., result);

    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);

    double propagator=1.;
    if (processProc==TVar::HSMHiggs){
      double mh, gah;
      RcdME->getHiggsMassWidth(mh, gah, 0);
      propagator = 1./(pow(pow(candMass, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
    }
    /* FIXME: MAY NEED TO ENABLE THIS
    else{
      double propagator;
      double mz = TUtil::GetMass(23);
      double gaz = TUtil::GetMass(23);
      if (fabs(candMass-mz)<=4.*gaz){
        double sh = pow(candMass, 2);
        double shdn = pow(mz-4.*gaz, 2);
        double shup = pow(mz+4.*gaz, 2);
        double prop_sh = 1./(pow(sh-pow(mz, 2), 2) + pow(mz*gaz, 2));
        double prop_shdn = 1./(pow(shdn-pow(mz, 2), 2) + pow(mz*gaz, 2));
        double prop_shup = 1./(pow(shup-pow(mz, 2), 2) + pow(mz*gaz, 2));
        double fsh = (sh-shdn)/(shup-shdn);
        propagator = prop_sh / (prop_shdn*(1.-fsh) + prop_shup*fsh);
      }
      else propagator=1.;
    }
    */
    double propagatorV=1;
    if (processProd==TVar::Had_WH || processProd==TVar::Had_ZH){
      // Get a simple event record, safest way to handle jet mass corrections
      int nRequested_AssociatedJets=2;
      int AssociationVCompatibility=0;
      int partIncCode=TVar::kUseAssociated_Jets;
      double mv=0, gav=0;
      if (processProd==TVar::Had_WH){
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
      if (mJJval>=0.) propagatorV = 1./(pow(pow(mJJval, 2)-pow(mv, 2), 2) + pow(mv*gav, 2));
    }

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "propagator=" << propagator << " "
      << "propagatorV=" << propagatorV << " "
      << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
      << endl;

    result *= propagator*propagatorV;
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);

  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::JJQCD
    &&
    processProc==TVar::bkgZZ
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();

    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);

    double propagator;
    double mz = TUtil::GetMass(23);
    double gaz = TUtil::GetMass(23);
    if (fabs(candMass-mz)<=4.*gaz){
      double sh = pow(candMass, 2);
      double shdn = pow(mz-4.*gaz, 2);
      double shup = pow(mz+4.*gaz, 2);
      double prop_sh = 1./(pow(sh-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double prop_shdn = 1./(pow(shdn-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double prop_shup = 1./(pow(shup-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double fsh = (sh-shdn)/(shup-shdn);
      propagator = prop_sh / (prop_shdn*(1.-fsh) + prop_shup*fsh);
    }
    else propagator=1.;

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "propagator=" << propagator << " "
      << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
      << "alphas(MZ)=" << alphasVal << "**2 "
      << endl;

    result *= propagator;
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);
    result *= pow(alphasVal, 2);

  }
  else{
    result = pow(10., result);
    if (verbosity>=TVar::DEBUG) cout << "MelaPConstant::Eval: Multiplying " << result << " by nothing." << endl;
  }

  if (verbosity>=TVar::DEBUG) cout << "End MelaPConstant::Eval with " << result << endl;
  return result;
}

double MelaPConstant::GetHPropagator(const MelaIO* RcdME, const TVar::VerbosityLevel& verbosity)const{
  double propagator=1.;
  double sh = pow(RcdME->melaCand->m(), 2);
  double mh, gah;
  RcdME->getHiggsMassWidth(mh, gah, 0);
  propagator = 1./(pow(sh-pow(mh, 2), 2) + pow(mh*gah, 2));
  if (verbosity>=TVar::DEBUG) cout
    << "MelaPConstant::GetHPropagator: "
    << "propagatorH=" << propagator << " "
    << endl;
  return propagator;
}
double MelaPConstant::GetZPropagator(const MelaIO* RcdME, const bool restricted, const TVar::VerbosityLevel& verbosity)const{
  double propagator=1.;
  double candMass = RcdME->melaCand->m();
  double sh = pow(candMass, 2);
  double mz = TUtil::GetMass(23);
  double gaz = TUtil::GetMass(23);
  double prop_sh = 1./(pow(sh-pow(mz, 2), 2) + pow(mz*gaz, 2));
  if (restricted){
    if (fabs(candMass-mz)<=4.*gaz){
      double shdn = pow(mz-4.*gaz, 2);
      double shup = pow(mz+4.*gaz, 2);
      double prop_shdn = 1./(pow(shdn-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double prop_shup = 1./(pow(shup-pow(mz, 2), 2) + pow(mz*gaz, 2));
      double fsh = (sh-shdn)/(shup-shdn);
      propagator = prop_sh / (prop_shdn*(1.-fsh) + prop_shup*fsh);
    }
  }
  else propagator = prop_sh;
  if (verbosity>=TVar::DEBUG) cout
    << "MelaPConstant::GetZPropagator: "
    << "propagatorZ=" << propagator << " "
    << endl;
  return propagator;
}
double MelaPConstant::GetAssociatedVjjPropagator(const MelaIO* RcdME, const TVar::VerbosityLevel& verbosity)const{
  double propagator=1.;
  if (processProd==TVar::Had_WH || processProd==TVar::Had_ZH){
    // Get a simple event record, safest way to handle jet mass corrections
    int nRequested_AssociatedJets=2;
    int AssociationVCompatibility=0;
    int partIncCode=TVar::kUseAssociated_Jets;
    double mv=0, gav=0;
    if (processProd==TVar::Had_WH){
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
  if (verbosity>=TVar::DEBUG) cout
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
  if (verbosity>=TVar::DEBUG) cout
    << "MelaPConstant::GetVDaughterCouplings: "
    << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
    << endl;
  if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
  if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);
  return result;
}
double MelaPConstant::GetAlphaSatMZ(const MelaIO* RcdME, const int p, const TVar::VerbosityLevel& verbosity)const{
  double result;
  double alphasVal = RcdME->getAlphaSatMZ();
  result = pow(alphasVal, p);
  if (verbosity>=TVar::DEBUG) cout
    << "MelaPConstant::GetAlphaSatMZ: "
    << "alphas(MZ)=" << alphasVal << " (**" << p << ") "
    << endl;
  return result;
}


