#include <iostream>
#include <cassert>
#include "MELAStreamHelpers.hh"
#include "MELACandidate.h"
#include "SuperDijetMela.h"


using namespace std;
using TVar::simple_event_record;
using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;


SuperDijetMela::SuperDijetMela(float sqrts_, TVar::VerbosityLevel verbosity_) :
sqrts(sqrts_),
verbosity(verbosity_)
{
  Build();
}
SuperDijetMela::SuperDijetMela(const SuperDijetMela& other) : 
sqrts(other.sqrts),
verbosity(other.verbosity)
{
  Build();
}
SuperDijetMela::~SuperDijetMela(){
  for (auto it=ResolutionModelMap.begin(); it!=ResolutionModelMap.end(); it++) delete it->second;
}

void SuperDijetMela::Build(){
  SetupResolutionModel(TVar::Had_ZH);
  SetupResolutionModel(TVar::Had_WH);
}

void SuperDijetMela::SetupResolutionModel(TVar::Production prod){
  // Setup file path
#ifdef _melapkgpathstr_
  const string MELAPKGPATH = _melapkgpathstr_;
#else
  MELAout << "SuperDijetMela::SetupResolutionModel: MELA package path is undefined! Please modify the makefle or the makefile-equivalent!" << endl;
  assert(0);
#endif
  TString path = TString(MELAPKGPATH.c_str()) + "data/resolution_mJJ_recoVStrue_";
  TString prodName;
  switch (prod){
  case TVar::Had_ZH:
    prodName = "ZH";
    break;
  case TVar::Had_WH:
    prodName = "WH";
    break;
  default:
    MELAout << "SuperDijetMela::SetupResolutionModel: Production " << TVar::ProductionName(prod) << " is unknown." << endl;
    return;
  }
  path += prodName;
  path += Form("_%.0fTeV%s", sqrts, ".root");

  TString appendName = Form("mJJReso_%.0fTeV", sqrts);
  MELADifermionResolutionModel* model = new MELADifermionResolutionModel(prod, sqrts, path, appendName);
  if (model->isValid()){
    int iprod = (int)prod;
    ResolutionModelMap[iprod] = model;
  }
  else MELAerr << "SuperDijetMela::SetupResolutionModel: Model for production " << TVar::ProductionName(prod) << " cannot be built." << endl;
}
float SuperDijetMela::GetConvBW(TVar::Production prod, MELACandidate* cand, bool useTrueBW){
  float result=-1;
  int iprod = (int)prod;
  if (cand!=0){
    float mJJval=-1;
    bool modelExists = (ResolutionModelMap.find(iprod)!=ResolutionModelMap.end());
    int AssociationVCompatibility=0;

    if (modelExists || useTrueBW){
      // Get a simple event record, safest way to handle jet mass corrections
      int nRequested_AssociatedJets=0;
      int partIncCode=TVar::kNoAssociated; // Just to avoid warnings
      if (prod == TVar::Had_ZH || prod == TVar::Had_WH){ // Only use associated partons
        partIncCode=TVar::kUseAssociated_Jets;
        nRequested_AssociatedJets=2;
      }
      if (prod==TVar::Had_WH) AssociationVCompatibility=24;
      else if (prod==TVar::Had_ZH) AssociationVCompatibility=23;
      simple_event_record mela_event;
      mela_event.AssociationCode=partIncCode;
      mela_event.AssociationVCompatibility=AssociationVCompatibility;
      mela_event.nRequested_AssociatedJets=nRequested_AssociatedJets;
      TUtil::GetBoostedParticleVectors(
        cand,
        mela_event,
        verbosity
        );

      vector<TLorentzVector> pJets;
      const SimpleParticleCollection_t& pAssociated = mela_event.pAssociated;
      for (auto& part : pAssociated){
        if (PDGHelpers::isAJet(part.first)) pJets.push_back(part.second);
        if (pJets.size()==2) break;
      }
      if (pJets.size()==2) mJJval = (pJets[0] + pJets[1]).M();
    }

    if (useTrueBW){
      float truePoleMass=-1;
      float truePoleWidth=-1;

      if (AssociationVCompatibility!=0){
        truePoleMass=TUtil::GetMass(AssociationVCompatibility);
        truePoleWidth=TUtil::GetDecayWidth(AssociationVCompatibility);
        if (mJJval>=0.) result = 1./(pow(pow(mJJval, 2)-pow(truePoleMass, 2), 2) + pow(truePoleMass*truePoleWidth, 2));
      }

      if (verbosity>=TVar::DEBUG) MELAout
        << "SuperDijetMela::GetConvBW[idV=" << AssociationVCompatibility << "]::"
        << "sqrt(s) = " << mJJval
        << ", true (m,Gamma) = ( " << truePoleMass << " , " << truePoleWidth << " )"
        << ", ideal BW = " << result
        << endl;
    }
    else if (modelExists){
      result = ResolutionModelMap[iprod]->getVal(mJJval);
      if (verbosity>=TVar::DEBUG) MELAout
        << "SuperDijetMela::GetConvBW[" << TVar::ProductionName(prod) << "]::"
        << "sqrt(s) = " << mJJval
        << ", reco BW = " << result
        << endl;
    }
    else{
      if (verbosity>=TVar::DEBUG) MELAout
        << "SuperDijetMela::GetConvBW:: No method for " << TVar::ProductionName(prod) << " is known. "
        << "sqrt(s) = " << mJJval
        << ", return value = " << result
        << endl;
    }
  }
  return result;
}
