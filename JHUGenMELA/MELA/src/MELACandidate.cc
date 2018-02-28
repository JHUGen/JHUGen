#include "MELACandidate.h"
#include "TUtilHelpers.hh"
#include "MELAStreamHelpers.hh"
#include "TMath.h"


using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;
using namespace PDGHelpers;


MELACandidate::MELACandidate() :
MELAParticle(),
associatedByHighestPt(false),
isShallowCopy(false),
selfDecayMode(TVar::CandidateDecay_Stable)
{}
MELACandidate::MELACandidate(int id_, bool associatedByHighestPt_) :
MELAParticle(id_),
associatedByHighestPt(associatedByHighestPt_),
isShallowCopy(false),
selfDecayMode(TVar::CandidateDecay_Stable)
{}
MELACandidate::MELACandidate(int id_, TLorentzVector p4_, bool associatedByHighestPt_) :
MELAParticle(id_, p4_),
associatedByHighestPt(associatedByHighestPt_),
isShallowCopy(false),
selfDecayMode(TVar::CandidateDecay_Stable)
{}
MELACandidate::MELACandidate(const MELACandidate& particle_) :
MELAParticle(particle_),
associatedByHighestPt(particle_.associatedByHighestPt),
isShallowCopy(false),
selfDecayMode(particle_.selfDecayMode),
associatedLeptons(particle_.associatedLeptons),
associatedNeutrinos(particle_.associatedNeutrinos),
associatedPhotons(particle_.associatedPhotons),
associatedJets(particle_.associatedJets),
associatedTops(particle_.associatedTops),
sortedDaughters(particle_.sortedDaughters)
{
  // Recreate the intermediate Vs since mother info needs to be reset
  createSortedVs();
  addAssociatedVs();
}
MELACandidate& MELACandidate::operator=(const MELACandidate& particle_){
  MELACandidate tmp(particle_);
  swap(tmp);
  return *this;
}
MELACandidate::~MELACandidate(){
  resetVs();
  sortedDaughters.clear();
  associatedTops.clear(); 
  associatedJets.clear();
  associatedLeptons.clear();
  associatedNeutrinos.clear();
  associatedPhotons.clear();
}

MELACandidate* MELACandidate::shallowCopy(){
  MELACandidate* cand = new MELACandidate(id, p4, associatedByHighestPt);

  // Copy particle content
  cand->setSelected(passSelection);
  cand->setGenStatus(genStatus);
  cand->setLifetime(lifetime);

  TUtilHelpers::copyVector(mothers, cand->mothers);
  TUtilHelpers::copyVector(daughters, cand->daughters);

  // Copy candidate content
  cand->setShallowCopy(true);
  TUtilHelpers::copyVector(associatedLeptons, cand->associatedLeptons);
  TUtilHelpers::copyVector(associatedNeutrinos, cand->associatedNeutrinos);
  TUtilHelpers::copyVector(associatedPhotons, cand->associatedPhotons);
  TUtilHelpers::copyVector(associatedJets, cand->associatedJets);
  TUtilHelpers::copyVector(associatedTops, cand->associatedTops);
  TUtilHelpers::copyVector(sortedDaughters, cand->sortedDaughters);
  TUtilHelpers::copyVector(sortedVs, cand->sortedVs);

  return cand;
}
void MELACandidate::swap(MELACandidate& particle_){
  MELAParticle::swap(particle_);
  std::swap(associatedByHighestPt, particle_.associatedByHighestPt);
  std::swap(isShallowCopy, particle_.isShallowCopy);
  std::swap(selfDecayMode, particle_.selfDecayMode);
  std::swap(associatedLeptons, particle_.associatedLeptons);
  std::swap(associatedNeutrinos, particle_.associatedNeutrinos);
  std::swap(associatedPhotons, particle_.associatedPhotons);
  std::swap(associatedJets, particle_.associatedJets);
  std::swap(associatedTops, particle_.associatedTops);
  std::swap(sortedDaughters, particle_.sortedDaughters);
  std::swap(sortedVs, particle_.sortedVs);
}

void MELACandidate::resetVs(){
  if (!isShallowCopy){ // Delete owned objects if not a shallow copy
    for (MELAParticle*& aV:sortedVs) delete aV;
  }
  sortedVs.clear();
  isShallowCopy=false; // Should delete sortedVs now
}
void MELACandidate::recreateVs(){
  resetVs();
  createSortedVs();
  addAssociatedVs();
}


void MELACandidate::setDecayMode(TVar::CandidateDecayMode flag){ selfDecayMode=flag; }
void MELACandidate::setAddAssociatedByHighestPt(bool associatedByHighestPt_){ associatedByHighestPt=associatedByHighestPt_; }
void MELACandidate::setShallowCopy(bool flag){ isShallowCopy=flag; }
bool MELACandidate::testShallowCopy(){ return isShallowCopy; }

void MELACandidate::sortDaughters(){
  setDecayMode(PDGHelpers::HDecayMode);
  if (debugVars::debugFlag) MELAout << "Starting MELACandidate::sortDaughters with self decay mode " << selfDecayMode << std::endl;
  if (debugVars::debugFlag) MELAout << "Starting MELACandidate::sortDaughtersInitial" << std::endl;
  sortDaughtersInitial();
  if (debugVars::debugFlag) MELAout << "Starting MELACandidate::sortDaughtersByBestZ1" << std::endl;
  sortDaughtersByBestZ1();
  if (debugVars::debugFlag) MELAout << "Starting MELACandidate::createSortedVs" << std::endl;
  createSortedVs();
}

std::vector<int> MELACandidate::getDaughterIds()const{
  std::vector<int> result;
  for (MELAParticle const* const& dau:sortedDaughters){
    if (dau) result.push_back(dau->id);
  }
  return result;
}
std::vector<int> MELACandidate::getAssociatedParticleIds()const{
  std::vector<int> result;
  std::vector<std::vector<MELAParticle*> const*> list; list.reserve(3);
  list.push_back(&associatedLeptons);
  list.push_back(&associatedPhotons);
  list.push_back(&associatedJets);
  for (std::vector<MELAParticle*> const*& ll:list){
    for (MELAParticle* const& part:(*ll)){
      if (part) result.push_back(part->id);
    }
  }
  return result;
}
MELAParticle* MELACandidate::getSortedDaughter(int index) const{
  if ((int)sortedDaughters.size()>index) return sortedDaughters.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getSortedV(int index) const{
  if ((int)sortedVs.size()>index) return sortedVs.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getAssociatedLepton(int index)const{
  if ((int)associatedLeptons.size()>index) return associatedLeptons.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getAssociatedNeutrino(int index)const{
  if ((int)associatedNeutrinos.size()>index) return associatedNeutrinos.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getAssociatedPhoton(int index)const{
  if ((int)associatedPhotons.size()>index) return associatedPhotons.at(index);
  else return 0;
}
MELAParticle* MELACandidate::getAssociatedJet(int index)const{
  if ((int)associatedJets.size()>index) return associatedJets.at(index);
  else return 0;
}
MELATopCandidate* MELACandidate::getAssociatedTop(int index)const{
  if ((int)associatedTops.size()>index) return associatedTops.at(index);
  else return 0;
}

std::vector<MELAParticle*>& MELACandidate::getSortedDaughters(){ return sortedDaughters; }
std::vector<MELAParticle*>& MELACandidate::getSortedVs(){ return sortedVs; }
std::vector<MELAParticle*>& MELACandidate::getAssociatedLeptons(){ return associatedLeptons; }
std::vector<MELAParticle*>& MELACandidate::getAssociatedNeutrinos(){ return associatedNeutrinos; }
std::vector<MELAParticle*>& MELACandidate::getAssociatedPhotons(){ return associatedPhotons; }
std::vector<MELAParticle*>& MELACandidate::getAssociatedJets(){ return associatedJets; }
std::vector<MELATopCandidate*>& MELACandidate::getAssociatedTops(){ return associatedTops; }

const std::vector<MELAParticle*>& MELACandidate::getSortedDaughters()const{ return sortedDaughters; }
const std::vector<MELAParticle*>& MELACandidate::getSortedVs()const{ return sortedVs; }
const std::vector<MELAParticle*>& MELACandidate::getAssociatedLeptons()const{ return associatedLeptons; }
const std::vector<MELAParticle*>& MELACandidate::getAssociatedNeutrinos()const{ return associatedNeutrinos; }
const std::vector<MELAParticle*>& MELACandidate::getAssociatedPhotons()const{ return associatedPhotons; }
const std::vector<MELAParticle*>& MELACandidate::getAssociatedJets()const{ return associatedJets; }
const std::vector<MELATopCandidate*>& MELACandidate::getAssociatedTops()const{ return associatedTops; }

std::vector<MELAParticle*> MELACandidate::getAssociatedSortedVs(){
  std::vector<MELAParticle*> res;
  std::vector<MELAParticle*>::iterator itEnd=sortedVs.end();
  std::vector<MELAParticle*>::iterator itBegin=itEnd;
  for (std::vector<MELAParticle*>::iterator it=sortedVs.begin(); it!=itEnd; it++){
    bool doSkip=false;
    for (auto const& dau:sortedDaughters){ if ((*it)->hasDaughter(dau)){ doSkip=true; break; } }
    if (!doSkip){ itBegin=it; break; }
  }
  if (itBegin!=itEnd) std::copy(itBegin, itEnd, std::back_inserter(res));
  return res;
}
std::vector<MELAParticle*> MELACandidate::getAssociatedSortedVs()const{
  std::vector<MELAParticle*> res;
  std::vector<MELAParticle*>::const_iterator itEnd=sortedVs.cend();
  std::vector<MELAParticle*>::const_iterator itBegin=itEnd;
  for (std::vector<MELAParticle*>::const_iterator it=sortedVs.cbegin(); it!=itEnd; it++){
    bool doSkip=false;
    for (auto const& dau:sortedDaughters){ if ((*it)->hasDaughter(dau)){ doSkip=true; break; } }
    if (!doSkip){ itBegin=it; break; }
  }
  if (itBegin!=itEnd) std::copy(itBegin, itEnd, std::back_inserter(res));
  return res;
}

void MELACandidate::sortDaughtersInitial(){
  bool beginWithIdPair = (
    selfDecayMode==TVar::CandidateDecay_ZZ
    || selfDecayMode==TVar::CandidateDecay_ZW
    || selfDecayMode==TVar::CandidateDecay_ZG
    || selfDecayMode==TVar::CandidateDecay_GG
    || selfDecayMode==TVar::CandidateDecay_ff
    );
  bool beginWithWPair = (selfDecayMode==TVar::CandidateDecay_WW || selfDecayMode==TVar::CandidateDecay_WG);
  int nDaughtersBooked=0;
  int tmpDindex[2]={ 0 };
  MELAParticle* df[2] ={ getDaughter(0), 0 };
  if (df[0]!=0) nDaughtersBooked++;
  for (int j=1; j<getNDaughters(); j++){
    MELAParticle* dtmp = getDaughter(j);
    if (
      (beginWithIdPair && std::abs(dtmp->charge())-std::abs(df[0]->charge())==0 && std::abs(dtmp->id)==std::abs(df[0]->id)) // First daughter in ZZ/ZW/ZG/GG/ff requires identical |Q| and |id|.
      ||
      ( // First daughter in WW/WG requires a pair with |sum(Q)|=1 and opposite id signs, or two unknown jets.
      beginWithWPair
      &&
      (
      (std::abs(dtmp->charge()+df[0]->charge())==1 && TMath::Sign(1, dtmp->id)==-TMath::Sign(1, df[0]->id))
      ||
      (PDGHelpers::isAnUnknownJet(dtmp->id) && PDGHelpers::isAnUnknownJet(df[0]->id))
      )
      )
      ){
      df[1] = dtmp;
      tmpDindex[1] = j;
      nDaughtersBooked++;
      break;
    }
  }
  MELAParticle* ds[2] ={ 0 };
  int sindex=0;
  for (int j=1; j<getNDaughters(); j++){
    if (j==tmpDindex[1]) continue;
    MELAParticle* dtmp = getDaughter(j);
    ds[sindex] = dtmp;
    nDaughtersBooked++;
    sindex++;
    if (sindex==2) break;
  }

  if (nDaughtersBooked!=getNDaughters()){
    if (getNDaughters()>4) MELAout << "MELACandidate::sortDaughtersInitial: Number of daughters passed " << getNDaughters() << ">4 is currently not supported." << std::endl;
    MELAout << "MELACandidate::sortDaughtersInitial: Number of daughters passed (" << getNDaughters() << ") is not the same as number of daughters booked (" << nDaughtersBooked << ")! Aborting, no daughters can be recorded." << std::endl;
    for (int idau=0; idau<getNDaughters(); idau++){
      MELAParticle* part = getDaughter(idau);
      if (part!=0) MELAout << "\t- Passed daughter " << idau << " (X, Y, Z, T, M, Id) = "
        << part->x() << " " << part->y() << " " << part->z() << " " << part->t() << " " << part->m() << " " << part->id << std::endl;
      else MELAout << "\t- Passed daughter " << idau << " = 0!" << std::endl;
    }
    for (unsigned int j=0; j<2; j++){ if (df[j]!=0) MELAout << "\t- df[" << j << "] (X, Y, Z, T, M, Id) = " << df[j]->x() << " " << df[j]->y() << " " << df[j]->z() << " " << df[j]->t() << " " << df[j]->m() << " " << df[j]->id << std::endl; }
    for (unsigned int j=0; j<2; j++){ if (ds[j]!=0) MELAout << "\t- ds[" << j << "] (X, Y, Z, T, M, Id) = " << ds[j]->x() << " " << ds[j]->y() << " " << ds[j]->z() << " " << ds[j]->t() << " " << ds[j]->m() << " " << ds[j]->id << std::endl; }
    return;
  }

  if (
    (df[0]!=0 && df[1]!=0)
    &&
    (
    // Order by ubar(0)v(1)
    (df[0]->id<df[1]->id && beginWithIdPair)
    ||
    (df[0]->id<df[1]->id && df[0]->id<0 && beginWithWPair)
    ||
    ((df[0]->id*df[1]->id>0 || (df[0]->id==0 && df[1]->id==0)) && df[0]->phi()<df[1]->phi())
    )
    ){
    MELAParticle* dtmp = df[0];
    df[0] = df[1];
    df[1] = dtmp;
  }

  if (
    (ds[0]!=0 && ds[1]!=0)
    &&
    (
    // Order by ubar(0)v(1)
    (ds[0]->id<ds[1]->id && beginWithIdPair)
    ||
    (ds[0]->id<ds[1]->id && ds[0]->id<0 && beginWithWPair)
    ||
    ((ds[0]->id*ds[1]->id>0 || (ds[0]->id==0 && ds[1]->id==0)) && ds[0]->phi()<ds[1]->phi())
    )
    ){
    MELAParticle* dtmp = ds[0];
    ds[0] = ds[1];
    ds[1] = dtmp;
  }
  if (df[1]==0 && df[0]!=0 && ds[0]!=0 && ds[1]!=0){ // Swap GZ to ZG
    for (int ip=0; ip<2; ip++){
      MELAParticle* dtmp = ds[ip];
      ds[ip] = df[ip];
      df[ip] = dtmp;
    }
  }
  for (int i=0; i<2; i++){
    if (df[i]!=0) sortedDaughters.push_back(df[i]);
  }
  for (int i=0; i<2; i++){
    if (ds[i]!=0) sortedDaughters.push_back(ds[i]);
  }
}
void MELACandidate::sortDaughtersByBestZ1(){
  bool beginWithZPair = (
    selfDecayMode==TVar::CandidateDecay_ZZ
    || selfDecayMode==TVar::CandidateDecay_ZW
    || selfDecayMode==TVar::CandidateDecay_ZG
    );
  bool beginWithMasslessPair = (
    selfDecayMode==TVar::CandidateDecay_GG
    || selfDecayMode==TVar::CandidateDecay_ff
    );
  bool beginWithIdPair = beginWithZPair || beginWithMasslessPair;
  bool beginWithWPair = (selfDecayMode==TVar::CandidateDecay_WW || selfDecayMode==TVar::CandidateDecay_WG);

  double HVVmass = PDGHelpers::Zeromass;
  if (beginWithZPair) HVVmass = PDGHelpers::Zmass;
  else if (beginWithWPair) HVVmass = PDGHelpers::Wmass;

  MELAParticle* orderedDs[2][2]={ { 0 } };
  TLorentzVector pZ1(0, 0, 0, 0);
  TLorentzVector pZ2(0, 0, 0, 0);
  if (sortedDaughters.size()>2){ // WW, ZZ, ZW, WG, ZG
    bool dauDiffType = true;
    if (debugVars::debugFlag) MELAout << "Ndaughters>2" << std::endl;

    for (int d=0; d<std::min(2, (int)sortedDaughters.size()); d++){
      if (sortedDaughters.at(d)!=0) pZ1 = pZ1 + sortedDaughters.at(d)->p4;
    }
    for (int d=std::min(2, (int)sortedDaughters.size()); d<std::min(4, (int)sortedDaughters.size()); d++){
      if (sortedDaughters.at(d)!=0) pZ2 = pZ2 + sortedDaughters.at(d)->p4;
    }

    if (debugVars::debugFlag) MELAout << "Preliminary pZ1 and pZ2 calculated!" << std::endl;

    if (sortedDaughters.size()>=4){
      if (
        beginWithIdPair &&
        (
        (isALepton(sortedDaughters.at(0)->id) && isALepton(sortedDaughters.at(1)->id) && isALepton(sortedDaughters.at(2)->id) && isALepton(sortedDaughters.at(3)->id))
        ||
        (isANeutrino(sortedDaughters.at(0)->id) && isANeutrino(sortedDaughters.at(1)->id) && isANeutrino(sortedDaughters.at(2)->id) && isANeutrino(sortedDaughters.at(3)->id))
        ||
        (isAPhoton(sortedDaughters.at(0)->id) && isAPhoton(sortedDaughters.at(1)->id) && isAPhoton(sortedDaughters.at(2)->id) && isAPhoton(sortedDaughters.at(3)->id))
        ||
        (isDownTypeQuark(sortedDaughters.at(0)->id) && isDownTypeQuark(sortedDaughters.at(1)->id) && isDownTypeQuark(sortedDaughters.at(2)->id) && isDownTypeQuark(sortedDaughters.at(3)->id))
        ||
        (isUpTypeQuark(sortedDaughters.at(0)->id) && isUpTypeQuark(sortedDaughters.at(1)->id) && isUpTypeQuark(sortedDaughters.at(2)->id) && isUpTypeQuark(sortedDaughters.at(3)->id))
        ||
        (isAGluon(sortedDaughters.at(0)->id) && isAGluon(sortedDaughters.at(1)->id) && isAGluon(sortedDaughters.at(2)->id) && isAGluon(sortedDaughters.at(3)->id))
        ||
        (isAnUnknownJet(sortedDaughters.at(0)->id) && isAnUnknownJet(sortedDaughters.at(1)->id) && isAnUnknownJet(sortedDaughters.at(2)->id) && isAnUnknownJet(sortedDaughters.at(3)->id))
        )
        ) dauDiffType=false;
    }

    if (
      (dauDiffType && beginWithIdPair && (
      sortedDaughters.size()<4 // WG, ZG
      ||
      (
      // WW, ZZ, ZW (mostly relevant for Z1-Z2 sorting)
      // Zll>Znn>?GG>Zdd>Zuu>?gg>Zjj. 0-2 swap may then give W+W-.
      sortedDaughters.size()>=4 && (
      isALepton(sortedDaughters.at(0)->id) ||
      (isANeutrino(sortedDaughters.at(0)->id) && !isALepton(sortedDaughters.at(2)->id)) ||
      (isAPhoton(sortedDaughters.at(0)->id) && !isALepton(sortedDaughters.at(2)->id) && !isANeutrino(sortedDaughters.at(2)->id)) ||
      (isDownTypeQuark(sortedDaughters.at(0)->id) && !isALepton(sortedDaughters.at(2)->id) && !isANeutrino(sortedDaughters.at(2)->id) && !isAPhoton(sortedDaughters.at(2)->id)) ||
      (isUpTypeQuark(sortedDaughters.at(0)->id) && !isALepton(sortedDaughters.at(2)->id) && !isANeutrino(sortedDaughters.at(2)->id) && !isAPhoton(sortedDaughters.at(2)->id) && !isDownTypeQuark(sortedDaughters.at(2)->id)) ||
      (isAGluon(sortedDaughters.at(0)->id) && !isALepton(sortedDaughters.at(2)->id) && !isANeutrino(sortedDaughters.at(2)->id) && !isAPhoton(sortedDaughters.at(2)->id) && !isDownTypeQuark(sortedDaughters.at(2)->id) && !isUpTypeQuark(sortedDaughters.at(2)->id)) ||
      (isAnUnknownJet(sortedDaughters.at(0)->id) && !isALepton(sortedDaughters.at(2)->id) && !isANeutrino(sortedDaughters.at(2)->id) && !isAPhoton(sortedDaughters.at(2)->id) && !isDownTypeQuark(sortedDaughters.at(2)->id) && !isUpTypeQuark(sortedDaughters.at(2)->id) && !isAGluon(sortedDaughters.at(2)->id))
      )
      )
      )
      )
      ||
      // ZZ->4f
      (
      !dauDiffType && beginWithIdPair && std::abs(pZ1.M() - HVVmass)<std::abs(pZ2.M() - HVVmass) // Z1 / Z2
      )
      ||
      (sortedDaughters.at(0)!=0 && sortedDaughters.at(1)!=0 && beginWithWPair && (sortedDaughters.at(0)->charge()+sortedDaughters.at(1)->charge())>0) // W+ / W-
      ){
      if (debugVars::debugFlag) MELAout << "pZ1 is closer to HVVmass " << HVVmass << std::endl;
      orderedDs[0][0]=sortedDaughters.at(0);
      orderedDs[0][1]=sortedDaughters.at(1);
      orderedDs[1][0]=((int)sortedDaughters.size()>2 ? sortedDaughters.at(2) : 0);
      orderedDs[1][1]=((int)sortedDaughters.size()>3 ? sortedDaughters.at(3) : 0);
    }
    else{
      if (debugVars::debugFlag) MELAout << "pZ2 is closer to HVVmass " << HVVmass << std::endl;
      orderedDs[0][0]=((int)sortedDaughters.size()>2 ? sortedDaughters.at(2) : 0);
      orderedDs[0][1]=((int)sortedDaughters.size()>3 ? sortedDaughters.at(3) : 0);
      orderedDs[1][0]=sortedDaughters.at(0);
      orderedDs[1][1]=sortedDaughters.at(1);
      TLorentzVector ptmp = pZ1;
      pZ1 = pZ2;
      pZ2 = ptmp;
    }
  }
  else if (sortedDaughters.size()==2){ // GG, ffbar
    if (debugVars::debugFlag) MELAout << "Ndaughters==2" << std::endl;

    if (sortedDaughters.at(0)!=0) pZ1 = pZ1 + sortedDaughters.at(0)->p4;
    if (sortedDaughters.at(1)!=0) pZ2 = pZ2 + sortedDaughters.at(1)->p4;
    orderedDs[0][0]=sortedDaughters.at(0);
    orderedDs[0][1]=0;
    orderedDs[1][0]=sortedDaughters.at(1);
    orderedDs[1][1]=0;
  }
  if (
    (
    (orderedDs[0][1]!=0 && orderedDs[1][1]!=0)
    &&
    (orderedDs[0][1]->id == orderedDs[1][1]->id)
    )
    ||
    (
    (orderedDs[1][0]!=0 && orderedDs[0][0]!=0)
    &&
    (orderedDs[1][0]->id == orderedDs[0][0]->id)
    )
    ){
    if (debugVars::debugFlag) MELAout << "Checking alternative pairings." << std::endl;

    MELAParticle* orderedDps[2][2]={ { 0 } };

    TLorentzVector pZ1p(0, 0, 0, 0);
    TLorentzVector pZ2p(0, 0, 0, 0);
    for (int d=0; d<2; d++){
      if (orderedDs[d][d]!=0) pZ1p = pZ1p + orderedDs[d][d]->p4;
    }
    for (int d=0; d<2; d++){
      if (orderedDs[1-d][d]!=0) pZ2p = pZ2p + orderedDs[1-d][d]->p4;
    }

    if (std::abs(pZ1p.M() - HVVmass)<std::abs(pZ2p.M() - HVVmass)){
      orderedDps[0][0]=orderedDs[0][0];
      orderedDps[0][1]=orderedDs[1][1];
      orderedDps[1][0]=orderedDs[1][0];
      orderedDps[1][1]=orderedDs[0][1];
    }
    else{
      orderedDps[0][0]=orderedDs[1][0];
      orderedDps[0][1]=orderedDs[0][1];
      orderedDps[1][0]=orderedDs[0][0];
      orderedDps[1][1]=orderedDs[1][1];
      TLorentzVector ptmp = pZ1p;
      pZ1p = pZ2p;
      pZ2p = ptmp;
    }
    if (std::abs(pZ1p.M() - HVVmass)<std::abs(pZ1.M() - HVVmass) || (std::abs(pZ1p.M() - HVVmass)==std::abs(pZ1.M() - HVVmass) && pZ2p.Pt()>pZ2.Pt())){
      for (int i=0; i<2; i++){
        for (int j=0; j<2; j++) orderedDs[i][j] = orderedDps[i][j];
      }
    }
  }
  sortedDaughters.clear();
  for (int i=0; i<2; i++){
    for (int j=0; j<2; j++){
      if (orderedDs[i][j]!=0) sortedDaughters.push_back(orderedDs[i][j]);
    }
  }
  if (debugVars::debugFlag) MELAout << "Final number of daughters in sortedDaughters: " << sortedDaughters.size() << std::endl;
}
void MELACandidate::createSortedVs(){
  TLorentzVector nullVector(0, 0, 0, 0);
  bool beginWithZPair = (
    selfDecayMode==TVar::CandidateDecay_ZZ
    || selfDecayMode==TVar::CandidateDecay_ZW
    || selfDecayMode==TVar::CandidateDecay_ZG
    );
  bool beginWithWPair = (selfDecayMode==TVar::CandidateDecay_WW || selfDecayMode==TVar::CandidateDecay_WG);
  int VID = 21; // Gamma/Gamma*
  if (beginWithWPair) VID = 24; // W+/W-
  else if (beginWithZPair) VID = 23; // Z

  int icutoff = (sortedDaughters.size()>2 ? 2 : std::min(1, (int)sortedDaughters.size()));
  int imax = std::min(4, (int)sortedDaughters.size());
  int V1id=-9000, V2id=-9000;
  if (icutoff==0) V1id=25;
  else{
    double Vcharge[2]={ 0 };
    for (int d=0; d<icutoff; d++){
      if (sortedDaughters.at(d)!=0){
        if (icutoff==2){
          V1id=VID;
          Vcharge[0] += sortedDaughters.at(d)->charge();
        }
        else if (icutoff==1) V1id=sortedDaughters.at(d)->id;
      }
    }
    for (int d=icutoff; d<imax; d++){
      if (sortedDaughters.at(d)!=0){
        if ((imax-icutoff)==2){
          V2id=(VID!=24 ? VID : -VID);
          Vcharge[1] += sortedDaughters.at(d)->charge();
        }
        else if ((imax-icutoff)==1) V2id=sortedDaughters.at(d)->id;
      }
    }
    // Override selfDecayMode if charges indicate some other final state
    if (fabs(Vcharge[0]-1.)<0.001) V1id=24;
    else if (fabs(Vcharge[0]+1.)<0.001) V1id=-24;
    if (fabs(Vcharge[1]-1.)<0.001) V2id=24;
    else if (fabs(Vcharge[1]+1.)<0.001) V2id=-24;
  }

  // If the number of Zs is less than 2, should still create empty particles
  MELAParticle* Z1 = new MELAParticle(V1id);
  MELAParticle* Z2 = new MELAParticle(V2id);

  // Add daughters and their 4-vectors
  if (icutoff==0) Z1->p4 = this->p4; // V1=H, V2=invalid
  else{
    for (int d=0; d<icutoff; d++){ if (sortedDaughters.at(d)!=0) *Z1 += sortedDaughters.at(d); }
    for (int d=icutoff; d<imax; d++){ if (sortedDaughters.at(d)!=0) *Z2 += sortedDaughters.at(d); }
  }

  // Configure relationships
  Z1->addMother(this);
  addSortedV(Z1);
  Z2->addMother(this);
  addSortedV(Z2);
}
TLorentzVector MELACandidate::getAlternativeVMomentum(int index)const{
  bool beginWithZPair = (
    selfDecayMode==TVar::CandidateDecay_ZZ
    || selfDecayMode==TVar::CandidateDecay_ZW
    || selfDecayMode==TVar::CandidateDecay_ZG
    );
  bool beginWithWPair = (selfDecayMode==TVar::CandidateDecay_WW || selfDecayMode==TVar::CandidateDecay_WG);

  double HVVmass = PDGHelpers::Zeromass;
  if (beginWithZPair) HVVmass = PDGHelpers::Zmass;
  else if (beginWithWPair) HVVmass = PDGHelpers::Wmass;

  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (sortedDaughters.size()>=3){
    TLorentzVector pZ1 = sortedDaughters.at(0)->p4;
    if (sortedDaughters.size()>3) pZ1 = pZ1 + sortedDaughters.at(3)->p4;
    TLorentzVector pZ2 = sortedDaughters.at(2)->p4+sortedDaughters.at(1)->p4;
    if (std::abs(pZ1.M() - HVVmass)>std::abs(pZ2.M() - HVVmass)){
      TLorentzVector pZtmp = pZ1;
      pZ1 = pZ2;
      pZ2 = pZtmp;
    }
    return (index==0 ? pZ1 : pZ2);
  }
  else return nullFourVector;
}
bool MELACandidate::daughtersInterfere()const{
  bool doInterfere = false;
  if (sortedDaughters.size()>3) doInterfere = (
    (sortedDaughters.at(0))->id == (sortedDaughters.at(2))->id
    &&
    (sortedDaughters.at(1))->id == (sortedDaughters.at(3))->id
    &&
    !PDGHelpers::isAnUnknownJet(sortedDaughters.at(0)->id) && !PDGHelpers::isAnUnknownJet(sortedDaughters.at(2)->id)
    );
  return doInterfere;
}

bool MELACandidate::checkDaughtership(MELAParticle const* myParticle)const{
  return hasDaughter(myParticle);
}

void MELACandidate::addAssociatedLeptons(MELAParticle* myParticle){
  addAssociatedParticleToArray(myParticle, associatedLeptons);
}
void MELACandidate::addAssociatedNeutrinos(MELAParticle* myParticle){
  addAssociatedParticleToArray(myParticle, associatedLeptons); // Neutrinos are leptons at the ZZ candidate level
  addAssociatedParticleToArray(myParticle, associatedNeutrinos);
}
void MELACandidate::addAssociatedPhotons(MELAParticle* myParticle){
  addAssociatedParticleToArray(myParticle, associatedPhotons);
}
void MELACandidate::addAssociatedJets(MELAParticle* myParticle){
  addAssociatedParticleToArray(myParticle, associatedJets);
}
void MELACandidate::addAssociatedTops(MELATopCandidate* myParticle){
  addAssociatedParticleToArray(myParticle, associatedTops);
}
void MELACandidate::addAssociatedParticleToArray(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray){
  if (checkDaughtership(myParticle)) return;
  if (associatedByHighestPt) MELACandidate::addByHighestPt(myParticle, particleArray);
  else MELACandidate::addUnordered(myParticle, particleArray);
}
void MELACandidate::addAssociatedParticleToArray(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray){
  if (checkDaughtership(myParticle)) return;
  if (associatedByHighestPt) MELACandidate::addByHighestPt(myParticle, particleArray);
  else MELACandidate::addUnordered(myParticle, particleArray);
}

void MELACandidate::addAssociatedVs(){
  createAssociatedVs(associatedJets);
  createAssociatedVs(associatedLeptons);
}
void MELACandidate::createAssociatedVs(std::vector<MELAParticle*>& particleArray){
  for (unsigned int i = 0; i<particleArray.size(); i++){
    double Qi = particleArray.at(i)->charge();
    const int& id_i = particleArray.at(i)->id;

    for (unsigned int j = i+1; j<particleArray.size(); j++){
      double Qj = particleArray.at(j)->charge();
      const int& id_j = particleArray.at(j)->id;

      int bosonId=-1;
      if ((Qi+Qj)==0 && (id_i+id_j)==0) bosonId = (id_i==0 ? 0 : 23); // Z->(f,fb) or 0->(0,0)
      else if (PDGHelpers::isAnUnknownJet(id_i) || PDGHelpers::isAnUnknownJet(id_j)) bosonId = 0; // 0->(q,0)/(0,qb)
      else if ( // W->f fb'
        abs(Qi+Qj)==1
        &&
        (
        ((PDGHelpers::isALepton(id_i) || PDGHelpers::isALepton(id_j)) && std::abs(id_i+id_j)==1) // Require SF lnu particle-antiparticle pairs
        ||
        (PDGHelpers::isUpTypeQuark(id_i) && PDGHelpers::isDownTypeQuark(id_j)) // Require ud- or du-type pairs, qqbar requirement is satisfied with charge.
        ||
        (PDGHelpers::isDownTypeQuark(id_i) && PDGHelpers::isUpTypeQuark(id_j))
        )
        ) bosonId=24*(Qi+Qj);

      if (bosonId!=-1){
        MELAParticle* boson = new MELAParticle(bosonId);
        // Order by f-f'b
        int firstdaughter = i, seconddaughter = j;
        if (
          (particleArray.at(firstdaughter)->id<particleArray.at(seconddaughter)->id)
          &&
          (
          !PDGHelpers::isAWBoson(bosonId)
          ||
          (particleArray.at(firstdaughter)->id<0 && PDGHelpers::isAWBoson(bosonId))
          )
          ) std::swap(firstdaughter, seconddaughter);
        *boson += particleArray.at(firstdaughter);
        *boson += particleArray.at(seconddaughter);

        boson->addMother(this);
        addSortedV(boson);
      }
    }
  }
}

void MELACandidate::testPreSelectedDaughters(){
  for (auto& dau:getDaughters()){
    if (!dau->passSelection){
      passSelection=false;
      break;
    }
  }
}

void MELACandidate::getRelatedParticles(std::vector<MELAParticle*>& particles){
  MELAParticle::getRelatedParticles(particles);
  for (auto& part:sortedDaughters) part->getRelatedParticles(particles); // Hopefully no particle gets added from here
  for (auto& part:associatedTops) part->getRelatedParticles(particles);
  for (auto& part:sortedVs) part->getRelatedParticles(particles);
  for (auto& part:associatedLeptons) part->getRelatedParticles(particles);
  for (auto& part:associatedPhotons) part->getRelatedParticles(particles);
  for (auto& part:associatedJets) part->getRelatedParticles(particles);
}

void MELACandidate::addUnordered(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray){
  bool inserted = MELAParticle::checkParticleExists(myParticle, particleArray); // Test if the particle is already in the vector
  if (!inserted) particleArray.push_back(myParticle);
}
void MELACandidate::addUnordered(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray){
  bool inserted = MELATopCandidate::checkTopCandidateExists(myParticle, particleArray); // Test if the particle is already in the vector
  if (!inserted) particleArray.push_back(myParticle);
}
void MELACandidate::addByHighestPt(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray){
  bool inserted = MELAParticle::checkParticleExists(myParticle, particleArray); // Test if the particle is already in the vector
  if (!inserted){
    for (std::vector<MELAParticle*>::iterator it = particleArray.begin(); it<particleArray.end(); it++){
      if ((*it)->pt()<myParticle->pt()){
        inserted=true;
        particleArray.insert(it, myParticle);
        break;
      }
    }
    if (!inserted) particleArray.push_back(myParticle);
  }
}
void MELACandidate::addByHighestPt(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray){
  bool inserted = MELATopCandidate::checkTopCandidateExists(myParticle, particleArray); // Test if the particle is already in the vector
  if (!inserted){
    for (std::vector<MELATopCandidate*>::iterator it = particleArray.begin(); it<particleArray.end(); it++){
      if ((*it)->pt()<myParticle->pt()){
        inserted=true;
        particleArray.insert(it, myParticle);
        break;
      }
    }
    if (!inserted) particleArray.push_back(myParticle);
  }
}
