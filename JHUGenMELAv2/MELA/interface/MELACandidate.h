#ifndef MELACANDIDATE_H
#define MELACANDIDATE_H

#include "MELAParticle.h"
#include "MELATopCandidate.h"

class MELACandidate : public MELAParticle{
public:
  MELACandidate(int id_, TLorentzVector p4_, bool associatedByHighestPt_=false);
  ~MELACandidate();
  MELACandidate* shallowCopy();

  // Member functions

  MELAParticle* getSortedDaughter(int index)const;
  MELAParticle* getSortedV(int index)const;
  MELAParticle* getAssociatedLepton(int index)const;
  MELAParticle* getAssociatedNeutrino(int index)const;
  MELAParticle* getAssociatedPhoton(int index)const;
  MELAParticle* getAssociatedJet(int index)const;
  MELATopCandidate* getAssociatedTop(int index)const;
  TLorentzVector getAlternativeVMomentum(int index)const;

  virtual std::vector<int> getDaughterIds()const;
  std::vector<int> getAssociatedParticleIds()const;

  TVar::CandidateDecayMode getDecayMode()const{ return selfDecayMode; }

  int getNAssociatedLeptons()const{ return associatedLeptons.size(); }
  int getNAssociatedNeutrinos()const{ return associatedNeutrinos.size(); }
  int getNAssociatedPhotons()const{ return associatedPhotons.size(); }
  int getNAssociatedJets()const{ return associatedJets.size(); }
  int getNAssociatedTops()const{ return associatedTops.size(); }
  int getNSortedVs()const{ return sortedVs.size(); }

  void addAssociatedLeptons(MELAParticle* myParticle);
  void addAssociatedNeutrinos(MELAParticle* myParticle);
  void addAssociatedPhotons(MELAParticle* myParticle);
  void addAssociatedJets(MELAParticle* myParticle);
  void addAssociatedTops(MELATopCandidate* myParticle);

  void addSortedV(MELAParticle* myParticle){ sortedVs.push_back(myParticle); }
  void addAssociatedVs();

  void sortDaughters();
  void testPreSelectedDaughters();
  bool testShallowCopy();

  bool daughtersInterfere()const;
  void setDecayMode(TVar::CandidateDecayMode flag);
  void setAddAssociatedByHighestPt(bool associatedByHighestPt_);
  void setShallowCopy(bool flag);

protected:
  bool associatedByHighestPt;
  bool isShallowCopy;

  TVar::CandidateDecayMode selfDecayMode;

  std::vector<MELAParticle*> associatedLeptons;
  std::vector<MELAParticle*> associatedNeutrinos;
  std::vector<MELAParticle*> associatedPhotons;
  std::vector<MELAParticle*> associatedJets;
  std::vector<MELATopCandidate*> associatedTops;

  std::vector<MELAParticle*> sortedDaughters;
  std::vector<MELAParticle*> sortedVs;

  void sortDaughtersInitial();
  void sortDaughtersByBestZ1();
  void createSortedVs();
  bool checkDaughtership(MELAParticle* myParticle)const;
  void createAssociatedVs(std::vector<MELAParticle*>& particleArray);
  void addByHighestPt(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray);
  void addByHighestPt(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray);

  bool checkTopCandidateExists(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray)const;
};



#endif
