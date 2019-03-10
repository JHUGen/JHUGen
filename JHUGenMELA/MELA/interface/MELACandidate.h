#ifndef MELACANDIDATE_H
#define MELACANDIDATE_H

#include "MELAParticle.h"
#include "MELATopCandidate.h"

class MELACandidate : public MELAParticle{
public:
  MELACandidate();
  MELACandidate(int id_, bool associatedByHighestPt_=false);
  MELACandidate(int id_, TLorentzVector p4_, bool associatedByHighestPt_=false);
  MELACandidate(const MELACandidate& particle_);
  MELACandidate& operator=(const MELACandidate& particle_);
  ~MELACandidate();
  MELACandidate* shallowCopy();
  void swap(MELACandidate& particle_);

  // Member functions

  MELAParticle* getSortedDaughter(int index)const;
  MELAParticle* getSortedV(int index)const;
  MELAParticle* getAssociatedLepton(int index)const;
  MELAParticle* getAssociatedNeutrino(int index)const;
  MELAParticle* getAssociatedPhoton(int index)const;
  MELAParticle* getAssociatedJet(int index)const;
  MELATopCandidate* getAssociatedTop(int index)const;

  std::vector<MELAParticle*>& getSortedDaughters();
  std::vector<MELAParticle*>& getSortedVs();
  std::vector<MELAParticle*>& getAssociatedLeptons();
  std::vector<MELAParticle*>& getAssociatedNeutrinos();
  std::vector<MELAParticle*>& getAssociatedPhotons();
  std::vector<MELAParticle*>& getAssociatedJets();
  std::vector<MELATopCandidate*>& getAssociatedTops();
  
  const std::vector<MELAParticle*>& getSortedDaughters()const;
  const std::vector<MELAParticle*>& getSortedVs()const;
  const std::vector<MELAParticle*>& getAssociatedLeptons()const;
  const std::vector<MELAParticle*>& getAssociatedNeutrinos()const;
  const std::vector<MELAParticle*>& getAssociatedPhotons()const;
  const std::vector<MELAParticle*>& getAssociatedJets()const;
  const std::vector<MELATopCandidate*>& getAssociatedTops()const;

  std::vector<MELAParticle*> getAssociatedSortedVs();
  std::vector<MELAParticle*> getAssociatedSortedVs()const;

  void getRelatedParticles(std::vector<MELAParticle*>& particles);

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

  void resetVs();
  void recreateVs();

  void sortDaughters();
  void testPreSelectedDaughters();
  bool testShallowCopy();

  bool daughtersInterfere()const;
  void setDecayMode(TVar::CandidateDecayMode flag);
  void setAddAssociatedByHighestPt(bool associatedByHighestPt_);
  void setShallowCopy(bool flag);

  static void addUnordered(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray);
  static void addUnordered(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray);
  static void addByHighestPt(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray);
  static void addByHighestPt(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray);

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
  bool checkDaughtership(MELAParticle const* myParticle)const;
  void createAssociatedVs(std::vector<MELAParticle*>& particleArray);

  void addAssociatedParticleToArray(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray);
  void addAssociatedParticleToArray(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray);

};



#endif
