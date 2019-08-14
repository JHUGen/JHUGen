#ifndef MELACANDIDATE_H
#define MELACANDIDATE_H

#include "MELAParticle.h"
#include "MELAThreeBodyDecayCandidate.h"

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
  MELATopCandidate_t* getAssociatedTop(int index)const;

  std::vector<MELAParticle*>& getSortedDaughters();
  std::vector<MELAParticle*>& getSortedVs();
  std::vector<MELAParticle*>& getAssociatedLeptons();
  std::vector<MELAParticle*>& getAssociatedNeutrinos();
  std::vector<MELAParticle*>& getAssociatedPhotons();
  std::vector<MELAParticle*>& getAssociatedJets();
  std::vector<MELATopCandidate_t*>& getAssociatedTops();
  
  const std::vector<MELAParticle*>& getSortedDaughters()const;
  const std::vector<MELAParticle*>& getSortedVs()const;
  const std::vector<MELAParticle*>& getAssociatedLeptons()const;
  const std::vector<MELAParticle*>& getAssociatedNeutrinos()const;
  const std::vector<MELAParticle*>& getAssociatedPhotons()const;
  const std::vector<MELAParticle*>& getAssociatedJets()const;
  const std::vector<MELATopCandidate_t*>& getAssociatedTops()const;

  std::vector<MELAParticle*> getAssociatedSortedVs();
  std::vector<MELAParticle*> getAssociatedSortedVs()const;

  void getRelatedParticles(std::vector<MELAParticle*>& particles) const;
  void getDaughterParticles(std::vector<MELAParticle*>& particles) const;

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

  void addAssociatedLepton(MELAParticle* myParticle);
  void addAssociatedNeutrino(MELAParticle* myParticle);
  void addAssociatedPhoton(MELAParticle* myParticle);
  void addAssociatedJet(MELAParticle* myParticle);
  void addAssociatedTop(MELATopCandidate_t* myParticle);

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
  static void addUnordered(MELAThreeBodyDecayCandidate* myParticle, std::vector<MELAThreeBodyDecayCandidate*>& particleArray);
  static void addByHighestPt(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray);
  static void addByHighestPt(MELAThreeBodyDecayCandidate* myParticle, std::vector<MELAThreeBodyDecayCandidate*>& particleArray);

protected:
  bool associatedByHighestPt;
  bool isShallowCopy;
  TVar::CandidateDecayMode selfDecayMode;

  std::vector<MELAParticle*> associatedLeptons;
  std::vector<MELAParticle*> associatedNeutrinos;
  std::vector<MELAParticle*> associatedPhotons;
  std::vector<MELAParticle*> associatedJets;
  std::vector<MELATopCandidate_t*> associatedTops;

  std::vector<MELAParticle*> sortedDaughters;
  std::vector<MELAParticle*> sortedVs;

  void sortDaughtersInitial();
  void sortDaughtersByBestZ1();
  void createSortedVs();
  bool checkDaughtership(MELAParticle const* myParticle)const;
  void createAssociatedVs(std::vector<MELAParticle*>& particleArray);

  void addAssociatedParticleToArray(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray);
  void addAssociatedParticleToArray(MELAThreeBodyDecayCandidate* myParticle, std::vector<MELAThreeBodyDecayCandidate*>& particleArray);

};



#endif
