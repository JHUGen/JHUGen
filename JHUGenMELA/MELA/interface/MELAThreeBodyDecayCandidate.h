#ifndef MELATHREEBODYDECAYCANDIDATE_H
#define MELATHREEBODYDECAYCANDIDATE_H

#include "MELAParticle.h"


class MELAThreeBodyDecayCandidate : public MELAParticle{
public:
  MELAThreeBodyDecayCandidate() : MELAParticle(), partnerParticle(0), Wferm(0), Wfermbar(0) {}
  MELAThreeBodyDecayCandidate(int id_, TLorentzVector p4_) : MELAParticle(id_, p4_), partnerParticle(0), Wferm(0), Wfermbar(0) {}
  MELAThreeBodyDecayCandidate(MELAParticle* partnerParticle_, MELAParticle* Wferm_, MELAParticle* Wfermbar_);
  MELAThreeBodyDecayCandidate(const MELAThreeBodyDecayCandidate& particle_) : MELAParticle(particle_), partnerParticle(particle_.partnerParticle), Wferm(particle_.Wferm), Wfermbar(particle_.Wfermbar) {}
  MELAThreeBodyDecayCandidate& operator=(const MELAThreeBodyDecayCandidate& particle_);
  ~MELAThreeBodyDecayCandidate(){}
  void swap(MELAThreeBodyDecayCandidate& particle_);

  void setPartnerParticle(MELAParticle* myParticle);
  void setWFermion(MELAParticle* myParticle);
  void setWAntifermion(MELAParticle* myParticle);

  MELAParticle* getPartnerParticle(){ return partnerParticle; }
  MELAParticle* getWFermion(){ return Wferm; }
  MELAParticle* getWAntifermion(){ return Wfermbar; }

  MELAParticle* getPartnerParticle()const{ return partnerParticle; }
  MELAParticle* getWFermion()const{ return Wferm; }
  MELAParticle* getWAntifermion()const{ return Wfermbar; }

  void testPreSelectedDaughters();

  double getWmass() const;

  static bool checkCandidateExists(MELAThreeBodyDecayCandidate const* myParticle, std::vector<MELAThreeBodyDecayCandidate*> const& particleArray);

protected:
  MELAParticle* partnerParticle;
  MELAParticle* Wferm;
  MELAParticle* Wfermbar;

};


typedef MELAThreeBodyDecayCandidate MELATopCandidate_t;
typedef MELAThreeBodyDecayCandidate MELATauCandidate_t;

#endif
