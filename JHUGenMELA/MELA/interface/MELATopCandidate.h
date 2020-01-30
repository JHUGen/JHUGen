#ifndef MELATOPCANDIDATE_H
#define MELATOPCANDIDATE_H

#include "MELAParticle.h"


class MELATopCandidate : public MELAParticle{
public:
  MELATopCandidate() : MELAParticle(), partnerParticle(0), Wferm(0), Wfermbar(0) {}
  MELATopCandidate(int id_, TLorentzVector p4_) : MELAParticle(id_, p4_), partnerParticle(0), Wferm(0), Wfermbar(0) {}
  MELATopCandidate(MELAParticle* partnerParticle_, MELAParticle* Wferm_, MELAParticle* Wfermbar_);
  MELATopCandidate(const MELATopCandidate& particle_) : MELAParticle(particle_), partnerParticle(particle_.partnerParticle), Wferm(particle_.Wferm), Wfermbar(particle_.Wfermbar) {}
  MELATopCandidate& operator=(const MELATopCandidate& particle_);
  ~MELATopCandidate(){}
  void swap(MELATopCandidate& particle_);

  void setPartnerParticle(MELAParticle* myParticle);
  void setWFermion(MELAParticle* myParticle);
  void setWAntifermion(MELAParticle* myParticle);

  MELAParticle* getPartnerParticle(){ return partnerParticle; }
  MELAParticle* getWFermion(){ return Wferm; }
  MELAParticle* getWAntifermion(){ return Wfermbar; }

  static bool checkTopCandidateExists(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray);

protected:
  MELAParticle* partnerParticle;
  MELAParticle* Wferm;
  MELAParticle* Wfermbar;

};


#endif
