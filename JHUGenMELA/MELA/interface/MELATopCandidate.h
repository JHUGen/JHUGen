#ifndef MELATOPCANDIDATE_H
#define MELATOPCANDIDATE_H

#include "MELAParticle.h"


class MELATopCandidate : public MELAParticle{
public:
  MELATopCandidate() : MELAParticle(), lightQuark(0), Wferm(0), Wfermbar(0) {}
  MELATopCandidate(int id_, TLorentzVector p4_) : MELAParticle(id_, p4_), lightQuark(0), Wferm(0), Wfermbar(0) {}
  MELATopCandidate(MELAParticle* lightQuark_, MELAParticle* Wferm_, MELAParticle* Wfermbar_);
  MELATopCandidate(const MELATopCandidate& particle_) : MELAParticle(particle_), lightQuark(particle_.lightQuark), Wferm(particle_.Wferm), Wfermbar(particle_.Wfermbar) {}
  MELATopCandidate& operator=(const MELATopCandidate& particle_);
  ~MELATopCandidate(){}
  void swap(MELATopCandidate& particle_);

  void setLightQuark(MELAParticle* myParticle);
  void setWFermion(MELAParticle* myParticle);
  void setWAntifermion(MELAParticle* myParticle);

  MELAParticle* getLightQuark(){ return lightQuark; }
  MELAParticle* getWFermion(){ return Wferm; }
  MELAParticle* getWAntifermion(){ return Wfermbar; }

  static bool checkTopCandidateExists(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray);

protected:
  MELAParticle* lightQuark;
  MELAParticle* Wferm;
  MELAParticle* Wfermbar;

};


#endif
