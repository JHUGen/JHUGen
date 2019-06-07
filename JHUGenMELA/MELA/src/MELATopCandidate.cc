#include "MELATopCandidate.h"
#include "TUtilHelpers.hh"


MELATopCandidate::MELATopCandidate(
  MELAParticle* partnerParticle_,
  MELAParticle* Wferm_,
  MELAParticle* Wfermbar_
  ) : MELAParticle(),
  partnerParticle(partnerParticle_),
  Wferm(Wferm_),
  Wfermbar(Wfermbar_)
{
  if (partnerParticle!=0){
    p4 = p4 + partnerParticle->p4;
    addDaughter(partnerParticle);
  }
  double Wcharge=0;
  if (Wferm!=0){
    p4 = p4 + Wferm->p4;
    Wcharge += Wferm->charge();
    addDaughter(Wferm);
  }
  if (Wfermbar!=0){
    p4 = p4 + Wfermbar->p4;
    Wcharge += Wfermbar->charge();
    addDaughter(Wfermbar);
  }
  if (
    (partnerParticle!=0 && PDGHelpers::isDownTypeQuark(partnerParticle->id) && partnerParticle->id>0)
    ||
    (Wcharge>0)
    ) id = 6;
  else if (
    (partnerParticle!=0 && PDGHelpers::isDownTypeQuark(partnerParticle->id) && partnerParticle->id<0)
    ||
    (Wcharge<0)
    ) id = -6;
  else id=0;
}
MELATopCandidate& MELATopCandidate::operator=(const MELATopCandidate& particle_){
  MELATopCandidate tmp(particle_);
  swap(tmp);
  return *this;
}

void MELATopCandidate::swap(MELATopCandidate& particle_){
  MELAParticle::swap(particle_);
  std::swap(partnerParticle, particle_.partnerParticle);
  std::swap(Wferm, particle_.Wferm);
  std::swap(Wfermbar, particle_.Wfermbar);
}

void MELATopCandidate::setPartnerParticle(MELAParticle* myParticle){ partnerParticle=myParticle; if (partnerParticle!=0) addDaughter(partnerParticle); }
void MELATopCandidate::setWFermion(MELAParticle* myParticle){ Wferm=myParticle; if (Wferm!=0) addDaughter(Wferm); }
void MELATopCandidate::setWAntifermion(MELAParticle* myParticle){ Wfermbar=myParticle; if (Wfermbar!=0) addDaughter(Wfermbar); }

bool MELATopCandidate::checkTopCandidateExists(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray){
  return TUtilHelpers::checkElementExists<MELATopCandidate const*, MELATopCandidate*>(myParticle, particleArray);
}

