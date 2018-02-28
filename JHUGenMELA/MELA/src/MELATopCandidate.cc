#include "MELATopCandidate.h"
#include "TUtilHelpers.hh"


MELATopCandidate::MELATopCandidate(
  MELAParticle* lightQuark_,
  MELAParticle* Wferm_,
  MELAParticle* Wfermbar_
  ) : MELAParticle(),
  lightQuark(lightQuark_),
  Wferm(Wferm_),
  Wfermbar(Wfermbar_)
{
  if (lightQuark!=0){
    p4 = p4 + lightQuark->p4;
    addDaughter(lightQuark);
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
    (lightQuark!=0 && PDGHelpers::isDownTypeQuark(lightQuark->id) && lightQuark->id>0)
    ||
    (Wcharge>0)
    ) id = 6;
  else if (
    (lightQuark!=0 && PDGHelpers::isDownTypeQuark(lightQuark->id) && lightQuark->id<0)
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
  std::swap(lightQuark, particle_.lightQuark);
  std::swap(Wferm, particle_.Wferm);
  std::swap(Wfermbar, particle_.Wfermbar);
}

void MELATopCandidate::setLightQuark(MELAParticle* myParticle){ lightQuark=myParticle; if (lightQuark!=0) addDaughter(lightQuark); }
void MELATopCandidate::setWFermion(MELAParticle* myParticle){ Wferm=myParticle; if (Wferm!=0) addDaughter(Wferm); }
void MELATopCandidate::setWAntifermion(MELAParticle* myParticle){ Wfermbar=myParticle; if (Wfermbar!=0) addDaughter(Wfermbar); }

bool MELATopCandidate::checkTopCandidateExists(MELATopCandidate* myParticle, std::vector<MELATopCandidate*>& particleArray){
  return TUtilHelpers::checkElementExists<MELATopCandidate const*, MELATopCandidate*>(myParticle, particleArray);
}

