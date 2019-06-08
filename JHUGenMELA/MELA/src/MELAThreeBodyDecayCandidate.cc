#include "MELAThreeBodyDecayCandidate.h"
#include "TUtilHelpers.hh"


MELAThreeBodyDecayCandidate::MELAThreeBodyDecayCandidate(
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
MELAThreeBodyDecayCandidate& MELAThreeBodyDecayCandidate::operator=(const MELAThreeBodyDecayCandidate& particle_){
  MELAThreeBodyDecayCandidate tmp(particle_);
  swap(tmp);
  return *this;
}

void MELAThreeBodyDecayCandidate::swap(MELAThreeBodyDecayCandidate& particle_){
  MELAParticle::swap(particle_);
  std::swap(partnerParticle, particle_.partnerParticle);
  std::swap(Wferm, particle_.Wferm);
  std::swap(Wfermbar, particle_.Wfermbar);
}

void MELAThreeBodyDecayCandidate::setPartnerParticle(MELAParticle* myParticle){ partnerParticle=myParticle; if (partnerParticle!=0) addDaughter(partnerParticle); }
void MELAThreeBodyDecayCandidate::setWFermion(MELAParticle* myParticle){ Wferm=myParticle; if (Wferm!=0) addDaughter(Wferm); }
void MELAThreeBodyDecayCandidate::setWAntifermion(MELAParticle* myParticle){ Wfermbar=myParticle; if (Wfermbar!=0) addDaughter(Wfermbar); }

double MELAThreeBodyDecayCandidate::getWmass() const{
  if (Wferm && Wfermbar) return (Wferm->p4+Wfermbar->p4).M();
  else return -1;
}

bool MELAThreeBodyDecayCandidate::checkCandidateExists(MELAThreeBodyDecayCandidate const* myParticle, std::vector<MELAThreeBodyDecayCandidate*> const& particleArray){
  return TUtilHelpers::checkElementExists<MELAThreeBodyDecayCandidate const*, MELAThreeBodyDecayCandidate*>(myParticle, particleArray);
}

void MELAThreeBodyDecayCandidate::testPreSelectedDaughters(){
  for (auto& dau:getDaughters()){
    if (!dau->passSelection){
      passSelection=false;
      break;
    }
  }
}

