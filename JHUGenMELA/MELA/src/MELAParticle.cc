#include "MELAParticle.h"
#include "TUtilHelpers.hh"
#include "MELAStreamHelpers.hh"


using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;
using namespace PDGHelpers;


namespace debugVars{
  bool debugFlag = false;
}


MELAParticle::MELAParticle():
id(-9000),
p4(0, 0, 0, 0),
passSelection(true),
genStatus(-2),
lifetime(0)
{}
MELAParticle::MELAParticle(int id_) :
id(id_),
p4(0, 0, 0, 0),
passSelection(true),
genStatus(-2),
lifetime(0)
{}
MELAParticle::MELAParticle(int id_, TLorentzVector p4_) :
id(id_),
p4(p4_),
passSelection(true),
genStatus(-2),
lifetime(0)
{}
MELAParticle::MELAParticle(const MELAParticle& particle_) :
id(particle_.id),
p4(particle_.p4),
passSelection(particle_.passSelection),
genStatus(particle_.genStatus),
lifetime(particle_.lifetime),
mothers(particle_.mothers),
daughters(particle_.daughters)
{}
MELAParticle& MELAParticle::operator=(const MELAParticle& particle_){
  MELAParticle tmp(particle_);
  swap(tmp);
  return *this;
}

void MELAParticle::swap(MELAParticle& particle_){
  std::swap(id, particle_.id);
  std::swap(p4, particle_.p4);
  std::swap(passSelection, particle_.passSelection);
  std::swap(genStatus, particle_.genStatus);
  std::swap(lifetime, particle_.lifetime);
  std::swap(mothers, particle_.mothers);
  std::swap(daughters, particle_.daughters);
}

void MELAParticle::addMother(MELAParticle* myParticle){ if (!checkParticleExists(myParticle, mothers)) mothers.push_back(myParticle); }
void MELAParticle::addDaughter(MELAParticle* myParticle){ if (!checkParticleExists(myParticle, daughters)) daughters.push_back(myParticle); }
MELAParticle* MELAParticle::getMother(int index)const{
  if ((int)mothers.size()>index) return mothers.at(index);
  else return 0;
}
MELAParticle* MELAParticle::getDaughter(int index)const{
  if ((int)daughters.size()>index) return daughters.at(index);
  else return 0;
}
std::vector<int> MELAParticle::getDaughterIds()const{
  std::vector<int> result;
  for (auto& dau:daughters){ if (dau) result.push_back(dau->id); }
  return result;
}
void MELAParticle::getRelatedParticles(std::vector<MELAParticle*>& particles) const{
  for (auto* part:mothers){ if (part) part->getRelatedParticles(particles); }
  for (auto* part:daughters){ if (part) part->getRelatedParticles(particles); }
  if (!checkParticleExists(this, particles)) particles.push_back((MELAParticle*) this);
}
void MELAParticle::getDaughterParticles(std::vector<MELAParticle*>& particles) const{
  for (auto* part:daughters){ if (part) part->getDaughterParticles(particles); }
  if (!checkParticleExists(this, particles)) particles.push_back((MELAParticle*) this); // A particle is its own daughter. This is to ensure that particle chains are reported properly.
}

bool MELAParticle::hasMother(MELAParticle const* part) const{ return MELAParticle::checkParticleExists(part, mothers); }
bool MELAParticle::hasDaughter(MELAParticle const* part) const{ return MELAParticle::checkParticleExists(part, daughters); }

double MELAParticle::charge()const{
  double cpos=0;
  if (isAWBoson(id) || abs(id)==37 || abs(id)==2212 || abs(id)==211 || abs(id)==321 || abs(id)==411 || abs(id)==521) cpos = 1.;
  else if (isALepton(id)) cpos = -1.;
  else if (isUpTypeQuark(id)) cpos = 2./3.;
  else if (isDownTypeQuark(id)) cpos = -1./3.;
  if (id<0) cpos *= -1.;
  return cpos;
}

void MELAParticle::boost(const TVector3& vec, bool boostAll){
  if (vec.Mag2()<1.){
    if (boostAll){
      std::vector<MELAParticle*> particles;
      this->getRelatedParticles(particles);
      for (std::vector<MELAParticle*>::iterator it = particles.begin(); it<particles.end(); it++) (*it)->boost(vec, false);
    }
    else p4.Boost(vec);
  }
  else MELAerr
    << "MELAParticle::boost: "
    << "|v|**2 = " << vec.Mag2() << " >=1 cannot be used to boost. "
    << "v = ( " << vec.X() << " , " << vec.Y() << " , " << vec.Z() << " )"
    << std::endl;
}

bool MELAParticle::checkParticleExists(MELAParticle const* myParticle, std::vector<MELAParticle*> const& particleArray){
  return TUtilHelpers::checkElementExists<MELAParticle const*, MELAParticle*>(myParticle, particleArray);
}
bool MELAParticle::checkDeepDaughtership(MELAParticle const* part1, MELAParticle const* part2){
  bool res = false;
  if (!part1 || !part2) return res;
  std::vector<MELAParticle*> daughters1; part1->getDaughterParticles(daughters1);
  std::vector<MELAParticle*> daughters2; part2->getDaughterParticles(daughters2);
  return TUtilHelpers::hasCommonElements(daughters1, daughters2);
}

TVector3 MELAParticle::calculateTotalDisplacement()const{
  TVector3 res;
  double const part_mass = this->m();
  TVector3 const part_vec = this->vect();

  // Calculate displacement
  if (part_mass>0.) res = part_vec * (lifetime/part_mass);
  //MELAout << "Displacement for id=" << this->id << " = " << res.X() << ", " << res.Y() << ", " << res.Z() << " [mom: " << this->p4 << ", lifetime: " << lifetime << "]" << std::endl;

  // Add mother displacement
  if (this->getNMothers()==1) res += this->getMother(0)->calculateTotalDisplacement();

  return res;
}

