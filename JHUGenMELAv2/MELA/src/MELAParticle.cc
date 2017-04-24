#include <algorithm>
#include <utility>
#include "MELAParticle.h"

using namespace PDGHelpers;

namespace debugVars{
  bool debugFlag = false;
}


MELAParticle::MELAParticle():
id(0),
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
mothers(particle.mothers),
daughters(particle.daughters)
{}
MELAParticle& MELAParticle::operator=(const MELAParticle& particle_){
  MELAParticle tmp(particle_);
  swap(tmp);
  return *this;
}

void swap(MELAParticle& particle_){
  std::swap(id, particle_.id);
  std::swap(p4, particle_.p4);
  std::swap(passSelection, particle_.passSelection);
  std::swap(genStatus, particle_.genStatus);
  std::swap(lifetime, particle_.lifetime);
  std::swap(mothers, particle_.mothers);
  std::swap(daughters, particle_.daughters);
}

bool MELAParticle::checkParticleExists(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray){
  for (std::vector<MELAParticle*>::iterator it = particleArray.begin(); it<particleArray.end(); it++){
    if ((*it)==myParticle) return true;
  }
  return false;
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
  for (unsigned int idau=0; idau<daughters.size(); idau++){
    if (daughters.at(idau)!=0) result.push_back(daughters.at(idau)->id);
  }
  return result;
}

double MELAParticle::charge()const{
  double cpos=0;
  if (isAWBoson(id) || abs(id)==37 || abs(id)==2212 || abs(id)==211 || abs(id)==321 || abs(id)==411 || abs(id)==521) cpos = 1.;
  else if (isALepton(id)) cpos = -1.;
  else if (isUpTypeQuark(id)) cpos = 2./3.;
  else if (isDownTypeQuark(id)) cpos = -1./3.;
  if (id<0) cpos *= -1.;
  return cpos;
}


