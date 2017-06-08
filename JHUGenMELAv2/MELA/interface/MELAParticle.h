#ifndef MELAPARTICLE_H
#define MELAPARTICLE_H

#include <vector>
#include "PDGHelpers.h"
#include "TVector3.h"
#include "TLorentzVector.h"

namespace debugVars{
  extern bool debugFlag;
}

class MELAParticle{
public:

  // Constructors

  MELAParticle();
  MELAParticle(int id_, TLorentzVector p4_);
  MELAParticle(const MELAParticle& particle_);
  MELAParticle& operator=(const MELAParticle& particle_);
  virtual ~MELAParticle(){};

  // Data

  int id;
  TLorentzVector p4;
  bool passSelection;
  int genStatus;
  double lifetime;

  // Member functions
  void setSelected(bool isSelected=true){ passSelection = isSelected; }
  void setGenStatus(int status_){ genStatus=status_; }
  void setLifetime(int life_){ lifetime=life_; }

  void addMother(MELAParticle* myParticle);
  void addDaughter(MELAParticle* myParticle);

  int getNMothers() const{ return mothers.size(); };
  int getNDaughters() const{ return daughters.size(); };
  virtual std::vector<int> getDaughterIds()const;

  MELAParticle* getMother(int index) const;
  MELAParticle* getDaughter(int index) const;

  double charge()const;
  double m()const{ return p4.M(); }
  double x()const{ return p4.X(); }
  double y()const{ return p4.Y(); }
  double z()const{ return p4.Z(); }
  double t()const{ return p4.T(); }
  double pt()const{ return p4.Pt(); }
  double eta()const{ return p4.Eta(); }
  double phi()const{ return p4.Phi(); }
  double rapidity()const{ return p4.Rapidity(); }
  double deltaR(const TLorentzVector& v)const{ return p4.DeltaR(v); }


protected:
  std::vector<MELAParticle*> mothers;
  std::vector<MELAParticle*> daughters;

  bool checkParticleExists(MELAParticle* myParticle, std::vector<MELAParticle*>& particleArray);
};


#endif
