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

  // Data members
public:
  int id;
  TLorentzVector p4;
  bool passSelection;
  int genStatus;
  double lifetime;

protected:
  std::vector<MELAParticle*> mothers;
  std::vector<MELAParticle*> daughters;


  // Member functions
public:

  // Constructors
  MELAParticle();
  MELAParticle(int id_);
  MELAParticle(int id_, TLorentzVector p4_);
  MELAParticle(const MELAParticle& particle_);
  MELAParticle& operator=(const MELAParticle& particle_);
  virtual ~MELAParticle(){};

  // Utility functions
  void swap(MELAParticle& particle_);

  void setSelected(bool isSelected=true){ passSelection = isSelected; }
  void setGenStatus(int status_){ genStatus=status_; }
  void setLifetime(double life_){ lifetime=life_; }

  void addMother(MELAParticle* myParticle);
  void addDaughter(MELAParticle* myParticle);

  int getNMothers() const{ return mothers.size(); };
  int getNDaughters() const{ return daughters.size(); };
  virtual std::vector<int> getDaughterIds()const;

  MELAParticle* getMother(int index) const;
  MELAParticle* getDaughter(int index) const;
  virtual void getRelatedParticles(std::vector<MELAParticle*>& particles) const;
  virtual void getDaughterParticles(std::vector<MELAParticle*>& particles) const;

  std::vector<MELAParticle*>& getMothers(){ return mothers; }
  std::vector<MELAParticle*>& getDaughters(){ return daughters; }
  const std::vector<MELAParticle*>& getMothers()const{ return mothers; }
  const std::vector<MELAParticle*>& getDaughters()const{ return daughters; }
  bool hasMother(MELAParticle const* part) const;
  bool hasDaughter(MELAParticle const* part) const;

  double charge()const;
  double m()const{ return p4.M(); }
  double x()const{ return p4.X(); }
  double y()const{ return p4.Y(); }
  double z()const{ return p4.Z(); }
  double t()const{ return p4.T(); }
  double p()const{ return p4.P(); }
  double pt()const{ return p4.Pt(); }
  double eta()const{ return p4.Eta(); }
  double phi()const{ return p4.Phi(); }
  double rapidity()const{ return p4.Rapidity(); }
  double dot(const TLorentzVector& v)const{ return p4.Dot(v); }
  double dot(const MELAParticle& part)const{ return dot(part.p4); }
  double dot(const MELAParticle* part)const{ if (part) return dot(*part); else return 0; }
  double euclidean_dot(const TLorentzVector& v)const{ return (p4.Vect().Dot(v.Vect())+p4.T()*v.T()); }
  double euclidean_dot(const MELAParticle& part)const{ return euclidean_dot(part.p4); }
  double euclidean_dot(const MELAParticle* part)const{ if (part) return euclidean_dot(*part); else return 0; }
  double deltaR(const TLorentzVector& v)const{ return p4.DeltaR(v); }
  double deltaR(const MELAParticle& part)const{ return deltaR(part.p4); }
  double deltaR(const MELAParticle* part)const{ if (part) return deltaR(*part); else return -1; }
  void boost(const TVector3& vec, bool boostAll=false);
  TVector3 vect()const{ return p4.Vect(); }

  // Function to calculate displacements from lifetimes
  TVector3 calculateTotalDisplacement()const;

  // Operators
  MELAParticle& operator+=(MELAParticle* part){ if (part){ p4 += part->p4; addDaughter(part); } return *this; }
  MELAParticle& operator+=(const TLorentzVector& mom){ p4 += mom; return *this; }

  // Helper functions
  static bool checkParticleExists(MELAParticle const* myParticle, std::vector<MELAParticle*> const& particleArray);
  static bool checkDeepDaughtership(MELAParticle const* part1, MELAParticle const* part2);

};


#endif
