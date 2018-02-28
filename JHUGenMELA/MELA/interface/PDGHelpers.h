#ifndef PDGHELPERS_H
#define PDGHELPERS_H

#include <iostream>
#include <cmath>
#include <vector>
// Mod_Parameters
#include "TModParameters.hh"
// TVar namespace
#include "TVar.hh"

namespace PDGHelpers{
  const double Wmass = 80.399;
  const double Zmass = 91.1876;
  const double Topmass = 173.2;
  const double Bottommass = 4.75;
  const double Zeromass = 0;
  const double m_el = 0.00051100;
  const double m_mu = 0.10566;
  const double m_tau = 1.7768;
  const double Wwidth = 2.085;
  const double Zwidth = 2.4952;
  const double Topwidth = 2.;

  extern TVar::CandidateDecayMode HDecayMode;

  bool isALepton(const int id);
  bool isANeutrino(const int id);
  bool isAJet(const int id);
  bool isAnUnknownJet(const int id);
  bool isInvalid(const int id);
  bool isAQuark(const int id);
  bool isUpTypeQuark(const int id);
  bool isDownTypeQuark(const int id);
  bool isAGluon(const int id);
  bool isAPhoton(const int id);
  bool isAZBoson(const int id);
  bool isAWBoson(const int id);
  bool isAHiggs(const int id);

  void orderParticles(
    const std::vector<int>& idlist,
    const std::vector<bool(*)(const int)>& testlist,
    std::vector<int>& ordering,
    bool allowUnknown=false
    );
  void groupIdenticalParticles(
    const std::vector<int>& ids,
    std::vector<std::vector<int>>& ordering,
    bool* hasUnknownParticles=0
    );
  void pairIdenticalParticles(
    const std::vector<int>& ids,
    std::vector<std::pair<int, int>>& ordering,
    bool allowUnknown=false
    );

  bool allEquivalent(std::vector<int> ids, bool allowUnknown=false);


  void setCandidateDecayMode(TVar::CandidateDecayMode mode);
  int getCoupledVertex(const int idfirst, const int idsecond, int* hel=0, int* useAHcoupl=0);

  int convertPythiaStatus(int pSt);
}

#endif
