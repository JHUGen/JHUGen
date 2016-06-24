#ifndef _TMODPARAMETERS_HH_
#define _TMODPARAMETERS_HH_

extern "C" {

  // From PDG ids to JHU ids
  int __modparameters_MOD_convertlhereverse(int* part);
#define convertLHEreverse __modparameters_MOD_convertlhereverse

  // From JHU ids to PDG ids
  int __modparameters_MOD_convertlhe(int* part);
#define convertLHE __modparameters_MOD_convertlhe

  double __modparameters_MOD_ckm(int* id1in, int* id2in); // |VCKM|*ScaleFactor
  double __modparameters_MOD_scalefactor(int* id1in, int* id2in); // ScaleFactor
}

#endif

