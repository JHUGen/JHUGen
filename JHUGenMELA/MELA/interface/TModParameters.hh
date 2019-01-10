#ifndef _TMODPARAMETERS_HH_
#define _TMODPARAMETERS_HH_

extern "C" {

  // Get the mass and width of the particles
  double __modparameters_MOD_getmass(int* part);
  double __modparameters_MOD_getdecaywidth(int* part);

  // From PDG ids to JHU ids
  int __modparameters_MOD_convertlhereverse(int* part);
#define convertLHEreverse __modparameters_MOD_convertlhereverse

  // From JHU ids to PDG ids
  int __modparameters_MOD_convertlhe(int* part);
#define convertLHE __modparameters_MOD_convertlhe

  int __modparameters_MOD_coupledvertex(int id[2], int* hel, int* useahcoupl);
#define CoupledVertex __modparameters_MOD_coupledvertex

  double __modparameters_MOD_scalefactor(int* id1in, int* id2in); // ScaleFactor
  double __modparameters_MOD_ckmbare(int* id1in, int* id2in); // VCKM
  double __modparameters_MOD_ckm(int* id1in, int* id2in); // VCKM*ScaleFactor

  void __modparameters_MOD_setmass(double* mass, int* ipart);
  void __modparameters_MOD_setdecaywidth(double* width, int* ipart);

  void __modparameters_MOD_computeckmelements(double* invckm_ud, double* invckm_us, double* invckm_cd, double* invckm_cs, double* invckm_ts, double* invckm_tb, double* invckm_ub=0, double* invckm_cb=0, double* invckm_td=0);
  void __modparameters_MOD_computeewvariables();
  void __modparameters_MOD_computeqcdvariables();

  void __modparameters_MOD_evalalphas();
}

#endif

