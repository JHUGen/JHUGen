#ifndef _TMODKINEMATICS_HH_
#define _TMODKINEMATICS_HH_

extern "C" {
  double __modkinematics_MOD_getbwpropagator(double* shat, int* scheme);
  double __modkinematics_MOD_reweightbwpropagator(double* shat);
  void __modkinematics_MOD_setrunningscales(double p[3][4], int id[4]);
  void __modkinematics_MOD_setpdfs(double* x1, double* x2, double pdf[2][13]);
}

#endif

