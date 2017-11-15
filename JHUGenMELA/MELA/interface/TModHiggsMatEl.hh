#ifndef _TMODHIGGSMATEL_HH_
#define _TMODHIGGSMATEL_HH_

extern "C" {
  void __modhiggs_MOD_evalamp_gg_h_vv(double P[6][4], int *MYIDUP, double *MatElSq);
  void __modhiggs_MOD_evalamp_h_vv(double P[6][4], int *MYIDUP, double *MatElSq);
  void __modhiggs_MOD_evalamp_h_ff(double P[2][4], double* mass_f, double *MatElSq);
  void __modhiggs_MOD_evalamp_h_tt_decay(double P[6][4], double* mass_f, double* ga_f, double *MatElSq);
}

#endif

