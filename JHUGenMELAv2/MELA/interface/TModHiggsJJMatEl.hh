#ifndef _TMODHIGGSJJMATEL_HH_
#define _TMODHIGGSJJMATEL_HH_

extern "C" {

  void __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa(double P[5][4], double MatElSq[11][11]);
  void __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select(double P[5][4], int* iSel, int* jSel, bool* zz_fusion, int* iflip, double MatElSq[11][11]);
  void __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(double P[5][4], int* iSel, int* jSel, int* rSel, int* sSel, double* MatElSq);
  void __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa(double P[5][4], double MatElSq[11][11]);
  void __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select(double P[5][4], int* iSel, int* jSel, int* flav_tag, int* iflip, double MatElSq[11][11]);
  void __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(double P[5][4], int* iSel, int* jSel, int* rSel, int* sSel, double* MatElSq);

  void __modhiggsjj_MOD_get_vbfchannelhash_nosplit(int ijSel[3][121], int* nijchannels);
  void __modhiggsjj_MOD_get_hjjchannelhash_nosplit(int ijSel[3][121], int* nijchannels);

}

#endif
