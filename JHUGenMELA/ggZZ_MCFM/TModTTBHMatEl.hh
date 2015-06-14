
#ifndef _TMODTTBHMATEL_HH_
#define _TMODTTBHMATEL_HH_

extern "C" {
  void __modttbh_MOD_initprocess_ttbh(double *MReso, double *MFerm);
//  void __modttbh_MOD_exitprocess_ttbh();
  void __modttbh_MOD_evalamp_gg_ttbh(double Ptth[13][4], double TTBHcoupl[2][2], int *TopDecays, double *MatElSq);
  void __modttbh_MOD_evalamp_qqb_ttbh(double Ptth[13][4], double TTBHcoupl[2][2], int *TopDecays, double *MatElSq);
  void __modttbh_MOD_evalxsec_pp_ttbh(double Ptth[13][4], double TTBHcoupl[2][2], int *TopDecays, int* SelectProcess, double *MatElSq);
  void  __modttbh_MOD_evalxsec_pp_bbbh(double Mom[13][4], double BBBHcoupl[2][2], int *SelectProcess, double *Res);
  void  __modttbh_MOD_evalamp_gg_bbbh(double Mom[13][4], double BBBHcoupl[2][2], double *SqAmp);
  void  __modttbh_MOD_evalamp_qqb_bbbh(double Mom[13][4], double BBBHcoupl[2][2], double *SqAmp);

}

#endif
