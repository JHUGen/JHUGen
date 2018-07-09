#ifndef _TMODJHUGENMELA_HH_
#define _TMODJHUGENMELA_HH_

// NOTE: LOGICAL==INT, LOGICAL*1==BOOL(!???): http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html

extern "C" {

  void __modjhugenmela_MOD_setdecaymodes(int idfirst[2], int idsecond[2]);
  void __modjhugenmela_MOD_setdistinguishwwcouplingsflag(int* doallow); // YES, THE ARGUMENT IS AN INT!
  void __modjhugenmela_MOD_sethdk(int* flag);
  void __modjhugenmela_MOD_sethiggsmasswidth(double *mass, double *width);
  void __modjhugenmela_MOD_setmurenfac(double* muren, double* mufac);
  void __modjhugenmela_MOD_resetmubarhgabarh();
  void __modjhugenmela_MOD_setmvgv();
  void __modjhugenmela_MOD_setmvgvfromvertex(int* idV);
  void __modjhugenmela_MOD_setspinonecouplings(double qqcoupl[2][2], double vvcoupl[2][2]);
  void __modjhugenmela_MOD_setspintwocouplings(double acoupl[5][2], double bcoupl[10][2], double qlr[2][2]);
  void __modjhugenmela_MOD_setspinzeroggcouplings(double ggcoupl[3][2]);
  void __modjhugenmela_MOD_setspinzeroqqcouplings(double qqcoupl[2][2]);
  void __modjhugenmela_MOD_setspinzerovvcouplings(double vvcoupl[39][2], int cqsq[3], double Lambda_qsq[4][3], int* usewwcoupl); // YES, THE LAST ARGUMENT IS AN INT!
  void __modjhugenmela_MOD_setspinzerovvcouplings_nogamma(double vvcoupl[32][2], int cqsq[3], double Lambda_qsq[4][3], int* usewwcoupl); // YES, THE LAST ARGUMENT IS AN INT!
  void __modjhugenmela_MOD_settopdecays(int* flag);

  void __modjhugenmela_MOD_getmvgv(double* mv, double* gv);
  void __modjhugenmela_MOD_getalphasalphasmz(double* val_as, double* val_asmz);
  void __modjhugenmela_MOD_getdecaycouplings(int* VVMode, int idordered[4], double* aL1, double* aR1, double* aL2, double* aR2);

}

#endif

