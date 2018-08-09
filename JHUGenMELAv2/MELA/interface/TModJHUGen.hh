#ifndef _TMODJHUGENMAIN_HH_
#define _TMODJHUGENMAIN_HH_

extern "C" {
  void __modjhugen_MOD_initfirsttime(char* gridfilename, int* lenfilename, int* irep);
  void __modjhugen_MOD_resetpdfs(char* gridfilename, int* lenfilename, int* irep, int* pdfid); // pdfid=1: CTEQ6L1   2: MRSW with best fit, 2xx: MSTW with eigenvector set xx=01..40    3: NNPDF3.0
}

#endif

