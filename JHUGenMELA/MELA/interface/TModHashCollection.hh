#ifndef _TMODHASHCOLLECTION_HH_
#define _TMODHASHCOLLECTION_HH_

extern "C" {

  void __modhashcollection_MOD_get_vbfchannelhash_nosplit(int ijSel[3][121], int* nijchannels);
  void __modhashcollection_MOD_get_hjjchannelhash_nosplit(int ijSel[3][121], int* nijchannels);
  void __modhashcollection_MOD_get_thchannelhash(int ijSel[3][121]);

}

#endif

