#ifndef MCFM_UTILITIES
#define MCFM_UTILITIES
#include <string>
#include <vector>
// PDGHelpers
#include "PDGHelpers.h"
// Numeric Util
#include "TNumericUtil.hh"
// ROOT includes
#include "TLorentzVector.h"


namespace TMCFMUtils{
  void AssociatedParticleOrdering_QQVVQQAny(int iSel, int jSel, int rSel, int sSel, int order[2]);
  std::vector<TNumericUtil::intQuad_t> Hash_QQVVQQAny();
  std::vector<TNumericUtil::intQuad_t> Hash_QQVVQQ();
  std::vector<TNumericUtil::intQuad_t> Hash_QQVVQQStrong();
}

#endif

