#ifndef JHUGEN_UTILITIES
#define JHUGEN_UTILITIES
#include <string>
#include <vector>
// PDGHelpers
#include "PDGHelpers.h"
// Numeric Util
#include "TNumericUtil.hh"
// ROOT includes
#include "TLorentzVector.h"


namespace TJHUGenUtils{
  std::vector<TNumericUtil::intTriplet_t> Hash_OnshellHJJHash();
  std::vector<TNumericUtil::intTriplet_t> Hash_OnshellVBFHash();

  const std::vector<TNumericUtil::intTriplet_t> JHUGenHash_OnshellHJJHash = Hash_OnshellHJJHash();
  const std::vector<TNumericUtil::intTriplet_t> JHUGenHash_OnshellVBFHash = Hash_OnshellVBFHash();
}

#endif

