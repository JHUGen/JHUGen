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
  std::vector<TNumericUtil::intTriplet_t> Hash_OnshellTQHHash();

  const std::vector<TNumericUtil::intTriplet_t>& Get_JHUGenHash_OnshellHJJHash();
  const std::vector<TNumericUtil::intTriplet_t>& Get_JHUGenHash_OnshellVBFHash();
  const std::vector<TNumericUtil::intTriplet_t>& Get_JHUGenHash_OnshellTQHHash();
}

#endif

