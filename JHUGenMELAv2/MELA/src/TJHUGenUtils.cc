#include "TModHiggsJJMatEl.hh"
#include "TJHUGenUtils.hh"
#include "TMath.h"

using namespace std;
using namespace TNumericUtil;

std::vector<intTriplet_t> TJHUGenUtils::Hash_OnshellHJJHash(){
  std::vector<intTriplet_t> pcfg;
  int ijsel[3][121];
  int nijchannels=77;
  __modhiggsjj_MOD_get_hjjchannelhash_nosplit(ijsel, &nijchannels);
  for (int ic=0; ic<nijchannels; ic++) pcfg.push_back(intTriplet_t(ijsel[0][ic], ijsel[1][ic], ijsel[2][ic]));
  return pcfg;
}
std::vector<intTriplet_t> TJHUGenUtils::Hash_OnshellVBFHash(){
  std::vector<intTriplet_t> pcfg;
  int ijsel[3][121];
  int nijchannels=68;
  __modhiggsjj_MOD_get_vbfchannelhash_nosplit(ijsel, &nijchannels);
  for (int ic=0; ic<nijchannels; ic++) pcfg.push_back(intTriplet_t(ijsel[0][ic], ijsel[1][ic], ijsel[2][ic]));
  return pcfg;
}
