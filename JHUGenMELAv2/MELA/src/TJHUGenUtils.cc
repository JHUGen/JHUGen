#include "TModHashCollection.hh"
#include "TJHUGenUtils.hh"
#include "TMath.h"

using namespace std;
using namespace TNumericUtil;

namespace TJHUGenUtils{
  const std::vector<TNumericUtil::intTriplet_t> JHUGenHash_OnshellHJJHash = Hash_OnshellHJJHash();
  const std::vector<TNumericUtil::intTriplet_t> JHUGenHash_OnshellVBFHash = Hash_OnshellVBFHash();
  const std::vector<TNumericUtil::intTriplet_t> JHUGenHash_OnshellTQHHash = Hash_OnshellTQHHash();
}

std::vector<intTriplet_t> TJHUGenUtils::Hash_OnshellHJJHash(){
  std::vector<intTriplet_t> pcfg;
  int ijsel[3][121];
  int nijchannels=77;
  __modhashcollection_MOD_get_hjjchannelhash_nosplit(ijsel, &nijchannels);
  for (int ic=0; ic<nijchannels; ic++) pcfg.push_back(intTriplet_t(ijsel[0][ic], ijsel[1][ic], ijsel[2][ic]));
  return pcfg;
}
std::vector<intTriplet_t> TJHUGenUtils::Hash_OnshellVBFHash(){
  std::vector<intTriplet_t> pcfg;
  int ijsel[3][121];
  int nijchannels=68;
  __modhashcollection_MOD_get_vbfchannelhash_nosplit(ijsel, &nijchannels);
  for (int ic=0; ic<nijchannels; ic++) pcfg.push_back(intTriplet_t(ijsel[0][ic], ijsel[1][ic], ijsel[2][ic]));
  return pcfg;
}
std::vector<intTriplet_t> TJHUGenUtils::Hash_OnshellTQHHash(){
  std::vector<intTriplet_t> pcfg;
  int ijsel[3][121];
  __modhashcollection_MOD_get_thchannelhash(ijsel);
  for (int ic=0; ic<121; ic++){
    if (ijsel[2][ic]>=0) pcfg.push_back(intTriplet_t(ijsel[0][ic], ijsel[1][ic], ijsel[2][ic]));
    else break;
  }
  return pcfg;
}

const std::vector<TNumericUtil::intTriplet_t>& TJHUGenUtils::Get_JHUGenHash_OnshellHJJHash(){ return JHUGenHash_OnshellHJJHash; }
const std::vector<TNumericUtil::intTriplet_t>& TJHUGenUtils::Get_JHUGenHash_OnshellVBFHash(){ return JHUGenHash_OnshellVBFHash; }
const std::vector<TNumericUtil::intTriplet_t>& TJHUGenUtils::Get_JHUGenHash_OnshellTQHHash(){ return JHUGenHash_OnshellTQHHash; }
