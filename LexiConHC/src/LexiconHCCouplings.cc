#include <cassert>
#include <iostream>
#include <string>
#include "LexiConHCCouplings.h"


using namespace std;


#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) \
case coupl_##PREFIX##_##NAME: \
  return #NAME;

std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::Amplitude_JHUGen_CouplingType type){
  switch (type){
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::EFT_JHUGen_CouplingType type){
  switch (type){
    EFT_JHUGEN_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::EFT_HiggsBasis_CouplingType type){
  switch (type){
    EFT_HIGGSBASIS_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}

#undef COUPLING_COMMAND
