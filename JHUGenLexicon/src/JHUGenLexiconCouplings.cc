#include <cassert>
#include <iostream>
#include <string>
#include "JHUGenLexiconCouplings.h"


using namespace std;


#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) \
case coupl_##PREFIX##_##NAME: \
  return #NAME;

std::string JHUGenLexiconCouplings::getCouplingName(JHUGenLexiconCouplings::Amplitude_JHUGen_CouplingType type){
  switch (type){
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS
  default:
    cerr << "JHUGenLexiconCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string JHUGenLexiconCouplings::getCouplingName(JHUGenLexiconCouplings::Amplitude_JHUGen_Include_Triple_CouplingType type){
  switch (type){
    AMPLITUDE_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS
  default:
    cerr << "JHUGenLexiconCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string JHUGenLexiconCouplings::getCouplingName(JHUGenLexiconCouplings::EFT_JHUGen_CouplingType type){
  switch (type){
    EFT_JHUGEN_COUPLING_COMMANDS
  default:
    cerr << "JHUGenLexiconCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string JHUGenLexiconCouplings::getCouplingName(JHUGenLexiconCouplings::EFT_JHUGen_Include_Triple_CouplingType type){
  switch (type){
    EFT_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS
  default:
    cerr << "JHUGenLexiconCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string JHUGenLexiconCouplings::getCouplingName(JHUGenLexiconCouplings::HiggsBasis_CouplingType type){
  switch (type){
    HIGGSBASIS_COUPLING_COMMANDS
  default:
    cerr << "JHUGenLexiconCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string JHUGenLexiconCouplings::getCouplingName(JHUGenLexiconCouplings::HiggsBasis_Include_Triple_CouplingType type){
  switch (type){
    HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS
  default:
    cerr << "JHUGenLexiconCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string JHUGenLexiconCouplings::getCouplingName(JHUGenLexiconCouplings::EFT_HiggsBasis_CouplingType type){
  switch (type){
    EFT_HIGGSBASIS_COUPLING_COMMANDS
  default:
    cerr << "JHUGenLexiconCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string JHUGenLexiconCouplings::getCouplingName(JHUGenLexiconCouplings::EFT_HiggsBasis_Include_Triple_CouplingType type){
  switch (type){
    EFT_HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS
  default:
    cerr << "JHUGenLexiconCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string JHUGenLexiconCouplings::getCouplingName(JHUGenLexiconCouplings::WarsawBasis_CouplingType type){
  switch (type){
    WARSAWBASIS_COUPLING_COMMANDS
  default:
    cerr << "JHUGenLexiconCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}



#undef COUPLING_COMMAND
