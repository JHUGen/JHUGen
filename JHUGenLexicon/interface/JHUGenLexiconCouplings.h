#ifndef JHUGenCouplings_COUPLINGS_H
#define JHUGenCouplings_COUPLINGS_H

#include <string>

// Default values of parameters as in mod_Parameters.F90
#define DEFVAL_MZ 91.1876
#define DEFVAL_MW 80.399
#define DEFVAL_SW 0.23119 //Actually is SW^2
#define DEFVAL_LAMBDA_VI 10000.
#define DEFVAL_ALPHA .0072973525693 //Used to set e^2
#define DEFVAL_ALPHA_S .11834 //Used to set gs^2
#define DEFVAL_VEV_LAM .0024622 //Vacuum Expectation Value GeV / Lambda Warsaw Basis


// Coupling definitions
// Added gluon couplings that appear in H->VV
#define AMPLITUDE_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS \
  COUPLING_COMMAND(ghz1, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghz1_prime2, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghz2, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghz4, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghw1, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghw1_prime2, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghw2, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghw4, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghzgs1_prime2, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghzgs2, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghzgs4, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghgsgs2, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghgsgs4, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghg2, ampjhutrip, 0.) \
  COUPLING_COMMAND(ghg4, ampjhutrip, 0.) \
  COUPLING_COMMAND(dV_Z, ampjhutrip, 0.) \
  COUPLING_COMMAND(dV_A, ampjhutrip, 0.) \
  COUPLING_COMMAND(dP_Z, ampjhutrip, 0.) \
  COUPLING_COMMAND(dP_A, ampjhutrip, 0.) \
  COUPLING_COMMAND(dM_Z, ampjhutrip, 0.) \
  COUPLING_COMMAND(dM_A, ampjhutrip, 0.) \
  COUPLING_COMMAND(dFour_Z, ampjhutrip, 0.) \
  COUPLING_COMMAND(dFour_A, ampjhutrip, 0.) \
  COUPLING_COMMAND(dZZWpWm, ampjhutrip, 0.) \
  COUPLING_COMMAND(dZAWpWm, ampjhutrip, 0.) \
  COUPLING_COMMAND(dAAWpWm, ampjhutrip, 0.) \

#define EFT_JHUGEN_COUPLING_COMMANDS \
  COUPLING_COMMAND(ghz1, eftjhu, 0.) \
  COUPLING_COMMAND(ghz1_prime2, eftjhu, 0.) \
  COUPLING_COMMAND(ghz2, eftjhu, 0.) \
  COUPLING_COMMAND(ghz4, eftjhu, 0.) \
  COUPLING_COMMAND(ghzgs2, eftjhu, 0.) \
  COUPLING_COMMAND(ghzgs4, eftjhu, 0.) \
  COUPLING_COMMAND(ghgsgs2, eftjhu, 0.) \
  COUPLING_COMMAND(ghgsgs4, eftjhu, 0.) \
  COUPLING_COMMAND(ghg2, eftjhu, 0.) \
  COUPLING_COMMAND(ghg4, eftjhu, 0.) \

#define HIGGSBASIS_COUPLING_COMMANDS \
  COUPLING_COMMAND(dCz, hbasis, 0.) \
  COUPLING_COMMAND(Czz, hbasis, 0.) \
  COUPLING_COMMAND(Czbx, hbasis, 0.) \
  COUPLING_COMMAND(tCzz, hbasis, 0.) \
  COUPLING_COMMAND(dCw, hbasis, 0.) \
  COUPLING_COMMAND(Cww, hbasis, 0.) \
  COUPLING_COMMAND(Cwbx, hbasis, 0.) \
  COUPLING_COMMAND(tCww, hbasis, 0.) \
  COUPLING_COMMAND(Cza, hbasis, 0.) \
  COUPLING_COMMAND(tCza, hbasis, 0.) \
  COUPLING_COMMAND(Cabx, hbasis, 0.) \
  COUPLING_COMMAND(Caa, hbasis, 0.) \
  COUPLING_COMMAND(tCaa, hbasis, 0.) \
  COUPLING_COMMAND(Cgg, hbasis, 0.) \
  COUPLING_COMMAND(tCgg, hbasis, 0.) \

#define EFT_HIGGSBASIS_COUPLING_COMMANDS \
  COUPLING_COMMAND(dCz, efthbasis, 0.) \
  COUPLING_COMMAND(Czbx, efthbasis, 0.) \
  COUPLING_COMMAND(Czz, efthbasis, 0.) \
  COUPLING_COMMAND(tCzz, efthbasis, 0.) \
  COUPLING_COMMAND(Cza, efthbasis, 0.) \
  COUPLING_COMMAND(tCza, efthbasis, 0.) \
  COUPLING_COMMAND(Caa, efthbasis, 0.) \
  COUPLING_COMMAND(tCaa, efthbasis, 0.) \
  COUPLING_COMMAND(Cgg, efthbasis, 0.) \
  COUPLING_COMMAND(tCgg, efthbasis, 0.) \

//This may be bad practice but here we define basis that use triple gauge or quartic gauge couplings for behind the scenes usage
// The first one is utilized is a user inputs eft_jhugen and want the eft_jhugen + the triple and quartic couplings returned back 

#define EFT_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS \
  COUPLING_COMMAND(ghz1, eftjhutrip, 0.) \
  COUPLING_COMMAND(ghz1_prime2, eftjhutrip, 0.) \
  COUPLING_COMMAND(ghz2, eftjhutrip, 0.) \
  COUPLING_COMMAND(ghz4, eftjhutrip, 0.) \
  COUPLING_COMMAND(ghzgs2, eftjhutrip, 0.) \
  COUPLING_COMMAND(ghzgs4, eftjhutrip, 0.) \
  COUPLING_COMMAND(ghgsgs2, eftjhutrip, 0.) \
  COUPLING_COMMAND(ghgsgs4, eftjhutrip, 0.) \
  COUPLING_COMMAND(ghg2, eftjhutrip, 0.) \
  COUPLING_COMMAND(ghg4, eftjhutrip, 0.) \
  COUPLING_COMMAND(dV_Z, eftjhutrip, 0.) \
  COUPLING_COMMAND(dV_A, eftjhutrip, 0.) \
  COUPLING_COMMAND(dP_Z, eftjhutrip, 0.) \
  COUPLING_COMMAND(dP_A, eftjhutrip, 0.) \
  COUPLING_COMMAND(dM_Z, eftjhutrip, 0.) \
  COUPLING_COMMAND(dM_A, eftjhutrip, 0.) \
  COUPLING_COMMAND(dFour_Z, eftjhutrip, 0.) \
  COUPLING_COMMAND(dFour_A, eftjhutrip, 0.) \
  COUPLING_COMMAND(dZZWpWm, eftjhutrip, 0.) \
  COUPLING_COMMAND(dZAWpWm, eftjhutrip, 0.) \
  COUPLING_COMMAND(dAAWpWm, eftjhutrip, 0.) \

// This new basis is for when eft_jhugen is input and the user only wants the jhugen amplitude couplings without quartic and gauge contributions

#define AMPLITUDE_JHUGEN_COUPLING_COMMANDS \
  COUPLING_COMMAND(ghz1, ampjhu, 0.) \
  COUPLING_COMMAND(ghz1_prime2, ampjhu, 0.) \
  COUPLING_COMMAND(ghz2, ampjhu, 0.) \
  COUPLING_COMMAND(ghz4, ampjhu, 0.) \
  COUPLING_COMMAND(ghw1, ampjhu, 0.) \
  COUPLING_COMMAND(ghw1_prime2, ampjhu, 0.) \
  COUPLING_COMMAND(ghw2, ampjhu, 0.) \
  COUPLING_COMMAND(ghw4, ampjhu, 0.) \
  COUPLING_COMMAND(ghzgs1_prime2, ampjhu, 0.) \
  COUPLING_COMMAND(ghzgs2, ampjhu, 0.) \
  COUPLING_COMMAND(ghzgs4, ampjhu, 0.) \
  COUPLING_COMMAND(ghgsgs2, ampjhu, 0.) \
  COUPLING_COMMAND(ghgsgs4, ampjhu, 0.) \
  COUPLING_COMMAND(ghg2, ampjhu, 0.) \
  COUPLING_COMMAND(ghg4, ampjhu, 0.) \

// This basis is for when user wants the output in higgs basis couplings with the triple gauge couplings included 
#define HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS \
  COUPLING_COMMAND(dCz, hbasistrip, 0.) \
  COUPLING_COMMAND(Czz, hbasistrip, 0.) \
  COUPLING_COMMAND(Czbx, hbasistrip, 0.) \
  COUPLING_COMMAND(tCzz, hbasistrip, 0.) \
  COUPLING_COMMAND(dCw, hbasistrip, 0.) \
  COUPLING_COMMAND(Cww, hbasistrip, 0.) \
  COUPLING_COMMAND(Cwbx, hbasistrip, 0.) \
  COUPLING_COMMAND(tCww, hbasistrip, 0.) \
  COUPLING_COMMAND(Cza, hbasistrip, 0.) \
  COUPLING_COMMAND(tCza, hbasistrip, 0.) \
  COUPLING_COMMAND(Cabx, hbasistrip, 0.) \
  COUPLING_COMMAND(Caa, hbasistrip, 0.) \
  COUPLING_COMMAND(tCaa, hbasistrip, 0.) \
  COUPLING_COMMAND(Cgg, hbasistrip, 0.) \
  COUPLING_COMMAND(tCgg, hbasistrip, 0.) \
  COUPLING_COMMAND(dKa, hbasistrip, 0.) \
  COUPLING_COMMAND(tKa, hbasistrip, 0.) \
  COUPLING_COMMAND(dKz, hbasistrip, 0.) \
  COUPLING_COMMAND(tKz, hbasistrip, 0.) \
  COUPLING_COMMAND(dg1z, hbasistrip, 0.) \

// This basis is for when the user input eft hbasis and wants the eft hbasis plus the new triple and quartic couplings
// This output is mixed by definition

#define EFT_HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS \
  COUPLING_COMMAND(dCz, efthbasistrip, 0.) \
  COUPLING_COMMAND(Czbx, efthbasistrip, 0.) \
  COUPLING_COMMAND(Czz, efthbasistrip, 0.) \
  COUPLING_COMMAND(tCzz, efthbasistrip, 0.) \
  COUPLING_COMMAND(Cza, efthbasistrip, 0.) \
  COUPLING_COMMAND(tCza, efthbasistrip, 0.) \
  COUPLING_COMMAND(Caa, efthbasistrip, 0.) \
  COUPLING_COMMAND(tCaa, efthbasistrip, 0.) \
  COUPLING_COMMAND(Cgg, efthbasistrip, 0.) \
  COUPLING_COMMAND(tCgg, efthbasistrip, 0.) \
  COUPLING_COMMAND(dKa, efthbasistrip, 0.) \
  COUPLING_COMMAND(tKa, efthbasistrip, 0.) \
  COUPLING_COMMAND(dKz, efthbasistrip, 0.) \
  COUPLING_COMMAND(tKz, efthbasistrip, 0.) \
  COUPLING_COMMAND(dg1z, efthbasistrip, 0.) \

// Internal coupling basis for Neutral current and charged current interactions 
// User must specify couplings for neutral or charged current interactions
// The default will be returning all of the couplings One must turn off charged current or neutral current manually
// This option is only useful if you have output as amp_jhugen for now

// Warsaw Basis Couplings Using cT instead of cHBx as in Rosetta etc
#define WARSAWBASIS_COUPLING_COMMANDS \
  COUPLING_COMMAND(cHbx, warsaw, 0.) \
  COUPLING_COMMAND(cHD, warsaw, 0.) \
  COUPLING_COMMAND(cHG, warsaw, 0.) \
  COUPLING_COMMAND(cHW, warsaw, 0.) \
  COUPLING_COMMAND(cHB, warsaw, 0.) \
  COUPLING_COMMAND(cHWB, warsaw, 0.) \
  COUPLING_COMMAND(tcHG, warsaw, 0.) \
  COUPLING_COMMAND(tcHW, warsaw, 0.) \
  COUPLING_COMMAND(tcHB, warsaw, 0.) \
  COUPLING_COMMAND(tcHWB, warsaw, 0.) \

namespace JHUGenLexiconCouplings{
#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) coupl_##PREFIX##_##NAME,

  enum Amplitude_JHUGen_CouplingType{
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS
    nAmplitude_JHUGen_CouplingTypes
  };
  enum Amplitude_JHUGen_Include_Triple_CouplingType{
    AMPLITUDE_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS
    nAmplitude_JHUGen_Include_Triple_CouplingTypes
  };
  enum EFT_JHUGen_CouplingType{
    EFT_JHUGEN_COUPLING_COMMANDS
    nEFT_JHUGen_CouplingTypes
  };
  enum EFT_JHUGen_Include_Triple_CouplingType{
    EFT_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS
    nEFT_JHUGen_Include_Triple_CouplingTypes
  };
  enum HiggsBasis_CouplingType{
    HIGGSBASIS_COUPLING_COMMANDS
    nHiggsBasis_CouplingTypes
  };
  enum HiggsBasis_Include_Triple_CouplingType {
    HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS
    nHiggsBasis_Include_Triple_CouplingTypes
  };
  enum EFT_HiggsBasis_CouplingType{
    EFT_HIGGSBASIS_COUPLING_COMMANDS
    nEFT_HiggsBasis_CouplingTypes
  };
  enum EFT_HiggsBasis_Include_Triple_CouplingType {
    EFT_HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS
    nEFT_HiggsBasis_Include_Triple_CouplingTypes
  };
  enum WarsawBasis_CouplingType{
    WARSAWBASIS_COUPLING_COMMANDS
    nWarsawBasis_CouplingTypes
  };

#undef COUPLING_COMMAND

  // Functions to get the coupling names from the indices
  std::string getCouplingName(JHUGenLexiconCouplings::Amplitude_JHUGen_CouplingType type);
  std::string getCouplingName(JHUGenLexiconCouplings::Amplitude_JHUGen_Include_Triple_CouplingType type);
  std::string getCouplingName(JHUGenLexiconCouplings::EFT_JHUGen_CouplingType type);
  std::string getCouplingName(JHUGenLexiconCouplings::EFT_JHUGen_Include_Triple_CouplingType type);
  std::string getCouplingName(JHUGenLexiconCouplings::HiggsBasis_CouplingType type);
  std::string getCouplingName(JHUGenLexiconCouplings::HiggsBasis_Include_Triple_CouplingType type);
  std::string getCouplingName(JHUGenLexiconCouplings::EFT_HiggsBasis_CouplingType type);
  std::string getCouplingName(JHUGenLexiconCouplings::EFT_HiggsBasis_Include_Triple_CouplingType type);
  std::string getCouplingName(JHUGenLexiconCouplings::WarsawBasis_CouplingType type);
}


#endif
