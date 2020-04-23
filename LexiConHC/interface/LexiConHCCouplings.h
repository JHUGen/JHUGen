#ifndef LEXICONHC_COUPLINGS_H
#define LEXICONHC_COUPLINGS_H


#include <string>

// Default values of parameters as in mod_Parameters.F90
#define DEFVAL_MZ 91.1876
#define DEFVAL_MW 80.399
#define DEFVAL_SW 0.23119
#define DEFVAL_LAMBDA_VI 10000.

// Coupling definitions
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
  COUPLING_COMMAND(dV_Z, ampjhu, 0.) \
  COUPLING_COMMAND(dV_A, ampjhu, 0.) \
  COUPLING_COMMAND(dP_Z, ampjhu, 0.) \
  COUPLING_COMMAND(dP_A, ampjhu, 0.) \
  COUPLING_COMMAND(dM_Z, ampjhu, 0.) \
  COUPLING_COMMAND(dM_A, ampjhu, 0.) \
  COUPLING_COMMAND(dFour_Z, ampjhu, 0.) \
  COUPLING_COMMAND(dFour_A, ampjhu, 0.) \
  COUPLING_COMMAND(dZZWpWm, ampjhu, 0.) \
  COUPLING_COMMAND(dZAWpWm, ampjhu, 0.) \
  COUPLING_COMMAND(dAAWpWm, ampjhu, 0.)

#define EFT_JHUGEN_COUPLING_COMMANDS \
  COUPLING_COMMAND(ghz1, eftjhu, 0.) \
  COUPLING_COMMAND(ghz1_prime2, eftjhu, 0.) \
  COUPLING_COMMAND(ghz2, eftjhu, 0.) \
  COUPLING_COMMAND(ghz4, eftjhu, 0.) \
  COUPLING_COMMAND(ghzgs2, eftjhu, 0.) \
  COUPLING_COMMAND(ghzgs4, eftjhu, 0.) \
  COUPLING_COMMAND(ghgsgs2, eftjhu, 0.) \
  COUPLING_COMMAND(ghgsgs4, eftjhu, 0.)

#define EFT_HIGGSBASIS_COUPLING_COMMANDS \
  COUPLING_COMMAND(d_cz, efthbasis, 0.) \
  COUPLING_COMMAND(czz, efthbasis, 0.) \
  COUPLING_COMMAND(cz_box, efthbasis, 0.) \
  COUPLING_COMMAND(czz_tilde, efthbasis, 0.) \
  COUPLING_COMMAND(d_cw, efthbasis, 0.) \
  COUPLING_COMMAND(cww, efthbasis, 0.) \
  COUPLING_COMMAND(cw_box, efthbasis, 0.) \
  COUPLING_COMMAND(cww_tilde, efthbasis, 0.) \
  COUPLING_COMMAND(czgs, efthbasis, 0.) \
  COUPLING_COMMAND(czgs_tilde, efthbasis, 0.) \
  COUPLING_COMMAND(cgs_box, efthbasis, 0.) \
  COUPLING_COMMAND(cgsgs, efthbasis, 0.) \
  COUPLING_COMMAND(cgsgs_tilde, efthbasis, 0.)


namespace LexiConHCCouplings{
#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) coupl_##PREFIX##_##NAME,

  enum Amplitude_JHUGen_CouplingType{
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS
    nAmplitude_JHUGen_CouplingTypes
  };
  enum EFT_JHUGen_CouplingType{
    EFT_JHUGEN_COUPLING_COMMANDS
    nEFT_JHUGen_CouplingTypes
  };
  enum EFT_HiggsBasis_CouplingType{
    EFT_HIGGSBASIS_COUPLING_COMMANDS
    nEFT_HiggsBasis_CouplingTypes
  };

#undef COUPLING_COMMAND

  // Functions to get the coupling names from the indices
  std::string getCouplingName(LexiConHCCouplings::Amplitude_JHUGen_CouplingType type);
  std::string getCouplingName(LexiConHCCouplings::EFT_JHUGen_CouplingType type);
  std::string getCouplingName(LexiConHCCouplings::EFT_HiggsBasis_CouplingType type);

}


#endif
