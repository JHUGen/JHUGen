#ifndef LEXICONHC_IOHELPERS_H
#define LEXICONHC_IOHELPERS_H

#include <string>

namespace LexiConHCIOHelpers{
  enum IOBasisType{
    bAmplitude_JHUGen = 0,

    bEFT_JHUGen,
    bEFT_HiggsBasis,

    nIOBases
  };

  LexiConHCIOHelpers::IOBasisType getIOBasisFromString(std::string const&);

}

#endif
