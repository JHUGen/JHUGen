#ifndef LEXICONDM_IOHELPERS_H
#define LEXICONDM_IOHELPERS_H

#include <string>

namespace LexiConDMIOHelpers{
  enum IOBasisType{
    bAmplitude_JHUGen = 0,

    bEFT_JHUGen,
    bEFT_HiggsBasis,

    nIOBases
  };

  LexiConDMIOHelpers::IOBasisType getIOBasisFromString(std::string const&);

}

#endif
