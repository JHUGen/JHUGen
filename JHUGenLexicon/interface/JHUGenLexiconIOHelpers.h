#ifndef JHUGENLEXICON_IOHELPERS_H
#define JHUGENLEXICON_IOHELPERS_H

#include <string>

namespace JHUGenLexiconIOHelpers{
  enum IOBasisType{
    bAmplitude_JHUGen = 0,

    bEFT_JHUGen,
    bEFT_HiggsBasis,
    bHiggsBasis,
    bWarsawBasis,

    nIOBases
  };

  JHUGenLexiconIOHelpers::IOBasisType getIOBasisFromString(std::string const&);

}

#endif
