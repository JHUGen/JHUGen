#ifndef JHUGENLEXICON_OPTIONPARSER_H
#define JHUGENLEXICON_OPTIONPARSER_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include "JHUGenLexiconIOHelpers.h"


class JHUGenLexiconOptionParser{
public:
  static bool globalHelpFlag;

protected:
  std::vector<std::string> rawOptions;

  JHUGenLexiconIOHelpers::IOBasisType basis_input;
  JHUGenLexiconIOHelpers::IOBasisType basis_output;

  std::unordered_map<std::string, bool> flags;
  std::unordered_map<std::string, double> parameters;
  std::unordered_map<std::string, std::pair<double, double> > couplings;

public:
  JHUGenLexiconOptionParser(int argc, char** argv);
  ~JHUGenLexiconOptionParser(){}

  void analyze();
  void interpretOption(std::string const& wish, std::string const& value, bool& invalidOption);
  void printOptionsHelp(bool command_fail) const;

  JHUGenLexiconIOHelpers::IOBasisType const& getInputBasis() const{ return basis_input; }
  JHUGenLexiconIOHelpers::IOBasisType const& getOutputBasis() const{ return basis_output; }

  std::unordered_map<std::string, bool>& getInputFlags(){ return flags; }
  std::unordered_map<std::string, bool> const& getInputFlags() const{ return flags; }

  std::unordered_map<std::string, double>& getInputParameters(){ return parameters; }
  std::unordered_map<std::string, double> const& getInputParameters() const{ return parameters; }

  std::unordered_map<std::string, std::pair<double, double> >& getInputCouplings(){ return couplings; }
  std::unordered_map<std::string, std::pair<double, double> > const& getInputCouplings() const{ return couplings; }

};

#endif
