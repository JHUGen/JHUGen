#ifndef LEXICONDM_OPTIONPARSER_H
#define LEXICONDM_OPTIONPARSER_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include "LexiConDMIOHelpers.h"


class LexiConDMOptionParser{
public:
  static bool globalHelpFlag;

protected:
  std::vector<std::string> rawOptions;

  LexiConDMIOHelpers::IOBasisType basis_input;
  LexiConDMIOHelpers::IOBasisType basis_output;

  std::unordered_map<std::string, bool> flags;
  std::unordered_map<std::string, double> parameters;
  std::unordered_map<std::string, std::pair<double, double> > couplings;

public:
  LexiConDMOptionParser(int argc, char** argv);
  ~LexiConDMOptionParser(){}

  void analyze();
  void interpretOption(std::string const& wish, std::string const& value, bool& invalidOption);
  void printOptionsHelp(bool command_fail) const;

  LexiConDMIOHelpers::IOBasisType const& getInputBasis() const{ return basis_input; }
  LexiConDMIOHelpers::IOBasisType const& getOutputBasis() const{ return basis_output; }

  std::unordered_map<std::string, bool>& getInputFlags(){ return flags; }
  std::unordered_map<std::string, bool> const& getInputFlags() const{ return flags; }

  std::unordered_map<std::string, double>& getInputParameters(){ return parameters; }
  std::unordered_map<std::string, double> const& getInputParameters() const{ return parameters; }

  std::unordered_map<std::string, std::pair<double, double> >& getInputCouplings(){ return couplings; }
  std::unordered_map<std::string, std::pair<double, double> > const& getInputCouplings() const{ return couplings; }

};

#endif
