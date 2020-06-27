#ifndef JHUGENLEXICON_TRANSLATOR_H
#define JHUGENLEXICON_TRANSLATOR_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include "JHUGenLexiconOptionParser.h"


class JHUGenLexiconTranslator{
protected:
  JHUGenLexiconOptionParser const& opts;

  std::unordered_map<std::string, double> result_parameters;
  std::unordered_map<std::string, std::pair<double, double> > result_couplings;

  std::vector<std::vector<double>> getTranslationMatrix(
    JHUGenLexiconIOHelpers::IOBasisType const& basis_input, JHUGenLexiconIOHelpers::IOBasisType const& basis_output,
    std::unordered_map<std::string, bool> const& input_flags,
    std::unordered_map<std::string, double> const& input_parameters
  ) const;

  std::vector< std::pair<double, double> > getOrderedInputCouplings(
    JHUGenLexiconIOHelpers::IOBasisType const& basis_input,
    std::unordered_map<std::string, bool> const& input_flags,
    std::unordered_map<std::string, double> const& input_parameters,
    std::unordered_map<std::string, std::pair<double, double> > const& input_couplings
  ) const;

  void interpretOutputCouplings(
    JHUGenLexiconIOHelpers::IOBasisType const& basis_output,
    std::unordered_map<std::string, bool> const& input_flags,
    std::unordered_map<std::string, double> const& input_parameters,
    std::vector< std::pair<double, double> >& output_vector
  );

  void translate();

public:
  JHUGenLexiconTranslator(JHUGenLexiconOptionParser const& opts_);
  ~JHUGenLexiconTranslator(){}

  std::unordered_map<std::string, double>& getParameters(){ return result_parameters; }
  std::unordered_map<std::string, double> const& getParameters() const{ return result_parameters; }

  std::unordered_map<std::string, std::pair<double, double> >& getCouplings(){ return result_couplings; }
  std::unordered_map<std::string, std::pair<double, double> > const& getCouplings() const{ return result_couplings; }

};

#endif
