#include <cassert>
#include <exception>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "LexiConDMCouplings.h"
#include "LexiConDMTranslator.h"
#include "LexiConDMHelperFunctions.h"


using namespace std;
using namespace LexiConDMHelperFunctions;
using namespace LexiConDMIOHelpers;
using namespace LexiConDMCouplings;


LexiConDMTranslator::LexiConDMTranslator(LexiConDMOptionParser const& opts_) :
  opts(opts_)
{
  translate();
}
void LexiConDMTranslator::translate(){
  auto const& input_flags = opts.getInputFlags();
  auto const& input_parameters = opts.getInputParameters();
  auto const& input_couplings = opts.getInputCouplings();
  auto const& basis_input = opts.getInputBasis();
  auto const& basis_output = opts.getOutputBasis();

  // Get the translation matrix
  std::vector<std::vector<double>> tmatrix = getTranslationMatrix(basis_input, basis_output, input_parameters);
  // Get the input vector
  std::vector< std::pair<double, double> > vinput = getOrderedInputCouplings(basis_input, input_flags, input_parameters, input_couplings);
  // Assign the output vector
  std::vector< std::pair<double, double> > voutput(tmatrix.size(), std::pair<double, double>(0, 0));

  // Do the matrix multiplication
  for (size_t i=0; i<tmatrix.size(); i++){
    auto const& trow = tmatrix.at(i);
    for (size_t j=0; j<trow.size(); j++){
      voutput.at(j).first += trow.at(j) * vinput.at(j).first;
      voutput.at(j).second += trow.at(j) * vinput.at(j).second;
    }
  }

  // Set the results
  interpretOutputCouplings(basis_output, input_flags, input_parameters, voutput);
}

std::vector<std::vector<double>> LexiConDMTranslator::getTranslationMatrix(
  LexiConDMIOHelpers::IOBasisType const& basis_input, LexiConDMIOHelpers::IOBasisType const& basis_output,
  std::unordered_map<std::string, double> const& input_parameters
) const{
  double sw; getValueWithDefault<std::string, double>(input_parameters, "sin2ThetaW", sw, DEFVAL_SW);
  std::vector<std::vector<double>> res;

  // This is where the translation matrixes should be coded carefully.
  // See also the getOrderedInputCouplings and interpretOutputCouplings functions for what the input/output vectors expect for JHUGen conventions (e.g. ghz1_prime2 scaled already by MZ^2/L1ZZ^2, so no need to put that for example).
  if (basis_input == bAmplitude_JHUGen){
    if (basis_output == bAmplitude_JHUGen){
      res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));
      for (size_t i=0; i<(size_t) nAmplitude_JHUGen_CouplingTypes; i++) res.at(i).at(i)=1;
    }
    else if (basis_output == bEFT_JHUGen){
      res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));
      /*
      res = std::vector<std::vector<double>>{
        {  },
        {  }
      };
      */
    }
    else if (basis_output == bEFT_HiggsBasis){
      res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nAmplitude_JHUGen_CouplingTypes, 0));
      /*
      res = std::vector<std::vector<double>>{
        {  },
        {  }
      };
      */
    }
  }
  else if (basis_input == bEFT_JHUGen){
    if (basis_output == bAmplitude_JHUGen){
      res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
      /*
      res = std::vector<std::vector<double>>{
        {  },
        {  }
      };
      */
    }
    else if (basis_output == bEFT_JHUGen){
      res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
      for (size_t i=0; i<(size_t) nEFT_JHUGen_CouplingTypes; i++) res.at(i).at(i)=1;
    }
    else if (basis_output == bEFT_HiggsBasis){
      res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nEFT_JHUGen_CouplingTypes, 0));
      /*
      res = std::vector<std::vector<double>>{
        {  },
        {  }
      };
      */
    }
  }
  else if (basis_input == bEFT_HiggsBasis){
    if (basis_output == bAmplitude_JHUGen){
      res.assign(nAmplitude_JHUGen_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
      /*
      res = std::vector<std::vector<double>>{
        {  },
        {  }
      };
      */
    }
    else if (basis_output == bEFT_JHUGen){
      res.assign(nEFT_JHUGen_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
      /*
      res = std::vector<std::vector<double>>{
        {  },
        {  }
      };
      */
    }
    else if (basis_output == bEFT_HiggsBasis){
      res.assign(nEFT_HiggsBasis_CouplingTypes, std::vector<double>(nEFT_HiggsBasis_CouplingTypes, 0));
      for (size_t i=0; i<(size_t) nEFT_HiggsBasis_CouplingTypes; i++) res.at(i).at(i)=1;
    }
  }

  if (res.empty()){
    cerr << "LexiConDMTranslator::getTranslationMatrix: Translation from input basis " << basis_input << " to output basis " << basis_output << " is not implemented." << endl;
    assert(0);
  }

  return res;
}

std::vector< std::pair<double, double> > LexiConDMTranslator::getOrderedInputCouplings(
  LexiConDMIOHelpers::IOBasisType const& basis_input,
  std::unordered_map<std::string, bool> const& input_flags,
  std::unordered_map<std::string, double> const& input_parameters,
  std::unordered_map<std::string, std::pair<double, double> > const& input_couplings
) const{
  bool useMCFMAtInput; getValueWithDefault<std::string, bool>(input_flags, "useMCFMAtInput", useMCFMAtInput, false);
  bool distinguish_HWWcouplings; getValueWithDefault<std::string, bool>(input_flags, "distinguish_HWWcouplings", distinguish_HWWcouplings, false);
  double MZ; getValueWithDefault<std::string, double>(input_parameters, "MZ", MZ, DEFVAL_MZ);
  double MW; getValueWithDefault<std::string, double>(input_parameters, "MW", MW, DEFVAL_MW);
  double Lambda_z1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_z1", Lambda_z1, DEFVAL_LAMBDA_VI);
  double Lambda_w1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_w1", Lambda_w1, DEFVAL_LAMBDA_VI);
  double Lambda_zgs1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_zgs1", Lambda_zgs1, DEFVAL_LAMBDA_VI);

  std::vector< std::pair<double, double> > res;
#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) \
  getValueWithDefault<std::string, std::pair<double, double>>(input_couplings, #NAME, res.at(coupl_##PREFIX##_##NAME), std::pair<double, double>(DEFVAL, 0)); \
  if (useMCFMAtInput && std::string(#NAME).find("ghz")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= 2.; res.at(coupl_##PREFIX##_##NAME).second *= 2.; } \
  if (std::string(#NAME).find("ghz")!=std::string::npos && std::string(#NAME).find("ghzgs")==std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= std::pow(MZ/Lambda_z1, 2); res.at(coupl_##PREFIX##_##NAME).second *= std::pow(MZ/Lambda_z1, 2); } \
  else if (std::string(#NAME).find("ghzgs")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= std::pow(MZ/Lambda_zgs1, 2); res.at(coupl_##PREFIX##_##NAME).second *= std::pow(MZ/Lambda_zgs1, 2); } \
  else if (std::string(#NAME).find("ghw")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ res.at(coupl_##PREFIX##_##NAME).first *= std::pow(MW/Lambda_w1, 2); res.at(coupl_##PREFIX##_##NAME).second *= std::pow(MW/Lambda_w1, 2); }

  switch (basis_input){
  case bAmplitude_JHUGen:
  {
    res.assign(nAmplitude_JHUGen_CouplingTypes, std::pair<double, double>(0, 0));
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS;

    if (!distinguish_HWWcouplings){
      static_assert((unsigned int) coupl_ampjhu_ghzgs1_prime2 - (unsigned int) coupl_ampjhu_ghw1 == (unsigned int) coupl_ampjhu_ghw1 - (unsigned int) coupl_ampjhu_ghz1);
      for (unsigned int i=(unsigned int) coupl_ampjhu_ghz1; i<(unsigned int) coupl_ampjhu_ghw1; i++){
        res.at(i+(unsigned int) coupl_ampjhu_ghw1 - (unsigned int) coupl_ampjhu_ghz1).first = res.at(i).first;
        res.at(i+(unsigned int) coupl_ampjhu_ghw1 - (unsigned int) coupl_ampjhu_ghz1).second = res.at(i).second;
      }
    }
    break;
  }
  case bEFT_JHUGen:
  {
    res.assign(nEFT_JHUGen_CouplingTypes, std::pair<double, double>(0, 0));
    EFT_JHUGEN_COUPLING_COMMANDS;
    break;
  }
  case bEFT_HiggsBasis:
  {
    res.assign(nEFT_HiggsBasis_CouplingTypes, std::pair<double, double>(0, 0));
    EFT_HIGGSBASIS_COUPLING_COMMANDS;
    break;
  }
  default:
    cerr << "LexiConDMTranslator::getOrderedInputCouplings: Input basis " << basis_input << " is not implemented." << endl;
    assert(0);
  }

#undef COUPLING_COMMAND

  return res;
}

void LexiConDMTranslator::interpretOutputCouplings(
  LexiConDMIOHelpers::IOBasisType const& basis_output,
  std::unordered_map<std::string, bool> const& input_flags,
  std::unordered_map<std::string, double> const& input_parameters,
  std::vector< std::pair<double, double> >& output_vector
){
  bool useMCFMAtOutput; getValueWithDefault<std::string, bool>(input_flags, "useMCFMAtOutput", useMCFMAtOutput, false);

  double MZ; getValueWithDefault<std::string, double>(input_parameters, "MZ", MZ, DEFVAL_MZ);
  double MW; getValueWithDefault<std::string, double>(input_parameters, "MW", MW, DEFVAL_MW);
  double Lambda_z1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_z1", Lambda_z1, DEFVAL_LAMBDA_VI);
  double Lambda_w1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_w1", Lambda_w1, DEFVAL_LAMBDA_VI);
  double Lambda_zgs1; getValueWithDefault<std::string, double>(input_parameters, "Lambda_zgs1", Lambda_zgs1, DEFVAL_LAMBDA_VI);

#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) \
  if (useMCFMAtOutput && std::string(#NAME).find("ghz")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= 2.; output_vector.at(coupl_##PREFIX##_##NAME).second /= 2.; } \
  if (std::string(#NAME).find("ghz")!=std::string::npos && std::string(#NAME).find("ghzgs")==std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= std::pow(MZ/Lambda_z1, 2); output_vector.at(coupl_##PREFIX##_##NAME).second /= std::pow(MZ/Lambda_z1, 2); } \
  else if (std::string(#NAME).find("ghzgs")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= std::pow(MZ/Lambda_zgs1, 2); output_vector.at(coupl_##PREFIX##_##NAME).second /= std::pow(MZ/Lambda_zgs1, 2); } \
  else if (std::string(#NAME).find("ghw")!=std::string::npos && std::string(#NAME).find("prime2")!=std::string::npos){ output_vector.at(coupl_##PREFIX##_##NAME).first /= std::pow(MW/Lambda_w1, 2); output_vector.at(coupl_##PREFIX##_##NAME).second /= std::pow(MW/Lambda_w1, 2); } \
  result_couplings[#NAME] = output_vector.at(coupl_##PREFIX##_##NAME);

  switch (basis_output){
  case bAmplitude_JHUGen:
  {
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS;
    break;
  }
  case bEFT_JHUGen:
  {
    EFT_JHUGEN_COUPLING_COMMANDS;
    break;
  }
  case bEFT_HiggsBasis:
  {
    EFT_HIGGSBASIS_COUPLING_COMMANDS;
    break;
  }
  default:
    cerr << "LexiConDMTranslator::interpretOutputCouplings: Output basis " << basis_output << " is not implemented." << endl;
    assert(0);
  }

#undef COUPLING_COMMAND
}
