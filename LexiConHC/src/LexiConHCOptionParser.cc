#include <exception>
#include <cstdlib>
#include <iostream>
#include "LexiConHCCouplings.h"
#include "LexiConHCOptionParser.h"
#include "LexiConHCHelperFunctions.h"


using namespace std;
using namespace LexiConHCHelperFunctions;
using namespace LexiConHCIOHelpers;


bool LexiConHCOptionParser::globalHelpFlag = false;

LexiConHCOptionParser::LexiConHCOptionParser(int argc, char** argv) :
  basis_input(nIOBases),
  basis_output(nIOBases)
{
  for (int a=0; a<argc; a++){
    string tmpArg(argv[a]);
    rawOptions.push_back(tmpArg);
  }

  analyze();
}
void LexiConHCOptionParser::analyze(){
  bool hasInvalidOption=false;
  char rawdelimiter = '=';
  for (unsigned int opt=1; opt<rawOptions.size(); opt++){
    string wish, value;
    splitOption(rawOptions.at(opt), wish, value, rawdelimiter);

    bool isInvalidOption=false;
    interpretOption(wish, value, isInvalidOption);
    hasInvalidOption |= isInvalidOption;
  }

  if (LexiConHCOptionParser::globalHelpFlag){ printOptionsHelp(false); return; }

  if (basis_input == nIOBases || basis_output == nIOBases){ cerr << "LexiConHCOptionParser::analyze: You have to specify the input and output basis types." << endl; if (!hasInvalidOption) hasInvalidOption=true; }

  bool distinguish_HWWcouplings; getValueWithDefault<std::string, bool>(flags, "distinguish_HWWcouplings", distinguish_HWWcouplings, false);
  if (!distinguish_HWWcouplings && basis_input == bAmplitude_JHUGen){
    using namespace LexiConHCCouplings;
    static_assert((unsigned int) coupl_ampjhu_ghzgs1_prime2 - (unsigned int) coupl_ampjhu_ghw1 == (unsigned int) coupl_ampjhu_ghw1 - (unsigned int) coupl_ampjhu_ghz1);
    for (unsigned int i=(unsigned int) coupl_ampjhu_ghz1; i<(unsigned int) coupl_ampjhu_ghw1; i++){
      couplings[getCouplingName((Amplitude_JHUGen_CouplingType) (i+(unsigned int) coupl_ampjhu_ghw1 - (unsigned int) coupl_ampjhu_ghz1))] =
        couplings[getCouplingName((Amplitude_JHUGen_CouplingType) i)];
    }
  }


  // Print help if needed and abort at this point, nowhere later
  if (hasInvalidOption) printOptionsHelp(hasInvalidOption);
}

void LexiConHCOptionParser::interpretOption(std::string const& wish, std::string const& value, bool& invalidOption){
  invalidOption=false;
  if (wish.empty()){
    if (value=="help") LexiConHCOptionParser::globalHelpFlag = true;
    else{
      cerr << "LexiConHCOptionParser::interpretOption: Unknown unspecified argument: " << value << endl;
      invalidOption=true;
    }
  }

  // Input and output bases
  else if (wish=="input_basis") basis_input = getIOBasisFromString(value);
  else if (wish=="output_basis") basis_output = getIOBasisFromString(value);

  // Flags come first because the formats of parameters and flags are the same
  else if (wish=="useMCFMAtInput"){ flags[wish] = false; castStringToValue(value, flags[wish]); }
  else if (wish=="useMCFMAtOutput"){ flags[wish] = false; castStringToValue(value, flags[wish]); }
  else if (wish=="distinguish_HWWcouplings"){ flags[wish] = false; castStringToValue(value, flags[wish]); }

  // Parameters come next, check if a comma is not found to assign them
  else if (value.find(",")==std::string::npos) parameters[wish] = stod(value);

  // Couplings last
  else if (value.find(",")!=std::string::npos){
    string vRe, vIm;
    splitOption(value, vRe, vIm, ',');
    couplings[wish] = std::pair<double, double>(stod(vRe), stod(vIm));
  }


  else cerr << "Unknown specified argument: " << value << " with specifier " << wish << endl;
}

void LexiConHCOptionParser::printOptionsHelp(bool command_fail)const{
  cout << endl;
  cout << "The options implemented in LexiConHC (format: specifier=value):\n\n";

  cout << "- input_basis: Input coupling conventions\n\n";
  cout << "- output_basis: Output coupling conventions\n\n";

  cout << "- useMCFMAtInput: Use MCFM conventions in the input JHUGen couplings. Assumes the ghv* couplings (ghv1=2 in SM) are divided by 2.\n\n";
  cout << "- useMCFMAtOutput: Use MCFM conventions in the output JHUGen couplings. Divides the ghv* couplings by 2.\n\n";
  cout << "- distinguish_HWWcouplings: Distinguish HZZ and HWW couplings in the JHUGen amplitude basis if it is the input. Default is false.\n\n";

  cout << "- The format to set any parameter is [specifier]=[value].\n\n";
  cout << "- The format to set any coupling [specifier] to the complex number ([vRe], [vIm]) is [specifier]=[vRe],[vIm].\n\n";
#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) cout << "\t" << #NAME << " (default = " << DEFVAL << ",0)\n";
  cout << "- Allowed couplings for the JHUGen amplitude formalism:\n";
  AMPLITUDE_JHUGEN_COUPLING_COMMANDS;
  cout << "- Allowed couplings for the JHUGen EFT formalism:\n";
  EFT_JHUGEN_COUPLING_COMMANDS;
  cout << "- Allowed couplings for the Higgs basis EFT formalism:\n";
  EFT_HIGGSBASIS_COUPLING_COMMANDS;
#undef COUPLING_COMMAND
  cout << "- Allowed parameters:\n";
  cout << "\t Lambda_z1 (default = " << DEFVAL_LAMBDA_VI << ")\n";
  cout << "\t Lambda_w1 (default = " << DEFVAL_LAMBDA_VI << ")\n";
  cout << "\t Lambda_zgs1 (default = " << DEFVAL_LAMBDA_VI << ")\n";
  cout << "\t MZ (default = " << DEFVAL_MZ << ")\n";
  cout << "\t MW (default = " << DEFVAL_MW << ")\n";
  cout << "\t sin2ThetaW (default = " << DEFVAL_SW << ")\n";

  cout << endl;
  if (command_fail) throw std::exception();
}
