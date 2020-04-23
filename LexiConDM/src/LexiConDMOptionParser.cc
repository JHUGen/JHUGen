#include <exception>
#include <cstdlib>
#include <iostream>
#include "LexiConDMCouplings.h"
#include "LexiConDMOptionParser.h"
#include "LexiConDMHelperFunctions.h"


using namespace std;
using namespace LexiConDMHelperFunctions;
using namespace LexiConDMIOHelpers;


bool LexiConDMOptionParser::globalHelpFlag = false;

LexiConDMOptionParser::LexiConDMOptionParser(int argc, char** argv) :
  basis_input(nIOBases),
  basis_output(nIOBases)
{
  for (int a=0; a<argc; a++){
    string tmpArg(argv[a]);
    rawOptions.push_back(tmpArg);
  }

  analyze();
}
void LexiConDMOptionParser::analyze(){
  bool hasInvalidOption=false;
  char rawdelimiter = '=';
  for (unsigned int opt=1; opt<rawOptions.size(); opt++){
    string wish, value;
    splitOption(rawOptions.at(opt), wish, value, rawdelimiter);

    bool isInvalidOption=false;
    interpretOption(wish, value, isInvalidOption);
    hasInvalidOption |= isInvalidOption;
  }

  if (LexiConDMOptionParser::globalHelpFlag){ printOptionsHelp(false); return; }

  if (basis_input == nIOBases || basis_output == nIOBases){ cerr << "LexiConDMOptionParser::analyze: You have to specify the input and output basis types." << endl; if (!hasInvalidOption) hasInvalidOption=true; }

  // Print help if needed and abort at this point, nowhere later
  if (hasInvalidOption) printOptionsHelp(hasInvalidOption);
}

void LexiConDMOptionParser::interpretOption(std::string const& wish, std::string const& value, bool& invalidOption){
  invalidOption=false;
  if (wish.empty()){
    if (value=="help") LexiConDMOptionParser::globalHelpFlag = true;
    else{
      cerr << "LexiConDMOptionParser::interpretOption: Unknown unspecified argument: " << value << endl;
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

void LexiConDMOptionParser::printOptionsHelp(bool command_fail)const{
  cout << endl;
  cout << "The options implemented in LexiConDM (format: specifier=value):\n\n";

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
