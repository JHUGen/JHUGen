#include <exception>
#include <iostream>
#include "LexiConHCOptionParser.h"
#include "LexiConHCTranslator.h"


int main(int argc, char** argv){
  using namespace std;

  try{
    LexiConHCOptionParser options(argc, argv);
    if (!LexiConHCOptionParser::globalHelpFlag){
      LexiConHCTranslator translator(options);
      bool firstEntry = true;
      for (auto const& it:translator.getParameters()){
        if (!firstEntry) cout << " ";
        cout << it.first << "=" << it.second;
        firstEntry = false;
      }
      for (auto const& it:translator.getCouplings()){
        if (!firstEntry) cout << " ";
        cout << it.first << "=" << it.second.first << "," << it.second.second;
        firstEntry = false;
      }
      cout << endl;
    }
    return 0;
  }
  catch (const std::exception&){
    return 1;
  }
}