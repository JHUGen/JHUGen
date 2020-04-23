#include <exception>
#include <iostream>
#include "LexiConDMOptionParser.h"
#include "LexiConDMTranslator.h"


int main(int argc, char** argv){
  using namespace std;

  try{
    LexiConDMOptionParser options(argc, argv);
    if (!LexiConDMOptionParser::globalHelpFlag){
      LexiConDMTranslator translator(options);
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