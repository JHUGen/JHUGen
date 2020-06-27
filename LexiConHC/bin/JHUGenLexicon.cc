#include <exception>
#include <iostream>
#include "JHUGenLexiconOptionParser.h"
#include "JHUGenLexiconTranslator.h"


int main(int argc, char** argv){
  using namespace std;

  try{
    JHUGenLexiconOptionParser options(argc, argv);
    if (!JHUGenLexiconOptionParser::globalHelpFlag){
      JHUGenLexiconTranslator translator(options);
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
