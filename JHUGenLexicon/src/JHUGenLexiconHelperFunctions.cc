#include <cassert>
#include <cctype>
#include <algorithm>
#include "JHUGenLexiconHelperFunctions.h"


using namespace std;


template<> void JHUGenLexiconHelperFunctions::lowercase(std::string const& name, std::string& val){
  val = name;
  std::transform(val.begin(), val.end(), val.begin(), [] (unsigned char c){ return std::tolower(c); });
}
template<> void JHUGenLexiconHelperFunctions::lowercase(const char* const& name, const char*& val){
  std::string strname = name;
  std::string strval;
  lowercase(strname, strval);
  val = strval.data();
}

template<> void JHUGenLexiconHelperFunctions::castStringToValue(std::string const& name, bool& val){
  std::string namelower=name;
  std::transform(namelower.begin(), namelower.end(), namelower.begin(), ::tolower);
  if (namelower=="true" || namelower=="t") val=true;
  else if (namelower=="false" || namelower=="f") val=false;
  else{ std::stringstream ss(name); ss >> val; }
}

void JHUGenLexiconHelperFunctions::splitOption(const std::string& rawoption, std::string& wish, std::string& value, char delimiter){
  size_t posEq = rawoption.find(delimiter);
  if (posEq!=string::npos){
    wish=rawoption;
    value=rawoption.substr(posEq+1);
    wish.erase(wish.begin()+posEq, wish.end());
  }
  else{
    wish="";
    value=rawoption;
  }
}
void JHUGenLexiconHelperFunctions::splitOptionRecursive(const std::string& rawoption, std::vector<std::string>& splitoptions, char delimiter, bool uniqueResults){
  string suboption=rawoption, result=rawoption;
  string remnant;
  while (result!=""){
    splitOption(suboption, result, remnant, delimiter);
    if (result!="" && (!uniqueResults || (uniqueResults && !checkListVariable(splitoptions, result)))) splitoptions.push_back(result);
    suboption = remnant;
    if (result=="" && suboption.find(delimiter)!=std::string::npos) result=suboption; // This can happen if the string starts with the delimiter.
  }
  if (remnant!="" && (!uniqueResults || (uniqueResults && !checkListVariable(splitoptions, remnant)))) splitoptions.push_back(remnant);
}
