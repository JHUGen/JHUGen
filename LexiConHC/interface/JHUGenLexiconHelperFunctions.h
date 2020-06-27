#ifndef JHUGENLEXICON_HELPERFUNCTIONS_H
#define JHUGENLEXICON_HELPERFUNCTIONS_H

#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>


namespace JHUGenLexiconHelperFunctions{
  template<typename T> void lowercase(T const& name, T& val);
  template<> void lowercase(std::string const& name, std::string& val);
  template<> void lowercase(const char* const& name, const char*& val);

  template<typename T> void castStringToValue(std::string const& name, T& val);
  template<> void castStringToValue(std::string const& name, bool& val);

  template<typename T> bool checkListVariable(std::vector<T> const& list, T const& var);
  template<typename T> bool hasCommonElements(std::vector<T> const& list1, std::vector<T> const& list2);

  template<typename T, typename U> bool getUnorderedMapIterator(T const& key, std::unordered_map<T, U> const& theMap, typename std::unordered_map<T, U>::const_iterator& it);
  template<typename T, typename U> bool getUnorderedMapIterator(T const& key, std::unordered_map<T, U>& theMap, typename std::unordered_map<T, U>::const_iterator& it);
  template<typename T, typename U> bool getUnorderedMapIterator(T const& key, std::unordered_map<T, U>& theMap, typename std::unordered_map<T, U>::iterator& it);

  template<typename T, typename U> bool hasElement(std::unordered_map<T, U> const& theMap, T const& key);
  template<typename T> bool hasElement(std::vector<T> const& list, T const& var);

  template<typename T, typename U> void getValueWithDefault(std::unordered_map<T, U> const& theMap, T const& key, U& val, U const& defval);


  void splitOption(const std::string& rawoption, std::string& wish, std::string& value, char delimiter);
  void splitOptionRecursive(const std::string& rawoption, std::vector<std::string>& splitoptions, char delimiter, bool uniqueResults=true);

}

template<typename T> void JHUGenLexiconHelperFunctions::castStringToValue(std::string const& name, T& val){ std::stringstream ss(name); ss >> val; }

template<typename T> bool JHUGenLexiconHelperFunctions::checkListVariable(std::vector<T> const& list, T const& var){ return (std::find(std::begin(list), std::end(list), var)!=std::end(list)); }
template<typename T> bool JHUGenLexiconHelperFunctions::hasCommonElements(std::vector<T> const& list1, std::vector<T> const& list2){
  for (T const& el1:list1){ if (checkListVariable(list2, el1)) return true; }
  return false;
}

template<typename T, typename U> bool JHUGenLexiconHelperFunctions::getUnorderedMapIterator(T const& key, std::unordered_map<T, U> const& theMap, typename std::unordered_map<T, U>::const_iterator& it){
  it = theMap.find(key);
  return (it!=theMap.cend());
}
template<typename T, typename U> bool JHUGenLexiconHelperFunctions::getUnorderedMapIterator(T const& key, std::unordered_map<T, U>& theMap, typename std::unordered_map<T, U>::const_iterator& it){
  it = theMap.find(key);
  return (it!=theMap.cend());
}
template<typename T, typename U> bool JHUGenLexiconHelperFunctions::getUnorderedMapIterator(T const& key, std::unordered_map<T, U>& theMap, typename std::unordered_map<T, U>::iterator& it){
  it = theMap.find(key);
  return (it!=theMap.end());
}

template<typename T, typename U> bool JHUGenLexiconHelperFunctions::hasElement(std::unordered_map<T, U> const& theMap, T const& key){
  typename std::unordered_map<T, U>::const_iterator it_dummy;
  return getUnorderedMapIterator(key, theMap, it_dummy);
}
template<typename T> bool JHUGenLexiconHelperFunctions::hasElement(std::vector<T> const& list, T const& var){ return checkListVariable(list, var); }

template<typename T, typename U> void JHUGenLexiconHelperFunctions::getValueWithDefault(std::unordered_map<T, U> const& theMap, T const& key, U& val, U const& defval){
  typename std::unordered_map<T, U>::const_iterator it;
  if (getUnorderedMapIterator(key, theMap, it)) val = it->second;
  else val = defval;
}


#endif
