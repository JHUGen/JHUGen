#ifndef MELAOUTPUTSTREAMER_H
#define MELAOUTPUTSTREAMER_H

#include <ios>
#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <streambuf>
#include <sstream>
#include <iomanip>
#include <string>
#include "TString.h"


class MELAOutputStreamer{
protected:
  std::ofstream theFile;
  std::ostream* stdout_ptr;

public:
  MELAOutputStreamer(const char* fname, std::ios_base::openmode fmode = std::ios_base::out, bool printError=false);
  ~MELAOutputStreamer();

  template<typename T> MELAOutputStreamer& operator<<(const T& val);
  MELAOutputStreamer& operator<<(std::ostream& (*fcn)(std::ostream&));
  MELAOutputStreamer& operator<<(std::ios& (*fcn)(std::ios&));
  MELAOutputStreamer& operator<<(std::ios_base& (*fcn)(std::ios_base&));

  std::streamsize width() const;
  std::streamsize width(std::streamsize wide);

  char fill() const;
  char fill(char fillch);

  void open(const char* fname, std::ios_base::openmode fmode = std::ios_base::out);
  void close();

  template<typename T> void writeCentered(const T& val, char fillch=' ', std::streamsize gapsize=0);

};

template<typename T> MELAOutputStreamer& MELAOutputStreamer::operator<<(const T& val){
  theFile << val;
  if (stdout_ptr) *stdout_ptr << val;
  return *this;
}
template MELAOutputStreamer& MELAOutputStreamer::operator<< <bool>(const bool& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <unsigned short>(const unsigned short& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <short>(const short& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <unsigned int>(const unsigned int& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <int>(const int& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <unsigned long>(const unsigned long& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <long>(const long& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <float>(const float& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <double>(const double& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <long double>(const long double& val);

template MELAOutputStreamer& MELAOutputStreamer::operator<< <char*>(char* const& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <char const*>(char const* const& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <char>(const char& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <signed char*>(signed char* const& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <signed char const*>(signed char const* const& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <signed char>(const signed char& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <unsigned char*>(unsigned char* const& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <unsigned char const*>(unsigned char const* const& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <unsigned char>(const unsigned char& val);

template MELAOutputStreamer& MELAOutputStreamer::operator<< <std::string>(std::string const& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <TString>(TString const& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <std::streambuf*>(std::streambuf* const& val);
template MELAOutputStreamer& MELAOutputStreamer::operator<< <void*>(void* const& val);

template<typename T> void MELAOutputStreamer::writeCentered(const T& val, char fillch, std::streamsize gapsize){
  char deffillch = this->fill(fillch);

  std::stringstream tmpss;
  tmpss << val;
  std::string tmpstr = tmpss.str();
  std::streamsize strlength = (std::streamsize) tmpstr.length();

  if (strlength>gapsize) *this << std::setw(gapsize) << "";
  else{
    std::streamsize leftgap = (gapsize+strlength)/2;
    std::streamsize rightgap = gapsize-leftgap;
    *this << std::setw(leftgap) << tmpstr << std::setw(rightgap) << "";
  }

  this->fill(deffillch);
}
template void MELAOutputStreamer::writeCentered<bool>(const bool& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<unsigned short>(const unsigned short& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<short>(const short& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<unsigned int>(const unsigned int& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<int>(const int& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<unsigned long>(const unsigned long& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<float>(const float& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<double>(const double& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<long double>(const long double& val, char fillch, std::streamsize gapsize);

template void MELAOutputStreamer::writeCentered<char*>(char* const& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<char const*>(char const* const& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<char>(char const& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<signed char*>(signed char* const& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<signed char const*>(signed char const* const& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<signed char>(signed char const& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<unsigned char*>(unsigned char* const& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<unsigned char const*>(unsigned char const* const& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<unsigned char>(unsigned char const& val, char fillch, std::streamsize gapsize);

template void MELAOutputStreamer::writeCentered<std::string>(std::string const& val, char fillch, std::streamsize gapsize);
template void MELAOutputStreamer::writeCentered<TString>(TString const& val, char fillch, std::streamsize gapsize);



#endif
