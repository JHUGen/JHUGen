#include "MELAOutputStreamer.h"


using namespace std;


MELAOutputStreamer::MELAOutputStreamer(const char* fname, std::ios_base::openmode fmode, bool printError) :
theFile(fname, fmode)
{
  if (!printError) stdout_ptr = &std::cout;
  else stdout_ptr = &std::cerr;
}
MELAOutputStreamer::~MELAOutputStreamer(){ this->close(); }

MELAOutputStreamer& MELAOutputStreamer::operator<<(std::ostream& (*fcn)(std::ostream&)){
  fcn(theFile);
  if (stdout_ptr) fcn(*stdout_ptr);
  return *this;
}
MELAOutputStreamer& MELAOutputStreamer::operator<<(std::ios& (*fcn)(std::ios&)){
  fcn(theFile);
  if (stdout_ptr) fcn(*stdout_ptr);
  return *this;
}
MELAOutputStreamer& MELAOutputStreamer::operator<<(std::ios_base& (*fcn)(std::ios_base&)){
  fcn(theFile);
  if (stdout_ptr) fcn(*stdout_ptr);
  return *this;
}

std::streamsize MELAOutputStreamer::width() const{ return theFile.width(); }
std::streamsize MELAOutputStreamer::width(std::streamsize wide){
  if (stdout_ptr) stdout_ptr->width(wide);
  return theFile.width(wide);
}

char MELAOutputStreamer::fill() const{ return theFile.fill(); }
char MELAOutputStreamer::fill(char fillch){
  if (stdout_ptr) stdout_ptr->fill(fillch);
  return theFile.fill(fillch);
}

void MELAOutputStreamer::close(){
  theFile.close();
}
void MELAOutputStreamer::open(const char* fname, std::ios_base::openmode fmode){
  this->close();
  theFile.open(fname, fmode);
}
