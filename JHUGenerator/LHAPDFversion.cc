#include <LHAPDF/Version.h>
#include <LHAPDF/PDFIndex.h>
#include <cassert>
#include <cstring>

extern "C" void lhapdfversion_(char *version) {
  std::string theversion = LHAPDF::version();
  if (theversion.size() > 5) assert(0);
  strcpy(version, theversion.data());
}
