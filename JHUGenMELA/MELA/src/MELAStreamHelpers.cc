#include "MELAStreamHelpers.hh"


namespace MELAStreamHelpers{
  MELAOutputStreamer MELAout("", std::ios_base::out, false);
  MELAOutputStreamer MELAerr("", std::ios_base::out, true);
}
