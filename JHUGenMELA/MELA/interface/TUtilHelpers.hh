#ifndef TUTILHELPERS
#define TUTILHELPERS

#include <vector>
#include <iterator>
#include <algorithm>
#include <utility>
#include "MELAParticle.h"
#include "MELACandidate.h"
#include "MELATopCandidate.h"
#include "TVar.hh"
#include "TNumericUtil.hh"

namespace TUtilHelpers{

  template<typename T> void copyVector(std::vector<T> const& input, std::vector<T>& target);
  template<typename T, typename U=T> bool checkElementExists(T const& element, std::vector<U> const& elementlist);

}

template<typename T> void TUtilHelpers::copyVector(std::vector<T> const& input, std::vector<T>& target){
  std::copy(input.cbegin(), input.cend(), std::back_inserter(target));
}
template void TUtilHelpers::copyVector<MELAParticle*>(std::vector<MELAParticle*> const& input, std::vector<MELAParticle*>& target);
template void TUtilHelpers::copyVector<MELACandidate*>(std::vector<MELACandidate*> const& input, std::vector<MELACandidate*>& target);
template void TUtilHelpers::copyVector<MELATopCandidate*>(std::vector<MELATopCandidate*> const& input, std::vector<MELATopCandidate*>& target);
template void TUtilHelpers::copyVector<SimpleParticle_t>(SimpleParticleCollection_t const& input, SimpleParticleCollection_t& target);
template void TUtilHelpers::copyVector<bool>(std::vector<bool> const& input, std::vector<bool>& target);
template void TUtilHelpers::copyVector<short>(std::vector<short> const& input, std::vector<short>& target);
template void TUtilHelpers::copyVector<unsigned int>(std::vector<unsigned int> const& input, std::vector<unsigned int>& target);
template void TUtilHelpers::copyVector<int>(std::vector<int> const& input, std::vector<int>& target);
template void TUtilHelpers::copyVector<float>(std::vector<float> const& input, std::vector<float>& target);
template void TUtilHelpers::copyVector<double>(std::vector<double> const& input, std::vector<double>& target);
template void TUtilHelpers::copyVector<TNumericUtil::intQuad_t>(std::vector<TNumericUtil::intQuad_t> const& input, std::vector<TNumericUtil::intQuad_t>& target);

template<typename T, typename U> bool TUtilHelpers::checkElementExists(T const& element, std::vector<U> const& elementlist){
  for (U const& el:elementlist){ if (el==element) return true; }
  return false;
}
template bool TUtilHelpers::checkElementExists<MELAParticle const*>(MELAParticle const* const& element, std::vector<MELAParticle const*> const& elementlist);
template bool TUtilHelpers::checkElementExists<MELAParticle const*, MELAParticle*>(MELAParticle const* const& element, std::vector<MELAParticle*> const& elementlist);
template bool TUtilHelpers::checkElementExists<MELACandidate const*>(MELACandidate const* const& element, std::vector<MELACandidate const*> const& elementlist);
template bool TUtilHelpers::checkElementExists<MELACandidate const*, MELACandidate*>(MELACandidate const* const& element, std::vector<MELACandidate*> const& elementlist);
template bool TUtilHelpers::checkElementExists<MELATopCandidate const*>(MELATopCandidate const* const& element, std::vector<MELATopCandidate const*> const& elementlist);
template bool TUtilHelpers::checkElementExists<MELATopCandidate const*, MELATopCandidate*>(MELATopCandidate const* const& element, std::vector<MELATopCandidate*> const& elementlist);
template bool TUtilHelpers::checkElementExists<bool>(bool const& element, std::vector<bool> const& elementlist);
template bool TUtilHelpers::checkElementExists<short>(short const& element, std::vector<short> const& elementlist);
template bool TUtilHelpers::checkElementExists<unsigned int>(unsigned int const& element, std::vector<unsigned int> const& elementlist);
template bool TUtilHelpers::checkElementExists<int>(int const& element, std::vector<int> const& elementlist);
template bool TUtilHelpers::checkElementExists<float>(float const& element, std::vector<float> const& elementlist);
template bool TUtilHelpers::checkElementExists<double>(double const& element, std::vector<double> const& elementlist);
template bool TUtilHelpers::checkElementExists<TNumericUtil::intQuad_t>(TNumericUtil::intQuad_t const& element, std::vector<TNumericUtil::intQuad_t> const& elementlist);

#endif

