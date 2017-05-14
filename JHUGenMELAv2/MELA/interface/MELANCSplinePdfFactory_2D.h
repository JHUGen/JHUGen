#ifndef MELANCSPLINEPDFFACTORY_2D
#define MELANCSPLINEPDFFACTORY_2D

#include <vector>
#include <utility>
#include <algorithm>
#include "TTree.h"
#include "TMCFMUtils.hh"
#include "RooConstVar.h"
#include "MELANCSplinePdf_2D_fast.h"

typedef TMCFMUtils::triplet<MELANCSplinePdfCore::T> splineTriplet_t;

class MELANCSplinePdfFactory_2D{
protected:
  TString appendName;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  MELANCSplinePdf_2D_fast* PDF;

  const std::vector<splineTriplet_t> getPoints(const std::vector<MELANCSplinePdfCore::T>& XList, const std::vector<MELANCSplinePdfCore::T>& YList, const std::vector<MELANCSplinePdfCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<splineTriplet_t>& pList);

  void addUnique(std::vector<MELANCSplinePdfCore::T>& list, MELANCSplinePdfCore::T val);

public:
  MELANCSplinePdfFactory_2D(RooAbsReal& XVar_, RooAbsReal& YVar_, TString appendName_="");
  ~MELANCSplinePdfFactory_2D();

  MELANCSplinePdf_2D_fast* getPDF();

  void setPoints(TTree* tree);
  void setPoints(const std::vector<splineTriplet_t>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& YList, const std::vector<inType>& FcnList){
    std::vector<MELANCSplinePdfCore::T> transXList;
    std::vector<MELANCSplinePdfCore::T> transYList;
    std::vector<MELANCSplinePdfCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((MELANCSplinePdfCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<YList.size(); ip++) transYList.push_back((MELANCSplinePdfCore::T)YList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((MELANCSplinePdfCore::T)FcnList.at(ip));
    const std::vector<splineTriplet_t> pList = getPoints(transXList, transYList, transFcnList);
    setPoints(pList);
  }

};

template void MELANCSplinePdfFactory_2D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& YList, const std::vector<Float_t>& FcnList);
template void MELANCSplinePdfFactory_2D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& YList, const std::vector<Double_t>& FcnList);

#endif



