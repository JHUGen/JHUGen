#ifndef MELANCSPLINEPDFFACTORY_1D
#define MELANCSPLINEPDFFACTORY_1D

#include <vector>
#include <utility>
#include <algorithm>
#include "MELANCSplinePdf_1D_fast.h"
#include "TGraph.h"

class MELANCSplinePdfFactory_1D{
protected:
  TString appendName;

  RooAbsReal* splineVar;
  MELANCSplinePdf_1D_fast* PDF;

  const std::vector<std::pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>> getPoints(const std::vector<MELANCSplinePdfCore::T>& XList, const std::vector<MELANCSplinePdfCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<std::pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>>& pList);

public:
  MELANCSplinePdfFactory_1D(RooAbsReal& splineVar_, TString appendName_="");
  ~MELANCSplinePdfFactory_1D();

  MELANCSplinePdf_1D_fast* getPDF();
  void setPoints(TGraph* tg);
  void setPoints(const std::vector<std::pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& FcnList){
    std::vector<MELANCSplinePdfCore::T> transXList;
    std::vector<MELANCSplinePdfCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((MELANCSplinePdfCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((MELANCSplinePdfCore::T)FcnList.at(ip));
    const std::vector<std::pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>> pList = getPoints(transXList, transFcnList);
    initPDF(pList);
  }

};

template void MELANCSplinePdfFactory_1D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& FcnList);
template void MELANCSplinePdfFactory_1D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& FcnList);

#endif



