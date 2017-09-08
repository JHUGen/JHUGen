#ifndef MELANCSPLINEFACTORY_1D
#define MELANCSPLINEFACTORY_1D

#include <vector>
#include <utility>
#include <algorithm>
#include "TGraph.h"
#include "TTree.h"
#include "MELANCSpline_1D_fast.h"
#include "MELAFuncPdf.h"

class MELANCSplineFactory_1D{
protected:
  TString appendName;

  RooAbsReal* splineVar;
  MELANCSpline_1D_fast* fcn;
  MELAFuncPdf* PDF;

  const std::vector<std::pair<MELANCSplineCore::T, MELANCSplineCore::T>> getPoints(const std::vector<MELANCSplineCore::T>& XList, const std::vector<MELANCSplineCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<std::pair<MELANCSplineCore::T, MELANCSplineCore::T>>& pList);

public:
  MELANCSplineFactory_1D(RooAbsReal& splineVar_, TString appendName_="");
  ~MELANCSplineFactory_1D();

  MELANCSpline_1D_fast* getFunc(){ return fcn; }
  MELAFuncPdf* getPDF(){ return PDF; }

  void setPoints(TTree* tree);
  void setPoints(TGraph* tg);
  void setPoints(const std::vector<std::pair<MELANCSplineCore::T, MELANCSplineCore::T>>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& FcnList){
    std::vector<MELANCSplineCore::T> transXList;
    std::vector<MELANCSplineCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((MELANCSplineCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((MELANCSplineCore::T)FcnList.at(ip));
    const std::vector<std::pair<MELANCSplineCore::T, MELANCSplineCore::T>> pList = getPoints(transXList, transFcnList);
    initPDF(pList);
  }

};

template void MELANCSplineFactory_1D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& FcnList);
template void MELANCSplineFactory_1D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& FcnList);

#endif



