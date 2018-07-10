#ifndef MELANCSPLINEFACTORY_2D
#define MELANCSPLINEFACTORY_2D

#include <vector>
#include <utility>
#include <algorithm>
#include "TTree.h"
#include "TNumericUtil.hh"
#include "MELANCSpline_2D_fast.h"
#include "MELAFuncPdf.h"

typedef TNumericUtil::triplet<MELANCSplineCore::T> splineTriplet_t;

class MELANCSplineFactory_2D{
protected:
  TString appendName;

  MELANCSplineCore::BoundaryCondition bcBeginX;
  MELANCSplineCore::BoundaryCondition bcEndX;
  MELANCSplineCore::BoundaryCondition bcBeginY;
  MELANCSplineCore::BoundaryCondition bcEndY;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  MELANCSpline_2D_fast* fcn;
  MELAFuncPdf* PDF;

  const std::vector<splineTriplet_t> getPoints(const std::vector<MELANCSplineCore::T>& XList, const std::vector<MELANCSplineCore::T>& YList, const std::vector<MELANCSplineCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<splineTriplet_t>& pList);

  void addUnique(std::vector<MELANCSplineCore::T>& list, MELANCSplineCore::T val);

public:
  MELANCSplineFactory_2D(
    RooAbsReal& XVar_, RooAbsReal& YVar_, TString appendName_="",
    MELANCSplineCore::BoundaryCondition const bcBeginX_=MELANCSplineCore::bcNaturalSpline,
    MELANCSplineCore::BoundaryCondition const bcEndX_=MELANCSplineCore::bcNaturalSpline,
    MELANCSplineCore::BoundaryCondition const bcBeginY_=MELANCSplineCore::bcNaturalSpline,
    MELANCSplineCore::BoundaryCondition const bcEndY_=MELANCSplineCore::bcNaturalSpline
  );
  ~MELANCSplineFactory_2D();

  MELANCSpline_2D_fast* getFunc(){ return fcn; }
  MELAFuncPdf* getPDF(){ return PDF; }

  void setEndConditions(
    MELANCSplineCore::BoundaryCondition const bcBegin,
    MELANCSplineCore::BoundaryCondition const bcEnd,
    const unsigned int direction
  );

  void setPoints(TTree* tree);
  void setPoints(const std::vector<splineTriplet_t>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& YList, const std::vector<inType>& FcnList){
    std::vector<MELANCSplineCore::T> transXList;
    std::vector<MELANCSplineCore::T> transYList;
    std::vector<MELANCSplineCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((MELANCSplineCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<YList.size(); ip++) transYList.push_back((MELANCSplineCore::T)YList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((MELANCSplineCore::T)FcnList.at(ip));
    const std::vector<splineTriplet_t> pList = getPoints(transXList, transYList, transFcnList);
    setPoints(pList);
  }

};

template void MELANCSplineFactory_2D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& YList, const std::vector<Float_t>& FcnList);
template void MELANCSplineFactory_2D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& YList, const std::vector<Double_t>& FcnList);

#endif



