#ifndef MELANCSPLINEFACTORY_3D
#define MELANCSPLINEFACTORY_3D

#include <vector>
#include <utility>
#include <algorithm>
#include "TTree.h"
#include "TNumericUtil.hh"
#include "MELANCSpline_3D_fast.h"
#include "MELAFuncPdf.h"

typedef TNumericUtil::quadruplet<MELANCSplineCore::T> splineQuadruplet_t;

class MELANCSplineFactory_3D{
protected:
  TString appendName;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  RooAbsReal* ZVar;
  MELANCSpline_3D_fast* fcn;
  MELAFuncPdf* PDF;

  const std::vector<splineQuadruplet_t> getPoints(const std::vector<MELANCSplineCore::T>& XList, const std::vector<MELANCSplineCore::T>& YList, const std::vector<MELANCSplineCore::T>& ZList, const std::vector<MELANCSplineCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<splineQuadruplet_t>& pList);

  void addUnique(std::vector<MELANCSplineCore::T>& list, MELANCSplineCore::T val);

public:
  MELANCSplineFactory_3D(RooAbsReal& XVar_, RooAbsReal& YVar_, RooAbsReal& ZVar_, TString appendName_="");
  ~MELANCSplineFactory_3D();

  MELANCSpline_3D_fast* getFunc(){ return fcn; }
  MELAFuncPdf* getPDF(){ return PDF; }

  void setPoints(TTree* tree);
  void setPoints(const std::vector<splineQuadruplet_t>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& YList, const std::vector<inType>& ZList, const std::vector<inType>& FcnList){
    std::vector<MELANCSplineCore::T> transXList;
    std::vector<MELANCSplineCore::T> transYList;
    std::vector<MELANCSplineCore::T> transZList;
    std::vector<MELANCSplineCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((MELANCSplineCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<YList.size(); ip++) transYList.push_back((MELANCSplineCore::T)YList.at(ip));
    for (unsigned int ip=0; ip<ZList.size(); ip++) transZList.push_back((MELANCSplineCore::T)ZList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((MELANCSplineCore::T)FcnList.at(ip));
    const std::vector<splineQuadruplet_t> pList = getPoints(transXList, transYList, transZList, transFcnList);
    setPoints(pList);
  }

};

template void MELANCSplineFactory_3D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& YList, const std::vector<Float_t>& ZList, const std::vector<Float_t>& FcnList);
template void MELANCSplineFactory_3D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& YList, const std::vector<Double_t>& ZList, const std::vector<Double_t>& FcnList);


#endif



