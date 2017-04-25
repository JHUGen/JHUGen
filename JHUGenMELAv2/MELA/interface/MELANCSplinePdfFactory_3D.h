#ifndef MELANCSPLINEPDFFACTORY_3D
#define MELANCSPLINEPDFFACTORY_3D

#include <vector>
#include <utility>
#include <algorithm>
#include "TMCFMUtils.hh"
#include "RooConstVar.h"
#include "MELANCSplinePdf_3D_fast.h"

typedef TMCFMUtils::quadruplet<MELANCSplinePdfCore::T> splineQuadruplet_t;

class MELANCSplinePdfFactory_3D{
protected:
  TString appendName;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  RooAbsReal* ZVar;
  MELANCSplinePdf_3D_fast* PDF;

  const std::vector<splineQuadruplet_t> getPoints(const std::vector<MELANCSplinePdfCore::T>& XList, const std::vector<MELANCSplinePdfCore::T>& YList, const std::vector<MELANCSplinePdfCore::T>& ZList, const std::vector<MELANCSplinePdfCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<splineQuadruplet_t>& pList);

  void addUnique(std::vector<MELANCSplinePdfCore::T>& list, MELANCSplinePdfCore::T val);

public:
  MELANCSplinePdfFactory_3D(RooAbsReal& XVar_, RooAbsReal& YVar_, RooAbsReal& ZVar_, TString appendName_="");
  ~MELANCSplinePdfFactory_3D();

  MELANCSplinePdf_3D_fast* getPDF();

  void setPoints(const std::vector<splineQuadruplet_t>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& YList, const std::vector<inType>& ZList, const std::vector<inType>& FcnList){
    std::vector<MELANCSplinePdfCore::T> transXList;
    std::vector<MELANCSplinePdfCore::T> transYList;
    std::vector<MELANCSplinePdfCore::T> transZList;
    std::vector<MELANCSplinePdfCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((MELANCSplinePdfCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<YList.size(); ip++) transYList.push_back((MELANCSplinePdfCore::T)YList.at(ip));
    for (unsigned int ip=0; ip<ZList.size(); ip++) transZList.push_back((MELANCSplinePdfCore::T)ZList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((MELANCSplinePdfCore::T)FcnList.at(ip));
    const std::vector<splineQuadruplet_t> pList = getPoints(transXList, transYList, transZList, transFcnList);
    setPoints(pList);
  }

};

template void MELANCSplinePdfFactory_3D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& YList, const std::vector<Float_t>& ZList, const std::vector<Float_t>& FcnList);
template void MELANCSplinePdfFactory_3D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& YList, const std::vector<Double_t>& ZList, const std::vector<Double_t>& FcnList);


#endif



