#ifndef MELANCSPLINEPDFFACTORY_2D
#define MELANCSPLINEPDFFACTORY_2D

#include <vector>
#include <utility>
#include <algorithm>
#include "RooConstVar.h"
#include "MELANCSplinePdf_2D_fast.h"

template<typename T> struct triplet{
  T value[3];
  triplet(T i1, T i2, T i3){
    value[0]=i1;
    value[1]=i2;
    value[2]=i3;
  }
  triplet(T i1){
    value[0]=i1;
    value[1]=i1;
    value[2]=i1;
  }
  triplet(){}
  T& operator[](std::size_t ipos){ return value[ipos]; } // Return by reference
  const T& operator[](std::size_t ipos)const{ return value[ipos]; } // Return by const reference
};
typedef triplet<MELANCSplinePdfCore::T> doubleTriplet_t;

class MELANCSplinePdfFactory_2D{
protected:
  TString appendName;

  RooAbsReal* XVar;
  RooAbsReal* YVar;
  MELANCSplinePdf_2D_fast* PDF;

  const std::vector<doubleTriplet_t> getPoints(const std::vector<MELANCSplinePdfCore::T>& XList, const std::vector<MELANCSplinePdfCore::T>& YList, const std::vector<MELANCSplinePdfCore::T>& FcnList);

  void destroyPDF();
  void initPDF(const std::vector<doubleTriplet_t>& pList);

  void addUnique(std::vector<MELANCSplinePdfCore::T>& list, MELANCSplinePdfCore::T val);

public:
  MELANCSplinePdfFactory_2D(RooAbsReal& XVar_, RooAbsReal& YVar_, TString appendName_="");
  ~MELANCSplinePdfFactory_2D();

  MELANCSplinePdf_2D_fast* getPDF();

  void setPoints(const std::vector<doubleTriplet_t>& pList){ initPDF(pList); }
  template<typename inType> void setPoints(const std::vector<inType>& XList, const std::vector<inType>& YList, const std::vector<inType>& FcnList){
    std::vector<MELANCSplinePdfCore::T> transXList;
    std::vector<MELANCSplinePdfCore::T> transYList;
    std::vector<MELANCSplinePdfCore::T> transFcnList;
    for (unsigned int ip=0; ip<XList.size(); ip++) transXList.push_back((MELANCSplinePdfCore::T)XList.at(ip));
    for (unsigned int ip=0; ip<YList.size(); ip++) transYList.push_back((MELANCSplinePdfCore::T)YList.at(ip));
    for (unsigned int ip=0; ip<FcnList.size(); ip++) transFcnList.push_back((MELANCSplinePdfCore::T)FcnList.at(ip));
    const std::vector<doubleTriplet_t> pList = getPoints(transXList, transYList, transFcnList);
    setPoints(pList);
  }

};

template void MELANCSplinePdfFactory_2D::setPoints<Float_t>(const std::vector<Float_t>& XList, const std::vector<Float_t>& YList, const std::vector<Float_t>& FcnList);
template void MELANCSplinePdfFactory_2D::setPoints<Double_t>(const std::vector<Double_t>& XList, const std::vector<Double_t>& YList, const std::vector<Double_t>& FcnList);

#endif



