#include "MELANCSplinePdfFactory_2D.h"

using namespace std;


MELANCSplinePdfFactory_2D::MELANCSplinePdfFactory_2D(RooAbsReal& XVar_, RooAbsReal& YVar_, TString appendName_) :
appendName(appendName_),
XVar(&XVar_), YVar(&YVar_),
PDF(0)
{}
MELANCSplinePdfFactory_2D::~MELANCSplinePdfFactory_2D(){
  destroyPDF();
}

void MELANCSplinePdfFactory_2D::addUnique(std::vector<MELANCSplinePdfCore::T>& list, MELANCSplinePdfCore::T val){
  for (unsigned int ip=0; ip<list.size(); ip++){ if (list.at(ip)==val) return; }
  list.push_back(val);
}
const std::vector<doubleTriplet_t> MELANCSplinePdfFactory_2D::getPoints(
  const std::vector<MELANCSplinePdfCore::T>& XList,
  const std::vector<MELANCSplinePdfCore::T>& YList,
  const std::vector<MELANCSplinePdfCore::T>& FcnList
  ){
  unsigned int nX = XList.size();
  unsigned int nY = YList.size();
  unsigned int n = FcnList.size();
  if (nX*nY!=n){
    cerr << "MELANCSplinePdfFactory_2D::getPoints: nX=" << nX << " x nY=" << nY << " != nFcn=" << n << endl;
    assert(0);
  }

  std::vector<doubleTriplet_t> pList; pList.reserve(n);
  for (unsigned int ix=0; ix<nX; ix++){
    MELANCSplinePdfCore::T xval = XList.at(ix);
    for (unsigned int iy=0; iy<nY; iy++){
      unsigned int ip = nY*ix + iy;
      MELANCSplinePdfCore::T yval = YList.at(iy);
      pList.push_back(doubleTriplet_t(xval, yval, FcnList.at(ip)));
    }
  }
  return pList;
}

void MELANCSplinePdfFactory_2D::destroyPDF(){ delete PDF; PDF=0; }
void MELANCSplinePdfFactory_2D::initPDF(const std::vector<doubleTriplet_t>& pList){
  destroyPDF();

  unsigned int n = pList.size();
  vector<MELANCSplinePdfCore::T> XList;
  vector<MELANCSplinePdfCore::T> YList;
  vector<vector<MELANCSplinePdfCore::T>> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    addUnique(XList, (pList.at(ip))[0]);
    addUnique(YList, (pList.at(ip))[1]);
  }
  FcnList.reserve(YList.size());
  for (unsigned int iy=0; iy<YList.size(); iy++){
    vector<MELANCSplinePdfCore::T> dum; dum.reserve(XList.size());
    for (unsigned int ix=0; ix<XList.size(); ix++){
      unsigned int ip = YList.size()*ix + iy;
      dum.push_back((pList.at(ip))[2]); // Do not use unique here
    }
    FcnList.push_back(dum);
  }

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new MELANCSplinePdf_2D_fast(
    name.Data(),
    title.Data(),
    *XVar, *YVar,
    XList, YList, FcnList
    );
}
MELANCSplinePdf_2D_fast* MELANCSplinePdfFactory_2D::getPDF(){ return PDF; }
