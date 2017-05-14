#include "MELANCSplinePdfFactory_3D.h"

using namespace std;


MELANCSplinePdfFactory_3D::MELANCSplinePdfFactory_3D(RooAbsReal& XVar_, RooAbsReal& YVar_, RooAbsReal& ZVar_, TString appendName_) :
appendName(appendName_),
XVar(&XVar_), YVar(&YVar_), ZVar(&ZVar_),
PDF(0)
{}
MELANCSplinePdfFactory_3D::~MELANCSplinePdfFactory_3D(){
  destroyPDF();
}

void MELANCSplinePdfFactory_3D::addUnique(std::vector<MELANCSplinePdfCore::T>& list, MELANCSplinePdfCore::T val){
  for (unsigned int ip=0; ip<list.size(); ip++){ if (list.at(ip)==val) return; }
  list.push_back(val);
}
const std::vector<splineQuadruplet_t> MELANCSplinePdfFactory_3D::getPoints(
  const std::vector<MELANCSplinePdfCore::T>& XList,
  const std::vector<MELANCSplinePdfCore::T>& YList,
  const std::vector<MELANCSplinePdfCore::T>& ZList,
  const std::vector<MELANCSplinePdfCore::T>& FcnList
  ){
  const unsigned int nX = XList.size();
  const unsigned int nY = YList.size();
  const unsigned int nZ = ZList.size();
  const unsigned int n = FcnList.size();
  if (nX*nY*nZ!=n){
    cerr << "MELANCSplinePdfFactory_3D::getPoints: nX=" << nX << " x nY=" << nY << " x nZ=" << nZ << " != nFcn=" << n << endl;
    assert(0);
  }

  std::vector<splineQuadruplet_t> pList; pList.reserve(n);
  for (unsigned int ix=0; ix<nX; ix++){
    MELANCSplinePdfCore::T xval = XList.at(ix);
    for (unsigned int iy=0; iy<nY; iy++){
      MELANCSplinePdfCore::T yval = YList.at(iy);
      for (unsigned int iz=0; iz<nZ; iz++){
        MELANCSplinePdfCore::T zval = ZList.at(iz);
        unsigned int ip = nZ*(nY*ix + iy) + iz;
        pList.push_back(splineQuadruplet_t(xval, yval, zval, FcnList.at(ip)));
      }
    }
  }
  return pList;
}

void MELANCSplinePdfFactory_3D::destroyPDF(){ delete PDF; PDF=0; }
void MELANCSplinePdfFactory_3D::initPDF(const std::vector<splineQuadruplet_t>& pList){
  destroyPDF();

  const unsigned int n = pList.size();
  vector<MELANCSplinePdfCore::T> XList;
  vector<MELANCSplinePdfCore::T> YList;
  vector<MELANCSplinePdfCore::T> ZList;
  vector<vector<vector<MELANCSplinePdfCore::T>>> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    addUnique(XList, (pList.at(ip))[0]);
    addUnique(YList, (pList.at(ip))[1]);
    addUnique(ZList, (pList.at(ip))[2]);
  }
  FcnList.reserve(ZList.size());
  for (unsigned int iz=0; iz<ZList.size(); iz++){
    vector<vector<MELANCSplinePdfCore::T>> dumz;
    dumz.reserve(YList.size());
    for (unsigned int iy=0; iy<YList.size(); iy++){
      vector<MELANCSplinePdfCore::T> dumy;
      dumy.reserve(XList.size());
      for (unsigned int ix=0; ix<XList.size(); ix++){
        unsigned int ip = ZList.size()*(YList.size()*ix + iy) + iz;
        dumy.push_back((pList.at(ip))[3]); // Do not use unique here
      }
      dumz.push_back(dumy);
    }
    FcnList.push_back(dumz);
  }

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new MELANCSplinePdf_3D_fast(
    name.Data(),
    title.Data(),
    *XVar, *YVar, *ZVar,
    XList, YList, ZList, FcnList
    );
}
MELANCSplinePdf_3D_fast* MELANCSplinePdfFactory_3D::getPDF(){ return PDF; }

void MELANCSplinePdfFactory_3D::setPoints(TTree* tree){
  vector<splineQuadruplet_t> pList;
  MELANCSplinePdfCore::T x, y, z, fcn;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Y", &y);
  tree->SetBranchAddress("Z", &z);
  tree->SetBranchAddress("Fcn", &fcn);
  int n = tree->GetEntries();
  for (int ip=0; ip<n; ip++){ tree->GetEntry(ip); pList.push_back(splineQuadruplet_t(x, y, z, fcn)); }
  setPoints(pList);
}
