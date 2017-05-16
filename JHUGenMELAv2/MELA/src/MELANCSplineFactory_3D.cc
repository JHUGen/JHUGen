#include "MELANCSplineFactory_3D.h"

using namespace std;


MELANCSplineFactory_3D::MELANCSplineFactory_3D(RooAbsReal& XVar_, RooAbsReal& YVar_, RooAbsReal& ZVar_, TString appendName_) :
appendName(appendName_),
XVar(&XVar_), YVar(&YVar_), ZVar(&ZVar_),
fcn(0),
PDF(0)
{}
MELANCSplineFactory_3D::~MELANCSplineFactory_3D(){
  destroyPDF();
}

void MELANCSplineFactory_3D::addUnique(std::vector<MELANCSplineCore::T>& list, MELANCSplineCore::T val){
  for (unsigned int ip=0; ip<list.size(); ip++){ if (list.at(ip)==val) return; }
  list.push_back(val);
}
const std::vector<splineQuadruplet_t> MELANCSplineFactory_3D::getPoints(
  const std::vector<MELANCSplineCore::T>& XList,
  const std::vector<MELANCSplineCore::T>& YList,
  const std::vector<MELANCSplineCore::T>& ZList,
  const std::vector<MELANCSplineCore::T>& FcnList
  ){
  const unsigned int nX = XList.size();
  const unsigned int nY = YList.size();
  const unsigned int nZ = ZList.size();
  const unsigned int n = FcnList.size();
  if (nX*nY*nZ!=n){
    cerr << "MELANCSplineFactory_3D::getPoints: nX=" << nX << " x nY=" << nY << " x nZ=" << nZ << " != nFcn=" << n << endl;
    assert(0);
  }

  std::vector<splineQuadruplet_t> pList; pList.reserve(n);
  for (unsigned int ix=0; ix<nX; ix++){
    MELANCSplineCore::T xval = XList.at(ix);
    for (unsigned int iy=0; iy<nY; iy++){
      MELANCSplineCore::T yval = YList.at(iy);
      for (unsigned int iz=0; iz<nZ; iz++){
        MELANCSplineCore::T zval = ZList.at(iz);
        unsigned int ip = nZ*(nY*ix + iy) + iz;
        pList.push_back(splineQuadruplet_t(xval, yval, zval, FcnList.at(ip)));
      }
    }
  }
  return pList;
}

void MELANCSplineFactory_3D::destroyPDF(){ delete PDF; PDF=0; delete fcn; fcn=0; }
void MELANCSplineFactory_3D::initPDF(const std::vector<splineQuadruplet_t>& pList){
  destroyPDF();

  const unsigned int n = pList.size();
  vector<MELANCSplineCore::T> XList;
  vector<MELANCSplineCore::T> YList;
  vector<MELANCSplineCore::T> ZList;
  vector<vector<vector<MELANCSplineCore::T>>> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    addUnique(XList, (pList.at(ip))[0]);
    addUnique(YList, (pList.at(ip))[1]);
    addUnique(ZList, (pList.at(ip))[2]);
  }
  FcnList.reserve(ZList.size());
  for (unsigned int iz=0; iz<ZList.size(); iz++){
    vector<vector<MELANCSplineCore::T>> dumz;
    dumz.reserve(YList.size());
    for (unsigned int iy=0; iy<YList.size(); iy++){
      vector<MELANCSplineCore::T> dumy;
      dumy.reserve(XList.size());
      for (unsigned int ix=0; ix<XList.size(); ix++){
        unsigned int ip = ZList.size()*(YList.size()*ix + iy) + iz;
        dumy.push_back((pList.at(ip))[3]); // Do not use unique here
      }
      dumz.push_back(dumy);
    }
    FcnList.push_back(dumz);
  }

  TString name = "Func";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  fcn = new MELANCSpline_3D_fast(
    name.Data(),
    title.Data(),
    *XVar, *YVar, *ZVar,
    XList, YList, ZList, FcnList
    );

  name.Prepend("PDF_"); title=name;
  PDF = new MELAFuncPdf(
    name.Data(),
    title.Data(),
    *fcn
    );
}

void MELANCSplineFactory_3D::setPoints(TTree* tree){
  vector<splineQuadruplet_t> pList;
  MELANCSplineCore::T x, y, z, fcn;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Y", &y);
  tree->SetBranchAddress("Z", &z);
  tree->SetBranchAddress("Fcn", &fcn);
  int n = tree->GetEntries();
  for (int ip=0; ip<n; ip++){ tree->GetEntry(ip); pList.push_back(splineQuadruplet_t(x, y, z, fcn)); }
  setPoints(pList);
}
