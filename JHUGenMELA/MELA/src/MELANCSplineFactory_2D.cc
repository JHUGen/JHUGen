#include "MELANCSplineFactory_2D.h"

using namespace std;


MELANCSplineFactory_2D::MELANCSplineFactory_2D(RooAbsReal& XVar_, RooAbsReal& YVar_, TString appendName_) :
appendName(appendName_),
XVar(&XVar_), YVar(&YVar_),
fcn(0),
PDF(0)
{}
MELANCSplineFactory_2D::~MELANCSplineFactory_2D(){
  destroyPDF();
}

void MELANCSplineFactory_2D::addUnique(std::vector<MELANCSplineCore::T>& list, MELANCSplineCore::T val){
  for (unsigned int ip=0; ip<list.size(); ip++){ if (list.at(ip)==val) return; }
  list.push_back(val);
}
const std::vector<splineTriplet_t> MELANCSplineFactory_2D::getPoints(
  const std::vector<MELANCSplineCore::T>& XList,
  const std::vector<MELANCSplineCore::T>& YList,
  const std::vector<MELANCSplineCore::T>& FcnList
  ){
  const unsigned int nX = XList.size();
  const unsigned int nY = YList.size();
  const unsigned int n = FcnList.size();
  if (nX*nY!=n){
    cerr << "MELANCSplineFactory_2D::getPoints: nX=" << nX << " x nY=" << nY << " != nFcn=" << n << endl;
    assert(0);
  }

  std::vector<splineTriplet_t> pList; pList.reserve(n);
  for (unsigned int ix=0; ix<nX; ix++){
    MELANCSplineCore::T xval = XList.at(ix);
    for (unsigned int iy=0; iy<nY; iy++){
      unsigned int ip = nY*ix + iy;
      MELANCSplineCore::T yval = YList.at(iy);
      pList.push_back(splineTriplet_t(xval, yval, FcnList.at(ip)));
    }
  }
  return pList;
}

void MELANCSplineFactory_2D::destroyPDF(){ delete PDF; PDF=0; delete fcn; fcn=0; }
void MELANCSplineFactory_2D::initPDF(const std::vector<splineTriplet_t>& pList){
  destroyPDF();

  const unsigned int n = pList.size();
  vector<MELANCSplineCore::T> XList;
  vector<MELANCSplineCore::T> YList;
  vector<vector<MELANCSplineCore::T>> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    addUnique(XList, (pList.at(ip))[0]);
    addUnique(YList, (pList.at(ip))[1]);
  }
  FcnList.reserve(YList.size());
  for (unsigned int iy=0; iy<YList.size(); iy++){
    vector<MELANCSplineCore::T> dum; dum.reserve(XList.size());
    for (unsigned int ix=0; ix<XList.size(); ix++){
      unsigned int ip = YList.size()*ix + iy;
      dum.push_back((pList.at(ip))[2]); // Do not use unique here
    }
    FcnList.push_back(dum);
  }

  TString name = "Func";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  fcn = new MELANCSpline_2D_fast(
    name.Data(),
    title.Data(),
    *XVar, *YVar,
    XList, YList, FcnList
    );

  name.Prepend("PDF_"); title=name;
  PDF = new MELAFuncPdf(
    name.Data(),
    title.Data(),
    *fcn
    );
}

void MELANCSplineFactory_2D::setPoints(TTree* tree){
  vector<splineTriplet_t> pList;
  MELANCSplineCore::T x, y, fcn;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Y", &y);
  tree->SetBranchAddress("Fcn", &fcn);
  int n = tree->GetEntries();
  for (int ip=0; ip<n; ip++){ tree->GetEntry(ip); pList.push_back(splineTriplet_t(x, y, fcn)); }
  setPoints(pList);
}
