#include "MELANCSplineFactory_1D.h"
#include <cassert>

using namespace std;


MELANCSplineFactory_1D::MELANCSplineFactory_1D(
  RooAbsReal& splineVar_, TString appendName_,
  MELANCSplineCore::BoundaryCondition const bcBeginX_,
  MELANCSplineCore::BoundaryCondition const bcEndX_
) :
  appendName(appendName_),
  bcBeginX(bcBeginX_), bcEndX(bcEndX_),
  splineVar(&splineVar_),
  fcn(0),
  PDF(0)
{}
MELANCSplineFactory_1D::~MELANCSplineFactory_1D(){
  destroyPDF();
}
void MELANCSplineFactory_1D::setPoints(TTree* tree){
  vector<pair<MELANCSplineCore::T, MELANCSplineCore::T>> pList;
  MELANCSplineCore::T x, fcn;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Fcn", &fcn);
  int n = tree->GetEntries();
  for (int ip=0; ip<n; ip++){ tree->GetEntry(ip); pList.push_back(pair<MELANCSplineCore::T, MELANCSplineCore::T>(x, fcn)); }
  setPoints(pList);
}
void MELANCSplineFactory_1D::setPoints(TGraph* tg){
  vector<pair<MELANCSplineCore::T, MELANCSplineCore::T>> pList;
  double* xx = tg->GetX();
  double* yy = tg->GetY();
  int n = tg->GetN();
  for (int ip=0; ip<n; ip++) pList.push_back(pair<MELANCSplineCore::T, MELANCSplineCore::T>(xx[ip], yy[ip]));
  setPoints(pList);
}
const std::vector<std::pair<MELANCSplineCore::T, MELANCSplineCore::T>> MELANCSplineFactory_1D::getPoints(const std::vector<MELANCSplineCore::T>& XList, const std::vector<MELANCSplineCore::T>& FcnList){
  const unsigned int nX = XList.size();
  const unsigned int n = FcnList.size();
  if (nX!=n){
    cerr << "MELANCSplineFactory_1D::getPoints: nX=" << nX << " != nFcn=" << n << endl;
    assert(0);
  }
  std::vector<std::pair<MELANCSplineCore::T, MELANCSplineCore::T>> pList; pList.reserve(n);
  for (unsigned int ip=0; ip<n; ip++) pList.push_back(pair<MELANCSplineCore::T, MELANCSplineCore::T>(XList.at(ip), FcnList.at(ip)));
  return pList;
}

void MELANCSplineFactory_1D::destroyPDF(){ delete PDF; PDF=0; delete fcn; fcn=0; }
void MELANCSplineFactory_1D::initPDF(const std::vector<std::pair<MELANCSplineCore::T, MELANCSplineCore::T>>& pList){
  destroyPDF();

  const unsigned int n = pList.size();
  std::vector<MELANCSplineCore::T> XList;
  std::vector<MELANCSplineCore::T> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    XList.push_back(pList.at(ip).first);
    FcnList.push_back(pList.at(ip).second);
  }

  TString name = "Func";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  fcn = new MELANCSpline_1D_fast(
    name.Data(),
    title.Data(),
    *splineVar,
    XList, FcnList,
    bcBeginX, bcEndX
    );

  name.Prepend("PDF_"); title=name;
  PDF = new MELAFuncPdf(
    name.Data(),
    title.Data(),
    *fcn
    );
}

void MELANCSplineFactory_1D::setEndConditions(
  MELANCSplineCore::BoundaryCondition const bcBegin,
  MELANCSplineCore::BoundaryCondition const bcEnd,
  const unsigned int /*direction*/
){
  bcBeginX=bcBegin;
  bcEndX=bcEnd;
}
