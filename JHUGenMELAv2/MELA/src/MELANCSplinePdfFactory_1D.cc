#include "MELANCSplinePdfFactory_1D.h"
#include <cassert>

using namespace std;


MELANCSplinePdfFactory_1D::MELANCSplinePdfFactory_1D(RooAbsReal& splineVar_, TString appendName_) :
appendName(appendName_),
splineVar(&splineVar_),
PDF(0)
{}
MELANCSplinePdfFactory_1D::~MELANCSplinePdfFactory_1D(){
  destroyPDF();
}
void MELANCSplinePdfFactory_1D::setPoints(TTree* tree){
  vector<pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>> pList;
  MELANCSplinePdfCore::T x, fcn;
  tree->SetBranchAddress("X", &x);
  tree->SetBranchAddress("Fcn", &fcn);
  int n = tree->GetEntries();
  for (int ip=0; ip<n; ip++){ tree->GetEntry(ip); pList.push_back(pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>(x, fcn)); }
  setPoints(pList);
}
void MELANCSplinePdfFactory_1D::setPoints(TGraph* tg){
  vector<pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>> pList;
  double* xx = tg->GetX();
  double* yy = tg->GetY();
  int n = tg->GetN();
  for (int ip=0; ip<n; ip++) pList.push_back(pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>(xx[ip], yy[ip]));
  setPoints(pList);
}
const std::vector<std::pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>> MELANCSplinePdfFactory_1D::getPoints(const std::vector<MELANCSplinePdfCore::T>& XList, const std::vector<MELANCSplinePdfCore::T>& FcnList){
  const unsigned int nX = XList.size();
  const unsigned int n = FcnList.size();
  if (nX!=n){
    cerr << "MELANCSplinePdfFactory_1D::getPoints: nX=" << nX << " != nFcn=" << n << endl;
    assert(0);
  }
  std::vector<std::pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>> pList; pList.reserve(n);
  for (unsigned int ip=0; ip<n; ip++) pList.push_back(pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>(XList.at(ip), FcnList.at(ip)));
  return pList;
}

void MELANCSplinePdfFactory_1D::destroyPDF(){ delete PDF; PDF=0; }
void MELANCSplinePdfFactory_1D::initPDF(const std::vector<std::pair<MELANCSplinePdfCore::T, MELANCSplinePdfCore::T>>& pList){
  destroyPDF();

  const unsigned int n = pList.size();
  std::vector<MELANCSplinePdfCore::T> XList;
  std::vector<MELANCSplinePdfCore::T> FcnList;
  for (unsigned int ip=0; ip<n; ip++){
    XList.push_back(pList.at(ip).first);
    FcnList.push_back(pList.at(ip).second);
  }

  TString name = "PDF";
  if (appendName!="") name = Form("%s_%s", name.Data(), appendName.Data());
  TString title=name;
  PDF = new MELANCSplinePdf_1D_fast(
    name.Data(),
    title.Data(),
    *splineVar,
    XList, FcnList
    );
}
MELANCSplinePdf_1D_fast* MELANCSplinePdfFactory_1D::getPDF(){ return PDF; }
