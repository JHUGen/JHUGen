#ifndef MELADIFERMIONRESOLUTIONMODEL_H
#define MELADIFERMIONRESOLUTIONMODEL_H


#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "Riostream.h"
#include "RooAbsReal.h"
#include "TUtil.hh"
#include "MELANCSplineFactory_1D.h"


struct SingleNCSpline1DResolution{
  typedef MELANCSplineCore::T T;
  RooAbsReal& varReco;
  MELANCSplineFactory_1D* spFactory;
  const T valTrue;

  SingleNCSpline1DResolution(RooAbsReal& var_, TGraphErrors* tg, const T valTrue_, TString appendName) : varReco(var_), valTrue(valTrue_) {
    TString spname = appendName;
    spname += Form("_At_%.4f", valTrue);
    spFactory = new MELANCSplineFactory_1D(varReco, spname);
    if (tg!=0) spFactory->setPoints(tg);
  }
  ~SingleNCSpline1DResolution(){ delete spFactory; }

  MELANCSpline_1D_fast* getPDF(){ return spFactory->getPDF(); }

};

struct SingleNCSpline2DResolution{
  RooAbsReal& varReco;
  RooAbsReal& varTrue;
  MELANCSplineFactory_2D* spFactory;

  SingleNCSpline2DResolution(RooAbsReal& varReco_, RooAbsReal& varTrue_, const std::vector<splineTriplet_t>& pList, TString appendName) : varReco(varReco_), varTrue(varTrue_) {
    TString spname = appendName;
    spname += "_varReco_VS_varTrue";
    spFactory = new MELANCSplineFactory_2D(varReco, varTrue, spname);
    spFactory->setPoints(pList);
  }
  ~SingleNCSpline2DResolution(){ delete spFactory; }

  MELANCSpline_2D_fast* getPDF(){ return spFactory->getPDF(); }

};

class MELADifermionResolutionModel{
protected:
  TVar::Production prod;
  float sqrts;
  float varTrueMin;
  float varTrueMax;
  RooRealVar* varReco;
  RooRealVar* varTrue;
  SingleNCSpline2DResolution* fullResolution;

  float interpolateResolution(float val, float mbw, float gbw){
    // Need to multiply by a true BW if needed

  }

public:
  MELADifermionResolutionModel(TVar::Production prod_, float sqrts_, TString filename) : prod(prod_), sqrts(sqrts_){
    TFile* finput = TFile::Open(filename, "read");
    TGraphErrors* tg_avg_reco = (TGraphErrors*)finput->Get("tg_avg_varReco");
    TH2F* h_true_reco = (TGraphErrors*)finput->Get("tg_avg_varReco");
    varTrueMin = h_true_reco->GetXaxis()->GetBinLowEdge(1);
    varTrueMax = h_true_reco->GetXaxis()->GetBinUpEdge(h_true_reco->GetNbinsX());
    int nbins = h_true_reco->GetNbinsX();
    if (nbins!=tg_avg_reco->GetN()){
      cerr
        << "Number of bins in 2D histogram not the same as number of points on the graph! Please check the input file."
        << endl;
      finput->Close();
      assert(0);
    }

    std::vector<splineTriplet_t> points;
    // LEFT HERE: Need to either evaluate points or convert resolution to 2D already in the input

    TString appendName = TVar::ProductionName(prod));
    varReco = new RooRealVar(Form("%s_varReco", appendName.Data(), "", 0., sqrts); // Could go from 0 to sqrts
    varTrue = new RooRealVar(Form("%s_varTrue", appendName.Data(), "", varTrueMin, varTrueMax); // The resolution function truth range is only valid within certain boundaries.
    RecoOverTrueMinusOne = new RooFormulaVar(Form("%s_varReco_Over_varTrue_Minus_One", appendName.Data(), "@1/@0-1.", RooArgList(*varReco, *varTrue));

    fullResolution = new SingleNCSpline2DResolution(RecoOverTrueMinusOne, varTrue, points, appendName);

    finput->Close();
  }
  ~MELADifermionResolutionModel(){ delete fullResolution; delete RecoOverTrueMinusOne; delete varTrue; delete varReco; }

  float getVal(float val, float mbw=-1., float gbw=0.){ return interpolateResolution(val, mbw, gbw); }
  bool isProduction(TVar::Production prod_){ return (prod_==prod); }

};



#endif
