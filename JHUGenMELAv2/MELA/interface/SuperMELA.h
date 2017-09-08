#ifndef MELA_SUPERMELA_H
#define MELA_SUPERMELA_H


#include <string>
#include <fstream>
#include <vector>
#include "Mela.h"
#include "MELACombinePdfs.h"
#include "Riostream.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "TLorentzVector.h"



class SuperMELA{

public:

  SuperMELA(double mH=125, string channel="4mu", double LHCsqrts=8.);
  ~SuperMELA();
  void init();

  double GetSigShapeSystematic(string parName);
  double GetSigShapeParameter(string parName);

  void SetVerbosity(bool verb=true){ verbose_=verb; }
  void SetDecayChannel(string myChan);
  void SetMH(double myMH){
    mHVal_=myMH;
    mH_rrv_->setVal(mHVal_);
    if (verbose_)std::cout<<"Setting MH to "<<mHVal_<<std::endl;
    init();
  }

  void SetPathToCards(string dirToCards){
    pathToCards_=dirToCards;
    if (verbose_)std::cout<<"New path to cards is "<<pathToCards_.c_str()<<std::endl;
  }

  std::pair<double, double> M4lProb(double m4l);
  std::pair<double, double> M4lProb(std::pair<double, double>);

protected:

  void splitLine(const string rawoption, vector<string>& splitoptions, char delimiter);
  void readSigParsFromFile(
    string& str_mean_CB,
    string& str_sigma_CB,
    string& str_n_CB,
    string& str_alpha_CB,
    string& str_n2_CB,
    string& str_alpha2_CB
    );
  void readBkgParsFromFile(std::vector<double>& apars);
  void readSigSystFromFile(
    double& str_mean_CB_err_e,
    double& str_mean_CB_err_m,
    double& str_sigma_CB_err_e,
    double& str_sigma_CB_err_m
    );

  void calc_mZZ_range(const double mHVal, double& low_M, double& high_M);
  bool checkChannel();

  double mHVal_;
  double sqrts_;
  double lowMH_, highMH_;
  string strChan_;
  int ch_;
  bool verbose_;
  string pathToCards_;

  RooRealVar* m4l_rrv_;
  RooRealVar* mH_rrv_;

  RooFormulaVar* n_CB_;
  RooFormulaVar* alpha_CB_;
  RooFormulaVar* n2_CB_;
  RooFormulaVar* alpha2_CB_;
  RooFormulaVar* mean_CB_;
  RooFormulaVar* sigma_CB_;
  RooFormulaVar* meanTOT_CB_;

  RooRealVar* mean_CB_err_;
  RooRealVar* sigma_CB_err_;

  RooRealVar* corr_mean_sig;
  RooRealVar* corr_sigma_sig;

  MELADoubleCB *sig_CB_;
  MELARelBWUFParam *sig_BW_;
  RooFFTConvPdf *sig_FFT_;
  RooRealVar* mean_BW_;
  RooRealVar* width_BW_;
  double norm_sig_CB_, norm_sig_FFT_;

  //qqZZ background m4l shape
  RooRealVar* a0_qqZZ_;
  RooRealVar* a1_qqZZ_;
  RooRealVar* a2_qqZZ_;
  RooRealVar* a3_qqZZ_;
  RooRealVar* a4_qqZZ_;
  RooRealVar* a5_qqZZ_;
  RooRealVar* a6_qqZZ_;
  RooRealVar* a7_qqZZ_;
  RooRealVar* a8_qqZZ_;
  RooRealVar* a9_qqZZ_;
  RooRealVar* a10_qqZZ_;
  RooRealVar* a11_qqZZ_;
  RooRealVar* a12_qqZZ_;
  RooRealVar* a13_qqZZ_;
  MELAqqZZPdf_v2 *qqZZ_pdf_;
  double norm_bkg_qqZZ_;

};


#endif
