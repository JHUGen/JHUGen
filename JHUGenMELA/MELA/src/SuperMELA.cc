#include <sstream>
#include <cmath>
#include <cassert>
#include "MELAStreamHelpers.hh"
#include "TVar.hh"
#include "SuperMELA.h"
#include "MELAHXSWidth.h"
#include "RooArgSet.h"
#include "RooArgList.h"


using namespace std;
using namespace RooFit;
using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;


SuperMELA::SuperMELA(
  double mH,
  string channel,
  double LHCsqrts
  ) :
  mHVal_(mH),
  sqrts_(LHCsqrts)
{
  sigma_CB_=nullptr;
  mean_CB_err_=nullptr;
  sigma_CB_err_=nullptr;

  n_CB_=nullptr;
  alpha_CB_=nullptr;
  n2_CB_=nullptr;
  alpha2_CB_=nullptr;
  mean_CB_=nullptr;
  meanTOT_CB_=nullptr;

  mean_BW_=nullptr;
  width_BW_=nullptr;

  corr_mean_sig=nullptr;
  corr_sigma_sig=nullptr;

  sig_CB_=nullptr;
  sig_BW_=nullptr;
  sig_FFT_=nullptr;

  a0_qqZZ_=nullptr;
  a1_qqZZ_=nullptr;
  a2_qqZZ_=nullptr;
  a3_qqZZ_=nullptr;
  a4_qqZZ_=nullptr;
  a5_qqZZ_=nullptr;
  a6_qqZZ_=nullptr;
  a7_qqZZ_=nullptr;
  a8_qqZZ_=nullptr;
  a9_qqZZ_=nullptr;
  a10_qqZZ_=nullptr;
  a11_qqZZ_=nullptr;
  a12_qqZZ_=nullptr;
  a13_qqZZ_=nullptr;

  qqZZ_pdf_=nullptr;

  verbose_=false;
  mH_rrv_=new RooRealVar("mH", "mH", mHVal_, 0., sqrts_*1000.);
  m4l_rrv_=nullptr;
  strChan_=channel;

  //pathToCards_="../../../HiggsAnalysis/HZZ4L_CombinationPy/CreateDatacards/SM_inputs_8TeV/";
  pathToCards_="../data/CombinationInputs/SM_inputs_8TeV/";
  //init();
}

SuperMELA::~SuperMELA(){
  // Delete in reverse order wrt SuperMELA::init()
  delete qqZZ_pdf_; qqZZ_pdf_=nullptr;

  delete a0_qqZZ_; delete a1_qqZZ_; delete a2_qqZZ_; delete a3_qqZZ_;
  delete a4_qqZZ_; delete a5_qqZZ_; delete a6_qqZZ_; delete a7_qqZZ_;
  delete a8_qqZZ_; delete a9_qqZZ_; delete a10_qqZZ_; delete a11_qqZZ_;
  delete a12_qqZZ_; delete a13_qqZZ_;
  a0_qqZZ_=nullptr; a1_qqZZ_=nullptr; a2_qqZZ_=nullptr; a3_qqZZ_=nullptr;
  a4_qqZZ_=nullptr; a5_qqZZ_=nullptr; a6_qqZZ_=nullptr; a7_qqZZ_=nullptr;
  a8_qqZZ_=nullptr; a9_qqZZ_=nullptr; a10_qqZZ_=nullptr; a11_qqZZ_=nullptr;
  a12_qqZZ_=nullptr; a13_qqZZ_=nullptr;

  delete sig_FFT_; sig_FFT_=nullptr;
  delete sig_BW_; sig_BW_=nullptr;
  delete sig_CB_; sig_CB_=nullptr;

  delete width_BW_; width_BW_=nullptr;
  delete mean_BW_; mean_BW_=nullptr;

  delete sigma_CB_; sigma_CB_=nullptr;
  delete meanTOT_CB_; meanTOT_CB_=nullptr;
  delete mean_CB_; mean_CB_=nullptr;

  delete corr_sigma_sig; corr_sigma_sig=nullptr;
  delete corr_mean_sig; corr_mean_sig=nullptr;

  delete alpha2_CB_; alpha2_CB_=nullptr;
  delete n2_CB_; n2_CB_=nullptr;
  delete alpha_CB_; alpha_CB_=nullptr;
  delete n_CB_; n_CB_=nullptr;

  delete sigma_CB_err_; sigma_CB_err_=nullptr;
  delete mean_CB_err_; mean_CB_err_=nullptr;

  delete m4l_rrv_; m4l_rrv_=nullptr;
  delete mH_rrv_; mH_rrv_=nullptr;
}

void SuperMELA::SetDecayChannel(string myChan){
  if (verbose_) MELAout << "SuperMELA::SetDecayChannel: Switching from " << strChan_ << " to " << myChan << endl;
  if (myChan != strChan_){ // do nothing if it's the same as before
    strChan_=myChan;
    bool newChanOK = checkChannel();
    if (verbose_) MELAout << "Setting decay channel of SuperMELA to " << strChan_.c_str() << " , re-initializing..." << endl;
    init();
    if (verbose_ && newChanOK) MELAout << "Decay channel set successfully to " << strChan_.c_str() << endl;
  }
}

double SuperMELA::GetSigShapeSystematic(string parName){
  if (parName=="meanCB") return mean_CB_err_->getVal();
  else if (parName=="sigmaCB") return sigma_CB_err_->getVal();
  else{
    MELAout << "Error from SuperMELA::GetSigShapeSystematic, unrecognized input: " << parName.c_str() << endl;
    try{
      throw 40;
    }
    catch (int e){
      if (e==40) MELAout << "Exception " << e << "  in SuperMELA::GetSigShapeSystematic! Unrecognized type of parameter requested: " << parName.c_str() << "  . Valid options are meanCB ; sigmaCB" << endl;
      return 0;
    }
  }
  return 0;
}

double SuperMELA::GetSigShapeParameter(string parName){
  if (parName=="meanCB") return meanTOT_CB_->getVal();
  else if (parName=="sigmaCB") return sigma_CB_->getVal();
  else if (parName=="alphaCB") return alpha_CB_->getVal();
  else if (parName=="nCB") return n_CB_->getVal();
  else{
    try{
      throw 40;
    }
    catch (int e){
      if (e==40) MELAout << "Exception " << e << "  in SuperMELA::GetSigShapeParameter! Unrecognized type of parameter requested: " << parName.c_str() << "  . Valid options are meanCB ; sigmaCB ; alphaCB ; nCB" << endl;
      return 0;
    }
  }
  return 0;
}


void SuperMELA::init(){
  if (verbose_)MELAout << "Begin SuperMELA::init..." << endl;

  // Calculate m4l ranges for the given mH, set range of rrv
  calc_mZZ_range(mHVal_, lowMH_, highMH_);
  if (verbose_) MELAout << "Range width=" << highMH_ - lowMH_ << endl;
  delete m4l_rrv_; m4l_rrv_=new RooRealVar("CMS_zz4l_mass", "CMS_zz4l_mass", mHVal_, lowMH_, highMH_);
  m4l_rrv_->setBins(2000, "fft");
  m4l_rrv_->setRange("shape", lowMH_, highMH_);

  //parameters for signal m4l shape systematics
  double str_mean_CB_err_e;
  double str_mean_CB_err_m;
  double str_sigma_CB_err_e;
  double str_sigma_CB_err_m;
  readSigSystFromFile(str_mean_CB_err_e, str_mean_CB_err_m, str_sigma_CB_err_e, str_sigma_CB_err_m);
  if (verbose_){
    MELAout << "mean_CB systematics (ele and mu): " << str_mean_CB_err_e << " / " << str_mean_CB_err_m << endl;
    MELAout << "sigma_CB systematics (ele and mu): " << str_sigma_CB_err_e << " / " << str_sigma_CB_err_m << endl;
  }

  //delete old stuff before reinitialization
  delete mean_CB_err_;
  delete sigma_CB_err_;
  if (strChan_=="4mu"){
    mean_CB_err_=new RooRealVar("mean_CB_err", "mean_CB_err", str_mean_CB_err_m);
    sigma_CB_err_=new RooRealVar("sigma_CB_err", "sigma_CB_err", str_sigma_CB_err_m);
  }
  else if (strChan_=="4e"){
    mean_CB_err_=new RooRealVar("mean_CB_err", "mean_CB_err", str_mean_CB_err_e);
    sigma_CB_err_=new RooRealVar("sigma_CB_err", "sigma_CB_err", str_sigma_CB_err_e);
  }
  else{//2e2mu, we should have already checked that the string of the channel is a sensible one
    mean_CB_err_=new RooRealVar("mean_CB_err", "mean_CB_err", (str_mean_CB_err_m + str_mean_CB_err_e));
    sigma_CB_err_=new RooRealVar("sigma_CB_err", "sigma_CB_err", sqrt(pow(str_sigma_CB_err_m, 2) + pow(str_sigma_CB_err_e, 2)));
  }
  if (mean_CB_err_->getVal()<0.){ MELAout << "Negative error on the m4l mean ! " << mean_CB_err_->getVal() << endl; }
  if (sigma_CB_err_->getVal()<0.){ MELAout << "Negative error on the m4l sigma ! " << sigma_CB_err_->getVal() << endl; }


  //set parameters for signal m4l shape and calculate normalization
  string str_n_CB;
  string str_alpha_CB;
  string str_n2_CB;
  string str_alpha2_CB;
  string str_mean_CB;
  string str_sigma_CB;
  if (verbose_)MELAout << "Reading signal shape formulas" << endl;
  readSigParsFromFile(str_mean_CB, str_sigma_CB, str_n_CB, str_alpha_CB, str_n2_CB, str_alpha2_CB);
  if (verbose_){
    MELAout << "Read from input card the following formulas: " << endl;
    MELAout << "Mean RooFormula (string): " << str_mean_CB.c_str() << endl;
    MELAout << "Sigma RooFormula (string): " << str_sigma_CB.c_str() << endl;
  }

  delete n_CB_;
  delete alpha_CB_;
  delete n2_CB_;
  delete alpha2_CB_;
  char rrvName[96];
  sprintf(rrvName, "CMS_zz4l_n_sig_%s_%d", strChan_.c_str(), int(sqrts_));
  if (verbose_) MELAout << "SuperMELA::init: Constructing n_CB_ from formula " << str_n_CB.c_str() << endl;
  n_CB_=new RooFormulaVar(rrvName, str_n_CB.c_str(), RooArgList(*mH_rrv_));
  sprintf(rrvName, "CMS_zz4l_alpha_sig_%s_%d", strChan_.c_str(), int(sqrts_));
  if (verbose_) MELAout << "SuperMELA::init: Constructing alpha_CB_ from formula " << str_alpha_CB.c_str() << endl;
  alpha_CB_=new RooFormulaVar(rrvName, str_alpha_CB.c_str(), RooArgList(*mH_rrv_));
  sprintf(rrvName, "CMS_zz4l_n2_sig_%s_%d", strChan_.c_str(), int(sqrts_));
  if (verbose_) MELAout << "SuperMELA::init: Constructing n2_CB_ from formula " << str_n2_CB.c_str() << endl;
  n2_CB_=new RooFormulaVar(rrvName, str_n2_CB.c_str(), RooArgList(*mH_rrv_));
  sprintf(rrvName, "CMS_zz4l_alpha2_sig_%s_%d", strChan_.c_str(), int(sqrts_));
  if (verbose_) MELAout << "SuperMELA::init: Constructing alpha2_CB_ from formula " << str_alpha2_CB.c_str() << endl;
  alpha2_CB_=new RooFormulaVar(rrvName, str_alpha2_CB.c_str(), RooArgList(*mH_rrv_));

  corr_mean_sig = new RooRealVar("CMS_zz4l_mean_sig_corrMH", "CMS_zz4l_mean_sig_corrMH", 0., -10., 10.);
  corr_sigma_sig = new RooRealVar("CMS_zz4l_sigma_sig_corrMH", "CMS_zz4l_sigma_sig_corrMH", 0., -10., 10.);
  corr_mean_sig->setConstant(true);
  corr_sigma_sig->setConstant(true);

  delete mean_CB_;
  delete meanTOT_CB_;
  delete sigma_CB_;
  mean_CB_=new RooFormulaVar("CMS_zz4l_mean_m_sig", Form("(%s)+@0*@1", str_mean_CB.c_str()), RooArgList(*mH_rrv_, *corr_mean_sig));//this is normalized by mHVal_
  meanTOT_CB_=new RooFormulaVar("CMS_zz4l_mean_sig", "(@0+@1)", RooArgList(*mH_rrv_, *mean_CB_));
  if (verbose_){ MELAout << "Signal Mean vals -> Correction: " << corr_mean_sig->getVal() << "  Mean: " << mean_CB_->getVal() << "  Total: " << meanTOT_CB_->getVal() << endl; }
  sigma_CB_=new RooFormulaVar("CMS_zz4l_sigma_m_sig", Form("(%s)*(1+@0*@1)", str_sigma_CB.c_str()), RooArgList(*mH_rrv_, *corr_sigma_sig));

  //for high-mass one also needs a gamma RooFormulaVar, 
  delete mean_BW_;
  delete width_BW_;
  sprintf(rrvName, "CMS_zz4l_mean_BW_sig_%s_%d", strChan_.c_str(), int(sqrts_));
  mean_BW_=new RooRealVar(rrvName, "CMS_zz4l_mean_BW", mHVal_, 100., 1000.);
  sprintf(rrvName, "CMS_zz4l_width_BW_sig_%s_%d", strChan_.c_str(), int(sqrts_));
  width_BW_=new RooRealVar(rrvName, "CMS_zz4l_width_BW", 1.);
  mean_BW_->setVal(mHVal_);
  mean_BW_->setConstant(true);
  width_BW_->setConstant(true);

  if (verbose_){
    MELAout << "Signal shape parameter values: " << endl;
    MELAout << "Mean (formula value) = " << meanTOT_CB_->getVal() << endl;
    MELAout << "Sigma (formula value) = " << sigma_CB_->getVal() << endl;
    MELAout << "n (formula value) = " << n_CB_->getVal() << endl;
    MELAout << "alpha (formula value) = " << alpha_CB_->getVal() << endl;
    MELAout << "n2 (formula value) = " << n2_CB_->getVal() << endl;
    MELAout << "alpha2 (formula value) = " << alpha2_CB_->getVal() << endl;
    MELAout << "Mean BW (realvar value) = " << mean_BW_->getVal() << endl;
    MELAout << "Width BW (realvar value) = " << width_BW_->getVal() << endl;
  }

  delete sig_FFT_;
  delete sig_CB_;
  delete sig_BW_;
  sig_CB_ =new MELADoubleCB("signalCB_ggH", "signalCB_ggH", *m4l_rrv_, *meanTOT_CB_, *sigma_CB_, *alpha_CB_, *n_CB_, *alpha2_CB_, *n2_CB_);
  sig_BW_ =new MELARelBWUFParam("signalBW_ggH", "signalBW_ggH", *m4l_rrv_, *mean_BW_, *width_BW_);
  //sig_FFT_=new RooFFTConvPdf("signal_ggH", "BW (X) CB", *m4l_rrv_, *sig_BW_, *sig_CB_, 2);
  //sig_FFT_->setBufferFraction(0.2);
  if (verbose_){
    MELAout << "Value of signal m4l CB shape is " << sig_CB_->getVal() << endl;
    MELAout << "Value of signal m4l BW shape is " << sig_BW_->getVal() << endl;
    //sig_FFT_->Print("v");
    //MELAout << "Value of signal m4l FFT shape is " << sig_FFT_->getVal() << endl;
  }

  RooAbsReal* tmpint;

  tmpint = sig_CB_->createIntegral(RooArgSet(*m4l_rrv_), RooFit::Range("shape"));
  norm_sig_CB_ =tmpint->getVal();
  delete tmpint;
  if (verbose_)MELAout << "Normalization of signal m4l CB shape is " << norm_sig_CB_ << endl;

  if (verbose_)MELAout << "\n---> Integrating Breit-Wigner:" << endl;
  tmpint = sig_BW_->createIntegral(RooArgSet(*m4l_rrv_), RooFit::Range("shape"));
  double norm_sig_BW_ =tmpint->getVal();
  delete tmpint;
  if (verbose_)MELAout << "Normalization of signal m4l BW shape is " << norm_sig_BW_ << endl;

  /*
  if (verbose_)MELAout << "\n---> Integrating full signal:" << endl;
  tmpint = sig_FFT_->createIntegral(RooArgSet(*m4l_rrv_), RooFit::Range("shape"));
  norm_sig_FFT_=tmpint->getVal();
  delete tmpint;
  if (verbose_)MELAout << "Normalization of signal m4l shape is " << norm_sig_FFT_ << endl;
  */

  if (verbose_)MELAout << "Reading background shape parameters" << endl;
  vector<double> v_apars;
  readBkgParsFromFile(v_apars);

  if (verbose_){
    MELAout << "Size of vector with bkg shape pars is " << v_apars.size() << endl;
    MELAout << "Param [0]=" << v_apars.at(0) << " [13]=" << v_apars.at(13) << endl;
  }

  delete a0_qqZZ_; a0_qqZZ_=new RooRealVar("CMS_zz4l_a0_qqZZ", "CMS_zz4l_a0_qqZZ", v_apars.at(0), 0., 200.);
  delete a1_qqZZ_; a1_qqZZ_=new RooRealVar("CMS_zz4l_a1_qqZZ", "CMS_zz4l_a1_qqZZ", v_apars.at(1), 0., 200.);
  delete a2_qqZZ_; a2_qqZZ_=new RooRealVar("CMS_zz4l_a2_qqZZ", "CMS_zz4l_a2_qqZZ", v_apars.at(2), 0., 200.);
  delete a3_qqZZ_; a3_qqZZ_=new RooRealVar("CMS_zz4l_a3_qqZZ", "CMS_zz4l_a3_qqZZ", v_apars.at(3), 0., 1.);
  delete a4_qqZZ_; a4_qqZZ_=new RooRealVar("CMS_zz4l_a4_qqZZ", "CMS_zz4l_a4_qqZZ", v_apars.at(4), 0., 200.);
  delete a5_qqZZ_; a5_qqZZ_=new RooRealVar("CMS_zz4l_a5_qqZZ", "CMS_zz4l_a5_qqZZ", v_apars.at(5), 0., 200.);
  delete a6_qqZZ_; a6_qqZZ_=new RooRealVar("CMS_zz4l_a6_qqZZ", "CMS_zz4l_a6_qqZZ", v_apars.at(6), 0., 100.);
  delete a7_qqZZ_; a7_qqZZ_=new RooRealVar("CMS_zz4l_a7_qqZZ", "CMS_zz4l_a7_qqZZ", v_apars.at(7), 0., 1.);
  delete a8_qqZZ_; a8_qqZZ_=new RooRealVar("CMS_zz4l_a8_qqZZ", "CMS_zz4l_a8_qqZZ", v_apars.at(8), 0., 200.);
  delete a9_qqZZ_; a9_qqZZ_=new RooRealVar("CMS_zz4l_a9_qqZZ", "CMS_zz4l_a9_qqZZ", v_apars.at(9), 0., 1.);
  delete a10_qqZZ_; a10_qqZZ_=new RooRealVar("CMS_zz4l_a10_qqZZ", "CMS_zz4l_a10_qqZZ", v_apars.at(10), 0., 200.);
  delete a11_qqZZ_; a11_qqZZ_=new RooRealVar("CMS_zz4l_a11_qqZZ", "CMS_zz4l_a11_qqZZ", v_apars.at(11), -100., 100.);
  delete a12_qqZZ_; a12_qqZZ_=new RooRealVar("CMS_zz4l_a12_qqZZ", "CMS_zz4l_a12_qqZZ", v_apars.at(12), 0., 10000.);
  delete a13_qqZZ_; a13_qqZZ_=new RooRealVar("CMS_zz4l_a13_qqZZ", "CMS_zz4l_a13_qqZZ", v_apars.at(13), 0., 1.);
  a0_qqZZ_->setConstant(kTRUE);
  a1_qqZZ_->setConstant(kTRUE);
  a2_qqZZ_->setConstant(kTRUE);
  a3_qqZZ_->setConstant(kTRUE);
  a4_qqZZ_->setConstant(kTRUE);
  a5_qqZZ_->setConstant(kTRUE);
  a6_qqZZ_->setConstant(kTRUE);
  a7_qqZZ_->setConstant(kTRUE);
  a8_qqZZ_->setConstant(kTRUE);
  a9_qqZZ_->setConstant(kTRUE);
  a10_qqZZ_->setConstant(kTRUE);
  a11_qqZZ_->setConstant(kTRUE);
  a12_qqZZ_->setConstant(kTRUE);
  a13_qqZZ_->setConstant(kTRUE);


  delete qqZZ_pdf_; qqZZ_pdf_ = new MELAqqZZPdf_v2("bkg_qqzz", "bkg_qqzz", *m4l_rrv_, *a0_qqZZ_, *a1_qqZZ_, *a2_qqZZ_, *a3_qqZZ_,
    *a4_qqZZ_, *a5_qqZZ_, *a6_qqZZ_, *a7_qqZZ_,
    *a8_qqZZ_, *a9_qqZZ_, *a10_qqZZ_, *a11_qqZZ_,
    *a12_qqZZ_, *a13_qqZZ_);

  tmpint = qqZZ_pdf_->createIntegral(RooArgSet(*m4l_rrv_), RooFit::Range("shape"));
  norm_bkg_qqZZ_=tmpint->getVal();
  delete tmpint;
}

void SuperMELA::SetMH(double myMH){
  mHVal_=myMH;
  mH_rrv_->setVal(mHVal_);
  if (verbose_) MELAout << "Setting MH to " << mHVal_ << endl;
  init();
}

void SuperMELA::SetPathToCards(string dirToCards){
  pathToCards_=dirToCards;
  if (verbose_) MELAout << "New path to cards is " << pathToCards_ << endl;
}

void SuperMELA::splitLine(const string rawoption, vector<string>& splitoptions, char delimiter){
  string suboption=rawoption, result=rawoption;
  string remnant;
  while (result!=""){
    size_t posEq = suboption.find(delimiter);
    if (posEq!=string::npos){
      result=suboption;
      remnant=suboption.substr(posEq+1);
      result.erase(result.begin()+posEq, result.end());
    }
    else{
      result="";
      remnant=suboption;
    }
    if (result!="") splitoptions.push_back(result);
    suboption = remnant;
  }
  if (remnant!="") splitoptions.push_back(remnant);
}

void SuperMELA::readSigSystFromFile(
  double& str_mean_CB_err_e,
  double& str_mean_CB_err_m,
  double& str_sigma_CB_err_e,
  double& str_sigma_CB_err_m
  ){
  bool mean_e_OK=false, sigma_e_OK=false;
  bool mean_m_OK=false, sigma_m_OK=false;

  //open text file
  string fCardName=pathToCards_+"inputs_"+strChan_+".txt";
  if (verbose_) MELAout << "SuperMELA::readSigSystFromFile: Parsing signal shape systematics from input card " << fCardName.c_str() << endl;
  ifstream card(fCardName.c_str(), ios::in);
  if (!card.good()){
    MELAerr << "SuperMELA::readSigSystFromFile: Input card " << fCardName << " is not good!" << endl;
    assert(0);
  }
  else if (verbose_) MELAout << "SuperMELA::readSigSystFromFile: Input card " << fCardName << " can be read!" << endl;

  string line;
  while (card.good()){
    getline(card, line);
    vector<string> fields;
    splitLine(line, fields, ' ');

    if (fields.size()<2 || !(fields[0]=="systematic"&&fields[1]=="param")) continue;

    if (fields.size()!=4){
      MELAerr << "Error in SuperMELA::readSigSystFromFile! Incorrect format of line " << line.c_str() << endl;
      break;
    }

    if (fields[2]=="CMS_zz4l_mean_m_sig"){
      str_mean_CB_err_m = atof(fields.at(3).c_str());
      mean_m_OK=true;
    }
    if (fields[2]=="CMS_zz4l_mean_e_sig"){
      str_mean_CB_err_e = atof(fields.at(3).c_str());
      mean_e_OK=true;
    }
    if (fields[2]=="CMS_zz4l_sigma_m_sig"){
      str_sigma_CB_err_m = atof(fields.at(3).c_str());
      sigma_m_OK=true;
    }
    if (fields[2]=="CMS_zz4l_sigma_e_sig"){
      str_sigma_CB_err_e = atof(fields.at(3).c_str());
      sigma_e_OK=true;
    }

    if (mean_e_OK && sigma_e_OK && mean_m_OK && sigma_m_OK) break;
  }

  try{
    if ((!(mean_e_OK&&sigma_e_OK))&&(!(mean_m_OK&&sigma_m_OK))) throw 20;
  }
  catch (int e){
    MELAout
      << "Exception " << e << " in SuperMELA::readSigSystFromFile! Not all signal shape formulas were read " << mean_e_OK << "  " << sigma_e_OK << "  " << mean_m_OK << "  " << sigma_m_OK
      << endl;
  }

  card.close();
}

void SuperMELA::readSigParsFromFile(
  string& str_mean_CB,
  string& str_sigma_CB,
  string& str_n_CB,
  string& str_alpha_CB,
  string& str_n2_CB,
  string& str_alpha2_CB
  ){
  bool meanOK=false, sigmaOK=false, nOK=false, alphaOK=false, n2OK=false, alpha2OK=false;
  //open text file
  string fCardName=pathToCards_+"inputs_"+strChan_+".txt";
  if (verbose_) MELAout << "Parsing input card " << fCardName.c_str() << endl;
  ifstream card(fCardName.c_str(), ios::in);
  string line;
  while (card.good()){
    getline(card, line);
    vector<string> fields;
    splitLine(line, fields, ' ');
    if (fields.size()==0 || fields[0]!="signalShape")continue;
    //ok, we found somethign interesting
    if (fields.size()<3){
      MELAout << "Error in SuperMELA::readSigParsFromFile! Incorrect format of line " << line.c_str() << " There should be at least three fields, I could find only " << fields.size() << " fields." << endl;
      break;
    }
    if (fields.at(1)=="n_CB"){ str_n_CB=fields.at(2); nOK=true; }
    if (fields.at(1)=="alpha_CB"){ str_alpha_CB=fields.at(2); alphaOK=true; }
    if (fields.at(1)=="n2_CB"){ str_n2_CB=fields.at(2); n2OK=true; }
    if (fields.at(1)=="alpha2_CB"){ str_alpha2_CB=fields.at(2); alpha2OK=true; }
    if (fields.at(1)=="mean_CB"){ str_mean_CB=fields.at(2); meanOK=true; }
    if (fields.at(1)=="sigma_CB"){ str_sigma_CB=fields.at(2); sigmaOK=true; }
    if (verbose_) MELAout << fields.at(1) << " == " << fields.at(2) << endl;

    if (meanOK && sigmaOK && alphaOK && nOK && alpha2OK && n2OK)break;
  }//end while loop on lines

  if (!(meanOK && sigmaOK && alphaOK && nOK && alpha2OK && n2OK)){
    MELAout << "Exception  in SuperMELA::readSigParsFromFile! Not all signal shape formulas were read " << meanOK << " " << sigmaOK << "  " << alphaOK << "  " << nOK << "  " << alpha2OK << "  " << n2OK << endl;
    throw 20;
  }

  card.close();
}

void SuperMELA::readBkgParsFromFile(vector<double>& apars){
  const int nPars=14;
  int nFound=0;
  apars.resize(14);
  string fCardName=pathToCards_+"inputs_"+strChan_+".txt";
  if (verbose_)MELAout << "Parsing input card " << fCardName.c_str() << endl;
  ifstream card(fCardName.c_str(), ios::in);
  string line;

  while (card.good()){
    getline(card, line);
    vector<string> fields;
    splitLine(line, fields, ' ');
    if (fields.size()==0 || fields[0]!="qqZZshape")continue;

    if (fields.size()<3){
      MELAout << "Error in SuperMELA::readSigParsFromFile! Incorrect format of line \'" << line.c_str() << "\' . It contains " << fields.size() << "fields (it should be 3)" << endl;
      break;
    }

    stringstream ssip;
    ssip << nFound;
    if (fields[1]=="a"+ssip.str()+"_bkgd"){
      apars[nFound]=atof(fields[2].c_str());
      nFound++;
    }
  }

  try{
    if (nFound!=nPars){
      throw 30;
    }
  }
  catch (int e){
    if (e==30){
      MELAerr << "Exception from void SuperMELA::readBkgParsFromFile(vector<double> apars ). Mismatched number of params of qqZZ shape read from card file " << fCardName.c_str() << " ---> " << nFound << " (it should be " << nPars << endl;
      for (unsigned int j=0; j<apars.size(); j++){ MELAerr << apars[j] << "   " << flush; }
      MELAerr << endl;
    }
  }
  card.close();
}

void SuperMELA::calc_mZZ_range(const double mHVal, double& low_M, double& high_M){
  //low_M=0.;
  //high_M=sqrts_*1000.;

#ifdef _melapkgpathstr_
  const string MELAPKGPATH = _melapkgpathstr_;
#else
  MELAout << "SuperMELA::calc_mZZ_range: MELA package path is undefined! Please modify the makefle or the makefile-equivalent!" << endl;
  assert(0);
#endif
  string path = MELAPKGPATH + "data/HiggsTotalWidth_YR3.txt";

  path.erase((find(path.rbegin(), path.rend(), '/').base()), path.end());
  MELAHXSWidth myCSW(path.c_str());
  double widthHVal =  myCSW.HiggsWidth(mHVal);
  double windowVal = max(widthHVal, 1.);
  double lowside = 100.;
  if (mHVal >= 275){ lowside = 180.; }
  else { lowside = 100.; }
  // Apply rounding
  low_M = int(max((mHVal - 20.*windowVal), lowside)+0.5);
  high_M = int(min((mHVal + 15.*windowVal), sqrts_*1000.)+0.5);
}

pair<double, double> SuperMELA::M4lProb(double m4l){
  double Psig =-1.;
  double Pbkg =-1.;
  if (m4l<lowMH_ || m4l>highMH_){
    if (verbose_) MELAerr << "SuperMELA::M4lProb: m4l = " << m4l << " outside range [" << lowMH_ << ", " << highMH_ << "]. Setting SuperMELA to dummy values." << endl;
  }
  else{
    m4l_rrv_->setVal(m4l);
    Psig=sig_CB_->getVal() / norm_sig_CB_;
    Pbkg=qqZZ_pdf_->getVal() / norm_bkg_qqZZ_;
    if (verbose_){
      MELAout << "SuperMELA::M4lProb: m4l = " << m4l
        << ", m4lPsig = " << Psig
        << ", m4lPbkg = " << Pbkg
        << endl;
    }
  }
  return make_pair(Psig, Pbkg);
}


pair<double, double> SuperMELA::M4lProb(pair<double, double> m4lPair){
  double Psig =-1.;
  double Pbkg =-1.;
  if (
    (m4lPair.first<lowMH_ || m4lPair.first>highMH_)
    ||
    (m4lPair.second<lowMH_ || m4lPair.second>highMH_)
    ){
    if (verbose_) MELAerr << "SuperMELA::M4lProb: m4l pairs " << m4lPair.first << " - " << m4lPair.second << " outside range [" << lowMH_ << ", " << highMH_ << "]. Setting SuperMELA to dummy values." << endl;
  }
  else{
    m4l_rrv_->setVal(m4lPair.first);
    Psig=sig_CB_->getVal() / norm_sig_CB_;
    m4l_rrv_->setVal(m4lPair.second);
    Pbkg=qqZZ_pdf_->getVal() / norm_bkg_qqZZ_;
    if (verbose_){
      MELAout << "SuperMELA::M4lProb: "
        << "m4lPsig(" << m4lPair.first << ") = " << Psig
        << ", m4lPbkg(" << m4lPair.second << ") = " << Pbkg
        << endl;
    }
  }
  return make_pair(Psig, Pbkg);
}

bool SuperMELA::checkChannel(){
  try{
    if (strChan_!="4mu" && strChan_!="4e" && strChan_!="2e2mu"  && strChan_!="2mu2e"){
      throw 10;
    }
  }
  catch (int e){
    MELAerr << "Exception " << e << " from SuperMELA::SetDecayChannel(string myChan). Unrecognized string for decay channel: " << strChan_.c_str() << endl;
    return false;
  }

  if (strChan_=="4mu") ch_=0;
  else if (strChan_=="4e") ch_=1;
  else ch_=2;

  if (strChan_=="2mu2e") strChan_="2e2mu";
  return true;
}

