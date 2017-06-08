#include <iostream>
#include <vector>
#include <cmath>
#include "MelaPConstant.h"
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"


using namespace std;


MelaPConstant::MelaPConstant(
  TVar::MatrixElement me_,
  TVar::Production prod_,
  TVar::Process proc_,
  const char* path,
  const char* spname
  ) :
  processME(me_),
  processProd(prod_),
  processProc(proc_),
  fcnLow(0),
  fcnHigh(0),
  fcnMid(0)
{
  GetFcnFromFile(path, spname);
}

MelaPConstant::~MelaPConstant(){}

void MelaPConstant::GetFcnFromFile(const char* path, const char* spname){
  TString spname_core = spname;
  spname_core.Append("_Smooth");
  spname_core.Prepend("sp_tg_");
  //cout << "MelaPConstant::GetFcnFromFile: Extracting " << spname_core << " from file " << path << endl;

  TFile* fin = TFile::Open(path, "read");
  gROOT->cd();

  if (fin!=0 && !fin->IsZombie() && fin->IsOpen()){
    TString spname_low = spname_core; spname_low.Prepend("lowFcn_");
    fcnLow = (TF1*)fin->Get(spname_low);

    TString spname_high = spname_core; spname_high.Prepend("highFcn_");
    fcnHigh = (TF1*)fin->Get(spname_high);

    TString spname_mid = spname_core;
    fcnMid = (TSpline3*)fin->Get(spname_mid);
  }
  else cerr << "MelaPConstant::GetFcnFromFile: Failed to open file in path " << path << endl;
  if (fin!=0 && fin->IsOpen()) fin->Close();
}

double MelaPConstant::Eval(MelaIO* RcdME, TVar::VerbosityLevel verbosity)const{
  if (verbosity>=TVar::DEBUG) cout << "Begin MelaPConstant::Eval" << endl;

  double result=1;
  if (RcdME->melaCand->genStatus==-1) return result;

  double candMass = RcdME->melaCand->m();
  if (verbosity>=TVar::DEBUG) cout << "MelaPConstant::Eval: Candidate mass is " << candMass << endl;
  if (candMass<=0.) return 0.;
  else if (fcnLow!=0 && candMass<fcnLow->GetXmax()) result = fcnLow->Eval(candMass);
  else if (fcnHigh!=0 && candMass>fcnHigh->GetXmin()) result = fcnHigh->Eval(candMass);
  else if (fcnMid!=0) result = fcnMid->Eval(candMass);
  else return result;

  if (verbosity>=TVar::DEBUG) cout << "MelaPConstant::Eval: Spline evaluated to " << result << endl;

  if (
    (processME==TVar::JHUGen || processME==TVar::MCFM)
    &&
    processProd==TVar::ZZGG
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();
    double propagator, mh, gah;
    RcdME->getHiggsMassWidth(mh, gah, 0);
    propagator = 1./(pow(pow(candMass, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "alphas(MZ)=" << alphasVal << "**2 "
      << "propagator=" << propagator << " "
      << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
      << endl;

    result *= pow(alphasVal, 2);
    result *= propagator;
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);

  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::ZZGG
    &&
    processProc==TVar::bkgZZ
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();
    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "alphas(MZ)=" << alphasVal << "**2 "
      << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
      << endl;

    result *= pow(alphasVal, 2);
    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);

  }
  else if (
    processME==TVar::MCFM
    &&
    processProd==TVar::ZZQQB
    &&
    processProc==TVar::bkgZZ
    ){
    result = pow(10., result);

    double aL1, aR1, aL2, aR2;
    RcdME->getVDaughterCouplings(aL1, aR1, 0);
    RcdME->getVDaughterCouplings(aL2, aR2, 1);

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "L**2+R**2 couplings=" << pow(aL1, 2)+pow(aR1, 2) << " " << pow(aL2, 2)+pow(aR2, 2) << " "
      << endl;

    if (fabs(aL1)>0. || fabs(aR1)>0.) result *= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) result *= pow(aL2, 2)+pow(aR2, 2);

  }
  else if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::JQCD
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "alphas(MZ)=" << alphasVal << "**3 "
      << endl;

    result *= pow(alphasVal, 3);
  }
  else if (
    processME==TVar::JHUGen
    &&
    processProd==TVar::JJQCD
    &&
    processProc==TVar::HSMHiggs
    ){
    result = pow(10., result);

    double alphasVal;
    alphasVal = RcdME->getAlphaSatMZ();

    if (verbosity>=TVar::DEBUG) cout
      << "MelaPConstant::Eval: Multiplying " << result << " by "
      << "alphas(MZ)=" << alphasVal << "**4 "
      << endl;

    result *= pow(alphasVal, 4);

  }
  else{
    result = pow(10., result);
    if (verbosity>=TVar::DEBUG) cout << "MelaPConstant::Eval: Multiplying " << result << " by nothing." << endl;
  }

  if (verbosity>=TVar::DEBUG) cout << "End MelaPConstant::Eval with " << result << endl;
  return result;
}



