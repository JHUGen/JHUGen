#ifndef MELAPCONSTANT_H
#define MELAPCONSTANT_H

#include "TUtil.hh"
#include "TF1.h"
#include "TSpline.h"


class MelaPConstant{
protected:

  TVar::MatrixElement processME;
  TVar::Production processProd;
  TVar::Process processProc;

  TF1* fcnLow;
  TF1* fcnHigh;
  TSpline3* fcnMid;
  TString fname;

  void GetFcnFromFile(const char* path, const char* spname);

  double GetHPropagator(const MelaIO* RcdME, const TVar::VerbosityLevel& verbosity)const;
  double GetZPropagator(const MelaIO* RcdME, const int& sigmaZ4lprop, const TVar::VerbosityLevel& verbosity)const;
  double GetAssociatedVjjPropagator(const MelaIO* RcdME, const TVar::VerbosityLevel& verbosity)const;
  double GetVDaughterCouplings(const MelaIO* RcdME, const TVar::VerbosityLevel& verbosity)const;
  double GetAlphaSatMZ(const MelaIO* RcdME, const unsigned int& powAlphaSatMZ, const TVar::VerbosityLevel& verbosity)const;

public:

  MelaPConstant(
    TVar::MatrixElement me_,
    TVar::Production prod_,
    TVar::Process proc_,
    const char* path,
    const char* spname
    );

  virtual ~MelaPConstant();

  double Eval(const MelaIO* RcdME, TVar::VerbosityLevel verbosity)const;

  bool IsValid(){ return bool(fcnMid!=0); }

  TString GetFileName(){ return fname; }
  TString GetSplineName(){ TString sname=""; if (IsValid()) sname=fcnMid->GetName(); return sname; }

};


#endif

