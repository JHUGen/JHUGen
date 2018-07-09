#ifndef MELAPCONSTANT_H
#define MELAPCONSTANT_H

#include "TUtil.hh"
#include "TF1.h"
#include "TSpline.h"


class MelaPConstant{
public:

  MelaPConstant(
    TVar::MatrixElement me_,
    TVar::Production prod_,
    TVar::Process proc_,
    const char* path,
    const char* spname
    );

  virtual ~MelaPConstant();

  double Eval(MelaIO* RcdME, TVar::VerbosityLevel verbosity)const;

private:

  TVar::MatrixElement processME;
  TVar::Production processProd;
  TVar::Process processProc;

  TF1* fcnLow;
  TF1* fcnHigh;
  TSpline3* fcnMid;

  void GetFcnFromFile(const char* path, const char* spname);

};


#endif

