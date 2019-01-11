#ifndef MELAINPUTOUTPUT_H
#define MELAINPUTOUTPUT_H

#include "TMCFM.hh"
#include "MELACandidate.h"


class MelaIO{
protected:

  double partonWeight[2][nmsq];
  double MEsq[nmsq][nmsq];
  double weightedMEsq[nmsq][nmsq];
  double sumME;
  double MEConst;

  double Qren;
  double Qfac;
  double alphas_mz;
  double alphas_Qren;

  double MH_GaH[nSupportedHiggses][2];

  double VDau1coupl[2]; // L/R
  double VDau2coupl[2]; // L/R

public:

  MELACandidate* melaCand; // Persistent container of the four-vectors, not owned by MelaIO

  void reset(){
    sumME=0;
    MEConst=1;
    for (int ix=0; ix<nmsq; ix++){
      for (int pp=0; pp<2; pp++) partonWeight[pp][ix]=0;
      for (int iy=0; iy<nmsq; iy++){
        MEsq[ix][iy]=0;
        weightedMEsq[ix][iy]=0;
      }
    }
    Qren=0;
    Qfac=0;
    alphas_mz=0;
    alphas_Qren=0;
    for (unsigned int jh=0; jh<nSupportedHiggses; jh++){ for (unsigned int img=0; img<2; img++) MH_GaH[jh][img]=0; }
    for (unsigned int ic=0; ic<2; ic++){
      VDau1coupl[ic]=0;
      VDau2coupl[ic]=0;
    }
  }
  MelaIO* getRef(){ return this; }

  // ME-related I/O
  void setPartonWeights(
    double partonOneWeight_[nmsq],
    double partonTwoWeight_[nmsq]
    ){
    for (int ix=0; ix<nmsq; ix++){
      partonWeight[0][ix]=partonOneWeight_[ix];
      partonWeight[1][ix]=partonTwoWeight_[ix];
    }
  }
  void setMEArray(double MEsq_[nmsq][nmsq], bool transpose=false){
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++){
        int jx=(transpose ? iy : ix);
        int jy=(transpose ? ix : iy);
        MEsq[ix][iy]=MEsq_[jx][jy];
      }
    }
  }
  void addMEArray(double MEsq_[nmsq][nmsq], double factor=1., bool transpose=false){
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++){
        int jx=(transpose ? iy : ix);
        int jy=(transpose ? ix : iy);
        MEsq[ix][iy]+=MEsq_[jx][jy]*factor;
      }
    }
    if (sumME!=0) computeWeightedMEArray();
  }
  void addMERecord(MelaIO* rcd, double factor=1., bool overwrite=false){
    double MEsq_[nmsq][nmsq]={ { 0 } };
    double partonOneWeight_[nmsq]={ 0 };
    double partonTwoWeight_[nmsq]={ 0 };
    rcd->getUnweightedMEArray(MEsq_);
    rcd->getPartonWeights(partonOneWeight_, partonTwoWeight_);

    if (overwrite) reset();
    setPartonWeights(partonOneWeight_, partonTwoWeight_);
    addMEArray(MEsq_, factor);
  }

  void computeWeightedMEArray(){
    sumME=0;
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++){
        weightedMEsq[ix][iy]=partonWeight[0][ix]*MEsq[ix][iy]*partonWeight[1][iy];
        sumME += weightedMEsq[ix][iy];
      }
    }
  }

  MelaIO(){ melaCand=0; reset(); }
  virtual ~MelaIO(){};

  double getSumME()const{ return sumME; }
  void setMEConst(const double& val){ MEConst=val; }
  void setMEConst(const float& val){ MEConst=(double)val; }
  double getMEConst()const{ return MEConst; }
  void getWeightedMEArray(double MEsq_[nmsq][nmsq]){
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++) MEsq_[ix][iy] = weightedMEsq[ix][iy];
    }
  }
  void getUnweightedMEArray(double MEsq_[nmsq][nmsq])const{
    for (int ix=0; ix<nmsq; ix++){
      for (int iy=0; iy<nmsq; iy++) MEsq_[ix][iy] = MEsq[ix][iy];
    }
  }
  void getPartonWeights(
    double partonOneWeight_[nmsq],
    double partonTwoWeight_[nmsq]
    )const{
    for (int ix=0; ix<nmsq; ix++){
      partonOneWeight_[ix] = partonWeight[0][ix];
      partonTwoWeight_[ix] = partonWeight[1][ix];
    }
  }


  // Scale-related I/O
  void setRenormalizationScale(const double& val){ Qren=val; }
  double getRenormalizationScale()const{ return Qren; }
  void setFactorizationScale(const double& val){ Qfac=val; }
  double getFactorizationScale()const{ return Qfac; }
  void setAlphaS(const double& val){ alphas_Qren=val; }
  double getAlphaS()const{ return alphas_Qren; }
  void setAlphaSatMZ(const double& val){ alphas_mz=val; }
  double getAlphaSatMZ()const{ return alphas_mz; }

  // Mass-related I/O
  void setHiggsMassWidth(const double& mass_, const double& width_, int jh){ if (jh<nSupportedHiggses){ MH_GaH[jh][0]=mass_; MH_GaH[jh][1]=width_; } }
  void getHiggsMassWidth(double& mass_, double& width_, int jh)const{ if (jh<nSupportedHiggses){ mass_=MH_GaH[jh][0]; width_=MH_GaH[jh][1]; } }

  // Z, W or gamma couplings to fermions
  void setVDaughterCouplings(const double& left, const double& right, int iv){
    if (iv==0){ VDau1coupl[0]=left; VDau1coupl[1]=right; }
    else if (iv==1){ VDau2coupl[0]=left; VDau2coupl[1]=right; }
  }
  void getVDaughterCouplings(double& left, double& right, int iv)const{
    if (iv==0){ left=VDau1coupl[0]; right=VDau1coupl[1]; }
    else if (iv==1){ left=VDau2coupl[0]; right=VDau2coupl[1]; }
  }


  ClassDef(MelaIO, 0)
};

ClassImp(MelaIO)

#endif
