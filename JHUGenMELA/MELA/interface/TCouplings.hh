#ifndef TCOUPLINGS_HH
#define TCOUPLINGS_HH

#include "TCouplingsBase.hh"


class SpinZeroCouplings{
public:
  SpinZeroCouplings();
  SpinZeroCouplings(SpinZeroCouplings const& other);
  virtual ~SpinZeroCouplings();

  void allow_WWZZSeparation(bool doAllow = true);
  void reset();
  void copy(SpinZeroCouplings const& other);
  SpinZeroCouplings* getRef();

  void SetHVVCouplings(unsigned int index, double c_real, double c_imag, bool setWW = false, int whichResonance=1);
  void SetHVVLambdaQ2(unsigned int gType, unsigned int index, double lambda, bool setWW = false, int whichResonance=1);
  void SetHVVSignCQ2(unsigned int index, int csign, bool setWW = false, int whichResonance=1);
  void SetHGGCouplings(unsigned int index, double c_real, double c_imag, int whichLoop=1, int whichResonance=1);
  void SetHQQCouplings(unsigned int index, double c_real, double c_imag, int qid=0, int whichResonance=1);

  void SetHVVpCouplings(unsigned int index, double c_real, double c_imag, bool setWWp = false, int whichResonance=1);
  void SetHVpVpCouplings(unsigned int index, double c_real, double c_imag, bool setWpWp = false, int whichResonance=1);

  double Hggcoupl[SIZE_HGG][2];
  double Hqqcoupl[SIZE_HQQ][2];
  double Httcoupl[SIZE_HQQ][2];
  double Hbbcoupl[SIZE_HQQ][2];
  double Hg4g4coupl[SIZE_HGG][2];
  double Ht4t4coupl[SIZE_HQQ][2];
  double Hb4b4coupl[SIZE_HQQ][2];

  double Hzzcoupl[SIZE_HVV][2];
  double Hwwcoupl[SIZE_HVV][2];
  double HzzLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ];
  double HwwLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ];
  int HzzCLambda_qsq[SIZE_HVV_CQSQ];
  int HwwCLambda_qsq[SIZE_HVV_CQSQ];

  double H2ggcoupl[SIZE_HGG][2];
  double H2qqcoupl[SIZE_HQQ][2];
  double H2ttcoupl[SIZE_HQQ][2];
  double H2bbcoupl[SIZE_HQQ][2];
  double H2g4g4coupl[SIZE_HGG][2];
  double H2t4t4coupl[SIZE_HQQ][2];
  double H2b4b4coupl[SIZE_HQQ][2];

  double H2zzcoupl[SIZE_HVV][2];
  double H2wwcoupl[SIZE_HVV][2];
  double H2zzLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ];
  double H2wwLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ];
  int H2zzCLambda_qsq[SIZE_HVV_CQSQ];
  int H2wwCLambda_qsq[SIZE_HVV_CQSQ];

  double Hzzpcoupl[SIZE_HVV][2];
  double Hzpzpcoupl[SIZE_HVV][2];
  double Hwwpcoupl[SIZE_HVV][2];
  double Hwpwpcoupl[SIZE_HVV][2];

  bool separateWWZZcouplings;
};

class SpinOneCouplings{
public:
  SpinOneCouplings();
  SpinOneCouplings(SpinOneCouplings const& other);
  virtual ~SpinOneCouplings();

  void reset();
  void copy(SpinOneCouplings const& other);
  SpinOneCouplings* getRef();

  void SetZVVCouplings(unsigned int index, double c_real, double c_imag);
  void SetZQQCouplings(unsigned int index, double c_real, double c_imag);

  double Zvvcoupl[SIZE_ZVV][2];
  double Zqqcoupl[SIZE_ZQQ][2];
};

class SpinTwoCouplings{
public:
  SpinTwoCouplings();
  SpinTwoCouplings(SpinTwoCouplings const& other);
  virtual ~SpinTwoCouplings();

  void reset();
  void copy(SpinTwoCouplings const& other);
  SpinTwoCouplings* getRef();

  void SetGVVCouplings(unsigned int index, double c_real, double c_imag);
  void SetGVVpCouplings(unsigned int index, double c_real, double c_imag);
  void SetGVpVpCouplings(unsigned int index, double c_real, double c_imag);
  void SetGQQCouplings(unsigned int index, double c_real, double c_imag);
  void SetGGGCouplings(unsigned int index, double c_real, double c_imag);

  double Gvvcoupl[SIZE_GVV][2];
  double Gvvpcoupl[SIZE_GVV][2];
  double Gvpvpcoupl[SIZE_GVV][2];
  double Gqqcoupl[SIZE_GQQ][2];
  double Gggcoupl[SIZE_GGG][2];
};

class VprimeCouplings{
public:
  VprimeCouplings();
  VprimeCouplings(VprimeCouplings const& other);
  virtual ~VprimeCouplings();

  void reset();
  void copy(VprimeCouplings const& other);
  VprimeCouplings* getRef();

  void SetVpffCouplings(unsigned int index, double c_real, double c_imag, bool setWpff = false, int whichResonance=1);
  void SetZPrimeMassWidth(double inmass, double inwidth);
  void SetWPrimeMassWidth(double inmass, double inwidth);

  double Zpffcoupl[SIZE_Vpff][2];
  double Wpffcoupl[SIZE_Vpff][2];

  double M_Zprime;
  double Ga_Zprime;
  double M_Wprime;
  double Ga_Wprime;
};

class aTQGCCouplings{
public:
  aTQGCCouplings();
  aTQGCCouplings(aTQGCCouplings const& other);
  virtual ~aTQGCCouplings();

  void reset();
  void copy(aTQGCCouplings const& other);
  aTQGCCouplings* getRef();

  void SetATQGCCouplings(unsigned int index, double c_real, double c_imag);

  double aTQGCcoupl[SIZE_ATQGC][2];
};

#endif
