#ifndef TCOUPLINGS_HH
#define TCOUPLINGS_HH

//---------------------------------
// Coupling array sizes
//---------------------------------
#define SIZE_HVV 39
#define SIZE_HGG 3
#define SIZE_ZQQ 2
#define SIZE_ZVV 2
#define SIZE_GQQ 2
#define SIZE_GGG 5
#define SIZE_GVV 10
#define SIZE_HVV_FREENORM 2
#define SIZE_HQQ 2

#include <iostream>
#include <vector>

using namespace std;

class SpinZeroCouplings{
public:
  SpinZeroCouplings(){ reset(); }

  void allow_WWZZSeparation(bool doAllow = true){ separateWWZZcouplings = doAllow; }
  void reset(){
    allow_WWZZSeparation(false);

    for (int im=0; im<2; im++){
      for (int ic=0; ic<SIZE_HVV; ic++){
        Hzzcoupl[ic][im] = 0;
        Hwwcoupl[ic][im] = 0;
        H2zzcoupl[ic][im] = 0;
        H2wwcoupl[ic][im] = 0;
      }
      for (int ic=0; ic<SIZE_HGG; ic++) Hggcoupl[ic][im]=0;
      for (int ic=0; ic<SIZE_HQQ; ic++) Hqqcoupl[ic][im]=0;
    }
    /*
    Hqqcoupl[0][0] = 1.0;
    Hggcoupl[0][0] = 1.0;
    Hzzcoupl[0][0] = 1.0;
    Hwwcoupl[0][0] = 1.0;
    */
    for (int ic=0; ic<SIZE_HVV_FREENORM; ic++) Hvvcoupl_freenorm[ic]=0;
    for (int ik=0; ik<3; ik++){
      HzzCLambda_qsq[ik]=0;
      HwwCLambda_qsq[ik]=0;
      H2zzCLambda_qsq[ik]=0;
      H2wwCLambda_qsq[ik]=0;
      for (int ic=0; ic<4; ic++){ // These default values do not matter as long as the c's are 0.
        HzzLambda_qsq[ic][ik] = 100.;
        HwwLambda_qsq[ic][ik] = 100.;
        H2zzLambda_qsq[ic][ik] = 100.;
        H2wwLambda_qsq[ic][ik] = 100.;
      }
    }
  };
  void copy(SpinZeroCouplings& other){
    allow_WWZZSeparation(other.separateWWZZcouplings);
    for (int im=0; im<2; im++){
      for (int ic=0; ic<SIZE_HVV; ic++){
        Hzzcoupl[ic][im] = (other.Hzzcoupl)[ic][im];
        Hwwcoupl[ic][im] = (other.Hwwcoupl)[ic][im];
      }
      for (int ic=0; ic<SIZE_HGG; ic++) Hggcoupl[ic][im]=(other.Hggcoupl)[ic][im];
      for (int ic=0; ic<SIZE_HQQ; ic++) Hqqcoupl[ic][im]=(other.Hqqcoupl)[ic][im];
    }
    for (int ic=0; ic<SIZE_HVV_FREENORM; ic++) Hvvcoupl_freenorm[ic]=(other.Hvvcoupl_freenorm)[ic];
    for (int ik=0; ik<3; ik++){
      HzzCLambda_qsq[ik]=(other.HzzCLambda_qsq)[ik];
      HwwCLambda_qsq[ik]=(other.HwwCLambda_qsq)[ik];;
      for (int ic=0; ic<4; ic++){
        HzzLambda_qsq[ic][ik] = (other.HzzLambda_qsq)[ic][ik];
        HwwLambda_qsq[ic][ik] = (other.HwwLambda_qsq)[ic][ik];
      }
    }
  };
  SpinZeroCouplings* getRef(){ return this; }

  void SetHVVFreeNormCouplings(unsigned int index, double cval){
    if (index>=SIZE_HVV_FREENORM) cerr << "Cannot set index " << index <<  " for the freenorm coupling, out of range." << endl;
    else Hvvcoupl_freenorm[index] = cval;
  };
  void SetHVVCouplings(unsigned int index, double c_real, double c_imag, bool setWW = false, int whichResonance=1){
    if (!separateWWZZcouplings && setWW) return;
    if (index>=SIZE_HVV){ cerr << "Cannot set index " << index << ", out of range for the type requested." << endl; }
    else if (whichResonance<0 || whichResonance>2) cerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << endl;
    else{
      if (whichResonance==1){ // First resonance
        if (setWW){
          Hwwcoupl[index][0] = c_real;
          Hwwcoupl[index][1] = c_imag;
        }
        else{
          Hzzcoupl[index][0] = c_real;
          Hzzcoupl[index][1] = c_imag;
        }
      }
      else{ // Second resonance
        if (setWW){
          H2wwcoupl[index][0] = c_real;
          H2wwcoupl[index][1] = c_imag;
        }
        else{
          H2zzcoupl[index][0] = c_real;
          H2zzcoupl[index][1] = c_imag;
        }
      }
    }
  };
  void SetHVVLambdaQ2(unsigned int gType, unsigned int index, double lambda, bool setWW = false, int whichResonance=1){
    if (!separateWWZZcouplings && setWW) return;
    if (index>=3 || gType>=4) cerr << "Cannot set index " << index <<  " for g" << (gType+1) << "_dyn, out of range." << endl;
    else if (whichResonance<0 || whichResonance>2) cerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << endl;
    else{
      if (whichResonance==1){
        if (setWW) HwwLambda_qsq[gType][index] = lambda;
        else HzzLambda_qsq[gType][index] = lambda;
      }
      else{
        if (setWW) H2wwLambda_qsq[gType][index] = lambda;
        else H2zzLambda_qsq[gType][index] = lambda;
      }
    }
  };
  void SetHVVSignCQ2(unsigned int index, int csign, bool setWW = false, int whichResonance=1){
    if (!separateWWZZcouplings && setWW) return;
    if (index>=3) cerr << "Cannot set index " << index << " for the c(z/w)qsq sign, out of range." << endl;
    else if (csign>1 || csign<-1) cerr << "Invalid csign argument. It has to be in the range [-1,1] with default to 0." << endl;
    else if (whichResonance<0 || whichResonance>2) cerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << endl;
    else{
      if (whichResonance==1){
        if (setWW) HwwCLambda_qsq[index] = csign;
        else HzzCLambda_qsq[index] = csign;
      }
      else{
        if (setWW) H2wwCLambda_qsq[index] = csign;
        else H2zzCLambda_qsq[index] = csign;
      }
    }
  };
  void SetHGGCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_HGG) cerr << "Cannot set index " << index << " for Hggcoupl, out of range for the type requested." << endl;
    else{
      Hggcoupl[index][0] = c_real;
      Hggcoupl[index][1] = c_imag;
    }
  };
  void SetHQQCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_HQQ) cerr << "Cannot set index " << index << " for Hqqcoupl, out of range for the type requested." << endl;
    else{
      Hqqcoupl[index][0] = c_real;
      Hqqcoupl[index][1] = c_imag;
    }
  };

  double Hqqcoupl[SIZE_HQQ][2];
  double Hggcoupl[SIZE_HGG][2];
  double Hvvcoupl_freenorm[SIZE_HVV_FREENORM];

  double Hzzcoupl[SIZE_HVV][2];
  double Hwwcoupl[SIZE_HVV][2];
  double HzzLambda_qsq[4][3];
  double HwwLambda_qsq[4][3];
  int HzzCLambda_qsq[3];
  int HwwCLambda_qsq[3];

  double H2zzcoupl[SIZE_HVV][2];
  double H2wwcoupl[SIZE_HVV][2];
  double H2zzLambda_qsq[4][3];
  double H2wwLambda_qsq[4][3];
  int H2zzCLambda_qsq[3];
  int H2wwCLambda_qsq[3];

  bool separateWWZZcouplings;

  inline virtual ~SpinZeroCouplings(){};
  ClassDef(SpinZeroCouplings, 0)
};

class SpinOneCouplings{
public:
  SpinOneCouplings(){ reset(); }

  void reset(){
    for (int im=0; im<2; im++){
      for (int ic=0; ic<SIZE_ZVV; ic++) Zvvcoupl[ic][im] = 0;
      for (int ic=0; ic<SIZE_ZQQ; ic++) Zqqcoupl[ic][im] = 0;
    }
    /*
    Zvvcoupl[0][0]=1.0;
    Zqqcoupl[0][0]=1.0;
    Zqqcoupl[1][0]=1.0;
    */
  };
  void copy(SpinOneCouplings& other){
    for (int im=0; im<2; im++){
      for (int ic=0; ic<SIZE_ZVV; ic++) Zvvcoupl[ic][im] = (other.Zvvcoupl)[ic][im];
      for (int ic=0; ic<SIZE_ZQQ; ic++) Zqqcoupl[ic][im] = (other.Zqqcoupl)[ic][im];
    }
  };
  SpinOneCouplings* getRef(){ return this; }

  void SetZVVCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_ZVV) cerr << "Cannot set index " << index << " for the Zvvcoupl, out of range." << endl;
    else{
      Zvvcoupl[index][0] = c_real;
      Zvvcoupl[index][1] = c_imag;
    }
  };
  void SetZQQCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_ZQQ) cerr << "Cannot set index " << index << " for the Zqqcoupl, out of range." << endl;
    else{
      Zqqcoupl[index][0] = c_real;
      Zqqcoupl[index][1] = c_imag;
    }
  };

  double Zvvcoupl[SIZE_ZVV][2];
  double Zqqcoupl[SIZE_ZQQ][2];

  inline virtual ~SpinOneCouplings(){};
  ClassDef(SpinOneCouplings, 0)
};

class SpinTwoCouplings{
public:
  SpinTwoCouplings(){ reset(); }

  void reset(){
    for (int im=0; im<2; im++){
      for (int ic=0; ic<SIZE_GVV; ic++) Gvvcoupl[ic][im] = 0;
      for (int ic=0; ic<SIZE_GGG; ic++) Gggcoupl[ic][im] = 0;
      for (int ic=0; ic<SIZE_GQQ; ic++) Gqqcoupl[ic][im] = 0;
    }
    /*
    Gggcoupl[0][0]=1.0;
    Gqqcoupl[0][0]=1.0;
    Gqqcoupl[1][0]=1.0;
    Gvvcoupl[0][0]=1.0;
    */
  };
  void copy(SpinTwoCouplings& other){
    for (int im=0; im<2; im++){
      for (int ic=0; ic<SIZE_GVV; ic++) Gvvcoupl[ic][im] = (other.Gvvcoupl)[ic][im];
      for (int ic=0; ic<SIZE_GGG; ic++) Gggcoupl[ic][im] = (other.Gggcoupl)[ic][im];
      for (int ic=0; ic<SIZE_GQQ; ic++) Gqqcoupl[ic][im] = (other.Gqqcoupl)[ic][im];
    }
  };
  SpinTwoCouplings* getRef(){ return this; }

  void SetGVVCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_GVV) cerr << "Cannot set index " << index << " for the Gvvcoupl, out of range." << endl;
    else{
      Gvvcoupl[index][0] = c_real;
      Gvvcoupl[index][1] = c_imag;
    }
  };
  void SetGQQCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_GQQ) cerr << "Cannot set index " << index << " for the Gqqcoupl, out of range." << endl;
    else{
      Gqqcoupl[index][0] = c_real;
      Gqqcoupl[index][1] = c_imag;
    }
  };
  void SetGGGCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_GGG) cerr << "Cannot set index " << index << " for the Gggcoupl, out of range." << endl;
    else{
      Gggcoupl[index][0] = c_real;
      Gggcoupl[index][1] = c_imag;
    }
  };

  double Gvvcoupl[SIZE_GVV][2];
  double Gqqcoupl[SIZE_GQQ][2];
  double Gggcoupl[SIZE_GGG][2];

  inline virtual ~SpinTwoCouplings(){};
  ClassDef(SpinTwoCouplings, 0)
};

ClassImp(SpinZeroCouplings)
ClassImp(SpinOneCouplings)
ClassImp(SpinTwoCouplings)


#endif
