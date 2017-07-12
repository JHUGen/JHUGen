#ifndef TCOUPLINGS_HH
#define TCOUPLINGS_HH

#include <iostream>
#include "TCouplingsBase.hh"


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
        Hzzpcoupl[ic][im] = 0;
        Hzpzpcoupl[ic][im] = 0;
        Hwwpcoupl[ic][im] = 0;
        Hwpwpcoupl[ic][im] = 0;
      }
      for (int ic=0; ic<SIZE_HGG; ic++){
        Hggcoupl[ic][im]=0;
        Hg4g4coupl[ic][im]=0;
        H2ggcoupl[ic][im]=0;
        H2g4g4coupl[ic][im]=0;
      }
      for (int ic=0; ic<SIZE_HQQ; ic++){
        Hqqcoupl[ic][im]=0;
        Httcoupl[ic][im]=0;
        Hbbcoupl[ic][im]=0;
        Ht4t4coupl[ic][im]=0;
        Hb4b4coupl[ic][im]=0;
        H2qqcoupl[ic][im]=0;
        H2ttcoupl[ic][im]=0;
        H2bbcoupl[ic][im]=0;
        H2t4t4coupl[ic][im]=0;
        H2b4b4coupl[ic][im]=0;
      }
      for (int ic = 0; ic < SIZE_Vp; ic++){
        Hzpcontact[ic][im]=0;
        Hwpcontact[ic][im]=0;
      }
    }
    for (int ik=0; ik<SIZE_HVV_CQSQ; ik++){
      HzzCLambda_qsq[ik]=0;
      HwwCLambda_qsq[ik]=0;
      H2zzCLambda_qsq[ik]=0;
      H2wwCLambda_qsq[ik]=0;
      for (int ic=0; ic<SIZE_HVV_LAMBDAQSQ; ic++){ // These default values do not matter as long as the c's are 0.
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
        Hzzpcoupl[ic][im] = (other.Hzzpcoupl)[ic][im];
        Hzpzpcoupl[ic][im] = (other.Hzpzpcoupl)[ic][im];
        Hwwpcoupl[ic][im] = (other.Hwwpcoupl)[ic][im];
        Hwpwpcoupl[ic][im] = (other.Hwpwpcoupl)[ic][im];
      }
      for (int ic=0; ic<SIZE_HGG; ic++){
        Hggcoupl[ic][im]=(other.Hggcoupl)[ic][im];
        Hg4g4coupl[ic][im]=(other.Hg4g4coupl)[ic][im];
        H2ggcoupl[ic][im]=(other.H2ggcoupl)[ic][im];
        H2g4g4coupl[ic][im]=(other.H2g4g4coupl)[ic][im];
      }
      for (int ic=0; ic<SIZE_HQQ; ic++){
        Hqqcoupl[ic][im]=(other.Hqqcoupl)[ic][im];
        Httcoupl[ic][im]=(other.Httcoupl)[ic][im];
        Hbbcoupl[ic][im]=(other.Hbbcoupl)[ic][im];
        Ht4t4coupl[ic][im]=(other.Ht4t4coupl)[ic][im];
        Hb4b4coupl[ic][im]=(other.Hb4b4coupl)[ic][im];
        H2qqcoupl[ic][im]=(other.H2qqcoupl)[ic][im];
        H2ttcoupl[ic][im]=(other.H2ttcoupl)[ic][im];
        H2bbcoupl[ic][im]=(other.H2bbcoupl)[ic][im];
        H2t4t4coupl[ic][im]=(other.H2t4t4coupl)[ic][im];
        H2b4b4coupl[ic][im]=(other.H2b4b4coupl)[ic][im];
      }
      for (int ic=0; ic<SIZE_Vp; ic++){
        Hzpcontact[ic][im] = (other.Hzpcontact)[ic][im];
        Hwpcontact[ic][im] = (other.Hwpcontact)[ic][im];
      }
    }
    for (int ik=0; ik<SIZE_HVV_CQSQ; ik++){
      HzzCLambda_qsq[ik]=(other.HzzCLambda_qsq)[ik];
      HwwCLambda_qsq[ik]=(other.HwwCLambda_qsq)[ik];;
      for (int ic=0; ic<SIZE_HVV_LAMBDAQSQ; ic++){
        HzzLambda_qsq[ic][ik] = (other.HzzLambda_qsq)[ic][ik];
        HwwLambda_qsq[ic][ik] = (other.HwwLambda_qsq)[ic][ik];
      }
    }
  };
  SpinZeroCouplings* getRef(){ return this; }

  void SetHVVCouplings(unsigned int index, double c_real, double c_imag, bool setWW = false, int whichResonance=1){
    if (!separateWWZZcouplings && setWW) return;
    if (index>=SIZE_HVV){ std::cerr << "Cannot set index " << index << ", out of range for the type requested." << std::endl; }
    else if (whichResonance<=0 || whichResonance>2) std::cerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << std::endl;
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
    if (index>=SIZE_HVV_CQSQ || gType>=SIZE_HVV_LAMBDAQSQ) std::cerr << "Cannot set index " << index <<  " for g" << (gType+1) << "_dyn, out of range." << std::endl;
    else if (whichResonance<=0 || whichResonance>2) std::cerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << std::endl;
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
    if (index>=SIZE_HVV_CQSQ) std::cerr << "Cannot set index " << index << " for the c(z/w)qsq sign, out of range." << std::endl;
    else if (csign>1 || csign<-1) std::cerr << "Invalid csign argument. It has to be in the range [-1,1] with default to 0." << std::endl;
    else if (whichResonance<=0 || whichResonance>2) std::cerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << std::endl;
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
  void SetHGGCouplings(unsigned int index, double c_real, double c_imag, int whichLoop=1, int whichResonance=1){
    if (index>=SIZE_HGG) std::cerr << "Cannot set index " << index << " for Hggcoupl, out of range for the type requested." << std::endl;
    else if (whichResonance<=0 || whichResonance>2) std::cerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << std::endl;
    else if (whichLoop<=0 || whichLoop>2) std::cerr << "gg loop " << whichLoop << " is not supported. Set it to 1 for the loop that corresponds to the top/bottom couplings, or 2 for the loop that corresponds to the tprime/bprime couplings." << std::endl;
    else{
      if (whichResonance==1){
        if (whichLoop==1){
          Hggcoupl[index][0] = c_real;
          Hggcoupl[index][1] = c_imag;
        }
        else{
          Hg4g4coupl[index][0] = c_real;
          Hg4g4coupl[index][1] = c_imag;
        }
      }
      else{
        if (whichLoop==1){
          H2ggcoupl[index][0] = c_real;
          H2ggcoupl[index][1] = c_imag;
        }
        else{
          H2g4g4coupl[index][0] = c_real;
          H2g4g4coupl[index][1] = c_imag;
        }
      }
    }
  };
  void SetHQQCouplings(unsigned int index, double c_real, double c_imag, int qid=0, int whichResonance=1){
    if (index>=SIZE_HQQ) std::cerr << "Cannot set index " << index << " for Hqqcoupl, out of range for the type requested." << std::endl;
    else if (whichResonance<=0 || whichResonance>2) std::cerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << std::endl;
    else if (abs(qid)>8) std::cerr << "Quark id=" << qid << ">8 is not supported. Please change the id to 0-4 (qq), or one of 5 bottom), 6 (top), 7 (bprime), 8 (tprime)." << std::endl;
    else{
      if (whichResonance==1){
        if (abs(qid)<5){
          Hqqcoupl[index][0] = c_real;
          Hqqcoupl[index][1] = c_imag;
        }
        else if (abs(qid)==5){
          Hbbcoupl[index][0] = c_real;
          Hbbcoupl[index][1] = c_imag;
        }
        else if (abs(qid)==6){
          Httcoupl[index][0] = c_real;
          Httcoupl[index][1] = c_imag;
        }
        else if (abs(qid)==7){
          Hb4b4coupl[index][0] = c_real;
          Hb4b4coupl[index][1] = c_imag;
        }
        else if (abs(qid)==8){
          Ht4t4coupl[index][0] = c_real;
          Ht4t4coupl[index][1] = c_imag;
        }
      }
      else{
        if (abs(qid)<5){
          H2qqcoupl[index][0] = c_real;
          H2qqcoupl[index][1] = c_imag;
        }
        else if (abs(qid)==5){
          H2bbcoupl[index][0] = c_real;
          H2bbcoupl[index][1] = c_imag;
        }
        else if (abs(qid)==6){
          H2ttcoupl[index][0] = c_real;
          H2ttcoupl[index][1] = c_imag;
        }
        else if (abs(qid)==7){
          H2b4b4coupl[index][0] = c_real;
          H2b4b4coupl[index][1] = c_imag;
        }
        else if (abs(qid)==8){
          H2t4t4coupl[index][0] = c_real;
          H2t4t4coupl[index][1] = c_imag;
        }
      }
    }
  };

  void SetHVVpCouplings(unsigned int index, double c_real, double c_imag, bool setWWp = false, int whichResonance=1){
    if (!separateWWZZcouplings && setWWp) return;
    if (index>=SIZE_HVV){ std::cerr << "Cannot set index " << index << ", out of range for the type requested." << std::endl; }
    else if (whichResonance!=1) {std::cerr << "Contact terms are only for the first resonance" << std::endl;}
    else{
      if (setWWp){
        Hwwpcoupl[index][0] = c_real;
        Hwwpcoupl[index][1] = c_imag;
      }
      else{
        Hzzpcoupl[index][0] = c_real;
        Hzzpcoupl[index][1] = c_imag;
      }
    }
  };

  void SetHVpVpCouplings(unsigned int index, double c_real, double c_imag, bool setWpWp = false, int whichResonance=1){
    if (!separateWWZZcouplings && setWpWp) return;
    if (index>=SIZE_HVV){ std::cerr << "Cannot set index " << index << ", out of range for the type requested." << std::endl; }
    else if (whichResonance!=1) {std::cerr << "Contact terms are only for the first resonance" << std::endl;}
    else{
      if (setWpWp){
        Hwwpcoupl[index][0] = c_real;
        Hwwpcoupl[index][1] = c_imag;
      }
      else{
        Hzzpcoupl[index][0] = c_real;
        Hzzpcoupl[index][1] = c_imag;
      }
    }
  };

  void SetZpcontactTerms(unsigned int index, double c_real, double c_imag, int whichResonance=1){
    if (whichResonance!=1) {std::cerr << "Contact terms are only for the first resonance" << std::endl;}
    else if (index > SIZE_Vp) {
      std::cerr << "index too big for SetZpcontactTerms: " << index << std::endl;
    }
    else {
      Hzpcontact[index][0] = c_real;
      Hzpcontact[index][1] = c_imag;
    }
  }

  void SetWpcontactTerms(unsigned int index, double c_real, double c_imag, int whichResonance=1){
    if (whichResonance!=1) {std::cerr << "Contact terms are only for the first resonance" << std::endl;}
    else if (index > SIZE_Vp) {
      std::cerr << "index too big for SetWpcontactTerms: " << index << std::endl;
    }
    else if (
             (   index == gHIGGS_Vp_L_N || index == gHIGGS_Vp_R_N
              || index == gHIGGS_Vp_L_D || index == gHIGGS_Vp_R_D
              || index == gHIGGS_Vp_L_S || index == gHIGGS_Vp_R_S
              || index == gHIGGS_Vp_L_B || index == gHIGGS_Vp_R_B
             ) && (c_real || c_imag)
            ) {
      std::cerr << "no W' contact terms for neutrino, down, strange, or bottom (use the lepton or up instead)" << std::endl;
    }
    else {
      Hwpcontact[index][0] = c_real;
      Hwpcontact[index][1] = c_imag;
    }
  }

  void SetUseVprime(bool useVp, double mass, double width) {
    UseVprime = useVp;
    M_Vprime = mass;
    Ga_Vprime = width;
  }

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
  double Hzpcontact[SIZE_Vp][2];
  double Hwwpcoupl[SIZE_HVV][2];
  double Hwpwpcoupl[SIZE_HVV][2];
  double Hwpcontact[SIZE_Vp][2];

  bool separateWWZZcouplings;
  bool UseVprime;
  double M_Vprime, Ga_Vprime;

  inline virtual ~SpinZeroCouplings(){};
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
    if (index>=SIZE_ZVV) std::cerr << "Cannot set index " << index << " for the Zvvcoupl, out of range." << std::endl;
    else{
      Zvvcoupl[index][0] = c_real;
      Zvvcoupl[index][1] = c_imag;
    }
  };
  void SetZQQCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_ZQQ) std::cerr << "Cannot set index " << index << " for the Zqqcoupl, out of range." << std::endl;
    else{
      Zqqcoupl[index][0] = c_real;
      Zqqcoupl[index][1] = c_imag;
    }
  };

  double Zvvcoupl[SIZE_ZVV][2];
  double Zqqcoupl[SIZE_ZQQ][2];

  inline virtual ~SpinOneCouplings(){};
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
    if (index>=SIZE_GVV) std::cerr << "Cannot set index " << index << " for the Gvvcoupl, out of range." << std::endl;
    else{
      Gvvcoupl[index][0] = c_real;
      Gvvcoupl[index][1] = c_imag;
    }
  };
  void SetGQQCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_GQQ) std::cerr << "Cannot set index " << index << " for the Gqqcoupl, out of range." << std::endl;
    else{
      Gqqcoupl[index][0] = c_real;
      Gqqcoupl[index][1] = c_imag;
    }
  };
  void SetGGGCouplings(unsigned int index, double c_real, double c_imag){
    if (index>=SIZE_GGG) std::cerr << "Cannot set index " << index << " for the Gggcoupl, out of range." << std::endl;
    else{
      Gggcoupl[index][0] = c_real;
      Gggcoupl[index][1] = c_imag;
    }
  };

  double Gvvcoupl[SIZE_GVV][2];
  double Gqqcoupl[SIZE_GQQ][2];
  double Gggcoupl[SIZE_GGG][2];

  inline virtual ~SpinTwoCouplings(){};
};


#endif
