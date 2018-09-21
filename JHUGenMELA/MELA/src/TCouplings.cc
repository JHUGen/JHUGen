#include <iostream>
#include "MELAStreamHelpers.hh"
#include "TCouplings.hh"


using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;
using namespace std;


/********** Spin-0 couplings **********/
SpinZeroCouplings::SpinZeroCouplings(){ reset(); }
SpinZeroCouplings::~SpinZeroCouplings(){}

void SpinZeroCouplings::allow_WWZZSeparation(bool doAllow){ separateWWZZcouplings = doAllow; }
void SpinZeroCouplings::reset(){
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
}
void SpinZeroCouplings::copy(SpinZeroCouplings& other){
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
  }
  for (int ik=0; ik<SIZE_HVV_CQSQ; ik++){
    HzzCLambda_qsq[ik]=(other.HzzCLambda_qsq)[ik];
    HwwCLambda_qsq[ik]=(other.HwwCLambda_qsq)[ik];;
    for (int ic=0; ic<SIZE_HVV_LAMBDAQSQ; ic++){
      HzzLambda_qsq[ic][ik] = (other.HzzLambda_qsq)[ic][ik];
      HwwLambda_qsq[ic][ik] = (other.HwwLambda_qsq)[ic][ik];
    }
  }
}
SpinZeroCouplings* SpinZeroCouplings::getRef(){ return this; }

void SpinZeroCouplings::SetHVVCouplings(unsigned int index, double c_real, double c_imag, bool setWW, int whichResonance){
  if (!separateWWZZcouplings && setWW) return;
  if (index>=SIZE_HVV){ MELAerr << "Cannot set index " << index << ", out of range for the type requested." << endl; }
  else if (whichResonance<=0 || whichResonance>2) MELAerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << endl;
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
}
void SpinZeroCouplings::SetHVVLambdaQ2(unsigned int gType, unsigned int index, double lambda, bool setWW, int whichResonance){
  if (!separateWWZZcouplings && setWW) return;
  if (index>=SIZE_HVV_CQSQ || gType>=SIZE_HVV_LAMBDAQSQ) MELAerr << "Cannot set index " << index <<  " for g" << (gType+1) << "_dyn, out of range." << endl;
  else if (whichResonance<=0 || whichResonance>2) MELAerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << endl;
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
}
void SpinZeroCouplings::SetHVVSignCQ2(unsigned int index, int csign, bool setWW, int whichResonance){
  if (!separateWWZZcouplings && setWW) return;
  if (index>=SIZE_HVV_CQSQ) MELAerr << "Cannot set index " << index << " for the c(z/w)qsq sign, out of range." << endl;
  else if (csign>1 || csign<-1) MELAerr << "Invalid csign argument. It has to be in the range [-1,1] with default to 0." << endl;
  else if (whichResonance<=0 || whichResonance>2) MELAerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << endl;
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
}
void SpinZeroCouplings::SetHGGCouplings(unsigned int index, double c_real, double c_imag, int whichLoop, int whichResonance){
  if (index>=SIZE_HGG) MELAerr << "Cannot set index " << index << " for Hggcoupl, out of range for the type requested." << endl;
  else if (whichResonance<=0 || whichResonance>2) MELAerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << endl;
  else if (whichLoop<=0 || whichLoop>2) MELAerr << "gg loop " << whichLoop << " is not supported. Set it to 1 for the loop that corresponds to the top/bottom couplings, or 2 for the loop that corresponds to the tprime/bprime couplings." << endl;
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
}
void SpinZeroCouplings::SetHQQCouplings(unsigned int index, double c_real, double c_imag, int qid, int whichResonance){
  if (index>=SIZE_HQQ) MELAerr << "Cannot set index " << index << " for Hqqcoupl, out of range for the type requested." << endl;
  else if (whichResonance<=0 || whichResonance>2) MELAerr << "Resonance " << whichResonance << " is not supported. Set it to 1 for the regular Higgs and 2 for the high-mass resonance." << endl;
  else if (abs(qid)>8) MELAerr << "Quark id=" << qid << ">8 is not supported. Please change the id to 0-4 (qq), or one of 5 bottom), 6 (top), 7 (bprime), 8 (tprime)." << endl;
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
}

void SpinZeroCouplings::SetHVVpCouplings(unsigned int index, double c_real, double c_imag, bool setWWp, int whichResonance){
  if (index>=SIZE_HVV){ MELAerr << "Cannot set index " << index << ", out of range for the type requested." << endl; }
  else if (whichResonance!=1) {MELAerr << "Contact terms are only for the first resonance" << endl;}
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
}
void SpinZeroCouplings::SetHVpVpCouplings(unsigned int index, double c_real, double c_imag, bool setWpWp, int whichResonance){
  if (whichResonance!=1) { MELAerr << "Contact terms are only for the first resonance" << endl; }
  else if (index>=SIZE_HVV){ MELAerr << "Cannot set index " << index << ", out of range for the type requested." << endl; }
  else{
    if (setWpWp){
      Hwpwpcoupl[index][0] = c_real;
      Hwpwpcoupl[index][1] = c_imag;
    }
    else{
      Hzpzpcoupl[index][0] = c_real;
      Hzpzpcoupl[index][1] = c_imag;
    }
  }
}


/********** Spin-1 couplings **********/
SpinOneCouplings::SpinOneCouplings(){ reset(); }
SpinOneCouplings::~SpinOneCouplings(){}

void SpinOneCouplings::reset(){
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_ZVV; ic++) Zvvcoupl[ic][im] = 0;
    for (int ic=0; ic<SIZE_ZQQ; ic++) Zqqcoupl[ic][im] = 0;
  }
  /*
  Zvvcoupl[0][0]=1.0;
  Zqqcoupl[0][0]=1.0;
  Zqqcoupl[1][0]=1.0;
  */
}
void SpinOneCouplings::copy(SpinOneCouplings& other){
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_ZVV; ic++) Zvvcoupl[ic][im] = (other.Zvvcoupl)[ic][im];
    for (int ic=0; ic<SIZE_ZQQ; ic++) Zqqcoupl[ic][im] = (other.Zqqcoupl)[ic][im];
  }
}
SpinOneCouplings* SpinOneCouplings::getRef(){ return this; }

void SpinOneCouplings::SetZVVCouplings(unsigned int index, double c_real, double c_imag){
  if (index>=SIZE_ZVV) MELAerr << "Cannot set index " << index << " for the Zvvcoupl, out of range." << endl;
  else{
    Zvvcoupl[index][0] = c_real;
    Zvvcoupl[index][1] = c_imag;
  }
}
void SpinOneCouplings::SetZQQCouplings(unsigned int index, double c_real, double c_imag){
  if (index>=SIZE_ZQQ) MELAerr << "Cannot set index " << index << " for the Zqqcoupl, out of range." << endl;
  else{
    Zqqcoupl[index][0] = c_real;
    Zqqcoupl[index][1] = c_imag;
  }
}


/********** Spin-2 couplings **********/
SpinTwoCouplings::SpinTwoCouplings(){ reset(); }
SpinTwoCouplings::~SpinTwoCouplings(){}

void SpinTwoCouplings::reset(){
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_GVV; ic++){
      Gvvcoupl[ic][im] = 0;
      Gvvpcoupl[ic][im] = 0;
      Gvpvpcoupl[ic][im] = 0;
    }
    for (int ic=0; ic<SIZE_GGG; ic++) Gggcoupl[ic][im] = 0;
    for (int ic=0; ic<SIZE_GQQ; ic++) Gqqcoupl[ic][im] = 0;
  }
  /*
  Gggcoupl[0][0]=1.0;
  Gqqcoupl[0][0]=1.0;
  Gqqcoupl[1][0]=1.0;
  Gvvcoupl[0][0]=1.0;
  */
}
void SpinTwoCouplings::copy(SpinTwoCouplings& other){
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_GVV; ic++) Gvvcoupl[ic][im] = (other.Gvvcoupl)[ic][im];
    for (int ic=0; ic<SIZE_GGG; ic++) Gggcoupl[ic][im] = (other.Gggcoupl)[ic][im];
    for (int ic=0; ic<SIZE_GQQ; ic++) Gqqcoupl[ic][im] = (other.Gqqcoupl)[ic][im];
  }
}
SpinTwoCouplings* SpinTwoCouplings::getRef(){ return this; }

void SpinTwoCouplings::SetGVVCouplings(unsigned int index, double c_real, double c_imag){
  if (index>=SIZE_GVV) MELAerr << "Cannot set index " << index << " for the Gvvcoupl, out of range." << endl;
  else{
    Gvvcoupl[index][0] = c_real;
    Gvvcoupl[index][1] = c_imag;
  }
}
void SpinTwoCouplings::SetGVVpCouplings(unsigned int index, double c_real, double c_imag){
  if (index>=SIZE_GVV) MELAerr << "Cannot set index " << index << " for the Gvvpcoupl, out of range." << endl;
  else{
    Gvvpcoupl[index][0] = c_real;
    Gvvpcoupl[index][1] = c_imag;
  }
}
void SpinTwoCouplings::SetGVpVpCouplings(unsigned int index, double c_real, double c_imag){
  if (index>=SIZE_GVV) MELAerr << "Cannot set index " << index << " for the Gvpvpcoupl, out of range." << endl;
  else{
    Gvpvpcoupl[index][0] = c_real;
    Gvpvpcoupl[index][1] = c_imag;
  }
}
void SpinTwoCouplings::SetGQQCouplings(unsigned int index, double c_real, double c_imag){
  if (index>=SIZE_GQQ) MELAerr << "Cannot set index " << index << " for the Gqqcoupl, out of range." << endl;
  else{
    Gqqcoupl[index][0] = c_real;
    Gqqcoupl[index][1] = c_imag;
  }
}
void SpinTwoCouplings::SetGGGCouplings(unsigned int index, double c_real, double c_imag){
  if (index>=SIZE_GGG) MELAerr << "Cannot set index " << index << " for the Gggcoupl, out of range." << endl;
  else{
    Gggcoupl[index][0] = c_real;
    Gggcoupl[index][1] = c_imag;
  }
}


/********** Vprime couplings **********/
VprimeCouplings::VprimeCouplings(){ reset(); }
VprimeCouplings::~VprimeCouplings(){}

void VprimeCouplings::reset(){
  for (int im=0; im<2; im++){
    for (int ic = 0; ic < SIZE_Vpff; ic++){
      Zpffcoupl[ic][im]=0;
      Wpffcoupl[ic][im]=0;
    }
  }

  SetZPrimeMassWidth(-1, 0);
  SetWPrimeMassWidth(-1, 0);
}
void VprimeCouplings::copy(VprimeCouplings& other){
  for (int im=0; im<2; im++){
    for (int ic=0; ic<SIZE_Vpff; ic++){
      Zpffcoupl[ic][im] = (other.Zpffcoupl)[ic][im];
      Wpffcoupl[ic][im] = (other.Wpffcoupl)[ic][im];
    }
  }

  M_Zprime = other.M_Zprime;
  Ga_Zprime = other.Ga_Zprime;
  M_Wprime = other.M_Wprime;
  Ga_Wprime = other.Ga_Wprime;
}
VprimeCouplings* VprimeCouplings::getRef(){ return this; }

void VprimeCouplings::SetVpffCouplings(unsigned int index, double c_real, double c_imag, bool setWpff, int whichResonance){
  if (whichResonance!=1){ MELAerr << "VprimeCouplings::SetVpffCouplings: Contact terms are only for the first resonance" << endl; }
  else if (index>=SIZE_Vpff){ MELAerr << "VprimeCouplings::SetVpffCouplings: Cannot set index " << index << ", out of range for the type requested." << endl; }
  else{
    if (!setWpff){
      Zpffcoupl[index][0] = c_real;
      Zpffcoupl[index][1] = c_imag;
    }
    else{
      if (
        (index == gHIGGS_Vp_NuE_left || index == gHIGGS_Vp_NuE_right
         || index == gHIGGS_Vp_Dn_left || index == gHIGGS_Vp_Dn_right
         || index == gHIGGS_Vp_Str_left || index == gHIGGS_Vp_Str_right
         || index == gHIGGS_Vp_Bot_left || index == gHIGGS_Vp_Bot_right
         ) && (c_real!=0. || c_imag!=0.)
        ) MELAerr << "No W' contact terms for neutrino, down, strange, or bottom; use the lepton or up-quark versions instead!" << endl;
      else{
        Wpffcoupl[index][0] = c_real;
        Wpffcoupl[index][1] = c_imag;
      }
    }
  }
}
void VprimeCouplings::SetZPrimeMassWidth(double inmass, double inwidth){ M_Zprime = inmass; Ga_Zprime = inwidth; }
void VprimeCouplings::SetWPrimeMassWidth(double inmass, double inwidth){ M_Wprime = inmass; Ga_Wprime = inwidth; }
