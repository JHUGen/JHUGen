#ifndef EvtProb_VAR
#define EvtProb_VAR

#include <TLorentzVector.h>
#include "TH2F.h"
#include "TH1F.h"

#define EBEAM 4000.00
#define fbGeV2 0.389379E12
#define smallnumber 1e-15
#define sixteen_2Pi_to_8 3.88650230418250561e+07
#define   eight_2Pi_to_5 7.83410393050320417e+04
#define    four_2Pi_to_2 39.478417604357432
class TVar{
public:
  enum VerbosityLevel {
    ERROR = 0,
    INFO = 1,
    DEBUG = 2
  };
  enum MatrixElement{
    MCFM = 0,
    MadGraph = 1,
    JHUGen = 2,
    ANALYTICAL = 3
  };
  enum Production{
    GG = 0,
    QQB = 1,
    INDEPENDENT=2
  };

  enum Process{
    ZZ_2e2m  =0, // eemm
    HZZ_4l   =1, // 0+
    PSHZZ_4l =2, // 0-
    TZZ_4l =3,   // spin 2 couplings have to set in TEvtProb.cc
    VZZ_4l =4,   // spin 1 couplings have to set in TEvtProb.cc
    ZZ_4e    =5,
    GGZZ_4l = 6,
    AVZZ_4l = 7,
    QQB_TZZ_4l = 8,
    HDHZZ_4l = 9, // 0h+
    TZZ_DECAY_4l = 10,
    VZZ_DECAY_4l = 11,
    AVZZ_DECAY_4l = 12,
    PTZZ_2hminus_4l = 13, // 2h-
    TZZ_2hplus_4l = 14, // 2h+
    TZZ_2bplus_4l = 15, // 2b+
    SummedBackgrounds = 16, // SuperMela Background standin process.
    HZZ_4l_MIXCP = 17,
    HJJNONVBF = 18,
    HJJVBF = 19,
    GGZZTOT_4l = 20,
    GGZZINT_4l = 21,
    Null
  };
  enum LeptonFlavor{
    Flavor_Dummy = 0, 
    Flavor_4e    = 1,
    Flavor_4mu   = 2,
    Flavor_2e2mu = 3    
  };

  //---------------------------------
  // Function
  //---------------------------------
  static TString ProcessName(int temp){ 
    if(temp==TVar::ZZ_2e2m   ) 
      return TString("ZZ_2e2m");
    if(temp==TVar::ZZ_4e   ) 
      return TString("ZZ_4e");
    if(temp==TVar::GGZZ_4l   ) 
      return TString("GGZZ_4l");
    else if(temp==TVar::HZZ_4l   ) 
      return TString("HZZ_4l");
    else if(temp==TVar::PSHZZ_4l   ) 
      return TString("PSHZZ_4l");
    else if(temp==TVar::HDHZZ_4l   ) 
      return TString("HDHZZ_4l");
    else if(temp==TVar::TZZ_4l   ) 
      return TString("TZZ_2mplus_4l");
    else if(temp==TVar::QQB_TZZ_4l   ) 
      return TString("QQB_TZZ_2mplus_4l");
    else if(temp==TVar::TZZ_DECAY_4l   ) 
      return TString("TZZ_decay_2mplus_4l");
    else if(temp==TVar::VZZ_4l   ) 
      return TString("VZZ_4l");
    else if(temp==TVar::AVZZ_4l   ) 
      return TString("AVZZ_4l");
    else if(temp==TVar::VZZ_DECAY_4l   ) 
      return TString("VZZ_decay_4l");
    else if(temp==TVar::AVZZ_DECAY_4l   ) 
      return TString("AVZZ_decay_4l");
    else if(temp==TVar::PTZZ_2hminus_4l   ) 
      return TString("PTZZ_2hminus_4l");
    else if(temp==TVar::TZZ_2hplus_4l   ) 
      return TString("TZZ_2hplus_4l");
    else if(temp==TVar::TZZ_2bplus_4l   ) 
      return TString("TZZ_2bplus_4l");
    else if(temp==TVar::SummedBackgrounds   ) 
      return TString("SummedBackgrounds");
    else if(temp==TVar::HZZ_4l_MIXCP   ) 
      return TString("HZZ_4l_MIXCP");
    else if(temp==TVar::HJJNONVBF ) 
      return TString("HJJNONVBF");
    else if(temp==TVar::HJJVBF ) 
      return TString("HJJVBF");
    else if(temp==TVar::GGZZTOT_4l   ) 
      return TString("GGZZ_total_4l");
   else if(temp==TVar::GGZZINT_4l   ) 
      return TString("GGZZ_interference_4l");
   else 
      return TString("UnKnown");
  };
  ClassDef(TVar,0)
};

struct branch_particle {
  int   PdgCode   ;
  int   Charge    ;
  double Px       ;
  double Py       ;
  double Pz       ;
  double E        ;
  double Eta      ;
  double Phi      ;

};
static const TString branch_format_particle =
 "PdgCode/I:"
 "Charge/I:"
 "Px/D:"
 "Py/D:"
 "Pz/D:"
 "E/D:"
 "Eta/D:"
 "Phi/D";

// in development
struct hzz4l_event_type{
  int PdgCode[4];
  TLorentzVector p[4];
  double Xsec   [10];
  double XsecErr[10];  
};
// in development
struct xjj_event_type{
  TLorentzVector p[3]; // p[0] and p[1] for the two jets and p[2] for the resonance X
  double Xsec   [10];
  double XsecErr[10];  
};

struct mcfm_event_type{
  int PdgCode[6];
  TLorentzVector p[6];
  double pswt;
};
struct event_type{
  TLorentzVector p1,p2,ep,em,nu,nb;
  double PSWeight;
};
struct anomcoup{
	   double delg1_z, delg1_g, lambda_g, lambda_z, delk_g, delk_z_,tevscale;
};

#endif
