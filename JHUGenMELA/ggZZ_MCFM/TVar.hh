#ifndef EvtProb_VAR
#define EvtProb_VAR

#include <cstring>
#include <string>
#include <TLorentzVector.h>
#include "TH2F.h"
#include "TH1F.h"

//#define EBEAM 4000.00
#define fbGeV2 0.389379E12
#define smallnumber 1e-15
#define sixteen_2Pi_to_8 3.88650230418250561e+07
#define eight_2Pi_to_5 7.83410393050320417e+04
#define four_2Pi_to_2 39.478417604357432
//---------------------------------
// Coupling array sizes
//---------------------------------
#define SIZE_HVV 39
#define SIZE_HVV_VBF 32
#define SIZE_HWW_VBF 32
#define SIZE_HGG 3
#define SIZE_ZQQ 2
#define SIZE_ZVV 2
#define SIZE_GQQ 2
#define SIZE_GGG 5
#define SIZE_GVV 10
#define SIZE_HVV_FREENORM 2


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
    ZZGG = 0,
    ZZQQB = 1,
    ZZQQB_STU = 2,
    ZZQQB_S = 3,
    ZZQQB_TU = 4,
    ZZINDEPENDENT= 5,
    JJGG = 6, // SBF
    JJVBF = 7, // WBF
    JH = 8, // H + 1 jet
    ZH = 9, // ZH
    WH = 10 // W(+/-)H
//
  };
  enum LeptonInterference{
    DefaultLeptonInterf = 0,
    InterfOn = 1,
    InterfOff=2
  };

  enum Process{

    HSMHiggs          = 0,    //0+, replacing HZZ_4l, when production is ZZGG, replacing HJJNONVBF/HJJVBF/HJJVH, when production is JJGG/JJVBF/JJVH. Call this for MCFM |H|**2.
    H0minus           = 1,    //0-, replacing PSHZZ_4l, when production is ZZGG, replacing PSHJJNONVBF/PSHJJVBF/PSHJJVH when production is JJGG/JJVBF/JJVH
    H0hplus           = 2,    //0h+, replacing HDHZZ_4l

    H1minus           = 3,    //1-, replacing VZZ_4l
    H1plus            = 4,    //1+, replacing AVZZ_4l

    H2_g8             = 5,    //2h-, replacing PTZZ_2hminus_4l	
    H2_g4             = 6,    //2h+, replacing TZZ_2hplus_4l
    H2_g5             = 7,    //2b+, replacing TZZ_2bplus_4l 
    H2_g1g5           = 8,    //2m+, replacing TZZ_4l 
		H2_g2					= 9, // 2h2+
		H2_g3					= 10, // 2h3+
		H2_g6					= 11, // 2h6+
		H2_g7					= 12, // 2h7+
		H2_g9					= 13, // 2h9-
		H2_g10				= 14, // 2h10-


    bkgZZ              = 15,    //qq->ZZ, replacing ZZ_2e2m & ZZ_4e, when production is ZZQQB, replacing GGZZ_4l when production is ggZZ, replacing SummedBackgrounds for superMela calculation 
	bkgZZ_SMHiggs      =16,    //ggZZ+SMHiggs, ggZZ always calculated by MCFM, ME stands for SMHiggs ME, JHUGen: MCFM ggZZ + JHUGen SMHiggs, MCFM: MCFM (ggZZ+ SMHiggs) 

    H0_g1prime2       = 17,   //g1=0, g1prime2=-12046.01, replacing H_g1prime2	

    /***For interaction terms **/
    D_g1g4            = 18,   //D_CP
    D_g1g4_pi_2       = 19,   //D_CP_T
    D_g1g2            = 20,   //D_int
    D_g1g2_pi_2       = 21,   //D_int_T
    D_g1g1prime2      = 22,   //D_int_lambda1	

    /***** Self Defined******/
    SelfDefine_spin0  = 23,
    SelfDefine_spin1  = 24,
    SelfDefine_spin2  = 25,
		/**** For width ***/
		D_gg10						= 26,
		
		H0_Zgs 						= 27,
		H0_gsgs 					= 28,
		D_zzzg						= 29,
		D_zzgg						= 30,

		H0_Zgs_PS 					= 31,
		H0_gsgs_PS 					= 32,
		D_zzzg_PS						= 33,
		D_zzgg_PS						= 34,
		
		H0_Zgsg1prime2						= 35,
		D_zzzg_g1prime2						= 36,
		D_zzzg_g1prime2_pi_2						= 37,

	/*** Are these ones still used? ***/
    //  QQB_TZZ_4l = 8,
    //  TZZ_DECAY_4l = 10,
    //  VZZ_DECAY_4l = 11,
    //  AVZZ_DECAY_4l = 12,
    //  HZZ_4l_MIXCP = 17,

    Null
  };
  enum LeptonFlavor{
    Flavor_Dummy = 0, // legacy code runs on 1/2/3
    Flavor_4e    = 1,
    Flavor_4mu   = 2,
    Flavor_2e2mu = 3    
  };
  enum SuperMelaSyst{
    SMSyst_None      = 0, // nominal value
    SMSyst_ScaleUp   = 1, //Scale Uncertaintie
    SMSyst_ScaleDown = 2,
    SMSyst_ResUp     = 3, // Resolution Uncertainty
    SMSyst_ResDown   = 4
  };

  //---------------------------------
  // Function
  //---------------------------------
  static TString ProcessName(int temp){
    if     (temp==TVar::HSMHiggs          ) return TString ("HSMHiggs");
    else if(temp==TVar::H0minus           ) return TString ("H0minus");           
    else if(temp==TVar::H0hplus           ) return TString ("H0hplus");          
                                                                               
    else if(temp==TVar::H1minus           ) return TString ("H1minus");           
    else if(temp==TVar::H1plus            ) return TString ("H1plus");            
                                                                               
    else if(temp==TVar::H2_g8             ) return TString ("H2_g8");             
    else if(temp==TVar::H2_g4             ) return TString ("H2_g4");             
    else if(temp==TVar::H2_g5             ) return TString ("H2_g5");             
    else if(temp==TVar::H2_g1g5           ) return TString ("H2_g1g5");           
    else if(temp == TVar::H2_g2 					) return TString ("H2_g2");
    else if(temp == TVar::H2_g3 					) return TString ("H2_g3");
    else if(temp == TVar::H2_g6 					) return TString ("H2_g6");
    else if(temp == TVar::H2_g7 					) return TString ("H2_g7");
    else if(temp == TVar::H2_g9 					) return TString ("H2_g9");
    else if(temp == TVar::H2_g10    		  ) return TString ("H2_g10");
                                                                               
                                                                               
    else if(temp==TVar::bkgZZ             ) return TString ("bkgZZ");             
	else if(temp==TVar::bkgZZ_SMHiggs     ) return TString ("bkgZZ_SMHiggs");     
                                                                               
    else if(temp==TVar::H0_g1prime2       ) return TString ("H0_g1prime2");       
                                                                               
    else if(temp==TVar::D_g1g4            ) return TString ("D_g1g4");            
    else if(temp==TVar::D_g1g4_pi_2       ) return TString ("D_g1g4_pi_2");       
    else if(temp==TVar::D_g1g2            ) return TString ("D_g1g2");            
    else if(temp==TVar::D_g1g2_pi_2       ) return TString ("D_g1g2_pi_2");       
    else if(temp==TVar::D_g1g1prime2      ) return TString ("D_g1g1prime2");      
                                                                               
    else if(temp==TVar::SelfDefine_spin0  ) return TString ("SelfDefine_spin0");  
    else if(temp==TVar::SelfDefine_spin1  ) return TString ("SelfDefine_spin1");  
    else if(temp==TVar::SelfDefine_spin2  ) return TString ("SelfDefine_spin2");  

	else if(temp==TVar::D_gg10  ) return TString ("D_gg10");  

	else if(temp==TVar::H0_Zgs  ) return TString ("H0_Zgs");  
	else if(temp==TVar::H0_gsgs  ) return TString ("H0_gsgs");  
	else if(temp==TVar::D_zzzg  ) return TString ("D_zzzg");  
	else if(temp==TVar::D_zzgg  ) return TString ("D_zzgg");  

	else if(temp==TVar::H0_Zgs_PS  ) return TString ("H0_Zgs_PS");  
	else if(temp==TVar::H0_gsgs_PS  ) return TString ("H0_gsgs_PS");  
	else if(temp==TVar::D_zzzg_PS  ) return TString ("D_zzzg_PS");  
	else if(temp==TVar::D_zzgg_PS  ) return TString ("D_zzgg_PS");  

	else if(temp==TVar::H0_Zgsg1prime2  ) return TString ("H0_Zgsg1prime2");  
	else if(temp==TVar::D_zzzg_g1prime2  ) return TString ("D_zzzg_g1prime2");  

    else return TString ("UnKnown");
  };

  inline virtual ~TVar(){};
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
struct vh_event_type{ // ME is 2 -> 3
  int PdgCode[3];
  TLorentzVector p[3];
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
struct EffHist{
  TH2F* els_eff_mc;
  TH2F* mus_eff_mc;
};


#endif
