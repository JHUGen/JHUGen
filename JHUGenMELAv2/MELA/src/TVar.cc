#include "TVar.hh"


TString TVar::ProcessName(TVar::Process temp){
  if (temp==TVar::HSMHiggs) return TString("HSMHiggs");
  else if (temp==TVar::H0minus) return TString("H0minus");
  else if (temp==TVar::H0hplus) return TString("H0hplus");
  else if (temp==TVar::H0_g1prime2) return TString("H0_g1prime2");
  else if (temp==TVar::H0_Zgs) return TString("H0_Zgs");
  else if (temp==TVar::H0_gsgs) return TString("H0_gsgs");
  else if (temp==TVar::H0_Zgs_PS) return TString("H0_Zgs_PS");
  else if (temp==TVar::H0_gsgs_PS) return TString("H0_gsgs_PS");
  else if (temp==TVar::H0_Zgsg1prime2) return TString("H0_Zgsg1prime2");

  else if (temp==TVar::D_g1g4) return TString("D_g1g4");
  else if (temp==TVar::D_g1g4_pi_2) return TString("D_g1g4_pi_2");
  else if (temp==TVar::D_g1g2) return TString("D_g1g2");
  else if (temp==TVar::D_g1g2_pi_2) return TString("D_g1g2_pi_2");
  else if (temp==TVar::D_g1g1prime2) return TString("D_g1g1prime2");
  else if (temp==TVar::D_zzzg) return TString("D_zzzg");
  else if (temp==TVar::D_zzgg) return TString("D_zzgg");
  else if (temp==TVar::D_zzzg_PS) return TString("D_zzzg_PS");
  else if (temp==TVar::D_zzgg_PS) return TString("D_zzgg_PS");
  else if (temp==TVar::D_zzzg_g1prime2) return TString("D_zzzg_g1prime2");
  else if (temp==TVar::D_zzzg_g1prime2_pi_2) return TString("D_zzzg_g1prime2_pi_2");

  else if (temp==TVar::H1minus) return TString("H1minus");
  else if (temp==TVar::H1plus) return TString("H1plus");

  else if (temp==TVar::H2_g1) return TString("H2_g1");
  else if (temp==TVar::H2_g2) return TString("H2_g2");
  else if (temp==TVar::H2_g3) return TString("H2_g3");
  else if (temp==TVar::H2_g4) return TString("H2_g4");
  else if (temp==TVar::H2_g5) return TString("H2_g5");
  else if (temp==TVar::H2_g1g5) return TString("H2_g1g5");
  else if (temp==TVar::H2_g6) return TString("H2_g6");
  else if (temp==TVar::H2_g7) return TString("H2_g7");
  else if (temp==TVar::H2_g8) return TString("H2_g8");
  else if (temp==TVar::H2_g9) return TString("H2_g9");
  else if (temp==TVar::H2_g10) return TString("H2_g10");

  else if (temp==TVar::bkgZGamma) return TString("bkgZGamma");
  else if (temp==TVar::bkgZJets) return TString("bkgZJets");
  else if (temp==TVar::bkgZZ) return TString("bkgZZ");
  else if (temp==TVar::bkgWW) return TString("bkgWW");
  else if (temp==TVar::bkgWWZZ) return TString("bkgWWZZ");
  else if (temp==TVar::bkgZZ_SMHiggs) return TString("bkgZZ_SMHiggs");
  else if (temp==TVar::bkgWW_SMHiggs) return TString("bkgWW_SMHiggs");
  else if (temp==TVar::bkgWWZZ_SMHiggs) return TString("bkgWWZZ_SMHiggs");
  else if (temp==TVar::HSMHiggs_WWZZ) return TString("HSMHiggs_WWZZ");

  else if (temp==TVar::D_gg10) return TString("D_gg10");

  else if (temp==TVar::SelfDefine_spin0) return TString("SelfDefine_spin0");
  else if (temp==TVar::SelfDefine_spin1) return TString("SelfDefine_spin1");
  else if (temp==TVar::SelfDefine_spin2) return TString("SelfDefine_spin2");

  else return TString("Unknown");
}

TString TVar::ProductionName(TVar::Production temp){
  if (temp==TVar::ZZGG) return TString("ZZGG");
  else if (temp==TVar::ZZQQB) return TString("ZZQQB");
  else if (temp==TVar::ZZQQB_STU) return TString("ZZQQB_STU");
  else if (temp==TVar::ZZINDEPENDENT) return TString("ZZINDEPENDENT");

  else if (temp==TVar::ttH) return TString("ttH");
  else if (temp==TVar::bbH) return TString("bbH");

  else if (temp==TVar::JQCD) return TString("JQCD");

  else if (temp==TVar::JJQCD) return TString("JJQCD");
  else if (temp==TVar::JJVBF) return TString("JJVBF");
  else if (temp==TVar::JJEW) return TString("JJEW");
  else if (temp==TVar::JJEWQCD) return TString("JJEWQCD");
  else if (temp==TVar::Had_ZH) return TString("Had_ZH");
  else if (temp==TVar::Had_WH) return TString("Had_WH");
  else if (temp==TVar::Lep_ZH) return TString("Lep_ZH");
  else if (temp==TVar::Lep_WH) return TString("Lep_WH");

  else if (temp==TVar::ZZQQB_S) return TString("ZZQQB_S");
  else if (temp==TVar::JJQCD_S) return TString("JJQCD_S");
  else if (temp==TVar::JJVBF_S) return TString("JJVBF_S");
  else if (temp==TVar::JJEW_S) return TString("JJEW_S");
  else if (temp==TVar::JJEWQCD_S) return TString("JJEWQCD_S");
  else if (temp==TVar::Had_ZH_S) return TString("Had_ZH_S");
  else if (temp==TVar::Had_WH_S) return TString("Had_WH_S");
  else if (temp==TVar::Lep_ZH_S) return TString("Lep_ZH_S");
  else if (temp==TVar::Lep_WH_S) return TString("Lep_WH_S");

  else if (temp==TVar::ZZQQB_TU) return TString("ZZQQB_TU");
  else if (temp==TVar::JJQCD_TU) return TString("JJQCD_TU");
  else if (temp==TVar::JJVBF_TU) return TString("JJVBF_TU");
  else if (temp==TVar::JJEW_TU) return TString("JJEW_TU");
  else if (temp==TVar::JJEWQCD_TU) return TString("JJEWQCD_TU");
  else if (temp==TVar::Had_ZH_TU) return TString("Had_ZH_TU");
  else if (temp==TVar::Had_WH_TU) return TString("Had_WH_TU");
  else if (temp==TVar::Lep_ZH_TU) return TString("Lep_ZH_TU");
  else if (temp==TVar::Lep_WH_TU) return TString("Lep_WH_TU");

  else if (temp==TVar::GammaH) return TString("GammaH");

  else return TString("Unknown");
}

TString TVar::MatrixElementName(TVar::MatrixElement temp){
  if (temp==TVar::MCFM) return TString("MCFM");
  else if (temp==TVar::JHUGen) return TString("JHUGen");
  else if (temp==TVar::ANALYTICAL) return TString("ANALYTICAL");

  else return TString("Unknown");
}

