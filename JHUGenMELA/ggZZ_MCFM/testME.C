// 
// Author Yanyan Gao (Yanyan.Gao@cern.ch)
// 

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TProfile.h"
#include <iostream>
#include "Math/LorentzVector.h"
#include "TLorentzRotation.h"
#include "Math/VectorUtil.h"
// ME related
#include "TVar.hh"
#include "TEvtProb.hh"
#include "math.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector; 

float ERRORthreshold=1.0;
using namespace std;

vector<TLorentzVector> Calculate4Momentum(float Mx,float M1,float M2,float theta,float theta1,float theta2,float Phi1,float Phi)
{
  float phi1,phi2;
  phi1=TMath::Pi()-Phi1;
  phi2=Phi1+Phi;


  float gamma1,gamma2,beta1,beta2;

  gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
  gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
  beta1=sqrt(1-1/(gamma1*gamma1));
  beta2=sqrt(1-1/(gamma2*gamma2));


  //gluon 4 vectors
  TLorentzVector p1CM(0,0,Mx/2,Mx/2);
  TLorentzVector p2CM(0,0,-Mx/2,Mx/2);

  //vector boson 4 vectors
  TLorentzVector kZ1(gamma1*M1*sin(theta)*beta1,0, gamma1*M1*cos(theta)*beta1,gamma1*M1*1);   
  TLorentzVector kZ2(-gamma2*M2*sin(theta)*beta2,0, -gamma2*M2*cos(theta)*beta2,gamma2*M2*1);

  //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.

  TLorentzRotation Z1ToZ,Z2ToZ;

  Z1ToZ.Boost(0,0,beta1);
  Z2ToZ.Boost(0,0,beta2);
  Z1ToZ.RotateY(theta);
  Z2ToZ.RotateY(TMath::Pi()+theta);


  //fermons 4 vectors in vector boson rest frame

  TLorentzVector p3Z1((M1/2)*sin(theta1)*cos(phi1),(M1/2)*sin(theta1)*sin(phi1),(M1/2)*cos(theta1),(M1/2)*1);

  TLorentzVector p4Z1(-(M1/2)*sin(theta1)*cos(phi1),-(M1/2)*sin(theta1)*sin(phi1),-(M1/2)*cos(theta1),(M1/2)*1);

  TLorentzVector p5Z2((M2/2)*sin(theta2)*cos(phi2),(M2/2)*sin(theta2)*sin(phi2),(M2/2)*cos(theta2),(M2/2)*1);

  TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2),-(M2/2)*sin(theta2)*sin(phi2),-(M2/2)*cos(theta2),(M2/2)*1);


  // fermions 4 vectors in CM frame

  TLorentzVector p3CM,p4CM,p5CM,p6CM;

  p3CM=Z1ToZ*p3Z1;
  p4CM=Z1ToZ*p4Z1;
  p5CM=Z2ToZ*p5Z2;
  p6CM=Z2ToZ*p6Z2;

  vector<TLorentzVector> p;

  p.push_back(p3CM);
  p.push_back(p4CM);
  p.push_back(p5CM);
  p.push_back(p6CM);

  return p;
}



//###################
//# main function
//###################
void testME() {

  TVar::VerbosityLevel verbosity = TVar::INFO;

  double dXsec_ZZ_MCFM = 0.;
  double dXsec_GGZZ_MCFM = 0.;
  double dXsec_GGZZTOT_MCFM = 0.;
  double dXsec_GGZZINT_MCFM = 0.;
  double dXsec_HZZ_MCFM = 0.;
  double dXsec_HZZ_JHU = 0.;
  double dXsec_PSHZZ_JHU = 0.;
  double dXsec_HDHZZ_JHU = 0.;
  double dXsec_TZZ_JHU = 0.;
  double dXsec_VZZ_JHU = 0.;
  double dXsec_AVZZ_JHU = 0.;
  double dXsec_QQB_TZZ_JHU = 0.;
  double dXsec_TZZ_DECAY_JHU = 0.;
  double dXsec_VZZ_DECAY_JHU = 0.;
  double dXsec_AVZZ_DECAY_JHU = 0.;
  double dXsec_PTZZ_2hminus_JHU = 0.;
  double dXsec_TZZ_2hplus_JHU = 0.;
  double dXsec_TZZ_2bplus_JHU = 0.;
  double dXsec_HZZ_MIXCP_JHU = 0.;
  double dXsec_HJJ_JHU = 0.;
  double dXsec_HJJVBF_JHU = 0.;

  float mzz = 126.; 
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  int mflavor = 3; 

	double couplingvals[SIZE_HVV_FREENORM] = {0};
	double selfDHggcoupl[SIZE_HVV][2] = {{0}};
	double selfDHvvcoupl_VBF[SIZE_HVV_VBF][2] = {{0}};
	double selfDHwwcoupl_VBF[SIZE_HWW_VBF][2] = {{0}};
	double selfDHvvcoupl[SIZE_HVV][2] = {{0}};
	double selfDZqqcoupl[SIZE_ZQQ][2] = {{0}};
	double selfDZvvcoupl[SIZE_ZVV][2] = {{0}};
	double selfDGqqcoupl[SIZE_GQQ][2] = {{0}};
	double selfDGggcoupl[SIZE_GGG][2] = {{0}};
	double selfDGvvcoupl[SIZE_GVV][2] = {{0}};

  // Create the instance of TEvtProb to calculate the differential cross-section
  TEvtProb Xcal2(4000);
  hzz4l_event_type hzz4l_event;
  Xcal2.ResetMCFM_EWKParameters(1.16639E-05,7.81751E-03,79.956049884402844,91.1876);

  // set four momenta
  vector<TLorentzVector> p;
  p=Calculate4Momentum(mzz,m1,m2,acos(hs),acos(h1),acos(h2),phi1,phi);

  TLorentzVector Z1_minus = p[0];
  TLorentzVector Z1_plus  = p[1];
  TLorentzVector Z2_minus = p[2];
  TLorentzVector Z2_plus  = p[3];

  hzz4l_event.p[0].SetXYZM(Z1_minus.Px(), Z1_minus.Py(), Z1_minus.Pz(), 0.);
  hzz4l_event.p[1].SetXYZM(Z1_plus.Px(), Z1_plus.Py(), Z1_plus.Pz(), 0.);
  hzz4l_event.p[2].SetXYZM(Z2_minus.Px(), Z2_minus.Py(), Z2_minus.Pz(), 0.);
  hzz4l_event.p[3].SetXYZM(Z2_plus.Px(), Z2_plus.Py(), Z2_plus.Pz(), 0.);


  // flavor 1 for 4e, 2 for 4m, 3 for 2e2mu  
  if ( mflavor == 1 ) {
    hzz4l_event.PdgCode[0] = 11;
    hzz4l_event.PdgCode[1] = -11;
    hzz4l_event.PdgCode[2] = 11;
    hzz4l_event.PdgCode[3] = -11;
  }
  if ( mflavor == 2 ) {
    hzz4l_event.PdgCode[0] = 13;
    hzz4l_event.PdgCode[1] = -13;
    hzz4l_event.PdgCode[2] = 13;
    hzz4l_event.PdgCode[3] = -13;
  }
  if ( mflavor == 3 ) {
    hzz4l_event.PdgCode[0] = 11;
    hzz4l_event.PdgCode[1] = -11;
    hzz4l_event.PdgCode[2] = 13;
    hzz4l_event.PdgCode[3] = -13;
  }


  float z1mass = (hzz4l_event.p[0]+hzz4l_event.p[1]).M();
  float z2mass = (hzz4l_event.p[2]+hzz4l_event.p[3]).M();
  float zzmass = (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).M();

  if (verbosity >= TVar::INFO) {
    std::cout << "Input: ==================================================" <<endl;
    printf("lep1 (Px, Py, Pz, E) = (%4.4f, %4.4f, %4.4f, %4.4f)\n",  p[0].Px(), p[0].Py(), p[0].Pz(), p[0].E());
    printf("lep2 (Px, Py, Pz, E) = (%4.4f, %4.4f, %4.4f, %4.4f)\n",  p[1].Px(), p[1].Py(), p[1].Pz(), p[1].E()); 
    printf("lep3 (Px, Py, Pz, E) = (%4.4f, %4.4f, %4.4f, %4.4f)\n",  p[2].Px(), p[2].Py(), p[2].Pz(), p[2].E());
    printf("lep4 (Px, Py, Pz, E) = (%4.4f, %4.4f, %4.4f, %4.4f)\n",  p[3].Px(), p[3].Py(), p[3].Pz(), p[3].E()); 
    std::cout << "ZZ system (pX, pY, pZ, E, mass) = ( " 
      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Px() << ", "
      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Py() << ", "
      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Pz() << ", "
      << (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Energy()  << ", "
      << zzmass << ")\n";
    std::cout << "Z1 mass = " << z1mass << "\tz2mass = " << z2mass << "\n";
    std::cout << "=========================================================\n";
  } 
  // finish event information

  // ==== Begin the differential cross-section calculation
  Xcal2.SetHiggsMass(zzmass,0.004);
  // calculate the ZZ using MCFM
  Xcal2.SetMatrixElement(TVar::MCFM);
  if ( mflavor < 3  )
    dXsec_ZZ_MCFM = Xcal2.XsecCalc(TVar::bkgZZ, TVar::ZZQQB, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);
  else 
    dXsec_ZZ_MCFM = Xcal2.XsecCalc(TVar::bkgZZ, TVar::ZZQQB, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);
  dXsec_GGZZ_MCFM = Xcal2.XsecCalc(TVar::bkgZZ, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);
  dXsec_HZZ_MCFM = Xcal2.XsecCalc(TVar::HSMHiggs, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);
  dXsec_GGZZTOT_MCFM = Xcal2.XsecCalc(TVar::bkgZZ_SMHiggs, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);
  dXsec_GGZZINT_MCFM = dXsec_GGZZTOT_MCFM - dXsec_HZZ_MCFM - dXsec_GGZZ_MCFM;

  if ( verbosity >= TVar::INFO ) {
    std::cout << "total HZZ+ggZZ (132 + 128 + 129) = " << dXsec_GGZZ_MCFM + dXsec_HZZ_MCFM + dXsec_GGZZINT_MCFM 
      << "\t, difference to process 131 " << dXsec_GGZZ_MCFM + dXsec_HZZ_MCFM + dXsec_GGZZINT_MCFM - dXsec_GGZZTOT_MCFM      
      << std::endl; 
    std::cout << "interference from (131 - 132 - 128) is " << dXsec_GGZZTOT_MCFM - dXsec_HZZ_MCFM - dXsec_GGZZ_MCFM << std::endl;
  }

  // calculate X->ZZ using JHUGen
  Xcal2.SetMatrixElement(TVar::JHUGen);

  // 0+ 
  dXsec_HZZ_JHU =  Xcal2.XsecCalc(TVar::HSMHiggs, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);

  // 0-
  dXsec_PSHZZ_JHU = Xcal2.XsecCalc(TVar::H0minus, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);

  // 0h+
  dXsec_HDHZZ_JHU = Xcal2.XsecCalc(TVar::H0hplus, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);

  // 0 mix cp
//  dXsec_HZZ_MIXCP_JHU = Xcal2.XsecCalc(TVar::D_g1g4, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);

  // 1-
  dXsec_VZZ_JHU = Xcal2.XsecCalc(TVar::H1minus, TVar::ZZQQB, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);    

  // 1+
  dXsec_AVZZ_JHU = Xcal2.XsecCalc(TVar::H1plus, TVar::ZZQQB, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);   

  // 1- decay only
  dXsec_VZZ_DECAY_JHU = Xcal2.XsecCalc(TVar::H1minus, TVar::ZZINDEPENDENT, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl); 

  // 1+ decay only 
  dXsec_AVZZ_DECAY_JHU = Xcal2.XsecCalc(TVar::H1plus, TVar::ZZINDEPENDENT, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl); 

  // 2m+ (gg production) 
  dXsec_TZZ_JHU = Xcal2.XsecCalc(TVar::H2_g1g5, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);
  
  // 2m+ (qqbar production)
  dXsec_QQB_TZZ_JHU = Xcal2.XsecCalc(TVar::H2_g1g5,TVar::ZZQQB, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);
  
  // 2m+ decay only
  dXsec_TZZ_DECAY_JHU = Xcal2.XsecCalc(TVar::H2_g1g5, TVar::ZZINDEPENDENT, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);

  // 2h-
  dXsec_PTZZ_2hminus_JHU = Xcal2.XsecCalc(TVar::H2_g8, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);

  // 2h+
  dXsec_TZZ_2hplus_JHU = Xcal2.XsecCalc(TVar::H2_g4, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);

  // 2b+
  dXsec_TZZ_2bplus_JHU = Xcal2.XsecCalc(TVar::H2_g5, TVar::ZZGG, hzz4l_event,verbosity,couplingvals,selfDHvvcoupl,selfDZqqcoupl,selfDZvvcoupl,selfDGqqcoupl,selfDGggcoupl,selfDGvvcoupl);

  // H+jj
  // calculate the p4 of the H + 2jets, boosted to have 0 pT
  TLorentzVector p4[3];
  TLorentzVector tot;
  //p4[0] for j1,  p4[1] for j2,  p4[2] for H
  p4[0].SetPxPyPzE( 335.75180561872444, 68.440486224675223, -4.90610220658908303, 342.69147198283938);   
  p4[1].SetPxPyPzE( 126.89634103814822, 60.252128195920029, -640.04704411152327, 655.28102291833976);      
  p4[2].SetPxPyPzE(-462.64814665687268, -128.69261442059525, -281.52865302540238, 570.73591010706917);      
  if ( verbosity >= TVar::INFO ) {
    std::cout << "========================================\n";
    std::cout << "Printing H+2j information " << "\n";
    std::cout << "========================================\n";
    std::cout << Form("Jet 1 (px,py,pz,m) = (%.5f, %.5f, %.5f, %.5f)\n", p4[0].Px(), p4[0].Py(), p4[0].Pz(), p4[0].M()); 
    std::cout << Form("Jet 2 (px,py,pz,m) = (%.5f, %.5f, %.5f, %.f)\n", p4[1].Px(), p4[1].Py(), p4[1].Pz(), p4[1].M()); 
    std::cout << Form("ZZ system (px,py,pz,m) = (%.5f, %.5f, %.5f, %.5f)\n", p4[2].Px(), p4[2].Py(), p4[2].Pz(), p4[2].M());
    std::cout << "========================================\n";
  }

  dXsec_HJJ_JHU = Xcal2.XsecCalcXJJ(TVar::HSMHiggs,TVar::JJGG, p4, verbosity,selfDHggcoupl,selfDHvvcoupl_VBF,selfDHwwcoupl_VBF);
  dXsec_HJJVBF_JHU = Xcal2.XsecCalcXJJ(TVar::HSMHiggs,TVar::JJVBF, p4, verbosity,selfDHggcoupl,selfDHvvcoupl_VBF,selfDHwwcoupl_VBF);

  if ( verbosity >= TVar::INFO ) { 
    FILE *output = fopen("output.txt", "w");
    fprintf(output, "==========Matrix Element outputs==============\n");
    fprintf(output, "==============================================\n");
    fprintf(output, "******From MCFM: ****\n");
    fprintf(output, "==============================================\n");
    fprintf(output, "gg->H->ZZ: %7.7e\n",  dXsec_HZZ_MCFM); 
    fprintf(output, "qqbar->ZZ: %7.7e\n",  dXsec_ZZ_MCFM); 
    fprintf(output, "gg->ZZ %7.7e\n", dXsec_GGZZ_MCFM);
    fprintf(output, "gg->ZZ(including HZZ): %7.7e\n", dXsec_GGZZTOT_MCFM);
    fprintf(output, "gg->ZZ and gg->H->ZZ interference:%7.7e\n", dXsec_GGZZINT_MCFM);
    fprintf(output, "==============================================\n");
    fprintf(output, "******From JHUGenMELA: ****\n");
    fprintf(output, "==============================================\n");
    fprintf(output, "gg->H->ZZ:%7.7e\n",  dXsec_HZZ_JHU);
    fprintf(output, "gg->X(0-)->ZZ:%7.7e\n",  dXsec_PSHZZ_JHU);
    fprintf(output, "gg->X(0h+)->ZZ:%7.7e\n",  dXsec_HDHZZ_JHU);
    fprintf(output, "gg->X(0 mix cp)->ZZ:%7.7e\n",  dXsec_HZZ_MIXCP_JHU);
    fprintf(output, "gg->X(1-)->ZZ:%7.7e\n",  dXsec_VZZ_JHU);
    fprintf(output, "gg->X(1+)->ZZ:%7.7e\n",  dXsec_AVZZ_JHU);
    fprintf(output, "X(1-)->ZZ production independent:%7.7e\n",  dXsec_VZZ_DECAY_JHU);
    fprintf(output, "X(1+)->ZZ production independent:%7.7e\n",  dXsec_AVZZ_DECAY_JHU);
    fprintf(output, "gg->X(2m+)->ZZ:%7.7e\n",  dXsec_TZZ_JHU);
    fprintf(output, "qqbar->X(2m+)->ZZ:%7.7e\n",  dXsec_QQB_TZZ_JHU);
    fprintf(output, "gg->X(2m-)->ZZ:%7.7e\n",  dXsec_PTZZ_2hminus_JHU);
    fprintf(output, "gg->X(2h+)->ZZ:%7.7e\n",  dXsec_TZZ_2hplus_JHU);
    fprintf(output, "gg->X(2b+)->ZZ:%7.7e\n",  dXsec_TZZ_2bplus_JHU);
    fprintf(output, "X(2m+)->ZZ production independent:%7.7e\n",  dXsec_TZZ_DECAY_JHU);
    fprintf(output, "H+JJ (SBF): %7.7e\n", dXsec_HJJ_JHU);
    fprintf(output, "H+JJ (WBF): %7.7e\n", dXsec_HJJVBF_JHU);
    fprintf(output, "===============================================\n" );
    fclose(output);
  }
}
