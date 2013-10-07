
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TUtil.hh"
#include "TCanvas.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "TProfile.h"

using namespace std;

void SetEwkCoupligParameters(){
  
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.aemmz_inp=7.81751E-03;
  ewinput_.wmass_inp=79.956049884402844;
  ewinput_.zmass_inp=91.1876;
  ewinput_.xw_inp=0.23116864;

}


void My_choose(TVar::Process process){

//ZZ_4l
if(process==TVar::ZZ_2e2m || process == TVar::GGZZ_4l  ){ 
 
    //81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
    //86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'
    //    nproc_.nproc=81;  
    //    chooser_();
  
  // these settings are identical to use the chooser_() function
    npart_.npart=4;
    nqcdjets_.nqcdjets=0;

    vsymfact_.vsymfact=1.0;                                                                                                               
    interference_.interference=false;

    nwz_.nwz=0;
    bveg1_mcfm_.ndim=10;
    masses_mcfm_.mb=0;
    breit_.n2=1;
    breit_.n3=1;

    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.0;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.0;
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;


 }  else if ( process == TVar::ZZ_4e) {

    // 90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))' 'L'
    // nproc_.nproc=90;
    //  chooser_();

    // these settings are  from 
    // ProdHep/chooser.f
    npart_.npart=4;
    nqcdjets_.nqcdjets=0;

    vsymfact_.vsymfact=0.25;                                                                                                               
    interference_.interference=true;

    nwz_.nwz=0;
    bveg1_mcfm_.ndim=10;
    masses_mcfm_.mb=0;
    breit_.n2=1;
    breit_.n3=1;

    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.0;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.0;
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;

    
 }  else if ( process == TVar::HZZ_4l) {

    // 114 '  f(p1)+f(p2) --> H(--> Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6))' 'N'
    // nproc_.nproc=114;
    // chooser_();
     
     npart_.npart=4;
     nqcdjets_.nqcdjets=0;

     bveg1_mcfm_.ndim=10;
     masses_mcfm_.mb=4.75; // necessary for the ggZZ and HZZ interference 

     breit_.n2=1;
     breit_.n3=1;

     breit_.mass2 =masses_mcfm_.zmass;
     breit_.width2=masses_mcfm_.zwidth;
     breit_.mass3 =masses_mcfm_.zmass;
     breit_.width3=masses_mcfm_.zwidth;
     
     zcouple_.l1=zcouple_.le;
     zcouple_.r1=zcouple_.re;
     
     zcouple_.l2=zcouple_.le;
     zcouple_.r2=zcouple_.re;

 } 
 else{
     std::cerr <<"[My_choose]: Can't identify Process: " << process <<endl;
 } 
}

bool My_masscuts(double s[][12],TVar::Process process){

 double minZmassSqr=10*10;

 if(process==TVar::ZZ_2e2m || process==TVar::ZZ_4e || process==TVar::GGZZ_4l){
   if(s[2][3]< minZmassSqr) return true;
   if(s[4][5]< minZmassSqr) return true;
 }
 return false;	 

}


bool My_smalls(double s[][12],int npart){

// Reject event if any s(i,j) is too small
// cutoff is defined in technical.Dat
	
      if ( 
       npart == 3 &&
       (
        (-s[5-1][1-1]< cutoff_.cutoff)  //gamma p1
     || (-s[5-1][2-1]< cutoff_.cutoff)  //gamma p2
     || (-s[4-1][1-1]< cutoff_.cutoff)  //e+    p1
     || (-s[4-1][2-1]< cutoff_.cutoff)  //e-    p2
     || (-s[3-1][1-1]< cutoff_.cutoff)  //nu    p1
     || (-s[3-1][2-1]< cutoff_.cutoff)  //nu    p2
     || (+s[5-1][4-1]< cutoff_.cutoff)  //gamma e+
     || (+s[5-1][3-1]< cutoff_.cutoff)  //gamma nu
     || (+s[4-1][3-1]< cutoff_.cutoff)  //e+    nu
	)	 
      ) 
        return true;
     
     else if (
       npart == 4 &&     
      (
        (-s[5-1][1-1]< cutoff_.cutoff)  //e-    p1
     || (-s[5-1][2-1]< cutoff_.cutoff)  //e-    p2
     || (-s[6-1][1-1]< cutoff_.cutoff)  //nb    p1
     || (-s[6-1][2-1]< cutoff_.cutoff)  //nb    p2
     || (+s[6-1][5-1]< cutoff_.cutoff)  //e-    nb
       )

     )
       
      return true;
     
     return false;
}




//Make sure
// 1. tot Energy Sum < 2EBEAM
// 2. PartonEnergy Fraction minimum<x0,x1<1
// 3. number of final state particle is defined
//
double SumMatrixElementPDF(TVar::Process process, mcfm_event_type* mcfm_event,double flavor_msq[nmsq][nmsq],double* flux){

  int NPart=npart_.npart+2;
  double p4[4][12];
  double msq[nmsq][nmsq];

  //Parton Density Function is always evalualted at pT=0 frame
  //Make sure parton Level Energy fraction is [0,1]
  //phase space function already makes sure the parton energy fraction between [min,1]
  //  x0 EBeam =>   <= -x1 EBeam
  
  double sysPz=mcfm_event->p[0].Pz()    +mcfm_event->p[1].Pz();
  double sysE =mcfm_event->p[0].Energy()+mcfm_event->p[1].Energy();
  
  //Ignore the Pt doesn't make significant effect
  //double sysPt_sqr=sysPx*sysPx+sysPy*sysPy;
  //if(sysPt_sqr>=1.0E-10)  sysE=TMath::Sqrt(sysE*sysE-sysPt_sqr);
  
  double xx[2]={(sysE+sysPz)/EBEAM/2,(sysE-sysPz)/EBEAM/2};
  if(xx[0] > 1.0 || xx[0]<=xmin_.xmin) return 0.0;
  if(xx[1] > 1.0 || xx[1]<=xmin_.xmin) return 0.0;
  
  //Convert TLorentzVector into 4x12 Matrix
  //reverse sign of incident partons back
  for(int ipar=0;ipar<2;ipar++){    
    if(mcfm_event->p[ipar].Energy()>0){
      p4[0][ipar] = -mcfm_event->p[ipar].Px();
      p4[1][ipar] = -mcfm_event->p[ipar].Py();
      p4[2][ipar] = -mcfm_event->p[ipar].Pz();
      p4[3][ipar] = -mcfm_event->p[ipar].Energy();
    }
  }
  //initialize decayed particles
  for(int ipar=2;ipar<NPart;ipar++){
    
    p4[0][ipar] = mcfm_event->p[ipar].Px();
    p4[1][ipar] = mcfm_event->p[ipar].Py();
    p4[2][ipar] = mcfm_event->p[ipar].Pz();
    p4[3][ipar] = mcfm_event->p[ipar].Energy();
    
  }
  
  //calculate invariant masses between partons/final state particles
  double s[12][12];
  for(int jdx=0;jdx< NPart ;jdx++){
    s[jdx][jdx]=0;
    for(int kdx=jdx+1;kdx<NPart;kdx++){
      s[jdx][kdx]=2*(p4[3][jdx]*p4[3][kdx]-p4[2][jdx]*p4[2][kdx]-p4[1][jdx]*p4[1][kdx]-p4[0][jdx]*p4[0][kdx]);
      s[kdx][jdx]=s[jdx][kdx];
    }
  }
  
  
  //remove events has small invariant mass
  // if(My_masscuts(s,process)) return 0.0;
  if(My_smalls(s,npart_.npart)) return 0.0;

  double msqjk=0;
  double msqgg=0;
    
  if( process==TVar::ZZ_2e2m || process==TVar::ZZ_4e )      qqb_zz_(p4[0],msq[0]);
  // MCFM-6.6 default gg->HZZ
  // if( process==TVar::HZZ_4l)     qqb_hzz_(p4[0],msq[0]);
  // MCFM-6.6 gg->HZZ + gg->ZZ including inteference
  if( process==TVar::HZZ_4l)  {
    double hcoupl[2] = {1.0, 0.0};
    gg_zz_int_freenorm_(p4[0], hcoupl, msq[0]);
  }
  if( process==TVar::GGZZ_4l)    gg_zz_  (p4[0],&msqgg);                     
  
  // by default assume only gg productions 
  // FOTRAN convention -5    -4   -3  -2    -1  0 1 2 3 4 5 
  //     parton flavor bbar cbar sbar ubar dbar g d u s c b
  // C++ convention     0     1   2    3    4   5 6 7 8 9 10

  msqjk=msq[5][5];
  if( process==TVar::ZZ_2e2m || process == TVar::ZZ_4e) msqjk=msq[3][7]+msq[7][3];
  // special for the GGZZ 
  if( process==TVar::GGZZ_4l) msqjk=msqgg;      
  
  (*flux)=fbGeV2/(8*xx[0]*xx[1]*EBEAM*EBEAM);
  
  if(msqjk != msqjk || flux!=flux ){
    cout << "SumMatrixPDF: "<< TVar::ProcessName(process) << " msqjk="  << msqjk << " flux="<< flux <<endl;
    msqjk=0;
    flux=0;
  }
  return msqjk;
  
}

//
// Test code from Markus to calculate the HZZ cross-section
// 
double JHUGenMatEl(TVar::Process process, TVar::Production production, mcfm_event_type* mcfm_event, double MReso, double GaReso, 
		   double Hggcoupl[3][2], double Hvvcoupl[4][2], double Zqqcoupl[2][2], double Zvvcoupl[2][2],
		   double Gqqcoupl[2][2], double Gggcoupl[5][2], double Gvvcoupl[10][2])
{
  // input unit = GeV/100 such that 125GeV is 1.25 in the code
  // this needs to be applied for all the p4
  MReso = MReso / 100.0;
  GaReso = GaReso /100.0;
  double p4[6][4];
  double MatElSq=0;
  int MYIDUP[4];

  int NPart = 6; 
  // p(i,0:3) = (E(i),px(i),py(i),pz(i))
  // i=0,1: glu1,glu2 (outgoing convention)
  // i=2,3: correspond to MY_IDUP(1),MY_IDUP(0)
  // i=4,5: correspond to MY_IDUP(3),MY_IDUP(2)
  for(int ipar=0;ipar<2;ipar++){   
    if(mcfm_event->p[ipar].Energy()>0){
      p4[ipar][0] = -mcfm_event->p[ipar].Energy()/100.;
      p4[ipar][1] = -mcfm_event->p[ipar].Px()/100.;
      p4[ipar][2] = -mcfm_event->p[ipar].Py()/100.;
      p4[ipar][3] = -mcfm_event->p[ipar].Pz()/100.;
    }
  }
  //initialize decayed particles
  for(int ipar=2;ipar<NPart;ipar++){
    p4[ipar][0] = mcfm_event->p[ipar].Energy()/100.;
    p4[ipar][1] = mcfm_event->p[ipar].Px()/100.;
    p4[ipar][2] = mcfm_event->p[ipar].Py()/100.;
    p4[ipar][3] = mcfm_event->p[ipar].Pz()/100.;
  }
  
  // particle ID: +7=e+,  -7=e-,  +8=mu+,  -8=mu-

  if ( TMath::Abs(mcfm_event->PdgCode[2]) == TMath::Abs(mcfm_event->PdgCode[3]) && 
       TMath::Abs(mcfm_event->PdgCode[3]) == TMath::Abs(mcfm_event->PdgCode[4]) && 
       TMath::Abs(mcfm_event->PdgCode[4]) == TMath::Abs(mcfm_event->PdgCode[5]) ) {
    if ( TMath::Abs(mcfm_event->PdgCode[2]) == 11  ) {
      MYIDUP[0]=+7;
      MYIDUP[1]=-7;
      MYIDUP[2]=+7;
      MYIDUP[3]=-7;
    } 
    if ( TMath::Abs(mcfm_event->PdgCode[2]) == 13  ) {
      MYIDUP[0]=+8;
      MYIDUP[1]=-8;
      MYIDUP[2]=+8;
      MYIDUP[3]=-8;
    } 
  } else {
      MYIDUP[0]=+7;
      MYIDUP[1]=-7;
      MYIDUP[2]=+8;
      MYIDUP[3]=-8;
  }
  if ( process == TVar::HZZ_4l || process == TVar::PSHZZ_4l || process == TVar::HDHZZ_4l || process == TVar::HZZ_4l_MIXCP ) {
    __modhiggs_MOD_evalamp_gg_h_vv(p4, &MReso,  &GaReso, Hggcoupl, Hvvcoupl, MYIDUP, &MatElSq);
  }
  if ( production == TVar::GG ) {
    if ( process == TVar::TZZ_4l || process == TVar::PTZZ_2hminus_4l || process == TVar::TZZ_2hplus_4l || process == TVar::TZZ_2bplus_4l ) {
      __modgraviton_MOD_evalamp_gg_g_vv(p4, &MReso,  &GaReso, Gggcoupl, Gvvcoupl, MYIDUP, &MatElSq);
    }
  }
  if ( production == TVar::INDEPENDENT ) {

    // special treatment of the 4-vectors
    // From Markus: 
    // Note that the momentum no.2, p(1:4,2), is a dummy which is not used. Momentum no.1,
    // p(1:4,1) = (/   -1.25d0,    0.00d0,    0.00d0,    0.00d0   /),
    // is the resonance momentum in its rest frame, which is crossed into the final state, 
    // i.e. the physical momentum is  -p(1:4,1) with a mass of 125GeV.
    
    double P[6][4];
    P[0][0]=-(mcfm_event->p[0]+mcfm_event->p[1]).M()/100.;
    P[0][1]=0.0;
    P[0][2]=0.0;
    P[0][3]=0.0;
    
    P[1][0]=0.0;
    P[1][1]=0.0;
    P[1][2]=0.0;
    P[1][3]=0.0;
    //initialize decayed particles
    for(int ipar=2;ipar<NPart;ipar++){
      P[ipar][0] = mcfm_event->p[ipar].Energy()/100.;
      P[ipar][1] = mcfm_event->p[ipar].Px()/100.;
      P[ipar][2] = mcfm_event->p[ipar].Py()/100.;
      P[ipar][3] = mcfm_event->p[ipar].Pz()/100.;
    }
    
    if ( process == TVar::TZZ_4l ) 
      __modgraviton_MOD_evalamp_g_vv(P, &MReso,  &GaReso, Gvvcoupl, MYIDUP, &MatElSq);

    if ( process == TVar::VZZ_4l || process == TVar::AVZZ_4l )
      __modzprime_MOD_evalamp_zprime_vv(P, &MReso,  &GaReso, Gvvcoupl, MYIDUP, &MatElSq);
  } 
  
  if ( production == TVar::QQB ) {
    if ( process == TVar::TZZ_4l ) {
      // -- YY: note that even if it is called xggcouplings, we are only testing xqq!
      __modgraviton_MOD_evalamp_qqb_g_vv(p4, &MReso,  &GaReso, Gggcoupl, Gvvcoupl, MYIDUP, &MatElSq);
    }
    if ( process == TVar::VZZ_4l || process == TVar::AVZZ_4l ) {
      // -- YY: note that even if it is called xggcouplings, we are only testing xqq!
      __modzprime_MOD_evalamp_qqb_zprime_vv(p4, &MReso,  &GaReso, Zqqcoupl, Zvvcoupl, MYIDUP, &MatElSq);
    }
  }
  
  /*
  printf("\n ");
  std::cout << "resoance = " << MReso *100. << ", width = " << GaReso*100. << "\n";
  for ( int i=0; i<NPart;i++) {
    std::cout << "p["<<i<<"] (E, Px, Py, Pz) = (" << p4[i][0] << ", " << p4[i][1] << ", " << p4[i][2] << ", " << p4[i][3] << ")\n";
  }
  printf("Matr.el. squared: %20.17e \n ",MatElSq);
  */
  // 
  // This constant is needed to account for the different units used in 
  // JHUGen compared to the MCFM
  // 
  double constant = 1.45/pow(10, 8);
  return MatElSq*constant;
}




//
// H+2j ME from Fabrizio Caola
// 

double  HJJMatEl(TVar::Process process, const TLorentzVector p[5], double Hggcoupl[3][2], double Hvvcoupl[3][2], TVar::VerbosityLevel verbosity)
{

  // by default assume only gg productions 
  // FOTRAN convention -5    -4   -3  -2    -1  0 1 2 3 4 5 
  //     parton flavor bbar cbar sbar ubar dbar g d u s c b
  // C++ convention     0     1   2    3    4   5 6 7 8 9 10
  //2-D matrix is reversed in fortran                                                                                                           
  // msq[ parton2 ] [ parton1 ]      
  //      flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];   
  double MatElsq[nmsq][nmsq];
  for ( int i = 0; i < nmsq; i++) {
    for ( int j = 0; j < nmsq; j++ ) {
      MatElsq[i][j] = 0;
    }
  }
  // input unit = GeV/100 such that 125GeV is 1.25 in the code
  // this needs to be applied for all the p4
  double p4[5][4];
  for (int i = 0; i < 5; i++) {
    p4[i][0] = p[i].Energy()/100.;
    p4[i][1] = p[i].Px()/100.;
    p4[i][2] = p[i].Py()/100.;
    p4[i][3] = p[i].Pz()/100.;

    // use out-going convention for the incoming particles
    if ( i < 2 ) {
      for ( int j = 0; j < 4; j++ ) {
	p4[i][j] = - p4[i][j];
      }
    }
  }      

  if ( verbosity >= TVar::DEBUG ) {
    std::cout << "p4[0] = "  << p4[0][0] << ", " <<  p4[0][1] << ", "  <<  p4[0][2] << ", "  <<  p4[0][3] << "\n";   
    std::cout << "p4[1] = "  << p4[1][0] << ", " <<  p4[1][1] << ", "  <<  p4[1][2] << ", "  <<  p4[1][3] << "\n";   
    std::cout << "p4[2] = "  << p4[2][0] << ", " <<  p4[2][1] << ", "  <<  p4[2][2] << ", "  <<  p4[2][3] << "\n"; 
    std::cout << "p4[3] = "  << p4[3][0] << ", " <<  p4[3][1] << ", "  <<  p4[3][2] << ", "  <<  p4[3][3] << "\n";   
    std::cout << "p4[4] = "  << p4[4][0] << ", " <<  p4[4][1] << ", "  <<  p4[4][2] << ", "  <<  p4[4][3] << "\n";   
  }

  
  if ( process == TVar::HJJNONVBF ) {
    __modhiggsjj_MOD_evalamp_gg_jjh(p4, Hggcoupl, MatElsq);
  }

  if ( process == TVar::HJJVBF ) {
    __modhiggsvbf_MOD_evalamp_vbfh(p4, Hvvcoupl, MatElsq);
  }

  //    FOTRAN convention    -5    -4   -3   -2   -1    0   1   2   3  4  5
  //     parton flavor      bbar  cbar  sbar ubar dbar  g   d   u   s  c  b
  //      C++ convention     0      1    2    3    4    5   6   7   8  9  10

  for(int ii = 0; ii < nmsq; ii++){
    for(int jj = 0; jj < nmsq; jj++){
      if ( verbosity >= TVar::DEBUG ) {
	std::cout<< "MatElsq: " << ii-5 << " " << jj-5 << " " << MatElsq[jj][ii] << "\n" ;
      }
    }
  }

  if ( process == TVar::HJJNONVBF ) {
    return SumMEPDF(p[0], p[1], MatElsq, verbosity);
    // return MatElsq[5][5];
  }
  
  if ( process == TVar::HJJVBF ) {
    // return MatElsq[6][7]+MatElsq[7][6];
    return SumMEPDF(p[0], p[1], MatElsq, verbosity);
  }

  return 0.;
}

// Below code sums over all production parton flavors according to PDF 
double SumMEPDF(const TLorentzVector p0, const TLorentzVector p1, double msq[nmsq][nmsq],  TVar::VerbosityLevel verbosity)
{
  //Calculate Pdf
  //Parton Density Function is always evalualted at pT=0 frame
  //Make sure parton Level Energy fraction is [0,1]
  //phase space function already makes sure the parton energy fraction between [min,1]
  //  x0 EBeam =>   <= -x1 EBeam
  double sysPz=p0.Pz()    + p1.Pz();
  double sysE =p0.Energy()+ p1.Energy();
  
  //Ignore the Pt doesn't make significant effect
  //double sysPt_sqr=sysPx*sysPx+sysPy*sysPy;
  //if(sysPt_sqr>=1.0E-10)  sysE=TMath::Sqrt(sysE*sysE-sysPt_sqr);
  double xx[2]={(sysE+sysPz)/EBEAM/2,(sysE-sysPz)/EBEAM/2};
  if ( verbosity >= TVar::DEBUG ) {
    std::cout << "xx[1]: " << xx[1] << "\t xx[2] = " << xx[2] << "\n";
  }

  if(xx[0] > 1.0 || xx[0]<=xmin_.xmin) return 0.0;
  if(xx[1] > 1.0 || xx[1]<=xmin_.xmin) return 0.0;
  double fx1[nmsq];
  double fx2[nmsq];

  //Always pass address through fortran function
  fdist_ (&density_.ih1, &xx[0], &scale_.scale, fx1); 
  fdist_ (&density_.ih2, &xx[1], &scale_.scale, fx2); 
  if ( verbosity >= TVar::DEBUG ) {
    for ( int i = 0; i < nmsq; i++ ) {
      std::cout << "fx1[" << i << "]: " <<  fx1[i] << "\t"
	"fx2[" << i << "]: " <<  fx2[i] << "\n";
    }
  }
  
  double msqjk(0.);
  double flavor_msq[nmsq][nmsq];
  
  for(int ii=0;ii<nmsq;ii++){
    for(int jj=0;jj<nmsq;jj++){
      flavor_msq[jj][ii] = 0.;
      //2-D matrix is reversed in fortran
      // msq[ parton2 ] [ parton1 ]
      flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];
      msqjk+=flavor_msq[jj][ii];
    }//ii
  }//jj
  
  return msqjk;

}
