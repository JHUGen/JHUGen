#include "TMath.h"
#include "TLorentzVector.h"
#include "TUtil.hh"
#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;

void SetEwkCouplingParameters(){

  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.aemmz_inp=7.81751E-03;
  ewinput_.wmass_inp=80.385;
  ewinput_.zmass_inp=91.1876;
  ewinput_.xw_inp=0.23116864; // Not used in the compiled MCFM ewcheme
/*
// SETTINGS TO MATCH JHUGen MEs:
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.aemmz_inp=7.81751E-03;
  ewinput_.wmass_inp=79.9549392; // Slightly differnt mW to set xW to match JHUGen
  ewinput_.zmass_inp=91.1876;
  ewinput_.xw_inp=0.23119;
*/
}

void SetAlphaS (double Q, int mynloop, int mynflav, string mypartons){
	if(Q<=1 || mynloop<=0 || mypartons.compare("Default")==0){
		if(Q<0) cout << "Invalid QCD scale for alpha_s, setting to mH/2..." << endl;
		Q = (masses_mcfm_.hmass)*0.5;
		mynloop = 1;
	};
	if(mypartons.compare("Default")!=0 && mypartons.compare("cteq6_l")!=0 && mypartons.compare("cteq6l1")!=0){
		cout << "Only default :: cteq6l1 or cteq6_l are supported. Modify mela.cc symlinks, put the pdf table into data/Pdfdata and retry. Setting to Default..." << endl;
		mypartons = "Default";
	};

	bool nflav_is_same = (nflav_.nflav == mynflav);
	scale_.scale = Q;
	scale_.musq = Q*Q;
	nlooprun_.nlooprun = mynloop;

// From pdfwrapper_linux.f:
	if(mypartons.compare("cteq6_l")==0) couple_.amz = 0.118;
	else if(mypartons.compare("cteq6l1")==0 || mypartons.compare("Default")==0) couple_.amz = 0.130;
	else couple_.amz = 0.118; // Add pdf as appropriate

// For proper pdfwrapper_linux.f execution (alpha_s computation does not use pdf but examines the pdf name to initialize amz.)
	if(!nflav_is_same){
		nflav_.nflav = mynflav;

		if(mypartons.compare("Default")!=0) sprintf(pdlabel_.pdlabel,"%s",mypartons.c_str());
		else sprintf(pdlabel_.pdlabel,"%s","cteq6l1"); // Default pdf is cteq6l1
		coupling2_();
	}
	else{
		qcdcouple_.as = alphas_(&(scale_.scale),&(couple_.amz),&(nlooprun_.nlooprun));
	};

	qcdcouple_.gsq = 4.0*TMath::Pi()*qcdcouple_.as;
	qcdcouple_.ason2pi = qcdcouple_.as/(2.0*TMath::Pi());
	qcdcouple_.ason4pi = qcdcouple_.as/(4.0*TMath::Pi());

// TEST RUNNING SCALE PER EVENT:
/*
	if(verbosity >= TVar::DEBUG){
		cout << "My pdf is: " << pdlabel_.pdlabel << endl;
		cout << "My Q: " << Q << " | My alpha_s: " << qcdcouple_.as << " at order " << nlooprun_.nlooprun << " with a(m_Z): " << couple_.amz << '\t'
			<< "Nflav: " << nflav_.nflav << endl;
*/
}


void SetMCFMHiggsDecayCouplings(bool useBSM, double Hvvcoupl[SIZE_HVV][2]){
	if (!useBSM){
		spinzerohiggs_anomcoupl_.AllowAnomalousCouplings = false;
		spinzerohiggs_anomcoupl_.ghz1[0] =  1; 
		spinzerohiggs_anomcoupl_.ghz2[0] =  0; 
		spinzerohiggs_anomcoupl_.ghz3[0] =  0;
		spinzerohiggs_anomcoupl_.ghz4[0] =  0;
		/*
		spinzerohiggs_anomcoupl_.ghzgs2[0] = 0; 
		spinzerohiggs_anomcoupl_.ghzgs3[0] = 0;  
		spinzerohiggs_anomcoupl_.ghzgs4[0] = 0;  
		spinzerohiggs_anomcoupl_.ghgsgs2[0] = 0;
		spinzerohiggs_anomcoupl_.ghgsgs3[0] = 0; 
		spinzerohiggs_anomcoupl_.ghgsgs4[0] = 0;       
		*/
		spinzerohiggs_anomcoupl_.ghz1_prime[0] = 0; 
		spinzerohiggs_anomcoupl_.ghz1_prime2[0] = 0; 
		spinzerohiggs_anomcoupl_.ghz1_prime3[0] = 0;
		spinzerohiggs_anomcoupl_.ghz1_prime4[0] = 0;
		spinzerohiggs_anomcoupl_.ghz1_prime5[0] = 0;
		spinzerohiggs_anomcoupl_.ghz1_prime6[0] = 0;
		spinzerohiggs_anomcoupl_.ghz1_prime7[0] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime[0] = 0; 
		spinzerohiggs_anomcoupl_.ghz2_prime2[0] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime3[0] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime4[0] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime5[0] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime6[0] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime7[0] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime[0] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime2[0] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime3[0] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime4[0] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime5[0] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime6[0] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime7[0] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime[0] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime2[0] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime3[0] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime4[0] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime5[0] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime6[0] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime7[0] = 0;
//		spinzerohiggs_anomcoupl_.ghzgs1_prime2[0] = 0;

		spinzerohiggs_anomcoupl_.ghz1[1] =  0; 
		spinzerohiggs_anomcoupl_.ghz2[1] =  0; 
		spinzerohiggs_anomcoupl_.ghz3[1] =  0;
		spinzerohiggs_anomcoupl_.ghz4[1] =  0;
		/*
		spinzerohiggs_anomcoupl_.ghzgs2[1] = 0; 
		spinzerohiggs_anomcoupl_.ghzgs3[1] = 0;  
		spinzerohiggs_anomcoupl_.ghzgs4[1] = 0;  
		spinzerohiggs_anomcoupl_.ghgsgs2[1] = 0;
		spinzerohiggs_anomcoupl_.ghgsgs3[1] = 0; 
		spinzerohiggs_anomcoupl_.ghgsgs4[1] = 0;       
		*/
		spinzerohiggs_anomcoupl_.ghz1_prime[1] = 0; 
		spinzerohiggs_anomcoupl_.ghz1_prime2[1] = 0; 
		spinzerohiggs_anomcoupl_.ghz1_prime3[1] = 0;
		spinzerohiggs_anomcoupl_.ghz1_prime4[1] = 0;
		spinzerohiggs_anomcoupl_.ghz1_prime5[1] = 0;
		spinzerohiggs_anomcoupl_.ghz1_prime6[1] = 0;
		spinzerohiggs_anomcoupl_.ghz1_prime7[1] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime[1] = 0; 
		spinzerohiggs_anomcoupl_.ghz2_prime2[1] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime3[1] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime4[1] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime5[1] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime6[1] = 0;
		spinzerohiggs_anomcoupl_.ghz2_prime7[1] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime[1] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime2[1] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime3[1] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime4[1] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime5[1] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime6[1] = 0;
		spinzerohiggs_anomcoupl_.ghz3_prime7[1] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime[1] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime2[1] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime3[1] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime4[1] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime5[1] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime6[1] = 0;
		spinzerohiggs_anomcoupl_.ghz4_prime7[1] = 0;
//		spinzerohiggs_anomcoupl_.ghzgs1_prime2[1] = 0;
//
	}
	else{
		spinzerohiggs_anomcoupl_.AllowAnomalousCouplings=true;
		spinzerohiggs_anomcoupl_.ghz1[0] =  Hvvcoupl[0][0]; 
		spinzerohiggs_anomcoupl_.ghz2[0] =  Hvvcoupl[1][0];
		spinzerohiggs_anomcoupl_.ghz3[0] =  Hvvcoupl[2][0];
		spinzerohiggs_anomcoupl_.ghz4[0] =  Hvvcoupl[3][0];
		/*
		spinzerohiggs_anomcoupl_.ghzgs2[0] = Hvvcoupl[4][0]; 
		spinzerohiggs_anomcoupl_.ghzgs3[0] = Hvvcoupl[5][0]; 
		spinzerohiggs_anomcoupl_.ghzgs4[0] = Hvvcoupl[6][0]; 
		spinzerohiggs_anomcoupl_.ghgsgs2[0] = Hvvcoupl[7][0];
		spinzerohiggs_anomcoupl_.ghgsgs3[0] = Hvvcoupl[8][0];
		spinzerohiggs_anomcoupl_.ghgsgs4[0] = Hvvcoupl[9][0];      
		*/
		spinzerohiggs_anomcoupl_.ghz1_prime[0] = Hvvcoupl[10][0];
		spinzerohiggs_anomcoupl_.ghz1_prime2[0] = Hvvcoupl[11][0];
		spinzerohiggs_anomcoupl_.ghz1_prime3[0] = Hvvcoupl[12][0];
		spinzerohiggs_anomcoupl_.ghz1_prime4[0] = Hvvcoupl[13][0];
		spinzerohiggs_anomcoupl_.ghz1_prime5[0] = Hvvcoupl[14][0];
		spinzerohiggs_anomcoupl_.ghz2_prime[0] = Hvvcoupl[15][0];
		spinzerohiggs_anomcoupl_.ghz2_prime2[0] = Hvvcoupl[16][0];
		spinzerohiggs_anomcoupl_.ghz2_prime3[0] = Hvvcoupl[17][0];
		spinzerohiggs_anomcoupl_.ghz2_prime4[0] = Hvvcoupl[18][0];
		spinzerohiggs_anomcoupl_.ghz2_prime5[0] = Hvvcoupl[19][0];
		spinzerohiggs_anomcoupl_.ghz3_prime[0] = Hvvcoupl[20][0];
		spinzerohiggs_anomcoupl_.ghz3_prime2[0] = Hvvcoupl[21][0];
		spinzerohiggs_anomcoupl_.ghz3_prime3[0] = Hvvcoupl[22][0];
		spinzerohiggs_anomcoupl_.ghz3_prime4[0] = Hvvcoupl[23][0];
		spinzerohiggs_anomcoupl_.ghz3_prime5[0] = Hvvcoupl[24][0];
		spinzerohiggs_anomcoupl_.ghz4_prime[0] = Hvvcoupl[25][0];
		spinzerohiggs_anomcoupl_.ghz4_prime2[0] = Hvvcoupl[26][0];
		spinzerohiggs_anomcoupl_.ghz4_prime3[0] = Hvvcoupl[27][0];
		spinzerohiggs_anomcoupl_.ghz4_prime4[0] = Hvvcoupl[28][0];
		spinzerohiggs_anomcoupl_.ghz4_prime5[0] = Hvvcoupl[29][0];
//		spinzerohiggs_anomcoupl_.ghzgs1_prime2[0] = Hvvcoupl[30][0];
		spinzerohiggs_anomcoupl_.ghz1_prime6[0] = Hvvcoupl[31][0];
		spinzerohiggs_anomcoupl_.ghz1_prime7[0] = Hvvcoupl[32][0];
		spinzerohiggs_anomcoupl_.ghz2_prime6[0] = Hvvcoupl[33][0];
		spinzerohiggs_anomcoupl_.ghz2_prime7[0] = Hvvcoupl[34][0];
		spinzerohiggs_anomcoupl_.ghz3_prime6[0] = Hvvcoupl[35][0];
		spinzerohiggs_anomcoupl_.ghz3_prime7[0] = Hvvcoupl[36][0];
		spinzerohiggs_anomcoupl_.ghz4_prime6[0] = Hvvcoupl[37][0];
		spinzerohiggs_anomcoupl_.ghz4_prime7[0] = Hvvcoupl[38][0];

		spinzerohiggs_anomcoupl_.ghz1[1] =  Hvvcoupl[0][1];
		spinzerohiggs_anomcoupl_.ghz2[1] =  Hvvcoupl[1][1];
		spinzerohiggs_anomcoupl_.ghz3[1] =  Hvvcoupl[2][1];
		spinzerohiggs_anomcoupl_.ghz4[1] =  Hvvcoupl[3][1];
		/*
		spinzerohiggs_anomcoupl_.ghzgs2[1] = Hvvcoupl[4][1]; 
		spinzerohiggs_anomcoupl_.ghzgs3[1] = Hvvcoupl[5][1]; 
		spinzerohiggs_anomcoupl_.ghzgs4[1] = Hvvcoupl[6][1];  
		spinzerohiggs_anomcoupl_.ghgsgs2[1] = Hvvcoupl[7][1];
		spinzerohiggs_anomcoupl_.ghgsgs3[1] = Hvvcoupl[8][1];
		spinzerohiggs_anomcoupl_.ghgsgs4[1] = Hvvcoupl[9][1];      
		*/
		spinzerohiggs_anomcoupl_.ghz1_prime[1] = Hvvcoupl[10][1];
		spinzerohiggs_anomcoupl_.ghz1_prime2[1] = Hvvcoupl[11][1];
		spinzerohiggs_anomcoupl_.ghz1_prime3[1] = Hvvcoupl[12][1];
		spinzerohiggs_anomcoupl_.ghz1_prime4[1] = Hvvcoupl[13][1];
		spinzerohiggs_anomcoupl_.ghz1_prime5[1] = Hvvcoupl[14][1];
		spinzerohiggs_anomcoupl_.ghz2_prime[1] = Hvvcoupl[15][1];
		spinzerohiggs_anomcoupl_.ghz2_prime2[1] = Hvvcoupl[16][1];
		spinzerohiggs_anomcoupl_.ghz2_prime3[1] = Hvvcoupl[17][1];
		spinzerohiggs_anomcoupl_.ghz2_prime4[1] = Hvvcoupl[18][1];
		spinzerohiggs_anomcoupl_.ghz2_prime5[1] = Hvvcoupl[19][1];
		spinzerohiggs_anomcoupl_.ghz3_prime[1] = Hvvcoupl[20][1];
		spinzerohiggs_anomcoupl_.ghz3_prime2[1] = Hvvcoupl[21][1];
		spinzerohiggs_anomcoupl_.ghz3_prime3[1] = Hvvcoupl[22][1];
		spinzerohiggs_anomcoupl_.ghz3_prime4[1] = Hvvcoupl[23][1];
		spinzerohiggs_anomcoupl_.ghz3_prime5[1] = Hvvcoupl[24][1];
		spinzerohiggs_anomcoupl_.ghz4_prime[1] = Hvvcoupl[25][1];
		spinzerohiggs_anomcoupl_.ghz4_prime2[1] = Hvvcoupl[26][1];
		spinzerohiggs_anomcoupl_.ghz4_prime3[1] = Hvvcoupl[27][1];
		spinzerohiggs_anomcoupl_.ghz4_prime4[1] = Hvvcoupl[28][1];
		spinzerohiggs_anomcoupl_.ghz4_prime5[1] = Hvvcoupl[29][1];
//		spinzerohiggs_anomcoupl_.ghzgs1_prime2[1] = Hvvcoupl[30][1];
		spinzerohiggs_anomcoupl_.ghz1_prime6[1] = Hvvcoupl[31][1];
		spinzerohiggs_anomcoupl_.ghz1_prime7[1] = Hvvcoupl[32][1];
		spinzerohiggs_anomcoupl_.ghz2_prime6[1] = Hvvcoupl[33][1];
		spinzerohiggs_anomcoupl_.ghz2_prime7[1] = Hvvcoupl[34][1];
		spinzerohiggs_anomcoupl_.ghz3_prime6[1] = Hvvcoupl[35][1];
		spinzerohiggs_anomcoupl_.ghz3_prime7[1] = Hvvcoupl[36][1];
		spinzerohiggs_anomcoupl_.ghz4_prime6[1] = Hvvcoupl[37][1];
		spinzerohiggs_anomcoupl_.ghz4_prime7[1] = Hvvcoupl[38][1];
//
	}
}



void My_choose(TVar::Process process, TVar::Production production, TVar::LeptonInterference leptonInterf, int flavor){
//void My_choose(TVar::Process process, TVar::Production production, int flavor){
 
//ZZ_4l
if(process==TVar::bkgZZ && (production == TVar::ZZQQB || production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU || production == TVar::ZZINDEPENDENT) ){ 
    //81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
    //86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'
    //90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))' 'L'
    //    nproc_.nproc=81;  
    //    chooser_();
  
  // these settings are identical to use the chooser_() function
    npart_.npart=4;
    nqcdjets_.nqcdjets=0;

    vsymfact_.vsymfact=1.0;                                                                                                               
    interference_.interference=false;

	if( (flavor == 1 || flavor == 0) && (leptonInterf==TVar::DefaultLeptonInterf || leptonInterf==TVar::InterfOn) ){
//	if( (flavor == 1 || flavor == 0) ){
    //90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))' 'L'
		vsymfact_.vsymfact=0.125; // MELA FACTOR (0.25 in MCFM 6.8)  --->   Result of just removing if(bw34_56) statements in FORTRAN code and not doing anything else
//		vsymfact_.vsymfact=0.25; // MELA FACTOR (Same 0.25 in MCFM 6.7)
		interference_.interference=true;
	}

    nwz_.nwz=0;
    bveg1_mcfm_.ndim=10;
    breit_.n2=1;
    breit_.n3=1;

    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.0; // Pretty important for MCFM 6.8+, does not matter for earlier versions
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.0; // Pretty important for MCFM 6.8+, does not matter for earlier versions
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;


} else if( ( (process==TVar::bkgZZ || process==TVar::HSMHiggs) && production == TVar::ZZGG) || process == TVar::bkgZZ_SMHiggs) {

    // 114 '  f(p1)+f(p2) --> H(--> Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6))' 'N'
/*
nprocs:
c--- 128 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [top, bottom loops, exact]' 'L'
c--- 129 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [only H, gg->ZZ intf.]' 'L' -> NOT IMPLEMENTED
c--- 130 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [H squared and H, gg->ZZ intf.]' 'L'
c--- 131 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6) [gg only, (H + gg->ZZ) squared]' 'L'
c--- 132 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6) [(gg->ZZ) squared]' 'L'
*/
     
     npart_.npart=4;
     nqcdjets_.nqcdjets=0;

     bveg1_mcfm_.ndim=10;

     vsymfact_.vsymfact=1.0;                                                                                                               
     interference_.interference=false;

	 if( (flavor == 1 || flavor == 0) && (leptonInterf==TVar::InterfOn) ){ // Notice, default lepton interf. is off for this calculation
	 	 vsymfact_.vsymfact=0.5;
		 interference_.interference=true;
 	 }

	 breit_.n2=1;
     breit_.n3=1;

     breit_.mass2 =masses_mcfm_.zmass;
     breit_.width2=masses_mcfm_.zwidth;
     breit_.mass3 =masses_mcfm_.zmass;
     breit_.width3=masses_mcfm_.zwidth;
     
    zcouple_.q1=-1.0;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.0;
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;

 } 
 else{
     std::cerr <<"[My_choose]: Can't identify Process: " << process <<endl;
 } 
}

bool My_masscuts(double s[][mxpart],TVar::Process process){

 double minZmassSqr=10*10;

 if(process==TVar::bkgZZ){
   if(s[2][3]< minZmassSqr) return true;
   if(s[4][5]< minZmassSqr) return true;
 }
 return false;	 

}


bool My_smalls(double s[][mxpart],int npart){

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
double SumMatrixElementPDF(TVar::Process process, TVar::Production production, TVar::MatrixElement myME, mcfm_event_type* mcfm_event,double flavor_msq[nmsq][nmsq],double* flux, double EBEAM, double coupling[SIZE_HVV_FREENORM]){

  int NPart=npart_.npart+2;
  double p4[4][mxpart];
  double s[mxpart][mxpart] = { { 0 } };
  double fx1[nmsq];
  double fx2[nmsq];
  double msq[nmsq][nmsq];
  int channeltoggle=0;
  
  
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
  double invariantP[5] = {0};
  //initialize decayed particles
  for(int ipar=2;ipar<NPart;ipar++){
    
    p4[0][ipar] = mcfm_event->p[ipar].Px();
    p4[1][ipar] = mcfm_event->p[ipar].Py();
    p4[2][ipar] = mcfm_event->p[ipar].Pz();
    p4[3][ipar] = mcfm_event->p[ipar].Energy();

	invariantP[1] += p4[3][ipar];
	invariantP[2] += p4[0][ipar];
	invariantP[3] += p4[1][ipar];
	invariantP[4] += p4[2][ipar];
  }

  invariantP[0] = pow(invariantP[1],2.0);
  for(int iq=2;iq<5;iq++) invariantP[0] -= pow(invariantP[iq],2.0);
  invariantP[0] = sqrt(fabs(invariantP[0]));

  double defaultScale = scale_.scale;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  SetAlphaS( invariantP[0]*0.5 , 1 , 5 , "cteq6_l"); // Set AlphaS(|Q|/2, mynloop, mynflav, mypartonPDF) for MCFM ME-related calculations

  //calculate invariant masses between partons/final state particles
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
  
  //Calculate Pdf
  //Always pass address through fortran function
  fdist_ (&density_.ih1, &xx[0], &scale_.scale, fx1); 
  fdist_ (&density_.ih2, &xx[1], &scale_.scale, fx2); 
/*
  if (process == TVar::bkgZZ && (production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU)){
	  if (production == TVar::ZZQQB_STU) cout << "STU" << endl;
	  if (production == TVar::ZZQQB_S) cout << "S" << endl;
	  if (production == TVar::ZZQQB_TU) cout << "TU" << endl;
  };
*/
  if( (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB) && process == TVar::bkgZZ)      qqb_zz_(p4[0],msq[0]);
  if( production == TVar::ZZQQB_STU && process == TVar::bkgZZ){
	  channeltoggle=0;
	  qqb_zz_stu_(p4[0], msq[0], &channeltoggle);
  };
  if( production == TVar::ZZQQB_S && process == TVar::bkgZZ){
	  channeltoggle=1;
	  qqb_zz_stu_(p4[0], msq[0], &channeltoggle);
  };
  if( production == TVar::ZZQQB_TU && process == TVar::bkgZZ){
	  channeltoggle=2;
	  qqb_zz_stu_(p4[0], msq[0], &channeltoggle);
  };
  //if( process==TVar::HZZ_4l)     qqb_hzz_(p4[0],msq[0]);
  // the subroutine for the calculations including the interfenrence             
  // ME =  sig + inter (sign, bkg)              
  // 1161 '  f(p1)+f(p2) --> H(--> Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6)) [including gg->ZZ intf.]' 'L'  
  if( process==TVar::bkgZZ_SMHiggs && myME==TVar::JHUGen)     gg_zz_int_freenorm_(p4[0],coupling,msq[0]); // |ggZZ + ggHZZ|**2 MCFM 6.6 version
  if( process==TVar::bkgZZ_SMHiggs && myME==TVar::MCFM)     gg_zz_all_ (p4[0],msq[0]); // |ggZZ + ggHZZ|**2
  if( process==TVar::HSMHiggs && production == TVar::ZZGG ) gg_hzz_tb_ (p4[0],msq[0]); // |ggHZZ|**2
  if( process==TVar::bkgZZ && production==TVar::ZZGG )    gg_zz_  (p4[0],&msqgg); // |ggZZ|**2

  /*
    // Below code sums over all production parton flavors according to PDF 
    // This is disabled as we are not using the intial production information
    // the below code is fine for the particle produced by single flavor of incoming partons

  for(int ii=0;ii<nmsq;ii++){
    for(int jj=0;jj<nmsq;jj++){
      
      //2-D matrix is transposed in fortran
      // msq[ parton2 ] [ parton1 ]
      //      flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];

      flavor_msq[jj][ii] = msq[jj][ii];
      //cout<<jj<<ii<<"="<<msq[jj][ii]<<"  ";
      msqjk+=flavor_msq[jj][ii];
    }//ii
    //    cout<<"\n";
  }//jj
  */
  // by default assume only gg productions 
  // FOTRAN convention -5    -4   -3  -2    -1  0 1 2 3 4 5 
  //     parton flavor bbar cbar sbar ubar dbar g d u s c b
  // C++ convention     0     1   2    3    4   5 6 7 8 9 10
  //
  msqjk=msq[5][5];
  if( process==TVar::bkgZZ && (production == TVar::ZZQQB || production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU || production ==TVar::ZZINDEPENDENT)) msqjk=msq[3][7]+msq[7][3];
/*
  if (process == TVar::bkgZZ && (production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU)){
	  for (int ix = 0; ix < 10; ix++){
		  for (int iy = 0; iy < 10; iy++) cout << msq[ix][iy] << '\t';
		  cout << endl;
	  }
  }
*/  // special for the GGZZ 
  if( process==TVar::bkgZZ && production == TVar::ZZGG ) msqjk=msqgg;      
  
  (*flux)=fbGeV2/(8*xx[0]*xx[1]*EBEAM*EBEAM);

  if(msqjk != msqjk || flux!=flux ){
    cout << "SumMatrixPDF: "<< TVar::ProcessName(process) << " msqjk="  << msqjk << " flux="<< *flux <<endl;
    msqjk=0;
    flux=0;
  }

  SetAlphaS( defaultScale , defaultNloop , defaultNflav , defaultPdflabel); // Protection for other probabilities
  return msqjk;
  
}

//
// Test code from Markus to calculate the HZZ cross-section
// 
double JHUGenMatEl(TVar::Process process, TVar::Production production, mcfm_event_type* mcfm_event, double MReso, double GaReso, 
		   double Hggcoupl[SIZE_HGG][2], double Hvvcoupl[SIZE_HVV][2], double Zqqcoupl[SIZE_ZQQ][2], double Zvvcoupl[SIZE_ZVV][2],
		   double Gqqcoupl[SIZE_GQQ][2], double Gggcoupl[SIZE_GGG][2], double Gvvcoupl[SIZE_GVV][2])
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

  //if ( process == TVar::HZZ_4l || process == TVar::PSHZZ_4l || process == TVar::HDHZZ_4l || process == TVar::HZZ_4l_MIXCP   || process == TVar::CPMixHZZ_4l || process == TVar::PSHZZ_g4star || process == TVar::HDHZZ_4l_g2star || process == TVar::HDMixHZZ_4l_pi_2|| process== TVar::HDMixHZZ_4l || process == TVar::CPMixHZZ_4l_pi_2 || process == TVar::SelfDefine || process == TVar::H_g1prime2) {
  if ( 
	  process == TVar::HSMHiggs
	  || process == TVar::H0minus
	  || process == TVar::H0hplus
	  || process == TVar::SelfDefine_spin0
	  || process == TVar::H0_g1prime2
	  || process== TVar::H0_Zgs
	  || process ==TVar::H0_gsgs
	  || process ==TVar::H0_Zgs_PS
	  || process ==TVar::H0_gsgs_PS
	  || process ==TVar::H0_Zgsg1prime2
	  ) {
    __modhiggs_MOD_evalamp_gg_h_vv(p4, &MReso,  &GaReso, Hggcoupl, Hvvcoupl, MYIDUP, &MatElSq);

  }
  if ( production == TVar::ZZGG ) {
    if ( process == TVar::H2_g1g5 || process == TVar::H2_g8 || process == TVar::H2_g4 || process == TVar::H2_g5 || process == TVar::H2_g2 || process == TVar::H2_g3 || process == TVar::H2_g6 || process == TVar::H2_g7 || process == TVar::H2_g9 || process == TVar::H2_g10 || process == TVar::SelfDefine_spin2) {
      __modgraviton_MOD_evalamp_gg_g_vv(p4, &MReso,  &GaReso, Gggcoupl, Gvvcoupl, MYIDUP, &MatElSq);
    }
  }
  if ( production == TVar::ZZINDEPENDENT ) {

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
    
    if ( process == TVar::H2_g1g5 || process == TVar::H2_g8 || process == TVar::H2_g4 || process == TVar::H2_g5 || process == TVar::H2_g2 || process == TVar::H2_g3 || process == TVar::H2_g6 || process == TVar::H2_g7 || process == TVar::H2_g9 || process == TVar::H2_g10 || process == TVar::SelfDefine_spin2) 
      __modgraviton_MOD_evalamp_g_vv(P, &MReso,  &GaReso, Gvvcoupl, MYIDUP, &MatElSq);

    if ( process == TVar::H1minus || process == TVar::H1plus || process == TVar::SelfDefine_spin1)
      __modzprime_MOD_evalamp_zprime_vv(P, &MReso,  &GaReso, Zvvcoupl, MYIDUP, &MatElSq);
  } 
  
  if ( production == TVar::ZZQQB ) {
    if ( process == TVar::H2_g1g5 || process == TVar::H2_g8 || process == TVar::H2_g4 || process == TVar::H2_g5 || process == TVar::H2_g2 || process == TVar::H2_g3 || process == TVar::H2_g6 || process == TVar::H2_g7 || process == TVar::H2_g9 || process == TVar::H2_g10 || process == TVar::SelfDefine_spin2) {
      __modgraviton_MOD_evalamp_qqb_g_vv(p4, &MReso,  &GaReso, Gqqcoupl, Gvvcoupl, MYIDUP, &MatElSq);
    }
    if ( process == TVar::H1minus || process == TVar::H1plus || process == TVar::SelfDefine_spin1) {
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
  double constant = 1.45e-8;
  return MatElSq*constant;

}


double HJJMatEl(TVar::Process process, TVar::Production production, const TLorentzVector p[5], double Hggcoupl[SIZE_HGG][2], double Hvvcoupl[SIZE_HVV_VBF][2], double Hwwcoupl[SIZE_HWW_VBF][2], TVar::VerbosityLevel verbosity, double EBEAM)
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

    // DO NOT Use out-going convention for the incoming particles for SumMEPDF
	// For HJJ, the subroutine already does its own calculation for p1, p2, so it does not matter if the sign remains flipped or not.
	// HJ exclusively takes lab-frame momenta, and p1 and p2 are used with a (-) sign. Thus, the sign would need to be flipped again.
/*
	if ( i < 2 ) {
		for ( int j = 0; j < 4; j++ ) {
			p4[i][j] = - p4[i][j];
		}
    }
*/
  }      
  if ( verbosity >= TVar::DEBUG ) {
    std::cout << "p4[0] = "  << p4[0][0] << ", " <<  p4[0][1] << ", "  <<  p4[0][2] << ", "  <<  p4[0][3] << "\n";   
    std::cout << "p4[1] = "  << p4[1][0] << ", " <<  p4[1][1] << ", "  <<  p4[1][2] << ", "  <<  p4[1][3] << "\n";   
    std::cout << "p4[2] = "  << p4[2][0] << ", " <<  p4[2][1] << ", "  <<  p4[2][2] << ", "  <<  p4[2][3] << "\n"; 
    std::cout << "p4[3] = "  << p4[3][0] << ", " <<  p4[3][1] << ", "  <<  p4[3][2] << ", "  <<  p4[3][3] << "\n";   
    std::cout << "p4[4] = "  << p4[4][0] << ", " <<  p4[4][1] << ", "  <<  p4[4][2] << ", "  <<  p4[4][3] << "\n";   
  }

  
  if ( production == TVar::JJGG ) {
    __modhiggsjj_MOD_evalamp_sbfh(p4, Hggcoupl, MatElsq);
  }
  if ( production == TVar::JJVBF) {
    __modhiggsjj_MOD_evalamp_wbfh(p4, Hvvcoupl, Hwwcoupl, MatElsq);
  }
  if ( production == TVar::JH) {
	double pOneJet[4][4] = { { 0 } };
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
//			if( i<2 ) pOneJet[i][j] = - p4[i][j]; // Revert back to lab-frame momenta
//			else pOneJet[i][j] = p4[i][j];
			pOneJet[i][j] = p4[i][j]; // Revert back to lab-frame momenta
		}
	}
	__modhiggsj_MOD_evalamp_hj(pOneJet, MatElsq);
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
  
  if ( production == TVar::JJGG || production == TVar::JJVBF || production == TVar::JH){
	return SumMEPDF(p[0], p[1], MatElsq, verbosity, EBEAM);
    //return MatElsq[5][5]; // jjgg
    //return MatElsq[6][7]+MatElsq[7][6]; // jjvbf
  }

  return 0.;
}

double VHiggsMatEl(TVar::Process process, TVar::Production production, TLorentzVector p[5], TLorentzVector pHdaughter[4], int Vdecay_id[6], double MReso, double GaReso, double Hvvcoupl[SIZE_HVV_VBF][2], TVar::VerbosityLevel verbosity, double EBEAM)
{
// Inputs to Fortran
// The dimensionality [9] should change to [11] once H->4l or 2l2nu is added to the amplitude.
  double p4[9][4] = { { 0 } };
  double masses[2][9] = { { 0 } };
  double helicities[9] = { 0 };
  int vh_ids[9] = { 0 };
  int n_HiggsFermions=0;

	MReso /=100.0;
	GaReso /= 100.0;

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
/*
	0 + 1 -> 2 (V*) -> 3 (V->5+6) + 4 (H->7+8)
	
	p[0]:=0
	p[1]:=1
	p[2]:=4 (H)
	p[3]:=5
	p[4]:=6
*/
  TLorentzVector pVH[9];
//  TLorentzVector pVH[11];
  TLorentzVector nullVector(0,0,0,0);
  for (int i = 0; i < 2; i++) pVH[i] = p[i];
  pVH[2] = p[0] + p[1]; // V*
  pVH[4] = p[2]; // H
  pVH[3] = pVH[2] - pVH[4]; // V
  pVH[5] = p[3]; // 5
  pVH[6] = p[4]; // 6
  if (
	  pHdaughter[0] != nullVector
	  && pHdaughter[1] != nullVector
	  && pHdaughter[2] == nullVector
	  && pHdaughter[3] == nullVector
	  ){
	  for (int i = 7; i < 9; i++) pVH[i] = pHdaughter[i-7]; // Higgs decay to 2 jets, variables.F90::H_DK=.true. should be set.
//	  for (int i = 9; i < 11; i++) pVH[i] = nullVector;
	  for (int i = 7; i < 9; i++) vh_ids[i] = Vdecay_id[i-5];
	  n_HiggsFermions=2;
  }
  else if (
	  pHdaughter[0] != nullVector
	  && pHdaughter[1] != nullVector
	  && pHdaughter[2] != nullVector
	  && pHdaughter[3] != nullVector
	  ){
	  for (int i = 7; i < 9; i++) pVH[i] = pHdaughter[2*(i-6)-2]+pHdaughter[2*(i-6)-1]; // Not yet supported fully
//	  for (int i = 9; i < 11; i++) pVH[i] = pHdaughter[i-7]; // Higgs decay to 4l, variables.F90::H_DK=.true. should be set.
	  for (int i = 7; i < 9; i++) vh_ids[i] = Vdecay_id[i-5];
//	  for (int i = 7; i < 11; i++) vh_ids[i] = Vdecay_id[i-5];
	  n_HiggsFermions=4;
  }
  else{
	  for (int i = 7; i < 9; i++) pVH[i] = pVH[4] * 0.5; // No Higgs decay is assumed, but conserve momentum in case of any mistake on variables.F90::H_DK=.false.
//	  for (int i = 7; i < 11; i++) pVH[i] = pVH[4] * 0.25; // No Higgs decay is assumed, but conserve momentum in case of any mistake on variables.F90::H_DK=.false.
  }

  // input unit = GeV/100 such that 125GeV is 1.25 in the code
  // this needs to be applied for all the p4
  for (int i = 0; i < 9; i++) {
    p4[i][0] = pVH[i].Energy()/100.;
    p4[i][1] = pVH[i].Px()/100.;
    p4[i][2] = pVH[i].Py()/100.;
    p4[i][3] = pVH[i].Pz()/100.;

    // DO NOT Use out-going convention for the incoming particles for SumMEPDF
	// VH exclusively takes lab-frame momenta, and p1 and p2 are used with a (-) sign.
/*
	if ( i < 2 ) {
		for ( int j = 0; j < 4; j++ ) {
			p4[i][j] = - p4[i][j];
		}
    }
*/
  }

// CAUTION: THESE HARDCODED NUMBERS HAVE TO BE THE SAME AS M/Ga_Z/W IN VARIABLES.F90
  if (production == TVar::ZH) {
	  masses[0][2] = 91.1876/100.;
	  masses[1][2] = 2.4952/100.;
	  masses[0][3] = masses[0][2];
	  masses[1][3] = masses[1][2];
  }
  if (production == TVar::WH) {
	  masses[0][2] = 80.399/100.;
	  masses[1][2] = 2.085/100.;
	  masses[0][3] = masses[0][2];
	  masses[1][3] = masses[1][2];
  }
  masses[0][4] = MReso; // Higgs
  masses[1][4] = GaReso; // Higgs

  vh_ids[4] = 25;
  vh_ids[5] = Vdecay_id[0]; // Handle jet-inclusive ME outside this function
  vh_ids[6] = Vdecay_id[1];
//  cout << "id5: " << vh_ids[5] << "\tid6: " << vh_ids[6] << endl;

  if ( verbosity >= TVar::DEBUG ) {
    for(int i=0;i<9;i++) std::cout << "p4[0] = "  << p4[i][0] << ", " <<  p4[i][1] << ", "  <<  p4[i][2] << ", "  <<  p4[i][3] << "\n";
    for(int i=0;i<9;i++) std::cout << "m(" << i << ") = "  << masses[0][i] << ", " <<  masses[0][i] << "\n";
//    for(int i=0;i<9;i++) std::cout << "id(" << i << ") = "  << vh_ids[i] << endl;
  }

  const double allowed_helicities[2] = { -1, 1 };
  double sumME=0;
  //    FOTRAN convention    -5    -4   -3   -2   -1    0   1   2   3  4  5
  //     parton flavor      bbar  cbar  sbar ubar dbar  g   d   u   s  c  b
  //      C++ convention     0      1    2    3    4    5   6   7   8  9  10
  for (int h0 = 0; h0 < 2; h0++){
	  helicities[0] = allowed_helicities[h0];
	  for (int h1 = 0; h1 < 2; h1++){
		  helicities[1] = allowed_helicities[h1];
		  for (int h5 = 0; h5 < 2; h5++){
			  helicities[5] = allowed_helicities[h5];
			  helicities[6] = -helicities[5];
			  for (int incoming1 = -nf; incoming1 <= nf; incoming1++){
				  if (production == TVar::ZH){
					  if(incoming1<=0) continue;
					  vh_ids[0] = incoming1;
					  vh_ids[1] = -incoming1;
					  vh_ids[2] = 23;
					  vh_ids[3] = 23;
					  double msq=0;
					  if(n_HiggsFermions==0) __modvhiggs_MOD_evalamp_vhiggs(vh_ids,helicities,p4,Hvvcoupl,masses,&msq);
					  else if(n_HiggsFermions==2){
						  for (int out1h = 0; out1h < 2; out1h++){
							  helicities[7] = allowed_helicities[out1h];
							  helicities[8] = helicities[7];
							  double msqtemp = 0;
							  __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, Hvvcoupl, masses, &msqtemp);
							  msq += msqtemp;
						  }
					  }
					  else if(n_HiggsFermions==4){
						  for (int out1h = 0; out1h < 2; out1h++){
							  for (int out2h = 0; out2h < 2; out2h++){
								  helicities[7] = allowed_helicities[out1h];
								  helicities[8] = allowed_helicities[out2h];
								  double msqtemp = 0;
								  __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, Hvvcoupl, masses, &msqtemp);
								  msq += msqtemp;
							  }
						  }
					  }
					  MatElsq[incoming1+5][-incoming1+5] += msq * 0.25; // Average over initial states with helicities +-1 only
					  MatElsq[-incoming1+5][incoming1+5] += msq * 0.25; // Average over initial states with helicities +-1 only
				  }
				  else if (production == TVar::WH){
					  if(incoming1==0) continue;

					  bool useWminus=false;
					  bool decaysToJets=false;

					  if( vh_ids[5] == -12 || vh_ids[5] == -14 || vh_ids[5] == -16 || vh_ids[5] == -2 || vh_ids[5] == -4 || vh_ids[6] == -2 || vh_ids[6] == -4 ) useWminus=true; // l- nu-bar or anti-up down -type quarks
					  if( abs(vh_ids[5])<=(nf+1) || abs(vh_ids[6])<=(nf+1) ) decaysToJets=true;

					  if (!useWminus){
						  if (incoming1 == 2 || incoming1 == 4){ // u or c to d-bar, b-bar or s-bar
							  for (int incoming2 = -nf; incoming2 < 0; incoming2++){
								  if( abs(incoming2)==abs(incoming1) ) continue;
								  vh_ids[0] = incoming1;
								  vh_ids[1] = incoming2;
								  vh_ids[2] = 24;
								  vh_ids[3] = 24;
								  double msq=0;
								  if(n_HiggsFermions==0) __modvhiggs_MOD_evalamp_vhiggs(vh_ids,helicities,p4,Hvvcoupl,masses,&msq);
								  else if(n_HiggsFermions==2){
									  for (int out1h = 0; out1h < 2; out1h++){
										  helicities[7] = allowed_helicities[out1h];
										  helicities[8] = helicities[7];
										  double msqtemp = 0;
										  __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, Hvvcoupl, masses, &msqtemp);
										  msq += msqtemp;
									  }
								  }
								  else if(n_HiggsFermions==4){
									  for (int out1h = 0; out1h < 2; out1h++){
										  for (int out2h = 0; out2h < 2; out2h++){
											  helicities[7] = allowed_helicities[out1h];
											  helicities[8] = allowed_helicities[out2h];
											  double msqtemp = 0;
											  __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, Hvvcoupl, masses, &msqtemp);
											  msq += msqtemp;
										  }
									  }
								  }
								  MatElsq[incoming1+5][incoming2+5] += msq * 0.25; // Average over initial states with helicities +-1 only
								  MatElsq[incoming2+5][incoming1+5] += msq * 0.25; // Average over initial states with helicities +-1 only
							  }
						  }
						  else continue;
					  }
					  else{
						  if (incoming1 == -2 || incoming1 == -4){ // u-bar or c-bar to d, b or s
							  for (int incoming2 = 1; incoming2 < nf+1; incoming2++){
								  if( abs(incoming2)==abs(incoming1) ) continue;
								  vh_ids[0] = incoming1;
								  vh_ids[1] = incoming2;
								  vh_ids[2] = 24;
								  vh_ids[3] = 24;
								  double msq=0;
								  if(n_HiggsFermions==0) __modvhiggs_MOD_evalamp_vhiggs(vh_ids,helicities,p4,Hvvcoupl,masses,&msq);
								  else if(n_HiggsFermions==2){
									  for (int out1h = 0; out1h < 2; out1h++){
										  helicities[7] = allowed_helicities[out1h];
										  helicities[8] = helicities[7];
										  double msqtemp = 0;
										  __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, Hvvcoupl, masses, &msqtemp);
										  msq += msqtemp;
									  }
								  }
								  else if(n_HiggsFermions==4){
									  for (int out1h = 0; out1h < 2; out1h++){
										  for (int out2h = 0; out2h < 2; out2h++){
											  helicities[7] = allowed_helicities[out1h];
											  helicities[8] = allowed_helicities[out2h];
											  double msqtemp = 0;
											  __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, Hvvcoupl, masses, &msqtemp);
											  msq += msqtemp;
										  }
									  }
								  }
								  MatElsq[incoming1+5][incoming2+5] += msq * 0.25; // Average over initial states with helicities +-1 only
								  MatElsq[incoming2+5][incoming1+5] += msq * 0.25; // Average over initial states with helicities +-1 only
							  }
						  }
						  else continue;
					  }
				  }
			  }
		  }
	  }
  }
	
  
  sumME = SumMEPDF(p[0], p[1], MatElsq, verbosity, EBEAM);
  return sumME;
}


// Below code sums over all production parton flavors according to PDF 
double SumMEPDF(const TLorentzVector p0, const TLorentzVector p1, double msq[nmsq][nmsq],  TVar::VerbosityLevel verbosity, double EBEAM)
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
    std::cout << "xx[0]: " << xx[0] << "\t xx[1] = " << xx[1] << "\n";
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
