//-----------------------------------------------------------------------------
//
// Class EventProb Module
//
//   EventProb Module
//
// March 21 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)
//-----------------------------------------------------------------------------

#include "TVar.hh"
#include "TEvtProb.hh"


ClassImp(TEvtProb)

using namespace std;

//-----------------------------------------------------------------------------
// Constructors and Destructor
//-----------------------------------------------------------------------------
TEvtProb::TEvtProb(double ebeam):EBEAM(ebeam){
  mcfm_init_((char *)"input.DAT",(char *)"./");
  SetEwkCouplingParameters();
  SetHiggsMass();
  energy_.sqrts = 2.*EBEAM;
  coupling_();
  SetLeptonInterf(TVar::DefaultLeptonInterf);
  spinzerohiggs_anomcoupl_.LambdaBSM=1000;
  spinzerohiggs_anomcoupl_.Lambda_z1=10000;
  spinzerohiggs_anomcoupl_.Lambda_z2=10000;
  spinzerohiggs_anomcoupl_.Lambda_z3=10000;
  spinzerohiggs_anomcoupl_.Lambda_z4=10000;
  spinzerohiggs_anomcoupl_.Lambda_Q=10000;
}


TEvtProb::~TEvtProb() {
}

/*
void TEvtProb::ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ){
  ewinput_.Gf_inp = ext_Gf;
  ewinput_.aemmz_inp = ext_aemmz;
  ewinput_.wmass_inp = ext_mW;
  ewinput_.zmass_inp = ext_mZ;
  ewinput_.xw_inp = 1.-pow(ext_mW/ext_mZ,2);
  coupling_();
}
*/

void TEvtProb::ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme){
  if (ext_ewscheme<-1 || ext_ewscheme>3) ext_ewscheme=3;
  ewinput_.Gf_inp = ext_Gf;
  ewinput_.aemmz_inp = ext_aemmz;
  ewinput_.wmass_inp = ext_mW;
  ewinput_.zmass_inp = ext_mZ;
  ewinput_.xw_inp = ext_xW;
  ewscheme_.ewscheme = ext_ewscheme;
  coupling_();
}

// Set NNPDF driver path
void TEvtProb::Set_LHAgrid(const char* path){
  char path_nnpdf_c[200];
  sprintf(path_nnpdf_c, "%s", path);
  int pathLength = strlen(path_nnpdf_c);
  nnpdfdriver_(path_nnpdf_c, &pathLength);
  pathLength=0;
  nninitpdf_(&pathLength);
}

//
// Directly calculate the ZZ->4l differential cross-section 
// WARNING: in development
// 
double TEvtProb::XsecCalc(TVar::Process proc, TVar::Production production, const hzz4l_event_type &hzz4l_event,
			TVar::VerbosityLevel verbosity,
			double couplingvals[SIZE_HVV_FREENORM],
			double selfDHvvcoupl[SIZE_HVV][2],
			double selfDZqqcoupl[SIZE_ZQQ][2],
			double selfDZvvcoupl[SIZE_ZVV][2],
			double selfDGqqcoupl[SIZE_GQQ][2],
			double selfDGggcoupl[SIZE_GGG][2],
			double selfDGvvcoupl[SIZE_GVV][2] ){

    //Initialize Process
    SetProcess(proc);
    SetProduction(production);
   
	int flavor = abs(hzz4l_event.PdgCode[0]) == abs(hzz4l_event.PdgCode[2]) ? 1 :3;
	bool needBSMHiggs=false;
	if (_matrixElement == TVar::MCFM || proc == TVar::bkgZZ_SMHiggs){ // Always uses MCFM
		for (int vv = 0; vv < SIZE_HVV; vv++){
//			if (selfDHvvcoupl[vv][1] != 0 || (vv != 0 && selfDHvvcoupl[vv][0] != 0)){ needBSMHiggs = true; break; }
			if ( selfDHvvcoupl[vv][1] != 0 || selfDHvvcoupl[vv][0] != 0 ){ needBSMHiggs = true; break; } // Only possible if selfDefine is called.
		}
		if (needBSMHiggs) SetLeptonInterf(TVar::InterfOn);
		SetMCFMHiggsDecayCouplings(needBSMHiggs, selfDHvvcoupl);
	}
    if ( _matrixElement == TVar::MCFM || proc == TVar::bkgZZ_SMHiggs) My_choose(proc, production, _leptonInterf, flavor);
    
    //constants
    double sqrts = 2.*EBEAM;
    double W=sqrts*sqrts;
    
    //Weight calculation
    double flux=1.;
    double dXsec=0.;
    
    mcfm_event_type mcfm_event; 
    // assign the right initial momentum
    // assumes the events are boosted to have 0 transverse momenta
    double sysPz= ( hzz4l_event.p[0] + hzz4l_event.p[1] + hzz4l_event.p[2] + hzz4l_event.p[3] ).Pz();
    double sysE = ( hzz4l_event.p[0] + hzz4l_event.p[1] + hzz4l_event.p[2] + hzz4l_event.p[3] ).Energy();
    double pz0 = (sysE+sysPz)/2.; 
    double pz1 = -(sysE-sysPz)/2.;
    mcfm_event.p[0].SetPxPyPzE   (0., 0., pz0, TMath::Abs(pz0));
    mcfm_event.p[1].SetPxPyPzE   (0., 0., pz1, TMath::Abs(pz1));
    mcfm_event.p[2].SetPxPyPzE   (hzz4l_event.p[0].Px(), hzz4l_event.p[0].Py(), hzz4l_event.p[0].Pz(), hzz4l_event.p[0].Energy());
    mcfm_event.p[3].SetPxPyPzE   (hzz4l_event.p[1].Px(), hzz4l_event.p[1].Py(), hzz4l_event.p[1].Pz(), hzz4l_event.p[1].Energy());
    mcfm_event.p[4].SetPxPyPzE   (hzz4l_event.p[2].Px(), hzz4l_event.p[2].Py(), hzz4l_event.p[2].Pz(), hzz4l_event.p[2].Energy());
    mcfm_event.p[5].SetPxPyPzE   (hzz4l_event.p[3].Px(), hzz4l_event.p[3].Py(), hzz4l_event.p[3].Pz(), hzz4l_event.p[3].Energy());
    
    mcfm_event.PdgCode[0] = 21;
    mcfm_event.PdgCode[1] = 21;
    mcfm_event.PdgCode[2] = hzz4l_event.PdgCode[0];
    mcfm_event.PdgCode[3] = hzz4l_event.PdgCode[1];
    mcfm_event.PdgCode[4] = hzz4l_event.PdgCode[2];
    mcfm_event.PdgCode[5] = hzz4l_event.PdgCode[3];

    //Matrix Element evaluation in qX=qY=0 frame
    //Evaluate f(x1)f(x2)|M(q)|/x1/x2 
    // 
    double qX=mcfm_event.p[0].Px()+mcfm_event.p[1].Px();
    double qY=mcfm_event.p[0].Py()+mcfm_event.p[1].Py();
    
    if((qX*qX+qY*qY)>0){
      double qE = mcfm_event.p[0].Energy()+mcfm_event.p[1].Energy();
      TVector3 boostV(qX/qE,qY/qE,0);
      for(int ipt=0;ipt<6;ipt++) mcfm_event.p[ipt].Boost(-boostV);
    }
    //event selections in Lab Frame
    double flavor_msq[nmsq][nmsq];
    double msqjk=0;
	if (_matrixElement == TVar::MCFM || proc == TVar::bkgZZ_SMHiggs){ // Always uses MCFM
		msqjk = SumMatrixElementPDF(_process, production, _matrixElement, &mcfm_event, flavor_msq, &flux, EBEAM, couplingvals);
		if(needBSMHiggs) SetLeptonInterf(TVar::DefaultLeptonInterf); // set defaults
		SetMCFMHiggsDecayCouplings(false, selfDHvvcoupl); // set defaults
	}
    else if ( _matrixElement == TVar::JHUGen ) {

	// all the possible couplings
      double Hggcoupl[SIZE_HGG][2] = { { 0 } };
      double Hvvcoupl[SIZE_HVV][2] = { { 0 } };
      double Zqqcoupl[SIZE_ZQQ][2] = { { 0 } };
      double Zvvcoupl[SIZE_ZVV][2] = { { 0 } };
      double Gqqcoupl[SIZE_GQQ][2] = { { 0 } };
      double Gggcoupl[SIZE_GGG][2] = { { 0 } };
      double Gvvcoupl[SIZE_GVV][2] = { { 0 } };
      
      // 
      // set spin 0 default numbers
      // 
      
      // By default set the Spin 0 couplings for SM case
      Hggcoupl[0][0]=1.0;  Hggcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
      Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
	  for (int i = 1; i < SIZE_HGG; i++){ for (int com = 0; com < 2; com++) Hggcoupl[i][com] = 0; }
      for (int i = 1; i<SIZE_HVV; i++){ for(int com=0; com<2; com++) Hvvcoupl[i][com] = 0; }

      // 
      // set spin 2 default numbers (2m+)
      // 
      Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0; // 2m+
      for (int i = 1; i<SIZE_GGG; i++){ for(int com=0; com<2; com++) Gggcoupl[i][com] = 0; }

      Gqqcoupl[0][0]=1.0;  Gqqcoupl[0][1]=0.0;  
      Gqqcoupl[1][0]=1.0;  Gqqcoupl[1][1]=0.0; 
      
      Gvvcoupl[0][0]=1.0;  Gvvcoupl[0][1]=0.0; // 2m+
      Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
      Gvvcoupl[4][0]=1.0;  Gvvcoupl[4][1]=0.0; // 2m+
      Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
      Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

      // 
      // set spin 1 default numbers (1-)
      // 
      Zqqcoupl[0][0]=1.0;  Zqqcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
      Zqqcoupl[1][0]=1.0;  Zqqcoupl[1][1]=0.0;
      // z->vv coupling constants
      Zvvcoupl[0][0]=1.0;  Zvvcoupl[0][1]=0.0; // 1-
      for (int i = 1; i<SIZE_ZVV; i++){ for(int com=0; com<2; com++) Zvvcoupl[i][com] = 0; }

      // 0-
      if ( proc == TVar::H0minus) {
		Hvvcoupl[0][0] = 0.0;
		Hvvcoupl[1][0] = 0.0;
		Hvvcoupl[2][0] = 0.0;
		Hvvcoupl[3][0] = 1.0;
      }

      if ( proc == TVar::H0_Zgs) {
				Hvvcoupl[4][0] = 0.0688;
				Hvvcoupl[0][0] = 0.;
      }
      if ( proc == TVar::H0_gsgs) {
				Hvvcoupl[7][0] = -0.0898;
				Hvvcoupl[0][0] = 0.;
      }
      if ( proc == TVar::H0_Zgs_PS) {
				Hvvcoupl[6][0] = 0.0855;
				Hvvcoupl[0][0] = 0.;
      }
      if ( proc == TVar::H0_gsgs_PS) {
				Hvvcoupl[9][0] = -0.0907;
				Hvvcoupl[0][0] = 0.;
      }
      if ( proc == TVar::H0_Zgsg1prime2) {
				Hvvcoupl[30][0] = -7591.914; // +- 6.613
				Hvvcoupl[0][0] = 0.;
      }

		if ( proc == TVar::SelfDefine_spin0){
			for(int i=0; i<SIZE_HVV; i++){
				for(int j=0;j<2;j++){
					Hvvcoupl [i][j] = selfDHvvcoupl[i][j];
				}
			}
		}
		if ( proc == TVar::SelfDefine_spin1){
			for(int i=0; i<SIZE_ZVV; i++){
				for(int j=0;j<2;j++){
					Zvvcoupl [i][j] = selfDZvvcoupl[i][j];
				}
			}
		}
		if ( proc == TVar::SelfDefine_spin2){
			for(int i=0; i<SIZE_GGG; i++){
				for(int j=0;j<2;j++){
					Gggcoupl [i][j] = selfDGggcoupl[i][j];
				}
			}
			for(int i=0; i<SIZE_GVV; i++){
				for(int j=0;j<2;j++){
					Gvvcoupl [i][j] = selfDGvvcoupl[i][j];
				}
			}
		}
		  // 0h+
		if ( proc == TVar::H0hplus) {
			Hvvcoupl[0][0] = 0.0;
			Hvvcoupl[1][0] = 1.0;
			Hvvcoupl[2][0] = 0.0;
			Hvvcoupl[3][0] = 0.0;
		}
		if( proc == TVar::H0_g1prime2){
			Hvvcoupl[0][0] = 0.;
			Hvvcoupl[11][0] = -12046.01;
		}
		// 2h-
		if ( proc == TVar::H2_g8 ){
			// gg production coupling constants
			Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
			Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
			Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
			Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
			Gggcoupl[4][0]=1.0;  Gggcoupl[4][1]=0.0; // 2h-
	
			// Graviton->ZZ coupling constants 
			Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
			Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
			Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
			Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
			Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
			Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
			Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
			Gvvcoupl[7][0]=1.0;  Gvvcoupl[7][1]=0.0; // 2h-
			Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
			Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0; 
		}
      
		  // 2h+
		if ( proc == TVar::H2_g4 ){
			// gg production coupling constants
			Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
			Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
			Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
			Gggcoupl[3][0]=1.0;  Gggcoupl[3][1]=0.0; // 2h+
			Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;
	
			// Graviton->ZZ coupling constants 
			Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
			Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
			Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
			Gvvcoupl[3][0]=1.0;  Gvvcoupl[3][1]=0.0; // 2h+
			Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
			Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
			Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
			Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
			Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
			Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;
		}
      
		  // 2b+
		if ( proc == TVar::H2_g5 ){
			// gg production coupling constants
			Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0;  // 2b+
			Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
			Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
			Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
			Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;
	
			// Graviton->ZZ coupling constants 
			Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
			Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
			Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
			Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
			Gvvcoupl[4][0]=1.0;  Gvvcoupl[4][1]=0.0; // 2b+
			Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
			Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
			Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
			Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
			Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;
		}

		if ( proc == TVar::H2_g2 ){			
			// gg production coupling constants
			Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
			Gggcoupl[1][0]=1.0;  Gggcoupl[1][1]=0.0;
			Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
			Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
			Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;
			
			// Graviton->ZZ coupling constants
			Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
			Gvvcoupl[1][0]=1.0;  Gvvcoupl[1][1]=0.0;
			Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
			Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
			Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
			Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
			Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
			Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
			Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
			Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;
		}

		// 2h3plus
		if ( proc == TVar::H2_g3 ){
			// gg production coupling constants
			Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
			Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
			Gggcoupl[2][0]=1.0;  Gggcoupl[2][1]=0.0;
			Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
			Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;
			
			// Graviton->ZZ coupling constants
			Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
			Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
			Gvvcoupl[2][0]=1.0;  Gvvcoupl[2][1]=0.0;
			Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
			Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
			Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
			Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
			Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
			Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
			Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;
		}

		// 2h6+
		if ( proc == TVar::H2_g6 ){
			// gg production coupling constants
			Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0;
			Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
			Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
			Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
			Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;
			
			// Graviton->ZZ coupling constants
			Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
			Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
			Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
			Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
			Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
			Gvvcoupl[5][0]=1.0;  Gvvcoupl[5][1]=0.0;
			Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
			Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
			Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
			Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;
		}
		
			// 2h7plus
		if ( proc == TVar::H2_g7 ){
			// gg production coupling constants
			Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0;
			Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
			Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
			Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
			Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;
			
			// Graviton->ZZ coupling constants
			Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
			Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
			Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
			Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
			Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
			Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
			Gvvcoupl[6][0]=1.0;  Gvvcoupl[6][1]=0.0;
			Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
			Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
			Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;
		}
		
		// 2h9minus
		if ( proc == TVar::H2_g9 ){
			// gg production coupling constants
			Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
			Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
			Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
			Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
			Gggcoupl[4][0]=1.0;  Gggcoupl[4][1]=0.0;
			
			// Graviton->ZZ coupling constants
			Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
			Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
			Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
			Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
			Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
			Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
			Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
			Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
			Gvvcoupl[8][0]=1.0;  Gvvcoupl[8][1]=0.0;
			Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;
		}
		
		// 2h10minus
		if ( proc == TVar::H2_g10 ){
			// gg production coupling constants
			Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
			Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
			Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
			Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
			Gggcoupl[4][0]=1.0;  Gggcoupl[4][1]=0.0;
			
			// Graviton->ZZ coupling constants
			Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
			Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
			Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
			Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
			Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
			Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
			Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
			Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
			Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
			Gvvcoupl[9][0]=1.0;  Gvvcoupl[9][1]=0.0;
		}
		if ( proc == TVar::H1plus) {
			// z->vv coupling constants
			Zvvcoupl[0][0]=0.0;  Zvvcoupl[0][1]=0.0;
			Zvvcoupl[1][0]=1.0;  Zvvcoupl[1][1]=0.0; // 1+
		}
		//cout<<_hwidth<<endl; 
		msqjk = JHUGenMatEl(proc, production, &mcfm_event, _hmass, _hwidth, Hggcoupl, Hvvcoupl, Zqqcoupl, Zvvcoupl, Gqqcoupl, Gggcoupl, Gvvcoupl);

	} // end of JHUGen matrix element calculations
    
    if(msqjk<=0){ mcfm_event.pswt=0; }
    
    flux=fbGeV2/(mcfm_event.p[0].Energy()*mcfm_event.p[1].Energy())	/(4*W);
    //    dXsec=msqjk*flux;
    dXsec=msqjk;
    
    
    if (verbosity >= TVar::DEBUG){
		cout <<"Process " << TVar::ProcessName(proc)
			<< " TEvtProb::XsecCalc(): dXsec=" << dXsec
			<< " Msq=" << msqjk 
			<< " flux=" << flux 
			<< endl;
	}

	return dXsec;
}

// Cross-section calculations for H + 2 jets
double TEvtProb::XsecCalcXJJ(TVar::Process proc, TVar::Production production, TLorentzVector p4[3], TVar::VerbosityLevel verbosity, double selfDHggcoupl[SIZE_HGG][2], double selfDHvvcoupl[SIZE_HVV_VBF][2], double selfDHwwcoupl[SIZE_HWW_VBF][2]){
  
  // Initialize Process
  SetProcess(proc);
  SetProduction(production);
  //constants
  //double sqrts = 2.*EBEAM;
  //double W=sqrts*sqrts;
  
  // first/second number is the real/imaginary part  

  double Hggcoupl[SIZE_HGG][2];
  Hggcoupl[0][0]=1.0;  Hggcoupl[0][1]=0.0; // g2 
    for(int i = 1; i<SIZE_HGG; i++){
       for(int j=0; j<2; j++){
         Hggcoupl[i][j]=0;
       }
    }

  double Hvvcoupl[SIZE_HVV_VBF][2];
  Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0; // g1
    for(int i = 1; i<SIZE_HVV_VBF; i++){
       for(int j=0; j<2; j++){
         Hvvcoupl[i][j]=0;
       }
    }
  double Hwwcoupl[SIZE_HWW_VBF][2];
    for(int i = 0; i<SIZE_HWW_VBF; i++){ // All are to be initialized at 0
       for(int j=0; j<2; j++){
         Hwwcoupl[i][j]=0;
       }
    }
// 0-
	if( proc == TVar::H0minus){
		Hggcoupl[0][0] = 0.0;
		Hggcoupl[1][0] = 0.0;
		Hggcoupl[2][0] = 1.0;
		Hvvcoupl[0][0] = 0.0;
		Hvvcoupl[1][0] = 0.0;
		Hvvcoupl[2][0] = 0.0;
		Hvvcoupl[3][0] = 1.0;
	}
// 0h+
	if ( proc == TVar::H0hplus) { // No need to re-set ggcoupl
		Hvvcoupl[0][0] = 0.0;
		Hvvcoupl[1][0] = 1.0;
		Hvvcoupl[2][0] = 0.0;
		Hvvcoupl[3][0] = 0.0;
	}
	if( proc == TVar::H0_g1prime2){ // No need to re-set ggcoupl
		Hvvcoupl[0][0] = 0.;
		Hvvcoupl[5][0] = -12046.01;
	}

	if ( proc == TVar::SelfDefine_spin0){
		for(int j=0;j<2;j++){
			for(int i=0; i<SIZE_HGG; i++){
				Hggcoupl [i][j] = selfDHggcoupl[i][j];
			}
			for(int i=0; i<SIZE_HVV_VBF; i++){
				Hvvcoupl [i][j] = selfDHvvcoupl[i][j];
			}
			for(int i=0; i<SIZE_HWW_VBF; i++){
				Hwwcoupl [i][j] = selfDHwwcoupl[i][j];
			}
		}
	}

  // input kinematics 
  //  !----- p1 and p2 used to get hadronic s
  //  !----- P(p1)+P(p2) -> j(p3) + j(p4) + H(p5)
  // p[0] -> p1
  // p[1] -> p2
  // p[2] -> p3
  // p[3] -> p4
  // p[4] -> p5
  TLorentzVector p[5];
  p[2].SetPxPyPzE ( p4[0].Px(), p4[0].Py(), p4[0].Pz(), p4[0].E() );
  p[3].SetPxPyPzE ( p4[1].Px(), p4[1].Py(), p4[1].Pz(), p4[1].E() );
  p[4].SetPxPyPzE ( p4[2].Px(), p4[2].Py(), p4[2].Pz(), p4[2].E() );

	TLorentzVector pCoM = p[2] + p[3] + p[4];
	double qX = pCoM.X();
	double qY = pCoM.Y();
	double qE = pCoM.E();
	double qPt = (qX*qX+qY*qY);
	if ( (qPt)>0 ){
		TVector3 boostV(-qX/qE,-qY/qE,0); // Unit boost vector
		for(int ipt=2;ipt<5;ipt++) p[ipt].Boost(boostV);
	}

  // assign the right initial momentum
  // assumes the events are boosted to have 0 transverse momenta
  double sysPz= ( pCoM ).Pz(); 
  double sysE = ( pCoM ).Energy(); 
  double pz0 = (sysE+sysPz)/2.; 
  double pz1 = -(sysE-sysPz)/2.;
  p[0].SetPxPyPzE   (0., 0., pz0, TMath::Abs(pz0));
  p[1].SetPxPyPzE   (0., 0., pz1, TMath::Abs(pz1));
    
  
  // calculate the matrix element squared
  double dXsec = 0;
  dXsec = HJJMatEl(proc, production, p, Hggcoupl, Hvvcoupl, Hwwcoupl, verbosity, EBEAM);
  if (verbosity >= TVar::DEBUG)
    {
      std::cout <<"Process " << TVar::ProcessName(proc) << 
	" TEvtProb::XsecCalc(): dXsec=" << dXsec << "\n";
    }
  return dXsec;
}


// Cross-section calculations for H (SM) + 1 jet
double TEvtProb::XsecCalcXJ(TVar::Process proc, TVar::Production production, TLorentzVector p4[2], TVar::VerbosityLevel verbosity){

	  double Hggcoupl[SIZE_HGG][2] = { { 0 } };
	  double Hvvcoupl[SIZE_HVV_VBF][2] = { { 0 } };
	  double Hwwcoupl[SIZE_HWW_VBF][2] = { { 0 } };

	// Initialize Process
	SetProcess(proc);
	SetProduction(production);
	//constants
	//double sqrts = 2.*EBEAM;
	//double W=sqrts*sqrts;

	// first/second number is the real/imaginary part  


  // input kinematics 
  //  !----- p1 and p2 used to get hadronic s
  //  !----- P(p1)+P(p2) -> H(p3) + j(p4)
  // p[0] -> p1
  // p[1] -> p2
  // p[2] -> p3
  // p[3] -> p4
  // p[4] -> 0
  TLorentzVector p[5];
  for (int mom = 0; mom < 2; mom++){
	  p[mom+2].SetPxPyPzE(p4[mom].Px(), p4[mom].Py(), p4[mom].Pz(), p4[mom].E());
  }
  p[4].SetPxPyPzE(0,0,0,0);

	TLorentzVector pCoM = p[2] + p[3];
	double qX = pCoM.X();
	double qY = pCoM.Y();
	double qE = pCoM.E();
	double qPt = (qX*qX+qY*qY);
	if ( (qPt)>0 ){
		TVector3 boostV(-qX/qE,-qY/qE,0); // Unit boost vector
		for(int ipt=2;ipt<4;ipt++) p[ipt].Boost(boostV);
	}

  // assign the right initial momentum
  // assumes the events are boosted to have 0 transverse momenta
  double sysPz = ( pCoM ).Pz(); 
  double sysE = ( pCoM ).Energy(); 
  double pz0 = (sysE+sysPz)/2.; 
  double pz1 = -(sysE-sysPz)/2.;
  p[0].SetPxPyPzE   (0., 0., pz0, TMath::Abs(pz0));
  p[1].SetPxPyPzE   (0., 0., pz1, TMath::Abs(pz1));
    
  
  // calculate the matrix element squared
  double dXsec = 0;
  dXsec = HJJMatEl(proc, production, p, Hggcoupl, Hvvcoupl, Hwwcoupl, verbosity, EBEAM);
  if (verbosity >= TVar::DEBUG)
    {
      std::cout <<"Process " << TVar::ProcessName(proc) << 
	" TEvtProb::XsecCalc(): dXsec=" << dXsec << "\n";
    }
  return dXsec;
}


double TEvtProb::XsecCalc_VX(TVar::Process proc, TVar::Production production, vh_event_type &vh_event,
			TVar::VerbosityLevel verbosity,
			double selfDHvvcoupl[SIZE_HVV_VBF][2]
			){

    //Initialize Process
    SetProcess(proc);
    SetProduction(production);

	// 0: Higgs; 1,2: V-daughters
	// PDG ID for V daughters could be passed as 0. While mothers cannot really be gluons, specifying 0 will mean averaging over u,d,c,s,b.
	bool inclusiveHadronicJets = (vh_event.PdgCode[1]==0 || vh_event.PdgCode[2]==0);
	const double N_Q=5.;
    
    //constants
    double sqrts = 2.*EBEAM;
    double W=sqrts*sqrts;
    
    //Weight calculation
    double msqjk=0;
    
	if (
		(abs(vh_event.PdgCode[2]) == 12 ||
		abs(vh_event.PdgCode[2]) == 14 ||
		abs(vh_event.PdgCode[2]) == 16) && production == TVar::WH
		){ // First daughter of W has to be a neutrino

		vh_event.p[1] = vh_event.p[1] + vh_event.p[2];
		vh_event.p[2] = vh_event.p[1] - vh_event.p[2];
		vh_event.p[1] = vh_event.p[1] - vh_event.p[2];
		vh_event.PdgCode[1] = vh_event.PdgCode[1] + vh_event.PdgCode[2];
		vh_event.PdgCode[2] = vh_event.PdgCode[1] - vh_event.PdgCode[2];
		vh_event.PdgCode[1] = vh_event.PdgCode[1] - vh_event.PdgCode[2];
	}

	int Vdecay_id[6] = {
		vh_event.PdgCode[1], vh_event.PdgCode[2],
		vh_event.PdgCode_Hdecay[0],
		vh_event.PdgCode_Hdecay[1],
		vh_event.PdgCode_Hdecay[2],
		vh_event.PdgCode_Hdecay[3]
	};

	TLorentzVector pCoM = vh_event.p[0] + vh_event.p[1] + vh_event.p[2];
	double qX = pCoM.X();
	double qY = pCoM.Y();
	double qE = pCoM.E();
	double qPt = (qX*qX+qY*qY);
	if ( (qPt)>0 ){
		TVector3 boostV(-qX/qE,-qY/qE,0); // Unit boost vector
		for(int ipt=0;ipt<3;ipt++) vh_event.p[ipt].Boost(boostV);
		for(int ipt=0;ipt<4;ipt++) vh_event.pHdecay[ipt].Boost(boostV);
	}

	// assign the right initial momentum
    // assumes the events are boosted to have 0 transverse momenta
    double sysPz= ( vh_event.p[0] + vh_event.p[1] + vh_event.p[2] ).Pz();
    double sysE = ( vh_event.p[0] + vh_event.p[1] + vh_event.p[2] ).Energy();
    double pz0 = (sysE+sysPz)/2.; 
    double pz1 = -(sysE-sysPz)/2.;
	TLorentzVector pVH[5];
	TLorentzVector pHdaughter[4];
	pVH[2].SetPxPyPzE   (vh_event.p[0].Px(), vh_event.p[0].Py(), vh_event.p[0].Pz(), vh_event.p[0].Energy());
	pVH[3].SetPxPyPzE   (vh_event.p[1].Px(), vh_event.p[1].Py(), vh_event.p[1].Pz(), vh_event.p[1].Energy());
	pVH[4].SetPxPyPzE   (vh_event.p[2].Px(), vh_event.p[2].Py(), vh_event.p[2].Pz(), vh_event.p[2].Energy());
	pVH[0].SetPxPyPzE   (0, 0, pz0, TMath::Abs(pz0));
	pVH[1].SetPxPyPzE   (0, 0, pz1, TMath::Abs(pz1));
	for(int ipt=0;ipt<4;ipt++) pHdaughter[ipt]=vh_event.pHdecay[ipt];

	if ( _matrixElement == TVar::JHUGen ) {
	// Set Couplings at the HVV* vertex
      double Hvvcoupl[SIZE_HVV_VBF][2] = { { 0 } };

	  // By default set the Spin 0 couplings for SM case
      Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
      for (int i = 1; i<SIZE_HVV_VBF; i++){ for(int com=0; com<2; com++) Hvvcoupl[i][com] = 0; }

      // 0-
      if ( proc == TVar::H0minus) {
		Hvvcoupl[0][0] = 0.0;
		Hvvcoupl[1][0] = 0.0;
		Hvvcoupl[2][0] = 0.0;
		Hvvcoupl[3][0] = 1.0;
      }

		if ( proc == TVar::SelfDefine_spin0){
			for(int i=0; i<SIZE_HVV_VBF; i++){
				for(int j=0;j<2;j++){
					Hvvcoupl [i][j] = selfDHvvcoupl[i][j];
				}
			}
		}
		// 0h+
		if ( proc == TVar::H0hplus) {
			Hvvcoupl[0][0] = 0.0;
			Hvvcoupl[1][0] = 1.0;
			Hvvcoupl[2][0] = 0.0;
			Hvvcoupl[3][0] = 0.0;
		}
		if( proc == TVar::H0_g1prime2){
			Hvvcoupl[0][0] = 0.;
			Hvvcoupl[5][0] = -12046.01;
		}

		if(!inclusiveHadronicJets) msqjk = VHiggsMatEl(proc, production, pVH, pHdaughter, Vdecay_id, _hmass, _hwidth, Hvvcoupl, verbosity, EBEAM);
		else{
			for (int outgoing1 = -nf; outgoing1 <= nf; outgoing1++){
				if (production == TVar::ZH){
					if (outgoing1 <= 0) continue;
					Vdecay_id[0] = outgoing1;
					Vdecay_id[1] = -outgoing1;
					msqjk += (VHiggsMatEl(proc, production, pVH, pHdaughter, Vdecay_id, _hmass, _hwidth, Hvvcoupl, verbosity, EBEAM)) / N_Q; // Average over quark flavors
				}
				else if (production == TVar::WH){
					if (outgoing1 == 0) continue;
					if (outgoing1 == 2 || outgoing1 == 4){ // u or c to d-bar, b-bar or s-bar
						for (int outgoing2 = -nf; outgoing2 < 0; outgoing2++){
							if (abs(outgoing2) == abs(outgoing1)) continue;
							Vdecay_id[0] = outgoing1;
							Vdecay_id[1] = outgoing2;
							msqjk += (VHiggsMatEl(proc, production, pVH, pHdaughter, Vdecay_id, _hmass, _hwidth, Hvvcoupl, verbosity, EBEAM)) / 12.; // Average over quark flavors; CAUTION about 12: Depends on nf, (nf+1)*(nf-1)/2 or nf**2/2
						}
					}
					if (outgoing1 == -2 || outgoing1 == -4){ // u-bar or c-bar to d, b or s
						for (int outgoing2 = 1; outgoing2 < nf + 1; outgoing2++){
							if (abs(outgoing2) == abs(outgoing1)) continue;
							Vdecay_id[0] = outgoing1;
							Vdecay_id[1] = outgoing2;
							msqjk += (VHiggsMatEl(proc, production, pVH, pHdaughter, Vdecay_id, _hmass, _hwidth, Hvvcoupl, verbosity, EBEAM)) / 12.; // Average over quark flavors; CAUTION about 12: Depends on nf
						}
					}
				}
			}
		}
	} // end of JHUGen matrix element calculations
	else{
		return 0.; // Analytical not implemented yet
	}
	return msqjk;
}

// Cross-section calculations for ttbar -> H
double TEvtProb::XsecCalc_TTX(
  TVar::Process proc, TVar::Production production, 
  tth_event_type &tth_event,
  int topDecay, int topProcess,
  TVar::VerbosityLevel verbosity,
  double selfDHvvcoupl[SIZE_TTH][2]
  ){

// Set Couplings at the TTH vertex
  double Hvvcoupl[SIZE_TTH][2] ={ { 0 } };

  //Initialize Process
  SetProcess(proc);
  SetProduction(production);

  double msq=0;

  TLorentzVector pCoM(0, 0, 0, 0);
  const int nV = 7;
  for(int vv=0;vv<nV;vv++) pCoM = pCoM + tth_event.p[vv];
  double qX = pCoM.X();
  double qY = pCoM.Y();
  double qE = pCoM.E();
  double qPt = (qX*qX+qY*qY);
  if ((qPt)>0){
    TVector3 boostV(-qX/qE, -qY/qE, 0); // Unit boost vector
    for (int vv=0; vv<nV; vv++) tth_event.p[vv].Boost(boostV);
  }

  TLorentzVector p_tbar(0, 0, 0, 0);
  TLorentzVector p_t(0, 0, 0, 0);
  
  bool unknownTopFlavor=false;
  int indexTTBAR=0;
  if (
    (tth_event.PdgCode_tdecay[0][0]==0 || tth_event.PdgCode_tdecay[1][0]==0) ||
    (tth_event.PdgCode_tdecay[0][0]>0 && tth_event.PdgCode_tdecay[1][0]>0) ||
    (tth_event.PdgCode_tdecay[0][0]<0 && tth_event.PdgCode_tdecay[1][0]<0)
    ) unknownTopFlavor=true;
  else if (tth_event.PdgCode_tdecay[0][0]>0) indexTTBAR=1;

  for (int vv=1; vv<4; vv++) p_tbar = p_tbar + tth_event.p[vv + 3*indexTTBAR];
  for (int vv=1; vv<4; vv++) p_t = p_t + tth_event.p[vv + 3*(1-indexTTBAR)];

  TLorentzVector pTTH[11];
  pTTH[0].SetXYZT(0, 0, EBEAM, EBEAM);
  pTTH[1].SetXYZT(0, 0, -EBEAM, EBEAM);
  pTTH[2].SetXYZT(tth_event.p[0].X(), tth_event.p[0].Y(), tth_event.p[0].Z(), tth_event.p[0].T());
  pTTH[3].SetXYZT(p_tbar.X(), p_tbar.Y(), p_tbar.Z(), p_tbar.T());
  pTTH[4].SetXYZT(p_t.X(), p_t.Y(), p_t.Z(), p_t.T());
  for (int vv=1; vv<4; vv++) pTTH[vv+4].SetXYZT(tth_event.p[vv + 3*indexTTBAR].X(), tth_event.p[vv + 3*indexTTBAR].Y(), tth_event.p[vv + 3*indexTTBAR].Z(), tth_event.p[vv + 3*indexTTBAR].T());
  for (int vv=1; vv<4; vv++) pTTH[vv+7].SetXYZT(tth_event.p[vv + 3*(1-indexTTBAR)].X(), tth_event.p[vv + 3*(1-indexTTBAR)].Y(), tth_event.p[vv + 3*(1-indexTTBAR)].Z(), tth_event.p[vv + 3*(1-indexTTBAR)].T());


  if (_matrixElement == TVar::JHUGen){
    // By default set the Spin 0 couplings for SM case
    Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
    for (int i = 1; i<SIZE_TTH; i++){ for (int com=0; com<2; com++) Hvvcoupl[i][com] = 0; }

    // 0-
    if (proc == TVar::H0minus) {
      Hvvcoupl[0][0] = 0.0;
      Hvvcoupl[1][0] = 1.0;
    }

    if (proc == TVar::SelfDefine_spin0){
      for (int i=0; i<SIZE_TTH; i++){
        for (int j=0; j<2; j++){
          Hvvcoupl[i][j] = selfDHvvcoupl[i][j];
        }
      }
    }

    double topMass = 173.2;
    double topWidth = 2.;
    if (production == TVar::bbH){
      topMass = 4.75;
      topWidth = 0;
    }
    msq = TTHiggsMatEl(production, pTTH, _hmass, _hwidth, topMass, topWidth, Hvvcoupl, topDecay, topProcess, verbosity);
    if (unknownTopFlavor){
//      cout << "Unknown top flavor" << endl;
      pTTH[4].SetXYZT(p_tbar.X(), p_tbar.Y(), p_tbar.Z(), p_tbar.T());
      pTTH[3].SetXYZT(p_t.X(), p_t.Y(), p_t.Z(), p_t.T());
      for (int vv=1; vv<4; vv++) pTTH[vv+7].SetXYZT(tth_event.p[vv + 3*indexTTBAR].X(), tth_event.p[vv + 3*indexTTBAR].Y(), tth_event.p[vv + 3*indexTTBAR].Z(), tth_event.p[vv + 3*indexTTBAR].T());
      for (int vv=1; vv<4; vv++) pTTH[vv+4].SetXYZT(tth_event.p[vv + 3*(1-indexTTBAR)].X(), tth_event.p[vv + 3*(1-indexTTBAR)].Y(), tth_event.p[vv + 3*(1-indexTTBAR)].Z(), tth_event.p[vv + 3*(1-indexTTBAR)].T());
      msq += TTHiggsMatEl(production, pTTH, _hmass, _hwidth, topMass, topWidth, Hvvcoupl, topDecay, topProcess, verbosity);
      msq *= 0.5;
    }
  }

  return msq;
}



// this appears to be some kind of 
// way of setting MCFM parameters through
// an interface defined in TMCFM.hh
void TEvtProb::SetHiggsMass(double mass, float wHiggs){
	if (wHiggs < 0){
		masses_mcfm_.hwidth=4.07e-3;
		_hwidth = 4.07e-3;
	}
	else{
		masses_mcfm_.hwidth = wHiggs;
    	_hwidth = wHiggs; 
	}
	if (mass < 0){
		masses_mcfm_.hmass=125.;
		_hmass = 125.;
	}
	else{
		masses_mcfm_.hmass=mass;
		_hmass = mass;
	}
//	cout << "Set JHUGen Higgs mass width to: " << _hmass << ", " << _hwidth << endl;
//	cout << "Set MCFM Higgs mass width to: " << masses_mcfm_.hmass << ", " << masses_mcfm_.hwidth << endl;
}
