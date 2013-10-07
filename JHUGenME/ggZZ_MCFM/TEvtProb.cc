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

#include "TEvtProb.hh"
#include "TVar.hh"


ClassImp(TEvtProb)

    using namespace std;

    //-----------------------------------------------------------------------------
    // Constructors and Destructor
    //-----------------------------------------------------------------------------
TEvtProb::TEvtProb() {
  mcfm_init_("input.DAT", "./");
  SetEwkCoupligParameters();
  coupling_();
}


TEvtProb::~TEvtProb() {
}

//
// Directly calculate the ZZ->4l differential cross-section 
// WARNING: in development
// 
double TEvtProb::XsecCalc(TVar::Process proc, TVar::Production production, const hzz4l_event_type &hzz4l_event,
			  TVar::VerbosityLevel verbosity){

  
    //Initialize Process

    SetProduction(production);

    TString meName = "MCFM";
    if (_matrixElement == TVar::JHUGen ) 
      meName = "JHUGen";

    if ( _matrixElement == TVar::MCFM) 
      My_choose(proc);
    
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

    /*
    for ( int i = 0; i < 6; i++ ) {
      std::cout << "Particle " << i << " (Px, Py, Pz, E): " <<  mcfm_event.p[i].Px() << ", " << mcfm_event.p[i].Py() 
		<< ", " << mcfm_event.p[i].Pz() << ", " << mcfm_event.p[i].Energy() <<  "\n";
    }
    */
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
    if ( _matrixElement == TVar::MCFM ) 
      msqjk = SumMatrixElementPDF(proc, &mcfm_event, flavor_msq, &flux);
    if ( _matrixElement == TVar::JHUGen ) {

      // all the possible couplings
      double Hggcoupl[3][2];
      double Hvvcoupl[4][2];
      double Zqqcoupl[2][2];
      double Zvvcoupl[2][2];
      double Gqqcoupl[2][2];
      double Gggcoupl[5][2];
      double Gvvcoupl[10][2];
      
      // 
      // set spin 0 default numbers
      // 
      
      // By default set the Spin 0 couplings for SM case
      Hggcoupl[0][0]=1.0;  Hggcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
      Hggcoupl[1][0]=0.0;  Hggcoupl[1][1]=0.0;  
      Hggcoupl[2][0]=0.0;  Hggcoupl[2][1]=0.0;  
      
      Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0;  
      Hvvcoupl[1][0]=0.0;  Hvvcoupl[1][1]=0.0;  
      Hvvcoupl[2][0]=0.0;  Hvvcoupl[2][1]=0.0;  
      Hvvcoupl[3][0]=0.0;  Hvvcoupl[3][1]=0.0;        
      
      // 
      // set spin 2 default numbers (2m+)
      // 
      Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0; // 2m+
      Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
      Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;

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
      Zvvcoupl[1][0]=0.0;  Zvvcoupl[1][1]=0.0; 
      

      // 0-
      if ( proc == TVar::PSHZZ_4l ) {
	Hvvcoupl[0][0] = 0.0;
	Hvvcoupl[1][0] = 0.0;
	Hvvcoupl[2][0] = 0.0;
	  Hvvcoupl[3][0] = 1.0;
      }
      // 0h+
      if ( proc == TVar::HDHZZ_4l ) {
	Hvvcoupl[0][0] = 0.0;
	Hvvcoupl[1][0] = 1.0;
	Hvvcoupl[2][0] = 0.0;
	Hvvcoupl[3][0] = 0.0;
      }
      // 0-mixture complex cp 
      // g1 = 1 
      // g2 = g4 = 1 + 2.5 i
      if ( proc == TVar::HZZ_4l_MIXCP ) {
	Hvvcoupl[0][0] = 1.0;
	Hvvcoupl[1][0] = 1.0;
	Hvvcoupl[2][0] = 0.0;
	Hvvcoupl[3][0] = 1.0;
 
	Hvvcoupl[0][1] = 0.0;
	Hvvcoupl[1][1] = 2.5;
	Hvvcoupl[2][1] = 0.0;
	Hvvcoupl[3][1] = 2.5;
      }
      // 2h-
      if ( proc == TVar::PTZZ_2hminus_4l ) {
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
      if ( proc == TVar::TZZ_2hplus_4l ) {
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
      if ( proc == TVar::TZZ_2bplus_4l ) {
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
      
      if ( proc == TVar::AVZZ_4l ) {
	// z->vv coupling constants
	Zvvcoupl[0][0]=0.0;  Zvvcoupl[0][1]=0.0;
	Zvvcoupl[1][0]=1.0;  Zvvcoupl[1][1]=0.0; // 1+
      }
      
      msqjk = JHUGenMatEl(proc, production, &mcfm_event, _hmass, _hwidth, Hggcoupl, Hvvcoupl, Zqqcoupl, Zvvcoupl, Gqqcoupl, Gggcoupl, Gvvcoupl);
      
    } // end of JHUGen matrix element calculations
    
    if(msqjk<=0){ mcfm_event.pswt=0; }
    
    flux=fbGeV2/(mcfm_event.p[0].Energy()*mcfm_event.p[1].Energy())	/(4*W);
    //    dXsec=msqjk*flux;
    dXsec=msqjk;
    
    if (verbosity >= TVar::INFO)
      {
	cout <<"Process " << TVar::ProcessName(proc) << 
	  " TEvtProb::XsecCalc(), using Method " << meName << ":  dXsec=" << dXsec
	     <<" Msq="<<msqjk 
	     <<endl;
      }

    return dXsec;

}

// Cross-section calculations for H+2j
double TEvtProb::XsecCalcXJJ(TVar::Process proc, TLorentzVector p4[3], TVar::VerbosityLevel verbosity){
  
  // Initialize Process
  SetProcess(proc);
  //constants
  //  double sqrts = 2.*EBEAM;
  //  double W=sqrts*sqrts;
  
  double Hggcoupl[3][2];
  // first/second number is the real/imaginary part	  
  Hggcoupl[0][0]=1.0;  Hggcoupl[0][1]=0.0; // g2 
  Hggcoupl[1][0]=0.0;  Hggcoupl[1][1]=0.0; // g3
  Hggcoupl[2][0]=0.0;  Hggcoupl[2][1]=0.0; // g4  

  double Hvvcoupl[4][2];
  // first/second number is the real/imaginary part
  Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0; // g1  
  Hvvcoupl[1][0]=0.0;  Hvvcoupl[1][1]=0.0; // g2 
  Hvvcoupl[2][0]=0.0;  Hvvcoupl[2][1]=0.0; // g3 
  Hvvcoupl[3][0]=0.0;  Hvvcoupl[3][1]=0.0; // g4   

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

  // assign the right initial momentum
  // assumes the events are boosted to have 0 transverse momenta
  double sysPz= ( p[2] + p[3] + p[4] ).Pz(); 
  double sysE = ( p[2] + p[3] + p[4] ).Energy(); 
  double pz0 = (sysE+sysPz)/2.; 
  double pz1 = -(sysE-sysPz)/2.;
  p[0].SetPxPyPzE   (0., 0., pz0, TMath::Abs(pz0));
  p[1].SetPxPyPzE   (0., 0., pz1, TMath::Abs(pz1));
    
  
  // calculate the matrix element squared
  double dXsec = 0;
  dXsec = HJJMatEl(proc, p, Hggcoupl, Hvvcoupl, verbosity);
  if (verbosity >= TVar::INFO)
    {
      std::cout <<"Process " << TVar::ProcessName(proc) << 
	" TEvtProb::XsecCalc(): dXsec=" << dXsec << "\n";
    }
  return dXsec;
}
// this appears to be some kind of 
// way of setting MCFM parameters through
// an interface defined in TMCFM.hh
void TEvtProb::SetHiggsMass(double mass){
    masses_mcfm_.hmass=mass;
    masses_mcfm_.hwidth=0.004; // use narrow width for 125
    _hmass = mass;
    _hwidth = 0.004;
}
