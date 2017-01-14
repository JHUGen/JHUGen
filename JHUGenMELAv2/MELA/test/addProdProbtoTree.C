#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include <vector>
#include "math.h"

using namespace RooFit;
bool includePathIsSet = false;

vector<TLorentzVector> Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi);

void addProdProbtoTree(){  
  char* inputFile = (char *)"/scratch0/hep/ianderso/CJLST/130720d/PRODFSR_8TeV/4mu/HZZ4lTree_ZZTo4mu";
  int flavor=1;
  int max=-1;
  int LHCsqrts=8;
  bool jhugen=false;

  gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/data/$SCRAM_ARCH/libmcfm.so");
  gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libZZMatrixElementMELA.so");

  // set up path for local cmssw includes
  // as well as roofit
  if (!includePathIsSet) {
    TString path = gSystem->GetIncludePath();
    path += "-I$CMSSW_BASE/src/ ";
    path += "-I$ROOFITSYS/include/ ";
    gSystem->SetIncludePath(path.Data());

    // this is awkward, but want to protect against 
    // multiple additions to the base include path
    // if this function is invoked multiple times per session
    includePathIsSet = true;
  }

  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MELA/interface/Mela.h+");

  Mela myMELA(LHCsqrts,flavor);

  char inputFileName[500];
  char outputFileName[500];
  sprintf(inputFileName,"%s.root",inputFile);
  sprintf(outputFileName,"%s_withProbabilities.root",inputFile);

  TFile* sigFile = new TFile(inputFileName);
  TTree* sigTree=0;
  if(sigFile)
    sigTree = (TTree*) sigFile->Get("SelectedTree");
  if(!sigTree){
    //2nd try with the name of data obs tree
    sigTree = (TTree*) sigFile->Get("data_obs");
    if(!sigTree){
      cout<<"ERROR could not find the tree!"<<endl;
      return;
    }
  }
  
  float m1,m2,h1,h2,hs,phi,phi1,mzz; 
  float jet1px,jet1py,jet1pz,jet1E;
  float jet2px,jet2py,jet2pz,jet2E;
  float ZZPx,ZZPy,ZZPz,ZZE,dR;
  vector<double> *JetPt=0;
  vector<double> *JetEta=0;
  vector<double> *JetPhi=0;
  vector<double> *JetMass=0;

  float dXsec_HJJ_JHU = 0.;
  float dXsec_HJJVBF_JHU = 0.;
  float dXsec_HJJVH_JHU = 0.;
  float dXsec_PSHJJ_JHU = 0.;
  float dXsec_PSHJJVBF_JHU = 0.;
  float dXsec_PSHJJVH_JHU = 0.;
  
  // -------- JHUGen TREES ---------------
  if ( sigTree->GetBranchStatus("Jet1Px") )
    sigTree->SetBranchAddress("Jet1Px",&jet1px);
  if ( sigTree->GetBranchStatus("Jet1Py") )
    sigTree->SetBranchAddress("Jet1Py",&jet1py);
  if ( sigTree->GetBranchStatus("Jet1Pz") )
    sigTree->SetBranchAddress("Jet1Pz",&jet1pz);
  if ( sigTree->GetBranchStatus("Jet1E") )
    sigTree->SetBranchAddress("Jet1E",&jet1E);
  if ( sigTree->GetBranchStatus("Jet2Px") )
    sigTree->SetBranchAddress("Jet2Px",&jet2px);
  if ( sigTree->GetBranchStatus("Jet2Py") )
    sigTree->SetBranchAddress("Jet2Py",&jet2py);
  if ( sigTree->GetBranchStatus("Jet2Pz") )
    sigTree->SetBranchAddress("Jet2Pz",&jet2pz);
  if ( sigTree->GetBranchStatus("Jet2E") )
    sigTree->SetBranchAddress("Jet2E",&jet2E);
  if ( sigTree->GetBranchStatus("ZZPx") )
    sigTree->SetBranchAddress("ZZPx",&ZZPx); 
  if ( sigTree->GetBranchStatus("ZZPy") )
    sigTree->SetBranchAddress("ZZPy",&ZZPy);
  if ( sigTree->GetBranchStatus("ZZPz") )
    sigTree->SetBranchAddress("ZZPz",&ZZPz);
  if ( sigTree->GetBranchStatus("ZZE") )
    sigTree->SetBranchAddress("ZZE",&ZZE); 
  if ( sigTree->GetBranchStatus("deltaR") )
    sigTree->SetBranchAddress("deltaR",&dR); 
  //----------- CJLST TREES ---------------
  if(!jhugen){
    if ( sigTree->GetBranchStatus("JetPt") ) 
      sigTree->SetBranchAddress( "JetPt"   , &JetPt);
    if ( sigTree->GetBranchStatus("JetEta") ) 
      sigTree->SetBranchAddress( "JetEta"   , &JetEta);
    if ( sigTree->GetBranchStatus("JetPhi") ) 
      sigTree->SetBranchAddress( "JetPhi"   , &JetPhi);
    if ( sigTree->GetBranchStatus("JetMass") ) 
      sigTree->SetBranchAddress( "JetMass"   , &JetMass);
    sigTree->SetBranchAddress( "Z1Mass"        , &m1      );   
    sigTree->SetBranchAddress( "Z2Mass"        , &m2      );   
    sigTree->SetBranchAddress( "helcosthetaZ1" , &h1      );   
    sigTree->SetBranchAddress( "helcosthetaZ2" , &h2      );   
    sigTree->SetBranchAddress( "costhetastar"  , &hs      );   
    sigTree->SetBranchAddress( "helphi"        , &phi     );   
    sigTree->SetBranchAddress( "phistarZ1"     , &phi1    );   
    sigTree->SetBranchAddress( "ZZMass"        , &mzz     );
  }
  //---------------------------------------
  
  TFile* newFile = new TFile(outputFileName,"RECREATE");
  TTree* newTree = sigTree->CloneTree(0);//new TTree("newTree","SelectedTree"); 
  
  newTree->Branch("dXsec_HJJ_JHU"   , &dXsec_HJJ_JHU   ,"dXsec_HJJ_JHU/F");
  newTree->Branch("dXsec_HJJVBF_JHU"   , &dXsec_HJJVBF_JHU   ,"dXsec_HJJVBF_JHU/F");
  newTree->Branch("dXsec_HJJVH_JHU"    , &dXsec_HJJVH_JHU    ,"dXsec_HJJVH_JHU/F");
  newTree->Branch("dXsec_PSHJJ_JHU"   , &dXsec_PSHJJ_JHU   ,"dXsec_PSHJJ_JHU/F");
  newTree->Branch("dXsec_PSHJJVBF_JHU"   , &dXsec_PSHJJVBF_JHU   ,"dXsec_PSHJJVBF_JHU/F");
  newTree->Branch("dXsec_PSHJJVH_JHU"    , &dXsec_PSHJJVH_JHU    ,"dXsec_PSHJJVH_JHU/F");

  if(max==-1) max=sigTree->GetEntries();

  cout<<sigTree->GetEntries()<<endl;

  for(int iEvt=0; iEvt<max; iEvt++){
    
    dXsec_HJJ_JHU = -99.;
    dXsec_HJJVBF_JHU = -99.;
    dXsec_HJJVH_JHU = -99.;
    dXsec_PSHJJ_JHU = -99.;
    dXsec_PSHJJVBF_JHU = -99.;
    dXsec_PSHJJVH_JHU = -99.;

    if(iEvt>=sigTree->GetEntries()) break;

    //if ( iEvt != 112442 ) continue;
    if(iEvt%1000==0) {
      cout<<"event: "<<iEvt<<endl;
    }
    sigTree->GetEntry(iEvt);

    TLorentzVector jet1,jet2,higgs;
    TLorentzVector p4[3];
    int NJets=0;
    double energy, p3sq, ratio, rdiff;
    if(jhugen){
      jet1.SetPxPyPzE(jet1px,jet1py,jet1pz,jet1E);
      jet2.SetPxPyPzE(jet2px,jet2py,jet2pz,jet2E);
      higgs.SetPxPyPzE(ZZPx,ZZPy,ZZPz,ZZE);
      p4[0]=jet1;
      p4[1]=jet2;
      for(int f=0;f<2;f++){
	energy = p4[f].Energy();
	p3sq = sqrt( p4[f].Px()*p4[f].Px() + p4[f].Py()*p4[f].Py() + p4[f].Pz()*p4[f].Pz()); 
	ratio = energy / p3sq; 
	p4[f].SetPxPyPzE ( p4[f].Px()*ratio, p4[f].Py()*ratio, p4[f].Pz()*ratio, energy);
	if(p4[f].Pt()>30. && fabs(p4[f].Eta())<4.7) NJets++;
      }
      rdiff=jet1.DeltaR(jet2);
      if(rdiff<0.5) NJets--;
      p4[2]=higgs;
    }
    if(!jhugen){
      //Gather jets
      double jetptc=0.;
      double jetetac=0.;
      vector<TLorentzVector> p;
      TLorentzVector jets[10];
      NJets=0;
      if ( JetPt != 0 ) {
	for (unsigned int k=0; k<JetPt->size();k++){
	  if (NJets==10) continue;
	  jetptc=JetPt->at(k);
	  jetetac=JetEta->at(k);
	  //cout<<jetptc<<" "<<jetetac<<endl;
	  jets[NJets].SetPtEtaPhiM(JetPt->at(k),JetEta->at(k),JetPhi->at(k),JetMass->at(k));
	  double energy = jets[NJets].Energy();
	  double p3sq = sqrt( jets[NJets].Px()*jets[NJets].Px() + jets[NJets].Py()*jets[NJets].Py() + jets[NJets].Pz()*jets[NJets].Pz()); 
	  double ratio = energy / p3sq; 
	  jets[NJets].SetPxPyPzE ( jets[NJets].Px()*ratio, jets[NJets].Py()*ratio, jets[NJets].Pz()*ratio, energy);
	  if (jetptc>30. && fabs(jetetac)<4.7){
	    NJets++;
	  }
	}
      }
      p4[0]=jets[0];
      p4[1]=jets[1];
      
      //Build Higgs object
      p=Calculate4Momentum(mzz,m1,m2,acos(hs),acos(h1),acos(h2),phi1,phi);
      TLorentzVector Z1_minus = p[0];
      TLorentzVector Z1_plus  = p[1];
      TLorentzVector Z2_minus = p[2];
      TLorentzVector Z2_plus  = p[3];
      p4[2] = p[0]+p[1]+p[2]+p[3];
    }

    if(NJets>1){
      myMELA.setProcess(TVar::HJJVBF,TVar::JHUGen,TVar::QQB);
      myMELA.computeProdP(p4[0],2,p4[1],2,p4[2],25,0.,0,dXsec_HJJVBF_JHU);
      myMELA.setProcess(TVar::PSHJJVBF,TVar::JHUGen,TVar::QQB);
      myMELA.computeProdP(p4[0],2,p4[1],2,p4[2],25,0.,0,dXsec_PSHJJVBF_JHU);
      myMELA.setProcess(TVar::HJJNONVBF,TVar::JHUGen,TVar::QQB);
      myMELA.computeProdP(p4[0],2,p4[1],2,p4[2],25,0.,0,dXsec_HJJ_JHU);
      myMELA.setProcess(TVar::PSHJJNONVBF,TVar::JHUGen,TVar::QQB);
      myMELA.computeProdP(p4[0],2,p4[1],2,p4[2],25,0.,0,dXsec_PSHJJ_JHU);
    }

    newTree->Fill();
    
  }
  
  newFile->cd();
  newTree->Write("SelectedTree"); 
  newFile->Close();
}

vector<TLorentzVector> Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi){
    double phi1,phi2;
    phi1=TMath::Pi()-Phi1;
    phi2=Phi1+Phi;

    double gamma1,gamma2,beta1,beta2;

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


    //fermion 4 vectors in vector boson rest frame

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
