#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/src/computeAngles.h>

using namespace RooFit;

using namespace std;

namespace {
  const float Zmass = 91.1876;
  const float gamZ = 2.5;
  const float M_muon = 0.105658389;
  const float M_electron = 0.00051099907;
  const float M_tau = 1.777;
  const int PDG_electron=11,PDG_muon=13,PDG_tau=15;
}



float getMCFMMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], double ggcoupl[2],int useConstant=0){
	float myprob=1.0;
	int myflavor=-1;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=1;
			else myflavor=2;
	}
	else myflavor=3;

	if(myflavor>=0) myMela.computeP(angularOrdered[0],angularOrdered[1],angularOrdered[2],angularOrdered[3],
		angularOrdered[4],angularOrdered[5],angularOrdered[6],angularOrdered[7],
	    myflavor,
	    myprob,
		useConstant
		);
	return myprob;
};
float getJHUGenMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], double selfDHvvcoupl[SIZE_HVV][2]){
	float myprob=1.0;
	int myflavor=-1;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=1;
			else myflavor=2;
	}
	else myflavor=3;

	if(myflavor>=0) myMela.computeP(angularOrdered[0],angularOrdered[1],angularOrdered[2],angularOrdered[3],
		angularOrdered[4],angularOrdered[5],angularOrdered[6],angularOrdered[7],
	    myflavor,
	    selfDHvvcoupl,
	    myprob
		);
	return myprob;
};
float getSuperMELA(Mela& myMela, int lepId[4], float mZZ, TVar::SuperMelaSyst syst){
	float myprob=1.0;
	TVar::LeptonFlavor myflavor=TVar::Flavor_Dummy;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=TVar::Flavor_4e;
			else myflavor=TVar::Flavor_4mu;
	}
	else myflavor=TVar::Flavor_2e2mu;

	if(myflavor!=TVar::Flavor_Dummy) myMela.computePM4l(mZZ,
	    myflavor,
	    syst,
	    myprob
		);
	return myprob;
};

double eta_to_costheta ( double eta ) { return cos(atan(exp(-eta))*2.); }

void test_JVBF(int erg_tev=8, float mSample=0, bool isggH=false){
  float mPOLE=mSample;
  if (mPOLE<=0) mPOLE=125.6;
  float wPOLE=4.15e-3;
  char TREE_NAME[] = "SelectedTree";

  //	TVar::VerbosityLevel verbosity = TVar::INFO;

  Mela mela(erg_tev, mPOLE);

  TFile* foutput;
  if (!isggH){
    if (mSample>0) foutput = new TFile(Form("HZZ4lTree_jvbfMELA_H%.0f_%iTeV.root", mSample, erg_tev), "recreate");
    else foutput = new TFile(Form("HZZ4lTree_jvbfMELA_HAll_%iTeV.root", erg_tev), "recreate");
  }
  else{
    if (mSample>0) foutput = new TFile(Form("HZZ4lTree_jvbfMELA_ggH%.0f_%iTeV.root", mSample, erg_tev), "recreate");
    else foutput = new TFile(Form("HZZ4lTree_jvbfMELA_ggHAll_%iTeV.root", erg_tev), "recreate");
  }
  TLorentzVector nullFourVector(0, 0, 0, 0);

  float MC_weight_noxsec;

  float pjvbf_VAJHU;
  float pjvbf_VAJHU_first;
  float pjvbf_VAJHU_second;
  float phj_VAJHU_first;
  float phj_VAJHU_second;
  float pAux_vbf;
  float pAux_vbf_first;
  float pAux_vbf_second;

  float jet1Pt, jet2Pt;
  float jet1Eta, jet2Eta;
  float jet1Phi, jet2Phi;
  float jet1E, jet2E;
  float jet1Pt_Fake, jet2Pt_Fake;
  float jet1Eta_Fake, jet2Eta_Fake;
  float jet1Phi_Fake, jet2Phi_Fake;
  float jet1E_Fake, jet2E_Fake;
  float jet1px, jet1py, jet1pz;
  float jet2px, jet2py, jet2pz;
  float ZZPx, ZZPy, ZZPz, ZZE, dR;
  short NJets30;
  std::vector<double> * JetPt=0;
  std::vector<double> * JetEta=0;
  std::vector<double> * JetPhi=0;
  std::vector<double> * JetMass=0;
  std::vector<double> myJetPt;
  std::vector<double> myJetCosTheta;
  std::vector<double> myJetEta;
  std::vector<double> myJetPhi;
  std::vector<double> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;

  float ZZMass, ZZPt, ZZPhi, ZZEta;

  int GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id;

  TChain* tree = new TChain(TREE_NAME);
  char* user_folder[3]={ "4mu", "4e", "2mu2e" };
  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR";
  if (erg_tev==8) cinput_main.Append("_8TeV");

//  TString cinput_main = "/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_";
//  cinput_main.Append(Form("%iTeV", erg_tev));

  for (int ff=0; ff<3; ff++){
    if (!isggH){
      if (mSample>0) tree->Add(Form("%s/%s/HZZ4lTree_VBFH%.0f.root", cinput_main.Data(), user_folder[ff], mSample));
      else tree->Add(Form("%s/%s/HZZ4lTree_VBFH*.root", cinput_main.Data(), user_folder[ff]));
    }
    else{
      if (mSample>0) tree->Add(Form("%s/%s/HZZ4lTree_minloH%.0f.root", cinput_main.Data(), user_folder[ff], mSample));
      else tree->Add(Form("%s/%s/HZZ4lTree_minloH*.root", cinput_main.Data(), user_folder[ff]));
    }
  }
  tree->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tree->SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
  tree->SetBranchAddress("JetMass", &JetMass, &bJetMass);
  tree->SetBranchAddress("ZZMass", &ZZMass);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("MC_weight_noxsec", &MC_weight_noxsec);
  newtree->Branch("ZZMass", &ZZMass);
  newtree->Branch("pAux_vbf", &pAux_vbf);
  newtree->Branch("pAux_vbf_first", &pAux_vbf_first);
  newtree->Branch("pAux_vbf_second", &pAux_vbf_second);
  newtree->Branch("pjvbf_VAJHU", &pjvbf_VAJHU);
  newtree->Branch("pjvbf_VAJHU_first", &pjvbf_VAJHU_first);
  newtree->Branch("pjvbf_VAJHU_second", &pjvbf_VAJHU_second);
  newtree->Branch("phj_VAJHU_first", &phj_VAJHU_first);
  newtree->Branch("phj_VAJHU_second", &phj_VAJHU_second);
  newtree->Branch("NJets30", &NJets30);
  newtree->Branch("jet1Pt", &jet1Pt);
  newtree->Branch("jet1Eta", &jet1Eta);
  newtree->Branch("jet1Phi", &jet1Phi);
  newtree->Branch("jet1E", &jet1E);
  newtree->Branch("jet2Pt", &jet2Pt);
  newtree->Branch("jet2Eta", &jet2Eta);
  newtree->Branch("jet2Phi", &jet2Phi);
  newtree->Branch("jet2E", &jet2E);
  newtree->Branch("jet1_Fake_Pt", &jet1Pt_Fake);
  newtree->Branch("jet1_Fake_Eta", &jet1Eta_Fake);
  newtree->Branch("jet1_Fake_Phi", &jet1Phi_Fake);
  newtree->Branch("jet1_Fake_E", &jet1E_Fake);
  newtree->Branch("jet2_Fake_Pt", &jet2Pt_Fake);
  newtree->Branch("jet2_Fake_Eta", &jet2Eta_Fake);
  newtree->Branch("jet2_Fake_Phi", &jet2Phi_Fake);
  newtree->Branch("jet2_Fake_E", &jet2E_Fake);
  newtree->Branch("JetPt", &myJetPt);
//  newtree->Branch("JetCosTheta", &myJetCosTheta);
  newtree->Branch("JetPhi", &myJetPhi);
  newtree->Branch("JetMass", &myJetMass);


  int nEntries = tree->GetEntries();
  double selfDHggcoupl[SIZE_HGG][2] ={ { 0 } };
  double selfDHvvcoupl[SIZE_HVV_VBF][2] ={ { 0 } };
  double selfDHwwcoupl[SIZE_HWW_VBF][2] ={ { 0 } };
  double ggvvcoupl[2]={ 0, 0 };
  mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
  int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    pjvbf_VAJHU=-1;
    pjvbf_VAJHU_first=-1;
    pjvbf_VAJHU_second=-1;
    jet1Pt=-1;
    jet2Pt=-1;
    jet1Eta=0;
    jet2Eta=0;

    tree->GetEntry(ev);
    GenLep1Id=11;
    GenLep2Id=-11;
    GenLep3Id=11;
    GenLep4Id=-11;

    myJetPt.clear();
    myJetEta.clear();
    myJetPhi.clear();
    myJetMass.clear();
    myJetCosTheta.clear();

    TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
    TLorentzVector p4[3], jets[2];
    higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, ZZMass);
    for (int i = 0; i < NJets30; i++){
      myJetPt.push_back(JetPt->at(i));
      myJetEta.push_back(JetEta->at(i));
      myJetPhi.push_back(JetPhi->at(i));
      myJetMass.push_back(JetMass->at(i));
      myJetCosTheta.push_back(eta_to_costheta(JetEta->at(i)));
    }

    int filled = 0;
    if (myJetPt.size()>=1){

      jets[0].SetPxPyPzE(0, 0, 0, 1);
      jets[1].SetPxPyPzE(0, 0, 0, 1);
      for (int i = 0; i < myJetPt.size(); i++){
        jets[filled].SetPtEtaPhiM(myJetPt[i], myJetEta[i], myJetPhi[i], myJetMass[i]);

        if (filled==0){
          double jetE = jets[filled].Energy();
          double jetP = jets[filled].P();
          double ratio = (jetP>0 ? jetE/jetP : 1);
          ratio = 1.;
          jet1.SetPxPyPzE(jets[filled].Px()*ratio, jets[filled].Py()*ratio, jets[filled].Pz()*ratio, jetE);
          filled++;
          jet1Pt = jet1.Pt();
          jet1Eta = jet1.Eta();
          jet1Phi = jet1.Phi();
          jet1E = jet1.E();
          jet2.SetXYZT(0, 0, 0, 0);
        }
        else if(filled==1){
          double jetE = jets[filled].Energy();
          double jetP = jets[filled].P();
          double ratio = (jetP>0 ? jetE/jetP : 1);
          ratio = 1.;
          jet2.SetXYZT(jets[filled].Px()*ratio, jets[filled].Py()*ratio, jets[filled].Pz()*ratio, jetE);
          filled++;
          jet2Pt = jet2.Pt();
          jet2Eta = jet2.Eta();
          jet2Phi = jet2.Phi();
          jet2E = jet2.E();
        }
        else continue;

/*
        if (filled == 0){
          if (jets[filled].Pt()>jet1.Pt()){
            jet1.SetPxPyPzE(jets[filled].Px(), jets[filled].Py(), jets[filled].Pz(), jetE);
            jet1Pt = jet1.Pt();
          }
          if (i == myJetPt.size() - 1){
            filled++;
            i = 0;
          }
        }
        else{
          if (jets[filled].Pt()<jet1.Pt() && jets[filled].Pt()>jet2.Pt()){
            jet2.SetPxPyPzE(jets[filled].Px(), jets[filled].Py(), jets[filled].Pz(), jetE);
            jet2Pt = jet2.Pt();
          }
          if (i == myJetPt.size() - 1){
            filled++;
          }
        }
        if (filled == 2) break;
*/
      }
//cos(atan(exp(-jet2Eta))*2)

      if (filled == 2){
        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
        mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, pjvbf_VAJHU);
        mela.get_PAux(pAux_vbf);

        TLorentzVector pTotal;
        pTotal.SetXYZT(0, 0, 0, 0);
        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
        mela.computeProdP(jet1, 2, pTotal, 2, higgs, 25, nullFourVector, 0, pjvbf_VAJHU_first);
        mela.get_PAux(pAux_vbf_first);
        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JH);
        mela.computeProdP(jet1, 2, pTotal, 2, higgs, 25, nullFourVector, 0, phj_VAJHU_first);
        mela::computeFakeJet(jet1, higgs, pTotal);
        jet2Pt_Fake = pTotal.Pt();
        jet2Eta_Fake = pTotal.Eta();
        jet2Phi_Fake = pTotal.Phi();
        jet2E_Fake = pTotal.E();

        pTotal.SetXYZT(0, 0, 0, 0);
        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
        mela.computeProdP(pTotal, 2, jet2, 2, higgs, 25, nullFourVector, 0, pjvbf_VAJHU_second);
        mela.get_PAux(pAux_vbf_second);
        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JH);
        mela.computeProdP(pTotal, 2, jet2, 2, higgs, 25, nullFourVector, 0, phj_VAJHU_second);
        mela::computeFakeJet(jet2, higgs, pTotal);
        jet1Pt_Fake = pTotal.Pt();
        jet1Eta_Fake = pTotal.Eta();
        jet1Phi_Fake = pTotal.Phi();
        jet1E_Fake = pTotal.E();

        newtree->Fill();
      }
      else if (filled == 1){
/*
        TLorentzVector pTotal = higgs+jet1;
        pTotal.SetVect(-pTotal.Vect());
        pTotal.SetE(pTotal.P());
        jet2 = pTotal;

        jet2Pt = jet2.Pt();
        jet2Eta = jet2.Eta();
        jet2Phi = jet2.Phi();
        jet2E = jet2.E();


        jet1Pt_Fake = jet1.Pt();
        jet1Eta_Fake = jet1.Eta();
        jet1Phi_Fake = jet1.Phi();
        jet1E_Fake = jet1.E();
        jet2Pt_Fake = jet2.Pt();
        jet2Eta_Fake = jet2.Eta();
        jet2Phi_Fake = jet2.Phi();
        jet2E_Fake = jet2.E();

        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
        mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, pjvbf_VAJHU);
*/
//
        jet2.SetXYZT(0,0,0,0);
        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
        mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, pjvbf_VAJHU);
        mela.get_PAux(pAux_vbf);
        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JH);
        mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, phj_VAJHU_first);
        mela::computeFakeJet(jet1, higgs, jet2);
/*
        cout << "TEST:"
          << " Higgs Pz: " <<  higgs.Pz()
          << " Higgs P: " <<  higgs.P()
          << " Higgs E: " <<  higgs.T()
          << " Jet 1 Pz: " <<  jet1.Pz()
          << " Jet 1 P: " <<  jet1.P()
          << " Jet 1 E: " <<  jet1.T()
          << " Jet 2 Pz: " <<  jet2.Pz()
          << " Jet 2 P: " <<  jet2.P()
          << " Jet 2 E: " <<  jet2.T() << '\n' << endl;
*/

        jet2Pt = jet2.Pt();
        jet2Eta = jet2.Eta();
        jet2Phi = jet2.Phi();
        jet2E = jet2.E();

        jet1Pt_Fake = jet1.Pt();
        jet1Eta_Fake = jet1.Eta();
        jet1Phi_Fake = jet1.Phi();
        jet1E_Fake = jet1.E();
        jet2Pt_Fake = jet2.Pt();
        jet2Eta_Fake = jet2.Eta();
        jet2Phi_Fake = jet2.Phi();
        jet2E_Fake = jet2.E();
//
        pjvbf_VAJHU_first = pjvbf_VAJHU;
        pjvbf_VAJHU_second = pjvbf_VAJHU;
        pAux_vbf_first = pAux_vbf;
        pAux_vbf_second = pAux_vbf;
        phj_VAJHU_second = phj_VAJHU_first;

        newtree->Fill();
      }
    }
  }


  foutput->WriteTObject(newtree);
  delete newtree;
  foutput->Close();
  delete tree;
}

void fit_JVBF(int erg_tev=8, float mSample=0){
  float mPOLE=mSample;
  if (mPOLE<=0) mPOLE=125.6;
  float wPOLE=4.15e-3;
  char TREE_NAME[] = "TestTree";

  TString cinput;
  if (mSample>0) cinput = Form("HZZ4lTree_jvbfMELA_H%.0f_%iTeV.root", mSample, erg_tev);
  else cinput = Form("HZZ4lTree_jvbfMELA_HAll_%iTeV.root", erg_tev);

  TFile* foutput;
  if (mSample>0) foutput = new TFile(Form("jvbfMELA_Fits_H%.0f_%iTeV.root", mSample, erg_tev), "recreate");
  else foutput = new TFile(Form("jvbfMELA_fits_wide_%iTeV.root", erg_tev), "recreate");

  TLorentzVector nullFourVector(0, 0, 0, 0);

  float MC_weight_noxsec;

  float pjvbf_VAJHU;
  float pjvbf_VAJHU_first;
  float pjvbf_VAJHU_second;

  float jet1Pt, jet2Pt;
  float jet1Eta, jet2Eta;
  float jet1Phi, jet2Phi;
  float jet1E, jet2E;
  float jet1Pt_Fake, jet2Pt_Fake;
  float jet1Eta_Fake, jet2Eta_Fake;
  float jet1Phi_Fake, jet2Phi_Fake;
  float jet1E_Fake, jet2E_Fake;
  float jet1px, jet1py, jet1pz;
  float jet2px, jet2py, jet2pz;
  float ZZPx, ZZPy, ZZPz, ZZE, dR;
  short NJets30;

  float ZZMass;

  TChain* tree = new TChain(TREE_NAME);
  tree->Add(cinput);
  tree->SetBranchAddress("MC_weight_noxsec", &MC_weight_noxsec);
  tree->SetBranchAddress("ZZMass", &ZZMass);
  tree->SetBranchAddress("pjvbf_VAJHU", &pjvbf_VAJHU);
  tree->SetBranchAddress("pjvbf_VAJHU_first", &pjvbf_VAJHU_first);
  tree->SetBranchAddress("pjvbf_VAJHU_second", &pjvbf_VAJHU_second);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("jet1Pt", &jet1Pt);
  tree->SetBranchAddress("jet1Eta", &jet1Eta);
  tree->SetBranchAddress("jet1Phi", &jet1Phi);
  tree->SetBranchAddress("jet1E", &jet1E);
  tree->SetBranchAddress("jet2Pt", &jet2Pt);
  tree->SetBranchAddress("jet2Eta", &jet2Eta);
  tree->SetBranchAddress("jet2Phi", &jet2Phi);
  tree->SetBranchAddress("jet2E", &jet2E);
  tree->SetBranchAddress("jet1_Fake_Pt", &jet1Pt_Fake);
  tree->SetBranchAddress("jet1_Fake_Eta", &jet1Eta_Fake);
  tree->SetBranchAddress("jet1_Fake_Phi", &jet1Phi_Fake);
  tree->SetBranchAddress("jet1_Fake_E", &jet1E_Fake);
  tree->SetBranchAddress("jet2_Fake_Pt", &jet2Pt_Fake);
  tree->SetBranchAddress("jet2_Fake_Eta", &jet2Eta_Fake);
  tree->SetBranchAddress("jet2_Fake_Phi", &jet2Phi_Fake);
  tree->SetBranchAddress("jet2_Fake_E", &jet2E_Fake);


  const int nbinsx=6;
  double bins_mzz[nbinsx+1]={70, 140, 190, 240, 340, 520, 20000};
  double sum_w[nbinsx]={ 0 };
  double sum_mzz[nbinsx][2]={ { 0 } };
  const int nbins_eta = 8;
  double bins_eta[nbins_eta+1]={ 0, 0.75, 1.5, 2.25, 3, 4, 5.5, 10.5, 12 };

  TProfile* hProb[nbinsx];
  TH1D* hProb_all;
  TProfile* hEta[nbinsx+1];
  TH1D* hDist[nbinsx+1];
  for (int bin=0; bin<nbinsx+1; bin++){
    if (bin<nbinsx) hProb[bin] = new TProfile(Form("njets1_pjvbf_bin_%i", bin+1), Form("m_{4l}: [%.0f, %.0f] GeV", bins_mzz[bin], bins_mzz[bin+1]), nbins_eta, bins_eta);
    else hProb_all = new TH1D("njets1_pjvbf_all", Form("m_{4l}: [%.0f, %.0f] GeV", bins_mzz[0], bins_mzz[nbinsx]), nbins_eta, bins_eta);
    if (bin<nbinsx) hEta[bin] = new TProfile(Form("njets1_eta_bin_%i", bin+1), Form("m_{4l}: [%.0f, %.0f] GeV", bins_mzz[bin], bins_mzz[bin+1]), nbins_eta, bins_eta);
    else hEta[bin] = new TProfile("njets1_eta_all", Form("m_{4l}: [%.0f, %.0f] GeV", bins_mzz[0], bins_mzz[nbinsx]), nbins_eta, bins_eta);
    if (bin<nbinsx) hDist[bin] = new TH1D(Form("njets1_dist_bin_%i", bin+1), Form("m_{4l}: [%.0f, %.0f] GeV", bins_mzz[bin], bins_mzz[bin+1]), nbins_eta, bins_eta);
    else hDist[bin] = new TH1D("njets1_dist_all", Form("m_{4l}: [%.0f, %.0f] GeV", bins_mzz[0], bins_mzz[nbinsx]), nbins_eta, bins_eta);
    if (bin<nbinsx) hProb[bin]->Sumw2();
    else hProb_all->Sumw2();
    hEta[bin]->Sumw2();
    hDist[bin]->Sumw2();
  }

  int nEntries = tree->GetEntries();
  cout << nEntries << endl;
  for (int ev = 0; ev < nEntries; ev++){
    tree->GetEntry(ev);

    int massbin=-1;
    for (int bin=0; bin<nbinsx; bin++){
      if (ZZMass>=bins_mzz[bin] && ZZMass<bins_mzz[bin+1]){
        massbin=bin;
      }
    }
    if (NJets30==1){
      hProb[massbin]->Fill(fabs(jet2Eta_Fake), pjvbf_VAJHU_second, MC_weight_noxsec);
      hEta[massbin]->Fill(fabs(jet2Eta_Fake), fabs(jet2Eta_Fake), MC_weight_noxsec);
      hDist[massbin]->Fill(fabs(jet2Eta_Fake), MC_weight_noxsec);

      hEta[nbinsx]->Fill(fabs(jet2Eta_Fake), fabs(jet2Eta_Fake), MC_weight_noxsec);
      hDist[nbinsx]->Fill(fabs(jet2Eta_Fake), MC_weight_noxsec);
    }
  }

  for (int bin=0; bin<nbinsx; bin++){
    hProb[bin]->Scale(1./hProb[bin]->Integral("width"));
  }
  for (int bin=0; bin<nbinsx+1; bin++){
    for (int bineta=1; bineta<=nbins_eta; bineta++){
      double bincontent = hDist[bin]->GetBinContent(bineta);
      double binerror = hDist[bin]->GetBinError(bineta);
      double binwidth = hDist[bin]->GetXaxis()->GetBinWidth(bineta);

      bincontent /= binwidth;
      binerror /= binwidth;

      hDist[bin]->SetBinContent(bineta, bincontent);
      hDist[bin]->SetBinError(bineta, binerror);
    }
    hDist[bin]->Scale(1./hDist[bin]->Integral("width"));
  }
  for (int bineta=1; bineta<=nbins_eta; bineta++){
    double sum_invwsq=0;
    double sum_prob_invwsq=0;
    for (int bin=0; bin<nbinsx; bin++){
      double bincontent = hProb[bin]->GetBinContent(bineta);
      double binerror = hProb[bin]->GetBinError(bineta);

      double invwsq = 0, prob_invwsq = 0;
      if (binerror!=0){
        invwsq = 1./pow(binerror, 2);
        prob_invwsq = bincontent*invwsq;
      }
      sum_invwsq += invwsq;
      sum_prob_invwsq += prob_invwsq;
    }

    if (sum_invwsq>0){
      sum_prob_invwsq /= sum_invwsq;
      sum_invwsq = 1./sqrt(sum_invwsq);
    }
    hProb_all->SetBinContent(bineta, sum_prob_invwsq);
    hProb_all->SetBinError(bineta, sum_invwsq);
  }
  hProb_all->Scale(1./hProb_all->Integral("width"));


  double eta_center_all[nbins_eta];
  double prob_center_all[nbins_eta];
  double eta_center_all_err[nbins_eta];
  double prob_center_all_err[nbins_eta];
  int nNonEmpty_all = 0;
  for (int bineta=1; bineta<=nbins_eta; bineta++){
    if (hProb_all->GetBinError(bineta)>0){
      eta_center_all[nNonEmpty_all] = hEta[nbinsx]->GetBinContent(bineta);
      eta_center_all_err[nNonEmpty_all] = hEta[nbinsx]->GetBinError(bineta);
      prob_center_all[nNonEmpty_all] = hProb_all->GetBinContent(bineta);
      prob_center_all_err[nNonEmpty_all] = hProb_all->GetBinError(bineta);
      cout << nNonEmpty_all << '\t' << prob_center_all[nNonEmpty_all] << '\t' << prob_center_all_err[nNonEmpty_all] << endl;
      nNonEmpty_all++;
    }
  }

  TGraphErrors* tg_all = new TGraphErrors(nNonEmpty_all, eta_center_all, prob_center_all, eta_center_all_err, prob_center_all_err);
  tg_all->SetNameTitle("tg_njets1_pjvbf_all", Form("m_{4l}: [%.0f, %.0f] GeV", bins_mzz[0], bins_mzz[nbinsx]));
  tg_all->GetXaxis()->SetTitle("#eta");
  tg_all->GetYaxis()->SetTitle("<P_{VBF}>");
  
  TF1* fit_tg_all = new TF1("fit_tg_njets1_pjvbf_all", "gaus*(1+gaus(3))", bins_eta[0], bins_eta[nbins_eta]);
  const int npars = 6;
  double mypars[npars]={ 0.155, 5.193, 1.46, 0.362, 3.524, 0.715 };
  fit_tg_all->SetParameters(mypars);
  fit_tg_all->SetParLimits(0, 0, 1);
  fit_tg_all->SetParLimits(3, 0, 1);
  fit_tg_all->SetParLimits(1, 0, 20);
  fit_tg_all->SetParLimits(4, 0, 20);
  fit_tg_all->SetParLimits(2, 0, 20);
  fit_tg_all->SetParLimits(5, 0, 20);
  tg_all->Fit(fit_tg_all);

  double* par_postfit = fit_tg_all->GetParameters();
  double* parerr_postfit = fit_tg_all->GetParErrors();
  double parIndex[npars];
  double parIndexErr[npars]={ 0 };
  for (int ip=0; ip<npars; ip++) parIndex[ip] = (double)ip;

  double integral = 0.5*par_postfit[0]*(1.+TMath::Erf(par_postfit[1]/par_postfit[2]/TMath::Sqrt(2)));
  double integral2 = 1./(2.*TMath::Sqrt(2.*TMath::Pi()*(TMath::Power(par_postfit[2], 2)+TMath::Power(par_postfit[5], 2))));
  integral2 *= par_postfit[0]*par_postfit[3]*TMath::Exp(-0.5*pow(par_postfit[1]-par_postfit[4], 2)/(TMath::Power(par_postfit[2], 2)+TMath::Power(par_postfit[5], 2)))*
    (1.+TMath::Erf((TMath::Power(par_postfit[2], 2)*par_postfit[4]+TMath::Power(par_postfit[5], 2)*par_postfit[1]) / (par_postfit[2]*par_postfit[5]*TMath::Sqrt(2.*(TMath::Power(par_postfit[2], 2)+TMath::Power(par_postfit[5], 2))))));
  integral += integral2;

  cout << integral << endl;
  par_postfit[0] /= integral;
  parerr_postfit[0] /= integral;

  TGraphErrors* tg_pars = new TGraphErrors(npars, parIndex, par_postfit, parIndexErr, parerr_postfit);
  tg_pars->SetNameTitle("njets1_pjvbf_pars", "gaus*(1+gaus(3))");
  tg_pars->GetXaxis()->SetTitle("Parameter index");
  tg_pars->GetYaxis()->SetTitle("Parameter");
  foutput->WriteTObject(tg_pars);
  delete tg_pars;

  foutput->WriteTObject(tg_all);
  delete tg_all;

  for (int bin=0; bin<nbinsx+1; bin++){
    if (bin<nbinsx) foutput->WriteTObject(hProb[bin]);
    else foutput->WriteTObject(hProb_all);
    foutput->WriteTObject(hEta[bin]);
    foutput->WriteTObject(hDist[bin]);
    delete hDist[bin];
    delete hEta[bin];
    if (bin<nbinsx) delete hProb[bin];
    else delete hProb_all;
  }
  delete tree;

  foutput->Close();
}

/*
gaus*(1+gaus(3))
->
Integrate[u*g[a,b,x]*(1 + v*h[c,d.x]), {x, 0, Infinity},Assumptions -> b > 0 && d > 0 && a > 0 && c > 0]
==
1/2 u (1 + Erf[a/(Sqrt[2] b)]) + (
E^(-((a - c)^2/(2 (b^2 + d^2))))
u v (1 + Erf[(b^2 c + a d^2)/(Sqrt[2] b d Sqrt[b^2 + d^2])]))/(
2 Sqrt[b^2 + d^2] Sqrt[2 \[Pi]])
*/
