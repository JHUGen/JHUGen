#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <memory>
#include "Mela.h"
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


using namespace RooFit;
using namespace std;

struct MELAwithOptions{
  Mela* mela;
  const double sqrts;
  const double mh;

  MELAwithOptions(double sqrts_, double mh_, TVar::VerbosityLevel verbosity_=TVar::ERROR) :
    sqrts(sqrts_), mh(mh_),
    mela(new Mela(sqrts_, mh_, verbosity_))
  { cout << "Called MELAwithOptions constructor!" << endl; }

  // FIXME: Not working yet, crashes after calling this destructor at delete mela
  ~MELAwithOptions(){ cout << "Called MELAwithOptions destructor!" << endl; delete mela; }
};

vector<unique_ptr<MELAwithOptions>> global_mela_list;

Mela* getMela(double sqrts, double mh, TVar::VerbosityLevel verbosity=TVar::ERROR){
  Mela* res=0;
  int it=0;
  for (auto& mwo : global_mela_list){
    if (mwo->sqrts==sqrts && mwo->mh==mh){
      res = mwo->mela;
      break;
    }
    it++;
  }
  if (res==0){
    unique_ptr<MELAwithOptions> tmp(new MELAwithOptions(sqrts, mh, verbosity));
    global_mela_list.push_back(std::move(tmp));
    res = global_mela_list.back()->mela;
  }
  res->setVerbosity(verbosity);
  return res;
}


shared_ptr<Mela> makemelaptr(int erg_tev, float mPOLE, TVar::VerbosityLevel verbosity){
  //function to make a shared_ptr in python, no idea how to do it directly
  return make_shared<Mela>(erg_tev, mPOLE, verbosity);
}


void testME_Dec_MCFM_Ping(int flavor=2, int useMothers=0, bool useConstants=false, shared_ptr<Mela> melaptr=nullptr){
  ofstream tout(TString("testME_Dec_MCFM_Ping_")+long(flavor)+"_"+long(useMothers)+"_"+long(useConstants)+".out");
  streambuf* coutbuf = cout.rdbuf();
  cout.rdbuf(tout.rdbuf());

  int erg_tev=13;
  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::DEBUG;
  if (!melaptr) {
    melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  }
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);

  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  //mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 79.9549392, 91.1876, 0.23119);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initializing" << endl;
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initialized" << endl;

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  /*
  const int nEntries = 3;
  double l1_array[nEntries][4] = {
  {1365.4973807340846,        10.289826593755228,        25.205694382277809,       -1365.2259480507332},
  {238.65751023078761,        9.2808858562825005,        15.827726043466324,       -237.95116187061188},
  {101.52463181523598,        27.359569630718468,      -0.90299073100241323,       -97.764458892691749}
  };
  double l2_array[nEntries][4] = {
  {22.786181013986834,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371},
  {101.67043553688544,        2.1169375132239789,       0.77953005873937187,       -101.64540506443268},
  {24.717634703436786,       -1.1722249478288802,       -5.9599387484197646,       -23.959684558009428}
  };
  double l3_array[nEntries][4] = {
  {1895.7562628816693,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620},
  {317.81904277258536,        2.5882005498984775,        21.352807448987718,       -317.09037005377883},
  {180.10885677707822,       -6.7240759244122792,        35.742176497019194,       -176.39865053838915}
  };
  double l4_array[nEntries][4] = {
  {471.71918486784784,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354},
  {95.655512770627581,       -13.986023919404957,       -37.960063551193414,       -86.679881365440792},
  {49.137252081251319,       -19.463268758477309,       -28.879247017597017,       -34.664676589120688}
  };
  */
  const int nEntries = 6;
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  for (int ev = 2; ev < 3; ev++){
    if (flavor == 2){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 1){
      GenLep1Id=11;
      GenLep2Id=-11;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 0){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 3){
      GenLep1Id=14;
      GenLep2Id=-14;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 4){
      GenLep1Id=0;
      GenLep2Id=-0;
      GenLep3Id=1;
      GenLep4Id=-1;
    }
    else if (flavor == 5){
      GenLep1Id=1;
      GenLep2Id=-1;
      GenLep3Id=2;
      GenLep4Id=-2;
    }

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    int idOrdered_WW[4] ={ 11, -12, -11, 12 };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters_WW;
    for (unsigned int idau=0; idau<4; idau++) daughters_WW.push_back(SimpleParticle_t(idOrdered_WW[idau], pOrdered[idau]));
    SimpleParticleCollection_t daughters_WWasZZ;
    for (unsigned int iv=0; iv<2; iv++){ for (int ivd=0; ivd<2; ivd++) daughters_WWasZZ.push_back(SimpleParticle_t(idOrdered_WW[iv+2*ivd], pOrdered[iv+2*ivd])); }
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    TLorentzVector pOrdered_ZG[3];
    pOrdered_ZG[0]=pOrdered[0];
    pOrdered_ZG[1]=pOrdered[1];
    pOrdered_ZG[2]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_ZG;
    for (unsigned int idau=0; idau<2; idau++) daughters_ZG.push_back(SimpleParticle_t(idOrdered[idau], pOrdered_ZG[idau]));
    for (unsigned int idau=2; idau<3; idau++) daughters_ZG.push_back(SimpleParticle_t(22, pOrdered_ZG[idau]));

    // Some gymnastics to get pseudo-mothers for these events
    TLorentzVector pTotal = pOrdered[0]+pOrdered[1]+pOrdered[2]+pOrdered[3];
    TLorentzVector pTotal_dummy=pTotal;
    TLorentzVector pTotal_perp(pTotal.X(), pTotal.Y(), 0, pTotal.T());
    pTotal_dummy.Boost(-pTotal_perp.BoostVector());
    TLorentzVector pMothers[2];
    pMothers[0].SetXYZT(0, 0, (pTotal_dummy.Z()+pTotal_dummy.T())/2., (pTotal_dummy.Z()+pTotal_dummy.T())/2.);
    pMothers[1].SetXYZT(0, 0, (pTotal_dummy.Z()-pTotal_dummy.T())/2., (-pTotal_dummy.Z()+pTotal_dummy.T())/2.);
    for (int im=0; im<2; im++) pMothers[im].Boost(pTotal_perp.BoostVector());
    SimpleParticleCollection_t mothers_QQB;
    for (unsigned int im=0; im<2; im++) mothers_QQB.push_back(SimpleParticle_t(1-2*im, pMothers[im]));
    SimpleParticleCollection_t mothers_GG;
    for (unsigned int im=0; im<2; im++) mothers_GG.push_back(SimpleParticle_t(21, pMothers[im]));
    SimpleParticleCollection_t* mothersPtr=0;
    if (useMothers==1) mothersPtr=&mothers_GG;
    else if (useMothers==2) mothersPtr=&mothers_QQB;

    TLorentzVector pOrdered_GG[3];
    pOrdered_GG[0]=pOrdered[0]+pOrdered[1];
    pOrdered_GG[1]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_GG;
    for (unsigned int idau=0; idau<2; idau++) daughters_GG.push_back(SimpleParticle_t(22, pOrdered_GG[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela.setInputEvent(&daughters_WW, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_WWasZZ, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));
    mela.setInputEvent(&daughters_ZG, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));
    mela.setCandidateDecayMode(TVar::CandidateDecay_GG);
    mela.setInputEvent(&daughters_GG, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));

    int cindex;

    /***** WW *****/
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    float pVAMCFM_qqWW_bkg;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(pVAMCFM_qqWW_bkg, useConstants);
    cout << "pVAMCFM_qqWW_bkg: " << pVAMCFM_qqWW_bkg << '\n' << endl;

    float pVAMCFM_ggVV_total;
    mela.setProcess(TVar::bkgWWZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_total, useConstants);
    cout << "pVAMCFM_ggVV_total: " << pVAMCFM_ggVV_total << '\n' << endl;
    float pVAMCFM_ggVV_bkg;
    mela.setProcess(TVar::bkgWWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_bkg, useConstants);
    cout << "pVAMCFM_ggVV_bkg: " << pVAMCFM_ggVV_bkg << '\n' << endl;
    float pVAMCFM_ggVV_sig;
    mela.setProcess(TVar::HSMHiggs_WWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_sig, useConstants);
    cout << "pVAMCFM_ggVV_sig: " << pVAMCFM_ggVV_sig << '\n' << endl;

    float pVAMCFM_ggWW_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_total, useConstants);
    cout << "pVAMCFM_ggWW_total: " << pVAMCFM_ggWW_total << '\n' << endl;
    float pVAMCFM_ggWW_bkg;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_bkg, useConstants);
    cout << "pVAMCFM_ggWW_bkg: " << pVAMCFM_ggWW_bkg << '\n' << endl;
    float pVAMCFM_ggWW_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_sig, useConstants);
    cout << "pVAMCFM_ggWW_sig: " << pVAMCFM_ggWW_sig << '\n' << endl;

    float pVAMCFM_ggWW_sig_selfDg1;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_selfDg1, useConstants);
    cout << "pVAMCFM_ggWW_sig_selfDg1: " << pVAMCFM_ggWW_sig_selfDg1 << '\n' << endl;

    float pVAMCFM_ggWW_sig_selfDg2;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_selfDg2, useConstants);
    cout << "pVAMCFM_ggWW_sig_selfDg2: " << pVAMCFM_ggWW_sig_selfDg2 << '\n' << endl;

    float pVAMCFM_ggWW_sig_selfDg4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_selfDg4, useConstants);
    cout << "pVAMCFM_ggWW_sig_selfDg4: " << pVAMCFM_ggWW_sig_selfDg4 << '\n' << endl;

    float pVAMCFM_ggWW_sig_selfDg1g2;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_selfDg1g2, useConstants);
    cout << "pVAMCFM_ggWW_sig_selfDg1g2: " << pVAMCFM_ggWW_sig_selfDg1g2 << '\n' << endl;

    float pVAMCFM_ggWW_sig_selfDg1g4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_selfDg1g4, useConstants);
    cout << "pVAMCFM_ggWW_sig_selfDg1g4: " << pVAMCFM_ggWW_sig_selfDg1g4 << '\n' << endl;

    float pVAMCFM_ggWW_sig_selfDg1g2im;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_selfDg1g2im, useConstants);
    cout << "pVAMCFM_ggWW_sig_selfDg1g2im: " << pVAMCFM_ggWW_sig_selfDg1g2im << '\n' << endl;

    float pVAMCFM_ggWW_sig_selfDg1g4im;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_selfDg1g4im, useConstants);
    cout << "pVAMCFM_ggWW_sig_selfDg1g4im: " << pVAMCFM_ggWW_sig_selfDg1g4im << '\n' << endl;

    float pVAMCFM_ggWW_sig_largemt4_kappat4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.resetMass(1e7, 8); // Already 1e5, but just to make sure
    mela.selfDHt4t4coupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_largemt4_kappat4, useConstants);
    cout << "pVAMCFM_ggWW_sig_largemt4_kappat4: " << pVAMCFM_ggWW_sig_largemt4_kappat4 << '\n' << endl;

    float pVAMCFM_ggWW_sig_largemt4_kappatildet4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.resetMass(1e7, 8); // Already 1e5, but just to make sure
    mela.selfDHt4t4coupl[0][gHIGGS_KAPPA_TILDE][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_largemt4_kappatildet4, useConstants);
    cout << "pVAMCFM_ggWW_sig_largemt4_kappatildet4: " << pVAMCFM_ggWW_sig_largemt4_kappatildet4 << '\n' << endl;

    float pVAMCFM_ggWW_sig_ghg2_gen4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHg4g4coupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_ghg2_gen4, useConstants);
    cout << "pVAMCFM_ggWW_sig_ghg2_gen4: " << pVAMCFM_ggWW_sig_ghg2_gen4 << '\n' << endl;

    float pVAMCFM_ggWW_sig_ghg4_gen4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHg4g4coupl[0][gHIGGS_GG_4][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_ghg4_gen4, useConstants);
    cout << "pVAMCFM_ggWW_sig_ghg4_gen4: " << pVAMCFM_ggWW_sig_ghg4_gen4 << '\n' << endl;

    float pVAMCFM_ggWW_sig_ghg2ghg4_gen4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHg4g4coupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHg4g4coupl[0][gHIGGS_GG_4][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_ghg2ghg4_gen4, useConstants);
    cout << "pVAMCFM_ggWW_sig_ghg2ghg4_gen4: " << pVAMCFM_ggWW_sig_ghg2ghg4_gen4 << '\n' << endl;

    float pVAMCFM_ggWW_sig_ghg2ghg4im_gen4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHg4g4coupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHg4g4coupl[0][gHIGGS_GG_4][1]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWW_sig_ghg2ghg4im_gen4, useConstants);
    cout << "pVAMCFM_ggWW_sig_ghg2ghg4im_gen4: " << pVAMCFM_ggWW_sig_ghg2ghg4im_gen4 << '\n' << endl;


    float pVAJHU_ggWW_sig_selfDg1;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_selfDg1, useConstants);
    cout << "pVAJHU_ggWW_sig_selfDg1: " << pVAJHU_ggWW_sig_selfDg1 << '\n' << endl;

    float pVAJHU_ggWW_sig_selfDg2;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_selfDg2, useConstants);
    cout << "pVAJHU_ggWW_sig_selfDg2: " << pVAJHU_ggWW_sig_selfDg2 << '\n' << endl;

    float pVAJHU_ggWW_sig_selfDg4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_selfDg4, useConstants);
    cout << "pVAJHU_ggWW_sig_selfDg4: " << pVAJHU_ggWW_sig_selfDg4 << '\n' << endl;

    float pVAJHU_ggWW_sig_selfDg1g2;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_selfDg1g2, useConstants);
    cout << "pVAJHU_ggWW_sig_selfDg1g2: " << pVAJHU_ggWW_sig_selfDg1g2 << '\n' << endl;

    float pVAJHU_ggWW_sig_selfDg1g4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_selfDg1g4, useConstants);
    cout << "pVAJHU_ggWW_sig_selfDg1g4: " << pVAJHU_ggWW_sig_selfDg1g4 << '\n' << endl;

    float pVAJHU_ggWW_sig_selfDg1g2im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_selfDg1g2im, useConstants);
    cout << "pVAJHU_ggWW_sig_selfDg1g2im: " << pVAJHU_ggWW_sig_selfDg1g2im << '\n' << endl;

    float pVAJHU_ggWW_sig_selfDg1g4im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_selfDg1g4im, useConstants);
    cout << "pVAJHU_ggWW_sig_selfDg1g4im: " << pVAJHU_ggWW_sig_selfDg1g4im << '\n' << endl;

    float pVAJHU_ggWW_sig_ghg4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_4][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_ghg4, useConstants);
    cout << "pVAJHU_ggWW_sig_ghg4: " << pVAJHU_ggWW_sig_ghg4 << '\n' << endl;

    float pVAJHU_ggWW_sig_ghg2ghg4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHggcoupl[0][gHIGGS_GG_4][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_ghg2ghg4, useConstants);
    cout << "pVAJHU_ggWW_sig_ghg2ghg4: " << pVAJHU_ggWW_sig_ghg2ghg4 << '\n' << endl;

    float pVAJHU_ggWW_sig_ghg2ghg4im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHggcoupl[0][gHIGGS_GG_4][1]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWW_sig_ghg2ghg4im, useConstants);
    cout << "pVAJHU_ggWW_sig_ghg2ghg4im: " << pVAJHU_ggWW_sig_ghg2ghg4im << '\n' << endl;

    cout << "MCFM vs JHUGen Ratio comparison:" << endl;
    cout << "ggWW_sig_selfDg2: " << pVAMCFM_ggWW_sig_selfDg2/pVAMCFM_ggWW_sig_selfDg1 << '\t' << pVAJHU_ggWW_sig_selfDg2/pVAJHU_ggWW_sig_selfDg1 << endl;
    cout << "ggWW_sig_selfDg1g2: " << pVAMCFM_ggWW_sig_selfDg1g2/pVAMCFM_ggWW_sig_selfDg1 << '\t' << pVAJHU_ggWW_sig_selfDg1g2/pVAJHU_ggWW_sig_selfDg1 << endl;
    cout << "ggWW_sig_selfDg1g2im: " << pVAMCFM_ggWW_sig_selfDg1g2im/pVAMCFM_ggWW_sig_selfDg1 << '\t' << pVAJHU_ggWW_sig_selfDg1g2im/pVAJHU_ggWW_sig_selfDg1 << endl;
    cout << "ggWW_sig_selfDg4: " << pVAMCFM_ggWW_sig_selfDg4/pVAMCFM_ggWW_sig_selfDg1 << '\t' << pVAJHU_ggWW_sig_selfDg4/pVAJHU_ggWW_sig_selfDg1 << endl;
    cout << "ggWW_sig_selfDg1g4: " << pVAMCFM_ggWW_sig_selfDg1g4/pVAMCFM_ggWW_sig_selfDg1 << '\t' << pVAJHU_ggWW_sig_selfDg1g4/pVAJHU_ggWW_sig_selfDg1 << endl;
    cout << "ggWW_sig_selfDg1g4im: " << pVAMCFM_ggWW_sig_selfDg1g4im/pVAMCFM_ggWW_sig_selfDg1 << '\t' << pVAJHU_ggWW_sig_selfDg1g4im/pVAJHU_ggWW_sig_selfDg1 << endl;
    cout << "***" << endl;
    cout << "MCFM_ggWW_sig_largemt4_kappatildet4: " << pVAMCFM_ggWW_sig_largemt4_kappatildet4/pVAMCFM_ggWW_sig_largemt4_kappat4 << endl;
    cout << "MCFM_ggWW_sig_ghg2_gen4: " << pVAMCFM_ggWW_sig_ghg2_gen4/pVAMCFM_ggWW_sig_largemt4_kappat4 << endl;
    cout << "ggWW_sig_ghg4_gen4: " << pVAMCFM_ggWW_sig_ghg4_gen4/pVAMCFM_ggWW_sig_ghg2_gen4 << '\t' << pVAJHU_ggWW_sig_ghg4/pVAJHU_ggWW_sig_selfDg1 << endl;
    cout << "ggWW_sig_ghg2ghg4_gen4: " << pVAMCFM_ggWW_sig_ghg2ghg4_gen4/pVAMCFM_ggWW_sig_ghg2_gen4 << '\t' << pVAJHU_ggWW_sig_ghg2ghg4/pVAJHU_ggWW_sig_selfDg1 << endl;
    cout << "ggWW_sig_ghg2ghg4im_gen4: " << pVAMCFM_ggWW_sig_ghg2ghg4im_gen4/pVAMCFM_ggWW_sig_ghg2_gen4 << '\t' << pVAJHU_ggWW_sig_ghg2ghg4im/pVAJHU_ggWW_sig_selfDg1 << endl;
    cout << endl;

    /***** WW (as ZZ) *****/
    cindex=1;
    mela.setCurrentCandidateFromIndex(cindex);

    float pVAMCFM_ggVV_fromZZ_total;
    mela.setProcess(TVar::bkgWWZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_fromZZ_total, useConstants);
    cout << "pVAMCFM_ggVV_fromZZ_total: " << pVAMCFM_ggVV_fromZZ_total << '\n' << endl;
    float pVAMCFM_ggVV_fromZZ_bkg;
    mela.setProcess(TVar::bkgWWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_fromZZ_bkg, useConstants);
    cout << "pVAMCFM_ggVV_fromZZ_bkg: " << pVAMCFM_ggVV_fromZZ_bkg << '\n' << endl;
    float pVAMCFM_ggVV_fromZZ_sig;
    mela.setProcess(TVar::HSMHiggs_WWZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggVV_fromZZ_sig, useConstants);
    cout << "pVAMCFM_ggVV_fromZZ_sig: " << pVAMCFM_ggVV_fromZZ_sig << '\n' << endl;

    float pVAMCFM_ggZZ_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_total, useConstants);
    cout << "pVAMCFM_ggZZ_total from WW as ZZ: " << pVAMCFM_ggZZ_total << '\n' << endl;
    float pVAMCFM_ggZZ_bkg;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_bkg, useConstants);
    cout << "pVAMCFM_ggZZ_bkg from WW as ZZ: " << pVAMCFM_ggZZ_bkg << '\n' << endl;
    float pVAMCFM_ggWWasZZ_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWWasZZ_sig, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig: " << pVAMCFM_ggWWasZZ_sig << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1: " << pVAMCFM_ggWWasZZ_sig_selfDg1 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg2;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg2, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg2: " << pVAMCFM_ggWWasZZ_sig_selfDg2 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg4: " << pVAMCFM_ggWWasZZ_sig_selfDg4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDgzgs2;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDgzgs2, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDgzgs2: " << pVAMCFM_ggWWasZZ_sig_selfDgzgs2 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDgzgs4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDgzgs4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDgzgs4: " << pVAMCFM_ggWWasZZ_sig_selfDgzgs4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDggsgs2;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDggsgs2, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDggsgs2: " << pVAMCFM_ggWWasZZ_sig_selfDggsgs2 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDggsgs4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDggsgs4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDggsgs4: " << pVAMCFM_ggWWasZZ_sig_selfDggsgs4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1g2;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1g2, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1g2: " << pVAMCFM_ggWWasZZ_sig_selfDg1g2 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1g4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1g4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1g4: " << pVAMCFM_ggWWasZZ_sig_selfDg1g4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2: " << pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4: " << pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2: " << pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4: " << pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1g2im;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1g2im, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1g2im: " << pVAMCFM_ggWWasZZ_sig_selfDg1g2im << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1g4im;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1g4im, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1g4im: " << pVAMCFM_ggWWasZZ_sig_selfDg1g4im << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2im;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_2][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2im, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2im: " << pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2im << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4im;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_4][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4im, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4im: " << pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4im << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2im;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_2][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2im, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2im: " << pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2im << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4im;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_4][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4im, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4im: " << pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4im << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_largemt4_kappat4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.resetMass(1e7, 8); // Already 1e5, but just to make sure
    mela.selfDHt4t4coupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_largemt4_kappat4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_largemt4_kappat4: " << pVAMCFM_ggWWasZZ_sig_largemt4_kappat4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_largemt4_kappatildet4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.resetMass(1e7, 8); // Already 1e5, but just to make sure
    mela.selfDHt4t4coupl[0][gHIGGS_KAPPA_TILDE][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_largemt4_kappatildet4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_largemt4_kappatildet4: " << pVAMCFM_ggWWasZZ_sig_largemt4_kappatildet4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_ghg2_gen4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHg4g4coupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_ghg2_gen4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_ghg2_gen4: " << pVAMCFM_ggWWasZZ_sig_ghg2_gen4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_ghg4_gen4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHg4g4coupl[0][gHIGGS_GG_4][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_ghg4_gen4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_ghg4_gen4: " << pVAMCFM_ggWWasZZ_sig_ghg4_gen4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_ghg2ghg4_gen4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHg4g4coupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHg4g4coupl[0][gHIGGS_GG_4][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_ghg2ghg4_gen4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_ghg2ghg4_gen4: " << pVAMCFM_ggWWasZZ_sig_ghg2ghg4_gen4 << '\n' << endl;

    float pVAMCFM_ggWWasZZ_sig_ghg2ghg4im_gen4;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHg4g4coupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHg4g4coupl[0][gHIGGS_GG_4][1]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggWWasZZ_sig_ghg2ghg4im_gen4, useConstants);
    cout << "pVAMCFM_ggWWasZZ_sig_ghg2ghg4im_gen4: " << pVAMCFM_ggWWasZZ_sig_ghg2ghg4im_gen4 << '\n' << endl;


    float pVAJHU_ggWWasZZ_sig_selfDg1;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1: " << pVAJHU_ggWWasZZ_sig_selfDg1 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg2;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg2, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg2: " << pVAJHU_ggWWasZZ_sig_selfDg2 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg4, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg4: " << pVAJHU_ggWWasZZ_sig_selfDg4 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDgzgs2;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDgzgs2, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDgzgs2: " << pVAJHU_ggWWasZZ_sig_selfDgzgs2 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDgzgs4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDgzgs4, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDgzgs4: " << pVAJHU_ggWWasZZ_sig_selfDgzgs4 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDggsgs2;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDggsgs2, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDggsgs2: " << pVAJHU_ggWWasZZ_sig_selfDggsgs2 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDggsgs4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDggsgs4, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDggsgs4: " << pVAJHU_ggWWasZZ_sig_selfDggsgs4 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1g2;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1g2, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1g2: " << pVAJHU_ggWWasZZ_sig_selfDg1g2 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1g4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1g4, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1g4: " << pVAJHU_ggWWasZZ_sig_selfDg1g4 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1gzgs2;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1gzgs2, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1gzgs2: " << pVAJHU_ggWWasZZ_sig_selfDg1gzgs2 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1gzgs4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1gzgs4, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1gzgs4: " << pVAJHU_ggWWasZZ_sig_selfDg1gzgs4 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_2][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2: " << pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4: " << pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1g2im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1g2im, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1g2im: " << pVAJHU_ggWWasZZ_sig_selfDg1g2im << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1g4im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1g4im, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1g4im: " << pVAJHU_ggWWasZZ_sig_selfDg1g4im << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1gzgs2im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_2][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1gzgs2im, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1gzgs2im: " << pVAJHU_ggWWasZZ_sig_selfDg1gzgs2im << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1gzgs4im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_4][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1gzgs4im, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1gzgs4im: " << pVAJHU_ggWWasZZ_sig_selfDg1gzgs4im << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_2][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2im, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2im: " << pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2im << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_AA_4][1]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4im, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4im: " << pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4im << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_ghg4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_4][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_ghg4, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_ghg4: " << pVAJHU_ggWWasZZ_sig_ghg4 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_ghg2ghg4;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHggcoupl[0][gHIGGS_GG_4][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_ghg2ghg4, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_ghg2ghg4: " << pVAJHU_ggWWasZZ_sig_ghg2ghg4 << '\n' << endl;

    float pVAJHU_ggWWasZZ_sig_ghg2ghg4im;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.0;
    mela.selfDHggcoupl[0][gHIGGS_GG_4][1]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAJHU_ggWWasZZ_sig_ghg2ghg4im, useConstants);
    cout << "pVAJHU_ggWWasZZ_sig_ghg2ghg4im: " << pVAJHU_ggWWasZZ_sig_ghg2ghg4im << '\n' << endl;


    cout << "MCFM vs JHUGen Ratio comparison:" << endl;
    cout << "ggWWasZZ_sig_selfDg2: " << pVAMCFM_ggWWasZZ_sig_selfDg2/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg2/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1g2: " << pVAMCFM_ggWWasZZ_sig_selfDg1g2/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1g2/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1g2im: " << pVAMCFM_ggWWasZZ_sig_selfDg1g2im/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1g2im/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg4: " << pVAMCFM_ggWWasZZ_sig_selfDg4/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg4/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1g4: " << pVAMCFM_ggWWasZZ_sig_selfDg1g4/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1g4/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1g4im: " << pVAMCFM_ggWWasZZ_sig_selfDg1g4im/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1g4im/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "***" << endl;
    cout << "ggWWasZZ_sig_selfDgzgs2: " << pVAMCFM_ggWWasZZ_sig_selfDgzgs2/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDgzgs2/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1gzgs2: " << pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1gzgs2/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1gzgs2im: " << pVAMCFM_ggWWasZZ_sig_selfDg1gzgs2im/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1gzgs2im/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDgzgs4: " << pVAMCFM_ggWWasZZ_sig_selfDgzgs4/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDgzgs4/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1gzgs4: " << pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1gzgs4/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1gzgs4im: " << pVAMCFM_ggWWasZZ_sig_selfDg1gzgs4im/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1gzgs4im/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "***" << endl;
    cout << "ggWWasZZ_sig_selfDggsgs2: " << pVAMCFM_ggWWasZZ_sig_selfDggsgs2/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDggsgs2/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1ggsgs2: " << pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1ggsgs2im: " << pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs2im/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1ggsgs2im/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDggsgs4: " << pVAMCFM_ggWWasZZ_sig_selfDggsgs4/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDggsgs4/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1ggsgs4: " << pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_selfDg1ggsgs4im: " << pVAMCFM_ggWWasZZ_sig_selfDg1ggsgs4im/pVAMCFM_ggWWasZZ_sig_selfDg1 << '\t' << pVAJHU_ggWWasZZ_sig_selfDg1ggsgs4im/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "***" << endl;
    cout << "MCFM_ggWWasZZ_sig_largemt4_kappatildet4: " << pVAMCFM_ggWWasZZ_sig_largemt4_kappatildet4/pVAMCFM_ggWWasZZ_sig_largemt4_kappat4 << endl;
    cout << "MCFM_ggWWasZZ_sig_ghg2_gen4: " << pVAMCFM_ggWWasZZ_sig_ghg2_gen4/pVAMCFM_ggWWasZZ_sig_largemt4_kappat4 << endl;
    cout << "ggWWasZZ_sig_ghg4_gen4: " << pVAMCFM_ggWWasZZ_sig_ghg4_gen4/pVAMCFM_ggWWasZZ_sig_ghg2_gen4 << '\t' << pVAJHU_ggWWasZZ_sig_ghg4/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_ghg2ghg4_gen4: " << pVAMCFM_ggWWasZZ_sig_ghg2ghg4_gen4/pVAMCFM_ggWWasZZ_sig_ghg2_gen4 << '\t' << pVAJHU_ggWWasZZ_sig_ghg2ghg4/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << "ggWWasZZ_sig_ghg2ghg4im_gen4: " << pVAMCFM_ggWWasZZ_sig_ghg2ghg4im_gen4/pVAMCFM_ggWWasZZ_sig_ghg2_gen4 << '\t' << pVAJHU_ggWWasZZ_sig_ghg2ghg4im/pVAJHU_ggWWasZZ_sig_selfDg1 << endl;
    cout << endl;

    /***** ZZ *****/
    cindex=2;
    mela.setCurrentCandidateFromIndex(cindex);

    float pVAMCFM_qqZZ_bkg;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(pVAMCFM_qqZZ_bkg, useConstants);
    cout << "pVAMCFM_qqZZ_bkg: " << pVAMCFM_qqZZ_bkg << '\n' << endl;

    float pVAMCFM_qqZJJ_bkg;
    mela.setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);
    mela.computeP(pVAMCFM_qqZJJ_bkg, useConstants);
    cout << "pVAMCFM_qqZJJ_bkg: " << pVAMCFM_qqZJJ_bkg << '\n' << endl;

    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_total, useConstants);
    cout << "pVAMCFM_ggZZ_total: " << pVAMCFM_ggZZ_total << '\n' << endl;

    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_bkg, useConstants);
    cout << "pVAMCFM_ggZZ_bkg: " << pVAMCFM_ggZZ_bkg << '\n' << endl;

    float pVAMCFM_ggZZ_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ggZZ_sig, useConstants);
    cout << "pVAMCFM_ggZZ_sig: " << pVAMCFM_ggZZ_sig << '\n' << endl;

    float pVAMCFM_ZZ_sig_selfDg1;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.0;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.0;
    mela.setMelaHiggsWidth(wPOLE);
    mela.setMelaLeptonInterference(TVar::InterfOn);
    mela.computeP(pVAMCFM_ZZ_sig_selfDg1, useConstants);
    cout << "pVAMCFM_ggZZ_sig_selfDg1: " << pVAMCFM_ZZ_sig_selfDg1 << '\n' << endl;

    /***** ZG *****/
    cindex=3;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAMCFM_ZG_bkg;
    mela.setProcess(TVar::bkgZGamma, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(pVAMCFM_ZG_bkg, useConstants);
    cout << "pVAMCFM_qqZG_bkg: " << pVAMCFM_ZG_bkg << '\n' << endl;

    if (verbosity>=TVar::DEBUG){
      cout << "Removing Mela candidates\nSummary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }
    mela.resetInputEvent();
    cout << "Removed..." << endl;
  }

  cout.rdbuf(coutbuf);
  tout.close();
  mela.setVerbosity(bkpverbosity);
}

void testME_VH_JHUGen_Ping(int erg_tev=13, bool useConstants=false, shared_ptr<Mela> melaptr=nullptr){
  TString strtout = Form("testME_VH_JHUGen_Ping_%iTeV_%i.out", erg_tev, (int)useConstants);
  ofstream tout(strtout.Data());
  streambuf* coutbuf = cout.rdbuf();
  cout.rdbuf(tout.rdbuf());

  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  if (verbosity>=TVar::DEBUG) cout << "Initializing Mela..." << endl;
  if (!melaptr) melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  //mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 79.9549392, 91.1876, 0.23119);

  int GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  const int nEntries = 6;
  double a1_array[nEntries][4] ={
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 },
    { 238.65751023078761, 9.2808858562825005, 15.827726043466324, -237.95116187061188 },
    { 101.52463181523598, 27.359569630718468, -0.90299073100241323, -97.764458892691749 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 101.67043553688544, 2.1169375132239789, 0.77953005873937187, -101.64540506443268 },
    { 24.717634703436786, -1.1722249478288802, -5.9599387484197646, -23.959684558009428 }
  };
  double a2_array[nEntries][4] ={
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 317.81904277258536, 2.5882005498984775, 21.352807448987718, -317.09037005377883 },
    { 180.10885677707822, -6.7240759244122792, 35.742176497019194, -176.39865053838915 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 95.655512770627581, -13.986023919404957, -37.960063551193414, -86.679881365440792 },
    { 49.137252081251319, -19.463268758477309, -28.879247017597017, -34.664676589120688 }
  };
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  // Decay mode does not matter, just use something
  GenLep1Id=13;
  GenLep2Id=-13;
  GenLep3Id=11;
  GenLep4Id=-11;

  for (int ev = 2; ev < 3; ev++){
    SimpleParticleCollection_t aparticles;
    TLorentzVector pAPart[2];
    pAPart[0].SetXYZT(a1_array[ev][1], a1_array[ev][2], a1_array[ev][3], a1_array[ev][0]);
    pAPart[1].SetXYZT(a2_array[ev][1], a2_array[ev][2], a2_array[ev][3], a2_array[ev][0]);
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t(0, pAPart[iap])); // q q'
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t((1-2*iap)*13, pAPart[iap])); // l- l+
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t((1-2*iap)*14, pAPart[iap])); // nub nu
    for (unsigned int iap=0; iap<1; iap++) aparticles.push_back(SimpleParticle_t(22, pAPart[iap])); // gamma

    int idOrdered[4] ={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters;
    for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters, &aparticles, (SimpleParticleCollection_t*)0, false);

    cout << "*******************************************************" << endl;
    for (int ic=0; ic<mela.getNCandidates(); ic++){
      cout << "Summary of candidate " << ic << ":" << endl;
      mela.setCurrentCandidateFromIndex(ic);
      TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      cout << "*******************************************************" << endl;
    }
    cout << "*******************************************************" << endl;
    cout << endl;

    int cindex;
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    TString strlh[3]={ "leptonic", "hadronic", "photonic" };
    TString strzwh[3]={ "ZH", "WH", "GammaH" };
    TVar::Production prod;
    for (unsigned int ilh=0; ilh<3; ilh++){
      TString lhapp = strlh[ilh];
      for (unsigned int izwa=0; izwa<3; izwa++){
        TString zwapp = strzwh[izwa];
        if (ilh==0 && izwa==0) prod=TVar::Lep_ZH;
        else if (ilh==0 && izwa==1) prod=TVar::Lep_WH;
        else if (ilh==1 && izwa==0) prod=TVar::Had_ZH;
        else if (ilh==1 && izwa==1) prod=TVar::Had_WH;
        else if (ilh==2 && izwa==2) prod=TVar::GammaH;
        else continue;

        cout << "*******************************************************" << endl;
        cout << "Computing MEs for " << strlh[ilh] << " " << strzwh[izwa] << endl;

        float p0mplus=0;
        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, prod);
        mela.computeProdP_VH(p0mplus, false, useConstants);
        cout << "p0mplus: " << p0mplus << '\n' << endl;

        float p0g1prime2=0;
        mela.setProcess(TVar::H0_g1prime2, TVar::JHUGen, prod);
        mela.computeProdP_VH(p0g1prime2, false, useConstants);
        cout << "p0g1prime2: " << p0g1prime2 << '\n' << endl;

        float p0hplus=0;
        mela.setProcess(TVar::H0hplus, TVar::JHUGen, prod);
        mela.computeProdP_VH(p0hplus, false, useConstants);
        cout << "p0hplus: " << p0hplus << '\n' << endl;

        float p0minus=0;
        mela.setProcess(TVar::H0minus, TVar::JHUGen, prod);
        mela.computeProdP_VH(p0minus, false, useConstants);
        cout << "p0minus: " << p0minus << '\n' << endl;

        float p0gzgs1prime2=0;
        mela.setProcess(TVar::H0_Zgsg1prime2, TVar::JHUGen, prod);
        mela.computeProdP_VH(p0gzgs1prime2, false, useConstants);
        cout << "p0gzgs1prime2: " << p0gzgs1prime2 << '\n' << endl;

        float p0hpluszgs=0;
        mela.setProcess(TVar::H0_Zgs, TVar::JHUGen, prod);
        mela.computeProdP_VH(p0hpluszgs, false, useConstants);
        cout << "p0hpluszgs: " << p0hpluszgs << '\n' << endl;

        // SelfD MEs
        float p0mplus_selfD=0;
        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
        mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        mela.computeProdP_VH(p0mplus_selfD, false, useConstants);
        cout << "p0mplus_selfD: " << p0mplus_selfD << '\n' << endl;

        if (izwa!=2){
          float p0mplus_selfD_CT=0;
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
          double xw = 0.23119;
          if (izwa==0){
            mela.selfDHzpzpcoupl[gHIGGS_VV_1][0]=1;
            mela.selfDZpffcoupl[gHIGGS_Vp_NuE_left][0]=0.5*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_El_left][0]=(-xw*(-1.)-0.5)*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Mu_left][0]=(-xw*(-1.)-0.5)*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Ta_left][0]=(-xw*(-1.)-0.5)*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Dn_left][0]=(-xw*(-1./3.)-0.5)*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Up_left][0]=(-xw*(2./3.)+0.5)*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Str_left][0]=(-xw*(-1./3.)-0.5)*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Chm_left][0]=(-xw*(2./3.)+0.5)*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Bot_left][0]=(-xw*(-1./3.)-0.5)*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_El_right][0]=(-xw*(-1.))*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Mu_right][0]=(-xw*(-1.))*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Ta_right][0]=(-xw*(-1.))*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Dn_right][0]=(-xw*(-1./3.))*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Up_right][0]=(-xw*(2./3.))*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Str_right][0]=(-xw*(-1./3.))*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Chm_right][0]=(-xw*(2./3.))*2.;
            mela.selfDZpffcoupl[gHIGGS_Vp_Bot_right][0]=(-xw*(-1./3.))*2.;
            mela.selfDM_Zprime=91.1876;
            mela.selfDGa_Zprime=2.4952;
          }
          else{
            mela.selfDHzpzpcoupl[gHIGGS_VV_1][0]=1;
            //mela.selfDHwpwpcoupl[gHIGGS_VV_1][0]=1;
            mela.selfDWpffcoupl[gHIGGS_Vp_El_left][0]=sqrt(2.*(1.-xw));
            mela.selfDWpffcoupl[gHIGGS_Vp_Mu_left][0]=sqrt(2.*(1.-xw));
            mela.selfDWpffcoupl[gHIGGS_Vp_Ta_left][0]=sqrt(2.*(1.-xw));
            mela.selfDWpffcoupl[gHIGGS_Vp_Up_left][0]=sqrt(2.*(1.-xw));
            mela.selfDWpffcoupl[gHIGGS_Vp_Chm_left][0]=sqrt(2.*(1.-xw));
            mela.selfDWpffcoupl[gHIGGS_Vp_Top_left][0]=sqrt(2.*(1.-xw));
            mela.selfDM_Wprime=80.399;
            mela.selfDGa_Wprime=2.085;
          }
          mela.computeProdP_VH(p0mplus_selfD_CT, false, useConstants);
          cout << "p0mplus_selfD_CT: " << p0mplus_selfD_CT << '\n' << endl;
        }

        float p0g1prime2_selfD=0;
        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
        mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=1;
        mela.computeProdP_VH(p0g1prime2_selfD, false, useConstants);
        cout << "p0g1prime2_selfD: " << p0g1prime2_selfD << '\n' << endl;

        float p0hplus_selfD=0;
        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
        mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1;
        mela.computeProdP_VH(p0hplus_selfD, false, useConstants);
        cout << "p0hplus_selfD: " << p0hplus_selfD << '\n' << endl;

        float p0minus_selfD=0;
        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
        mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
        mela.computeProdP_VH(p0minus_selfD, false, useConstants);
        cout << "p0minus_selfD: " << p0minus_selfD << '\n' << endl;

        float p0gzgs1prime2_selfD=0;
        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
        mela.selfDHzzcoupl[0][gHIGGS_ZA_1_PRIME2][0]=1;
        mela.computeProdP_VH(p0gzgs1prime2_selfD, false, useConstants);
        cout << "p0gzgs1prime2_selfD: " << p0gzgs1prime2_selfD << '\n' << endl;

        float p0hpluszgs_selfD=0;
        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
        mela.selfDHzzcoupl[0][gHIGGS_ZA_2][0]=1;
        mela.computeProdP_VH(p0hpluszgs_selfD, false, useConstants);
        cout << "p0hpluszgs_selfD: " << p0hpluszgs_selfD << '\n' << endl;

        float recoBW=-99;
        mela.computeDijetConvBW(recoBW);
        cout << "Reco BW: " << recoBW << '\n' << endl;

        float idealBW=-99;
        mela.computeDijetConvBW(idealBW, true);
        cout << "Ideal BW: " << idealBW << '\n' << endl;

        cout << "*******************************************************" << endl;

        float costhetastar = 0, costheta1 = 0, costheta2 = 0, Phi = 0, Phi1 = 0;
        if (prod != TVar::GammaH) mela.computeVHAngles(
          costheta1,
          costheta2,
          Phi,
          costhetastar,
          Phi1
        );
        cout << "VH (" << TVar::ProductionName(prod) << ") angles: " << costheta1 << " " << costheta2 << " " << Phi << " " << costhetastar << " " << Phi1 << endl;

      }
    }

    mela.resetInputEvent();
  }

  cout.rdbuf(coutbuf);
  tout.close();
  mela.setVerbosity(bkpverbosity);
}

void testME_Dec_JHUGen_Ping(int erg_tev=13, bool useConstants=false, shared_ptr<Mela> melaptr=nullptr){
  TString strtout = Form("testME_Dec_JHUGen_Ping_%iTeV_%i.out", erg_tev, (int)useConstants);
  ofstream tout(strtout.Data());
  streambuf* coutbuf = cout.rdbuf();
  cout.rdbuf(tout.rdbuf());

  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  if (verbosity>=TVar::DEBUG) cout << "Initializing Mela..." << endl;
  if (!melaptr) melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  auto bkpprecision = cout.precision(10);
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  //mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 79.9549392, 91.1876, 0.23119);

  int GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  const int nEntries = 6;
  double a1_array[nEntries][4] ={
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 },
    { 238.65751023078761, 9.2808858562825005, 15.827726043466324, -237.95116187061188 },
    { 101.52463181523598, 27.359569630718468, -0.90299073100241323, -97.764458892691749 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 101.67043553688544, 2.1169375132239789, 0.77953005873937187, -101.64540506443268 },
    { 24.717634703436786, -1.1722249478288802, -5.9599387484197646, -23.959684558009428 }
  };
  double a2_array[nEntries][4] ={
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 317.81904277258536, 2.5882005498984775, 21.352807448987718, -317.09037005377883 },
    { 180.10885677707822, -6.7240759244122792, 35.742176497019194, -176.39865053838915 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 95.655512770627581, -13.986023919404957, -37.960063551193414, -86.679881365440792 },
    { 49.137252081251319, -19.463268758477309, -28.879247017597017, -34.664676589120688 }
  };
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  // Decay mode does not matter, just use something
  GenLep1Id=13;
  GenLep2Id=-13;
  GenLep3Id=11;
  GenLep4Id=-11;

  for (int ev = 2; ev < 3; ev++){
    SimpleParticleCollection_t aparticles;
    TLorentzVector pAPart[2];
    pAPart[0].SetXYZT(a1_array[ev][1], a1_array[ev][2], a1_array[ev][3], a1_array[ev][0]);
    pAPart[1].SetXYZT(a2_array[ev][1], a2_array[ev][2], a2_array[ev][3], a2_array[ev][0]);
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t(0, pAPart[iap])); // q q'
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t((1-2*iap)*13, pAPart[iap])); // l- l+
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t((1-2*iap)*14, pAPart[iap])); // nub nu
    for (unsigned int iap=0; iap<1; iap++) aparticles.push_back(SimpleParticle_t(22, pAPart[iap])); // gamma

    int idOrdered[4] ={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters;
    for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters, &aparticles, (SimpleParticleCollection_t*)0, false);

    cout << "*******************************************************" << endl;
    for (int ic=0; ic<mela.getNCandidates(); ic++){
      cout << "Summary of candidate " << ic << ":" << endl;
      mela.setCurrentCandidateFromIndex(ic);
      TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      cout << "*******************************************************" << endl;
    }
    cout << "*******************************************************" << endl;
    cout << endl;

    int cindex;
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    TString strlh[3]={ "leptonic", "hadronic", "photonic" };
    TString strzwh[3]={ "ZH", "WH", "GammaH" };
    TVar::Production prod;

    cout << "*******************************************************" << endl;
    cout << "Computing MEs for ZZ decay" << endl;

    float p0mplus=0;
    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(p0mplus, useConstants);
    cout << "p0mplus: " << p0mplus << '\n' << endl;

    float p0g1prime2=0;
    mela.setProcess(TVar::H0_g1prime2, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(p0g1prime2, useConstants);
    cout << "p0g1prime2: " << p0g1prime2 << '\n' << endl;

    float p0hplus=0;
    mela.setProcess(TVar::H0hplus, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(p0hplus, useConstants);
    cout << "p0hplus: " << p0hplus << '\n' << endl;

    float p0minus=0;
    mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(p0minus, useConstants);
    cout << "p0minus: " << p0minus << '\n' << endl;

    float p0gzgs1prime2=0;
    mela.setProcess(TVar::H0_Zgsg1prime2, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(p0gzgs1prime2, useConstants);
    cout << "p0gzgs1prime2: " << p0gzgs1prime2 << '\n' << endl;

    float p0hpluszgs=0;
    mela.setProcess(TVar::H0_Zgs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(p0hpluszgs, useConstants);
    cout << "p0hpluszgs: " << p0hpluszgs << '\n' << endl;

    // SelfD MEs
    float p0mplus_selfD=0;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p0mplus_selfD, useConstants);
    cout << "p0mplus_selfD: " << p0mplus_selfD << '\n' << endl;

    float p0mplus_selfD_CT=0;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    double xw = 0.23119;
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzpzpcoupl[gHIGGS_VV_1][0]=1;
    mela.selfDZpffcoupl[gHIGGS_Vp_NuE_left][0]=0.5*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_El_left][0]=(-xw*(-1.)-0.5)*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Mu_left][0]=(-xw*(-1.)-0.5)*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Ta_left][0]=(-xw*(-1.)-0.5)*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Dn_left][0]=(-xw*(-1./3.)-0.5)*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Up_left][0]=(-xw*(2./3.)+0.5)*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Str_left][0]=(-xw*(-1./3.)-0.5)*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Chm_left][0]=(-xw*(2./3.)+0.5)*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Bot_left][0]=(-xw*(-1./3.)-0.5)*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_El_right][0]=(-xw*(-1.))*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Mu_right][0]=(-xw*(-1.))*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Ta_right][0]=(-xw*(-1.))*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Dn_right][0]=(-xw*(-1./3.))*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Up_right][0]=(-xw*(2./3.))*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Str_right][0]=(-xw*(-1./3.))*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Chm_right][0]=(-xw*(2./3.))*2.;
    mela.selfDZpffcoupl[gHIGGS_Vp_Bot_right][0]=(-xw*(-1./3.))*2.;
    mela.selfDM_Zprime=91.1876;
    mela.selfDGa_Zprime=2.4952;
    mela.computeP(p0mplus_selfD_CT, useConstants);
    cout << "p0mplus_selfD_CT: " << p0mplus_selfD_CT << '\n' << endl;

    float p0g1prime2_selfD=0;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=1;
    mela.computeP(p0g1prime2_selfD, useConstants);
    cout << "p0g1prime2_selfD: " << p0g1prime2_selfD << '\n' << endl;

    float p0hplus_selfD=0;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1;
    mela.computeP(p0hplus_selfD, useConstants);
    cout << "p0hplus_selfD: " << p0hplus_selfD << '\n' << endl;

    float p0minus_selfD=0;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
    mela.computeP(p0minus_selfD, useConstants);
    cout << "p0minus_selfD: " << p0minus_selfD << '\n' << endl;

    float p0gzgs1prime2_selfD=0;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_1_PRIME2][0]=1;
    mela.computeP(p0gzgs1prime2_selfD, useConstants);
    cout << "p0gzgs1prime2_selfD: " << p0gzgs1prime2_selfD << '\n' << endl;

    float p0hpluszgs_selfD=0;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_2][0]=1;
    mela.computeP(p0hpluszgs_selfD, useConstants);
    cout << "p0hpluszgs_selfD: " << p0hpluszgs_selfD << '\n' << endl;

    float costhetastar = 0, costheta1 = 0, costheta2 = 0, Phi = 0, Phi1 = 0;
    float costhetastarMELA = 0, costheta1MELA = 0, costheta2MELA = 0, PhiMELA = 0, Phi1MELA = 0, m1MELA = 0, m2MELA = 0, m4lMELA = 0;
    TUtil::computeAngles(
      costhetastar,
      costheta1,
      costheta2,
      Phi,
      Phi1,

      daughters.at(2).second, daughters.at(2).first,
      daughters.at(3).second, daughters.at(3).first,
      daughters.at(0).second, daughters.at(0).first,
      daughters.at(1).second, daughters.at(1).first
    );
    mela.computeDecayAngles(
      m4lMELA,
      m1MELA,
      m2MELA,
      costheta1MELA,
      costheta2MELA,
      PhiMELA,
      costhetastarMELA,
      Phi1MELA
    );

    cout << "TUtil decay angles: " << costheta1 << " " << costheta2 << " " << Phi << " " << costhetastar << " " << Phi1 << endl;
    cout << "MELA decay angles: " << costheta1MELA << " " << costheta2MELA << " " << PhiMELA << " " << costhetastarMELA << " " << Phi1MELA << endl;
    cout << "MELA masses: " << m4lMELA << " " << m1MELA << " " << m2MELA << endl;

    mela.resetInputEvent();
  }

  cout.rdbuf(coutbuf);
  tout.close();
  mela.setVerbosity(bkpverbosity);
  cout.precision(bkpprecision);
}

void testME_VBF_JHUGen_Ping(int erg_tev=13, bool useConstants=false, shared_ptr<Mela> melaptr=nullptr){
  TString strtout = Form("testME_VBF_JHUGen_Ping_%iTeV_%i.out", erg_tev, (int)useConstants);
  ofstream tout(strtout.Data());
  streambuf* coutbuf = cout.rdbuf();
  cout.rdbuf(tout.rdbuf());

  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  if (verbosity>=TVar::DEBUG) cout << "Initializing Mela..." << endl;
  if (!melaptr) melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  //mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 79.9549392, 91.1876, 0.23119);

  int GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  const int nEntries = 6;
  double a1_array[nEntries][4] ={
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 },
    { 238.65751023078761, 9.2808858562825005, 15.827726043466324, -237.95116187061188 },
    { 101.52463181523598, 27.359569630718468, -0.90299073100241323, -97.764458892691749 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 101.67043553688544, 2.1169375132239789, 0.77953005873937187, -101.64540506443268 },
    { 24.717634703436786, -1.1722249478288802, -5.9599387484197646, -23.959684558009428 }
  };
  double a2_array[nEntries][4] ={
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 317.81904277258536, 2.5882005498984775, 21.352807448987718, -317.09037005377883 },
    { 180.10885677707822, -6.7240759244122792, 35.742176497019194, -176.39865053838915 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 95.655512770627581, -13.986023919404957, -37.960063551193414, -86.679881365440792 },
    { 49.137252081251319, -19.463268758477309, -28.879247017597017, -34.664676589120688 }
  };
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  // Decay mode does not matter, just use something
  GenLep1Id=13;
  GenLep2Id=-13;
  GenLep3Id=11;
  GenLep4Id=-11;

  for (int ev = 2; ev < 3; ev++){
    SimpleParticleCollection_t aparticles;
    TLorentzVector pAPart[2];
    pAPart[0].SetXYZT(a1_array[ev][1], a1_array[ev][2], a1_array[ev][3], a1_array[ev][0]);
    pAPart[1].SetXYZT(a2_array[ev][1], a2_array[ev][2], a2_array[ev][3], a2_array[ev][0]);
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t(0, pAPart[iap])); // q q'
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t((1-2*iap)*13, pAPart[iap])); // l- l+
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t((1-2*iap)*14, pAPart[iap])); // nub nu
    for (unsigned int iap=0; iap<1; iap++) aparticles.push_back(SimpleParticle_t(22, pAPart[iap])); // gamma

    int idOrdered[4] ={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters;
    for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters, &aparticles, (SimpleParticleCollection_t*)0, false);

    if (verbosity>=TVar::DEBUG){
      cout << "*******************************************************" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Summary of candidate " << ic << ":" << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
        cout << "*******************************************************" << endl;
      }
      cout << "*******************************************************" << endl;
      cout << endl;
    }

    int cindex;
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    TString strlh[1]={ "hadronic" };
    TString strzwh[3]={ "VBF", "Z-BF", "W-BF" };
    TVar::Production prod;
    for (unsigned int ilh=0; ilh<1; ilh++){
      TString lhapp = strlh[ilh];
      for (unsigned int izwa=0; izwa<3; izwa++){
        TString zwapp = strzwh[izwa];
        prod=TVar::JJVBF;

        cout << "*******************************************************" << endl;
        cout << "Computing MEs for " << strlh[ilh] << " " << strzwh[izwa] << endl;

        if (izwa==0){
          float p0mplus=0;
          mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, prod);
          mela.computeProdP(p0mplus, useConstants);
          cout << "p0mplus: " << p0mplus << '\n' << endl;

          float p0g1prime2=0;
          mela.setProcess(TVar::H0_g1prime2, TVar::JHUGen, prod);
          mela.computeProdP(p0g1prime2, useConstants);
          cout << "p0g1prime2: " << p0g1prime2 << '\n' << endl;

          float p0hplus=0;
          mela.setProcess(TVar::H0hplus, TVar::JHUGen, prod);
          mela.computeProdP(p0hplus, useConstants);
          cout << "p0hplus: " << p0hplus << '\n' << endl;

          float p0minus=0;
          mela.setProcess(TVar::H0minus, TVar::JHUGen, prod);
          mela.computeProdP(p0minus, useConstants);
          cout << "p0minus: " << p0minus << '\n' << endl;

          float p0gzgs1prime2=0;
          mela.setProcess(TVar::H0_Zgsg1prime2, TVar::JHUGen, prod);
          mela.computeProdP(p0gzgs1prime2, useConstants);
          cout << "p0gzgs1prime2: " << p0gzgs1prime2 << '\n' << endl;

          float p0hpluszgs=0;
          mela.setProcess(TVar::H0_Zgs, TVar::JHUGen, prod);
          mela.computeProdP(p0hpluszgs, useConstants);
          cout << "p0hpluszgs: " << p0hpluszgs << '\n' << endl;
        }

        // SelfD MEs
        float p0mplus_selfD=0;
        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
        if (izwa==1 || izwa==0) mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        else mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
        if (izwa!=0) mela.differentiate_HWW_HZZ=true;
        mela.computeProdP(p0mplus_selfD, useConstants);
        cout << "p0mplus_selfD: " << p0mplus_selfD << '\n' << endl;

        float p0mplus_selfD_CT=0;
        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
        double xw = 0.23119;
        if (izwa==1 || izwa==0){
          mela.selfDHzpzpcoupl[gHIGGS_VV_1][0]=1;
          mela.selfDZpffcoupl[gHIGGS_Vp_NuE_left][0]=0.5*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_El_left][0]=(-xw*(-1.)-0.5)*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Mu_left][0]=(-xw*(-1.)-0.5)*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Ta_left][0]=(-xw*(-1.)-0.5)*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Dn_left][0]=(-xw*(-1./3.)-0.5)*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Up_left][0]=(-xw*(2./3.)+0.5)*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Str_left][0]=(-xw*(-1./3.)-0.5)*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Chm_left][0]=(-xw*(2./3.)+0.5)*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Bot_left][0]=(-xw*(-1./3.)-0.5)*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_El_right][0]=(-xw*(-1.))*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Mu_right][0]=(-xw*(-1.))*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Ta_right][0]=(-xw*(-1.))*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Dn_right][0]=(-xw*(-1./3.))*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Up_right][0]=(-xw*(2./3.))*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Str_right][0]=(-xw*(-1./3.))*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Chm_right][0]=(-xw*(2./3.))*2.;
          mela.selfDZpffcoupl[gHIGGS_Vp_Bot_right][0]=(-xw*(-1./3.))*2.;
          mela.selfDM_Zprime=91.1876;
          mela.selfDGa_Zprime=2.4952;
        }
        if (izwa==2 || izwa==0){
          mela.selfDHwpwpcoupl[gHIGGS_VV_1][0]=1;
          mela.selfDWpffcoupl[gHIGGS_Vp_El_left][0]=sqrt(2.*(1.-xw));
          mela.selfDWpffcoupl[gHIGGS_Vp_Mu_left][0]=sqrt(2.*(1.-xw));
          mela.selfDWpffcoupl[gHIGGS_Vp_Ta_left][0]=sqrt(2.*(1.-xw));
          mela.selfDWpffcoupl[gHIGGS_Vp_Up_left][0]=sqrt(2.*(1.-xw));
          mela.selfDWpffcoupl[gHIGGS_Vp_Chm_left][0]=sqrt(2.*(1.-xw));
          mela.selfDWpffcoupl[gHIGGS_Vp_Top_left][0]=sqrt(2.*(1.-xw));
          /*
          mela.selfDWpffcoupl[gHIGGS_Vp_El_left][0]=1;
          mela.selfDWpffcoupl[gHIGGS_Vp_Mu_left][0]=1;
          mela.selfDWpffcoupl[gHIGGS_Vp_Ta_left][0]=1;
          mela.selfDWpffcoupl[gHIGGS_Vp_Up_left][0]=1;
          mela.selfDWpffcoupl[gHIGGS_Vp_Chm_left][0]=1;
          mela.selfDWpffcoupl[gHIGGS_Vp_Top_left][0]=1;
          */
          mela.selfDM_Wprime=80.399;
          mela.selfDGa_Wprime=2.085;
        }
        mela.differentiate_HWW_HZZ=true;
        mela.computeProdP(p0mplus_selfD_CT, useConstants);
        cout << "p0mplus_selfD_CT: " << p0mplus_selfD_CT << '\n' << endl;

        if (izwa==0){
          float p0g1prime2_selfD=0;
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
          mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=1;
          mela.computeProdP(p0g1prime2_selfD, useConstants);
          cout << "p0g1prime2_selfD: " << p0g1prime2_selfD << '\n' << endl;

          float p0hplus_selfD=0;
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
          mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1;
          mela.computeProdP(p0hplus_selfD, useConstants);
          cout << "p0hplus_selfD: " << p0hplus_selfD << '\n' << endl;

          float p0minus_selfD=0;
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
          mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
          mela.computeProdP(p0minus_selfD, useConstants);
          cout << "p0minus_selfD: " << p0minus_selfD << '\n' << endl;

          float p0gzgs1prime2_selfD=0;
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
          mela.selfDHzzcoupl[0][gHIGGS_ZA_1_PRIME2][0]=1;
          mela.computeProdP(p0gzgs1prime2_selfD, useConstants);
          cout << "p0gzgs1prime2_selfD: " << p0gzgs1prime2_selfD << '\n' << endl;

          float p0hpluszgs_selfD=0;
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, prod);
          mela.selfDHzzcoupl[0][gHIGGS_ZA_2][0]=1;
          mela.computeProdP(p0hpluszgs_selfD, useConstants);
          cout << "p0hpluszgs_selfD: " << p0hpluszgs_selfD << '\n' << endl;
        }

        cout << "*******************************************************" << endl;
      }
    }

    float costhetastar, costheta1re, costheta1im, costheta2re, costheta2im, Phi, Phi1, Q2V1, Q2V2;
    mela.computeVBFAngles(
      Q2V1,
      Q2V2,
      costheta1re,
      costheta2re,
      Phi,
      costhetastar,
      Phi1
    );
    cout << "Lab-frame VBF angles: " << costheta1re << " " << costheta2re << " " << Phi << " " << costhetastar << " " << Phi1 << " " << Q2V1 << " " << Q2V2 << endl;
    mela.computeVBFAngles_ComplexBoost(
      Q2V1,
      Q2V2,
      costheta1re, costheta1im,
      costheta2re, costheta2im,
      Phi,
      costhetastar,
      Phi1
    );
    cout << "Complex VBF angles: " << costheta1re << " " << costheta1im << " " << costheta2re << " " << costheta2im << " " << Phi << " " << costhetastar << " " << Phi1 << " " << Q2V1 << " " << Q2V2 << endl;

    mela.resetInputEvent();
    cout << "Removed..." << endl;
  }

  cout.rdbuf(coutbuf);
  tout.close();
  mela.setVerbosity(bkpverbosity);
}


void testME_Prop_Ping(int useMothers=0, shared_ptr<Mela> melaptr=nullptr){
  ofstream tout(TString("testME_Prop_Ping_")+(long)useMothers+".out");
  streambuf* coutbuf = cout.rdbuf();
  cout.rdbuf(tout.rdbuf());

  int erg_tev=13;
  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::DEBUG;
  if (!melaptr) {
    melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  }
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  //mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 79.9549392, 91.1876, 0.23119);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initializing" << endl;
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initialized" << endl;

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  const int nEntries = 6;
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  GenLep1Id=13;
  GenLep2Id=-13;
  GenLep3Id=11;
  GenLep4Id=-11;

  for (int ev = 2; ev < 3; ev++){
    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    // Some gymnastics to get pseudo-mothers for these events
    TLorentzVector pTotal = pOrdered[0]+pOrdered[1]+pOrdered[2]+pOrdered[3];
    TLorentzVector pTotal_dummy=pTotal;
    TLorentzVector pTotal_perp(pTotal.X(), pTotal.Y(), 0, pTotal.T());
    pTotal_dummy.Boost(-pTotal_perp.BoostVector());
    TLorentzVector pMothers[2];
    pMothers[0].SetXYZT(0, 0, (pTotal_dummy.Z()+pTotal_dummy.T())/2., (pTotal_dummy.Z()+pTotal_dummy.T())/2.);
    pMothers[1].SetXYZT(0, 0, (pTotal_dummy.Z()-pTotal_dummy.T())/2., (-pTotal_dummy.Z()+pTotal_dummy.T())/2.);
    for (int im=0; im<2; im++) pMothers[im].Boost(pTotal_perp.BoostVector());
    SimpleParticleCollection_t mothers_QQB;
    for (unsigned int im=0; im<2; im++) mothers_QQB.push_back(SimpleParticle_t(1-2*im, pMothers[im]));
    SimpleParticleCollection_t mothers_GG;
    for (unsigned int im=0; im<2; im++) mothers_GG.push_back(SimpleParticle_t(21, pMothers[im]));
    SimpleParticleCollection_t* mothersPtr=0;
    if (useMothers==1) mothersPtr=&mothers_GG;
    else if (useMothers==2) mothersPtr=&mothers_QQB;

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*)0, mothersPtr, (mothersPtr!=0));

    int cindex;
    /***** ZZ *****/
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    double mH=500.;
    cout << "Computing Prop(500 GeV) at mass=" << mela.getCurrentCandidate()->m() << endl;
    for (int ps=0; ps<=3; ps++){
      /*
      NoPropagator=0,
      RunningWidth=1,
      FixedWidth=2,
      CPS=3
      */
      float prop;
      mela.setMelaHiggsMassWidth(mH, -1, 0);
      mela.getXPropagator((TVar::ResonancePropagatorScheme)ps, prop);
      cout << "prop[" << ps << "]: " << prop << endl;
    }

    if (verbosity>=TVar::DEBUG){
      cout << "Removing Mela candidates\nSummary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }
    mela.resetInputEvent();
    cout << "Removed..." << endl;
  }

  cout.rdbuf(coutbuf);
  tout.close();
  mela.setVerbosity(bkpverbosity);
}


void testME_ProdDec_MCFM_Ordering(int iSel, int jSel, int rSel, int sSel){
  int order[2];
  TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(iSel, jSel, rSel, sSel, order);
  if (order[0]!=-1) cout << "testME_ProdDec_MCFM_Ordering::Order of particles should be "
    << order[0] << " "
    << order[1] << endl;
}

void testME_ProdDec_MCFM_JHUGen_WBFZZWW_Comparison_Ping(int motherflavor=0, int isZZWW=0 /*1==ZZ, 2==WW*/, int vbfvhchannel=0 /*0==VBF, 1==HadVH, 2==LepVH*/, int decZZWW=1 /*1==ZZ, 2==WW*/, int hasInterf=0 /*0==2l2l, 1==4l*/, shared_ptr<Mela> melaptr=nullptr){
  TString outname;
  int ZZWWdec_onevertexflag;
  if (decZZWW==1){
    outname = Form("testME_ProdDec_MCFM_JHUGen_WBFZZ_Comparison_Ping_%i_%i_%i_%s.out", motherflavor, isZZWW, vbfvhchannel, (hasInterf ? "4l" : "2l2l"));
    ZZWWdec_onevertexflag=2;
  }
  else if(decZZWW==2){
    outname = Form("testME_ProdDec_MCFM_JHUGen_WBFWW_Comparison_Ping_%i_%i_%i.out", motherflavor, isZZWW, vbfvhchannel);
    ZZWWdec_onevertexflag=1;
  }
  else return;
  if (hasInterf==1 && decZZWW!=1) return;

  struct mcfmme{
    float proddecme;
    double mearray[nmsq][nmsq];
    mcfmme(){
      proddecme=0;
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]=0; }
    }
    void add(const mcfmme& other){
      proddecme+=other.proddecme;
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]+=(other.mearray)[ii][jj]; }
    }
    void multiplyarray(const float val){
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]*=val; }
    }
    void printarray(){
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) cout << '\t' << mearray[ii][jj]; cout << endl; }
    }
  };
  struct jhume{
    float prodme;
    float decme;
    float proddecme;
    double mearray[nmsq][nmsq];
    jhume(){
      prodme=0;
      decme=0;
      proddecme=0;
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]=0; }
    }
    void multiplyarray(const float val){
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]*=val; }
    }
    void printarray(){
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) cout << '\t' << mearray[ii][jj]; cout << endl; }
    }
  };

  int erg_tev=13;
  float mPOLE=125.0;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::ERROR;

  bool doEval=false;
  int idMother[2]={ 0 };
  if (vbfvhchannel<=1){
    // VBF ZZ(+)WW
    if (motherflavor==0){
      doEval = true;
    }
    else if (motherflavor==1){
      idMother[0]=2; idMother[1]=1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0);
    }
    else if (motherflavor==2){
      idMother[0]=-2; idMother[1]=-1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0);
    }
    // VBF ZZ-only(+)WH
    else if (motherflavor==3){
      idMother[0]=2; idMother[1]=-1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel==0) ||
        (isZZWW==0 && vbfvhchannel>=1);
    }
    else if (motherflavor==4){
      idMother[0]=-2; idMother[1]=1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel==0) ||
        (isZZWW==0 && vbfvhchannel==1);
    }
    // VBF ZZ(+)ZH or WW(+)ZH
    else if (motherflavor==5){
      idMother[0]=2; idMother[1]=-2;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==6){
      idMother[0]=-2; idMother[1]=2;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==7){
      idMother[0]=1; idMother[1]=-1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==8){
      idMother[0]=-1; idMother[1]=1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    // Extra pieces
    // VBF ZZ(+)WW
    else if (motherflavor==9){
      idMother[0]=4; idMother[1]=3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0);
    }
    else if (motherflavor==10){
      idMother[0]=4; idMother[1]=1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0);
    }
    else if (motherflavor==11){
      idMother[0]=2; idMother[1]=3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0);
    }
    else if (motherflavor==12){
      idMother[0]=-4; idMother[1]=-3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0);
    }
    else if (motherflavor==13){
      idMother[0]=-4; idMother[1]=-1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0);
    }
    else if (motherflavor==14){
      idMother[0]=-2; idMother[1]=-3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0);
    }
    // VBF ZZ-only(+)WH
    else if (motherflavor==15){
      idMother[0]=4; idMother[1]=-3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel==0) ||
        (isZZWW==0 && vbfvhchannel>=1);
    }
    else if (motherflavor==16){
      idMother[0]=4; idMother[1]=-1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel==0) ||
        (isZZWW==0 && vbfvhchannel>=1);
    }
    else if (motherflavor==17){
      idMother[0]=2; idMother[1]=-3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel==0) ||
        (isZZWW==0 && vbfvhchannel>=1);
    }
    else if (motherflavor==18){
      idMother[0]=-4; idMother[1]=3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel==0) ||
        (isZZWW==0 && vbfvhchannel==1);
    }
    else if (motherflavor==19){
      idMother[0]=-4; idMother[1]=1;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel==0) ||
        (isZZWW==0 && vbfvhchannel==1);
    }
    else if (motherflavor==20){
      idMother[0]=-2; idMother[1]=3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel==0) ||
        (isZZWW==0 && vbfvhchannel==1);
    }
    // VBF ZZ(+)ZH or WW(+)ZH
    else if (motherflavor==21){
      idMother[0]=4; idMother[1]=-4;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==22){
      idMother[0]=-4; idMother[1]=4;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==23){
      idMother[0]=3; idMother[1]=-3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==24){
      idMother[0]=-3; idMother[1]=3;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==25){
      idMother[0]=5; idMother[1]=-5;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==26){
      idMother[0]=-5; idMother[1]=5;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor<=31){
      idMother[0]=-(motherflavor-26); idMother[1]=0;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1) ||
        (isZZWW==2 && vbfvhchannel>=1);
    }
    else if (motherflavor<=36){
      idMother[0]=(motherflavor-31); idMother[1]=0;
      doEval =
        (isZZWW==1 && vbfvhchannel==0) ||
        (isZZWW==2 && vbfvhchannel==0) ||
        (isZZWW==1 && vbfvhchannel>=1) ||
        (isZZWW==2 && vbfvhchannel>=1);
    }
  }
  else if (vbfvhchannel==2){
    // VBF ZZ-only(+)WH
    if (motherflavor==3){
      idMother[0]=14; idMother[1]=-13;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel>=1);
    }
    else if (motherflavor==4){
      idMother[0]=-14; idMother[1]=13;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1) ||
        (isZZWW==0 && vbfvhchannel>=1);
    }
    // VBF ZZ(+)ZH or WW(+)ZH
    else if (motherflavor==5){
      idMother[0]=14; idMother[1]=-14;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==6){
      idMother[0]=-14; idMother[1]=14;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==7){
      idMother[0]=13; idMother[1]=-13;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==8){
      idMother[0]=-13; idMother[1]=13;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
  }

  TVar::Production prod;
  if (vbfvhchannel==0) prod=TVar::JJVBF_S;
  else if (vbfvhchannel==1){
    if (
      (idMother[0]==-idMother[1] && idMother[0]!=0)
      ||
      (idMother[0]==-idMother[1] && idMother[0]==0 && isZZWW==1)
      ||
      (idMother[0]!=-idMother[1] && (idMother[0]==0 || idMother[1]==0) && isZZWW==1)
      ) prod=TVar::Had_ZH_S;
    else if (
      (TMath::Sign(1, idMother[0])==-TMath::Sign(1, idMother[1]) && abs(idMother[0])%2!=abs(idMother[1])%2)
      ||
      (idMother[0]==-idMother[1] && idMother[0]==0 && isZZWW==2)
      ||
      ((idMother[0]==0 || idMother[1]==0) && isZZWW==2)
      ) prod=TVar::Had_WH_S;
    else doEval=false;
  }
  else if (vbfvhchannel==2){
    if (idMother[0]==-idMother[1] && idMother[0]!=0) prod=TVar::Lep_ZH_S;
    else if (TMath::Sign(1, idMother[0])==-TMath::Sign(1, idMother[1]) && abs(idMother[0])%2!=abs(idMother[1])%2) prod=TVar::Lep_WH_S;
    else doEval=false;
  }
  else doEval=false;

  if (doEval){
    ofstream tout(outname.Data());
    streambuf* coutbuf = cout.rdbuf();
    cout.rdbuf(tout.rdbuf());

    if (!melaptr) {
      melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
    }
    Mela& mela = *melaptr;
    TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
    mela.setVerbosity(verbosity);

    jhume p_prod_0mplus_dec_0mplus_VAJHU;
    jhume p_prod_0minus_dec_0minus_VAJHU;
    jhume p_prod_fa3_dec_fa3_VAJHU;
    jhume p_prod_0mplusL2_dec_0mplusL2_VAJHU;
    jhume p_prod_fL2_dec_fL2_VAJHU;
    jhume p_prod_0minusZA_dec_0minusZA_VAJHU;
    jhume p_prod_fa3ZA_dec_fa3ZA_VAJHU;
    jhume p_prod_0minusAA_dec_0minusAA_VAJHU;
    jhume p_prod_fa3AA_dec_fa3AA_VAJHU;

    mcfmme p_bkg_VAMCFM, p_bkg_VAMCFM_rssum;
    mcfmme p_prod_0mplus_dec_0mplus_VAMCFM;
    mcfmme p_prod_0minus_dec_0minus_VAMCFM;
    mcfmme p_prod_0mplusL2_dec_0mplusL2_VAMCFM;
    mcfmme p_prod_fL2_dec_fL2_VAMCFM;
    mcfmme p_prod_fa3_dec_fa3_VAMCFM;
    mcfmme p_prod_0minusZA_dec_0minusZA_VAMCFM;
    mcfmme p_prod_fa3ZA_dec_fa3ZA_VAMCFM;
    mcfmme p_prod_0minusAA_dec_0minusAA_VAMCFM;
    mcfmme p_prod_fa3AA_dec_fa3AA_VAMCFM;

    float mzz = 0;
    float mjj = 0;
    float mjjzz = 0;
    float sysZ = 0.;

    float costhetastar=0;
    float costheta1=0;
    float costheta2=0;
    float Phi=0;
    float Phi1=0;
    float Q2V1=0;
    float Q2V2=0;

    float pingMom[8][4]={
      { 0, 0, 865.37881546721542, 865.37881546721542 },
      { 0, 0, -624.03396598421773, 624.03396598421773 },
      { 7.6145299215002638, -17.259247740062808, 9.4660586470659975, 21.106135714241464 },
      { 90.901719112641416, -69.683681833050798, 32.066319224729980, 118.94194752090492 },
      { 78.476352131782917, -35.264818847819797, -8.8615639484695272, 86.490881645951262 },
      { 191.68369742375290, -197.85205601463366, 100.99437243828194, 293.40746273989180 },
      { -131.59521398083137, 330.56000090294270, 437.01695094737875, 563.53440884737279 },
      { -237.08108460884614, -10.500196467375645, -329.33728782598945, 405.93194498307093 }
    };

    int GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
    if (decZZWW==1){
      GenLep1Id=13;
      GenLep2Id=-13;
      if (hasInterf==0){
        GenLep3Id=11;
        GenLep4Id=-11;
      }
      else{
        GenLep3Id=GenLep1Id;
        GenLep4Id=GenLep2Id;
      }
      mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    }
    else{
      GenLep1Id=13;
      GenLep2Id=-14;
      GenLep3Id=12;
      GenLep4Id=-11;
      mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    }
    int idOrdered[4] ={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };

    SimpleParticleCollection_t mothers;
    for (unsigned int ip=0; ip<2; ip++){
      mothers.push_back(
        SimpleParticle_t(
        0,
        TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
        )
        );
    };
    SimpleParticleCollection_t daughters;
    for (unsigned int ip=2; ip<6; ip++){
      daughters.push_back(
        SimpleParticle_t(
        idOrdered[ip-2],
        TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
        )
        );
    };
    SimpleParticleCollection_t associated;
    for (unsigned int ip=6; ip<8; ip++){
      associated.push_back(
        SimpleParticle_t(
        idMother[ip-6], // Is this wrong? No, not really. We want to check all initial particles.
        TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
        )
        );
    };
    mjj = (associated.at(0).second+associated.at(1).second).M();
    mzz = (daughters.at(0).second+daughters.at(1).second+daughters.at(2).second+daughters.at(3).second).M();

    mela.setInputEvent(&daughters, &associated, &mothers, true);

    if (mothers.size()>1){
      if (vbfvhchannel==0) TUtil::computeVBFAngles(
        costhetastar,
        costheta1,
        costheta2,
        Phi,
        Phi1,
        Q2V1,
        Q2V2,

        daughters.at(0).second, daughters.at(0).first,
        daughters.at(1).second, daughters.at(1).first,
        daughters.at(2).second, daughters.at(2).first,
        daughters.at(3).second, daughters.at(3).first,

        associated.at(0).second, associated.at(0).first,
        associated.at(1).second, associated.at(1).first,

        &(mothers.at(0).second), mothers.at(0).first,
        &(mothers.at(1).second), mothers.at(1).first
        );
      else TUtil::computeVHAngles(
        costhetastar,
        costheta1,
        costheta2,
        Phi,
        Phi1,

        daughters.at(0).second, daughters.at(0).first,
        daughters.at(1).second, daughters.at(1).first,
        daughters.at(2).second, daughters.at(2).first,
        daughters.at(3).second, daughters.at(3).first,

        associated.at(0).second, associated.at(0).first,
        associated.at(1).second, associated.at(1).first,

        &(mothers.at(0).second), mothers.at(0).first,
        &(mothers.at(1).second), mothers.at(1).first
        );
    }
    else{
      if (vbfvhchannel==0) TUtil::computeVBFAngles(
        costhetastar,
        costheta1,
        costheta2,
        Phi,
        Phi1,
        Q2V1,
        Q2V2,

        daughters.at(0).second, daughters.at(0).first,
        daughters.at(1).second, daughters.at(1).first,
        daughters.at(2).second, daughters.at(2).first,
        daughters.at(3).second, daughters.at(3).first,

        associated.at(0).second, associated.at(0).first,
        associated.at(1).second, associated.at(1).first
        );
      else TUtil::computeVHAngles(
        costhetastar,
        costheta1,
        costheta2,
        Phi,
        Phi1,

        daughters.at(0).second, daughters.at(0).first,
        daughters.at(1).second, daughters.at(1).first,
        daughters.at(2).second, daughters.at(2).first,
        daughters.at(3).second, daughters.at(3).first,

        associated.at(0).second, associated.at(0).first,
        associated.at(1).second, associated.at(1).first
        );
    }

    // Back-up CKM
    double bkpvckm_ud=TUtil::GetCKMElement(2, 1);
    double bkpvckm_us=TUtil::GetCKMElement(2, 3);
    double bkpvckm_ub=TUtil::GetCKMElement(2, 5);
    double bkpvckm_cd=TUtil::GetCKMElement(4, 1);
    double bkpvckm_cs=TUtil::GetCKMElement(4, 3);
    double bkpvckm_cb=TUtil::GetCKMElement(4, 5);
    double bkpvckm_ts=TUtil::GetCKMElement(6, 1);
    double bkpvckm_tb=TUtil::GetCKMElement(6, 3);
    double bkpvckm_td=TUtil::GetCKMElement(6, 5);
    // Set CKM to be diagonal for this test
    double invckm_ud=1, invckm_us=0, invckm_cd=0, invckm_cs=1, invckm_ts=0, invckm_tb=1, invckm_ub=0, invckm_cb=0, invckm_td=0;
    TUtil::SetCKMElements(&invckm_ud, &invckm_us, &invckm_cd, &invckm_cs, &invckm_ts, &invckm_tb, &invckm_ub, &invckm_cb, &invckm_td);

    /***** JHUGen *****/
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZINDEPENDENT);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p_prod_0mplus_dec_0mplus_VAJHU.decme, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
    mela.computeP(p_prod_0minus_dec_0minus_VAJHU.decme, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
    mela.computeP(p_prod_fa3_dec_fa3_VAJHU.decme, false);

    if (decZZWW==1){
      mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
      mela.computeP(p_prod_0minusZA_dec_0minusZA_VAJHU.decme, false);

      mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
      mela.computeP(p_prod_0minusAA_dec_0minusAA_VAJHU.decme, false);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
      mela.computeP(p_prod_fa3ZA_dec_fa3ZA_VAJHU.decme, false);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
      mela.computeP(p_prod_fa3AA_dec_fa3AA_VAJHU.decme, false);

      p_prod_0mplusL2_dec_0mplusL2_VAJHU.decme=p_prod_0mplus_dec_0mplus_VAJHU.decme;
      p_prod_fL2_dec_fL2_VAJHU.decme=p_prod_0mplus_dec_0mplus_VAJHU.decme;
    }
    else{
      mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME3][0]=10000;
      mela.computeP(p_prod_0mplusL2_dec_0mplusL2_VAJHU.decme, false);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME3][0]=10000;
      mela.computeP(p_prod_fL2_dec_fL2_VAJHU.decme, false);

      p_prod_0minusZA_dec_0minusZA_VAJHU.decme=p_prod_0mplus_dec_0mplus_VAJHU.decme;
      p_prod_0minusAA_dec_0minusAA_VAJHU.decme=p_prod_0mplus_dec_0mplus_VAJHU.decme;
      p_prod_fa3ZA_dec_fa3ZA_VAJHU.decme=p_prod_0mplus_dec_0mplus_VAJHU.decme;
      p_prod_fa3AA_dec_fa3AA_VAJHU.decme=p_prod_0mplus_dec_0mplus_VAJHU.decme;
    }

    bool computeL2WWprod=(isZZWW==2 || (decZZWW==2 && vbfvhchannel==0));
    bool computeJHUZA=true;
    if (vbfvhchannel>=1){
      if (isZZWW==2 && decZZWW==1){ // WH->ZZ
      }
      else if (isZZWW==1 && decZZWW==2){ // ZH->WW
      }
      else if (isZZWW==1 && decZZWW==1){ // ZH->ZZ
      }
      else computeJHUZA=false;
    }
    else{
      if (isZZWW==2 && decZZWW==1){ // WW->ZZ
      }
      else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
      }
      else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
      }
      else computeJHUZA=false;
    }
    if (vbfvhchannel==0){
      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
      cout << "JHUGen production chosen: " << TVar::ProductionName(TVar::JJVBF) << endl;

      if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
      else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
      mela.computeProdP(p_prod_0mplus_dec_0mplus_VAJHU.prodme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_0mplus_dec_0mplus_VAJHU.mearray);

      if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
      else{ mela.selfDHwwcoupl[0][gHIGGS_VV_4][0]=1; mela.differentiate_HWW_HZZ=true; }
      mela.computeProdP(p_prod_0minus_dec_0minus_VAJHU.prodme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_0minus_dec_0minus_VAJHU.mearray);

      if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
      else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHwwcoupl[0][gHIGGS_VV_4][0]=1; mela.differentiate_HWW_HZZ=true; }
      mela.computeProdP(p_prod_fa3_dec_fa3_VAJHU.prodme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_fa3_dec_fa3_VAJHU.mearray);

      if (computeL2WWprod){
        if (decZZWW==2 && isZZWW==1) mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        mela.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME3][0]=10000; mela.differentiate_HWW_HZZ=true;
        mela.computeProdP(p_prod_0mplusL2_dec_0mplusL2_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0mplusL2_dec_0mplusL2_VAJHU.mearray);

        if (decZZWW==2 && isZZWW==1) mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME3][0]=10000; mela.differentiate_HWW_HZZ=true;
        mela.computeProdP(p_prod_fL2_dec_fL2_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fL2_dec_fL2_VAJHU.mearray);
      }
      else{
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_0mplusL2_dec_0mplusL2_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0mplusL2_dec_0mplusL2_VAJHU.mearray);

        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_fL2_dec_fL2_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fL2_dec_fL2_VAJHU.mearray);
      }
      if (computeJHUZA){
        if (isZZWW==2 && decZZWW==1){ // WW->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_0minusZA_dec_0minusZA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0minusZA_dec_0minusZA_VAJHU.mearray);

        if (isZZWW==2 && decZZWW==1){ // WW->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_0minusAA_dec_0minusAA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0minusAA_dec_0minusAA_VAJHU.mearray);

        if (isZZWW==2 && decZZWW==1){ // WW->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_fa3ZA_dec_fa3ZA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fa3ZA_dec_fa3ZA_VAJHU.mearray);

        if (isZZWW==2 && decZZWW==1){ // WW->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_fa3AA_dec_fa3AA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fa3AA_dec_fa3AA_VAJHU.mearray);
      }
      else{
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_0minusZA_dec_0minusZA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0minusZA_dec_0minusZA_VAJHU.mearray);

        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_0minusAA_dec_0minusAA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0minusAA_dec_0minusAA_VAJHU.mearray);

        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_fa3ZA_dec_fa3ZA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fa3ZA_dec_fa3ZA_VAJHU.mearray);

        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_fa3AA_dec_fa3AA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fa3AA_dec_fa3AA_VAJHU.mearray);
      }
    }
    else{
      if (prod==TVar::Had_ZH_S || prod==TVar::Lep_ZH_S){
        if (prod==TVar::Had_ZH_S){
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_ZH);
          cout << "JHUGen production chosen: " << TVar::ProductionName(TVar::Had_ZH) << endl;
        }
        else{
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_ZH);
          cout << "JHUGen production chosen: " << TVar::ProductionName(TVar::Lep_ZH) << endl;
        }

        if (isZZWW!=2) mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        else mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=0;
        mela.computeProdP(p_prod_0mplus_dec_0mplus_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0mplus_dec_0mplus_VAJHU.mearray);

        if (isZZWW!=2) mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
        else mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=0;
        mela.computeProdP(p_prod_0minus_dec_0minus_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0minus_dec_0minus_VAJHU.mearray);

        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; }
        else{ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=0; mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=0; }
        mela.computeProdP(p_prod_fa3_dec_fa3_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fa3_dec_fa3_VAJHU.mearray);

        {
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.computeProdP(p_prod_0mplusL2_dec_0mplusL2_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_0mplusL2_dec_0mplusL2_VAJHU.mearray);

          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.computeProdP(p_prod_fL2_dec_fL2_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_fL2_dec_fL2_VAJHU.mearray);
        }
        if (computeJHUZA){
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.computeProdP(p_prod_0minusZA_dec_0minusZA_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_0minusZA_dec_0minusZA_VAJHU.mearray);

          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.computeProdP(p_prod_0minusAA_dec_0minusAA_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_0minusAA_dec_0minusAA_VAJHU.mearray);

          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.computeProdP(p_prod_fa3ZA_dec_fa3ZA_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_fa3ZA_dec_fa3ZA_VAJHU.mearray);

          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.computeProdP(p_prod_fa3AA_dec_fa3AA_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_fa3AA_dec_fa3AA_VAJHU.mearray);
        }
        else{
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.computeProdP(p_prod_0minusZA_dec_0minusZA_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_0minusZA_dec_0minusZA_VAJHU.mearray);

          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.computeProdP(p_prod_0minusAA_dec_0minusAA_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_0minusAA_dec_0minusAA_VAJHU.mearray);

          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.computeProdP(p_prod_fa3ZA_dec_fa3ZA_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_fa3ZA_dec_fa3ZA_VAJHU.mearray);

          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.computeProdP(p_prod_fa3AA_dec_fa3AA_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_fa3AA_dec_fa3AA_VAJHU.mearray);
        }
      }
      else if (prod==TVar::Had_WH_S || prod==TVar::Lep_WH_S){
        if (prod==TVar::Had_WH_S){
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_WH);
          cout << "JHUGen production chosen: " << TVar::ProductionName(TVar::Had_WH) << endl;
        }
        else{
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_WH);
          cout << "JHUGen production chosen: " << TVar::ProductionName(TVar::Lep_WH) << endl;
        }

        if (isZZWW!=1) mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        mela.computeProdP(p_prod_0mplus_dec_0mplus_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0mplus_dec_0mplus_VAJHU.mearray);

        if (isZZWW!=1) mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
        mela.computeProdP(p_prod_0minus_dec_0minus_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0minus_dec_0minus_VAJHU.mearray);

        if (isZZWW!=1){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; }
        else{ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=0; mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=0; }
        mela.computeProdP(p_prod_fa3_dec_fa3_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fa3_dec_fa3_VAJHU.mearray);

        mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        mela.computeProdP(p_prod_0minusZA_dec_0minusZA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0minusZA_dec_0minusZA_VAJHU.mearray);

        mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        mela.computeProdP(p_prod_0minusAA_dec_0minusAA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_0minusAA_dec_0minusAA_VAJHU.mearray);

        mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        mela.computeProdP(p_prod_fa3ZA_dec_fa3ZA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fa3ZA_dec_fa3ZA_VAJHU.mearray);

        mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
        mela.computeProdP(p_prod_fa3AA_dec_fa3AA_VAJHU.prodme, false);
        mela.getIORecord()->getUnweightedMEArray(p_prod_fa3AA_dec_fa3AA_VAJHU.mearray);

        {
          mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME3][0]=10000;
          mela.computeProdP(p_prod_0mplusL2_dec_0mplusL2_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_0mplusL2_dec_0mplusL2_VAJHU.mearray);

          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME3][0]=10000;
          mela.computeProdP(p_prod_fL2_dec_fL2_VAJHU.prodme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_fL2_dec_fL2_VAJHU.mearray);
        }
      }
    }

    double mh=mPOLE;
    double gah=wPOLE;
    double propagator = 1./(pow(pow(mzz, 2)-pow(mPOLE, 2), 2)+pow(mPOLE*wPOLE, 2));
    if (vbfvhchannel>=1){ // JHUGen VH pseudo-propagator
      mela.getIORecord()->getHiggsMassWidth(mh, gah, 0);
      propagator /= 1./(pow(pow(mzz, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
    }
    p_prod_0mplus_dec_0mplus_VAJHU.prodme *= propagator; p_prod_0mplus_dec_0mplus_VAJHU.multiplyarray(propagator);
    p_prod_0minus_dec_0minus_VAJHU.prodme *= propagator; p_prod_0minus_dec_0minus_VAJHU.multiplyarray(propagator);
    p_prod_0mplusL2_dec_0mplusL2_VAJHU.prodme *= propagator; p_prod_0mplusL2_dec_0mplusL2_VAJHU.multiplyarray(propagator);
    p_prod_0minusZA_dec_0minusZA_VAJHU.prodme *= propagator; p_prod_0minusZA_dec_0minusZA_VAJHU.multiplyarray(propagator);
    p_prod_0minusAA_dec_0minusAA_VAJHU.prodme *= propagator; p_prod_0minusAA_dec_0minusAA_VAJHU.multiplyarray(propagator);
    p_prod_fa3_dec_fa3_VAJHU.prodme *= propagator; p_prod_fa3_dec_fa3_VAJHU.multiplyarray(propagator);
    p_prod_fL2_dec_fL2_VAJHU.prodme *= propagator; p_prod_fL2_dec_fL2_VAJHU.multiplyarray(propagator);
    p_prod_fa3ZA_dec_fa3ZA_VAJHU.prodme *= propagator; p_prod_fa3ZA_dec_fa3ZA_VAJHU.multiplyarray(propagator);
    p_prod_fa3AA_dec_fa3AA_VAJHU.prodme *= propagator; p_prod_fa3AA_dec_fa3AA_VAJHU.multiplyarray(propagator);

    p_prod_0mplus_dec_0mplus_VAJHU.proddecme = p_prod_0mplus_dec_0mplus_VAJHU.prodme*p_prod_0mplus_dec_0mplus_VAJHU.decme; p_prod_0mplus_dec_0mplus_VAJHU.multiplyarray(p_prod_0mplus_dec_0mplus_VAJHU.decme);
    if (isZZWW==ZZWWdec_onevertexflag){ // MCFM setting below turns off ZZ anomalous couplings if WH with WW couplings tested, so we should do the same here.
      p_prod_0minus_dec_0minus_VAJHU.proddecme = p_prod_0minus_dec_0minus_VAJHU.prodme*p_prod_0mplus_dec_0mplus_VAJHU.decme; p_prod_0minus_dec_0minus_VAJHU.multiplyarray(p_prod_0mplus_dec_0mplus_VAJHU.decme);
      p_prod_fa3_dec_fa3_VAJHU.proddecme = p_prod_fa3_dec_fa3_VAJHU.prodme*p_prod_0mplus_dec_0mplus_VAJHU.decme; p_prod_fa3_dec_fa3_VAJHU.multiplyarray(p_prod_0mplus_dec_0mplus_VAJHU.decme);
    }
    else{
      p_prod_0minus_dec_0minus_VAJHU.proddecme = p_prod_0minus_dec_0minus_VAJHU.prodme*p_prod_0minus_dec_0minus_VAJHU.decme; p_prod_0minus_dec_0minus_VAJHU.multiplyarray(p_prod_0minus_dec_0minus_VAJHU.decme);
      p_prod_fa3_dec_fa3_VAJHU.proddecme = p_prod_fa3_dec_fa3_VAJHU.prodme*p_prod_fa3_dec_fa3_VAJHU.decme; p_prod_fa3_dec_fa3_VAJHU.multiplyarray(p_prod_fa3_dec_fa3_VAJHU.decme);
    }
    // These MEs only test WW anomalous couplings, so everything is fine here.
    p_prod_0mplusL2_dec_0mplusL2_VAJHU.proddecme = p_prod_0mplusL2_dec_0mplusL2_VAJHU.prodme*p_prod_0mplusL2_dec_0mplusL2_VAJHU.decme; p_prod_0mplusL2_dec_0mplusL2_VAJHU.multiplyarray(p_prod_0mplusL2_dec_0mplusL2_VAJHU.decme);
    p_prod_fL2_dec_fL2_VAJHU.proddecme = p_prod_fL2_dec_fL2_VAJHU.prodme*p_prod_fL2_dec_fL2_VAJHU.decme; p_prod_fL2_dec_fL2_VAJHU.multiplyarray(p_prod_fL2_dec_fL2_VAJHU.decme);
    // These MEs are a test of WWZZ ZZ anomalous coupling, so test here separately
    p_prod_0minusZA_dec_0minusZA_VAJHU.proddecme = p_prod_0minusZA_dec_0minusZA_VAJHU.prodme*p_prod_0minusZA_dec_0minusZA_VAJHU.decme; p_prod_0minusZA_dec_0minusZA_VAJHU.multiplyarray(p_prod_0minusZA_dec_0minusZA_VAJHU.decme);
    p_prod_0minusAA_dec_0minusAA_VAJHU.proddecme = p_prod_0minusAA_dec_0minusAA_VAJHU.prodme*p_prod_0minusAA_dec_0minusAA_VAJHU.decme; p_prod_0minusAA_dec_0minusAA_VAJHU.multiplyarray(p_prod_0minusAA_dec_0minusAA_VAJHU.decme);
    p_prod_fa3ZA_dec_fa3ZA_VAJHU.proddecme = p_prod_fa3ZA_dec_fa3ZA_VAJHU.prodme*p_prod_fa3ZA_dec_fa3ZA_VAJHU.decme; p_prod_fa3ZA_dec_fa3ZA_VAJHU.multiplyarray(p_prod_fa3ZA_dec_fa3ZA_VAJHU.decme);
    p_prod_fa3AA_dec_fa3AA_VAJHU.proddecme = p_prod_fa3AA_dec_fa3AA_VAJHU.prodme*p_prod_fa3AA_dec_fa3AA_VAJHU.decme; p_prod_fa3AA_dec_fa3AA_VAJHU.multiplyarray(p_prod_fa3AA_dec_fa3AA_VAJHU.decme);

    /***** MCFM *****/
    // Reset these in case the function needs to be repeated
    spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
    spinzerohiggs_anomcoupl_.AnomalCouplPR=1;

    if (isZZWW==ZZWWdec_onevertexflag) spinzerohiggs_anomcoupl_.AnomalCouplDK=0; // Test WW couplings in ZZ decay or ZZ couplings in WW decay
    else spinzerohiggs_anomcoupl_.AnomalCouplDK=1; // Test prod*decay couplings

    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, prod);
    cout << "MCFM production chosen: " << TVar::ProductionName(prod) << endl;
    cout << "spinzerohiggs_anomcoupl_.AnomalCouplDK=" << spinzerohiggs_anomcoupl_.AnomalCouplDK << endl;

    if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
    else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
    mela.computeProdDecP(p_prod_0mplus_dec_0mplus_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_0mplus_dec_0mplus_VAMCFM.mearray);

    if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
    else{ mela.selfDHwwcoupl[0][gHIGGS_VV_4][0]=1; mela.differentiate_HWW_HZZ=true; }
    mela.computeProdDecP(p_prod_0minus_dec_0minus_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_0minus_dec_0minus_VAMCFM.mearray);

    if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
    else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHwwcoupl[0][gHIGGS_VV_4][0]=1; mela.differentiate_HWW_HZZ=true; }
    mela.computeProdDecP(p_prod_fa3_dec_fa3_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_fa3_dec_fa3_VAMCFM.mearray);

    // Test L2 asymmetric coupling in WW as well
    bool computeL2WWproddec=(isZZWW==2 || decZZWW==2);
    if (computeL2WWproddec){
      mela.differentiate_HWW_HZZ=true;
      if (isZZWW==2 && decZZWW!=2){
        if (vbfvhchannel==0){ // If WW fusion and ZZ decay, turn off ZZ couplings, and AC in DK
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=0;
          spinzerohiggs_anomcoupl_.AnomalCouplDK=0;
          spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
        }
        else{ // If WH and ZZ decay, again do the same
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=0;
          spinzerohiggs_anomcoupl_.AnomalCouplDK=0;
          spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
        }
      }
      else if (isZZWW!=2 && decZZWW==2){
        if (vbfvhchannel==0){ // If ZZ fusion and WW decay, turn off ZZ couplings, and AC in PR
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
          spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
        }
        else{ // If ZH and WW decay, keep ZZ coupling on
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
          spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
        }
      }
      else{ // If WW in both production and decay, only turn off ZZ
        mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=0;
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      }
      mela.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME3][0]=10000;
      /*
      cerr << "mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=" << mela.selfDHzzcoupl[0][gHIGGS_VV_1][0] << endl;
      cerr << "mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=" << mela.selfDHwwcoupl[0][gHIGGS_VV_1][0] << endl;
      cerr << "mela.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME3][0]=" << mela.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME3][0] << endl;
      cerr << "mela.differentiate_HWW_HZZ=" << mela.differentiate_HWW_HZZ << endl;
      cerr << "spinzerohiggs_anomcoupl_.AnomalCouplPR=" << spinzerohiggs_anomcoupl_.AnomalCouplPR << endl;
      cerr << "spinzerohiggs_anomcoupl_.AnomalCouplDK=" << spinzerohiggs_anomcoupl_.AnomalCouplDK << endl;
      */
      mela.computeProdDecP(p_prod_0mplusL2_dec_0mplusL2_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_0mplusL2_dec_0mplusL2_VAMCFM.mearray);

      mela.differentiate_HWW_HZZ=true;
      if (isZZWW==2 && decZZWW!=2){
        if (vbfvhchannel==0){ // If WW fusion and ZZ decay, turn off ZZ couplings, and AC in DK
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=0;
          spinzerohiggs_anomcoupl_.AnomalCouplDK=0;
          spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
        }
        else{ // If WH and ZZ decay, again do the same
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=0;
          spinzerohiggs_anomcoupl_.AnomalCouplDK=0;
          spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
        }
      }
      else if (isZZWW!=2 && decZZWW==2){
        if (vbfvhchannel==0){ // If ZZ fusion and WW decay, turn off ZZ couplings, and AC in PR
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
          spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
        }
        else{ // If ZH and WW decay, keep ZZ coupling on
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
          spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
        }
      }
      else{ // If WW in both production and decay, only turn off ZZ
        mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=0;
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      }
      mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHwwcoupl[0][gHIGGS_VV_1_PRIME3][0]=10000;
      mela.computeProdDecP(p_prod_fL2_dec_fL2_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_fL2_dec_fL2_VAMCFM.mearray);

      // Reset
      spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
      spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
    }
    else{
      for (int ii=0; ii<nmsq; ii++){
        for (int jj=0; jj<nmsq; jj++){
          p_prod_0mplusL2_dec_0mplusL2_VAMCFM.mearray[ii][jj]=p_prod_0mplus_dec_0mplus_VAMCFM.mearray[ii][jj]; p_prod_0mplusL2_dec_0mplusL2_VAMCFM.proddecme=p_prod_0mplus_dec_0mplus_VAMCFM.proddecme;
          p_prod_fL2_dec_fL2_VAMCFM.mearray[ii][jj]=p_prod_0mplus_dec_0mplus_VAMCFM.mearray[ii][jj]; p_prod_fL2_dec_fL2_VAMCFM.proddecme=p_prod_0mplus_dec_0mplus_VAMCFM.proddecme;
        }
      }
    }

    // Test ZA and AA couplings in WWZZ in WBFZZ or ampvbf in WBFWW as well
    bool testZAcoupl=true;
    if (vbfvhchannel>=1){
      if (isZZWW==2 && decZZWW==1){ // WH->ZZ
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=0;
      }
      else if (isZZWW==1 && decZZWW==2){ // ZH->WW
        spinzerohiggs_anomcoupl_.AnomalCouplDK=0;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      }
      else if (isZZWW==1 && decZZWW==1){ // ZH->ZZ
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      }
      else testZAcoupl=false;
    }
    else{
      if (isZZWW==2 && decZZWW==1){ // WW->ZZ
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      }
      else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      }
      else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      }
      else testZAcoupl=false;
    }

    if (testZAcoupl){
      if (vbfvhchannel>=1){
        if (isZZWW==2 && decZZWW==1){ // WH->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZH->WW
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZH->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
      }
      else{
        if (isZZWW==2 && decZZWW==1){ // WW->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
      }
      mela.computeProdDecP(p_prod_0minusZA_dec_0minusZA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_0minusZA_dec_0minusZA_VAMCFM.mearray);

      if (vbfvhchannel>=1){
        if (isZZWW==2 && decZZWW==1){ // WH->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZH->WW
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZH->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
      }
      else{
        if (isZZWW==2 && decZZWW==1){ // WW->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
      }
      mela.computeProdDecP(p_prod_0minusAA_dec_0minusAA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_0minusAA_dec_0minusAA_VAMCFM.mearray);

      if (vbfvhchannel>=1){
        if (isZZWW==2 && decZZWW==1){ // WH->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZH->WW
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZH->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
      }
      else{
        if (isZZWW==2 && decZZWW==1){ // WW->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
      }
      mela.computeProdDecP(p_prod_fa3ZA_dec_fa3ZA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_fa3ZA_dec_fa3ZA_VAMCFM.mearray);

      if (vbfvhchannel>=1){
        if (isZZWW==2 && decZZWW==1){ // WH->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZH->WW
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZH->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
      }
      else{
        if (isZZWW==2 && decZZWW==1){ // WW->ZZ
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==2){ // ZZ->WW
          mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
        else if (isZZWW==1 && decZZWW==1){ // ZZ->ZZ
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
          mela.differentiate_HWW_HZZ=true;
        }
      }
      mela.computeProdDecP(p_prod_fa3AA_dec_fa3AA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_fa3AA_dec_fa3AA_VAMCFM.mearray);
    }
    else{
      for (int ii=0; ii<nmsq; ii++){
        for (int jj=0; jj<nmsq; jj++){
          p_prod_0minusZA_dec_0minusZA_VAMCFM.mearray[ii][jj]=p_prod_0mplus_dec_0mplus_VAMCFM.mearray[ii][jj]; p_prod_0minusZA_dec_0minusZA_VAMCFM.proddecme=p_prod_0mplus_dec_0mplus_VAMCFM.proddecme;
          p_prod_0minusAA_dec_0minusAA_VAMCFM.mearray[ii][jj]=p_prod_0mplus_dec_0mplus_VAMCFM.mearray[ii][jj]; p_prod_0minusAA_dec_0minusAA_VAMCFM.proddecme=p_prod_0mplus_dec_0mplus_VAMCFM.proddecme;
          p_prod_fa3ZA_dec_fa3ZA_VAMCFM.mearray[ii][jj]=p_prod_0mplus_dec_0mplus_VAMCFM.mearray[ii][jj]; p_prod_fa3ZA_dec_fa3ZA_VAMCFM.proddecme=p_prod_0mplus_dec_0mplus_VAMCFM.proddecme;
          p_prod_fa3AA_dec_fa3AA_VAMCFM.mearray[ii][jj]=p_prod_0mplus_dec_0mplus_VAMCFM.mearray[ii][jj]; p_prod_fa3AA_dec_fa3AA_VAMCFM.proddecme=p_prod_0mplus_dec_0mplus_VAMCFM.proddecme;
        }
      }
    }

    if (decZZWW==1) mela.setProcess(TVar::bkgZZ, TVar::MCFM, prod);
    else mela.setProcess(TVar::bkgWW, TVar::MCFM, prod);
    mela.computeProdDecP(p_bkg_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_bkg_VAMCFM.mearray);
    if (motherflavor==0 && vbfvhchannel<2){
      for (int r=-5; r<=5; r++){
        for (int s=-5; s<=5; s++){
          MELACandidate* cand = mela.getCurrentCandidate();
          int idj[2] ={
            cand->getAssociatedJet(0)->id,
            cand->getAssociatedJet(1)->id
          };

          cand->getAssociatedJet(0)->id=r;
          cand->getAssociatedJet(1)->id=s;
          mcfmme p_bkg_VAMCFM_rsindiv;
          mela.computeProdDecP(p_bkg_VAMCFM_rsindiv.proddecme, false);
          mela.getIORecord()->getUnweightedMEArray(p_bkg_VAMCFM_rsindiv.mearray);
          p_bkg_VAMCFM_rssum.add(p_bkg_VAMCFM_rsindiv);
          cand->getAssociatedJet(0)->id=idj[0];
          cand->getAssociatedJet(1)->id=idj[1];
        }
      }
    }

    cout << "Production variables:\n";
    cout << "\tmJJ = " << mjj << endl;
    cout << "\tPhi = " << Phi << endl;
    cout << "\tJHUGen (mass, width): (" << mh << ", " << gah << ")" << endl;
    cout << "\tJHUGen propagator: " << propagator << endl;
    cout << "Bkg" << endl;
    cout << "\tMCFM ME: " << p_bkg_VAMCFM.proddecme << endl;
    cout << "0mplus" << endl;
    cout << "\tJHUGen decay-alone: " << p_prod_0mplus_dec_0mplus_VAJHU.decme << endl;
    cout << "\tJHUGen prod.-alone: " << p_prod_0mplus_dec_0mplus_VAJHU.prodme << endl;
    cout << "\tJHUGen ME: " << p_prod_0mplus_dec_0mplus_VAJHU.proddecme << endl;
    cout << "\tMCFM ME: " << p_prod_0mplus_dec_0mplus_VAMCFM.proddecme << endl;
    cout << "0minus" << endl;
    cout << "\tJHUGen decay-alone: " << p_prod_0minus_dec_0minus_VAJHU.decme << endl;
    cout << "\tJHUGen prod.-alone: " << p_prod_0minus_dec_0minus_VAJHU.prodme << endl;
    cout << "\tJHUGen ME: " << p_prod_0minus_dec_0minus_VAJHU.proddecme << endl;
    cout << "\tMCFM ME: " << p_prod_0minus_dec_0minus_VAMCFM.proddecme << endl;
    cout << "fa3" << endl;
    cout << "\tJHUGen decay-alone: " << p_prod_fa3_dec_fa3_VAJHU.decme << endl;
    cout << "\tJHUGen prod.-alone: " << p_prod_fa3_dec_fa3_VAJHU.prodme << endl;
    cout << "\tJHUGen ME: " << p_prod_fa3_dec_fa3_VAJHU.proddecme << endl;
    cout << "\tMCFM ME: " << p_prod_fa3_dec_fa3_VAMCFM.proddecme << endl;
    cout << "0mplusL2" << endl;
    cout << "\tJHUGen decay-alone: " << p_prod_0mplusL2_dec_0mplusL2_VAJHU.decme << endl;
    cout << "\tJHUGen prod.-alone: " << p_prod_0mplusL2_dec_0mplusL2_VAJHU.prodme << endl;
    cout << "\tJHUGen ME: " << p_prod_0mplusL2_dec_0mplusL2_VAJHU.proddecme << endl;
    cout << "\tMCFM ME: " << p_prod_0mplusL2_dec_0mplusL2_VAMCFM.proddecme << endl;
    cout << "fL2" << endl;
    cout << "\tJHUGen decay-alone: " << p_prod_fL2_dec_fL2_VAJHU.decme << endl;
    cout << "\tJHUGen prod.-alone: " << p_prod_fL2_dec_fL2_VAJHU.prodme << endl;
    cout << "\tJHUGen ME: " << p_prod_fL2_dec_fL2_VAJHU.proddecme << endl;
    cout << "\tMCFM ME: " << p_prod_fL2_dec_fL2_VAMCFM.proddecme << endl;
    cout << "0minusZA" << endl;
    cout << "\tJHUGen decay-alone: " << p_prod_0minusZA_dec_0minusZA_VAJHU.decme << endl;
    cout << "\tJHUGen prod.-alone: " << p_prod_0minusZA_dec_0minusZA_VAJHU.prodme << endl;
    cout << "\tJHUGen ME: " << p_prod_0minusZA_dec_0minusZA_VAJHU.proddecme << endl;
    cout << "\tMCFM ME: " << p_prod_0minusZA_dec_0minusZA_VAMCFM.proddecme << endl;
    cout << "fa3ZA" << endl;
    cout << "\tJHUGen decay-alone: " << p_prod_fa3ZA_dec_fa3ZA_VAJHU.decme << endl;
    cout << "\tJHUGen prod.-alone: " << p_prod_fa3ZA_dec_fa3ZA_VAJHU.prodme << endl;
    cout << "\tJHUGen ME: " << p_prod_fa3ZA_dec_fa3ZA_VAJHU.proddecme << endl;
    cout << "\tMCFM ME: " << p_prod_fa3ZA_dec_fa3ZA_VAMCFM.proddecme << endl;
    cout << "0minusAA" << endl;
    cout << "\tJHUGen decay-alone: " << p_prod_0minusAA_dec_0minusAA_VAJHU.decme << endl;
    cout << "\tJHUGen prod.-alone: " << p_prod_0minusAA_dec_0minusAA_VAJHU.prodme << endl;
    cout << "\tJHUGen ME: " << p_prod_0minusAA_dec_0minusAA_VAJHU.proddecme << endl;
    cout << "\tMCFM ME: " << p_prod_0minusAA_dec_0minusAA_VAMCFM.proddecme << endl;
    cout << "fa3AA" << endl;
    cout << "\tJHUGen decay-alone: " << p_prod_fa3AA_dec_fa3AA_VAJHU.decme << endl;
    cout << "\tJHUGen prod.-alone: " << p_prod_fa3AA_dec_fa3AA_VAJHU.prodme << endl;
    cout << "\tJHUGen ME: " << p_prod_fa3AA_dec_fa3AA_VAJHU.proddecme << endl;
    cout << "\tMCFM ME: " << p_prod_fa3AA_dec_fa3AA_VAMCFM.proddecme << endl;

    cout << "Arrays:" << endl;
    cout << "0mplus" << endl;
    cout << "\tJHUGen" << endl;
    p_prod_0mplus_dec_0mplus_VAJHU.printarray();
    cout << "\tMCFM" << endl;
    p_prod_0mplus_dec_0mplus_VAMCFM.printarray();
    cout << "\tJHUGen/MCFM Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_prod_0mplus_dec_0mplus_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_0mplus_dec_0mplus_VAJHU.mearray[ii][jj]/p_prod_0mplus_dec_0mplus_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    cout << "0minus" << endl;
    cout << "\tJHUGen" << endl;
    p_prod_0minus_dec_0minus_VAJHU.printarray();
    cout << "\tMCFM" << endl;
    p_prod_0minus_dec_0minus_VAMCFM.printarray();
    cout << "\tJHUGen/MCFM Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_prod_0minus_dec_0minus_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_0minus_dec_0minus_VAJHU.mearray[ii][jj]/p_prod_0minus_dec_0minus_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    cout << "fa3" << endl;
    cout << "\tJHUGen" << endl;
    p_prod_fa3_dec_fa3_VAJHU.printarray();
    cout << "\tMCFM" << endl;
    p_prod_fa3_dec_fa3_VAMCFM.printarray();
    cout << "\tJHUGen/MCFM Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_prod_fa3_dec_fa3_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_fa3_dec_fa3_VAJHU.mearray[ii][jj]/p_prod_fa3_dec_fa3_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    cout << "0mplusL2" << endl;
    cout << "\tJHUGen" << endl;
    p_prod_0mplusL2_dec_0mplusL2_VAJHU.printarray();
    cout << "\tMCFM" << endl;
    p_prod_0mplusL2_dec_0mplusL2_VAMCFM.printarray();
    cout << "\tJHUGen/MCFM Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_prod_0mplusL2_dec_0mplusL2_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_0mplusL2_dec_0mplusL2_VAJHU.mearray[ii][jj]/p_prod_0mplusL2_dec_0mplusL2_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    cout << "fL2" << endl;
    cout << "\tJHUGen" << endl;
    p_prod_fL2_dec_fL2_VAJHU.printarray();
    cout << "\tMCFM" << endl;
    p_prod_fL2_dec_fL2_VAMCFM.printarray();
    cout << "\tJHUGen/MCFM Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_prod_fL2_dec_fL2_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_fL2_dec_fL2_VAJHU.mearray[ii][jj]/p_prod_fL2_dec_fL2_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    cout << "0minusZA" << endl;
    cout << "\tJHUGen" << endl;
    p_prod_0minusZA_dec_0minusZA_VAJHU.printarray();
    cout << "\tMCFM" << endl;
    p_prod_0minusZA_dec_0minusZA_VAMCFM.printarray();
    cout << "\tJHUGen/MCFM Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_prod_0minusZA_dec_0minusZA_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_0minusZA_dec_0minusZA_VAJHU.mearray[ii][jj]/p_prod_0minusZA_dec_0minusZA_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    cout << "fa3ZA" << endl;
    cout << "\tJHUGen" << endl;
    p_prod_fa3ZA_dec_fa3ZA_VAJHU.printarray();
    cout << "\tMCFM" << endl;
    p_prod_fa3ZA_dec_fa3ZA_VAMCFM.printarray();
    cout << "\tJHUGen/MCFM Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_prod_fa3ZA_dec_fa3ZA_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_fa3ZA_dec_fa3ZA_VAJHU.mearray[ii][jj]/p_prod_fa3ZA_dec_fa3ZA_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    cout << "0minusAA" << endl;
    cout << "\tJHUGen" << endl;
    p_prod_0minusAA_dec_0minusAA_VAJHU.printarray();
    cout << "\tMCFM" << endl;
    p_prod_0minusAA_dec_0minusAA_VAMCFM.printarray();
    cout << "\tJHUGen/MCFM Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_prod_0minusAA_dec_0minusAA_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_0minusAA_dec_0minusAA_VAJHU.mearray[ii][jj]/p_prod_0minusAA_dec_0minusAA_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    cout << "fa3AA" << endl;
    cout << "\tJHUGen" << endl;
    p_prod_fa3AA_dec_fa3AA_VAJHU.printarray();
    cout << "\tMCFM" << endl;
    p_prod_fa3AA_dec_fa3AA_VAMCFM.printarray();
    cout << "\tJHUGen/MCFM Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_prod_fa3AA_dec_fa3AA_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_fa3AA_dec_fa3AA_VAJHU.mearray[ii][jj]/p_prod_fa3AA_dec_fa3AA_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    cout << "Bkg" << endl;
    cout << "\tMCFM" << endl;
    p_bkg_VAMCFM.printarray();
    cout << "\tMCFM 0mplus/Bkg Ratio" << endl;
    for (int ii=0; ii<nmsq; ii++){
      for (int jj=0; jj<nmsq; jj++){
        cout << '\t';
        if (p_bkg_VAMCFM.mearray[ii][jj]!=0.) cout << p_prod_0mplus_dec_0mplus_VAJHU.mearray[ii][jj]/p_bkg_VAMCFM.mearray[ii][jj];
        else cout << 0;
      }
      cout << endl;
    }
    if (motherflavor==0 && vbfvhchannel<2){
      cout << "Bkg manual sum" << endl;
      cout << "\tMCFM" << endl;
      p_bkg_VAMCFM_rssum.printarray();
      cout << "\tMCFM Bkg (re-sum)/Bkg Ratio" << endl;
      bool wrongRatio=false;
      for (int ii=0; ii<nmsq; ii++){
        for (int jj=0; jj<nmsq; jj++){
          cout << '\t';
          float rr=0;
          if (p_bkg_VAMCFM.mearray[ii][jj]!=0.) rr = p_bkg_VAMCFM_rssum.mearray[ii][jj]/p_bkg_VAMCFM.mearray[ii][jj];
          else if (p_bkg_VAMCFM_rssum.mearray[ii][jj]!=0.) rr = -9999;
          cout << rr;
          if (rr!=0. && rr!=4.) wrongRatio=true;
        }
        cout << endl;
      }

      if (wrongRatio){ // Print all non-zero contributions in the manual sum
        for (int r=-5; r<=5; r++){
          for (int s=-5; s<=5; s++){
            MELACandidate* cand = mela.getCurrentCandidate();
            int idj[2] ={
              cand->getAssociatedJet(0)->id,
              cand->getAssociatedJet(1)->id
            };

            cand->getAssociatedJet(0)->id=r;
            cand->getAssociatedJet(1)->id=s;
            mcfmme p_bkg_VAMCFM_rsindiv;
            mela.computeProdDecP(p_bkg_VAMCFM_rsindiv.proddecme, false);
            mela.getIORecord()->getUnweightedMEArray(p_bkg_VAMCFM_rsindiv.mearray);
            if (p_bkg_VAMCFM_rsindiv.proddecme>0.){
              mela.setVerbosity(TVar::DEBUG_VERBOSE);
              mela.computeProdDecP(p_bkg_VAMCFM_rsindiv.proddecme, false);
              cout << endl;
              cout << "Outgoing id1, id2 = " << r << " , " << s << endl;

              for (int ii=0; ii<nmsq; ii++){
                for (int jj=0; jj<nmsq; jj++){
                  cout << '\t';
                  cout << p_bkg_VAMCFM_rsindiv.mearray[ii][jj];
                }
                cout << endl;
              }

              cout << endl;
              mela.setVerbosity(verbosity);
            }
            cand->getAssociatedJet(0)->id=idj[0];
            cand->getAssociatedJet(1)->id=idj[1];
          }
        }
      }

    }

    TUtil::PrintCandidateSummary(mela.getCurrentCandidate());

    mela.resetInputEvent();

    // Reset the CKM elements so that they don't transfer to the next function call
    TUtil::SetCKMElements(&bkpvckm_ud, &bkpvckm_us, &bkpvckm_cd, &bkpvckm_cs, &bkpvckm_ts, &bkpvckm_tb, &bkpvckm_ub, &bkpvckm_cb, &bkpvckm_td);

    // Reset the buffer
    cout.rdbuf(coutbuf);
    tout.close();
    mela.setVerbosity(bkpverbosity);
  }
}

void testME_ProdDec_MCFM_JHUGen_JJQCDZZWW_Comparison_Ping(int motherflavor=0, int decZZWW=1 /*1==ZZ, 2==WW*/, int hasInterf=0 /*0==2l2l, 1==4l*/, shared_ptr<Mela> melaptr=nullptr){
  if (hasInterf==1 && decZZWW==2) return;
  TString outname;
  if (decZZWW==1) outname = Form("testME_ProdDec_MCFM_JHUGen_JJQCDZZ_Comparison_Ping_%i_%s.out", motherflavor, (hasInterf ? "4l" : "2l2l"));
  else if (decZZWW==2) outname = Form("testME_ProdDec_MCFM_JHUGen_JJQCDWW_Comparison_Ping_%i.out", motherflavor);
  else return;

  struct mcfmme{
    float proddecme;
    double mearray[nmsq][nmsq];
    mcfmme(){
      proddecme=0;
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]=0; }
    }
    void add(const mcfmme& other){
      proddecme+=other.proddecme;
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]+=(other.mearray)[ii][jj]; }
    }
    void multiplyarray(const float val){
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]*=val; }
    }
    void printarray(){
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) cout << '\t' << mearray[ii][jj]; cout << endl; }
    }
  };

  int erg_tev=13;
  float mPOLE=125.0;
  float wPOLE=4.07e-3;

  bool doEval=true;
  int idMother[2]={ 0 };
  if (motherflavor==0){}
  else if (motherflavor==1){ idMother[0]=2; idMother[1]=-2; }
  else if (motherflavor==2){ idMother[0]=4; idMother[1]=-4; }
  else if (motherflavor==3){ idMother[0]=-4; idMother[1]=2; }
  else if (motherflavor==4){ idMother[0]=-2; idMother[1]=4; }
  else if (motherflavor==5){ idMother[0]=1; idMother[1]=-1; }
  else if (motherflavor==6){ idMother[0]=3; idMother[1]=-3; }
  else if (motherflavor==7){ idMother[0]=5; idMother[1]=-5; }
  else if (motherflavor==8){ idMother[0]=-3; idMother[1]=1; }
  else if (motherflavor==9){ idMother[0]=-5; idMother[1]=1; }
  else if (motherflavor==10){ idMother[0]=-1; idMother[1]=3; }
  else if (motherflavor==11){ idMother[0]=-1; idMother[1]=5; }

  else if (motherflavor==12){ idMother[0]=1; idMother[1]=1; }
  else if (motherflavor==13){ idMother[0]=2; idMother[1]=2; }
  else if (motherflavor==14){ idMother[0]=3; idMother[1]=3; }
  else if (motherflavor==15){ idMother[0]=4; idMother[1]=4; }
  else if (motherflavor==16){ idMother[0]=5; idMother[1]=5; }

  else if (motherflavor==17){ idMother[0]=1; idMother[1]=21; }
  else if (motherflavor==18){ idMother[0]=21; idMother[1]=3; }
  else if (motherflavor==19){ idMother[0]=21; idMother[1]=5; }
  else if (motherflavor==20){ idMother[0]=-1; idMother[1]=21; }
  else if (motherflavor==21){ idMother[0]=21; idMother[1]=-3; }
  else if (motherflavor==22){ idMother[0]=21; idMother[1]=-5; }
  else if (motherflavor==23){ idMother[0]=2; idMother[1]=21; }
  else if (motherflavor==24){ idMother[0]=21; idMother[1]=4; }
  else if (motherflavor==25){ idMother[0]=-2; idMother[1]=21; }
  else if (motherflavor==26){ idMother[0]=21; idMother[1]=-4; }
  else if (motherflavor==27){ idMother[0]=21; idMother[1]=21; }

  else if (motherflavor==28){ idMother[0]=1; idMother[1]=0; }
  else if (motherflavor==29){ idMother[0]=-1; idMother[1]=0; }
  else if (motherflavor==30){ idMother[0]=2; idMother[1]=0; }
  else if (motherflavor==31){ idMother[0]=-2; idMother[1]=0; }
  else if (motherflavor==32){ idMother[0]=3; idMother[1]=0; }
  else if (motherflavor==33){ idMother[0]=-3; idMother[1]=0; }
  else if (motherflavor==34){ idMother[0]=4; idMother[1]=0; }
  else if (motherflavor==35){ idMother[0]=-4; idMother[1]=0; }
  else if (motherflavor==36){ idMother[0]=5; idMother[1]=0; }
  else if (motherflavor==37){ idMother[0]=-5; idMother[1]=0; }
  else if (motherflavor==38){ idMother[0]=21; idMother[1]=0; }

  else doEval=false;

  if (doEval){
    ofstream tout(outname.Data());
    streambuf* coutbuf = cout.rdbuf();
    cout.rdbuf(tout.rdbuf());

    TVar::VerbosityLevel verbosity = TVar::ERROR;
    TVar::Production prod = TVar::JJQCD;

    if (!melaptr) {
      melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
    }
    Mela& mela = *melaptr;
    TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
    mela.setVerbosity(verbosity);

    mcfmme p_prod_JJQCD_VAMCFM, p_prod_JJQCD_VAMCFM_rssum;

    float mzz = 0;
    float mjj = 0;

    float pingMom[8][4]={
      { 0, 0, 865.37881546721542, 865.37881546721542 },
      { 0, 0, -624.03396598421773, 624.03396598421773 },
      { 7.6145299215002638, -17.259247740062808, 9.4660586470659975, 21.106135714241464 },
      { 90.901719112641416, -69.683681833050798, 32.066319224729980, 118.94194752090492 },
      { 78.476352131782917, -35.264818847819797, -8.8615639484695272, 86.490881645951262 },
      { 191.68369742375290, -197.85205601463366, 100.99437243828194, 293.40746273989180 },
      { -131.59521398083137, 330.56000090294270, 437.01695094737875, 563.53440884737279 },
      { -237.08108460884614, -10.500196467375645, -329.33728782598945, 405.93194498307093 }
    };

    int GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
    if (decZZWW==1){
      GenLep1Id=13;
      GenLep2Id=-13;
      if (hasInterf==0){
        GenLep3Id=11;
        GenLep4Id=-11;
      }
      else{
        GenLep3Id=GenLep1Id;
        GenLep4Id=GenLep2Id;
      }
      mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    }
    else{
      GenLep1Id=13;
      GenLep2Id=-14;
      GenLep3Id=12;
      GenLep4Id=-11;
      mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    }

    int idOrdered[4] ={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };

    SimpleParticleCollection_t mothers;
    for (unsigned int ip=0; ip<2; ip++){
      mothers.push_back(
        SimpleParticle_t(
        0,
        TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
        )
        );
    };
    SimpleParticleCollection_t daughters;
    for (unsigned int ip=2; ip<6; ip++){
      daughters.push_back(
        SimpleParticle_t(
        idOrdered[ip-2],
        TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
        )
        );
    };
    SimpleParticleCollection_t associated;
    for (unsigned int ip=6; ip<8; ip++){
      associated.push_back(
        SimpleParticle_t(
        idMother[ip-6], // Is this wrong? No, not really. We want to check all initial particles.
        TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
        )
        );
    };
    mjj = (associated.at(0).second+associated.at(1).second).M();
    mzz = (daughters.at(0).second+daughters.at(1).second+daughters.at(2).second+daughters.at(3).second).M();

    mela.setInputEvent(&daughters, &associated, &mothers, true);

    /***** MCFM *****/

    mela.setProcess((decZZWW==1 ? TVar::bkgZZ : TVar::bkgWW), TVar::MCFM, prod);

    mela.computeProdDecP(p_prod_JJQCD_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_JJQCD_VAMCFM.mearray);
    if (motherflavor==0){
      for (int r=-5; r<=21; r++){
        if (r>5 && r<21) continue;
        for (int s=-5; s<=21; s++){
          if (s>5 && s<21) continue;
          MELACandidate* cand = mela.getCurrentCandidate();
          int idj[2] ={
            cand->getAssociatedJet(0)->id,
            cand->getAssociatedJet(1)->id
          };

          cand->getAssociatedJet(0)->id=r;
          cand->getAssociatedJet(1)->id=s;
          mcfmme p_prod_JJQCD_VAMCFM_rsindiv;
          mela.computeProdDecP(p_prod_JJQCD_VAMCFM_rsindiv.proddecme, false);
          mela.getIORecord()->getUnweightedMEArray(p_prod_JJQCD_VAMCFM_rsindiv.mearray);
          p_prod_JJQCD_VAMCFM_rssum.add(p_prod_JJQCD_VAMCFM_rsindiv);
          cand->getAssociatedJet(0)->id=idj[0];
          cand->getAssociatedJet(1)->id=idj[1];
        }
      }
    }

    cout << "Production variables:\n";
    cout << "\tmZZ = " << mzz << endl;
    cout << "\tmJJ = " << mjj << endl;
    cout << "MEsq:" << endl;
    cout << "\tMCFM ME: " << p_prod_JJQCD_VAMCFM.proddecme << endl;
    cout << "\tArray:" << endl;
    p_prod_JJQCD_VAMCFM.printarray();
    if (motherflavor==0){
      cout << "Bkg manual sum" << endl;
      cout << "\tMCFM" << endl;
      p_prod_JJQCD_VAMCFM_rssum.printarray();
      cout << "\tMCFM Bkg (re-sum)/Bkg Ratio" << endl;
      bool wrongRatio=false;
      vector<pair<int, int>> wrongRatioArray;
      for (int ii=0; ii<nmsq; ii++){
        for (int jj=0; jj<nmsq; jj++){
          cout << '\t';
          float rr=0;
          if (p_prod_JJQCD_VAMCFM.mearray[ii][jj]!=0.) rr = p_prod_JJQCD_VAMCFM_rssum.mearray[ii][jj]/p_prod_JJQCD_VAMCFM.mearray[ii][jj];
          else if (p_prod_JJQCD_VAMCFM_rssum.mearray[ii][jj]!=0.) rr = -9999;
          cout << rr;
          if (rr!=0. && rr!=4.){
            wrongRatio=true;
            wrongRatioArray.push_back(pair<int, int>(ii, jj));
          }
        }
        cout << endl;
      }
      if (wrongRatio){ // Print all non-zero contributions in the manual sum
        for (int r=-5; r<=21; r++){
          if (r>5 && r<21) continue;
          for (int s=-5; s<=21; s++){
            if (s>5 && s<21) continue;
            MELACandidate* cand = mela.getCurrentCandidate();
            int idj[2] ={
              cand->getAssociatedJet(0)->id,
              cand->getAssociatedJet(1)->id
            };

            cand->getAssociatedJet(0)->id=r;
            cand->getAssociatedJet(1)->id=s;
            mcfmme p_prod_JJQCD_VAMCFM_rsindiv;
            mela.computeProdDecP(p_prod_JJQCD_VAMCFM_rsindiv.proddecme, false);
            mela.getIORecord()->getUnweightedMEArray(p_prod_JJQCD_VAMCFM_rsindiv.mearray);
            if (p_prod_JJQCD_VAMCFM_rsindiv.proddecme>0.){
              mela.computeProdDecP(p_prod_JJQCD_VAMCFM_rsindiv.proddecme, false);
              bool hasCorrespondance=false;
              for (auto& p:wrongRatioArray){
                for (int ii=0; ii<nmsq; ii++){
                  for (int jj=0; jj<nmsq; jj++){
                    if (p.first==ii && p.second==jj && p_prod_JJQCD_VAMCFM_rsindiv.mearray[ii][jj]!=0.) hasCorrespondance=true;
                    if (hasCorrespondance) break;
                  }
                  if (hasCorrespondance) break;
                }
                if (hasCorrespondance) break;
              }
              if (hasCorrespondance){
                //mela.setVerbosity(TVar::DEBUG_VERBOSE);
                //mela.computeProdDecP(p_prod_JJQCD_VAMCFM_rsindiv.proddecme, false);

                cout << endl;
                cout << "Outgoing id1, id2 = " << r << " , " << s << endl;

                for (int ii=0; ii<nmsq; ii++){
                  for (int jj=0; jj<nmsq; jj++){
                    cout << '\t';
                    cout << p_prod_JJQCD_VAMCFM_rsindiv.mearray[ii][jj];
                  }
                  cout << endl;
                }

                cout << endl;
              }
              mela.setVerbosity(verbosity);
            }
            cand->getAssociatedJet(0)->id=idj[0];
            cand->getAssociatedJet(1)->id=idj[1];
          }
        }
      }

    }

    TUtil::PrintCandidateSummary(mela.getCurrentCandidate());

    mela.resetInputEvent();

    cout.rdbuf(coutbuf);
    tout.close();
    mela.setVerbosity(bkpverbosity);
  }
}

void testME_ProdDec_MCFM_JHUGen_WBFZZWW_TU_Comparison_Ping(int motherflavor=0, int isZZWW=0 /*1==ZZ, 2==WW*/, int vbfvhchannel=1 /*0==VBF, 1==HadVH, 2==LepVH*/, shared_ptr<Mela> melaptr=nullptr){
  if (vbfvhchannel<1 || vbfvhchannel>2) return;
  TString outname;
  if (isZZWW==1) outname = Form("testME_ProdDec_MCFM_JHUGen_WBFZZ_TU_Comparison_Ping_%i_%i_%i.out", motherflavor, isZZWW, vbfvhchannel);
  else if (isZZWW==2) outname = Form("testME_ProdDec_MCFM_JHUGen_WBFWW_TU_Comparison_Ping_%i_%i_%i.out", motherflavor, isZZWW, vbfvhchannel);
  else return;

  struct mcfmme{
    float proddecme;
    double mearray[nmsq][nmsq];
    mcfmme(){
      proddecme=0;
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]=0; }
    }
    void add(const mcfmme& other){
      proddecme+=other.proddecme;
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]+=(other.mearray)[ii][jj]; }
    }
    void multiplyarray(const float val){
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) mearray[ii][jj]*=val; }
    }
    void printarray(){
      for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) cout << '\t' << mearray[ii][jj]; cout << endl; }
    }
  };

  int erg_tev=13;
  float mPOLE=125.0;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::ERROR;

  bool doEval=false;
  int idMother[2]={ 0 };
  if (vbfvhchannel<=1){
    // VBF ZZ-only(+)WH
    if (motherflavor==1){
      idMother[0]=2; idMother[1]=-1;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    else if (motherflavor==2){
      idMother[0]=-2; idMother[1]=1;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    // VBF ZZ(+)ZH or WW(+)ZH
    else if (motherflavor==3){
      idMother[0]=2; idMother[1]=-2;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==4){
      idMother[0]=-2; idMother[1]=2;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==5){
      idMother[0]=1; idMother[1]=-1;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==6){
      idMother[0]=-1; idMother[1]=1;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    // VBF ZZ-only(+)WH
    else if (motherflavor==7){
      idMother[0]=4; idMother[1]=-3;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    else if (motherflavor==8){
      idMother[0]=4; idMother[1]=-1;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    else if (motherflavor==9){
      idMother[0]=2; idMother[1]=-3;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    else if (motherflavor==10){
      idMother[0]=-4; idMother[1]=3;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    else if (motherflavor==11){
      idMother[0]=-4; idMother[1]=1;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    else if (motherflavor==12){
      idMother[0]=-2; idMother[1]=3;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    // VBF ZZ(+)ZH or WW(+)ZH
    else if (motherflavor==13){
      idMother[0]=4; idMother[1]=-4;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==14){
      idMother[0]=-4; idMother[1]=4;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==15){
      idMother[0]=3; idMother[1]=-3;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==16){
      idMother[0]=-3; idMother[1]=3;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==17){
      idMother[0]=5; idMother[1]=-5;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==18){
      idMother[0]=-5; idMother[1]=5;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
  }
  else if (vbfvhchannel==2){
    // VBF ZZ-only(+)WH
    if (motherflavor==1){
      idMother[0]=14; idMother[1]=-13;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    else if (motherflavor==2){
      idMother[0]=-14; idMother[1]=13;
      doEval =
        (isZZWW==2 && vbfvhchannel>=1);
    }
    // VBF ZZ(+)ZH or WW(+)ZH
    else if (motherflavor==3){
      idMother[0]=14; idMother[1]=-14;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==4){
      idMother[0]=-14; idMother[1]=14;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==5){
      idMother[0]=13; idMother[1]=-13;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
    else if (motherflavor==6){
      idMother[0]=-13; idMother[1]=13;
      doEval =
        (isZZWW==1 && vbfvhchannel>=1);
    }
  }

  TVar::Production prod, prod_tu, prod_stu;
  if (vbfvhchannel==1){
    if (
      (idMother[0]==-idMother[1] && idMother[0]!=0)
      ) prod=TVar::Had_ZH_S;
    else if (
      (TMath::Sign(1, idMother[0])==-TMath::Sign(1, idMother[1]) && abs(idMother[0])%2!=abs(idMother[1])%2)
      ) prod=TVar::Had_WH_S;
    else doEval=false;
  }
  else if (vbfvhchannel==2){
    if (idMother[0]==-idMother[1] && idMother[0]!=0) prod=TVar::Lep_ZH_S;
    else if (TMath::Sign(1, idMother[0])==-TMath::Sign(1, idMother[1]) && abs(idMother[0])%2!=abs(idMother[1])%2) prod=TVar::Lep_WH_S;
    else doEval=false;
  }
  else doEval=false;

  if (prod==TVar::Had_ZH_S){
    prod_tu=TVar::Had_ZH_TU;
    prod_stu=TVar::Had_ZH;
  }
  else if (prod==TVar::Had_WH_S){
    prod_tu=TVar::Had_WH_TU;
    prod_stu=TVar::Had_WH;
  }
  else if (prod==TVar::Lep_ZH_S){
    prod_tu=TVar::Lep_ZH_TU;
    prod_stu=TVar::Lep_ZH;
  }
  else if (prod==TVar::Lep_WH_S){
    prod_tu=TVar::Lep_WH_TU;
    prod_stu=TVar::Lep_WH;
  }
  else doEval=false;

  if (doEval){
    ofstream tout(outname.Data());
    streambuf* coutbuf = cout.rdbuf();
    cout.rdbuf(tout.rdbuf());

    if (!melaptr) {
      melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
    }
    Mela& mela = *melaptr;
    TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
    mela.setVerbosity(verbosity);

    mcfmme p_bkg_VAMCFM;
    mcfmme p_prod_0mplus_dec_0mplus_VAMCFM;
    mcfmme p_prod_0minus_dec_0minus_VAMCFM;
    mcfmme p_prod_fa3_dec_fa3_VAMCFM;
    mcfmme p_prod_0minusZA_dec_0minusZA_VAMCFM;
    mcfmme p_prod_fa3ZA_dec_fa3ZA_VAMCFM;
    mcfmme p_prod_0minusAA_dec_0minusAA_VAMCFM;
    mcfmme p_prod_fa3AA_dec_fa3AA_VAMCFM;

    mcfmme p_bkg_tu_VAMCFM;
    mcfmme p_prod_tu_0mplus_dec_0mplus_VAMCFM;
    mcfmme p_prod_tu_0minus_dec_0minus_VAMCFM;
    mcfmme p_prod_tu_fa3_dec_fa3_VAMCFM;
    mcfmme p_prod_tu_0minusZA_dec_0minusZA_VAMCFM;
    mcfmme p_prod_tu_fa3ZA_dec_fa3ZA_VAMCFM;
    mcfmme p_prod_tu_0minusAA_dec_0minusAA_VAMCFM;
    mcfmme p_prod_tu_fa3AA_dec_fa3AA_VAMCFM;

    mcfmme p_bkg_stu_VAMCFM;
    mcfmme p_prod_stu_0mplus_dec_0mplus_VAMCFM;
    mcfmme p_prod_stu_0minus_dec_0minus_VAMCFM;
    mcfmme p_prod_stu_fa3_dec_fa3_VAMCFM;
    mcfmme p_prod_stu_0minusZA_dec_0minusZA_VAMCFM;
    mcfmme p_prod_stu_fa3ZA_dec_fa3ZA_VAMCFM;
    mcfmme p_prod_stu_0minusAA_dec_0minusAA_VAMCFM;
    mcfmme p_prod_stu_fa3AA_dec_fa3AA_VAMCFM;

    float mzz = 0;
    float mjj = 0;
    float mjjzz = 0;
    float sysZ = 0.;

    float costhetastar=0;
    float costheta1=0;
    float costheta2=0;
    float Phi=0;
    float Phi1=0;
    float Q2V1=0;
    float Q2V2=0;

    float pingMom[8][4]={
      { 0, 0, 865.37881546721542, 865.37881546721542 },
      { 0, 0, -624.03396598421773, 624.03396598421773 },
      { 7.6145299215002638, -17.259247740062808, 9.4660586470659975, 21.106135714241464 },
      { 90.901719112641416, -69.683681833050798, 32.066319224729980, 118.94194752090492 },
      { 78.476352131782917, -35.264818847819797, -8.8615639484695272, 86.490881645951262 },
      { 191.68369742375290, -197.85205601463366, 100.99437243828194, 293.40746273989180 },
      { -131.59521398083137, 330.56000090294270, 437.01695094737875, 563.53440884737279 },
      { -237.08108460884614, -10.500196467375645, -329.33728782598945, 405.93194498307093 }
    };

    int GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
    if (isZZWW==1){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
      mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    }
    else{
      GenLep1Id=13;
      GenLep2Id=-14;
      GenLep3Id=12;
      GenLep4Id=-11;
      mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    }
    int idOrdered[4] ={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };

    SimpleParticleCollection_t mothers;
    for (unsigned int ip=0; ip<2; ip++){
      mothers.push_back(
        SimpleParticle_t(
        0,
        TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
        )
        );
    };
    SimpleParticleCollection_t daughters;
    for (unsigned int ip=2; ip<6; ip++){
      daughters.push_back(
        SimpleParticle_t(
        idOrdered[ip-2],
        TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
        )
        );
    };
    SimpleParticleCollection_t associated;
    for (unsigned int ip=6; ip<8; ip++){
      associated.push_back(
        SimpleParticle_t(
        idMother[ip-6], // Is this wrong? No, not really. We want to check all initial particles.
        TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
        )
        );
    };
    mjj = (associated.at(0).second+associated.at(1).second).M();
    mzz = (daughters.at(0).second+daughters.at(1).second+daughters.at(2).second+daughters.at(3).second).M();

    mela.setInputEvent(&daughters, &associated, &mothers, true);

    SimpleParticleCollection_t daughters_tu;
    SimpleParticleCollection_t associated_tu;
    if (isZZWW==1 || (isZZWW==2 && (idMother[0]+idMother[1])>0)){
      for (unsigned int ip=2; ip<4; ip++){
        daughters_tu.push_back(
          SimpleParticle_t(
          idOrdered[ip-2],
          TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
          )
          );
      };
      for (unsigned int ip=4; ip<6; ip++){
        associated_tu.push_back(
          SimpleParticle_t(
          idOrdered[ip-2],
          TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
          )
          );
      };
      for (unsigned int ip=6; ip<8; ip++){
        daughters_tu.push_back(
          SimpleParticle_t(
          idMother[ip-6], // Is this wrong? No, not really. We want to check all initial particles.
          TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
          )
          );
      };
    }
    else{
      for (unsigned int ip=2; ip<4; ip++){
        associated_tu.push_back(
          SimpleParticle_t(
          idOrdered[ip-2],
          TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
          )
          );
      };
      for (unsigned int ip=4; ip<6; ip++){
        daughters_tu.push_back(
          SimpleParticle_t(
          idOrdered[ip-2],
          TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
          )
          );
      };
      for (unsigned int ip=6; ip<8; ip++){
        daughters_tu.push_back(
          SimpleParticle_t(
          idMother[ip-6], // Is this wrong? No, not really. We want to check all initial particles.
          TLorentzVector(pingMom[ip][0], pingMom[ip][1], pingMom[ip][2], pingMom[ip][3])
          )
          );
      };
    }
    mela.setInputEvent(&daughters_tu, &associated_tu, &mothers, true);

    if (mothers.size()>1){
      if (vbfvhchannel==0) TUtil::computeVBFAngles(
        costhetastar,
        costheta1,
        costheta2,
        Phi,
        Phi1,
        Q2V1,
        Q2V2,

        daughters.at(0).second, daughters.at(0).first,
        daughters.at(1).second, daughters.at(1).first,
        daughters.at(2).second, daughters.at(2).first,
        daughters.at(3).second, daughters.at(3).first,

        associated.at(0).second, associated.at(0).first,
        associated.at(1).second, associated.at(1).first,

        &(mothers.at(0).second), mothers.at(0).first,
        &(mothers.at(1).second), mothers.at(1).first
        );
      else TUtil::computeVHAngles(
        costhetastar,
        costheta1,
        costheta2,
        Phi,
        Phi1,

        daughters.at(0).second, daughters.at(0).first,
        daughters.at(1).second, daughters.at(1).first,
        daughters.at(2).second, daughters.at(2).first,
        daughters.at(3).second, daughters.at(3).first,

        associated.at(0).second, associated.at(0).first,
        associated.at(1).second, associated.at(1).first,

        &(mothers.at(0).second), mothers.at(0).first,
        &(mothers.at(1).second), mothers.at(1).first
        );
    }
    else{
      if (vbfvhchannel==0) TUtil::computeVBFAngles(
        costhetastar,
        costheta1,
        costheta2,
        Phi,
        Phi1,
        Q2V1,
        Q2V2,

        daughters.at(0).second, daughters.at(0).first,
        daughters.at(1).second, daughters.at(1).first,
        daughters.at(2).second, daughters.at(2).first,
        daughters.at(3).second, daughters.at(3).first,

        associated.at(0).second, associated.at(0).first,
        associated.at(1).second, associated.at(1).first
        );
      else TUtil::computeVHAngles(
        costhetastar,
        costheta1,
        costheta2,
        Phi,
        Phi1,

        daughters.at(0).second, daughters.at(0).first,
        daughters.at(1).second, daughters.at(1).first,
        daughters.at(2).second, daughters.at(2).first,
        daughters.at(3).second, daughters.at(3).first,

        associated.at(0).second, associated.at(0).first,
        associated.at(1).second, associated.at(1).first
        );
    }

    // Back-up CKM
    double bkpvckm_ud=TUtil::GetCKMElement(2, 1);
    double bkpvckm_us=TUtil::GetCKMElement(2, 3);
    double bkpvckm_ub=TUtil::GetCKMElement(2, 5);
    double bkpvckm_cd=TUtil::GetCKMElement(4, 1);
    double bkpvckm_cs=TUtil::GetCKMElement(4, 3);
    double bkpvckm_cb=TUtil::GetCKMElement(4, 5);
    double bkpvckm_ts=TUtil::GetCKMElement(6, 1);
    double bkpvckm_tb=TUtil::GetCKMElement(6, 3);
    double bkpvckm_td=TUtil::GetCKMElement(6, 5);
    // Set CKM to be diagonal for this test
    double invckm_ud=1, invckm_us=0, invckm_cd=0, invckm_cs=1, invckm_ts=0, invckm_tb=1, invckm_ub=0, invckm_cb=0, invckm_td=0;
    TUtil::SetCKMElements(&invckm_ud, &invckm_us, &invckm_cd, &invckm_cs, &invckm_ts, &invckm_tb, &invckm_ub, &invckm_cb, &invckm_td);
    MELACandidate* cand;

    /***** MCFM *****/
    mela.setCurrentCandidateFromIndex(1);
    cand = mela.getCurrentCandidate();
    TUtil::PrintCandidateSummary(cand);
    if ((prod==TVar::Had_WH_S || prod==TVar::Had_ZH_S) && (PDGHelpers::isALepton(cand->getSortedV(2)->getDaughter(0)->id) || PDGHelpers::isANeutrino(cand->getSortedV(2)->getDaughter(0)->id))){
      if (prod==TVar::Had_WH_S) prod=TVar::Lep_WH_S;
      else prod=TVar::Lep_ZH_S;
    }
    else if ((prod==TVar::Lep_WH_S || prod==TVar::Lep_ZH_S) && PDGHelpers::isAJet(cand->getSortedV(2)->getDaughter(0)->id)){
      if (prod==TVar::Lep_WH_S) prod=TVar::Had_WH_S;
      else prod=TVar::Had_ZH_S;
    }

    // Reset these in case the function needs to be repeated
    spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
    spinzerohiggs_anomcoupl_.AnomalCouplPR=1;

    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, prod);
    mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
    if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
    else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
    mela.computeProdDecP(p_prod_0mplus_dec_0mplus_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_0mplus_dec_0mplus_VAMCFM.mearray);

    mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
    if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
    else{ mela.selfDHwwcoupl[0][gHIGGS_VV_4][0]=1; mela.differentiate_HWW_HZZ=true; }
    mela.computeProdDecP(p_prod_0minus_dec_0minus_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_0minus_dec_0minus_VAMCFM.mearray);

    mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
    if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
    else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHwwcoupl[0][gHIGGS_VV_4][0]=1; mela.differentiate_HWW_HZZ=true; }
    mela.computeProdDecP(p_prod_fa3_dec_fa3_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_fa3_dec_fa3_VAMCFM.mearray);

    // Test ZA and AA couplings in WWZZ in WBFZZ or ampvbf in WBFWW as well
    bool testZAcoupl=true;
    if (vbfvhchannel>=1){
      if (isZZWW==1){ // ZH->ZZ
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      }
      else testZAcoupl=false;
    }

    if (testZAcoupl){
      mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
      mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
      mela.differentiate_HWW_HZZ=true;
      mela.computeProdDecP(p_prod_0minusZA_dec_0minusZA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_0minusZA_dec_0minusZA_VAMCFM.mearray);

      mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
      mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
      mela.differentiate_HWW_HZZ=true;
      mela.computeProdDecP(p_prod_0minusAA_dec_0minusAA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_0minusAA_dec_0minusAA_VAMCFM.mearray);

      mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
      mela.differentiate_HWW_HZZ=true;
      mela.computeProdDecP(p_prod_fa3ZA_dec_fa3ZA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_fa3ZA_dec_fa3ZA_VAMCFM.mearray);

      mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
      mela.differentiate_HWW_HZZ=true;
      mela.computeProdDecP(p_prod_fa3AA_dec_fa3AA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_fa3AA_dec_fa3AA_VAMCFM.mearray);
    }

    //mela.setVerbosity(TVar::DEBUG_VERBOSE);
    mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
    if (isZZWW==1) mela.setProcess(TVar::bkgZZ, TVar::MCFM, prod);
    else mela.setProcess(TVar::bkgWW, TVar::MCFM, prod);
    mela.computeProdDecP(p_bkg_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_bkg_VAMCFM.mearray);
    //mela.setVerbosity(verbosity);


    mela.setCurrentCandidateFromIndex(0);
    cand = mela.getCurrentCandidate();
    TUtil::PrintCandidateSummary(cand);
    if ((prod_tu==TVar::Had_WH_TU || prod_tu==TVar::Had_ZH_TU) && (PDGHelpers::isALepton(cand->getSortedV(2)->getDaughter(0)->id) || PDGHelpers::isANeutrino(cand->getSortedV(2)->getDaughter(0)->id))){
      if (prod_tu==TVar::Had_WH_TU) prod_tu=TVar::Lep_WH_TU;
      else prod_tu=TVar::Lep_ZH_TU;
    }
    else if ((prod_tu==TVar::Lep_WH_TU || prod_tu==TVar::Lep_ZH_TU) && PDGHelpers::isAJet(cand->getSortedV(2)->getDaughter(0)->id)){
      if (prod_tu==TVar::Lep_WH_TU) prod_tu=TVar::Had_WH_TU;
      else prod_tu=TVar::Had_ZH_TU;
    }

    // Reset these in case the function needs to be repeated
    spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
    spinzerohiggs_anomcoupl_.AnomalCouplPR=1;

    cout << "MCFM production chosen: " << TVar::ProductionName(prod_tu) << endl;
    cout << "spinzerohiggs_anomcoupl_.AnomalCouplDK=" << spinzerohiggs_anomcoupl_.AnomalCouplDK << endl;

    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, prod_tu);
    mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
    if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
    else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.differentiate_HWW_HZZ=true; }
    mela.computeProdDecP(p_prod_tu_0mplus_dec_0mplus_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_tu_0mplus_dec_0mplus_VAMCFM.mearray);

    mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
    if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
    else{ mela.selfDHwwcoupl[0][gHIGGS_VV_4][0]=1; mela.differentiate_HWW_HZZ=true; }
    mela.computeProdDecP(p_prod_tu_0minus_dec_0minus_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_tu_0minus_dec_0minus_VAMCFM.mearray);

    mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
    if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
    else{ mela.selfDHwwcoupl[0][gHIGGS_VV_1][0]=1; mela.selfDHwwcoupl[0][gHIGGS_VV_4][0]=1; mela.differentiate_HWW_HZZ=true; }
    mela.computeProdDecP(p_prod_tu_fa3_dec_fa3_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_prod_tu_fa3_dec_fa3_VAMCFM.mearray);

    // Test ZA and AA couplings in WWZZ in WBFZZ or ampvbf in WBFWW as well
    testZAcoupl=true;
    if (vbfvhchannel>=1){
      if (isZZWW==1){ // ZH->ZZ
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      }
      else testZAcoupl=false;
    }

    if (testZAcoupl){
      mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
      mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
      mela.differentiate_HWW_HZZ=true;
      mela.computeProdDecP(p_prod_tu_0minusZA_dec_0minusZA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_tu_0minusZA_dec_0minusZA_VAMCFM.mearray);

      mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
      mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
      mela.differentiate_HWW_HZZ=true;
      mela.computeProdDecP(p_prod_tu_0minusAA_dec_0minusAA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_tu_0minusAA_dec_0minusAA_VAMCFM.mearray);

      mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1;
      mela.differentiate_HWW_HZZ=true;
      mela.computeProdDecP(p_prod_tu_fa3ZA_dec_fa3ZA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_tu_fa3ZA_dec_fa3ZA_VAMCFM.mearray);

      mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1;
      mela.differentiate_HWW_HZZ=true;
      mela.computeProdDecP(p_prod_tu_fa3AA_dec_fa3AA_VAMCFM.proddecme, false);
      mela.getIORecord()->getUnweightedMEArray(p_prod_tu_fa3AA_dec_fa3AA_VAMCFM.mearray);
    }

    //mela.setVerbosity(TVar::DEBUG_VERBOSE);
    mela.setRenFacScaleMode(TVar::Fixed_mH, TVar::Fixed_mH, 1, 1);
    if (isZZWW==1) mela.setProcess(TVar::bkgZZ, TVar::MCFM, prod_tu);
    else mela.setProcess(TVar::bkgWW, TVar::MCFM, prod_tu);
    mela.computeProdDecP(p_bkg_tu_VAMCFM.proddecme, false);
    mela.getIORecord()->getUnweightedMEArray(p_bkg_tu_VAMCFM.mearray);
    //mela.setVerbosity(verbosity);

    // Reset these in case the function needs to be repeated
    spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
    spinzerohiggs_anomcoupl_.AnomalCouplPR=1;

    cout << "Production variables:\n";
    cout << "\tmJJ = " << mjj << endl;
    cout << "\tPhi = " << Phi << endl;
    cout << "Bkg" << endl;
    cout << "\tMCFM s ME: " << p_bkg_VAMCFM.proddecme << endl;
    cout << "\tMCFM tu ME: " << p_bkg_tu_VAMCFM.proddecme << endl;
    cout << "0mplus" << endl;
    cout << "\tMCFM s ME: " << p_prod_0mplus_dec_0mplus_VAMCFM.proddecme << endl;
    cout << "\tMCFM tu ME: " << p_prod_tu_0mplus_dec_0mplus_VAMCFM.proddecme << endl;
    cout << "0minus" << endl;
    cout << "\tMCFM s ME: " << p_prod_0minus_dec_0minus_VAMCFM.proddecme << endl;
    cout << "\tMCFM tu ME: " << p_prod_tu_0minus_dec_0minus_VAMCFM.proddecme << endl;
    cout << "fa3" << endl;
    cout << "\tMCFM s ME: " << p_prod_fa3_dec_fa3_VAMCFM.proddecme << endl;
    cout << "\tMCFM tu ME: " << p_prod_tu_fa3_dec_fa3_VAMCFM.proddecme << endl;
    cout << "0minusZA" << endl;
    cout << "\tMCFM s ME: " << p_prod_0minusZA_dec_0minusZA_VAMCFM.proddecme << endl;
    cout << "\tMCFM tu ME: " << p_prod_tu_0minusZA_dec_0minusZA_VAMCFM.proddecme << endl;
    cout << "fa3ZA" << endl;
    cout << "\tMCFM s ME: " << p_prod_fa3ZA_dec_fa3ZA_VAMCFM.proddecme << endl;
    cout << "\tMCFM tu ME: " << p_prod_tu_fa3ZA_dec_fa3ZA_VAMCFM.proddecme << endl;
    cout << "0minusAA" << endl;
    cout << "\tMCFM s ME: " << p_prod_0minusAA_dec_0minusAA_VAMCFM.proddecme << endl;
    cout << "\tMCFM tu ME: " << p_prod_tu_0minusAA_dec_0minusAA_VAMCFM.proddecme << endl;
    cout << "fa3AA" << endl;
    cout << "\tMCFM s ME: " << p_prod_fa3AA_dec_fa3AA_VAMCFM.proddecme << endl;
    cout << "\tMCFM tu ME: " << p_prod_tu_fa3AA_dec_fa3AA_VAMCFM.proddecme << endl;

    cout << "Arrays:" << endl;
    cout << "0mplus" << endl;
    cout << "\tMCFM s" << endl;
    p_prod_0mplus_dec_0mplus_VAMCFM.printarray();
    cout << "\tMCFM tu" << endl;
    p_prod_tu_0mplus_dec_0mplus_VAMCFM.printarray();
    cout << "0minus" << endl;
    cout << "\tMCFM s" << endl;
    p_prod_0minus_dec_0minus_VAMCFM.printarray();
    cout << "\tMCFM tu" << endl;
    p_prod_tu_0minus_dec_0minus_VAMCFM.printarray();
    cout << "fa3" << endl;
    cout << "\tMCFM s" << endl;
    p_prod_fa3_dec_fa3_VAMCFM.printarray();
    cout << "\tMCFM tu" << endl;
    p_prod_tu_fa3_dec_fa3_VAMCFM.printarray();
    cout << "0minusZA" << endl;
    cout << "\tMCFM s" << endl;
    p_prod_0minusZA_dec_0minusZA_VAMCFM.printarray();
    cout << "\tMCFM tu" << endl;
    p_prod_tu_0minusZA_dec_0minusZA_VAMCFM.printarray();
    cout << "fa3ZA" << endl;
    cout << "\tMCFM s" << endl;
    p_prod_fa3ZA_dec_fa3ZA_VAMCFM.printarray();
    cout << "\tMCFM tu" << endl;
    p_prod_tu_fa3ZA_dec_fa3ZA_VAMCFM.printarray();
    cout << "0minusAA" << endl;
    cout << "\tMCFM s" << endl;
    p_prod_0minusAA_dec_0minusAA_VAMCFM.printarray();
    cout << "\tMCFM tu" << endl;
    p_prod_tu_0minusAA_dec_0minusAA_VAMCFM.printarray();
    cout << "fa3AA" << endl;
    cout << "\tMCFM s" << endl;
    p_prod_fa3AA_dec_fa3AA_VAMCFM.printarray();
    cout << "\tMCFM tu" << endl;
    p_prod_tu_fa3AA_dec_fa3AA_VAMCFM.printarray();
    cout << "Bkg" << endl;
    cout << "\tMCFM s" << endl;
    p_bkg_VAMCFM.printarray();
    cout << "\tMCFM tu" << endl;
    p_bkg_tu_VAMCFM.printarray();

    mela.resetInputEvent();

    // Reset the CKM elements so that they don't transfer to the next function call
    TUtil::SetCKMElements(&bkpvckm_ud, &bkpvckm_us, &bkpvckm_cd, &bkpvckm_cs, &bkpvckm_ts, &bkpvckm_tb, &bkpvckm_ub, &bkpvckm_cb, &bkpvckm_td);

    // Reset the buffer
    cout.rdbuf(coutbuf);
    tout.close();
    mela.setVerbosity(bkpverbosity);
  }
}


void testME_ProdDec_MCFM_Ping(int flavor=2, shared_ptr<Mela> melaptr=nullptr){
  ofstream tout(TString("testME_ProdDec_MCFM_Ping_")+long(flavor)+".out");
  streambuf* coutbuf = cout.rdbuf();
  cout.rdbuf(tout.rdbuf());

  int erg_tev=13;
  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::DEBUG;
  if (!melaptr) melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  //mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 79.9549392, 91.1876, 0.23119);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initializing" << endl;
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  if (verbosity>=TVar::DEBUG) cout << "Mela candidate decay mode initialized" << endl;

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  const int nEntries = 6;
  double a1_array[nEntries][4] ={
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 },
    { 238.65751023078761, 9.2808858562825005, 15.827726043466324, -237.95116187061188 },
    { 101.52463181523598, 27.359569630718468, -0.90299073100241323, -97.764458892691749 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 101.67043553688544, 2.1169375132239789, 0.77953005873937187, -101.64540506443268 },
    { 24.717634703436786, -1.1722249478288802, -5.9599387484197646, -23.959684558009428 }
  };
  double a2_array[nEntries][4] ={
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 317.81904277258536, 2.5882005498984775, 21.352807448987718, -317.09037005377883 },
    { 180.10885677707822, -6.7240759244122792, 35.742176497019194, -176.39865053838915 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 95.655512770627581, -13.986023919404957, -37.960063551193414, -86.679881365440792 },
    { 49.137252081251319, -19.463268758477309, -28.879247017597017, -34.664676589120688 }
  };
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  if (flavor == 2){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 1){
    GenLep1Id=11;
    GenLep2Id=-11;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 0){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 3){
    GenLep1Id=14;
    GenLep2Id=-14;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 4){
    GenLep1Id=0;
    GenLep2Id=-0;
    GenLep3Id=1;
    GenLep4Id=-1;
  }
  else if (flavor == 5){
    GenLep1Id=1;
    GenLep2Id=-1;
    GenLep3Id=2;
    GenLep4Id=-2;
  }
  else if (flavor == 6){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=2;
    GenLep4Id=-2;
  }

  for (int ev = 0; ev < 1; ev++){
    SimpleParticleCollection_t aparticles;
    TLorentzVector pAPart[4];
    pAPart[0].SetXYZT(a1_array[ev][1], a1_array[ev][2], a1_array[ev][3], a1_array[ev][0]);
    pAPart[1].SetXYZT(a2_array[ev][1], a2_array[ev][2], a2_array[ev][3], a2_array[ev][0]);
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t(1-2*iap, pAPart[iap]));
    for (unsigned int iap=1; iap<2; iap++) aparticles.push_back(SimpleParticle_t(-2*iap, pAPart[iap]));
    for (unsigned int iap=0; iap<2; iap++) aparticles.push_back(SimpleParticle_t((1-2*iap)*13, pAPart[iap]));

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    int idOrdered_WW[4] ={ 11, -12, -11, 12 };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters_WW;
    for (unsigned int idau=0; idau<4; idau++) daughters_WW.push_back(SimpleParticle_t(idOrdered_WW[idau], pOrdered[idau]));
    SimpleParticleCollection_t daughters_WWasZZ;
    for (unsigned int iv=0; iv<2; iv++){ for (int ivd=0; ivd<2; ivd++) daughters_WWasZZ.push_back(SimpleParticle_t(idOrdered_WW[iv+2*ivd], pOrdered[iv+2*ivd])); }
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela.setInputEvent(&daughters_WW, &aparticles, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_WWasZZ, &aparticles, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_ZZ, &aparticles, (SimpleParticleCollection_t*)0, false);

    const unsigned int nprods=12;
    TVar::Production theProds[12]={
      TVar::JJVBF, TVar::JJVBF_S, TVar::JJVBF_TU,
      TVar::Had_ZH, TVar::Had_ZH_S, TVar::Had_ZH_TU,
      TVar::Had_WH, TVar::Had_WH_S, TVar::Had_WH_TU,
      TVar::JJEW, TVar::JJEW_S, TVar::JJEW_TU
    };
    int cindex;
    for (unsigned int iprod=0; iprod<nprods; iprod++){
      /***** WW *****/
      cindex=0;
      mela.setCurrentCandidateFromIndex(cindex);

      TString theProdName=TVar::ProductionName(theProds[iprod]);

      float pVAMCFM_wbfVV_total;
      mela.setProcess(TVar::bkgWWZZ_SMHiggs, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfVV_total, true);
      cout << "pVAMCFM_wbfVV_" << theProdName << "_total: " << pVAMCFM_wbfVV_total << '\n' << endl;
      float pVAMCFM_wbfVV_bkg;
      mela.setProcess(TVar::bkgWWZZ, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfVV_bkg, true);
      cout << "pVAMCFM_wbfVV_" << theProdName << "_bkg: " << pVAMCFM_wbfVV_bkg << '\n' << endl;
      float pVAMCFM_wbfVV_sig;
      mela.setProcess(TVar::HSMHiggs_WWZZ, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfVV_sig, true);
      cout << "pVAMCFM_wbfVV_" << theProdName << "_sig: " << pVAMCFM_wbfVV_sig << '\n' << endl;

      float pVAMCFM_wbfWW_total;
      mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfWW_total, true);
      cout << "pVAMCFM_wbfWW_" << theProdName << "_total: " << pVAMCFM_wbfWW_total << '\n' << endl;
      float pVAMCFM_wbfWW_bkg;
      mela.setProcess(TVar::bkgWW, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfWW_bkg, true);
      cout << "pVAMCFM_wbfWW_" << theProdName << "_bkg: " << pVAMCFM_wbfWW_bkg << '\n' << endl;
      float pVAMCFM_wbfWW_sig;
      mela.setProcess(TVar::HSMHiggs, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfWW_sig, true);
      cout << "pVAMCFM_wbfWW_" << theProdName << "_sig: " << pVAMCFM_wbfWW_sig << '\n' << endl;

      /***** WW (as ZZ) *****/
      cindex=1;
      mela.setCurrentCandidateFromIndex(cindex);
      float pVAMCFM_wbfWWasZZ_total;
      mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfWWasZZ_total, true);
      cout << "pVAMCFM_wbfWWasZZ_" << theProdName << "_total: " << pVAMCFM_wbfWWasZZ_total << '\n' << endl;
      float pVAMCFM_wbfWWasZZ_bkg;
      mela.setProcess(TVar::bkgZZ, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfWWasZZ_bkg, true);
      cout << "pVAMCFM_wbfWWasZZ_" << theProdName << "_bkg: " << pVAMCFM_wbfWWasZZ_bkg << '\n' << endl;
      float pVAMCFM_wbfWWasZZ_sig;
      mela.setProcess(TVar::HSMHiggs, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfWWasZZ_sig, true);
      cout << "pVAMCFM_wbfWWasZZ_" << theProdName << "_sig: " << pVAMCFM_wbfWWasZZ_sig << '\n' << endl;

      /***** ZZ *****/
      cindex=2;
      mela.setCurrentCandidateFromIndex(cindex);
      float pVAMCFM_wbfZZ_total;
      mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfZZ_total, true);
      cout << "pVAMCFM_wbfZZ_" << theProdName << "_total: " << pVAMCFM_wbfZZ_total << '\n' << endl;

      float pVAMCFM_wbfZZ_bkg;
      mela.setProcess(TVar::bkgZZ, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfZZ_bkg, true);
      cout << "pVAMCFM_wbfZZ_" << theProdName << "_bkg: " << pVAMCFM_wbfZZ_bkg << '\n' << endl;

      float pVAMCFM_wbfZZ_sig;
      mela.setProcess(TVar::HSMHiggs, TVar::MCFM, theProds[iprod]);
      mela.computeProdDecP(pVAMCFM_wbfZZ_sig, true);
      cout << "pVAMCFM_wbfZZ_" << theProdName << "_sig: " << pVAMCFM_wbfZZ_sig << '\n' << endl;
    }
    if (verbosity>=TVar::DEBUG){
      cout << "Removing Mela candidates\nSummary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }
    mela.resetInputEvent();
    cout << "Removed..." << endl;
  }

  cout.rdbuf(coutbuf);
  tout.close();

  mela.setVerbosity(bkpverbosity);
}


void testME_ProdDec_MCFM_JHUGen_Comparison(int flavor=2, bool useBkgSample=false, int motherflavor=0, int isZZWW=0 /*1==ZZ, 2==WW*/, int vbfvhchannel=1 /*1==VBF, 2==VH*/, shared_ptr<Mela> melaptr=nullptr){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.15e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;

  int idMother[2]={ 0 };
  // VBF ZZ(+)WW
  if (motherflavor==1){ idMother[0]=2; idMother[1]=1; }
  else if (motherflavor==2){ idMother[0]=-2; idMother[1]=-1; }
  // VBF ZZ-only(+)WH
  else if (motherflavor==3){ idMother[0]=2; idMother[1]=-1; }
  else if (motherflavor==4){ idMother[0]=-2; idMother[1]=1; }
  // VBF ZZ(+)ZH or WW(+)ZH
  else if (motherflavor==5){ idMother[0]=2; idMother[1]=-2; }
  else if (motherflavor==6){ idMother[0]=-2; idMother[1]=2; }
  else if (motherflavor==7){ idMother[0]=1; idMother[1]=-1; }
  else if (motherflavor==8){ idMother[0]=-1; idMother[1]=1; }

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TString coutput;
  if (useBkgSample){
    if (vbfvhchannel==1) coutput = Form("HZZ4lTree_ZZTo%s_vbfMELA_MCFMJHUGenComparison", (flavor>=2 ? "2e2mu" : "4e"));
    else coutput = Form("HZZ4lTree_ZZTo%s_vhMELA_MCFMJHUGenComparison", (flavor>=2 ? "2e2mu" : "4e"));
  }
  else{
    if (vbfvhchannel==1){
      coutput = "HZZ4lTree_VBF0P_H125.6_vbfMELA_MCFMJHUGenComparison";
    }
    else{
      coutput = "HZZ4lTree_WH125_vhMELA_MCFMJHUGenComparison";
      mPOLE=125.; wPOLE=4.07e-3;
    }
  }
  if (isZZWW==1) coutput.Append("_ZZCouplOnly");
  else if (isZZWW==2) coutput.Append("_WWCouplOnly");
  if (idMother[0]!=0 || idMother[1]!=0) coutput.Append(Form("_MotherId_%i_%i", idMother[0], idMother[1]));
  coutput.Append(".root");

  if (!melaptr) {
    melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  }
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);

  TFile* finput;
  TFile* foutput;
  if (useBkgSample){
    finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor>=2 ? "2mu2e" : "4e"), (flavor==2 ? "2e2mu" : "4e")), "read");
  }
  else{
    if (vbfvhchannel==1){
      finput = new TFile(Form("%s/%s/HZZ4lTree_VBF0P_H125.6.root", cinput_main.Data(), (flavor==2 ? "2mu2e" : "4e")), "read");
    }
    else{
      finput = new TFile(Form("%s/%s/HZZ4lTree_WH125.root", cinput_main.Data(), (flavor==2 ? "2mu2e" : "4e")), "read");
    }
  }
  foutput = new TFile(coutput, "recreate");


  float p_prod_0plus_dec_0plus_VAJHU;
  float p_prod_0plus_dec_0hplus_VAJHU;
  float p_prod_0plus_dec_0minus_VAJHU;
  float p_prod_0plus_dec_0_g1prime2_VAJHU;
  float p_prod_0hplus_dec_0hplus_VAJHU;
  float p_prod_0minus_dec_0minus_VAJHU;
  float p_prod_0_g1prime2_dec_0_g1prime2_VAJHU;
  float p_prod_0hplus_dec_0plus_VAJHU;
  float p_prod_0minus_dec_0plus_VAJHU;
  float p_prod_0_g1prime2_dec_0plus_VAJHU;

  float p_prod_0plus_dec_0plus_VAMCFM;
  float p_prod_0plus_dec_0hplus_VAMCFM;
  float p_prod_0plus_dec_0minus_VAMCFM;
  float p_prod_0plus_dec_0_g1prime2_VAMCFM;
  float p_prod_0hplus_dec_0hplus_VAMCFM;
  float p_prod_0minus_dec_0minus_VAMCFM;
  float p_prod_0_g1prime2_dec_0_g1prime2_VAMCFM;
  float p_prod_0hplus_dec_0plus_VAMCFM;
  float p_prod_0minus_dec_0plus_VAMCFM;
  float p_prod_0_g1prime2_dec_0plus_VAMCFM;


  short NJets30;
  std::vector<double>* JetPt=0;
  std::vector<double>* JetEta=0;
  std::vector<double>* JetPhi=0;
  std::vector<double>* JetMass=0;
  std::vector<double> myJetPt;
  std::vector<double> myJetEta;
  std::vector<double> myJetPhi;
  std::vector<double> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;

  float mjj = 0;
  float mjjzz = 0;
  float sysZ = 0.;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;

  float costhetastar=0;
  float costheta1=0;
  float costheta2=0;
  float Phi=0;
  float Phi1=0;
  float Q2V1=0;
  float Q2V2=0;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tree->SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
  tree->SetBranchAddress("JetMass", &JetMass, &bJetMass);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");

  newtree->Branch("p_prod_0plus_dec_0plus_VAJHU", &p_prod_0plus_dec_0plus_VAJHU);
  newtree->Branch("p_prod_0plus_dec_0hplus_VAJHU", &p_prod_0plus_dec_0hplus_VAJHU);
  newtree->Branch("p_prod_0plus_dec_0minus_VAJHU", &p_prod_0plus_dec_0minus_VAJHU);
  newtree->Branch("p_prod_0plus_dec_0_g1prime2_VAJHU", &p_prod_0plus_dec_0_g1prime2_VAJHU);
  newtree->Branch("p_prod_0hplus_dec_0hplus_VAJHU", &p_prod_0hplus_dec_0hplus_VAJHU);
  newtree->Branch("p_prod_0minus_dec_0minus_VAJHU", &p_prod_0minus_dec_0minus_VAJHU);
  newtree->Branch("p_prod_0_g1prime2_dec_0_g1prime2_VAJHU", &p_prod_0_g1prime2_dec_0_g1prime2_VAJHU);
  newtree->Branch("p_prod_0hplus_dec_0plus_VAJHU", &p_prod_0hplus_dec_0plus_VAJHU);
  newtree->Branch("p_prod_0minus_dec_0plus_VAJHU", &p_prod_0minus_dec_0plus_VAJHU);
  newtree->Branch("p_prod_0_g1prime2_dec_0plus_VAJHU", &p_prod_0_g1prime2_dec_0plus_VAJHU);

  newtree->Branch("p_prod_0plus_dec_0plus_VAMCFM", &p_prod_0plus_dec_0plus_VAMCFM);
  newtree->Branch("p_prod_0plus_dec_0hplus_VAMCFM", &p_prod_0plus_dec_0hplus_VAMCFM);
  newtree->Branch("p_prod_0plus_dec_0minus_VAMCFM", &p_prod_0plus_dec_0minus_VAMCFM);
  newtree->Branch("p_prod_0plus_dec_0_g1prime2_VAMCFM", &p_prod_0plus_dec_0_g1prime2_VAMCFM);
  newtree->Branch("p_prod_0hplus_dec_0hplus_VAMCFM", &p_prod_0hplus_dec_0hplus_VAMCFM);
  newtree->Branch("p_prod_0minus_dec_0minus_VAMCFM", &p_prod_0minus_dec_0minus_VAMCFM);
  newtree->Branch("p_prod_0_g1prime2_dec_0_g1prime2_VAMCFM", &p_prod_0_g1prime2_dec_0_g1prime2_VAMCFM);
  newtree->Branch("p_prod_0hplus_dec_0plus_VAMCFM", &p_prod_0hplus_dec_0plus_VAMCFM);
  newtree->Branch("p_prod_0minus_dec_0plus_VAMCFM", &p_prod_0minus_dec_0plus_VAMCFM);
  newtree->Branch("p_prod_0_g1prime2_dec_0plus_VAMCFM", &p_prod_0_g1prime2_dec_0plus_VAMCFM);

  newtree->Branch("ZZMass", &mzz);
  newtree->Branch("DiJetMass", &mjj);
  newtree->Branch("DiJetVVMass", &mjjzz);
  newtree->Branch("DiJetVVPz", &sysZ);
  newtree->Branch("costhetastar", &costhetastar);
  newtree->Branch("costheta1", &costheta1);
  newtree->Branch("costheta2", &costheta2);
  newtree->Branch("Phi", &Phi);
  newtree->Branch("Phi1", &Phi1);
  newtree->Branch("Q2V1", &Q2V1);
  newtree->Branch("Q2V2", &Q2V2);


  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  if (flavor == 2){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 1){
    GenLep1Id=11;
    GenLep2Id=-11;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 0){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 3){
    GenLep1Id=14;
    GenLep2Id=-14;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 4){
    GenLep1Id=0;
    GenLep2Id=-0;
    GenLep3Id=1;
    GenLep4Id=-1;
  }
  else if (flavor == 5){
    GenLep1Id=1;
    GenLep2Id=-1;
    GenLep3Id=2;
    GenLep4Id=-2;
  }
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };

  int nEntries = tree->GetEntries();
  int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=(verbosity>=TVar::DEBUG ? 4 : (!useBkgSample ? 1000 : nEntries))) break;
    //if (recorded>=1) break;
    tree->GetEntry(ev);
    if (recorded%1000==0) cout << "Nrecorded = " << recorded << "..." << endl;

    if (JetPt->size()>=2 && NJets30>=2){
      TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet1));
      associated.push_back(SimpleParticle_t(0, jet2));
      mjj = (jet1+jet2).M();

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters_ZZ;
      for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pDaughters[idau]));

      SimpleParticleCollection_t mothers;
      TLorentzVector pTotal(0, 0, 0, 0);
      for (unsigned int ip=0; ip<associated.size(); ip++) pTotal = pTotal + associated.at(ip).second;
      for (unsigned int ip=0; ip<daughters_ZZ.size(); ip++) pTotal = pTotal + daughters_ZZ.at(ip).second;
      mjjzz = pTotal.M();
      sysZ = pTotal.Z();
      TVector3 boostT(-pTotal.X()/pTotal.T(), -pTotal.Y()/pTotal.T(), 0);
      pTotal.Boost(boostT);
      TLorentzVector pMother[2];
      pMother[0].SetXYZT(0, 0, (pTotal.T()+pTotal.Z())/2., (pTotal.T()+pTotal.Z())/2.);
      pMother[1].SetXYZT(0, 0, (-pTotal.T()+pTotal.Z())/2., (pTotal.T()-pTotal.Z())/2.);
      for (unsigned int imot=0; imot<2; imot++){
        pMother[imot].Boost(-boostT);
        mothers.push_back(SimpleParticle_t(idMother[imot], pMother[imot]));
      }
      mela.setInputEvent(&daughters_ZZ, &associated, &mothers, true);

      if (mothers.size()>1){
        if (vbfvhchannel==1) TUtil::computeVBFAngles(
          costhetastar,
          costheta1,
          costheta2,
          Phi,
          Phi1,
          Q2V1,
          Q2V2,

          daughters_ZZ.at(0).second, daughters_ZZ.at(0).first,
          daughters_ZZ.at(1).second, daughters_ZZ.at(1).first,
          daughters_ZZ.at(2).second, daughters_ZZ.at(2).first,
          daughters_ZZ.at(3).second, daughters_ZZ.at(3).first,

          associated.at(0).second, associated.at(0).first,
          associated.at(1).second, associated.at(1).first,

          &(mothers.at(0).second), mothers.at(0).first,
          &(mothers.at(1).second), mothers.at(1).first
          );
        else TUtil::computeVHAngles(
          costhetastar,
          costheta1,
          costheta2,
          Phi,
          Phi1,

          daughters_ZZ.at(0).second, daughters_ZZ.at(0).first,
          daughters_ZZ.at(1).second, daughters_ZZ.at(1).first,
          daughters_ZZ.at(2).second, daughters_ZZ.at(2).first,
          daughters_ZZ.at(3).second, daughters_ZZ.at(3).first,

          associated.at(0).second, associated.at(0).first,
          associated.at(1).second, associated.at(1).first,

          &(mothers.at(0).second), mothers.at(0).first,
          &(mothers.at(1).second), mothers.at(1).first
          );
      }
      else{
        if (vbfvhchannel==1) TUtil::computeVBFAngles(
          costhetastar,
          costheta1,
          costheta2,
          Phi,
          Phi1,
          Q2V1,
          Q2V2,

          daughters_ZZ.at(0).second, daughters_ZZ.at(0).first,
          daughters_ZZ.at(1).second, daughters_ZZ.at(1).first,
          daughters_ZZ.at(2).second, daughters_ZZ.at(2).first,
          daughters_ZZ.at(3).second, daughters_ZZ.at(3).first,

          associated.at(0).second, associated.at(0).first,
          associated.at(1).second, associated.at(1).first
          );
        else TUtil::computeVHAngles(
          costhetastar,
          costheta1,
          costheta2,
          Phi,
          Phi1,

          daughters_ZZ.at(0).second, daughters_ZZ.at(0).first,
          daughters_ZZ.at(1).second, daughters_ZZ.at(1).first,
          daughters_ZZ.at(2).second, daughters_ZZ.at(2).first,
          daughters_ZZ.at(3).second, daughters_ZZ.at(3).first,

          associated.at(0).second, associated.at(0).first,
          associated.at(1).second, associated.at(1).first
          );
      }


      /***** JHUGen *****/

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZINDEPENDENT);
      float p_dec_0plus_VAJHU;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.computeP(p_dec_0plus_VAJHU, false);
      float p_dec_0hplus_VAJHU=p_dec_0plus_VAJHU;
      if (verbosity<=TVar::ERROR && isZZWW!=2){
        mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1;
        mela.computeP(p_dec_0hplus_VAJHU, false);
      }
      float p_dec_0minus_VAJHU=p_dec_0plus_VAJHU;
      if (verbosity<=TVar::ERROR && isZZWW!=2){
        mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
        mela.computeP(p_dec_0minus_VAJHU, false);
      }
      float p_dec_0_g1prime2_VAJHU=p_dec_0plus_VAJHU;
      if (verbosity<=TVar::ERROR && isZZWW!=2){
        mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-10000;
        mela.computeP(p_dec_0_g1prime2_VAJHU, false);
      }

      float p_prod_0plus_VAJHU=0;
      float p_prod_0hplus_VAJHU=0;
      float p_prod_0minus_VAJHU=0;
      float p_prod_0_g1prime2_VAJHU=0;

      if (vbfvhchannel==1){
        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][0][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdP(p_prod_0plus_VAJHU, false);
        if (verbosity<=TVar::ERROR){
          if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
          else{ mela.selfDHwwcoupl[0][1][0]=1; mela.differentiate_HWW_HZZ=true; }
          mela.computeProdP(p_prod_0hplus_VAJHU, false);

          if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
          else{ mela.selfDHwwcoupl[0][3][0]=1; mela.differentiate_HWW_HZZ=true; }
          mela.computeProdP(p_prod_0minus_VAJHU, false);

          if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-10000; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
          else{ mela.selfDHwwcoupl[0][11][0]=-10000; mela.differentiate_HWW_HZZ=true; }
          mela.computeProdP(p_prod_0_g1prime2_VAJHU, false);
        }
      }
      else{
        float p_wh_0plus_VAJHU=0;
        float p_wh_0hplus_VAJHU=0;
        float p_wh_0minus_VAJHU=0;
        float p_wh_0_g1prime2_VAJHU=0;

        if (isZZWW!=1){
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_WH);
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.computeProdP(p_wh_0plus_VAJHU, false);
          if (verbosity<=TVar::ERROR){
            mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1;
            mela.computeProdP(p_wh_0hplus_VAJHU, false);

            mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
            mela.computeProdP(p_wh_0minus_VAJHU, false);

            mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-10000;
            mela.computeProdP(p_wh_0_g1prime2_VAJHU, false);
          }
        }

        float p_zh_0plus_VAJHU=0;
        float p_zh_0hplus_VAJHU=0;
        float p_zh_0minus_VAJHU=0;
        float p_zh_0_g1prime2_VAJHU=0;

        if (isZZWW!=2){
          mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_ZH);
          mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
          mela.computeProdP(p_zh_0plus_VAJHU, false);
          if (verbosity<=TVar::ERROR){
            mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1;
            mela.computeProdP(p_zh_0hplus_VAJHU, false);

            mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
            mela.computeProdP(p_zh_0minus_VAJHU, false);

            mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-10000;
            mela.computeProdP(p_zh_0_g1prime2_VAJHU, false);
          }
        }

        p_prod_0plus_VAJHU = p_wh_0plus_VAJHU + p_zh_0plus_VAJHU;
        p_prod_0hplus_VAJHU = p_wh_0hplus_VAJHU + p_zh_0hplus_VAJHU;
        p_prod_0minus_VAJHU = p_wh_0minus_VAJHU + p_zh_0minus_VAJHU;
        p_prod_0_g1prime2_VAJHU = p_wh_0_g1prime2_VAJHU + p_zh_0_g1prime2_VAJHU;
      }

      p_prod_0plus_dec_0plus_VAJHU=p_prod_0plus_VAJHU*p_dec_0plus_VAJHU;
      p_prod_0plus_dec_0hplus_VAJHU=p_prod_0plus_VAJHU*p_dec_0hplus_VAJHU;
      p_prod_0plus_dec_0minus_VAJHU=p_prod_0plus_VAJHU*p_dec_0minus_VAJHU;
      p_prod_0plus_dec_0_g1prime2_VAJHU=p_prod_0plus_VAJHU*p_dec_0_g1prime2_VAJHU;
      p_prod_0hplus_dec_0hplus_VAJHU=p_prod_0hplus_VAJHU*p_dec_0hplus_VAJHU;
      p_prod_0minus_dec_0minus_VAJHU=p_prod_0minus_VAJHU*p_dec_0minus_VAJHU;
      p_prod_0_g1prime2_dec_0_g1prime2_VAJHU=p_prod_0_g1prime2_VAJHU*p_dec_0_g1prime2_VAJHU;
      p_prod_0hplus_dec_0plus_VAJHU=p_prod_0hplus_VAJHU*p_dec_0plus_VAJHU;
      p_prod_0minus_dec_0plus_VAJHU=p_prod_0minus_VAJHU*p_dec_0plus_VAJHU;
      p_prod_0_g1prime2_dec_0plus_VAJHU=p_prod_0_g1prime2_VAJHU*p_dec_0plus_VAJHU;

      double propagator = 1./(pow(pow(mzz, 2)-pow(mPOLE, 2), 2)+pow(mPOLE*wPOLE, 2));
      if (vbfvhchannel!=1){ // JHUGen VH pseudo-propagator
        double mh, gah;
        mela.getIORecord()->getHiggsMassWidth(mh, gah, 0);
        propagator /= 1./(pow(pow(mzz, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
      }

      /***** MCFM *****/

      mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF);

      spinzerohiggs_anomcoupl_.channeltoggle_stu=0;
      spinzerohiggs_anomcoupl_.vvhvvtoggle_vbfvh=(vbfvhchannel==1 ? 0 : 1);

      spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
      spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
      if (isZZWW==2) spinzerohiggs_anomcoupl_.AnomalCouplDK=0;
      if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
      else{ mela.selfDHwwcoupl[0][0][0]=1; mela.differentiate_HWW_HZZ=true; }
      mela.computeProdDecP(p_prod_0plus_dec_0plus_VAMCFM, false);
      /*
      if (verbosity<=TVar::ERROR){
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdDecP(p_prod_0hplus_dec_0hplus_VAMCFM, false);
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][3][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdDecP(p_prod_0minus_dec_0minus_VAMCFM, false);
        mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-10000;
        mela.computeProdDecP(p_prod_0_g1prime2_dec_0_g1prime2_VAMCFM, false);

        spinzerohiggs_anomcoupl_.AnomalCouplPR=0;
        spinzerohiggs_anomcoupl_.AnomalCouplDK=1;
        if (isZZWW==2) spinzerohiggs_anomcoupl_.AnomalCouplDK=0;
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdDecP(p_prod_0plus_dec_0hplus_VAMCFM, false);
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][3][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdDecP(p_prod_0plus_dec_0minus_VAMCFM, false);
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-10000; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][11][0]=-10000; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdDecP(p_prod_0plus_dec_0_g1prime2_VAMCFM, false);

        spinzerohiggs_anomcoupl_.AnomalCouplPR=1;
        spinzerohiggs_anomcoupl_.AnomalCouplDK=0;
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][1][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdDecP(p_prod_0hplus_dec_0plus_VAMCFM, false);
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][3][0]=1; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdDecP(p_prod_0minus_dec_0plus_VAMCFM, false);
        if (isZZWW!=2){ mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-10000; if (isZZWW==1) mela.differentiate_HWW_HZZ=true; }
        else{ mela.selfDHwwcoupl[0][11][0]=-10000; mela.differentiate_HWW_HZZ=true; }
        mela.computeProdDecP(p_prod_0_g1prime2_dec_0plus_VAMCFM, false);
      }
      */

      if (verbosity>=TVar::DEBUG) TUtil::PrintCandidateSummary(mela.getCurrentCandidate());

      mela.resetInputEvent();

      if (p_prod_0plus_dec_0plus_VAJHU>0.){
        p_prod_0plus_dec_0hplus_VAJHU /= p_prod_0plus_dec_0plus_VAJHU;
        p_prod_0plus_dec_0minus_VAJHU /= p_prod_0plus_dec_0plus_VAJHU;
        p_prod_0plus_dec_0_g1prime2_VAJHU /= p_prod_0plus_dec_0plus_VAJHU;
        p_prod_0hplus_dec_0hplus_VAJHU /= p_prod_0plus_dec_0plus_VAJHU;
        p_prod_0minus_dec_0minus_VAJHU /= p_prod_0plus_dec_0plus_VAJHU;
        p_prod_0_g1prime2_dec_0_g1prime2_VAJHU /= p_prod_0plus_dec_0plus_VAJHU;
        p_prod_0hplus_dec_0plus_VAJHU /= p_prod_0plus_dec_0plus_VAJHU;
        p_prod_0minus_dec_0plus_VAJHU /= p_prod_0plus_dec_0plus_VAJHU;
        p_prod_0_g1prime2_dec_0plus_VAJHU /= p_prod_0plus_dec_0plus_VAJHU;
        p_prod_0plus_dec_0plus_VAJHU *= propagator;
      }
      if (p_prod_0plus_dec_0plus_VAMCFM>0.){
        p_prod_0plus_dec_0hplus_VAMCFM /= p_prod_0plus_dec_0plus_VAMCFM;
        p_prod_0plus_dec_0minus_VAMCFM /= p_prod_0plus_dec_0plus_VAMCFM;
        p_prod_0plus_dec_0_g1prime2_VAMCFM /= p_prod_0plus_dec_0plus_VAMCFM;
        p_prod_0hplus_dec_0hplus_VAMCFM /= p_prod_0plus_dec_0plus_VAMCFM;
        p_prod_0minus_dec_0minus_VAMCFM /= p_prod_0plus_dec_0plus_VAMCFM;
        p_prod_0_g1prime2_dec_0_g1prime2_VAMCFM /= p_prod_0plus_dec_0plus_VAMCFM;
        p_prod_0hplus_dec_0plus_VAMCFM /= p_prod_0plus_dec_0plus_VAMCFM;
        p_prod_0minus_dec_0plus_VAMCFM /= p_prod_0plus_dec_0plus_VAMCFM;
        p_prod_0_g1prime2_dec_0plus_VAMCFM /= p_prod_0plus_dec_0plus_VAMCFM;
      }

      newtree->Fill();
      recorded++;
    }
  }
  foutput->WriteTObject(newtree);
  delete newtree;
  foutput->Close();
  finput->Close();
  mela.setVerbosity(bkpverbosity);
}

void testME_Dec_JHUGenMCFM_Ping(int flavor=2, shared_ptr<Mela> melaptr=nullptr){
  ofstream tout(TString("testME_Dec_JHUGenMCFM_Ping_")+long(flavor)+".out");
  streambuf* coutbuf = cout.rdbuf();
  cout.rdbuf(tout.rdbuf());

  int erg_tev=13;
  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  if (!melaptr) melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);

  double aL1, aL2, aR1, aR2;

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  /*
  const int nEntries = 3;
  double l1_array[nEntries][4] = {
  {1365.4973807340846,        10.289826593755228,        25.205694382277809,       -1365.2259480507332},
  {238.65751023078761,        9.2808858562825005,        15.827726043466324,       -237.95116187061188},
  {101.52463181523598,        27.359569630718468,      -0.90299073100241323,       -97.764458892691749}
  };
  double l2_array[nEntries][4] = {
  {22.786181013986834,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371},
  {101.67043553688544,        2.1169375132239789,       0.77953005873937187,       -101.64540506443268},
  {24.717634703436786,       -1.1722249478288802,       -5.9599387484197646,       -23.959684558009428}
  };
  double l3_array[nEntries][4] = {
  {1895.7562628816693,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620},
  {317.81904277258536,        2.5882005498984775,        21.352807448987718,       -317.09037005377883},
  {180.10885677707822,       -6.7240759244122792,        35.742176497019194,       -176.39865053838915}
  };
  double l4_array[nEntries][4] = {
  {471.71918486784784,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354},
  {95.655512770627581,       -13.986023919404957,       -37.960063551193414,       -86.679881365440792},
  {49.137252081251319,       -19.463268758477309,       -28.879247017597017,       -34.664676589120688}
  };
  */
  const int nEntries = 6;
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  for (int ev = 0; ev < 1; ev++){
    if (flavor == 2){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 1){
      GenLep1Id=11;
      GenLep2Id=-11;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 0){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 3){
      GenLep1Id=14;
      GenLep2Id=-14;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 4){
      GenLep1Id=0;
      GenLep2Id=-0;
      GenLep3Id=1;
      GenLep4Id=-1;
    }
    else if (flavor == 5){
      GenLep1Id=1;
      GenLep2Id=-1;
      GenLep3Id=2;
      GenLep4Id=-2;
    }

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    int idOrdered_WW[4] ={ 11, -12, -11, 12 };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters_WW;
    for (unsigned int idau=0; idau<4; idau++) daughters_WW.push_back(SimpleParticle_t(idOrdered_WW[idau], pOrdered[idau]));
    SimpleParticleCollection_t daughters_WWasZZ;
    for (unsigned int iv=0; iv<2; iv++){ for (int ivd=0; ivd<2; ivd++) daughters_WWasZZ.push_back(SimpleParticle_t(idOrdered_WW[iv+2*ivd], pOrdered[iv+2*ivd])); }
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    TLorentzVector pOrdered_ZG[3];
    pOrdered_ZG[0]=pOrdered[0];
    pOrdered_ZG[1]=pOrdered[1];
    pOrdered_ZG[2]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_ZG;
    for (unsigned int idau=0; idau<2; idau++) daughters_ZG.push_back(SimpleParticle_t(idOrdered[idau], pOrdered_ZG[idau]));
    for (unsigned int idau=2; idau<3; idau++) daughters_ZG.push_back(SimpleParticle_t(22, pOrdered_ZG[idau]));

    TLorentzVector pOrdered_GG[3];
    pOrdered_GG[0]=pOrdered[0]+pOrdered[1];
    pOrdered_GG[1]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_GG;
    for (unsigned int idau=0; idau<2; idau++) daughters_GG.push_back(SimpleParticle_t(22, pOrdered_GG[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela.setInputEvent(&daughters_WW, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setInputEvent(&daughters_ZG, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_GG);
    mela.setInputEvent(&daughters_GG, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

    int cindex;

    /***** WW *****/
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    // JHUGen ME
    float pVAJHUGen_ggWW_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggWW_SM_sig, true);
    cout << "pVAJHUGen_ggWW_SM_sig: " << pVAJHUGen_ggWW_SM_sig << '\n' << endl;

    float pVAJHUGen_ggWW_0pm_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.;
    mela.computeP(pVAJHUGen_ggWW_0pm_sig, true);
    cout << "pVAJHUGen_ggWW_0pm_sig: " << pVAJHUGen_ggWW_0pm_sig << '\n' << endl;
    mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
    mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);

    float pVAJHUGen_ggWW_0pm_CT_sig;
    //mela.setVerbosity(TVar::DEBUG_VERBOSE);
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzpzpcoupl[gHIGGS_VV_1][0]=1.;
    //mela.selfDHwpwpcoupl[gHIGGS_VV_1][0]=1.;
    mela.selfDWpffcoupl[gHIGGS_Vp_El_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Mu_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Ta_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Up_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Chm_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Top_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_El_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Mu_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Ta_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Up_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Chm_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Top_right][0]=aR1;
    mela.selfDM_Wprime=80.399;
    mela.selfDGa_Wprime=2.085;
    mela.computeP(pVAJHUGen_ggWW_0pm_CT_sig, true);
    cout << "pVAJHUGen_ggWW_0pm_CT_sig: " << pVAJHUGen_ggWW_0pm_CT_sig << '\n' << endl;
    //mela.setVerbosity(verbosity);

    float pVAJHUGen_ggWW_PS_sig;
    mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggWW_PS_sig, true);
    cout << "pVAJHUGen_ggWW_PS_sig: " << pVAJHUGen_ggWW_PS_sig << '\n' << endl;

    float pVAJHUGen_ggWW_0m_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.;
    mela.computeP(pVAJHUGen_ggWW_0m_sig, true);
    cout << "pVAJHUGen_ggWW_0m_sig: " << pVAJHUGen_ggWW_0m_sig << '\n' << endl;

    // MCFM ME
    float pVAMCFM_ggWW_SM_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_SM_total, true);
    cout << "pVAMCFM_ggWW_SM_total: " << pVAMCFM_ggWW_SM_total << '\n' << endl;

    float pVAMCFM_ggWW_0pm_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.;
    mela.computeP(pVAMCFM_ggWW_0pm_total, true);
    cout << "pVAMCFM_ggWW_0pm_total: " << pVAMCFM_ggWW_0pm_total << '\n' << endl;

    float pVAMCFM_ggWW_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggWW_SM_sig, true);
    cout << "pVAMCFM_ggWW_SM_sig: " << pVAMCFM_ggWW_SM_sig << '\n' << endl;

    float pVAMCFM_ggWW_0pm_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.;
    mela.computeP(pVAMCFM_ggWW_0pm_sig, true);
    cout << "pVAMCFM_ggWW_0pm_sig: " << pVAMCFM_ggWW_0pm_sig << '\n' << endl;

    float pVAMCFM_ggWW_0m_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.;
    mela.computeP(pVAMCFM_ggWW_0m_sig, true);
    cout << "pVAMCFM_ggWW_0m_sig: " << pVAMCFM_ggWW_0m_sig << '\n' << endl;

    // MCFM HM ME
    float pVAMCFM_HM_ggWW_SM_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.computeP(pVAMCFM_HM_ggWW_SM_total, true);
    cout << "pVAMCFM_HM_ggWW_SM_total: " << pVAMCFM_HM_ggWW_SM_total << '\n' << endl;

    float pVAMCFM_HM_ggWW_0pm_total;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[1][0][0]=1.;
    mela.computeP(pVAMCFM_HM_ggWW_0pm_total, true);
    cout << "pVAMCFM_HM_ggWW_0pm_total: " << pVAMCFM_HM_ggWW_0pm_total << '\n' << endl;

    float pVAMCFM_HM_ggWW_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.computeP(pVAMCFM_HM_ggWW_SM_sig, true);
    cout << "pVAMCFM_HM_ggWW_SM_sig: " << pVAMCFM_HM_ggWW_SM_sig << '\n' << endl;

    float pVAMCFM_HM_ggWW_0pm_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[1][0][0]=1.;
    mela.computeP(pVAMCFM_HM_ggWW_0pm_sig, true);
    cout << "pVAMCFM_HM_ggWW_0pm_sig: " << pVAMCFM_HM_ggWW_0pm_sig << '\n' << endl;

    float pVAMCFM_HM_ggWW_0m_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[1][3][0]=1.;
    mela.computeP(pVAMCFM_HM_ggWW_0m_sig, true);
    cout << "pVAMCFM_HM_ggWW_0m_sig: " << pVAMCFM_HM_ggWW_0m_sig << '\n' << endl;

    /***** ZZ *****/
    cindex=1;
    mela.setCurrentCandidateFromIndex(cindex);

    // JHUGen ME
    float pVAJHUGen_ggZZ_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggZZ_SM_sig, true);
    cout << "pVAJHUGen_ggZZ_SM_sig: " << pVAJHUGen_ggZZ_SM_sig << '\n' << endl;

    float pVAJHUGen_ggZZ_0pm_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.;
    mela.computeP(pVAJHUGen_ggZZ_0pm_sig, true);
    cout << "pVAJHUGen_ggZZ_0pm_sig: " << pVAJHUGen_ggZZ_0pm_sig << '\n' << endl;
    mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
    mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);

    float pVAJHUGen_ggZZ_0pm_CT_sig;
    //mela.setVerbosity(TVar::DEBUG_VERBOSE);
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzpzpcoupl[gHIGGS_VV_1][0]=1.;
    mela.selfDZpffcoupl[gHIGGS_Vp_El_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Mu_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Ta_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_NuE_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Up_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Chm_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Top_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Dn_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Str_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Bot_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_El_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Mu_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Ta_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_NuE_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Up_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Chm_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Top_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Dn_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Str_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Bot_right][0]=aR1;
    mela.selfDM_Zprime=91.1876;
    mela.selfDGa_Zprime=2.4952;
    mela.computeP(pVAJHUGen_ggZZ_0pm_CT_sig, true);
    cout << "pVAJHUGen_ggZZ_0pm_CT_sig: " << pVAJHUGen_ggZZ_0pm_CT_sig << '\n' << endl;
    //mela.setVerbosity(verbosity);

    float pVAJHUGen_ggZZ_PS_sig;
    mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggZZ_PS_sig, true);
    cout << "pVAJHUGen_ggZZ_PS_sig: " << pVAJHUGen_ggZZ_PS_sig << '\n' << endl;

    float pVAJHUGen_ggZZ_0m_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.;
    mela.computeP(pVAJHUGen_ggZZ_0m_sig, true);
    cout << "pVAJHUGen_ggZZ_0m_sig: " << pVAJHUGen_ggZZ_0m_sig << '\n' << endl;

    // MCFM ME
    float pVAMCFM_ggZZ_SM_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_SM_total, true);
    cout << "pVAMCFM_ggZZ_SM_total: " << pVAMCFM_ggZZ_SM_total << '\n' << endl;

    float pVAMCFM_ggZZ_0pm_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.;
    mela.computeP(pVAMCFM_ggZZ_0pm_total, true);
    cout << "pVAMCFM_ggZZ_0pm_total: " << pVAMCFM_ggZZ_0pm_total << '\n' << endl;

    float pVAMCFM_ggZZ_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(pVAMCFM_ggZZ_SM_sig, true);
    cout << "pVAMCFM_ggZZ_SM_sig: " << pVAMCFM_ggZZ_SM_sig << '\n' << endl;

    float pVAMCFM_ggZZ_0pm_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1.;
    mela.computeP(pVAMCFM_ggZZ_0pm_sig, true);
    cout << "pVAMCFM_ggZZ_0pm_sig: " << pVAMCFM_ggZZ_0pm_sig << '\n' << endl;

    float pVAMCFM_ggZZ_0m_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHttcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[0][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1.;
    mela.computeP(pVAMCFM_ggZZ_0m_sig, true);
    cout << "pVAMCFM_ggZZ_0m_sig: " << pVAMCFM_ggZZ_0m_sig << '\n' << endl;

    // MCFM HM ME
    float pVAMCFM_HM_ggZZ_SM_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.computeP(pVAMCFM_HM_ggZZ_SM_total, true);
    cout << "pVAMCFM_HM_ggZZ_SM_total: " << pVAMCFM_HM_ggZZ_SM_total << '\n' << endl;

    float pVAMCFM_HM_ggZZ_0pm_total;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[1][0][0]=1.;
    mela.computeP(pVAMCFM_HM_ggZZ_0pm_total, true);
    cout << "pVAMCFM_HM_ggZZ_0pm_total: " << pVAMCFM_HM_ggZZ_0pm_total << '\n' << endl;

    float pVAMCFM_HM_ggZZ_SM_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.computeP(pVAMCFM_HM_ggZZ_SM_sig, true);
    cout << "pVAMCFM_HM_ggZZ_SM_sig: " << pVAMCFM_HM_ggZZ_SM_sig << '\n' << endl;

    float pVAMCFM_HM_ggZZ_0pm_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[1][0][0]=1.;
    mela.computeP(pVAMCFM_HM_ggZZ_0pm_sig, true);
    cout << "pVAMCFM_HM_ggZZ_0pm_sig: " << pVAMCFM_HM_ggZZ_0pm_sig << '\n' << endl;

    float pVAMCFM_HM_ggZZ_0m_sig;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.setMelaHiggsMassWidth(-1, 0, 0); mela.setMelaHiggsMassWidth(mPOLE, -1, 1);
    mela.selfDHttcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHbbcoupl[1][gHIGGS_KAPPA][0]=1.;
    mela.selfDHzzcoupl[1][3][0]=1.;
    mela.computeP(pVAMCFM_HM_ggZZ_0m_sig, true);
    cout << "pVAMCFM_HM_ggZZ_0m_sig: " << pVAMCFM_HM_ggZZ_0m_sig << '\n' << endl;

    /***** ZG *****/
    cindex=2;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAJHUGen_ggZG_SM_sig;
    mela.setProcess(TVar::H0_Zgs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggZG_SM_sig, true);
    cout << "pVAJHUGen_ggZG_SM_sig: " << pVAJHUGen_ggZG_SM_sig << '\n' << endl;

    float pVAJHUGen_ggZG_0pm_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_2][0]=1.;
    mela.computeP(pVAJHUGen_ggZG_0pm_sig, true);
    cout << "pVAJHUGen_ggZG_0pm_sig: " << pVAJHUGen_ggZG_0pm_sig << '\n' << endl;

    float pVAJHUGen_ggZG_PS_sig;
    mela.setProcess(TVar::H0_Zgs_PS, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggZG_PS_sig, true);
    cout << "pVAJHUGen_ggZG_PS_sig: " << pVAJHUGen_ggZG_PS_sig << '\n' << endl;

    float pVAJHUGen_ggZG_0m_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_ZA_4][0]=1.;
    mela.computeP(pVAJHUGen_ggZG_0m_sig, true);
    cout << "pVAJHUGen_ggZG_0m_sig: " << pVAJHUGen_ggZG_0m_sig << '\n' << endl;

    /***** GG *****/
    cindex=3;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAJHUGen_ggGG_SM_sig;
    mela.setProcess(TVar::H0_gsgs, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggGG_SM_sig, true);
    cout << "pVAJHUGen_ggGG_SM_sig: " << pVAJHUGen_ggGG_SM_sig << '\n' << endl;

    float pVAJHUGen_ggGG_0pm_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_AA_2][0]=1.;
    mela.computeP(pVAJHUGen_ggGG_0pm_sig, true);
    cout << "pVAJHUGen_ggGG_0pm_sig: " << pVAJHUGen_ggGG_0pm_sig << '\n' << endl;

    float pVAJHUGen_ggGG_PS_sig;
    mela.setProcess(TVar::H0_gsgs_PS, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(pVAJHUGen_ggGG_PS_sig, true);
    cout << "pVAJHUGen_ggGG_PS_sig: " << pVAJHUGen_ggGG_PS_sig << '\n' << endl;

    float pVAJHUGen_ggGG_0m_sig;
    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1.;
    mela.selfDHzzcoupl[0][gHIGGS_AA_4][0]=1.;
    mela.computeP(pVAJHUGen_ggGG_0m_sig, true);
    cout << "pVAJHUGen_ggGG_0m_sig: " << pVAJHUGen_ggGG_0m_sig << '\n' << endl;


    if (verbosity>=TVar::DEBUG){
      cout << "Removing Mela candidates\nSummary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }
    mela.resetInputEvent();
    cout << "Removed..." << endl;
  }

  cout.rdbuf(coutbuf);
  tout.close();
  mela.setVerbosity(bkpverbosity);
}

void testME_ProdDec_JHUGen_SpinTwo_Ping(int flavor=2, shared_ptr<Mela> melaptr=nullptr){
  ofstream tout(TString("testME_ProdDec_JHUGen_SpinTwo_Ping_")+long(flavor)+".out");
  streambuf* coutbuf = cout.rdbuf();
  cout.rdbuf(tout.rdbuf());

  int erg_tev=13;
  float mPOLE=125.;
  float wPOLE=4.07e-3;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  if (!melaptr) melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);

  double aL1, aL2, aR1, aR2;

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  /*
  const int nEntries = 3;
  double l1_array[nEntries][4] = {
  {1365.4973807340846,        10.289826593755228,        25.205694382277809,       -1365.2259480507332},
  {238.65751023078761,        9.2808858562825005,        15.827726043466324,       -237.95116187061188},
  {101.52463181523598,        27.359569630718468,      -0.90299073100241323,       -97.764458892691749}
  };
  double l2_array[nEntries][4] = {
  {22.786181013986834,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371},
  {101.67043553688544,        2.1169375132239789,       0.77953005873937187,       -101.64540506443268},
  {24.717634703436786,       -1.1722249478288802,       -5.9599387484197646,       -23.959684558009428}
  };
  double l3_array[nEntries][4] = {
  {1895.7562628816693,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620},
  {317.81904277258536,        2.5882005498984775,        21.352807448987718,       -317.09037005377883},
  {180.10885677707822,       -6.7240759244122792,        35.742176497019194,       -176.39865053838915}
  };
  double l4_array[nEntries][4] = {
  {471.71918486784784,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354},
  {95.655512770627581,       -13.986023919404957,       -37.960063551193414,       -86.679881365440792},
  {49.137252081251319,       -19.463268758477309,       -28.879247017597017,       -34.664676589120688}
  };
  */
  const int nEntries = 6;
  double l1_array[nEntries][4] ={
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 51.374202, 25.924766, 12.290178, 42.616376 },
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973807340846, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Massless
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 }, // Muon via E
    { 1365.4973848483, 10.289826593755228, 25.205694382277809, -1365.2259480507332 } // Muon via E
  };
  double l2_array[nEntries][4] ={
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 }
  };
  double l3_array[nEntries][4] ={
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 75.823478, -16.640412, 23.246999, 70.227220 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562628816693, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 },
    { 1895.7562658451, 25.837804322120054, -28.821887970086259, -1895.3610513294620 }
  };
  double l4_array[nEntries][4] ={
    { 21.481452, 9.489680, -9.336587, 16.858699 },
    { 271.875752, 70.427173, -11.138146, 261.769598 },
    { 471.71918486784784, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.786181013986834, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 },
    { 471.7191967775, -35.976267906053060, 4.5169691019519895, -470.32360615864354 },
    { 22.7864275656, -0.15136300982222117, -0.90077551414353962, -22.767866345236371 }
  };

  for (int ev = 0; ev < 1; ev++){
    if (flavor == 2){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 1){
      GenLep1Id=11;
      GenLep2Id=-11;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else if (flavor == 0){
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 3){
      GenLep1Id=14;
      GenLep2Id=-14;
      GenLep3Id=13;
      GenLep4Id=-13;
    }
    else if (flavor == 4){
      GenLep1Id=0;
      GenLep2Id=-0;
      GenLep3Id=1;
      GenLep4Id=-1;
    }
    else if (flavor == 5){
      GenLep1Id=1;
      GenLep2Id=-1;
      GenLep3Id=2;
      GenLep4Id=-2;
    }

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    int idOrdered_WW[4] ={ 11, -12, -11, 12 };
    TLorentzVector pOrdered[4];
    pOrdered[0].SetXYZT(l1_array[ev][1], l1_array[ev][2], l1_array[ev][3], l1_array[ev][0]);
    pOrdered[1].SetXYZT(l2_array[ev][1], l2_array[ev][2], l2_array[ev][3], l2_array[ev][0]);
    pOrdered[2].SetXYZT(l3_array[ev][1], l3_array[ev][2], l3_array[ev][3], l3_array[ev][0]);
    pOrdered[3].SetXYZT(l4_array[ev][1], l4_array[ev][2], l4_array[ev][3], l4_array[ev][0]);
    SimpleParticleCollection_t daughters_WW;
    for (unsigned int idau=0; idau<4; idau++) daughters_WW.push_back(SimpleParticle_t(idOrdered_WW[idau], pOrdered[idau]));
    SimpleParticleCollection_t daughters_WWasZZ;
    for (unsigned int iv=0; iv<2; iv++){ for (int ivd=0; ivd<2; ivd++) daughters_WWasZZ.push_back(SimpleParticle_t(idOrdered_WW[iv+2*ivd], pOrdered[iv+2*ivd])); }
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));

    TLorentzVector pOrdered_ZG[3];
    pOrdered_ZG[0]=pOrdered[0];
    pOrdered_ZG[1]=pOrdered[1];
    pOrdered_ZG[2]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_ZG;
    for (unsigned int idau=0; idau<2; idau++) daughters_ZG.push_back(SimpleParticle_t(idOrdered[idau], pOrdered_ZG[idau]));
    for (unsigned int idau=2; idau<3; idau++) daughters_ZG.push_back(SimpleParticle_t(22, pOrdered_ZG[idau]));

    TLorentzVector pOrdered_GG[3];
    pOrdered_GG[0]=pOrdered[0]+pOrdered[1];
    pOrdered_GG[1]=pOrdered[2]+pOrdered[3];
    SimpleParticleCollection_t daughters_GG;
    for (unsigned int idau=0; idau<2; idau++) daughters_GG.push_back(SimpleParticle_t(22, pOrdered_GG[idau]));

    mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela.setInputEvent(&daughters_WW, (SimpleParticleCollection_t*) 0, (SimpleParticleCollection_t*) 0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*) 0, (SimpleParticleCollection_t*) 0, false);
    mela.setInputEvent(&daughters_ZG, (SimpleParticleCollection_t*) 0, (SimpleParticleCollection_t*) 0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_GG);
    mela.setInputEvent(&daughters_GG, (SimpleParticleCollection_t*) 0, (SimpleParticleCollection_t*) 0, false);

    int cindex;

    /***** WW *****/
    cindex=0;
    mela.setCurrentCandidateFromIndex(cindex);

    // JHUGen ME
    float pVAJHUGen_ggWW_a1_b1_sig;
    //mela.setVerbosity(TVar::DEBUG_VERBOSE);
    mela.setProcess(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG);
    mela.selfDGggcoupl[gGRAVITON_GG_1][0]=1.;
    mela.selfDGvvcoupl[gGRAVITON_VV_1][0]=1.;
    mela.computeP(pVAJHUGen_ggWW_a1_b1_sig, true);
    cout << "pVAJHUGen_ggWW_a1_b1_sig: " << pVAJHUGen_ggWW_a1_b1_sig << '\n' << endl;
    //mela.setVerbosity(verbosity);
    mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
    mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
    //cout << "aL/R 1, 2 = " << aL1 << ", " << aR1 << ", " << aL2 << ", " << aR2 << endl;

    float pVAJHUGen_ggWW_a1_b1_CT_sig;
    //mela.setVerbosity(TVar::DEBUG_VERBOSE);
    mela.setProcess(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG);
    mela.selfDGggcoupl[gGRAVITON_GG_1][0]=1.;
    mela.selfDGvpvpcoupl[gGRAVITON_VV_1][0]=1.;
    mela.selfDWpffcoupl[gHIGGS_Vp_El_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Mu_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Ta_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Up_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Chm_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Top_left][0]=aL1;
    mela.selfDWpffcoupl[gHIGGS_Vp_El_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Mu_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Ta_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Up_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Chm_right][0]=aR1;
    mela.selfDWpffcoupl[gHIGGS_Vp_Top_right][0]=aR1;
    mela.selfDM_Wprime=80.399;
    mela.selfDGa_Wprime=2.085;
    mela.computeP(pVAJHUGen_ggWW_a1_b1_CT_sig, true);
    cout << "pVAJHUGen_ggWW_a1_b1_CT_sig: " << pVAJHUGen_ggWW_a1_b1_CT_sig << '\n' << endl;
    //mela.setVerbosity(verbosity);


    /***** ZZ *****/
    cindex=1;
    mela.setCurrentCandidateFromIndex(cindex);

    // JHUGen ME
    float pVAJHUGen_ggZZ_a1_b1_sig;
    mela.setProcess(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG);
    mela.selfDGggcoupl[gGRAVITON_GG_1][0]=1.;
    mela.selfDGvvcoupl[gGRAVITON_VV_1][0]=1.;
    mela.computeP(pVAJHUGen_ggZZ_a1_b1_sig, true);
    cout << "pVAJHUGen_ggZZ_a1_b1_sig: " << pVAJHUGen_ggZZ_a1_b1_sig << '\n' << endl;
    mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
    mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
    //cout << "aL/R 1, 2 = " << aL1 << ", " << aR1 << ", " << aL2 << ", " << aR2 << endl;

    float pVAJHUGen_ggZZ_a1_b1_CT_sig;
    mela.setProcess(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG);
    mela.selfDGggcoupl[gGRAVITON_GG_1][0]=1.;
    mela.selfDGvpvpcoupl[gGRAVITON_VV_1][0]=1.;
    mela.selfDZpffcoupl[gHIGGS_Vp_El_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Mu_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Ta_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_NuE_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Up_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Chm_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Top_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Dn_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Str_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Bot_left][0]=aL1;
    mela.selfDZpffcoupl[gHIGGS_Vp_El_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Mu_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Ta_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_NuE_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Up_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Chm_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Top_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Dn_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Str_right][0]=aR1;
    mela.selfDZpffcoupl[gHIGGS_Vp_Bot_right][0]=aR1;
    mela.selfDM_Zprime=91.1876;
    mela.selfDGa_Zprime=2.4952;
    mela.computeP(pVAJHUGen_ggZZ_a1_b1_CT_sig, true);
    cout << "pVAJHUGen_ggZZ_a1_b1_CT_sig: " << pVAJHUGen_ggZZ_a1_b1_CT_sig << '\n' << endl;

    /***** ZG *****/
    cindex=2;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAJHUGen_ggZG_a1_b1_sig;
    mela.setProcess(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG);
    mela.selfDGggcoupl[gGRAVITON_GG_1][0]=1.;
    mela.selfDGvvcoupl[gGRAVITON_ZA_1][0]=1.;
    mela.computeP(pVAJHUGen_ggZG_a1_b1_sig, true);
    cout << "pVAJHUGen_ggZG_a1_b1_sig: " << pVAJHUGen_ggZG_a1_b1_sig << '\n' << endl;

    /***** GG *****/
    cindex=3;
    mela.setCurrentCandidateFromIndex(cindex);
    float pVAJHUGen_ggGG_a1_b1_sig;
    mela.setProcess(TVar::SelfDefine_spin2, TVar::JHUGen, TVar::ZZGG);
    mela.selfDGggcoupl[gGRAVITON_GG_1][0]=1.;
    mela.selfDGvvcoupl[gGRAVITON_AA_1][0]=1.;
    mela.computeP(pVAJHUGen_ggGG_a1_b1_sig, true);
    cout << "pVAJHUGen_ggGG_a1_b1_sig: " << pVAJHUGen_ggGG_a1_b1_sig << '\n' << endl;

    if (verbosity>=TVar::DEBUG){
      cout << "Removing Mela candidates\nSummary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }
    mela.resetInputEvent();
    cout << "Removed..." << endl;
  }

  cout.rdbuf(coutbuf);
  tout.close();
  mela.setVerbosity(bkpverbosity);
}

void testME_Dec_FullSim(int flavor=2, bool useConstants=false, bool useBkgSample=false, shared_ptr<Mela> melaptr=nullptr){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  if (!melaptr) {
    melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  }
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput;
  TFile* foutput;
  if (!useBkgSample){
    finput = new TFile(Form("%s/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root", cinput_main.Data(), (flavor>=2 ? "2mu2e" : "4e")), "read");
    foutput = new TFile(Form("HZZ4lTree_powheg15jhuGenV3-0PMH125.6_%s_OriginalMEv2ValidationTestOnly.root", (flavor>=2 ? "2mu2e" : "4e")), "recreate");
  }
  else{
    finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor>=2 ? "2mu2e" : "4e"), (flavor==2 ? "2e2mu" : "4e")), "read");
    foutput = new TFile(Form("HZZ4lTree_ZZTo%s_OriginalMEv2ValidationTestOnly.root", (flavor>=2 ? "2e2mu" : "4e")), "recreate");
  }

  float p0plus_VAJHU;
  float p0hplus_VAJHU;
  float p0minus_VAJHU;
  float p0_g1prime2_VAJHU;
  float pg1g1prime2_VAJHU;
  //float pg1g1prime2_pi2_VAJHU;
  float pg1g2_VAJHU;
  float pg1g2_pi2_VAJHU;
  float pg1g4_VAJHU;
  float pg1g4_pi2_VAJHU;

  float p0plus_VAJHU_NEW;
  float p0hplus_VAJHU_NEW;
  float p0minus_VAJHU_NEW;
  float p0_g1prime2_VAJHU_NEW;
  float pg1g1prime2_VAJHU_NEW;
  float pg1g1prime2_pi2_VAJHU_NEW;
  float pg1g2_VAJHU_NEW;
  float pg1g2_pi2_VAJHU_NEW;
  float pg1g4_VAJHU_NEW;
  float pg1g4_pi2_VAJHU_NEW;

  float p0hplus_VAMCFM_NEW;
  float p0minus_VAMCFM_NEW;
  float p0_g1prime2_VAMCFM_NEW;
  float pg1g1prime2_VAMCFM_NEW;
  float pg1g1prime2_pi2_VAMCFM_NEW;
  float pg1g2_VAMCFM_NEW;
  float pg1g2_pi2_VAMCFM_NEW;
  float pg1g4_VAMCFM_NEW;
  float pg1g4_pi2_VAMCFM_NEW;

  float p0hplus_VAMCFM_ratio_NEW;
  float p0minus_VAMCFM_ratio_NEW;
  float p0_g1prime2_VAMCFM_ratio_NEW;
  float pg1g1prime2_VAMCFM_ratio_NEW;
  float pg1g1prime2_pi2_VAMCFM_ratio_NEW;
  float pg1g2_VAMCFM_ratio_NEW;
  float pg1g2_pi2_VAMCFM_ratio_NEW;
  float pg1g4_VAMCFM_ratio_NEW;
  float pg1g4_pi2_VAMCFM_ratio_NEW;

  float p0hplus_VAJHU_ratio_NEW;
  float p0minus_VAJHU_ratio_NEW;
  float p0_g1prime2_VAJHU_ratio_NEW;
  float pg1g1prime2_VAJHU_ratio_NEW;
  float pg1g1prime2_pi2_VAJHU_ratio_NEW;
  float pg1g2_VAJHU_ratio_NEW;
  float pg1g2_pi2_VAJHU_ratio_NEW;
  float pg1g4_VAJHU_ratio_NEW;
  float pg1g4_pi2_VAJHU_ratio_NEW;

  float p0plus_VAMCFM, ggzz_p0plus_VAMCFM;
  float bkg_VAMCFM, ggzz_VAMCFM;
  float p0plus_VAMCFM_NEW, ggzz_p0plus_VAMCFM_NEW, p0plus_VAMCFM_NEW_BSMOn;
  float bkg_VAMCFM_NEW, ggzz_VAMCFM_NEW;
  //float bkg_VAMCFM_STU, bkg_VAMCFM_TU, bkg_VAMCFM_S;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);
  tree->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
  tree->SetBranchAddress("p0hplus_VAJHU", &p0hplus_VAJHU);
  tree->SetBranchAddress("p0minus_VAJHU", &p0minus_VAJHU);
  tree->SetBranchAddress("p0_g1prime2_VAJHU", &p0_g1prime2_VAJHU);
  tree->SetBranchAddress("pg1g2_VAJHU", &pg1g2_VAJHU);
  tree->SetBranchAddress("pg1g2_pi2_VAJHU", &pg1g2_pi2_VAJHU);
  tree->SetBranchAddress("pg1g4_VAJHU", &pg1g4_VAJHU);
  tree->SetBranchAddress("pg1g4_pi2_VAJHU", &pg1g4_pi2_VAJHU);
  tree->SetBranchAddress("pg1g1prime2_VAJHU", &pg1g1prime2_VAJHU);
  tree->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
  tree->SetBranchAddress("ggzz_VAMCFM", &ggzz_VAMCFM);
  tree->SetBranchAddress("p0plus_VAMCFM", &p0plus_VAMCFM);
  tree->SetBranchAddress("ggzz_p0plus_VAMCFM", &ggzz_p0plus_VAMCFM);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("ZZMass", &mzz);
  newtree->Branch("Z1Mass", &m1);
  newtree->Branch("Z2Mass", &m2);
  newtree->Branch("helcosthetaZ1", &h1);
  newtree->Branch("helcosthetaZ2", &h2);
  newtree->Branch("helphi", &phi);
  newtree->Branch("costhetastar", &hs);
  newtree->Branch("phistarZ1", &phi1);

  newtree->Branch("p0plus_VAJHU", &p0plus_VAJHU);
  newtree->Branch("p0hplus_VAJHU", &p0hplus_VAJHU);
  newtree->Branch("p0minus_VAJHU", &p0minus_VAJHU);
  newtree->Branch("p0_g1prime2_VAJHU", &p0_g1prime2_VAJHU);
  newtree->Branch("pg1g2_VAJHU", &pg1g2_VAJHU);
  newtree->Branch("pg1g2_pi2_VAJHU", &pg1g2_pi2_VAJHU);
  newtree->Branch("pg1g4_VAJHU", &pg1g4_VAJHU);
  newtree->Branch("pg1g4_pi2_VAJHU", &pg1g4_pi2_VAJHU);
  newtree->Branch("pg1g1prime2_VAJHU", &pg1g1prime2_VAJHU);

  newtree->Branch("p0plus_VAJHU_NEW", &p0plus_VAJHU_NEW);
  newtree->Branch("p0hplus_VAJHU_NEW", &p0hplus_VAJHU_NEW);
  newtree->Branch("p0minus_VAJHU_NEW", &p0minus_VAJHU_NEW);
  newtree->Branch("p0_g1prime2_VAJHU_NEW", &p0_g1prime2_VAJHU_NEW);
  newtree->Branch("pg1g2_VAJHU_NEW", &pg1g2_VAJHU_NEW);
  newtree->Branch("pg1g2_pi2_VAJHU_NEW", &pg1g2_pi2_VAJHU_NEW);
  newtree->Branch("pg1g4_VAJHU_NEW", &pg1g4_VAJHU_NEW);
  newtree->Branch("pg1g4_pi2_VAJHU_NEW", &pg1g4_pi2_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_VAJHU_NEW", &pg1g1prime2_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAJHU_NEW", &pg1g1prime2_pi2_VAJHU_NEW);


  newtree->Branch("bkg_VAMCFM", &bkg_VAMCFM);
  newtree->Branch("ggzz_VAMCFM", &ggzz_VAMCFM);
  newtree->Branch("p0plus_VAMCFM", &p0plus_VAMCFM);
  newtree->Branch("ggzz_p0plus_VAMCFM", &ggzz_p0plus_VAMCFM);
  newtree->Branch("bkg_VAMCFM_NEW", &bkg_VAMCFM_NEW);
  newtree->Branch("ggzz_VAMCFM_NEW", &ggzz_VAMCFM_NEW);
  newtree->Branch("p0plus_VAMCFM_NEW", &p0plus_VAMCFM_NEW);
  newtree->Branch("ggzz_p0plus_VAMCFM_NEW", &ggzz_p0plus_VAMCFM_NEW);
  newtree->Branch("p0plus_VAMCFM_NEW_BSMOn", &p0plus_VAMCFM_NEW_BSMOn);

  newtree->Branch("p0hplus_VAMCFM_NEW", &p0hplus_VAMCFM_NEW);
  newtree->Branch("p0minus_VAMCFM_NEW", &p0minus_VAMCFM_NEW);
  newtree->Branch("p0_g1prime2_VAMCFM_NEW", &p0_g1prime2_VAMCFM_NEW);
  newtree->Branch("pg1g2_VAMCFM_NEW", &pg1g2_VAMCFM_NEW);
  newtree->Branch("pg1g2_pi2_VAMCFM_NEW", &pg1g2_pi2_VAMCFM_NEW);
  newtree->Branch("pg1g4_VAMCFM_NEW", &pg1g4_VAMCFM_NEW);
  newtree->Branch("pg1g4_pi2_VAMCFM_NEW", &pg1g4_pi2_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_VAMCFM_NEW", &pg1g1prime2_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAMCFM_NEW", &pg1g1prime2_pi2_VAMCFM_NEW);

  newtree->Branch("p0hplus_VAMCFM_ratio_NEW", &p0hplus_VAMCFM_ratio_NEW);
  newtree->Branch("p0minus_VAMCFM_ratio_NEW", &p0minus_VAMCFM_ratio_NEW);
  newtree->Branch("p0_g1prime2_VAMCFM_ratio_NEW", &p0_g1prime2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_VAMCFM_ratio_NEW", &pg1g2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_pi2_VAMCFM_ratio_NEW", &pg1g2_pi2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_VAMCFM_ratio_NEW", &pg1g4_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_pi2_VAMCFM_ratio_NEW", &pg1g4_pi2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_VAMCFM_ratio_NEW", &pg1g1prime2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAMCFM_ratio_NEW", &pg1g1prime2_pi2_VAMCFM_ratio_NEW);

  newtree->Branch("p0hplus_VAJHU_ratio_NEW", &p0hplus_VAJHU_ratio_NEW);
  newtree->Branch("p0minus_VAJHU_ratio_NEW", &p0minus_VAJHU_ratio_NEW);
  newtree->Branch("p0_g1prime2_VAJHU_ratio_NEW", &p0_g1prime2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_VAJHU_ratio_NEW", &pg1g2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_pi2_VAJHU_ratio_NEW", &pg1g2_pi2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_VAJHU_ratio_NEW", &pg1g4_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_pi2_VAJHU_ratio_NEW", &pg1g4_pi2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_VAJHU_ratio_NEW", &pg1g1prime2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAJHU_ratio_NEW", &pg1g1prime2_pi2_VAJHU_ratio_NEW);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  if (flavor == 2){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 1){
    GenLep1Id=11;
    GenLep2Id=-11;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 0){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 3){
    GenLep1Id=14;
    GenLep2Id=-14;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 4){
    GenLep1Id=0;
    GenLep2Id=-0;
    GenLep3Id=1;
    GenLep4Id=-1;
  }
  else if (flavor == 5){
    GenLep1Id=1;
    GenLep2Id=-1;
    GenLep3Id=2;
    GenLep4Id=-2;
  }
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (int ev = 0; ev < (!useBkgSample ? min(1000, (int)tree->GetEntries()) : (int)tree->GetEntries()); ev++){
    tree->GetEntry(ev);

    if (ev%10000 == 0) cout << "Processing event " << ev << endl;

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    TLorentzVector pOrdered[4];
    std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
    for (int ip=0; ip<min(4, (int)daus.size()); ip++) pOrdered[ip]=daus.at(ip);
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

    /***** ZZ *****/
    if (verbosity>=TVar::DEBUG) cout << "Computing bkg_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(bkg_VAMCFM_NEW, useConstants);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    //mela.setMelaLeptonInterference(TVar::InterfOff);
    mela.computeP(ggzz_VAMCFM_NEW, useConstants);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    //mela.setMelaLeptonInterference(TVar::InterfOff);
    mela.computeP(p0plus_VAMCFM_NEW, useConstants);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    //mela.setMelaLeptonInterference(TVar::InterfOff);
    mela.computeP(ggzz_p0plus_VAMCFM_NEW, useConstants);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p0plus_VAMCFM_NEW_BSMOn, useConstants);

    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p0plus_VAJHU_NEW, useConstants);

    if (!useBkgSample){
      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
      mela.computeP(p0hplus_VAJHU_NEW, useConstants);

      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
      mela.computeP(p0minus_VAJHU_NEW, useConstants);

      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-12046.01;
      mela.computeP(p0_g1prime2_VAJHU_NEW, useConstants);

      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
      mela.computeP(pg1g2_VAJHU_NEW, useConstants);
      pg1g2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.638;
      mela.computeP(pg1g2_pi2_VAJHU_NEW, useConstants);
      pg1g2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
      mela.computeP(pg1g4_VAJHU_NEW, useConstants);
      pg1g4_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=2.521;
      mela.computeP(pg1g4_pi2_VAJHU_NEW, useConstants);
      pg1g4_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=12046.01;
      mela.computeP(pg1g1prime2_VAJHU_NEW, useConstants);
      pg1g1prime2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);

      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][1]=12046.01;
      mela.computeP(pg1g1prime2_pi2_VAJHU_NEW, useConstants);
      pg1g1prime2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);


      mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.computeP(p0plus_VAMCFM_NEW, useConstants);

      mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
      mela.computeP(p0hplus_VAMCFM_NEW, useConstants);

      mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
      mela.computeP(p0minus_VAMCFM_NEW, useConstants);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-12046.01;
      mela.computeP(p0_g1prime2_VAMCFM_NEW, useConstants);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
      mela.computeP(pg1g2_VAMCFM_NEW, useConstants);
      pg1g2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0hplus_VAMCFM_NEW);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.638;
      mela.computeP(pg1g2_pi2_VAMCFM_NEW, useConstants);
      pg1g2_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0hplus_VAMCFM_NEW);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
      mela.computeP(pg1g4_VAMCFM_NEW, useConstants);
      pg1g4_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0minus_VAMCFM_NEW);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=2.521;
      mela.computeP(pg1g4_pi2_VAMCFM_NEW, useConstants);
      pg1g4_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0minus_VAMCFM_NEW);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=12046.01;
      mela.computeP(pg1g1prime2_VAMCFM_NEW, useConstants);
      pg1g1prime2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0_g1prime2_VAMCFM_NEW);

      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][1]=12046.01;
      mela.computeP(pg1g1prime2_pi2_VAMCFM_NEW, useConstants);
      pg1g1prime2_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0_g1prime2_VAMCFM_NEW);
    }

    p0hplus_VAJHU_ratio_NEW = p0hplus_VAJHU_NEW/p0plus_VAJHU_NEW;
    p0minus_VAJHU_ratio_NEW = p0minus_VAJHU_NEW/p0plus_VAJHU_NEW;
    p0_g1prime2_VAJHU_ratio_NEW = p0_g1prime2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g2_VAJHU_ratio_NEW = pg1g2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g4_VAJHU_ratio_NEW = pg1g4_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g1prime2_VAJHU_ratio_NEW = pg1g1prime2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g2_pi2_VAJHU_ratio_NEW = pg1g2_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g4_pi2_VAJHU_ratio_NEW = pg1g4_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g1prime2_pi2_VAJHU_ratio_NEW = pg1g1prime2_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;

    p0hplus_VAMCFM_ratio_NEW = p0hplus_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    p0minus_VAMCFM_ratio_NEW = p0minus_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    p0_g1prime2_VAMCFM_ratio_NEW = p0_g1prime2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g2_VAMCFM_ratio_NEW = pg1g2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g4_VAMCFM_ratio_NEW = pg1g4_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g1prime2_VAMCFM_ratio_NEW = pg1g1prime2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g2_pi2_VAMCFM_ratio_NEW = pg1g2_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g4_pi2_VAMCFM_ratio_NEW = pg1g4_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g1prime2_pi2_VAMCFM_ratio_NEW = pg1g1prime2_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;


    mela.resetInputEvent();
    newtree->Fill();
  }


  foutput->WriteTObject(newtree);
  foutput->Close();
  finput->Close();
  mela.setVerbosity(bkpverbosity);
}

void testME_ProdP_VBFHJJ_FullSim(int flavor=2, bool useConstants=false, bool useBkgSample=false, shared_ptr<Mela> melaptr=nullptr){
  int erg_tev=13;
  float mPOLE=125.6;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  if (!melaptr) {
    melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  }
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput;
  TFile* foutput;
  if (useBkgSample){
    finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor>=2 ? "2mu2e" : "4e"), (flavor==2 ? "2e2mu" : "4e")), "read");
    foutput = new TFile(Form("HZZ4lTree_ZZTo%s_vbfMELATest.root", (flavor>=2 ? "2e2mu" : "4e")), "recreate");
  }
  else{
    finput = new TFile(Form("%s/%s/HZZ4lTree_VBF0P_H125.6.root", cinput_main.Data(), (flavor==2 ? "2mu2e" : "4e")), "read");
    foutput = new TFile("HZZ4lTree_VBF0P_H125.6_vbfMELATest.root", "recreate");
  }

  float phjj_VAJHU_old;
  float pvbf_VAJHU_old;
  float phjj_VAJHU_old_NEW;
  float pvbf_VAJHU_old_NEW;
  float phjj_VAJHU_old_NEW_selfD;
  float pvbf_VAJHU_old_NEW_selfD;
  float phjj0minus_VAJHU_old_NEW;
  float pvbf0minus_VAJHU_old_NEW;
  float phjj0minus_VAJHU_old_NEW_selfD;
  float pvbf0minus_VAJHU_old_NEW_selfD;

  //float jet1Pt, jet2Pt;
  //float jet1px, jet1py, jet1pz, jet1E;
  //float jet2px, jet2py, jet2pz, jet2E;
  //float ZZPx, ZZPy, ZZPz, ZZE, dR;
  short NJets30;
  std::vector<double>* JetPt=0;
  std::vector<double>* JetEta=0;
  std::vector<double>* JetPhi=0;
  std::vector<double>* JetMass=0;
  std::vector<double> myJetPt;
  std::vector<double> myJetEta;
  std::vector<double> myJetPhi;
  std::vector<double> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tree->SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
  tree->SetBranchAddress("JetMass", &JetMass, &bJetMass);
  tree->SetBranchAddress("phjj_VAJHU_old", &phjj_VAJHU_old);
  tree->SetBranchAddress("pvbf_VAJHU_old", &pvbf_VAJHU_old);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("phjj_VAJHU_old", &phjj_VAJHU_old);
  newtree->Branch("pvbf_VAJHU_old", &pvbf_VAJHU_old);
  newtree->Branch("phjj_VAJHU_old_NEW", &phjj_VAJHU_old_NEW);
  newtree->Branch("pvbf_VAJHU_old_NEW", &pvbf_VAJHU_old_NEW);
  newtree->Branch("phjj0minus_VAJHU_old_NEW", &phjj0minus_VAJHU_old_NEW);
  newtree->Branch("pvbf0minus_VAJHU_old_NEW", &pvbf0minus_VAJHU_old_NEW);
  newtree->Branch("phjj_VAJHU_old_NEW_selfD", &phjj_VAJHU_old_NEW_selfD);
  newtree->Branch("pvbf_VAJHU_old_NEW_selfD", &pvbf_VAJHU_old_NEW_selfD);
  newtree->Branch("phjj0minus_VAJHU_old_NEW_selfD", &phjj0minus_VAJHU_old_NEW_selfD);
  newtree->Branch("pvbf0minus_VAJHU_old_NEW_selfD", &pvbf0minus_VAJHU_old_NEW_selfD);
  newtree->Branch("ZZMass", &mzz);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  if (flavor == 2){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 1){
    GenLep1Id=11;
    GenLep2Id=-11;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 0){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 3){
    GenLep1Id=14;
    GenLep2Id=-14;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 4){
    GenLep1Id=0;
    GenLep2Id=-0;
    GenLep3Id=1;
    GenLep4Id=-1;
  }
  else if (flavor == 5){
    GenLep1Id=1;
    GenLep2Id=-1;
    GenLep3Id=2;
    GenLep4Id=-2;
  }
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };

  int nEntries = tree->GetEntries();
  int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=(verbosity>=TVar::DEBUG ? 1 : (!useBkgSample ? 1000 : nEntries))) break;
    tree->GetEntry(ev);

    if (JetPt->size()>=2 && NJets30>=2){
      TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet1));
      associated.push_back(SimpleParticle_t(0, jet2));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters_ZZ;
      for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters_ZZ, &associated, (SimpleParticleCollection_t*)0, false);


      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
      mela.computeProdP(pvbf_VAJHU_old_NEW, useConstants);

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
      mela.computeProdP(phjj_VAJHU_old_NEW, useConstants);
      /*
      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJVBF);
      mela.computeProdP(pvbf0minus_VAJHU_old_NEW, useConstants);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJQCD);
      mela.computeProdP(phjj0minus_VAJHU_old_NEW, useConstants);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.computeProdP(pvbf_VAJHU_old_NEW_selfD, useConstants);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJQCD);
      mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
      mela.computeProdP(phjj_VAJHU_old_NEW_selfD, useConstants);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
      mela.computeProdP(pvbf0minus_VAJHU_old_NEW_selfD, useConstants);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJQCD);
      mela.selfDHggcoupl[2][0]=1;
      mela.computeProdP(phjj0minus_VAJHU_old_NEW_selfD, useConstants);
      */
      newtree->Fill();
      recorded++;
      mela.resetInputEvent();
    }
  }


  foutput->WriteTObject(newtree);
  delete newtree;
  foutput->Close();
  finput->Close();
  mela.setVerbosity(bkpverbosity);
}

void testME_ProdP_VH_FullSim(shared_ptr<Mela> melaptr=nullptr){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::DEBUG;
  if (!melaptr) {
    melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  }
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), "2mu2e", "2e2mu"), "read");
  TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_vhMELATest.root", "2e2mu"), "recreate");

  float pwh_leptonic_VAJHU_old_NEW;
  float pzh_leptonic_VAJHU_old_NEW;
  float pwh_leptonic_VAJHU_old_NEW_selfD;
  float pzh_leptonic_VAJHU_old_NEW_selfD;
  float pwh0minus_leptonic_VAJHU_old_NEW;
  float pzh0minus_leptonic_VAJHU_old_NEW;
  float pwh0minus_leptonic_VAJHU_old_NEW_selfD;
  float pzh0minus_leptonic_VAJHU_old_NEW_selfD;
  float pwh_hadronic_VAJHU_old_NEW;
  float pzh_hadronic_VAJHU_old_NEW;
  float pwh_hadronic_VAJHU_old_NEW_selfD;
  float pzh_hadronic_VAJHU_old_NEW_selfD;
  float pwh0minus_hadronic_VAJHU_old_NEW;
  float pzh0minus_hadronic_VAJHU_old_NEW;
  float pwh0minus_hadronic_VAJHU_old_NEW_selfD;
  float pzh0minus_hadronic_VAJHU_old_NEW_selfD;

  //float jet1Pt, jet2Pt;
  //float jet1px, jet1py, jet1pz, jet1E;
  //float jet2px, jet2py, jet2pz, jet2E;
  //float ZZPx, ZZPy, ZZPz, ZZE, dR;
  short NJets30;
  std::vector<double>* JetPt=0;
  std::vector<double>* JetEta=0;
  std::vector<double>* JetPhi=0;
  std::vector<double>* JetMass=0;
  std::vector<double> myJetPt;
  std::vector<double> myJetEta;
  std::vector<double> myJetPhi;
  std::vector<double> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tree->SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
  tree->SetBranchAddress("JetMass", &JetMass, &bJetMass);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("pwh_leptonic_VAJHU_old_NEW", &pwh_leptonic_VAJHU_old_NEW);
  newtree->Branch("pzh_leptonic_VAJHU_old_NEW", &pzh_leptonic_VAJHU_old_NEW);
  newtree->Branch("pwh0minus_leptonic_VAJHU_old_NEW", &pwh0minus_leptonic_VAJHU_old_NEW);
  newtree->Branch("pzh0minus_leptonic_VAJHU_old_NEW", &pzh0minus_leptonic_VAJHU_old_NEW);
  newtree->Branch("pwh_leptonic_VAJHU_old_NEW_selfD", &pwh_leptonic_VAJHU_old_NEW_selfD);
  newtree->Branch("pzh_leptonic_VAJHU_old_NEW_selfD", &pzh_leptonic_VAJHU_old_NEW_selfD);
  newtree->Branch("pwh0minus_leptonic_VAJHU_old_NEW_selfD", &pwh0minus_leptonic_VAJHU_old_NEW_selfD);
  newtree->Branch("pzh0minus_leptonic_VAJHU_old_NEW_selfD", &pzh0minus_leptonic_VAJHU_old_NEW_selfD);
  newtree->Branch("pwh_hadronic_VAJHU_old_NEW", &pwh_hadronic_VAJHU_old_NEW);
  newtree->Branch("pzh_hadronic_VAJHU_old_NEW", &pzh_hadronic_VAJHU_old_NEW);
  newtree->Branch("pwh0minus_hadronic_VAJHU_old_NEW", &pwh0minus_hadronic_VAJHU_old_NEW);
  newtree->Branch("pzh0minus_hadronic_VAJHU_old_NEW", &pzh0minus_hadronic_VAJHU_old_NEW);
  newtree->Branch("pwh_hadronic_VAJHU_old_NEW_selfD", &pwh_hadronic_VAJHU_old_NEW_selfD);
  newtree->Branch("pzh_hadronic_VAJHU_old_NEW_selfD", &pzh_hadronic_VAJHU_old_NEW_selfD);
  newtree->Branch("pwh0minus_hadronic_VAJHU_old_NEW_selfD", &pwh0minus_hadronic_VAJHU_old_NEW_selfD);
  newtree->Branch("pzh0minus_hadronic_VAJHU_old_NEW_selfD", &pzh0minus_hadronic_VAJHU_old_NEW_selfD);
  newtree->Branch("ZZMass", &mzz);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  GenLep1Id=13;
  GenLep2Id=-13;
  GenLep3Id=14;
  GenLep4Id=-14;
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };

  int nEntries = tree->GetEntries();
  int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=(verbosity>=TVar::DEBUG ? 1 : 1000)) break;
    tree->GetEntry(ev);

    if (JetPt->size()>=2 && NJets30>=2){
      TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet1));
      associated.push_back(SimpleParticle_t(0, jet2));
      for (unsigned int idau=0; idau<4; idau++) associated.push_back(SimpleParticle_t(idOrdered[idau], pDaughters[idau])); // Let's put some extra leptons

      SimpleParticleCollection_t daughters_ZZ;
      for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters_ZZ, &associated, (SimpleParticleCollection_t*)0, false);

      if (verbosity>=TVar::DEBUG){
        cout << "Mela candidates summary:" << endl;
        for (int ic=0; ic<mela.getNCandidates(); ic++){
          cout << "Candidate " << ic << endl;
          mela.setCurrentCandidateFromIndex(ic);
          TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
        }
        cout << endl;
      }

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Lep_ZH);
      mela.computeProdP_VH(pzh_leptonic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Lep_WH);
      mela.computeProdP_VH(pwh_leptonic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::Lep_ZH);
      mela.computeProdP_VH(pzh0minus_leptonic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::Lep_WH);
      mela.computeProdP_VH(pwh0minus_leptonic_VAJHU_old_NEW, false, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_ZH);
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.computeProdP_VH(pzh_leptonic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_WH);
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.computeProdP_VH(pwh_leptonic_VAJHU_old_NEW_selfD, false, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_ZH);
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
      mela.computeProdP_VH(pzh0minus_leptonic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Lep_WH);
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
      mela.computeProdP_VH(pwh0minus_leptonic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_ZH);
      mela.computeProdP_VH(pzh_hadronic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::Had_WH);
      mela.computeProdP_VH(pwh_hadronic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::Had_ZH);
      mela.computeProdP_VH(pzh0minus_hadronic_VAJHU_old_NEW, false, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::Had_WH);
      mela.computeProdP_VH(pwh0minus_hadronic_VAJHU_old_NEW, false, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_ZH);
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.computeProdP_VH(pzh_hadronic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_WH);
      mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
      mela.computeProdP_VH(pwh_hadronic_VAJHU_old_NEW_selfD, false, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_ZH);
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
      mela.computeProdP_VH(pzh0minus_hadronic_VAJHU_old_NEW_selfD, false, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::Had_WH);
      mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=1;
      mela.computeProdP_VH(pwh0minus_hadronic_VAJHU_old_NEW_selfD, false, false);

      newtree->Fill();
      recorded++;
      mela.resetInputEvent();
    }
  }


  foutput->WriteTObject(newtree);
  delete newtree;
  foutput->Close();
  finput->Close();
  mela.setVerbosity(bkpverbosity);
}

void testME_ProdP_TTHBBH_FullSim(shared_ptr<Mela> melaptr=nullptr){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  if (!melaptr) {
    melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  }
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), "2mu2e", "2e2mu"), "read");
  TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_ttHbbHMELATest.root", "2e2mu"), "recreate");

  float pbbh_VAJHU_old_NEW;
  float ptth_VAJHU_old_NEW;
  float pbbh_VAJHU_old_NEW_selfD;
  float ptth_VAJHU_old_NEW_selfD;
  float pbbh0minus_VAJHU_old_NEW;
  float ptth0minus_VAJHU_old_NEW;
  float pbbh0minus_VAJHU_old_NEW_selfD;
  float ptth0minus_VAJHU_old_NEW_selfD;

  short NJets30;
  std::vector<double>* JetPt=0;
  std::vector<double>* JetEta=0;
  std::vector<double>* JetPhi=0;
  std::vector<double>* JetMass=0;
  std::vector<double> myJetPt;
  std::vector<double> myJetEta;
  std::vector<double> myJetPhi;
  std::vector<double> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tree->SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
  tree->SetBranchAddress("JetMass", &JetMass, &bJetMass);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("pbbh_VAJHU_old_NEW", &pbbh_VAJHU_old_NEW);
  newtree->Branch("ptth_VAJHU_old_NEW", &ptth_VAJHU_old_NEW);
  newtree->Branch("pbbh0minus_VAJHU_old_NEW", &pbbh0minus_VAJHU_old_NEW);
  newtree->Branch("ptth0minus_VAJHU_old_NEW", &ptth0minus_VAJHU_old_NEW);
  newtree->Branch("pbbh_VAJHU_old_NEW_selfD", &pbbh_VAJHU_old_NEW_selfD);
  newtree->Branch("ptth_VAJHU_old_NEW_selfD", &ptth_VAJHU_old_NEW_selfD);
  newtree->Branch("pbbh0minus_VAJHU_old_NEW_selfD", &pbbh0minus_VAJHU_old_NEW_selfD);
  newtree->Branch("ptth0minus_VAJHU_old_NEW_selfD", &ptth0minus_VAJHU_old_NEW_selfD);
  newtree->Branch("ZZMass", &mzz);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  GenLep1Id=13;
  GenLep2Id=-13;
  GenLep3Id=11;
  GenLep4Id=-11;
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
  int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };

  int nEntries = tree->GetEntries();
  int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=1000) break;
    tree->GetEntry(ev);

    if (JetPt->size()>=2 && NJets30>=2){
      TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet1));
      associated.push_back(SimpleParticle_t(0, jet2));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(idOrdered[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      // This is some hack to get top daughter four-vectors
      SimpleParticleCollection_t topDaughters;
      topDaughters.push_back(SimpleParticle_t(0, jet1));
      for (unsigned int idau=0; idau<2; idau++) topDaughters.push_back(SimpleParticle_t(0, pDaughters[idau]));
      mela.appendTopCandidate(&topDaughters);
      SimpleParticleCollection_t antitopDaughters;
      antitopDaughters.push_back(SimpleParticle_t(0, jet2));
      for (unsigned int idau=2; idau<4; idau++) antitopDaughters.push_back(SimpleParticle_t(0, pDaughters[idau]));
      mela.appendTopCandidate(&antitopDaughters);

      if (verbosity>=TVar::DEBUG){
        cout << "Mela candidates summary:" << endl;
        for (int ic=0; ic<mela.getNCandidates(); ic++){
          cout << "Candidate " << ic << endl;
          mela.setCurrentCandidateFromIndex(ic);
          TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
        }
        cout << endl;
      }

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ttH);
      mela.computeProdP_ttH(ptth_VAJHU_old_NEW, 2, 0, false);

      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::bbH);
      mela.computeProdP_ttH(pbbh_VAJHU_old_NEW, 2, 0, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::ttH);
      mela.computeProdP_ttH(ptth0minus_VAJHU_old_NEW, 2, 0, false);

      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::bbH);
      mela.computeProdP_ttH(pbbh0minus_VAJHU_old_NEW, 2, 0, false);


      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ttH);
      mela.selfDHqqcoupl[0][gHIGGS_KAPPA][0]=1;
      mela.computeProdP_ttH(ptth_VAJHU_old_NEW_selfD, 2, 0, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::bbH);
      mela.selfDHqqcoupl[0][gHIGGS_KAPPA][0]=1;
      mela.computeProdP_ttH(pbbh_VAJHU_old_NEW_selfD, 2, 0, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ttH);
      mela.selfDHqqcoupl[0][gHIGGS_KAPPA_TILDE][0]=1;
      mela.computeProdP_ttH(ptth0minus_VAJHU_old_NEW_selfD, 2, 0, false);

      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::bbH);
      mela.selfDHqqcoupl[0][gHIGGS_KAPPA_TILDE][0]=1;
      mela.computeProdP_ttH(pbbh0minus_VAJHU_old_NEW_selfD, 2, 0, false);

      newtree->Fill();
      recorded++;
      mela.resetInputEvent();
    }
  }


  foutput->WriteTObject(newtree);
  delete newtree;
  foutput->Close();
  finput->Close();
  mela.setVerbosity(bkpverbosity);
}


void testME_Dec_ZZWWComparison_FullSim(shared_ptr<Mela> melaptr=nullptr){
  int erg_tev=8;
  float mPOLE=125.;
  float wPOLE=4.07e-3;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  if (!melaptr) {
    melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  }
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;
  mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root", cinput_main.Data(), "2mu2e"), "read");
  TFile* foutput = new TFile(Form("HZZ4lTree_powheg15jhuGenV3-0PMH125.6_%s_ZZWWComparison.root", "2l2nu"), "recreate");

  float p0plus_VAJHU_NEW;
  float p0hplus_VAJHU_NEW;
  float p0minus_VAJHU_NEW;
  float p0_g1prime2_VAJHU_NEW;
  float pg1g1prime2_VAJHU_NEW;
  float pg1g1prime2_pi2_VAJHU_NEW;
  float pg1g2_VAJHU_NEW;
  float pg1g2_pi2_VAJHU_NEW;
  float pg1g4_VAJHU_NEW;
  float pg1g4_pi2_VAJHU_NEW;

  float p0hplus_VAMCFM_NEW;
  float p0minus_VAMCFM_NEW;
  float p0_g1prime2_VAMCFM_NEW;
  float pg1g1prime2_VAMCFM_NEW;
  float pg1g1prime2_pi2_VAMCFM_NEW;
  float pg1g2_VAMCFM_NEW;
  float pg1g2_pi2_VAMCFM_NEW;
  float pg1g4_VAMCFM_NEW;
  float pg1g4_pi2_VAMCFM_NEW;

  float p0hplus_VAMCFM_ratio_NEW;
  float p0minus_VAMCFM_ratio_NEW;
  float p0_g1prime2_VAMCFM_ratio_NEW;
  float pg1g1prime2_VAMCFM_ratio_NEW;
  float pg1g1prime2_pi2_VAMCFM_ratio_NEW;
  float pg1g2_VAMCFM_ratio_NEW;
  float pg1g2_pi2_VAMCFM_ratio_NEW;
  float pg1g4_VAMCFM_ratio_NEW;
  float pg1g4_pi2_VAMCFM_ratio_NEW;

  float p0hplus_VAJHU_ratio_NEW;
  float p0minus_VAJHU_ratio_NEW;
  float p0_g1prime2_VAJHU_ratio_NEW;
  float pg1g1prime2_VAJHU_ratio_NEW;
  float pg1g1prime2_pi2_VAJHU_ratio_NEW;
  float pg1g2_VAJHU_ratio_NEW;
  float pg1g2_pi2_VAJHU_ratio_NEW;
  float pg1g4_VAJHU_ratio_NEW;
  float pg1g4_pi2_VAJHU_ratio_NEW;

  float p0plus_VAMCFM_NEW, ggzz_p0plus_VAMCFM_NEW, p0plus_VAMCFM_NEW_BSMOn, ggzz_p0plus_VAMCFM_NEW_BSMOn;
  float bkg_VAMCFM_NEW, ggzz_VAMCFM_NEW;

  float p0plus_WW_VAJHU_NEW;
  float p0hplus_WW_VAJHU_NEW;
  float p0minus_WW_VAJHU_NEW;
  float p0_g1prime2_WW_VAJHU_NEW;
  float pg1g1prime2_WW_VAJHU_NEW;
  float pg1g1prime2_pi2_WW_VAJHU_NEW;
  float pg1g2_WW_VAJHU_NEW;
  float pg1g2_pi2_WW_VAJHU_NEW;
  float pg1g4_WW_VAJHU_NEW;
  float pg1g4_pi2_WW_VAJHU_NEW;

  float p0hplus_WW_VAMCFM_NEW;
  float p0minus_WW_VAMCFM_NEW;
  float p0_g1prime2_WW_VAMCFM_NEW;
  float pg1g1prime2_WW_VAMCFM_NEW;
  float pg1g1prime2_pi2_WW_VAMCFM_NEW;
  float pg1g2_WW_VAMCFM_NEW;
  float pg1g2_pi2_WW_VAMCFM_NEW;
  float pg1g4_WW_VAMCFM_NEW;
  float pg1g4_pi2_WW_VAMCFM_NEW;

  float p0hplus_WW_VAMCFM_ratio_NEW;
  float p0minus_WW_VAMCFM_ratio_NEW;
  float p0_g1prime2_WW_VAMCFM_ratio_NEW;
  float pg1g1prime2_WW_VAMCFM_ratio_NEW;
  float pg1g1prime2_pi2_WW_VAMCFM_ratio_NEW;
  float pg1g2_WW_VAMCFM_ratio_NEW;
  float pg1g2_pi2_WW_VAMCFM_ratio_NEW;
  float pg1g4_WW_VAMCFM_ratio_NEW;
  float pg1g4_pi2_WW_VAMCFM_ratio_NEW;

  float p0hplus_WW_VAJHU_ratio_NEW;
  float p0minus_WW_VAJHU_ratio_NEW;
  float p0_g1prime2_WW_VAJHU_ratio_NEW;
  float pg1g1prime2_WW_VAJHU_ratio_NEW;
  float pg1g1prime2_pi2_WW_VAJHU_ratio_NEW;
  float pg1g2_WW_VAJHU_ratio_NEW;
  float pg1g2_pi2_WW_VAJHU_ratio_NEW;
  float pg1g4_WW_VAJHU_ratio_NEW;
  float pg1g4_pi2_WW_VAJHU_ratio_NEW;

  float p0plus_WW_VAMCFM_NEW, ggzz_p0plus_WW_VAMCFM_NEW, p0plus_WW_VAMCFM_NEW_BSMOn, ggzz_p0plus_WW_VAMCFM_NEW_BSMOn;
  float bkg_WW_VAMCFM_NEW, ggzz_WW_VAMCFM_NEW;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("ZZMass", &mzz);
  newtree->Branch("Z1Mass", &m1);
  newtree->Branch("Z2Mass", &m2);
  newtree->Branch("helcosthetaZ1", &h1);
  newtree->Branch("helcosthetaZ2", &h2);
  newtree->Branch("helphi", &phi);
  newtree->Branch("costhetastar", &hs);
  newtree->Branch("phistarZ1", &phi1);

  newtree->Branch("p0plus_VAJHU_NEW", &p0plus_VAJHU_NEW);
  newtree->Branch("p0hplus_VAJHU_NEW", &p0hplus_VAJHU_NEW);
  newtree->Branch("p0minus_VAJHU_NEW", &p0minus_VAJHU_NEW);
  newtree->Branch("p0_g1prime2_VAJHU_NEW", &p0_g1prime2_VAJHU_NEW);
  newtree->Branch("pg1g2_VAJHU_NEW", &pg1g2_VAJHU_NEW);
  newtree->Branch("pg1g2_pi2_VAJHU_NEW", &pg1g2_pi2_VAJHU_NEW);
  newtree->Branch("pg1g4_VAJHU_NEW", &pg1g4_VAJHU_NEW);
  newtree->Branch("pg1g4_pi2_VAJHU_NEW", &pg1g4_pi2_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_VAJHU_NEW", &pg1g1prime2_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAJHU_NEW", &pg1g1prime2_pi2_VAJHU_NEW);

  newtree->Branch("bkg_VAMCFM_NEW", &bkg_VAMCFM_NEW);
  newtree->Branch("ggzz_VAMCFM_NEW", &ggzz_VAMCFM_NEW);
  newtree->Branch("p0plus_VAMCFM_NEW", &p0plus_VAMCFM_NEW);
  newtree->Branch("ggzz_p0plus_VAMCFM_NEW", &ggzz_p0plus_VAMCFM_NEW);
  newtree->Branch("p0plus_VAMCFM_NEW_BSMOn", &p0plus_VAMCFM_NEW_BSMOn);
  newtree->Branch("ggzz_p0plus_VAMCFM_NEW_BSMOn", &ggzz_p0plus_VAMCFM_NEW_BSMOn);

  newtree->Branch("p0hplus_VAMCFM_NEW", &p0hplus_VAMCFM_NEW);
  newtree->Branch("p0minus_VAMCFM_NEW", &p0minus_VAMCFM_NEW);
  newtree->Branch("p0_g1prime2_VAMCFM_NEW", &p0_g1prime2_VAMCFM_NEW);
  newtree->Branch("pg1g2_VAMCFM_NEW", &pg1g2_VAMCFM_NEW);
  newtree->Branch("pg1g2_pi2_VAMCFM_NEW", &pg1g2_pi2_VAMCFM_NEW);
  newtree->Branch("pg1g4_VAMCFM_NEW", &pg1g4_VAMCFM_NEW);
  newtree->Branch("pg1g4_pi2_VAMCFM_NEW", &pg1g4_pi2_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_VAMCFM_NEW", &pg1g1prime2_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAMCFM_NEW", &pg1g1prime2_pi2_VAMCFM_NEW);

  newtree->Branch("p0hplus_VAMCFM_ratio_NEW", &p0hplus_VAMCFM_ratio_NEW);
  newtree->Branch("p0minus_VAMCFM_ratio_NEW", &p0minus_VAMCFM_ratio_NEW);
  newtree->Branch("p0_g1prime2_VAMCFM_ratio_NEW", &p0_g1prime2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_VAMCFM_ratio_NEW", &pg1g2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_pi2_VAMCFM_ratio_NEW", &pg1g2_pi2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_VAMCFM_ratio_NEW", &pg1g4_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_pi2_VAMCFM_ratio_NEW", &pg1g4_pi2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_VAMCFM_ratio_NEW", &pg1g1prime2_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAMCFM_ratio_NEW", &pg1g1prime2_pi2_VAMCFM_ratio_NEW);

  newtree->Branch("p0hplus_VAJHU_ratio_NEW", &p0hplus_VAJHU_ratio_NEW);
  newtree->Branch("p0minus_VAJHU_ratio_NEW", &p0minus_VAJHU_ratio_NEW);
  newtree->Branch("p0_g1prime2_VAJHU_ratio_NEW", &p0_g1prime2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_VAJHU_ratio_NEW", &pg1g2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_pi2_VAJHU_ratio_NEW", &pg1g2_pi2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_VAJHU_ratio_NEW", &pg1g4_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_pi2_VAJHU_ratio_NEW", &pg1g4_pi2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_VAJHU_ratio_NEW", &pg1g1prime2_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_VAJHU_ratio_NEW", &pg1g1prime2_pi2_VAJHU_ratio_NEW);

  newtree->Branch("p0plus_WW_VAJHU_NEW", &p0plus_WW_VAJHU_NEW);
  newtree->Branch("p0hplus_WW_VAJHU_NEW", &p0hplus_WW_VAJHU_NEW);
  newtree->Branch("p0minus_WW_VAJHU_NEW", &p0minus_WW_VAJHU_NEW);
  newtree->Branch("p0_g1prime2_WW_VAJHU_NEW", &p0_g1prime2_WW_VAJHU_NEW);
  newtree->Branch("pg1g2_WW_VAJHU_NEW", &pg1g2_WW_VAJHU_NEW);
  newtree->Branch("pg1g2_pi2_WW_VAJHU_NEW", &pg1g2_pi2_WW_VAJHU_NEW);
  newtree->Branch("pg1g4_WW_VAJHU_NEW", &pg1g4_WW_VAJHU_NEW);
  newtree->Branch("pg1g4_pi2_WW_VAJHU_NEW", &pg1g4_pi2_WW_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_WW_VAJHU_NEW", &pg1g1prime2_WW_VAJHU_NEW);
  newtree->Branch("pg1g1prime2_pi2_WW_VAJHU_NEW", &pg1g1prime2_pi2_WW_VAJHU_NEW);

  newtree->Branch("bkg_WW_VAMCFM_NEW", &bkg_WW_VAMCFM_NEW);
  newtree->Branch("ggzz_WW_VAMCFM_NEW", &ggzz_WW_VAMCFM_NEW);
  newtree->Branch("p0plus_WW_VAMCFM_NEW", &p0plus_WW_VAMCFM_NEW);
  newtree->Branch("ggzz_p0plus_WW_VAMCFM_NEW", &ggzz_p0plus_WW_VAMCFM_NEW);
  newtree->Branch("p0plus_WW_VAMCFM_NEW_BSMOn", &p0plus_WW_VAMCFM_NEW_BSMOn);
  newtree->Branch("ggzz_p0plus_WW_VAMCFM_NEW_BSMOn", &ggzz_p0plus_WW_VAMCFM_NEW_BSMOn);

  newtree->Branch("p0hplus_WW_VAMCFM_NEW", &p0hplus_WW_VAMCFM_NEW);
  newtree->Branch("p0minus_WW_VAMCFM_NEW", &p0minus_WW_VAMCFM_NEW);
  newtree->Branch("p0_g1prime2_WW_VAMCFM_NEW", &p0_g1prime2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g2_WW_VAMCFM_NEW", &pg1g2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g2_pi2_WW_VAMCFM_NEW", &pg1g2_pi2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g4_WW_VAMCFM_NEW", &pg1g4_WW_VAMCFM_NEW);
  newtree->Branch("pg1g4_pi2_WW_VAMCFM_NEW", &pg1g4_pi2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_WW_VAMCFM_NEW", &pg1g1prime2_WW_VAMCFM_NEW);
  newtree->Branch("pg1g1prime2_pi2_WW_VAMCFM_NEW", &pg1g1prime2_pi2_WW_VAMCFM_NEW);

  newtree->Branch("p0hplus_WW_VAMCFM_ratio_NEW", &p0hplus_WW_VAMCFM_ratio_NEW);
  newtree->Branch("p0minus_WW_VAMCFM_ratio_NEW", &p0minus_WW_VAMCFM_ratio_NEW);
  newtree->Branch("p0_g1prime2_WW_VAMCFM_ratio_NEW", &p0_g1prime2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_WW_VAMCFM_ratio_NEW", &pg1g2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g2_pi2_WW_VAMCFM_ratio_NEW", &pg1g2_pi2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_WW_VAMCFM_ratio_NEW", &pg1g4_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g4_pi2_WW_VAMCFM_ratio_NEW", &pg1g4_pi2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_WW_VAMCFM_ratio_NEW", &pg1g1prime2_WW_VAMCFM_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_WW_VAMCFM_ratio_NEW", &pg1g1prime2_pi2_WW_VAMCFM_ratio_NEW);

  newtree->Branch("p0hplus_WW_VAJHU_ratio_NEW", &p0hplus_WW_VAJHU_ratio_NEW);
  newtree->Branch("p0minus_WW_VAJHU_ratio_NEW", &p0minus_WW_VAJHU_ratio_NEW);
  newtree->Branch("p0_g1prime2_WW_VAJHU_ratio_NEW", &p0_g1prime2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_WW_VAJHU_ratio_NEW", &pg1g2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g2_pi2_WW_VAJHU_ratio_NEW", &pg1g2_pi2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_WW_VAJHU_ratio_NEW", &pg1g4_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g4_pi2_WW_VAJHU_ratio_NEW", &pg1g4_pi2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_WW_VAJHU_ratio_NEW", &pg1g1prime2_WW_VAJHU_ratio_NEW);
  newtree->Branch("pg1g1prime2_pi2_WW_VAJHU_ratio_NEW", &pg1g1prime2_pi2_WW_VAJHU_ratio_NEW);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  GenLep1Id=14;
  GenLep2Id=-14;
  GenLep3Id=13;
  GenLep4Id=-13;

  for (int ev = 0; ev < min(1000, (int)tree->GetEntries()); ev++){
    tree->GetEntry(ev);

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    TLorentzVector pOrdered[4];
    std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
    for (int ip=0; ip<min(4, (int)daus.size()); ip++) pOrdered[ip]=daus.at(ip);
    SimpleParticleCollection_t daughters;
    for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));
    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);
    mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    mela.setCandidateDecayMode(TVar::CandidateDecay_WW);
    mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
    if (verbosity>=TVar::DEBUG){
      cout << "Mela candidate summary:" << endl;
      for (int ic=0; ic<mela.getNCandidates(); ic++){
        cout << "Candidate " << ic << endl;
        mela.setCurrentCandidateFromIndex(ic);
        TUtil::PrintCandidateSummary(mela.getCurrentCandidate());
      }
      cout << endl;
    }

    /***** ZZ *****/
    mela.setCurrentCandidateFromIndex(0);

    if (verbosity>=TVar::DEBUG) cout << "Computing bkg_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(bkg_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    mela.computeP(ggzz_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(p0plus_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(ggzz_p0plus_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p0plus_VAMCFM_NEW_BSMOn, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(ggzz_p0plus_VAMCFM_NEW_BSMOn, false);

    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p0plus_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
    mela.computeP(p0hplus_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
    mela.computeP(p0minus_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-12046.01;
    mela.computeP(p0_g1prime2_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
    mela.computeP(pg1g2_VAJHU_NEW, false);
    pg1g2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.638;
    mela.computeP(pg1g2_pi2_VAJHU_NEW, false);
    pg1g2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
    mela.computeP(pg1g4_VAJHU_NEW, false);
    pg1g4_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=2.521;
    mela.computeP(pg1g4_pi2_VAJHU_NEW, false);
    pg1g4_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=12046.01;
    mela.computeP(pg1g1prime2_VAJHU_NEW, false);
    pg1g1prime2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_VAJHU_NEW, false);
    pg1g1prime2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);


    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p0plus_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
    mela.computeP(p0hplus_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
    mela.computeP(p0minus_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-12046.01;
    mela.computeP(p0_g1prime2_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
    mela.computeP(pg1g2_VAMCFM_NEW, false);
    pg1g2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0hplus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.638;
    mela.computeP(pg1g2_pi2_VAMCFM_NEW, false);
    pg1g2_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0hplus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
    mela.computeP(pg1g4_VAMCFM_NEW, false);
    pg1g4_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0minus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=2.521;
    mela.computeP(pg1g4_pi2_VAMCFM_NEW, false);
    pg1g4_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0minus_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=12046.01;
    mela.computeP(pg1g1prime2_VAMCFM_NEW, false);
    pg1g1prime2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0_g1prime2_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_VAMCFM_NEW, false);
    pg1g1prime2_pi2_VAMCFM_NEW -= (p0plus_VAMCFM_NEW + p0_g1prime2_VAMCFM_NEW);


    p0hplus_VAJHU_ratio_NEW = p0hplus_VAJHU_NEW/p0plus_VAJHU_NEW;
    p0minus_VAJHU_ratio_NEW = p0minus_VAJHU_NEW/p0plus_VAJHU_NEW;
    p0_g1prime2_VAJHU_ratio_NEW = p0_g1prime2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g2_VAJHU_ratio_NEW = pg1g2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g4_VAJHU_ratio_NEW = pg1g4_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g1prime2_VAJHU_ratio_NEW = pg1g1prime2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g2_pi2_VAJHU_ratio_NEW = pg1g2_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g4_pi2_VAJHU_ratio_NEW = pg1g4_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;
    pg1g1prime2_pi2_VAJHU_ratio_NEW = pg1g1prime2_pi2_VAJHU_NEW/p0plus_VAJHU_NEW;

    p0hplus_VAMCFM_ratio_NEW = p0hplus_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    p0minus_VAMCFM_ratio_NEW = p0minus_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    p0_g1prime2_VAMCFM_ratio_NEW = p0_g1prime2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g2_VAMCFM_ratio_NEW = pg1g2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g4_VAMCFM_ratio_NEW = pg1g4_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g1prime2_VAMCFM_ratio_NEW = pg1g1prime2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g2_pi2_VAMCFM_ratio_NEW = pg1g2_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g4_pi2_VAMCFM_ratio_NEW = pg1g4_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;
    pg1g1prime2_pi2_VAMCFM_ratio_NEW = pg1g1prime2_pi2_VAMCFM_NEW/p0plus_VAMCFM_NEW;



    /***** WW *****/
    mela.setCurrentCandidateFromIndex(1);

    if (verbosity>=TVar::DEBUG) cout << "Computing bkg_WW_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZQQB);
    mela.computeP(bkg_WW_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_WW_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgWW, TVar::MCFM, TVar::ZZGG);
    mela.computeP(ggzz_WW_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_WW_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(p0plus_WW_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_WW_VAMCFM_NEW" << endl;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(ggzz_p0plus_WW_VAMCFM_NEW, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing p0plus_WW_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p0plus_WW_VAMCFM_NEW_BSMOn, false);

    if (verbosity>=TVar::DEBUG) cout << "Computing ggzz_p0plus_WW_VAMCFM_NEW_BSMOn" << endl;
    mela.setProcess(TVar::bkgWW_SMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(ggzz_p0plus_WW_VAMCFM_NEW_BSMOn, false);

    mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p0plus_WW_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
    mela.computeP(p0hplus_WW_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
    mela.computeP(p0minus_WW_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-12046.01;
    mela.computeP(p0_g1prime2_WW_VAJHU_NEW, false);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
    mela.computeP(pg1g2_WW_VAJHU_NEW, false);
    pg1g2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0hplus_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.638;
    mela.computeP(pg1g2_pi2_WW_VAJHU_NEW, false);
    pg1g2_pi2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0hplus_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
    mela.computeP(pg1g4_WW_VAJHU_NEW, false);
    pg1g4_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0minus_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=2.521;
    mela.computeP(pg1g4_pi2_WW_VAJHU_NEW, false);
    pg1g4_pi2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0minus_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=12046.01;
    mela.computeP(pg1g1prime2_WW_VAJHU_NEW, false);
    pg1g1prime2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0_g1prime2_WW_VAJHU_NEW);

    mela.selfDHggcoupl[0][gHIGGS_GG_2][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_WW_VAJHU_NEW, false);
    pg1g1prime2_pi2_WW_VAJHU_NEW -= (p0plus_WW_VAJHU_NEW + p0_g1prime2_WW_VAJHU_NEW);


    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.computeP(p0plus_WW_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
    mela.computeP(p0hplus_WW_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
    mela.computeP(p0minus_WW_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=-12046.01;
    mela.computeP(p0_g1prime2_WW_VAMCFM_NEW, false);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][0]=1.638;
    mela.computeP(pg1g2_WW_VAMCFM_NEW, false);
    pg1g2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0hplus_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_2][1]=1.638;
    mela.computeP(pg1g2_pi2_WW_VAMCFM_NEW, false);
    pg1g2_pi2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0hplus_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][0]=2.521;
    mela.computeP(pg1g4_WW_VAMCFM_NEW, false);
    pg1g4_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0minus_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_4][1]=2.521;
    mela.computeP(pg1g4_pi2_WW_VAMCFM_NEW, false);
    pg1g4_pi2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0minus_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][0]=12046.01;
    mela.computeP(pg1g1prime2_WW_VAMCFM_NEW, false);
    pg1g1prime2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0_g1prime2_WW_VAMCFM_NEW);

    mela.selfDHzzcoupl[0][gHIGGS_VV_1][0]=1;
    mela.selfDHzzcoupl[0][gHIGGS_VV_1_PRIME2][1]=12046.01;
    mela.computeP(pg1g1prime2_pi2_WW_VAMCFM_NEW, false);
    pg1g1prime2_pi2_WW_VAMCFM_NEW -= (p0plus_WW_VAMCFM_NEW + p0_g1prime2_WW_VAMCFM_NEW);


    p0hplus_WW_VAJHU_ratio_NEW = p0hplus_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    p0minus_WW_VAJHU_ratio_NEW = p0minus_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    p0_g1prime2_WW_VAJHU_ratio_NEW = p0_g1prime2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g2_WW_VAJHU_ratio_NEW = pg1g2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g4_WW_VAJHU_ratio_NEW = pg1g4_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g1prime2_WW_VAJHU_ratio_NEW = pg1g1prime2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g2_pi2_WW_VAJHU_ratio_NEW = pg1g2_pi2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g4_pi2_WW_VAJHU_ratio_NEW = pg1g4_pi2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;
    pg1g1prime2_pi2_WW_VAJHU_ratio_NEW = pg1g1prime2_pi2_WW_VAJHU_NEW/p0plus_WW_VAJHU_NEW;

    p0hplus_WW_VAMCFM_ratio_NEW = p0hplus_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    p0minus_WW_VAMCFM_ratio_NEW = p0minus_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    p0_g1prime2_WW_VAMCFM_ratio_NEW = p0_g1prime2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g2_WW_VAMCFM_ratio_NEW = pg1g2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g4_WW_VAMCFM_ratio_NEW = pg1g4_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g1prime2_WW_VAMCFM_ratio_NEW = pg1g1prime2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g2_pi2_WW_VAMCFM_ratio_NEW = pg1g2_pi2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g4_pi2_WW_VAMCFM_ratio_NEW = pg1g4_pi2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;
    pg1g1prime2_pi2_WW_VAMCFM_ratio_NEW = pg1g1prime2_pi2_WW_VAMCFM_NEW/p0plus_WW_VAMCFM_NEW;


    mela.resetInputEvent();
    newtree->Fill();
  }


  foutput->WriteTObject(newtree);
  foutput->Close();
  finput->Close();
  mela.setVerbosity(bkpverbosity);
}

void testME_SuperMela_FullSim(int flavor=2, bool useBkgSample=false, bool debug=false, shared_ptr<Mela> melaptr=nullptr){
  int erg_tev=8;
  float mPOLE=125.6;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = (debug ? TVar::DEBUG : TVar::ERROR);
  if (!melaptr) {
    melaptr.reset(new Mela(erg_tev, mPOLE, verbosity));
  }
  Mela& mela = *melaptr;
  TVar::VerbosityLevel bkpverbosity = mela.getVerbosity();
  mela.setVerbosity(verbosity);
  if (verbosity>=TVar::DEBUG) cout << "Mela is initialized" << endl;

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput;
  TFile* foutput;
  if (!useBkgSample){
    finput = new TFile(Form("%s/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root", cinput_main.Data(), (flavor>=2 ? "2mu2e" : "4e")), "read");
    foutput = new TFile(Form("HZZ4lTree_powheg15jhuGenV3-0PMH125.6_%s_OriginalMEv2ValidationTestOnly.root", (flavor>=2 ? "2mu2e" : "4e")), "recreate");
  }
  else{
    finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor>=2 ? "2mu2e" : "4e"), (flavor==2 ? "2e2mu" : "4e")), "read");
    foutput = new TFile(Form("HZZ4lTree_ZZTo%s_OriginalMEv2ValidationTestOnly.root", (flavor>=2 ? "2e2mu" : "4e")), "recreate");
  }

  float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
  float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

  float p0plus_m4l_NEW, p0plus_m4l_ScaleUp_NEW, p0plus_m4l_ScaleDown_NEW, p0plus_m4l_ResUp_NEW, p0plus_m4l_ResDown_NEW;
  float bkg_m4l_NEW, bkg_m4l_ScaleUp_NEW, bkg_m4l_ScaleDown_NEW, bkg_m4l_ResUp_NEW, bkg_m4l_ResDown_NEW;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);
  tree->SetBranchAddress("p0plus_m4l", &p0plus_m4l);
  tree->SetBranchAddress("p0plus_m4l_ScaleUp", &p0plus_m4l_ScaleUp);
  tree->SetBranchAddress("p0plus_m4l_ScaleDown", &p0plus_m4l_ScaleDown);
  tree->SetBranchAddress("p0plus_m4l_ResUp", &p0plus_m4l_ResUp);
  tree->SetBranchAddress("p0plus_m4l_ResDown", &p0plus_m4l_ResDown);
  tree->SetBranchAddress("bkg_m4l", &bkg_m4l);
  tree->SetBranchAddress("bkg_m4l_ScaleUp", &bkg_m4l_ScaleUp);
  tree->SetBranchAddress("bkg_m4l_ScaleDown", &bkg_m4l_ScaleDown);
  tree->SetBranchAddress("bkg_m4l_ResUp", &bkg_m4l_ResUp);
  tree->SetBranchAddress("bkg_m4l_ResDown", &bkg_m4l_ResDown);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("ZZMass", &mzz);
  newtree->Branch("Z1Mass", &m1);
  newtree->Branch("Z2Mass", &m2);
  newtree->Branch("helcosthetaZ1", &h1);
  newtree->Branch("helcosthetaZ2", &h2);
  newtree->Branch("helphi", &phi);
  newtree->Branch("costhetastar", &hs);
  newtree->Branch("phistarZ1", &phi1);

  newtree->Branch("p0plus_m4l", &p0plus_m4l);
  newtree->Branch("p0plus_m4l_ScaleUp", &p0plus_m4l_ScaleUp);
  newtree->Branch("p0plus_m4l_ScaleDown", &p0plus_m4l_ScaleDown);
  newtree->Branch("p0plus_m4l_ResUp", &p0plus_m4l_ResUp);
  newtree->Branch("p0plus_m4l_ResDown", &p0plus_m4l_ResDown);
  newtree->Branch("bkg_m4l", &bkg_m4l);
  newtree->Branch("bkg_m4l_ScaleUp", &bkg_m4l_ScaleUp);
  newtree->Branch("bkg_m4l_ScaleDown", &bkg_m4l_ScaleDown);
  newtree->Branch("bkg_m4l_ResUp", &bkg_m4l_ResUp);
  newtree->Branch("bkg_m4l_ResDown", &bkg_m4l_ResDown);

  newtree->Branch("p0plus_m4l_NEW", &p0plus_m4l_NEW);
  newtree->Branch("p0plus_m4l_ScaleUp_NEW", &p0plus_m4l_ScaleUp_NEW);
  newtree->Branch("p0plus_m4l_ScaleDown_NEW", &p0plus_m4l_ScaleDown_NEW);
  newtree->Branch("p0plus_m4l_ResUp_NEW", &p0plus_m4l_ResUp_NEW);
  newtree->Branch("p0plus_m4l_ResDown_NEW", &p0plus_m4l_ResDown_NEW);
  newtree->Branch("bkg_m4l_NEW", &bkg_m4l_NEW);
  newtree->Branch("bkg_m4l_ScaleUp_NEW", &bkg_m4l_ScaleUp_NEW);
  newtree->Branch("bkg_m4l_ScaleDown_NEW", &bkg_m4l_ScaleDown_NEW);
  newtree->Branch("bkg_m4l_ResUp_NEW", &bkg_m4l_ResUp_NEW);
  newtree->Branch("bkg_m4l_ResDown_NEW", &bkg_m4l_ResDown_NEW);

  float GenLep1Id=0, GenLep2Id=0, GenLep3Id=0, GenLep4Id=0;
  if (flavor == 2){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 1){
    GenLep1Id=11;
    GenLep2Id=-11;
    GenLep3Id=11;
    GenLep4Id=-11;
  }
  else if (flavor == 0){
    GenLep1Id=13;
    GenLep2Id=-13;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 3){
    GenLep1Id=14;
    GenLep2Id=-14;
    GenLep3Id=13;
    GenLep4Id=-13;
  }
  else if (flavor == 4){
    GenLep1Id=0;
    GenLep2Id=-0;
    GenLep3Id=1;
    GenLep4Id=-1;
  }
  else if (flavor == 5){
    GenLep1Id=1;
    GenLep2Id=-1;
    GenLep3Id=2;
    GenLep4Id=-2;
  }
  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (int ev = 0; ev < (verbosity>=TVar::DEBUG ? 1 : (!useBkgSample ? min(1000, (int)tree->GetEntries()) : (int)tree->GetEntries())); ev++){
    tree->GetEntry(ev);
    if (ev%10000 == 0) cout << "Processing event " << ev << endl;

    int idOrdered[4] ={ static_cast<int>(GenLep1Id), static_cast<int>(GenLep2Id), static_cast<int>(GenLep3Id), static_cast<int>(GenLep4Id) };
    TLorentzVector pOrdered[4];
    std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
    for (int ip=0; ip<min(4, (int)daus.size()); ip++) pOrdered[ip]=daus.at(ip);
    SimpleParticleCollection_t daughters_ZZ;
    for (unsigned int idau=0; idau<4; idau++) daughters_ZZ.push_back(SimpleParticle_t(idOrdered[idau], pOrdered[idau]));
    mela.setInputEvent(&daughters_ZZ, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

    /***** ZZ *****/

    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela.computePM4l(TVar::SMSyst_None, p0plus_m4l_NEW);
    mela.computePM4l(TVar::SMSyst_ScaleUp, p0plus_m4l_ScaleUp_NEW);
    mela.computePM4l(TVar::SMSyst_ScaleDown, p0plus_m4l_ScaleDown_NEW);
    mela.computePM4l(TVar::SMSyst_ResUp, p0plus_m4l_ResUp_NEW);
    mela.computePM4l(TVar::SMSyst_ResDown, p0plus_m4l_ResDown_NEW);

    mela.setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
    mela.computePM4l(TVar::SMSyst_None, bkg_m4l_NEW);
    mela.computePM4l(TVar::SMSyst_ScaleUp, bkg_m4l_ScaleUp_NEW);
    mela.computePM4l(TVar::SMSyst_ScaleDown, bkg_m4l_ScaleDown_NEW);
    mela.computePM4l(TVar::SMSyst_ResUp, bkg_m4l_ResUp_NEW);
    mela.computePM4l(TVar::SMSyst_ResDown, bkg_m4l_ResDown_NEW);

    mela.resetInputEvent();
    newtree->Fill();
  }

  foutput->WriteTObject(newtree);
  foutput->Close();
  finput->Close();
  mela.setVerbosity(bkpverbosity);
}
