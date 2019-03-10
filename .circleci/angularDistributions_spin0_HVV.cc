#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "Mela.h"

#include "RooBinning.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDecay.h"
#include "RooGaussModel.h"
#include "RooPlot.h"
#include "RooRealVar.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>


using namespace RooFit;
using namespace std;


bool checkListVariable(const vector<string>& list, const string& var);
void splitOption(const string rawoption, string& wish, string& value, char delimiter);
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter);
void extractCoupling(string opt, ScalarPdfFactory_HVV& factory);

void angularDistributions_spin0_HVV(string cinput, string strcouplings, string decaytype, string coutput, double mPOLE = 125., int nbins=80, bool useTaus=false){
  vector<string> strcouplingsList;
  splitOptionRecursive(strcouplings, strcouplingsList, ';');

  RooRealVar* mzz = new RooRealVar("GenHMass", "M_{ZZ} (GeV)", mPOLE, mPOLE-0.02, mPOLE+0.02);
  RooRealVar* z1mass = new RooRealVar("GenZ1Mass", "m_{Z1} (GeV)", 0.0, min(120., mPOLE));
  RooRealVar* z2mass = new RooRealVar("GenZ2Mass", "m_{Z2} (GeV)", 0.0, min(120., (mPOLE-90.)*20./35.+55.));
  RooRealVar* hs = new RooRealVar("Gencosthetastar", "cos#theta^{*}", -1, 1);
  RooRealVar* h1 = new RooRealVar("GenhelcosthetaZ1", "cos#theta_{Z1}", -1, 1);
  RooRealVar* h2 = new RooRealVar("GenhelcosthetaZ2", "cos#theta_{Z2}", -1, 1);
  RooRealVar* Phi = new RooRealVar("Genhelphi", "#Phi", -TMath::Pi(), TMath::Pi());
  RooRealVar* Phi1 = new RooRealVar("GenphistarZ1", "#Phi_{Z1}", -TMath::Pi(), TMath::Pi());
  RooRealVar* Y = new RooRealVar("GenY", "Y", 0);

  RooSpin::modelMeasurables measurables_;
  measurables_.h1 = h1;
  measurables_.h2 = h2;
  measurables_.Phi = Phi;
  measurables_.m1 = z1mass;
  measurables_.m2 = z2mass;
  measurables_.m12 = mzz;
  measurables_.hs = hs;
  measurables_.Phi1 = Phi1;
  measurables_.Y = Y;

  const int nVars = 8;
  const int nPlots=nVars-1;
  RooArgSet treeargs;
  RooRealVar* measurables[nPlots]={ hs, Phi1, z1mass, h1, z2mass, h2, Phi };

  float kd_vars[nVars];
  TString strKDs[nVars]={
    "GenHMass", "GenZ1Mass", "GenZ2Mass",
    "GenhelcosthetaZ1", "GenhelcosthetaZ2", "GenphistarZ1",
    "Gencosthetastar", "Genhelphi"
  };
  TString strGenLepId[4]={ "GenLep1Id", "GenLep2Id", "GenLep3Id", "GenLep4Id" };
  int GenLepId[4]={ 0 };
  RooSpin::VdecayType Vdecay1=(decaytype.find("WW")!=string::npos ? RooSpin::kVdecayType_Wany : RooSpin::kVdecayType_Zll);
  RooSpin::VdecayType Vdecay2=(decaytype.find("WW")!=string::npos ? RooSpin::kVdecayType_Wany : RooSpin::kVdecayType_Zll);
  int nToPlot=nPlots;
  if (decaytype.find("ZG")!=string::npos){
    nToPlot-=3;
    Vdecay2=RooSpin::kVdecayType_GammaOnshell;
    if (decaytype.find("2l")!=string::npos) Vdecay1=RooSpin::kVdecayType_Zll;
    else if (decaytype.find("2nu")!=string::npos) Vdecay1=RooSpin::kVdecayType_Znn;
    else if (decaytype.find("2q")!=string::npos) Vdecay1=RooSpin::kVdecayType_Zud;
    else{
      cerr << "Could not find the Z decays! Exiting." << endl;
      return;
    }
  }
  else if (decaytype.find("GG")!=string::npos && decaytype.find("GGto")==string::npos){ nToPlot-=6; Vdecay1=RooSpin::kVdecayType_GammaOnshell; Vdecay2=RooSpin::kVdecayType_GammaOnshell; }
  else if (decaytype.find("ZZ")!=string::npos){
    if (decaytype.find("4l")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Zll; }
    else if (decaytype.find("4nu")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Znn; Vdecay2=RooSpin::kVdecayType_Znn; }
    else if (decaytype.find("4q")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zud; Vdecay2=RooSpin::kVdecayType_Zud; }
    else if (decaytype.find("2l2q")!=string::npos || decaytype.find("2q2l")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Zud; }
    else if (decaytype.find("2l2nu")!=string::npos || decaytype.find("2nu2l")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Znn; }
    else if (decaytype.find("2q2nu")!=string::npos || decaytype.find("2nu2q")!=string::npos){ Vdecay1=RooSpin::kVdecayType_Zll; Vdecay2=RooSpin::kVdecayType_Zud; }
    else{
      cerr << "Could not find the Z decays! Exiting." << endl;
      return;
    }
  }
  else if (decaytype.find("WW")==string::npos){
    cerr << "Could not find the V decays! Exiting." << endl;
    return;
  }
  for (int r=0; r<nToPlot; r++) treeargs.add(*(measurables[r]));
  if (Vdecay1==RooSpin::kVdecayType_GammaOnshell){ z1mass->setVal(0); z1mass->setConstant(true); }
  if (Vdecay2==RooSpin::kVdecayType_GammaOnshell){ z2mass->setVal(0); z2mass->setConstant(true); }

  cout << "Decay modes found are " << Vdecay1 << '\t' << Vdecay2 << endl;
  ScalarPdfFactory_HVV* someHiggs = new ScalarPdfFactory_HVV(measurables_, false, Vdecay1, Vdecay2);
  someHiggs->makeParamsConst(false);
  someHiggs->makeCouplingsConst(false);
  for (auto& opt:strcouplingsList) extractCoupling(opt, *someHiggs);
  someHiggs->makeParamsConst(true);
  someHiggs->makeCouplingsConst(true);
  someHiggs->getPDF()->printParameters();

  if (coutput==""){
    size_t lastSlash = cinput.find_last_of("/\\");
    string finName;
    finName=cinput;
    finName.resize(lastSlash+1);
    coutput = finName + "Validation/";
  }
  cout << "Output folder is " << coutput << endl;
  string strCmd = "mkdir -p ";
  strCmd.append(coutput);
  gSystem->Exec(strCmd.c_str());

  TChain* tree = new TChain("SelectedTree");
  tree->Add(cinput.c_str());
  TTree* reducedTree = new TTree("ReducedTree", "");
  for (int v=0; v<8; v++){
    tree->SetBranchAddress(strKDs[v], (kd_vars+v));
    reducedTree->Branch(strKDs[v], (kd_vars+v));
  }
  for (int v=0; v<4; v++) tree->SetBranchAddress(strGenLepId[v], (GenLepId+v));

  int nRecorded=0;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if ((kd_vars[0]-mPOLE)>=0.02) continue;
    if (
      ((GenLepId[0]==11 || GenLepId[0]==13 || GenLepId[0]==15) && (GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Zll))
      ||
      ((GenLepId[0]==12 || GenLepId[0]==14 || GenLepId[0]==16) && (GenLepId[2]==12 || GenLepId[2]==14 || GenLepId[2]==16) && (Vdecay1==RooSpin::kVdecayType_Znn && Vdecay2==RooSpin::kVdecayType_Znn))
      ||
      ((GenLepId[0]>0 && GenLepId[0]<6) && (GenLepId[2]>0 && GenLepId[2]<6) && (Vdecay1==RooSpin::kVdecayType_Zud && Vdecay2==RooSpin::kVdecayType_Zud))
      ||
      ((((GenLepId[0]==11 || GenLepId[0]==13 || GenLepId[0]==15) && (GenLepId[2]==12 || GenLepId[2]==14 || GenLepId[2]==16)) || ((GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (GenLepId[0]==12 || GenLepId[0]==14 || GenLepId[0]==16))) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Znn))
      ||
      ((((GenLepId[0]==11 || GenLepId[0]==13 || GenLepId[0]==15) && (GenLepId[2]>0 && GenLepId[2]<6)) || ((GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (GenLepId[0]>0 && GenLepId[0]<6))) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Zud))
      ||
      ((((GenLepId[0]==12 || GenLepId[0]==14 || GenLepId[0]==16) && (GenLepId[2]>0 && GenLepId[2]<6)) || ((GenLepId[2]==12 || GenLepId[2]==14 || GenLepId[2]==16) && (GenLepId[0]>0 && GenLepId[0]<6))) && (Vdecay1==RooSpin::kVdecayType_Znn && Vdecay2==RooSpin::kVdecayType_Zud))
      ||
      (Vdecay1==RooSpin::kVdecayType_Wany && Vdecay2==RooSpin::kVdecayType_Wany)
      ||
      (Vdecay1==RooSpin::kVdecayType_GammaOnshell || Vdecay2==RooSpin::kVdecayType_GammaOnshell)
      ){
      if (
        ((GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (GenLepId[0]==12 || GenLepId[0]==14 || GenLepId[0]==16) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Znn))
        ||
        ((GenLepId[2]==11 || GenLepId[2]==13 || GenLepId[2]==15) && (GenLepId[0]>0 && GenLepId[0]<6) && (Vdecay1==RooSpin::kVdecayType_Zll && Vdecay2==RooSpin::kVdecayType_Zud))
        ||
        ((GenLepId[2]==12 || GenLepId[2]==14 || GenLepId[2]==16) && (GenLepId[0]>0 && GenLepId[0]<6) && (Vdecay1==RooSpin::kVdecayType_Znn && Vdecay2==RooSpin::kVdecayType_Zud))
        ){
        //swap(strKDs[2], strKDs[2]);
        //swap(strKDs[3], strKDs[4]);
        continue; // Don't bother to recalculate angles etc.
      }
      if (GenLepId[0]==GenLepId[2] && GenLepId[1]==GenLepId[3] && Vdecay1!=RooSpin::kVdecayType_GammaOnshell && Vdecay2!=RooSpin::kVdecayType_GammaOnshell) continue;
      if (!useTaus && (GenLepId[0]==15 || GenLepId[1]==-15 || GenLepId[2]==15 || GenLepId[3]==-15)) continue;
      reducedTree->Fill();
      nRecorded++;
    }
  }
  cout << "Number of entries recorded: " << nRecorded << '/' << tree->GetEntries() << endl;

  RooDataSet* dataSM = new RooDataSet("data", "data", reducedTree, treeargs);
  for (int plotIndex=0; plotIndex<(int)treeargs.getSize(); plotIndex++){
    cout << plotIndex << endl;

    string m_name = measurables[plotIndex]->GetName();
    int nbins_plot=nbins;
    if (m_name.find("Mass")!=string::npos) nbins_plot*=2;
    RooPlot* plot = measurables[plotIndex]->frame(nbins_plot);
    plot->GetXaxis()->CenterTitle();
    plot->GetYaxis()->SetTitleOffset(1.2);
    plot->GetYaxis()->CenterTitle();
    plot->GetYaxis()->SetTitle("Number of Events");
    plot->GetXaxis()->SetNdivisions(-505);
    plot->SetTitle(m_name.c_str());

    dataSM->plotOn(plot, MarkerColor(kRed), MarkerStyle(3), MarkerSize(1.2), LineWidth(0), XErrorSize(0), DataError(RooAbsData::Poisson));
    RooSpinZero_7DComplex_withAccep_HVV* pdf = (RooSpinZero_7DComplex_withAccep_HVV*)someHiggs->getPDF();
    pdf->plotOn(plot, LineColor(kRed), LineWidth(2));

    TGaxis::SetMaxDigits(3);

    TCanvas* can = new TCanvas("can", "can", 600, 600);

    plot->Draw();

    char cname[200];
    sprintf(cname, "%s", m_name.c_str());
    string cname_pdf=cname;
    string cname_eps=cname;
    string cname_png=cname;
    cname_pdf = (coutput + cname_pdf) + ".pdf";
    cname_eps = (coutput + cname_eps) + ".eps";
    cname_png = (coutput + cname_png) + ".png";

    can->SaveAs(cname_pdf.c_str());
    can->SaveAs(cname_eps.c_str());
    can->SaveAs(cname_png.c_str());
    can->Close();
  }

  delete dataSM;
  delete reducedTree;
  delete tree;
  delete someHiggs;
  delete mzz;
  delete z1mass;
  delete z2mass;
  delete hs;
  delete h1;
  delete h2;
  delete Phi;
  delete Phi1;
  delete Y;
}

void extractCoupling(string opt, ScalarPdfFactory_HVV& factory){
  RooSpinZero::modelCouplings& couplings = factory.couplings;
  RooSpin::modelParameters& parameters = factory.parameters;

  string wish, strVal, strValRe, strValIm;
  // Use double precision for couplings
  Double_t valRe=0;
  Double_t valIm=0;
  splitOption(opt, wish, strVal, '=');
  // Lambda and cz/cw couplings have no imaginary components, so do not expect to parse them with ','.
  if (
    wish.find("Lambda")==string::npos
    &&
    wish.find("q1sq")==string::npos
    &&
    wish.find("q2sq")==string::npos
    &&
    wish.find("q12sq")==string::npos
    &&
    wish.find("mW")==string::npos
    &&
    wish.find("mZ")==string::npos
    &&
    wish.find("gVprimeff_decay")==string::npos
    ){
    splitOption(strVal, strValRe, strValIm, ',');
    valRe = atof(strValRe.c_str());
    valIm = atof(strValIm.c_str());
  }
  else valRe = atof(strVal.c_str());

  bool setanomalous, setghv1;

  // Here we go again, sillions of couplings
  if (wish=="cv_q1sq"){ ((RooRealVar*)couplings.cLambda_qsq[cLambdaHIGGS_VV_QSQ1])->setVal((int) valRe); }
  else if (wish=="cv_q2sq"){ ((RooRealVar*)couplings.cLambda_qsq[cLambdaHIGGS_VV_QSQ2])->setVal((int) valRe); }
  else if (wish=="cv_q12sq"){ ((RooRealVar*)couplings.cLambda_qsq[cLambdaHIGGS_VV_QSQ12])->setVal((int) valRe); }

  else if (wish=="Lambda_v11"){ ((RooRealVar*) couplings.Lambda_z1qsq[cLambdaHIGGS_VV_QSQ1])->setVal(valRe); }
  else if (wish=="Lambda_v12"){ ((RooRealVar*) couplings.Lambda_z1qsq[cLambdaHIGGS_VV_QSQ2])->setVal(valRe); }
  else if (wish=="Lambda_v10"){ ((RooRealVar*) couplings.Lambda_z1qsq[cLambdaHIGGS_VV_QSQ12])->setVal(valRe); }

  else if (wish=="Lambda_v21"){ ((RooRealVar*) couplings.Lambda_z2qsq[cLambdaHIGGS_VV_QSQ1])->setVal(valRe); }
  else if (wish=="Lambda_v22"){ ((RooRealVar*) couplings.Lambda_z2qsq[cLambdaHIGGS_VV_QSQ2])->setVal(valRe); }
  else if (wish=="Lambda_v20"){ ((RooRealVar*) couplings.Lambda_z2qsq[cLambdaHIGGS_VV_QSQ12])->setVal(valRe); }

  else if (wish=="Lambda_v31"){ ((RooRealVar*) couplings.Lambda_z3qsq[cLambdaHIGGS_VV_QSQ1])->setVal(valRe); }
  else if (wish=="Lambda_v32"){ ((RooRealVar*) couplings.Lambda_z3qsq[cLambdaHIGGS_VV_QSQ2])->setVal(valRe); }
  else if (wish=="Lambda_v30"){ ((RooRealVar*) couplings.Lambda_z3qsq[cLambdaHIGGS_VV_QSQ12])->setVal(valRe); }

  else if (wish=="Lambda_v41"){ ((RooRealVar*) couplings.Lambda_z4qsq[cLambdaHIGGS_VV_QSQ1])->setVal(valRe); }
  else if (wish=="Lambda_v42"){ ((RooRealVar*) couplings.Lambda_z4qsq[cLambdaHIGGS_VV_QSQ2])->setVal(valRe); }
  else if (wish=="Lambda_v40"){ ((RooRealVar*) couplings.Lambda_z4qsq[cLambdaHIGGS_VV_QSQ12])->setVal(valRe); }

  else if (wish=="ghv1"){ ((RooRealVar*) couplings.g1List[0][0])->setVal(valRe); ((RooRealVar*) couplings.g1List[0][1])->setVal(valIm); setghv1 = true; }
  else if (wish=="ghv1_prime"){ ((RooRealVar*) couplings.g1List[1][0])->setVal(valRe); ((RooRealVar*) couplings.g1List[1][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv1_prime2"){ ((RooRealVar*) couplings.g1List[2][0])->setVal(valRe); ((RooRealVar*) couplings.g1List[2][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv1_prime3"){ ((RooRealVar*) couplings.g1List[3][0])->setVal(valRe); ((RooRealVar*) couplings.g1List[3][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv1_prime4"){ ((RooRealVar*) couplings.g1List[4][0])->setVal(valRe); ((RooRealVar*) couplings.g1List[4][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv1_prime5"){ ((RooRealVar*) couplings.g1List[5][0])->setVal(valRe); ((RooRealVar*) couplings.g1List[5][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv1_prime6"){ ((RooRealVar*) couplings.g1List[6][0])->setVal(valRe); ((RooRealVar*) couplings.g1List[6][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv1_prime7"){ ((RooRealVar*) couplings.g1List[7][0])->setVal(valRe); ((RooRealVar*) couplings.g1List[7][1])->setVal(valIm); setanomalous = true; }

  else if (wish=="ghv2"){ ((RooRealVar*) couplings.g2List[0][0])->setVal(valRe); ((RooRealVar*) couplings.g2List[0][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv2_prime"){ ((RooRealVar*) couplings.g2List[1][0])->setVal(valRe); ((RooRealVar*) couplings.g2List[1][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv2_prime2"){ ((RooRealVar*) couplings.g2List[2][0])->setVal(valRe); ((RooRealVar*) couplings.g2List[2][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv2_prime3"){ ((RooRealVar*) couplings.g2List[3][0])->setVal(valRe); ((RooRealVar*) couplings.g2List[3][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv2_prime4"){ ((RooRealVar*) couplings.g2List[4][0])->setVal(valRe); ((RooRealVar*) couplings.g2List[4][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv2_prime5"){ ((RooRealVar*) couplings.g2List[5][0])->setVal(valRe); ((RooRealVar*) couplings.g2List[5][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv2_prime6"){ ((RooRealVar*) couplings.g2List[6][0])->setVal(valRe); ((RooRealVar*) couplings.g2List[6][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv2_prime7"){ ((RooRealVar*) couplings.g2List[7][0])->setVal(valRe); ((RooRealVar*) couplings.g2List[7][1])->setVal(valIm); setanomalous = true; }

  else if (wish=="ghv3"){ ((RooRealVar*) couplings.g3List[0][0])->setVal(valRe); ((RooRealVar*) couplings.g3List[0][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv3_prime"){ ((RooRealVar*) couplings.g3List[1][0])->setVal(valRe); ((RooRealVar*) couplings.g3List[1][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv3_prime2"){ ((RooRealVar*) couplings.g3List[2][0])->setVal(valRe); ((RooRealVar*) couplings.g3List[2][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv3_prime3"){ ((RooRealVar*) couplings.g3List[3][0])->setVal(valRe); ((RooRealVar*) couplings.g3List[3][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv3_prime4"){ ((RooRealVar*) couplings.g3List[4][0])->setVal(valRe); ((RooRealVar*) couplings.g3List[4][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv3_prime5"){ ((RooRealVar*) couplings.g3List[5][0])->setVal(valRe); ((RooRealVar*) couplings.g3List[5][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv3_prime6"){ ((RooRealVar*) couplings.g3List[6][0])->setVal(valRe); ((RooRealVar*) couplings.g3List[6][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv3_prime7"){ ((RooRealVar*) couplings.g3List[7][0])->setVal(valRe); ((RooRealVar*) couplings.g3List[7][1])->setVal(valIm); setanomalous = true; }

  else if (wish=="ghv4"){ ((RooRealVar*) couplings.g4List[0][0])->setVal(valRe); ((RooRealVar*) couplings.g4List[0][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv4_prime"){ ((RooRealVar*) couplings.g4List[1][0])->setVal(valRe); ((RooRealVar*) couplings.g4List[1][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv4_prime2"){ ((RooRealVar*) couplings.g4List[2][0])->setVal(valRe); ((RooRealVar*) couplings.g4List[2][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv4_prime3"){ ((RooRealVar*) couplings.g4List[3][0])->setVal(valRe); ((RooRealVar*) couplings.g4List[3][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv4_prime4"){ ((RooRealVar*) couplings.g4List[4][0])->setVal(valRe); ((RooRealVar*) couplings.g4List[4][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv4_prime5"){ ((RooRealVar*) couplings.g4List[5][0])->setVal(valRe); ((RooRealVar*) couplings.g4List[5][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv4_prime6"){ ((RooRealVar*) couplings.g4List[6][0])->setVal(valRe); ((RooRealVar*) couplings.g4List[6][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghv4_prime7"){ ((RooRealVar*) couplings.g4List[7][0])->setVal(valRe); ((RooRealVar*) couplings.g4List[7][1])->setVal(valIm); setanomalous = true; }

  else if (wish=="ghzgs1_prime2"){ ((RooRealVar*) couplings.gzgs1List[0][0])->setVal(valRe); ((RooRealVar*) couplings.gzgs1List[0][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghzgs2"){ ((RooRealVar*) couplings.gzgs2List[0][0])->setVal(valRe); ((RooRealVar*) couplings.gzgs2List[0][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghzgs3"){ ((RooRealVar*) couplings.gzgs3List[0][0])->setVal(valRe); ((RooRealVar*) couplings.gzgs3List[0][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghzgs4"){ ((RooRealVar*) couplings.gzgs4List[0][0])->setVal(valRe); ((RooRealVar*) couplings.gzgs4List[0][1])->setVal(valIm); setanomalous = true; }

  else if (wish=="ghgsgs2"){ ((RooRealVar*) couplings.ggsgs2List[0][0])->setVal(valRe); ((RooRealVar*) couplings.ggsgs2List[0][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghgsgs3"){ ((RooRealVar*) couplings.ggsgs3List[0][0])->setVal(valRe); ((RooRealVar*) couplings.ggsgs3List[0][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghgsgs4"){ ((RooRealVar*) couplings.ggsgs4List[0][0])->setVal(valRe); ((RooRealVar*) couplings.ggsgs4List[0][1])->setVal(valIm); setanomalous = true; }

  else if (wish=="ghvvp1"){ ((RooRealVar*) couplings.gvvp1List[0][0])->setVal(valRe); ((RooRealVar*) couplings.gvvp1List[0][1])->setVal(valIm); setanomalous = true; }
  else if (wish=="ghvpvp1"){ ((RooRealVar*) couplings.gvpvp1List[0][0])->setVal(valRe); ((RooRealVar*) couplings.gvpvp1List[0][1])->setVal(valIm); setanomalous = true; }

  else if (wish=="mW") ((RooRealVar*) parameters.mW)->setVal(valRe);
  else if (wish=="mZ") ((RooRealVar*) parameters.mZ)->setVal(valRe);
  else if (wish=="mWprime") ((RooRealVar*) parameters.mWprime)->setVal(valRe);
  else if (wish=="mZprime") ((RooRealVar*) parameters.mZprime)->setVal(valRe);
  else if (wish=="gamW") ((RooRealVar*) parameters.gamW)->setVal(valRe);
  else if (wish=="gamZ") ((RooRealVar*) parameters.gamZ)->setVal(valRe);
  else if (wish=="gamWprime") ((RooRealVar*) parameters.gamWprime)->setVal(valRe);
  else if (wish=="gamZprime") ((RooRealVar*) parameters.gamZprime)->setVal(valRe);

  else if (wish=="gVprimeff_decay1_left") ((RooRealVar*) parameters.gVprimeff_decay1_left)->setVal(valRe);
  else if (wish=="gVprimeff_decay1_right") ((RooRealVar*) parameters.gVprimeff_decay1_right)->setVal(valRe);
  else if (wish=="gVprimeff_decay2_left") ((RooRealVar*) parameters.gVprimeff_decay2_left)->setVal(valRe);
  else if (wish=="gVprimeff_decay2_right") ((RooRealVar*) parameters.gVprimeff_decay2_right)->setVal(valRe);

  else cerr << "extractCoupling: Coupling " << wish << " is not supported!" << endl;

  if (setanomalous && !setghv1){
    cerr << "set an anomalous HVV coupling, but didn't set ghv1!" << endl;
    assert(0);
  }
}

void splitOption(const string rawoption, string& wish, string& value, char delimiter){
  size_t posEq = rawoption.find(delimiter);
  if (posEq!=string::npos){
    wish=rawoption;
    value=rawoption.substr(posEq+1);
    wish.erase(wish.begin()+posEq, wish.end());
  }
  else{
    wish="";
    value=rawoption;
  }
}
void splitOptionRecursive(const string rawoption, vector<string>& splitoptions, char delimiter){
  string suboption=rawoption, result=rawoption;
  string remnant;
  while (result!=""){
    splitOption(suboption, result, remnant, delimiter);
    if (result!="" && !checkListVariable(splitoptions, result)) splitoptions.push_back(result);
    suboption = remnant;
  }
  if (remnant!="") splitoptions.push_back(remnant);
}
bool checkListVariable(const vector<string>& list, const string& var){
  for (unsigned int v=0; v<list.size(); v++){
    if (list.at(v)==var) return true; // Look for exact match
  }
  return false;
}
