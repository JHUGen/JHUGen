#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include "TROOT.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TF1.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TRandom.h"
#include "RooNumIntConfig.h"
#include "RooRealIntegral.h"
#include "Mela.h"
#include "ScalarPdfFactory_VH.h"


using namespace std;
using namespace RooFit;


TString inputdir_7TeV = "/work-zfs/lhc/ianderso/hep/CJLST/140519/PRODFSR";
TString inputdir_8TeV = "/work-zfs/lhc/ianderso/hep/CJLST/140519b/PRODFSR_8TeV";
TString inputdir_13TeV = "/work-zfs/lhc/CJLSTtrees/180214_2016MC";


template<typename T> void addByLowest(std::vector<T>& valArray, T val, bool unique){
  bool inserted = false;
  if (unique){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it==val){
        inserted=true;
        break;
      }
    }
  }
  if (!inserted){
    for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if (*it>=val){
        inserted=true;
        valArray.insert(it, val);
        break;
      }
    }
  }
  if (!inserted) valArray.push_back(val);
}

template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, T val, U index){
  bool inserted = false;
  for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if ((*it).first>=val){
      inserted=true;
      if ((*it).second!=index) valArray.insert(it, std::pair<T, U>(val, index));
      break;
    }
  }
  if (!inserted) valArray.push_back(std::pair<T, U>(val, index));
}

template<typename T> void appendVector(std::vector<T>& a, std::vector<T>& b){ a.insert(a.end(), b.begin(), b.end()); }

template<typename T> bool checkListVariable(const vector<T>& list, const T& var){
  for (unsigned int v=0; v<list.size(); v++){
    if (list.at(v)==var) return true; // Look for exact match
  }
  return false;
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

struct ExtBin{
  double binlow;
  double binhigh;
  vector<int> events;
  vector<float> masses;
  vector<float> mevals;
  vector<float> me2vals;
  vector<float> weights;

  void addEvent(float mass, float me, float me2, float weight=1){
    masses.push_back(mass);
    mevals.push_back(me);
    me2vals.push_back(me2);
    weights.push_back(weight);
  }
  void sift(){
    if (masses.size()==0) return;

    vector<int> take_out;
    vector<pair<float, int>> me_entry[2];
    for (unsigned int ev=0; ev<masses.size(); ev++){
      addByLowest<float, int>(me_entry[0], mevals.at(ev), ev);
      addByLowest<float, int>(me_entry[1], me2vals.at(ev), ev);
    }
    for (unsigned int im=0; im<2; im++){
      int at99p8ev = (float(me_entry[im].size()))*0.998;
      if (at99p8ev==(int)me_entry[im].size()) at99p8ev--;
      int bin=me_entry[im].size()-1;
      while (bin>at99p8ev){
        if (
          me_entry[im].at(at99p8ev).first*2.<me_entry[im].at(bin).first
          ) addByLowest<int>(take_out, me_entry[im].at(bin).second, true);
        bin--;
      }
    }
    for (int bin=take_out.size()-1; bin>=0; bin--){
      int t_ev = take_out.at(bin);
      cout << "Discarding event " << events.at(t_ev) << endl;
      vector<int>::iterator itev;
      vector<float>::iterator it;
      itev=events.begin()+t_ev; events.erase(itev);
      it=masses.begin()+t_ev; masses.erase(it);
      it=mevals.begin()+t_ev; mevals.erase(it);
      it=me2vals.begin()+t_ev; me2vals.erase(it);
      it=weights.begin()+t_ev; weights.erase(it);
    }
  }
  void adjustWeights(float refth=3.){
    if (masses.size()==0) return;

    vector<pair<float, int>> weight_entry;
    for (unsigned int ev=0; ev<masses.size(); ev++) addByLowest<float, int>(weight_entry, weights.at(ev), ev);
    int bin=weight_entry.size()-1;

    int at99p0ev = (float(weight_entry.size()))*0.99;
    if (at99p0ev==(int)weight_entry.size()) at99p0ev--;
    float threshold = weight_entry.at(at99p0ev).first;
    while (bin>at99p0ev){
      if (
        threshold*2.<weight_entry.at(bin).first
        ){
        weights.at(weight_entry.at(bin).second) = pow(threshold, 2)/weight_entry.at(bin).first;
        cout << "Adjusting weight for event " << bin << endl;
      }
      bin--;
    }
    float sum[2]={ 0 };
    for (unsigned int ev=0; ev<weights.size(); ev++){
      sum[0] += weights.at(ev);
      sum[1] += 1;
    }
    for (unsigned int ev=0; ev<weights.size(); ev++) weights.at(ev) *= sum[1]/sum[0];
    if (refth>0.){
      for (unsigned int ev=0; ev<weights.size(); ev++){
        if (weights.at(ev)<=refth) continue;
        weights.at(ev) = pow(refth, 2)/weights.at(ev);
      }
    }
  }
  void mergeBin(const ExtBin& other){
    if (
      other.events.size()!=other.masses.size()
      ||
      other.events.size()!=other.mevals.size()
      ||
      other.events.size()!=other.me2vals.size()
      ||
      other.events.size()!=other.weights.size()
      ){
      cerr << "Size of the following elements are not the same:" << endl;
      cerr << "\t- Events: " << other.events.size() << endl;
      cerr << "\t- Masses: " << other.masses.size() << endl;
      cerr << "\t- MEs: " << other.mevals.size() << endl;
      cerr << "\t- MEs (2): " << other.me2vals.size() << endl;
      cerr << "\t- Weights: " << other.weights.size() << endl;
    }
    for (unsigned int ev=0; ev<other.events.size(); ev++){
      events.push_back(other.events.at(ev));
      addEvent(other.masses.at(ev), other.mevals.at(ev), other.me2vals.at(ev), other.weights.at(ev));
    }
  }
};

float findPoleMass(TString samplename){
  float mass = -1;
  string strtmp = samplename.Data();
  std::size_t extpos = strtmp.find(".root");
  if (extpos!=string::npos) strtmp.erase(extpos, 5);
  vector<string> strsplit;
  splitOptionRecursive(strtmp, strsplit, 'H');
  if (strsplit.size()>1){
    string strmass = strsplit.at(1);
    strsplit.clear();
    splitOptionRecursive(strmass, strsplit, '_');
    strmass = strsplit.at(0);
    mass = std::stod(strmass);
  }
  return mass;
}
TTree* findTree(vector<TTree*> treeList, int evid){
  int ev = evid;
  for (unsigned int t=0; t<treeList.size(); t++){
    TTree* tree = treeList.at(t);
    int nevts = tree->GetEntries();
    if (ev<nevts) return tree;
    else ev -= nevts;
    if (ev<0) cerr << "findTree::ERROR: Could not find the event " << evid << endl;
  }
  return 0;
}
void getEntry(vector<TTree*> treeList, int evid){
  int ev = evid;
  for (unsigned int t=0; t<treeList.size(); t++){
    TTree* tree = treeList.at(t);
    if (tree==0) cerr << "Something went wrong in getEntry. Tree = " << tree << endl;
    int nevts = tree->GetEntries();
    if (ev<nevts){ tree->GetEntry(ev); break; }
    else ev -= nevts;
    if (ev<0) cerr << "findTree::ERROR: Could not find the event " << evid << endl;
  }
}

vector<TString> constructSamplesList(TString strsample, float sqrts){
  vector<TString> samples;
  if (strsample=="JJVBF"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_VBFH116.root");
      samples.push_back("HZZ4lTree_VBFH117.root");
      samples.push_back("HZZ4lTree_VBFH118.root");
      samples.push_back("HZZ4lTree_VBFH119.root");
      samples.push_back("HZZ4lTree_VBFH120.root");
      samples.push_back("HZZ4lTree_VBFH121.root");
      samples.push_back("HZZ4lTree_VBFH122.root");
      samples.push_back("HZZ4lTree_VBFH123.root");
      samples.push_back("HZZ4lTree_VBFH124.root");
      samples.push_back("HZZ4lTree_VBFH125.root");
      samples.push_back("HZZ4lTree_VBFH126.root");
      samples.push_back("HZZ4lTree_VBFH127.root");
      samples.push_back("HZZ4lTree_VBFH128.root");
      samples.push_back("HZZ4lTree_VBFH129.root");
      samples.push_back("HZZ4lTree_VBFH130.root");
      samples.push_back("HZZ4lTree_VBFH135.root");
      samples.push_back("HZZ4lTree_VBFH140.root");
      samples.push_back("HZZ4lTree_VBFH145.root");
      samples.push_back("HZZ4lTree_VBFH150.root");
      samples.push_back("HZZ4lTree_VBFH160.root");
      samples.push_back("HZZ4lTree_VBFH170.root");
      samples.push_back("HZZ4lTree_VBFH180.root");
      samples.push_back("HZZ4lTree_VBFH190.root");
      samples.push_back("HZZ4lTree_powheg15VBFH200.root");
      samples.push_back("HZZ4lTree_powheg15VBFH225.root");
      samples.push_back("HZZ4lTree_powheg15VBFH250.root");
      samples.push_back("HZZ4lTree_powheg15VBFH275.root");
      samples.push_back("HZZ4lTree_powheg15VBFH300.root");
      samples.push_back("HZZ4lTree_powheg15VBFH350.root");
      samples.push_back("HZZ4lTree_powheg15VBFH400.root");
      samples.push_back("HZZ4lTree_powheg15VBFH450.root");
      samples.push_back("HZZ4lTree_powheg15VBFH500.root");
      samples.push_back("HZZ4lTree_powheg15VBFH550.root");
      samples.push_back("HZZ4lTree_powheg15VBFH600.root");
      samples.push_back("HZZ4lTree_powheg15VBFH650.root");
      samples.push_back("HZZ4lTree_powheg15VBFH700.root");
      samples.push_back("HZZ4lTree_powheg15VBFH750.root");
      samples.push_back("HZZ4lTree_powheg15VBFH800.root");
      samples.push_back("HZZ4lTree_powheg15VBFH850.root");
      samples.push_back("HZZ4lTree_powheg15VBFH900.root");
      samples.push_back("HZZ4lTree_powheg15VBFH950.root");
      samples.push_back("HZZ4lTree_powheg15VBFH1000.root");
    }
    else{
      samples.push_back("VBFH115");
      samples.push_back("VBFH120");
      samples.push_back("VBFH124");
      samples.push_back("VBFH125");
      samples.push_back("VBFH126");
      samples.push_back("VBFH130");
      samples.push_back("VBFH135");
      samples.push_back("VBFH140");
      samples.push_back("VBFH145");
      samples.push_back("VBFH150");
      samples.push_back("VBFH155");
      samples.push_back("VBFH160");
      samples.push_back("VBFH165");
      samples.push_back("VBFH170");
      samples.push_back("VBFH175");
      samples.push_back("VBFH180");
      samples.push_back("VBFH190");
      samples.push_back("VBFH200");
      samples.push_back("VBFH210");
      samples.push_back("VBFH230");
      samples.push_back("VBFH250");
      samples.push_back("VBFH270");
      samples.push_back("VBFH300");
      samples.push_back("VBFH350");
      samples.push_back("VBFH400");
      samples.push_back("VBFH450");
      samples.push_back("VBFH500");
      samples.push_back("VBFH550");
      samples.push_back("VBFH600");
      samples.push_back("VBFH700");
      samples.push_back("VBFH750");
      samples.push_back("VBFH800");
      samples.push_back("VBFH900");
      samples.push_back("VBFH1000");
      samples.push_back("VBFH1500");
      samples.push_back("VBFH2000");
      samples.push_back("VBFH2500");
      samples.push_back("VBFH3000");
    }
  }
  else if (strsample=="WH"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_WH110.root");
      samples.push_back("HZZ4lTree_WH115.root");
      samples.push_back("HZZ4lTree_WH120.root");
      samples.push_back("HZZ4lTree_WH125.root");
      samples.push_back("HZZ4lTree_WH126.root");
      samples.push_back("HZZ4lTree_WH130.root");
      samples.push_back("HZZ4lTree_WH140.root");
      samples.push_back("HZZ4lTree_WH150.root");
      samples.push_back("HZZ4lTree_WH160.root");
      samples.push_back("HZZ4lTree_WH180.root");
      samples.push_back("HZZ4lTree_WH200.root");
    }
    else{
      samples.push_back("WminusH115");
      samples.push_back("WminusH120");
      samples.push_back("WminusH124");
      samples.push_back("WminusH125");
      samples.push_back("WminusH126");
      samples.push_back("WminusH130");
      samples.push_back("WminusH135");
      samples.push_back("WminusH140");
      samples.push_back("WminusH145");
      samples.push_back("WminusH150");
      samples.push_back("WminusH155");
      samples.push_back("WminusH160");
      samples.push_back("WminusH165");
      samples.push_back("WminusH170");
      samples.push_back("WminusH175");
      samples.push_back("WminusH180");
      samples.push_back("WminusH190");
      samples.push_back("WminusH200");
      samples.push_back("WminusH210");
      samples.push_back("WminusH230");
      samples.push_back("WminusH250");
      samples.push_back("WminusH270");
      samples.push_back("WminusH300");
      samples.push_back("WminusH350");
      samples.push_back("WminusH400");
      samples.push_back("WminusH450");
      samples.push_back("WminusH500");
      samples.push_back("WminusH550");
      samples.push_back("WminusH600");
      samples.push_back("WminusH700");
      samples.push_back("WminusH750");
      samples.push_back("WminusH800");
      samples.push_back("WminusH900");
      samples.push_back("WminusH1000");
      samples.push_back("WminusH1500");
      samples.push_back("WminusH2000");
      samples.push_back("WminusH2500");
      samples.push_back("WminusH3000");
      samples.push_back("WplusH115");
      samples.push_back("WplusH120");
      samples.push_back("WplusH124");
      samples.push_back("WplusH125");
      samples.push_back("WplusH126");
      samples.push_back("WplusH130");
      samples.push_back("WplusH135");
      samples.push_back("WplusH140");
      samples.push_back("WplusH145");
      samples.push_back("WplusH150");
      samples.push_back("WplusH155");
      samples.push_back("WplusH160");
      samples.push_back("WplusH165");
      samples.push_back("WplusH170");
      samples.push_back("WplusH175");
      samples.push_back("WplusH180");
      samples.push_back("WplusH190");
      samples.push_back("WplusH200");
      samples.push_back("WplusH210");
      samples.push_back("WplusH230");
      samples.push_back("WplusH250");
      samples.push_back("WplusH270");
      samples.push_back("WplusH300");
      samples.push_back("WplusH350");
      samples.push_back("WplusH400");
      samples.push_back("WplusH450");
      samples.push_back("WplusH500");
      samples.push_back("WplusH550");
      samples.push_back("WplusH600");
      samples.push_back("WplusH700");
      samples.push_back("WplusH750");
      samples.push_back("WplusH800");
      samples.push_back("WplusH900");
      samples.push_back("WplusH1000");
      samples.push_back("WplusH1500");
      samples.push_back("WplusH2000");
      samples.push_back("WplusH2500");
      samples.push_back("WplusH3000");
    }
  }
  else if (strsample=="ZH"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_ZH110.root");
      samples.push_back("HZZ4lTree_ZH115.root");
      samples.push_back("HZZ4lTree_ZH120.root");
      samples.push_back("HZZ4lTree_ZH125.root");
      samples.push_back("HZZ4lTree_ZH126.root");
      samples.push_back("HZZ4lTree_ZH130.root");
      samples.push_back("HZZ4lTree_ZH140.root");
      samples.push_back("HZZ4lTree_ZH150.root");
      samples.push_back("HZZ4lTree_ZH160.root");
      samples.push_back("HZZ4lTree_ZH180.root");
      samples.push_back("HZZ4lTree_ZH200.root");
    }
    else{
      samples.push_back("ZH115");
      samples.push_back("ZH120");
      samples.push_back("ZH124");
      samples.push_back("ZH125");
      samples.push_back("ZH126");
      samples.push_back("ZH130");
      samples.push_back("ZH135");
      samples.push_back("ZH140");
      samples.push_back("ZH145");
      samples.push_back("ZH150");
      samples.push_back("ZH155");
      samples.push_back("ZH160");
      samples.push_back("ZH165");
      samples.push_back("ZH170");
      samples.push_back("ZH175");
      samples.push_back("ZH180");
      samples.push_back("ZH190");
      samples.push_back("ZH200");
      samples.push_back("ZH210");
      samples.push_back("ZH230");
      samples.push_back("ZH250");
      samples.push_back("ZH270");
      samples.push_back("ZH300");
      samples.push_back("ZH350");
      samples.push_back("ZH400");
      samples.push_back("ZH450");
      samples.push_back("ZH500");
      samples.push_back("ZH550");
      samples.push_back("ZH600");
      samples.push_back("ZH700");
      samples.push_back("ZH750");
      samples.push_back("ZH800");
      samples.push_back("ZH900");
      samples.push_back("ZH1000");
      samples.push_back("ZH1500");
      samples.push_back("ZH2000");
      samples.push_back("ZH2500");
      samples.push_back("ZH3000");
    }
  }
  else if (strsample=="JJQCD"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_minloH90.root");
      samples.push_back("HZZ4lTree_minloH95.root");
      samples.push_back("HZZ4lTree_minloH100.root");
      samples.push_back("HZZ4lTree_minloH105.root");
      samples.push_back("HZZ4lTree_minloH110.root");
      samples.push_back("HZZ4lTree_minloH115.root");
      samples.push_back("HZZ4lTree_minloH120.root");
      samples.push_back("HZZ4lTree_minloH124.root");
      samples.push_back("HZZ4lTree_minloH125.root");
      samples.push_back("HZZ4lTree_minloH126.root");
      samples.push_back("HZZ4lTree_minloH130.root");
      samples.push_back("HZZ4lTree_minloH135.root");
      samples.push_back("HZZ4lTree_minloH140.root");
      samples.push_back("HZZ4lTree_minloH145.root");
      samples.push_back("HZZ4lTree_minloH150.root");
      samples.push_back("HZZ4lTree_minloH155.root");
      samples.push_back("HZZ4lTree_minloH160.root");
      samples.push_back("HZZ4lTree_minloH170.root");
      samples.push_back("HZZ4lTree_minloH180.root");
      samples.push_back("HZZ4lTree_minloH190.root");
      samples.push_back("HZZ4lTree_minloH200.root");
      samples.push_back("HZZ4lTree_minloH250.root");
      samples.push_back("HZZ4lTree_minloH300.root");
      samples.push_back("HZZ4lTree_minloH350.root");
      samples.push_back("HZZ4lTree_minloH400.root");
      samples.push_back("HZZ4lTree_minloH450.root");
      samples.push_back("HZZ4lTree_minloH500.root");
      samples.push_back("HZZ4lTree_minloH550.root");
      samples.push_back("HZZ4lTree_minloH600.root");
      samples.push_back("HZZ4lTree_minloH650.root");
      samples.push_back("HZZ4lTree_minloH700.root");
      samples.push_back("HZZ4lTree_minloH750.root");
      samples.push_back("HZZ4lTree_minloH800.root");
      samples.push_back("HZZ4lTree_minloH850.root");
      samples.push_back("HZZ4lTree_minloH900.root");
      samples.push_back("HZZ4lTree_minloH950.root");
      samples.push_back("HZZ4lTree_minloH1000.root");
    }
    else{
      samples.push_back("ggH115");
      samples.push_back("ggH120");
      samples.push_back("ggH124");
      samples.push_back("ggH125");
      samples.push_back("ggH126");
      samples.push_back("ggH130");
      samples.push_back("ggH135");
      samples.push_back("ggH140");
      samples.push_back("ggH145");
      samples.push_back("ggH150");
      samples.push_back("ggH155");
      samples.push_back("ggH160");
      samples.push_back("ggH165");
      samples.push_back("ggH170");
      samples.push_back("ggH175");
      samples.push_back("ggH180");
      samples.push_back("ggH190");
      samples.push_back("ggH200");
      samples.push_back("ggH210");
      samples.push_back("ggH230");
      samples.push_back("ggH250");
      samples.push_back("ggH270");
      samples.push_back("ggH300");
      samples.push_back("ggH350");
      samples.push_back("ggH400");
      samples.push_back("ggH450");
      samples.push_back("ggH500");
      samples.push_back("ggH550");
      samples.push_back("ggH600");
      samples.push_back("ggH700");
      samples.push_back("ggH750");
      samples.push_back("ggH800");
      samples.push_back("ggH900");
      samples.push_back("ggH1000");
      samples.push_back("ggH1500");
      samples.push_back("ggH2000");
      samples.push_back("ggH2500");
      samples.push_back("ggH3000");
    }
  }
  else if (strsample=="gg_Sig_JHUGen"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_jhuGenV4-H91.2.root");
      samples.push_back("HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root");
    }
    else{
    }
  }
  else if (strsample=="gg_Sig_MCFM"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_ggTo4mu_SMH-MCFM67_H125.6.root");
      samples.push_back("HZZ4lTree_ggTo4e_SMH-MCFM67_H125.6.root");
      samples.push_back("HZZ4lTree_ggTo2e2mu_SMH-MCFM67_H125.6.root");
    }
    else{
      samples.push_back("ggTo4mu_0PMH125_MCFM701");
      samples.push_back("ggTo4e_0PMH125_MCFM701");
      samples.push_back("ggTo2e2mu_0PMH125_MCFM701");
      samples.push_back("ggTo2e2tau_0PMH125_MCFM701");
      samples.push_back("ggTo2mu2tau_0PMH125_MCFM701");
      samples.push_back("ggTo4tau_0PMH125_MCFM701");
    }
  }
  else if (strsample=="gg_Bkg_MCFM"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_ggTo2e2mu_Contin-MCFM67.root");
      samples.push_back("HZZ4lTree_ggTo4mu_Contin-MCFM67.root");
      samples.push_back("HZZ4lTree_ggTo4e_Contin-MCFM67.root");
    }
    else{
      samples.push_back("ggTo4mu_Contin_MCFM701");
      samples.push_back("ggTo4e_Contin_MCFM701");
      samples.push_back("ggTo2e2mu_Contin_MCFM701");
      samples.push_back("ggTo2e2tau_Contin_MCFM701");
      samples.push_back("ggTo2mu2tau_Contin_MCFM701");
      samples.push_back("ggTo4tau_Contin_MCFM701");
    }
  }
  else if (strsample=="gg_Sig_ggVV"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_ggTo2l2l_H125.6.root");
      samples.push_back("HZZ4lTree_ggTo4l_H125.6.root");
    }
    else{
    }
  }
  else if (strsample=="gg_Bkg_ggVV"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_ggTo4l_Continuum.root");
      samples.push_back("HZZ4lTree_ggZZ4l.root");
      samples.push_back("HZZ4lTree_ggTo2l2l_Continuum.root");
      samples.push_back("HZZ4lTree_ggZZ2l2l.root");
    }
    else{
    }
  }
  else if (strsample=="VV_Sig_MCFM"){
    if (sqrts<10.){
    }
    else{
      samples.push_back("VBFTo2e2muJJ_0PMH125_phantom128");
      samples.push_back("VBFTo4muJJ_0PMH125_phantom128");
      samples.push_back("VBFTo4eJJ_0PMH125_phantom128");
    }
  }
  else if (strsample=="VV_Bkg_MCFM"){
    if (sqrts<10.){
    }
    else{
      samples.push_back("VBFTo2e2muJJ_Contin_phantom128");
      samples.push_back("VBFTo4muJJ_Contin_phantom128");
      samples.push_back("VBFTo4eJJ_Contin_phantom128");
    }
  }
  else if (strsample=="qq_Bkg"){
    if (sqrts<10.){
      samples.push_back("HZZ4lTree_ZZTo2e2mu.root");
      samples.push_back("HZZ4lTree_ZZTo2e2tau.root");
      samples.push_back("HZZ4lTree_ZZTo2mu2tau.root");
      samples.push_back("HZZ4lTree_ZZTo4mu.root");
      samples.push_back("HZZ4lTree_ZZTo4e.root");
      samples.push_back("HZZ4lTree_ZZTo4tau.root");
      samples.push_back("HZZ4lTree_ZZ95-160To2e2mu.root");
      samples.push_back("HZZ4lTree_ZZ95-160To2e2tau.root");
      samples.push_back("HZZ4lTree_ZZ95-160To2mu2tau.root");
      samples.push_back("HZZ4lTree_ZZ95-160To4mu.root");
      samples.push_back("HZZ4lTree_ZZ95-160To4e.root");
      samples.push_back("HZZ4lTree_ZZ95-160To4tau.root");
    }
    else{
      samples.push_back("ZZTo4l");
      samples.push_back("ZZTo4l_ext");
    }
  }
  return samples;
}

/*
GENERAL COMMENTS
- ALL MES ARE TRANSFORMED TO LOG10(ME).
*/

/* SPECIFIC COMMENT: NONE */
void get_PAvgProfile_JHUGen_JJVBF_HSMHiggs_7or8TeV(int sqrts=8){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "SelectedTree";
  bool writeFinalTree=true;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::JJVBF;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

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
  float jetptetaphimass[2][4];

  float mesq_conserveDifermMass=0;
  float mesq_jetPtoEScale=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  int LepID[4]={ 13, -13, 11, -11 };
  short Z1Flav, Z2Flav;

  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TString cinput_main;
  if (sqrts==8) cinput_main = inputdir_8TeV;
  else if (sqrts==7) cinput_main = inputdir_7TeV;
  else return;

  vector<TString> strSamples = constructSamplesList("JJVBF", sqrts);

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;

  for (int is=0; is<(int)strSamples.size(); is++){
    for (int ic=0; ic<3; ic++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
            tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
            tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("NJets30", 1); tree->SetBranchAddress("NJets30", &NJets30);
            tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
            tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
            tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
            tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);

            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }

  const int nSamples = treeList.size();

  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;

  vector<pair<float, int>> index;
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    if (NJets30<2) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    addByLowest(index, mzz, ev);
    ev_acc++;
  }

  float firstVal=index.at(0).first;
  float lastVal=index.at(index.size()-1).first;
  float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
  float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
  cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

  float divisor=85000;
  int nbins = index.size()/divisor;
  const int nbins_th=10/*50*/;
  while (nbins<nbins_th){
    if (divisor>1000) divisor -= 1000;
    else if (divisor>100) divisor -= 100;
    else break;
    nbins=index.size()/divisor;
  }
  cout << "nbins=" << nbins << endl;
  if (nbins<3) cerr << "Not enough bins!" << endl;
  vector<ExtBin> binList;
  float* binning = new float[nbins+1];
  binning[0]=infimum;
  binning[nbins]=supremum;
  int ev_stepsize = index.size()/nbins;
  cout << "Event step size: " << ev_stepsize << endl;
  cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
  for (int ix=1; ix<nbins; ix++){
    binning[ix]=(index[ix*ev_stepsize-1].first+index[ix*ev_stepsize].first)*0.5;
    ExtBin tmpbin;
    tmpbin.binlow = binning[ix-1];
    tmpbin.binhigh = binning[ix];
    for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[(ix-1)*ev_stepsize+bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
  }
  ExtBin tmpbin;
  tmpbin.binlow = binning[nbins-1];
  tmpbin.binhigh = binning[nbins];
  for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index.size(); bin++) tmpbin.events.push_back(index[bin].second);
  binList.push_back(tmpbin);
  cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
  cout << "Bin list has the following bins:" << endl;
  for (unsigned int ib=0; ib<binList.size(); ib++){
    cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
  }

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  TProfile* hvar = new TProfile("candMass", "", nbins, binning); hvar->Sumw2();
  TProfile* hmesq_conserveDifermMass = new TProfile("P_ConserveDifermionMass", "", nbins, binning); hmesq_conserveDifermMass->Sumw2();
  TProfile* hmesq_jetPtoEScale = new TProfile("P_MomentumToEnergy", "", nbins, binning); hmesq_jetPtoEScale->Sumw2();
  
  TTree* newtree=0;
  if (writeFinalTree){
    newtree = new TTree("FinalTree", "");
    newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
    newtree->Branch("mesq_jetPtoEScale", &mesq_jetPtoEScale);
    newtree->Branch("ZZMass", &mzz);
  }

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (unsigned int bin=0; bin<binList.size(); bin++){
    cout << "Bin " << bin << " is now being scrutinized..." << endl;
    for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
      int getEv = binList.at(bin).events.at(ev);
      getEntry(treeList, getEv);
      if (ev%1000==0) cout << "Doing event " << getEv << endl;
      if (JetPt->size()<2 || JetEta->size()<2 || JetPhi->size()<2 || JetMass->size()<2){
        cerr << "Jet array sizes are less than Njets!" << endl;
        continue;
      }

      TLorentzVector jet[2], higgs;
      for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(jetptetaphimass[ij][0], jetptetaphimass[ij][1], jetptetaphimass[ij][2], jetptetaphimass[ij][3]);
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));
      associated.push_back(SimpleParticle_t(0, jet[1]));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      mela.setProcess(proc, me, prod);

      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_conserveDifermMass, false);
      TUtil::setJetMassScheme(TVar::MomentumToEnergy);
      mela.computeProdP(mesq_jetPtoEScale, false);

      bool doFill = !(
        isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
        ||
        isnan(mesq_jetPtoEScale) || isinf(mesq_jetPtoEScale)
        );

      if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, mesq_jetPtoEScale);

      mela.resetInputEvent();
    }

    binList.at(bin).sift(); binList.at(bin).adjustWeights();

    for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
      mzz = binList.at(bin).masses.at(ev);
      mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
      mesq_jetPtoEScale = binList.at(bin).me2vals.at(ev);
      hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
      hmesq_jetPtoEScale->Fill(mzz, mesq_jetPtoEScale);
      hvar->Fill(mzz, mzz);
      if (writeFinalTree) newtree->Fill();
    }
  }

  double* xexyey[2][4];
  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) xexyey[inorm][ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[inorm][0][bin] = hvar->GetBinContent(bin+1);
      xexyey[inorm][1][bin] = hvar->GetBinError(bin+1);

      if (inorm==0) cout << "Bin " << bin << " x-center: " << xexyey[inorm][0][bin] << " +- " << xexyey[inorm][1][bin] << endl;

      if (inorm==0){
        xexyey[inorm][2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      }
      else{
        xexyey[inorm][2][bin] = hmesq_jetPtoEScale->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_jetPtoEScale->GetBinError(bin+1);
      }
      xexyey[inorm][3][bin] = log10(xexyey[inorm][3][bin])/xexyey[inorm][2][bin];
      xexyey[inorm][2][bin] = log10(xexyey[inorm][2][bin]);
    }
  }

  for (int inorm=0; inorm<2; inorm++){
    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[inorm][0], xexyey[inorm][2], xexyey[inorm][1], xexyey[inorm][3]);
    if (inorm==0) tg->SetName("tg_P_ConserveDifermionMass");
    else tg->SetName("tg_P_MomentumToEnergy");
    foutput->WriteTObject(tg);
    delete tg;
  }

  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) delete[] xexyey[inorm][ix];
  }
  foutput->WriteTObject(hmesq_jetPtoEScale);
  foutput->WriteTObject(hmesq_conserveDifermMass);
  foutput->WriteTObject(hvar);
  if (writeFinalTree) foutput->WriteTObject(newtree);
  if (writeFinalTree) delete newtree;
  delete hmesq_conserveDifermMass;
  delete hmesq_jetPtoEScale;
  delete hvar;
  foutput->Close();
  delete[] binning;
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}
/*
Function
[0]*exp(-x/[1])*(1+[2]*exp(-pow(x/[3],2)))
with parameters
0.0187935
489.335
0.0870576
256.215
fits well.
*/


/*
SPECIFIC COMMENT: OUTPUT ME DIVIDED BY 
- ALPHAS(MZ)**4 TO REMAIN INDEPENDENT OF PDF CHOICE TO FIRST APPROXIMATION
*/
void get_PAvgProfile_JHUGen_JJQCD_HSMHiggs_7or8TeV(int sqrts=8){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "SelectedTree";
  bool writeFinalTree=true;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::JJQCD;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

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
  float jetptetaphimass[2][4];

  float mesq_conserveDifermMass=0;
  float mesq_jetPtoEScale=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  int LepID[4]={ 13, -13, 11, -11 };
  short Z1Flav, Z2Flav;

  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TString cinput_main;
  if (sqrts==8) cinput_main = inputdir_8TeV;
  else if (sqrts==7) cinput_main = inputdir_7TeV;
  else return;

  vector<TString> strSamples = constructSamplesList("JJQCD", sqrts);

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;

  for (int is=0; is<(int)strSamples.size(); is++){
    for (int ic=0; ic<3; ic++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
            tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
            tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("NJets30", 1); tree->SetBranchAddress("NJets30", &NJets30);
            tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
            tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
            tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
            tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);

            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }

  const int nSamples = treeList.size();

  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;

  vector<pair<float, int>> index;
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    if (NJets30<2) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    addByLowest(index, mzz, ev);
    ev_acc++;
  }

  float firstVal=index.at(0).first;
  float lastVal=index.at(index.size()-1).first;
  float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
  float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
  cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

  float divisor=85000;
  int nbins = index.size()/divisor;
  const int nbins_th=10/*50*/;
  while (nbins<nbins_th){
    if (divisor>1000) divisor -= 1000;
    else if (divisor>100) divisor -= 100;
    else break;
    nbins=index.size()/divisor;
  }
  cout << "nbins=" << nbins << endl;
  if (nbins<3) cerr << "Not enough bins!" << endl;
  vector<ExtBin> binList;
  float* binning = new float[nbins+1];
  binning[0]=infimum;
  binning[nbins]=supremum;
  int ev_stepsize = index.size()/nbins;
  cout << "Event step size: " << ev_stepsize << endl;
  cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
  for (int ix=1; ix<nbins; ix++){
    binning[ix]=(index[ix*ev_stepsize-1].first+index[ix*ev_stepsize].first)*0.5;
    ExtBin tmpbin;
    tmpbin.binlow = binning[ix-1];
    tmpbin.binhigh = binning[ix];
    for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[(ix-1)*ev_stepsize+bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
  }
  ExtBin tmpbin;
  tmpbin.binlow = binning[nbins-1];
  tmpbin.binhigh = binning[nbins];
  for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index.size(); bin++) tmpbin.events.push_back(index[bin].second);
  binList.push_back(tmpbin);
  cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
  cout << "Bin list has the following bins:" << endl;
  for (unsigned int ib=0; ib<binList.size(); ib++){
    cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
  }

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  TProfile* hvar = new TProfile("candMass", "", nbins, binning); hvar->Sumw2();
  TProfile* hmesq_conserveDifermMass = new TProfile("P_ConserveDifermionMass", "", nbins, binning); hmesq_conserveDifermMass->Sumw2();
  TProfile* hmesq_jetPtoEScale = new TProfile("P_MomentumToEnergy", "", nbins, binning); hmesq_jetPtoEScale->Sumw2();

  TTree* newtree=0;
  if (writeFinalTree){
    newtree = new TTree("FinalTree", "");
    newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
    newtree->Branch("mesq_jetPtoEScale", &mesq_jetPtoEScale);
    newtree->Branch("ZZMass", &mzz);
  }

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (unsigned int bin=0; bin<binList.size(); bin++){
    cout << "Bin " << bin << " is now being scrutinized..." << endl;
    for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
      int getEv = binList.at(bin).events.at(ev);
      getEntry(treeList, getEv);
      if (ev%1000==0) cout << "Doing event " << getEv << endl;
      if (JetPt->size()<2 || JetEta->size()<2 || JetPhi->size()<2 || JetMass->size()<2){
        cerr << "Jet array sizes are less than Njets!" << endl;
        continue;
      }

      TLorentzVector jet[2], higgs;
      for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(jetptetaphimass[ij][0], jetptetaphimass[ij][1], jetptetaphimass[ij][2], jetptetaphimass[ij][3]);
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));
      associated.push_back(SimpleParticle_t(0, jet[1]));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      mela.setProcess(proc, me, prod);

      double alphasVal;

      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_conserveDifermMass, false);
      alphasVal = mela.getIORecord()->getAlphaSatMZ();
      mesq_conserveDifermMass /= pow(alphasVal, 4);

      TUtil::setJetMassScheme(TVar::MomentumToEnergy);
      mela.computeProdP(mesq_jetPtoEScale, false);
      alphasVal = mela.getIORecord()->getAlphaSatMZ();
      mesq_jetPtoEScale /= pow(alphasVal, 4);

      bool doFill = !(
        isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
        ||
        isnan(mesq_jetPtoEScale) || isinf(mesq_jetPtoEScale)
        );

      if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, mesq_jetPtoEScale);

      mela.resetInputEvent();
    }

    binList.at(bin).sift(); binList.at(bin).adjustWeights();

    for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
      mzz = binList.at(bin).masses.at(ev);
      mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
      mesq_jetPtoEScale = binList.at(bin).me2vals.at(ev);
      hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
      hmesq_jetPtoEScale->Fill(mzz, mesq_jetPtoEScale);
      hvar->Fill(mzz, mzz);
      if (writeFinalTree) newtree->Fill();
    }
  }

  double* xexyey[2][4];
  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) xexyey[inorm][ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[inorm][0][bin] = hvar->GetBinContent(bin+1);
      xexyey[inorm][1][bin] = hvar->GetBinError(bin+1);

      if (inorm==0) cout << "Bin " << bin << " x-center: " << xexyey[inorm][0][bin] << " +- " << xexyey[inorm][1][bin] << endl;

      if (inorm==0){
        xexyey[inorm][2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      }
      else{
        xexyey[inorm][2][bin] = hmesq_jetPtoEScale->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_jetPtoEScale->GetBinError(bin+1);
      }
      xexyey[inorm][3][bin] = log10(xexyey[inorm][3][bin])/xexyey[inorm][2][bin];
      xexyey[inorm][2][bin] = log10(xexyey[inorm][2][bin]);
    }
  }

  for (int inorm=0; inorm<2; inorm++){
    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[inorm][0], xexyey[inorm][2], xexyey[inorm][1], xexyey[inorm][3]);
    if (inorm==0) tg->SetName("tg_P_ConserveDifermionMass");
    else tg->SetName("tg_P_MomentumToEnergy");
    foutput->WriteTObject(tg);
    delete tg;
  }

  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) delete[] xexyey[inorm][ix];
  }
  foutput->WriteTObject(hmesq_jetPtoEScale);
  foutput->WriteTObject(hmesq_conserveDifermMass);
  foutput->WriteTObject(hvar);
  if (writeFinalTree) foutput->WriteTObject(newtree);
  if (writeFinalTree) delete newtree;
  delete hmesq_conserveDifermMass;
  delete hmesq_jetPtoEScale;
  delete hvar;
  foutput->Close();
  delete[] binning;
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}
/*
Function
[0]*exp(-x/[1])*(1+[2]*exp(-pow(x/[3],2)) + [4]*exp(-x/[5]))
with parameters
152.197
496.507
0.812036
349.521
9.00697
55.4923
fits well.
*/

/* SPECIFIC COMMENT: OUTPUT ME DIVIDED BY ALPHAS(MZ)**3 TO REMAIN INDEPENDENT OF PDF CHOICE TO FIRST APPROXIMATION */
void get_PAvgProfile_JHUGen_JQCD_HSMHiggs_7or8TeV(int sqrts=8){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "SelectedTree";
  bool writeFinalTree=true;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::JQCD;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

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
  float jetptetaphimass[2][4];

  float mesq_conserveDifermMass=0;
  float mesq_jetPtoEScale=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  int LepID[4]={ 13, -13, 11, -11 };
  short Z1Flav, Z2Flav;

  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TString cinput_main;
  if (sqrts==8) cinput_main = inputdir_8TeV;
  else if (sqrts==7) cinput_main = inputdir_7TeV;
  else return;

  vector<TString> strSamples = constructSamplesList("JJQCD", sqrts);

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;

  for (int is=0; is<(int)strSamples.size(); is++){
    for (int ic=0; ic<3; ic++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
            tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
            tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("NJets30", 1); tree->SetBranchAddress("NJets30", &NJets30);
            tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
            tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
            tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
            tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);

            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }

  const int nSamples = treeList.size();

  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;

  vector<pair<float, int>> index;
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    if (NJets30!=1) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    addByLowest(index, mzz, ev);
    ev_acc++;
  }

  float firstVal=index.at(0).first;
  float lastVal=index.at(index.size()-1).first;
  float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
  float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
  cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

  float divisor=85000;
  int nbins = index.size()/divisor;
  const int nbins_th=10/*50*/;
  while (nbins<nbins_th){
    if (divisor>1000) divisor -= 1000;
    else if (divisor>100) divisor -= 100;
    else break;
    nbins=index.size()/divisor;
  }
  cout << "nbins=" << nbins << endl;
  if (nbins<3) cerr << "Not enough bins!" << endl;
  vector<ExtBin> binList;
  float* binning = new float[nbins+1];
  binning[0]=infimum;
  binning[nbins]=supremum;
  int ev_stepsize = index.size()/nbins;
  cout << "Event step size: " << ev_stepsize << endl;
  cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
  for (int ix=1; ix<nbins; ix++){
    binning[ix]=(index[ix*ev_stepsize-1].first+index[ix*ev_stepsize].first)*0.5;
    ExtBin tmpbin;
    tmpbin.binlow = binning[ix-1];
    tmpbin.binhigh = binning[ix];
    for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[(ix-1)*ev_stepsize+bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
  }
  ExtBin tmpbin;
  tmpbin.binlow = binning[nbins-1];
  tmpbin.binhigh = binning[nbins];
  for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index.size(); bin++) tmpbin.events.push_back(index[bin].second);
  binList.push_back(tmpbin);
  cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
  cout << "Bin list has the following bins:" << endl;
  for (unsigned int ib=0; ib<binList.size(); ib++){
    cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
  }

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  TProfile* hvar = new TProfile("candMass", "", nbins, binning); hvar->Sumw2();
  TProfile* hmesq_conserveDifermMass = new TProfile("P_ConserveDifermionMass", "", nbins, binning); hmesq_conserveDifermMass->Sumw2();
  TProfile* hmesq_jetPtoEScale = new TProfile("P_MomentumToEnergy", "", nbins, binning); hmesq_jetPtoEScale->Sumw2();

  TTree* newtree=0;
  if (writeFinalTree){
    newtree = new TTree("FinalTree", "");
    newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
    newtree->Branch("mesq_jetPtoEScale", &mesq_jetPtoEScale);
    newtree->Branch("ZZMass", &mzz);
  }

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (unsigned int bin=0; bin<binList.size(); bin++){
    cout << "Bin " << bin << " is now being scrutinized..." << endl;
    for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
      int getEv = binList.at(bin).events.at(ev);
      getEntry(treeList, getEv);
      if (ev%1000==0) cout << "Doing event " << getEv << endl;
      if (JetPt->size()<1 || JetEta->size()<1 || JetPhi->size()<1 || JetMass->size()<1){
        cerr << "Jet array sizes are less than Njets!" << endl;
        continue;
      }

      TLorentzVector jet, higgs;
      jet.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      mela.setProcess(proc, me, prod);

      double alphasVal;

      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_conserveDifermMass, false);
      alphasVal = mela.getIORecord()->getAlphaSatMZ();
      mesq_conserveDifermMass /= pow(alphasVal, 3);

      TUtil::setJetMassScheme(TVar::MomentumToEnergy);
      mela.computeProdP(mesq_jetPtoEScale, false);
      alphasVal = mela.getIORecord()->getAlphaSatMZ();
      mesq_jetPtoEScale /= pow(alphasVal, 3);

      bool doFill = !(
        isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
        ||
        isnan(mesq_jetPtoEScale) || isinf(mesq_jetPtoEScale)
        );

      if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, mesq_jetPtoEScale);

      mela.resetInputEvent();
    }

    binList.at(bin).sift(); binList.at(bin).adjustWeights();

    for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
      mzz = binList.at(bin).masses.at(ev);
      mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
      mesq_jetPtoEScale = binList.at(bin).me2vals.at(ev);
      hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
      hmesq_jetPtoEScale->Fill(mzz, mesq_jetPtoEScale);
      hvar->Fill(mzz, mzz);
      if (writeFinalTree) newtree->Fill();
    }
  }

  double* xexyey[2][4];
  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) xexyey[inorm][ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[inorm][0][bin] = hvar->GetBinContent(bin+1);
      xexyey[inorm][1][bin] = hvar->GetBinError(bin+1);

      if (inorm==0) cout << "Bin " << bin << " x-center: " << xexyey[inorm][0][bin] << " +- " << xexyey[inorm][1][bin] << endl;

      if (inorm==0){
        xexyey[inorm][2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      }
      else{
        xexyey[inorm][2][bin] = hmesq_jetPtoEScale->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_jetPtoEScale->GetBinError(bin+1);
      }
      xexyey[inorm][3][bin] = log10(xexyey[inorm][3][bin])/xexyey[inorm][2][bin];
      xexyey[inorm][2][bin] = log10(xexyey[inorm][2][bin]);
    }
  }

  for (int inorm=0; inorm<2; inorm++){
    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[inorm][0], xexyey[inorm][2], xexyey[inorm][1], xexyey[inorm][3]);
    if (inorm==0) tg->SetName("tg_P_ConserveDifermionMass");
    else tg->SetName("tg_P_MomentumToEnergy");
    foutput->WriteTObject(tg);
    delete tg;
  }

  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) delete[] xexyey[inorm][ix];
  }
  foutput->WriteTObject(hmesq_jetPtoEScale);
  foutput->WriteTObject(hmesq_conserveDifermMass);
  foutput->WriteTObject(hvar);
  if (writeFinalTree) foutput->WriteTObject(newtree);
  if (writeFinalTree) delete newtree;
  delete hmesq_conserveDifermMass;
  delete hmesq_jetPtoEScale;
  delete hvar;
  foutput->Close();
  delete[] binning;
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}
/*
Function
[0]*exp(-x/[1])*(1+[2]*exp(-pow(x/[3],2)) + [4]*exp(-x/[5]))
with parameters
196.358
291.176
14.6094
92.3443
13.2622
133.669
fits well.
*/


/***** 13 TeV Samples *****/
/*
SPECIFIC COMMENT:
- NJets30 -> nCleanedJetsPt30
- vector<double> -> vector<float>
*/

/* SPECIFIC COMMENT: NONE */
void get_PAvgProfile_JHUGen_JJVBF_HSMHiggs_13TeV(int sqrts=13, bool recalculate = true){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";
  bool writeFinalTree=true;

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::JJVBF;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  short NJets30;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float> myJetPt;
  std::vector<float> myJetEta;
  std::vector<float> myJetPhi;
  std::vector<float> myJetMass;

  float mesq_calc=0., cconst_calc=1.;
  float mesq_conserveDifermMass=0;
  float mesq_jetPtoEScale=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  float genHEPMCweight;
  float wgt=1;

  TString cinput_main;
  if (sqrts==13) cinput_main = inputdir_13TeV;
  else return;

  vector<TString> strSamples = constructSamplesList("JJVBF", sqrts);

  unordered_map<TTree*, pair<float, float>> nGenMap;
  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  for (unsigned int is=0; is<strSamples.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          if (!recalculate && tree->GetBranchStatus("pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal")==0) recalculate=true;
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("genHEPMCweight", 1); tree->SetBranchAddress("genHEPMCweight", &genHEPMCweight);
          tree->SetBranchStatus("nCleanedJetsPt30", 1); tree->SetBranchAddress("nCleanedJetsPt30", &NJets30);
          tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
          tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
          tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
          tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
          tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
          tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("LepLepId", 1); tree->SetBranchAddress("LepLepId", &LepLepId);
          tree->SetBranchStatus("LepPt", 1); tree->SetBranchAddress("LepPt", &LepPt);
          tree->SetBranchStatus("LepEta", 1); tree->SetBranchAddress("LepEta", &LepEta);
          tree->SetBranchStatus("LepPhi", 1); tree->SetBranchAddress("LepPhi", &LepPhi);
          if (!recalculate){
            tree->SetBranchStatus("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", 1);
            tree->SetBranchAddress("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &mesq_calc);
            tree->SetBranchStatus("pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", 1);
            tree->SetBranchAddress("pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal", &cconst_calc);
          }

          TH1F* htmp = (TH1F*)finput->Get(COUNTERS_NAME);
          pair<float, float> nsum(htmp->GetBinContent(1), htmp->GetBinContent(41)); // No PU reweighting
          nGenMap[tree]=nsum;

          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }
  cout << "NEntries = " << nEntries << endl;

  vector<pair<float, int>> index;
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    if (NJets30<2) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    addByLowest(index, mzz, ev);
    ev_acc++;
  }

  float firstVal=index.at(0).first;
  float lastVal=index.at(index.size()-1).first;
  float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
  float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
  cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

  float divisor=85000;
  int nbins = index.size()/divisor;
  const int nbins_th=10/*50*/;
  while (nbins<nbins_th){
    if (divisor>1000) divisor -= 1000;
    else if (divisor>100) divisor -= 100;
    else break;
    nbins=index.size()/divisor;
  }
  cout << "nbins=" << nbins << endl;
  if (nbins<3) cerr << "Not enough bins!" << endl;
  vector<ExtBin> binList;
  float* binning = new float[nbins+1];
  binning[0]=infimum;
  binning[nbins]=supremum;
  int ev_stepsize = index.size()/nbins;
  cout << "Event step size: " << ev_stepsize << endl;
  cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
  for (int ix=1; ix<nbins; ix++){
    binning[ix]=(index[ix*ev_stepsize-1].first+index[ix*ev_stepsize].first)*0.5;
    ExtBin tmpbin;
    tmpbin.binlow = binning[ix-1];
    tmpbin.binhigh = binning[ix];
    for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[(ix-1)*ev_stepsize+bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
  }
  ExtBin tmpbin;
  tmpbin.binlow = binning[nbins-1];
  tmpbin.binhigh = binning[nbins];
  for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index.size(); bin++) tmpbin.events.push_back(index[bin].second);
  binList.push_back(tmpbin);
  cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
  cout << "Bin list has the following bins:" << endl;
  for (unsigned int ib=0; ib<binList.size(); ib++){
    cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    /*
    cout << "\tEvents: ";
    for (unsigned int ev=0; ev<binList.at(ib).events.size(); ev++){
      cout << binList.at(ib).events.at(ev) << " ";
    }
    cout << endl;
    */
  }
  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  TProfile* hvar = new TProfile("candMass", "", nbins, binning); hvar->Sumw2();
  TProfile* hmesq_conserveDifermMass = new TProfile("P_ConserveDifermionMass", "", nbins, binning); hmesq_conserveDifermMass->Sumw2();
  TProfile* hmesq_jetPtoEScale = new TProfile("P_MomentumToEnergy", "", nbins, binning); hmesq_jetPtoEScale->Sumw2();


  TTree* newtree=0;
  if (writeFinalTree){
    newtree = new TTree("FinalTree", "");
    newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
    newtree->Branch("mesq_jetPtoEScale", &mesq_jetPtoEScale);
    newtree->Branch("ZZMass", &mzz);
  }

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (unsigned int bin=0; bin<binList.size(); bin++){
    cout << "Bin " << bin << " is now being scrutinized..." << endl;
    for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
      int getEv = binList.at(bin).events.at(ev);
      getEntry(treeList, getEv);
      if (ev%1000==0) cout << "Doing event " << getEv << endl;
      if (JetPt->size()<2 || JetEta->size()<2 || JetPhi->size()<2 || JetMass->size()<2){
        cerr << "Jet array sizes are less than Njets!" << endl;
        continue;
      }
      TTree* tree = findTree(treeList, getEv);
      wgt = fabs(genHEPMCweight*nGenMap[tree].first/nGenMap[tree].second);

      TLorentzVector jet[2], higgs;
      for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(JetPt->at(ij), JetEta->at(ij), JetPhi->at(ij), JetMass->at(ij));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));
      associated.push_back(SimpleParticle_t(0, jet[1]));

      SimpleParticleCollection_t daughters;
      for (int id=0; id<4; id++){
        double mass=0;
        if (abs(LepLepId->at(id))==13) mass = 0.105658;
        else if (abs(LepLepId->at(id))==11) mass = 0.000511;
        TLorentzVector pDaughter;
        pDaughter.SetPtEtaPhiM(LepPt->at(id), LepEta->at(id), LepPhi->at(id), mass);
        daughters.push_back(SimpleParticle_t(LepLepId->at(id), pDaughter));
      }
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      mela.setProcess(proc, me, prod);

      if (recalculate){
        TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
        mela.computeProdP(mesq_conserveDifermMass, false);
      }
      else{
        mesq_calc /= cconst_calc;
        mesq_conserveDifermMass=mesq_calc;
      }
      TUtil::setJetMassScheme(TVar::MomentumToEnergy);
      mela.computeProdP(mesq_jetPtoEScale, false);

      bool doFill = !(
        isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
        ||
        isnan(mesq_jetPtoEScale) || isinf(mesq_jetPtoEScale)
        );

      if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, mesq_jetPtoEScale, wgt);

      mela.resetInputEvent();
    }

    binList.at(bin).sift(); binList.at(bin).adjustWeights();

    for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
      mzz = binList.at(bin).masses.at(ev);
      mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
      mesq_jetPtoEScale = binList.at(bin).me2vals.at(ev);
      wgt = binList.at(bin).weights.at(ev);
      hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass, wgt);
      hmesq_jetPtoEScale->Fill(mzz, mesq_jetPtoEScale, wgt);
      hvar->Fill(mzz, mzz, wgt);
      if (writeFinalTree) newtree->Fill();
    }
  }

  double* xexyey[2][4];
  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) xexyey[inorm][ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[inorm][0][bin] = hvar->GetBinContent(bin+1);
      xexyey[inorm][1][bin] = hvar->GetBinError(bin+1);

      if (inorm==0) cout << "Bin " << bin << " x-center: " << xexyey[inorm][0][bin] << " +- " << xexyey[inorm][1][bin] << endl;

      if (inorm==0){
        xexyey[inorm][2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      }
      else{
        xexyey[inorm][2][bin] = hmesq_jetPtoEScale->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_jetPtoEScale->GetBinError(bin+1);
      }
      xexyey[inorm][3][bin] = log10(xexyey[inorm][3][bin])/xexyey[inorm][2][bin];
      xexyey[inorm][2][bin] = log10(xexyey[inorm][2][bin]);
    }
  }

  for (int inorm=0; inorm<2; inorm++){
    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[inorm][0], xexyey[inorm][2], xexyey[inorm][1], xexyey[inorm][3]);
    if (inorm==0) tg->SetName("tg_P_ConserveDifermionMass");
    else tg->SetName("tg_P_MomentumToEnergy");
    foutput->WriteTObject(tg);
    delete tg;
  }

  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) delete[] xexyey[inorm][ix];
  }
  foutput->WriteTObject(hmesq_jetPtoEScale);
  foutput->WriteTObject(hmesq_conserveDifermMass);
  foutput->WriteTObject(hvar);
  if (writeFinalTree) foutput->WriteTObject(newtree);
  if (writeFinalTree) delete newtree;
  delete hmesq_conserveDifermMass;
  delete hmesq_jetPtoEScale;
  delete hvar;
  foutput->Close();
  delete[] binning;
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}


/* SPECIFIC COMMENT: OUTPUT ME DIVIDED BY ALPHAS(MZ)**4 TO REMAIN INDEPENDENT OF PDF CHOICE TO FIRST APPROXIMATION */
void get_PAvgProfile_JHUGen_JJQCD_HSMHiggs_13TeV(int sqrts=13, bool recalculate = true){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";
  bool writeFinalTree=true;

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::JJQCD;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  short NJets30;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float> myJetPt;
  std::vector<float> myJetEta;
  std::vector<float> myJetPhi;
  std::vector<float> myJetMass;

  float mesq_calc=0., cconst_calc=1.;
  float mesq_conserveDifermMass=0;
  float mesq_jetPtoEScale=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  float genHEPMCweight;
  float wgt=1;

  TString cinput_main;
  if (sqrts==13) cinput_main = inputdir_13TeV;
  else return;

  vector<TString> strSamples = constructSamplesList("JJQCD", sqrts);

  unordered_map<TTree*, pair<float, float>> nGenMap;
  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  for (unsigned int is=0; is<strSamples.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          if (!recalculate && tree->GetBranchStatus("pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal")==0) recalculate=true;
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("genHEPMCweight", 1); tree->SetBranchAddress("genHEPMCweight", &genHEPMCweight);
          tree->SetBranchStatus("nCleanedJetsPt30", 1); tree->SetBranchAddress("nCleanedJetsPt30", &NJets30);
          tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
          tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
          tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
          tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
          tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
          tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("LepLepId", 1); tree->SetBranchAddress("LepLepId", &LepLepId);
          tree->SetBranchStatus("LepPt", 1); tree->SetBranchAddress("LepPt", &LepPt);
          tree->SetBranchStatus("LepEta", 1); tree->SetBranchAddress("LepEta", &LepEta);
          tree->SetBranchStatus("LepPhi", 1); tree->SetBranchAddress("LepPhi", &LepPhi);
          if (!recalculate){
            tree->SetBranchStatus("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", 1);
            tree->SetBranchAddress("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &mesq_calc);
            tree->SetBranchStatus("pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", 1);
            tree->SetBranchAddress("pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal", &cconst_calc);
          }

          TH1F* htmp = (TH1F*)finput->Get(COUNTERS_NAME);
          pair<float, float> nsum(htmp->GetBinContent(1), htmp->GetBinContent(41)); // No PU reweighting
          nGenMap[tree]=nsum;

          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }
  cout << "NEntries = " << nEntries << endl;

  vector<pair<float, int>> index;
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    if (NJets30<2) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    addByLowest(index, mzz, ev);
    ev_acc++;
  }

  float firstVal=index.at(0).first;
  float lastVal=index.at(index.size()-1).first;
  float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
  float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
  cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

  float divisor=85000;
  int nbins = index.size()/divisor;
  const int nbins_th=10/*50*/;
  while (nbins<nbins_th){
    if (divisor>1000) divisor -= 1000;
    else if (divisor>100) divisor -= 100;
    else break;
    nbins=index.size()/divisor;
  }
  cout << "nbins=" << nbins << endl;
  if (nbins<3) cerr << "Not enough bins!" << endl;
  vector<ExtBin> binList;
  float* binning = new float[nbins+1];
  binning[0]=infimum;
  binning[nbins]=supremum;
  int ev_stepsize = index.size()/nbins;
  cout << "Event step size: " << ev_stepsize << endl;
  cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
  for (int ix=1; ix<nbins; ix++){
    binning[ix]=(index[ix*ev_stepsize-1].first+index[ix*ev_stepsize].first)*0.5;
    ExtBin tmpbin;
    tmpbin.binlow = binning[ix-1];
    tmpbin.binhigh = binning[ix];
    for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[(ix-1)*ev_stepsize+bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
  }
  ExtBin tmpbin;
  tmpbin.binlow = binning[nbins-1];
  tmpbin.binhigh = binning[nbins];
  for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index.size(); bin++) tmpbin.events.push_back(index[bin].second);
  binList.push_back(tmpbin);
  cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
  cout << "Bin list has the following bins:" << endl;
  for (unsigned int ib=0; ib<binList.size(); ib++){
    cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    /*
    cout << "\tEvents: ";
    for (unsigned int ev=0; ev<binList.at(ib).events.size(); ev++){
    cout << binList.at(ib).events.at(ev) << " ";
    }
    cout << endl;
    */
  }
  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  TProfile* hvar = new TProfile("candMass", "", nbins, binning); hvar->Sumw2();
  TProfile* hmesq_conserveDifermMass = new TProfile("P_ConserveDifermionMass", "", nbins, binning); hmesq_conserveDifermMass->Sumw2();
  TProfile* hmesq_jetPtoEScale = new TProfile("P_MomentumToEnergy", "", nbins, binning); hmesq_jetPtoEScale->Sumw2();


  TTree* newtree=0;
  if (writeFinalTree){
    newtree = new TTree("FinalTree", "");
    newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
    newtree->Branch("mesq_jetPtoEScale", &mesq_jetPtoEScale);
    newtree->Branch("ZZMass", &mzz);
  }

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (unsigned int bin=0; bin<binList.size(); bin++){
    cout << "Bin " << bin << " is now being scrutinized..." << endl;
    for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
      int getEv = binList.at(bin).events.at(ev);
      getEntry(treeList, getEv);
      if (ev%1000==0) cout << "Doing event " << getEv << endl;
      if (JetPt->size()<2 || JetEta->size()<2 || JetPhi->size()<2 || JetMass->size()<2){
        cerr << "Jet array sizes are less than Njets!" << endl;
        continue;
      }
      TTree* tree = findTree(treeList, getEv);
      wgt = fabs(genHEPMCweight*nGenMap[tree].first/nGenMap[tree].second);

      TLorentzVector jet[2], higgs;
      for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(JetPt->at(ij), JetEta->at(ij), JetPhi->at(ij), JetMass->at(ij));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));
      associated.push_back(SimpleParticle_t(0, jet[1]));

      SimpleParticleCollection_t daughters;
      for (int id=0; id<4; id++){
        double mass=0;
        if (abs(LepLepId->at(id))==13) mass = 0.105658;
        else if (abs(LepLepId->at(id))==11) mass = 0.000511;
        TLorentzVector pDaughter;
        pDaughter.SetPtEtaPhiM(LepPt->at(id), LepEta->at(id), LepPhi->at(id), mass);
        daughters.push_back(SimpleParticle_t(LepLepId->at(id), pDaughter));
      }
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      mela.setProcess(proc, me, prod);

      double alphasVal;
      if (recalculate){
        TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
        mela.computeProdP(mesq_conserveDifermMass, false);
        alphasVal = mela.getIORecord()->getAlphaSatMZ();
        mesq_conserveDifermMass /= pow(alphasVal, 4);
      }
      else{
        mesq_calc /= cconst_calc;
        mesq_conserveDifermMass=mesq_calc;
      }
      TUtil::setJetMassScheme(TVar::MomentumToEnergy);
      mela.computeProdP(mesq_jetPtoEScale, false);
      alphasVal = mela.getIORecord()->getAlphaSatMZ();
      mesq_jetPtoEScale /= pow(alphasVal, 4);

      bool doFill = !(
        isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
        ||
        isnan(mesq_jetPtoEScale) || isinf(mesq_jetPtoEScale)
        );

      if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, mesq_jetPtoEScale, wgt);

      mela.resetInputEvent();
    }

    binList.at(bin).sift(); binList.at(bin).adjustWeights();

    for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
      mzz = binList.at(bin).masses.at(ev);
      mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
      mesq_jetPtoEScale = binList.at(bin).me2vals.at(ev);
      wgt = binList.at(bin).weights.at(ev);
      hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass, wgt);
      hmesq_jetPtoEScale->Fill(mzz, mesq_jetPtoEScale, wgt);
      hvar->Fill(mzz, mzz, wgt);
      if (writeFinalTree) newtree->Fill();
    }
  }

  double* xexyey[2][4];
  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) xexyey[inorm][ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[inorm][0][bin] = hvar->GetBinContent(bin+1);
      xexyey[inorm][1][bin] = hvar->GetBinError(bin+1);

      if (inorm==0) cout << "Bin " << bin << " x-center: " << xexyey[inorm][0][bin] << " +- " << xexyey[inorm][1][bin] << endl;

      if (inorm==0){
        xexyey[inorm][2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      }
      else{
        xexyey[inorm][2][bin] = hmesq_jetPtoEScale->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_jetPtoEScale->GetBinError(bin+1);
      }
      xexyey[inorm][3][bin] = log10(xexyey[inorm][3][bin])/xexyey[inorm][2][bin];
      xexyey[inorm][2][bin] = log10(xexyey[inorm][2][bin]);
    }
  }

  for (int inorm=0; inorm<2; inorm++){
    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[inorm][0], xexyey[inorm][2], xexyey[inorm][1], xexyey[inorm][3]);
    if (inorm==0) tg->SetName("tg_P_ConserveDifermionMass");
    else tg->SetName("tg_P_MomentumToEnergy");
    foutput->WriteTObject(tg);
    delete tg;
  }

  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) delete[] xexyey[inorm][ix];
  }
  foutput->WriteTObject(hmesq_jetPtoEScale);
  foutput->WriteTObject(hmesq_conserveDifermMass);
  foutput->WriteTObject(hvar);
  if (writeFinalTree) foutput->WriteTObject(newtree);
  if (writeFinalTree) delete newtree;
  delete hmesq_conserveDifermMass;
  delete hmesq_jetPtoEScale;
  delete hvar;
  foutput->Close();
  delete[] binning;
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}


/* SPECIFIC COMMENT: OUTPUT ME DIVIDED BY ALPHAS(MZ)**3 TO REMAIN INDEPENDENT OF PDF CHOICE TO FIRST APPROXIMATION */
void get_PAvgProfile_JHUGen_JQCD_HSMHiggs_13TeV(int sqrts=13, bool recalculate = true){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";
  bool writeFinalTree=true;

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::JQCD;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  short NJets30;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float> myJetPt;
  std::vector<float> myJetEta;
  std::vector<float> myJetPhi;
  std::vector<float> myJetMass;

  float mesq_calc=0., cconst_calc=1.;
  float mesq_conserveDifermMass=0;
  float mesq_jetPtoEScale=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  float genHEPMCweight;
  float wgt=1;

  TString cinput_main;
  if (sqrts==13) cinput_main = inputdir_13TeV;
  else return;

  vector<TString> strSamples = constructSamplesList("JJQCD", sqrts);

  unordered_map<TTree*, pair<float, float>> nGenMap;
  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  for (unsigned int is=0; is<strSamples.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          if (!recalculate && tree->GetBranchStatus("pConst_JQCD_SIG_ghg2_1_JHUGen_JECNominal")==0) recalculate=true;
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("genHEPMCweight", 1); tree->SetBranchAddress("genHEPMCweight", &genHEPMCweight);
          tree->SetBranchStatus("nCleanedJetsPt30", 1); tree->SetBranchAddress("nCleanedJetsPt30", &NJets30);
          tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
          tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
          tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
          tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
          tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
          tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("LepLepId", 1); tree->SetBranchAddress("LepLepId", &LepLepId);
          tree->SetBranchStatus("LepPt", 1); tree->SetBranchAddress("LepPt", &LepPt);
          tree->SetBranchStatus("LepEta", 1); tree->SetBranchAddress("LepEta", &LepEta);
          tree->SetBranchStatus("LepPhi", 1); tree->SetBranchAddress("LepPhi", &LepPhi);
          if (!recalculate){
            tree->SetBranchStatus("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", 1);
            tree->SetBranchAddress("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &mesq_calc);
            tree->SetBranchStatus("pConst_JQCD_SIG_ghg2_1_JHUGen_JECNominal", 1);
            tree->SetBranchAddress("pConst_JQCD_SIG_ghg2_1_JHUGen_JECNominal", &cconst_calc);
          }

          TH1F* htmp = (TH1F*)finput->Get(COUNTERS_NAME);
          pair<float, float> nsum(htmp->GetBinContent(1), htmp->GetBinContent(41)); // No PU reweighting
          nGenMap[tree]=nsum;

          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }
  cout << "NEntries = " << nEntries << endl;

  vector<pair<float, int>> index;
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    if (NJets30!=1) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    addByLowest(index, mzz, ev);
    ev_acc++;
  }

  float firstVal=index.at(0).first;
  float lastVal=index.at(index.size()-1).first;
  float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
  float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
  cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

  float divisor=85000;
  int nbins = index.size()/divisor;
  const int nbins_th=10/*50*/;
  while (nbins<nbins_th){
    if (divisor>1000) divisor -= 1000;
    else if (divisor>100) divisor -= 100;
    else break;
    nbins=index.size()/divisor;
  }
  cout << "nbins=" << nbins << endl;
  if (nbins<3) cerr << "Not enough bins!" << endl;
  vector<ExtBin> binList;
  float* binning = new float[nbins+1];
  binning[0]=infimum;
  binning[nbins]=supremum;
  int ev_stepsize = index.size()/nbins;
  cout << "Event step size: " << ev_stepsize << endl;
  cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
  for (int ix=1; ix<nbins; ix++){
    binning[ix]=(index[ix*ev_stepsize-1].first+index[ix*ev_stepsize].first)*0.5;
    ExtBin tmpbin;
    tmpbin.binlow = binning[ix-1];
    tmpbin.binhigh = binning[ix];
    for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[(ix-1)*ev_stepsize+bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
  }
  ExtBin tmpbin;
  tmpbin.binlow = binning[nbins-1];
  tmpbin.binhigh = binning[nbins];
  for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index.size(); bin++) tmpbin.events.push_back(index[bin].second);
  binList.push_back(tmpbin);
  cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
  cout << "Bin list has the following bins:" << endl;
  for (unsigned int ib=0; ib<binList.size(); ib++){
    cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
  }
  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  TProfile* hvar = new TProfile("candMass", "", nbins, binning); hvar->Sumw2();
  TProfile* hmesq_conserveDifermMass = new TProfile("P_ConserveDifermionMass", "", nbins, binning); hmesq_conserveDifermMass->Sumw2();
  TProfile* hmesq_jetPtoEScale = new TProfile("P_MomentumToEnergy", "", nbins, binning); hmesq_jetPtoEScale->Sumw2();


  TTree* newtree=0;
  if (writeFinalTree){
    newtree = new TTree("FinalTree", "");
    newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
    newtree->Branch("mesq_jetPtoEScale", &mesq_jetPtoEScale);
    newtree->Branch("ZZMass", &mzz);
  }

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (unsigned int bin=0; bin<binList.size(); bin++){
    cout << "Bin " << bin << " is now being scrutinized..." << endl;
    for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
      int getEv = binList.at(bin).events.at(ev);
      getEntry(treeList, getEv);
      if (ev%1000==0) cout << "Doing event " << getEv << endl;
      if (JetPt->size()<1 || JetEta->size()<1 || JetPhi->size()<1 || JetMass->size()<1){
        cerr << "Jet array sizes are less than Njets!" << endl;
        continue;
      }
      TTree* tree = findTree(treeList, getEv);
      wgt = fabs(genHEPMCweight*nGenMap[tree].first/nGenMap[tree].second);

      TLorentzVector jet, higgs;
      jet.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet));

      SimpleParticleCollection_t daughters;
      for (int id=0; id<4; id++){
        double mass=0;
        if (abs(LepLepId->at(id))==13) mass = 0.105658;
        else if (abs(LepLepId->at(id))==11) mass = 0.000511;
        TLorentzVector pDaughter;
        pDaughter.SetPtEtaPhiM(LepPt->at(id), LepEta->at(id), LepPhi->at(id), mass);
        daughters.push_back(SimpleParticle_t(LepLepId->at(id), pDaughter));
      }
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      mela.setProcess(proc, me, prod);

      double alphasVal;
      if (recalculate){
        TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
        mela.computeProdP(mesq_conserveDifermMass, false);
        alphasVal = mela.getIORecord()->getAlphaSatMZ();
        mesq_conserveDifermMass /= pow(alphasVal, 3);
      }
      else{
        mesq_calc /= cconst_calc;
        mesq_conserveDifermMass=mesq_calc;
      }
      TUtil::setJetMassScheme(TVar::MomentumToEnergy);
      mela.computeProdP(mesq_jetPtoEScale, false);
      alphasVal = mela.getIORecord()->getAlphaSatMZ();
      mesq_jetPtoEScale /= pow(alphasVal, 3);

      bool doFill = !(
        isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
        ||
        isnan(mesq_jetPtoEScale) || isinf(mesq_jetPtoEScale)
        );

      if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, mesq_jetPtoEScale, wgt);

      mela.resetInputEvent();
    }

    binList.at(bin).sift(); binList.at(bin).adjustWeights();

    for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
      mzz = binList.at(bin).masses.at(ev);
      mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
      mesq_jetPtoEScale = binList.at(bin).me2vals.at(ev);
      wgt = binList.at(bin).weights.at(ev);
      hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass, wgt);
      hmesq_jetPtoEScale->Fill(mzz, mesq_jetPtoEScale, wgt);
      hvar->Fill(mzz, mzz, wgt);
      if (writeFinalTree) newtree->Fill();
    }
  }

  double* xexyey[2][4];
  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) xexyey[inorm][ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[inorm][0][bin] = hvar->GetBinContent(bin+1);
      xexyey[inorm][1][bin] = hvar->GetBinError(bin+1);

      if (inorm==0) cout << "Bin " << bin << " x-center: " << xexyey[inorm][0][bin] << " +- " << xexyey[inorm][1][bin] << endl;

      if (inorm==0){
        xexyey[inorm][2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      }
      else{
        xexyey[inorm][2][bin] = hmesq_jetPtoEScale->GetBinContent(bin+1);
        xexyey[inorm][3][bin] = hmesq_jetPtoEScale->GetBinError(bin+1);
      }
      xexyey[inorm][3][bin] = log10(xexyey[inorm][3][bin])/xexyey[inorm][2][bin];
      xexyey[inorm][2][bin] = log10(xexyey[inorm][2][bin]);
    }
  }

  for (int inorm=0; inorm<2; inorm++){
    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[inorm][0], xexyey[inorm][2], xexyey[inorm][1], xexyey[inorm][3]);
    if (inorm==0) tg->SetName("tg_P_ConserveDifermionMass");
    else tg->SetName("tg_P_MomentumToEnergy");
    foutput->WriteTObject(tg);
    delete tg;
  }

  for (int inorm=0; inorm<2; inorm++){
    for (int ix=0; ix<4; ix++) delete[] xexyey[inorm][ix];
  }
  foutput->WriteTObject(hmesq_jetPtoEScale);
  foutput->WriteTObject(hmesq_conserveDifermMass);
  foutput->WriteTObject(hvar);
  if (writeFinalTree) foutput->WriteTObject(newtree);
  if (writeFinalTree) delete newtree;
  delete hmesq_conserveDifermMass;
  delete hmesq_jetPtoEScale;
  delete hvar;
  foutput->Close();
  delete[] binning;
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}

/*
SPECIFIC COMMENTS:
OUTPUT ME DIVIDED BY
- H(1) PROPAGATOR
- Ideal mJJ propagator is taken out in WH or ZH.
- Reco mJJ propagator is multiplied.
*/
void get_PAvgProfile_JHUGen_HadVH_HSMHiggs_13TeV(TString strprod, int sqrts=13, bool recalculate = true){
  if (!(strprod == "Had_ZH" || strprod == "Had_WH")) return;
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";
  bool writeFinalTree=true;

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod;
  for (int iprod=(int)TVar::JJVBF; iprod<(int)TVar::nProductions; iprod++){
    prod = (TVar::Production)iprod;
    if (TVar::ProductionName(prod)==strprod) break;
  }
  TString strsamples="";
  //TString strsamples2="";
  if (strprod == "Had_ZH") strsamples="ZH";
  else if (strprod == "Had_WH") strsamples="WH";
  //if (strprod == "Had_ZH") strsamples2="WH";
  //else if (strprod == "Had_WH") strsamples2="ZH";

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  short NJets30;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float> myJetPt;
  std::vector<float> myJetEta;
  std::vector<float> myJetPhi;
  std::vector<float> myJetMass;

  float mesq_calc=0., cconst_calc=1., pmavjj=1.;
  float mesq_conserveDifermMass=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  short Z1Flav, Z2Flav;
  float genxsec=0, genBR=0, ngen=0, genHEPMCweight=0;
  float wgt=1;
  float rewgt=1;
  bool doRewgt=false;
  bool recomputePmavjj=false;
  TString strrewgtbranch="";
  TString strrecalcbranch="";
  TString strrecalcconstbranch="";
  TString strpmavjjbranch="";

  TString cinput_main;
  if (sqrts==13) cinput_main = inputdir_13TeV;
  else return;

  vector<TString> strSamples = constructSamplesList(strsamples, sqrts);
  //vector<TString> strSamples2 = constructSamplesList(strsamples2, sqrts);
  //appendVector(strSamples, strSamples2);

  unordered_map<TTree*, pair<float, float>> nGenMap;
  unordered_map<TTree*, float> mass_map;
  unordered_map<float, vector<TTree*>> mass_sample_map;

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  for (unsigned int is=0; is<strSamples.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          if (!recalculate){
            if (prod==TVar::Had_ZH && tree->GetBranchStatus("pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal")!=0){
              strrecalcbranch = "p_HadZH_SIG_ghz1_1_JHUGen_JECNominal";
              strrecalcconstbranch = "pConst_HadZH_SIG_ghz1_1_JHUGen_JECNominal";
            }
            else if (prod==TVar::Had_WH && tree->GetBranchStatus("pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal")!=0){
              strrecalcbranch = "p_HadWH_SIG_ghw1_1_JHUGen_JECNominal";
              strrecalcconstbranch = "pConst_HadWH_SIG_ghw1_1_JHUGen_JECNominal";
            }
            else recalculate=true;
            if (prod==TVar::Had_ZH && tree->GetBranchStatus("p_HadZH_mavjj_JECNominal")!=0) strpmavjjbranch="p_HadZH_mavjj_JECNominal";
            else if (prod==TVar::Had_WH && tree->GetBranchStatus("p_HadWH_mavjj_JECNominal")!=0) strpmavjjbranch="p_HadWH_mavjj_JECNominal";
            else recomputePmavjj=true;
          }
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("genxsec", 1); tree->SetBranchAddress("genxsec", &genxsec); tree->SetBranchStatus("genBR", 1); tree->SetBranchAddress("genBR", &genBR);
          tree->SetBranchStatus("genHEPMCweight", 1); tree->SetBranchAddress("genHEPMCweight", &genHEPMCweight);
          tree->SetBranchStatus("nCleanedJetsPt30", 1); tree->SetBranchAddress("nCleanedJetsPt30", &NJets30);
          tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
          tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
          tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
          tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
          tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
          tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("LepLepId", 1); tree->SetBranchAddress("LepLepId", &LepLepId);
          tree->SetBranchStatus("LepPt", 1); tree->SetBranchAddress("LepPt", &LepPt);
          tree->SetBranchStatus("LepEta", 1); tree->SetBranchAddress("LepEta", &LepEta);
          tree->SetBranchStatus("LepPhi", 1); tree->SetBranchAddress("LepPhi", &LepPhi);
          tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
          tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
          if (!recalculate){
            tree->SetBranchStatus(strrecalcbranch, 1); tree->SetBranchAddress(strrecalcbranch, &mesq_calc);
            tree->SetBranchStatus(strrecalcconstbranch, 1); tree->SetBranchAddress(strrecalcconstbranch, &cconst_calc);
            cout << "Extracting ME from " << strrecalcbranch << " and const from " << strrecalcconstbranch << endl;
          }
          else cout << "Recalculating the ME" << endl;
          if (!recomputePmavjj && strpmavjjbranch!=""){
            tree->SetBranchStatus(strpmavjjbranch, 1);
            tree->SetBranchAddress(strpmavjjbranch, &pmavjj);
            cout << "Extracting P_V(mJJ) from " << strpmavjjbranch << endl;
          }
          else if (recomputePmavjj) cout << "Recomputing P_V(mJJ)" << endl;
          else cout << "No valid P_V(mJJ)" << endl;
          tree->GetEntry(0);
          cout << "Cross section = " << genxsec*genBR << endl;

          TH1F* htmp = (TH1F*)finput->Get(COUNTERS_NAME);
          pair<float, float> nsum(htmp->GetBinContent(1), htmp->GetBinContent(41)); // No PU reweighting
          nGenMap[tree]=nsum;

          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);

          float polemass = findPoleMass(strSamples[is]);
          cout << "Pole mass = " << polemass << endl;
          mass_map[tree]=polemass;
          if (mass_sample_map.find(polemass)==mass_sample_map.end()){
            cout << "Inserting new pole mass sample array" << endl;
            vector<TTree*> dumarr;
            mass_sample_map[polemass] = dumarr;
          }
          mass_sample_map[polemass].push_back(tree);
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }
  cout << "NEntries = " << nEntries << endl;

  for (auto it = mass_sample_map.begin(); it != mass_sample_map.end(); ++it){
    float sum_ngen=0;
    float sum_xsec=0;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      it->second.at(ix)->GetEntry(0);
      sum_ngen += nGenMap[it->second.at(ix)].first;
      sum_xsec += genxsec*genBR;
    }
    float overallWeight = sum_ngen/sum_xsec;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      cout << "Sum Hep MC weights in tree " << ix << " / " << it->second.size() << " was " << nGenMap[it->second.at(ix)].second << " over " << nGenMap[it->second.at(ix)].first << " total gen. events." << endl;
      nGenMap[it->second.at(ix)].first = overallWeight/nGenMap[it->second.at(ix)].second;
      cout << "Event scale for tree " << ix << " / " << it->second.size() << " at pole mass " << it->first << " = " << nGenMap[it->second.at(ix)].first << endl;
    }
  }

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  vector<pair<float, int>> index;
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=NJets30>=2;
    if (!doProcess) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    unsigned int ic = (Z1Flav*Z2Flav==pow(13, 4))*0 + (Z1Flav*Z2Flav==pow(11, 4))*1 + (Z1Flav*Z2Flav==pow(11*13, 2))*2;
    addByLowest(index, mzz, ev);
    ev_acc++;
  }

  float firstVal=index.at(0).first;
  float lastVal=index.at(index.size()-1).first;
  float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
  float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
  cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

  float divisor=30000;
  //float divisor=60000;
  if (prod==TVar::Had_WH) divisor=45000;
  int nbins = index.size()/divisor;
  const int nbins_th=8/*50*/;
  while (nbins<nbins_th){
    if (divisor>1000) divisor -= 1000;
    else if (divisor>100) divisor -= 100;
    else break;
    nbins=index.size()/divisor;
  }
  cout << "nbins=" << nbins << endl;
  if (nbins<3) cerr << "Not enough bins!" << endl;
  vector<ExtBin> binList;
  float* binning = new float[nbins+1];
  binning[0]=infimum;
  binning[nbins]=supremum;
  int ev_stepsize = index.size()/nbins;
  cout << "Event step size: " << ev_stepsize << endl;
  cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
  for (int ix=1; ix<nbins; ix++){
    binning[ix]=(index[ix*ev_stepsize-1].first+index[ix*ev_stepsize].first)*0.5;
    ExtBin tmpbin;
    tmpbin.binlow = binning[ix-1];
    tmpbin.binhigh = binning[ix];
    for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[(ix-1)*ev_stepsize+bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
  }
  ExtBin tmpbin;
  tmpbin.binlow = binning[nbins-1];
  tmpbin.binhigh = binning[nbins];
  for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index.size(); bin++) tmpbin.events.push_back(index[bin].second);
  binList.push_back(tmpbin);
  cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
  cout << "Bin list has the following bins:" << endl;
  for (unsigned int ib=0; ib<binList.size(); ib++){
    cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
  }

  TProfile* hvar = new TProfile("candMass", "", nbins, binning); hvar->Sumw2();
  TProfile* hmesq_conserveDifermMass = new TProfile("P_ConserveDifermionMass", "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

  TTree* newtree=0;
  if (writeFinalTree){
    newtree = new TTree("FinalTree", "");
    newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
    newtree->Branch("ZZMass", &mzz);
    newtree->Branch("weight", &wgt);
  }

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  double aL1=0, aR1=0, aL2=0, aR2=0;
  double mh=0, gah=0;
  double mv=0, gav=0;
  ev_acc=0;
  for (unsigned int bin=0; bin<binList.size(); bin++){
    cout << "Bin " << bin << " is now being scrutinized..." << endl;
    for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
      int getEv = binList.at(bin).events.at(ev);
      getEntry(treeList, getEv);
      if (ev%1000==0) cout << "Doing event " << getEv << endl;
      if (JetPt->size()<2 || JetEta->size()<2 || JetPhi->size()<2 || JetMass->size()<2){
        cerr << "Jet array sizes are less than Njets!" << endl;
        continue;
      }
      TTree* tree = findTree(treeList, getEv);
      wgt = fabs(genxsec*genBR*genHEPMCweight*nGenMap[tree].first*rewgt);

      TLorentzVector jet[2], higgs;
      for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(JetPt->at(ij), JetEta->at(ij), JetPhi->at(ij), JetMass->at(ij));
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));
      associated.push_back(SimpleParticle_t(0, jet[1]));

      SimpleParticleCollection_t daughters;
      for (int id=0; id<4; id++){
        double mass=0;
        if (abs(LepLepId->at(id))==13) mass = 0.105658;
        else if (abs(LepLepId->at(id))==11) mass = 0.000511;
        TLorentzVector pDaughter;
        pDaughter.SetPtEtaPhiM(LepPt->at(id), LepEta->at(id), LepPhi->at(id), mass);
        daughters.push_back(SimpleParticle_t(LepLepId->at(id), pDaughter));
      }
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      mela.setProcess(proc, me, prod);

      double propagator=1;
      double propagatorV=1;
      if (recalculate || ev_acc==0){
        TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
        mela.computeProdP(mesq_conserveDifermMass, false);
        if (!recalculate){
          mesq_calc /= cconst_calc;
          mesq_conserveDifermMass=mesq_calc;

          mh = mass_map[tree];
          gah = mela.getHiggsWidthAtPoleMass(mh);
        }
        else{
          mela.getIORecord()->getHiggsMassWidth(mh, gah, 0);
        }
        mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
        mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
        if (prod==TVar::Had_ZH){
          mv = mela.getPrimaryMass(23);
          gav = mela.getPrimaryWidth(23);
        }
        else if (prod==TVar::Had_WH){
          mv = mela.getPrimaryMass(24);
          gav = mela.getPrimaryWidth(24);
        }
      }
      else{
        mesq_calc /= cconst_calc;
        mesq_conserveDifermMass=mesq_calc;
      }
      propagator = 1./(pow(pow(mzz, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
      mesq_conserveDifermMass /= propagator;
      if (prod==TVar::Had_ZH || prod==TVar::Had_WH){
        double mjj = (jet[0] + jet[1]).M();
        propagatorV = 1./(pow(pow(mjj, 2)-pow(mv, 2), 2) + pow(mv*gav, 2));
        mesq_conserveDifermMass /= propagatorV;
      }
      if (recomputePmavjj) mela.computeDijetConvBW(pmavjj, false);
      mesq_conserveDifermMass *= pmavjj;

      bool doFill = !(
        isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
        );

      if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0, wgt);

      mela.resetInputEvent();
      ev_acc++;
    }

    binList.at(bin).sift(); binList.at(bin).adjustWeights(10);

    for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
      mzz = binList.at(bin).masses.at(ev);
      mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
      wgt = binList.at(bin).weights.at(ev);
      hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
      hvar->Fill(mzz, mzz);
      if (writeFinalTree) newtree->Fill();
    }
  }
  cout << "Total events accumulated: " << ev_acc << endl;

  double* xexyey[4];
  for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
  for (int bin=0; bin<nbins; bin++){
    xexyey[0][bin] = hvar->GetBinContent(bin+1);
    xexyey[1][bin] = hvar->GetBinError(bin+1);

    cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;

    xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
    xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
    xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
    xexyey[2][bin] = log10(xexyey[2][bin]);
  }

  TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
  tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
  foutput->WriteTObject(tg);
  delete tg;

  for (int ix=0; ix<4; ix++) delete[] xexyey[ix];

  foutput->WriteTObject(hmesq_conserveDifermMass);
  foutput->WriteTObject(hvar);
  if (writeFinalTree) foutput->WriteTObject(newtree);
  if (writeFinalTree) delete newtree;
  delete hmesq_conserveDifermMass;
  delete hvar;
  delete[] binning;

  foutput->Close();
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}


/*
SPECIFIC COMMENT: OUTPUT ME DIVIDED BY
- H(1) PROPAGATOR
- (aL1**2+aR1**2)*(aL2**2+aR2**2) TO REMAIN INDEPENDENT OF CHANNEL
- Ideal mJJ propagator is taken out in WH or ZH.
- Reco mJJ propagator is multiplied.
*/
void get_PAvgProfile_MCFM_JJPROD_S_HSMHiggs_13TeV(TString strprod, int sqrts=13, bool recalculate = true){
  if (!(strprod == "Had_ZH_S" || strprod == "Had_WH_S" || strprod == "JJVBF_S")) return;
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";
  bool writeFinalTree=true;

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod;
  for (int iprod=(int) TVar::ZZGG; iprod<(int) TVar::nProductions; iprod++){
    prod = (TVar::Production)iprod;
    if (TVar::ProductionName(prod)==strprod) break;
  }
  TString strsamples;
  if (strprod == "Had_ZH_S") strsamples="ZH";
  else if (strprod == "Had_WH_S") strsamples="WH";
  else if (strprod == "JJVBF_S") strsamples="JJVBF";

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  short NJets30;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float> myJetPt;
  std::vector<float> myJetEta;
  std::vector<float> myJetPhi;
  std::vector<float> myJetMass;

  float mesq_calc=0., cconst_calc=1., pmavjj=1.;
  float mesq_conserveDifermMass=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  short Z1Flav, Z2Flav;
  float genxsec=0, genBR=0, ngen=0, genHEPMCweight=0;
  float wgt=1;
  float rewgt=1;
  bool doRewgt=false;
  bool recomputePmavjj=false;
  TString strrewgtbranch="";
  TString strrecalcbranch="";
  TString strrecalcconstbranch="";
  TString strpmavjjbranch="";

  TString cinput_main;
  if (sqrts==13) cinput_main = inputdir_13TeV;
  else return;

  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  vector<TString> strSamples = constructSamplesList(strsamples, sqrts);

  unordered_map<TTree*, pair<float, float>> nGenMap;
  unordered_map<TTree*, float> mass_map;
  unordered_map<float, vector<TTree*>> mass_sample_map;

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  for (unsigned int is=0; is<strSamples.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          if (!recalculate){
            if (prod==TVar::Had_ZH_S && tree->GetBranchStatus("pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal")!=0){
              strrecalcbranch = "p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal";
              strrecalcconstbranch = "pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal";
            }
            else if (prod==TVar::Had_WH_S && tree->GetBranchStatus("pConst_HadWH_S_SIG_ghv1_1_MCFM_JECNominal")!=0){
              strrecalcbranch = "p_HadWH_S_SIG_ghv1_1_MCFM_JECNominal";
              strrecalcconstbranch = "pConst_HadWH_S_SIG_ghv1_1_MCFM_JECNominal";
            }
            else if (prod==TVar::JJVBF_S && tree->GetBranchStatus("pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal")!=0){
              strrecalcbranch = "p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal";
              strrecalcconstbranch = "pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal";
            }
            else recalculate=true;
            if (prod==TVar::Had_ZH_S && tree->GetBranchStatus("p_HadZH_mavjj_JECNominal")!=0) strpmavjjbranch="p_HadZH_mavjj_JECNominal";
            else if (prod==TVar::Had_WH_S && tree->GetBranchStatus("p_HadWH_mavjj_JECNominal")!=0) strpmavjjbranch="p_HadWH_mavjj_JECNominal";
            else if (prod!=TVar::JJVBF_S) recomputePmavjj=true;
          }
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("genxsec", 1); tree->SetBranchAddress("genxsec", &genxsec); tree->SetBranchStatus("genBR", 1); tree->SetBranchAddress("genBR", &genBR);
          tree->SetBranchStatus("genHEPMCweight", 1); tree->SetBranchAddress("genHEPMCweight", &genHEPMCweight);
          tree->SetBranchStatus("nCleanedJetsPt30", 1); tree->SetBranchAddress("nCleanedJetsPt30", &NJets30);
          tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
          tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
          tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
          tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
          tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
          tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("LepLepId", 1); tree->SetBranchAddress("LepLepId", &LepLepId);
          tree->SetBranchStatus("LepPt", 1); tree->SetBranchAddress("LepPt", &LepPt);
          tree->SetBranchStatus("LepEta", 1); tree->SetBranchAddress("LepEta", &LepEta);
          tree->SetBranchStatus("LepPhi", 1); tree->SetBranchAddress("LepPhi", &LepPhi);
          tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
          tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
          if (!recalculate){
            tree->SetBranchStatus(strrecalcbranch, 1); tree->SetBranchAddress(strrecalcbranch, &mesq_calc);
            tree->SetBranchStatus(strrecalcconstbranch, 1); tree->SetBranchAddress(strrecalcconstbranch, &cconst_calc);
            cout << "Extracting ME from " << strrecalcbranch << " and const from " << strrecalcconstbranch << endl;
          }
          else cout << "Recalculating the ME" << endl;
          if (!recomputePmavjj && strpmavjjbranch!=""){
            tree->SetBranchStatus(strpmavjjbranch, 1);
            tree->SetBranchAddress(strpmavjjbranch, &pmavjj);
            cout << "Extracting P_V(mJJ) from " << strpmavjjbranch << endl;
          }
          else if (recomputePmavjj) cout << "Recomputing P_V(mJJ)" << endl;
          else cout << "No valid P_V(mJJ)" << endl;

          tree->GetEntry(0);
          cout << "Cross section = " << genxsec*genBR << endl;

          TH1F* htmp = (TH1F*)finput->Get(COUNTERS_NAME);
          pair<float, float> nsum(htmp->GetBinContent(1), htmp->GetBinContent(41)); // No PU reweighting
          nGenMap[tree]=nsum;

          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);

          float polemass = findPoleMass(strSamples[is]);
          cout << "Pole mass = " << polemass << endl;
          mass_map[tree]=polemass;
          if (mass_sample_map.find(polemass)==mass_sample_map.end()){
            cout << "Inserting new pole mass sample array" << endl;
            vector<TTree*> dumarr;
            mass_sample_map[polemass] = dumarr;
          }
          mass_sample_map[polemass].push_back(tree);
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }
  cout << "NEntries = " << nEntries << endl;

  for (auto it = mass_sample_map.begin(); it != mass_sample_map.end(); ++it){
    float sum_ngen=0;
    float sum_xsec=0;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      it->second.at(ix)->GetEntry(0);
      sum_ngen += nGenMap[it->second.at(ix)].first;
      sum_xsec += genxsec*genBR;
    }
    float overallWeight = sum_ngen/sum_xsec;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      cout << "Sum Hep MC weights in tree " << ix << " / " << it->second.size() << " was " << nGenMap[it->second.at(ix)].second << " over " << nGenMap[it->second.at(ix)].first << " total gen. events." << endl;
      nGenMap[it->second.at(ix)].first = overallWeight/nGenMap[it->second.at(ix)].second;
      cout << "Event scale for tree " << ix << " / " << it->second.size() << " at pole mass " << it->first << " = " << nGenMap[it->second.at(ix)].first << endl;
    }
  }

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  vector<pair<float, int>> index[3];
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=
      NJets30>=2 && (
      Z1Flav*Z2Flav==pow(13, 4)
      ||
      Z1Flav*Z2Flav==pow(11, 4)
      ||
      Z1Flav*Z2Flav==pow(11*13, 2)
      )
      ;
    if (!doProcess) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    unsigned int ic = (Z1Flav*Z2Flav==pow(13, 4))*0 + (Z1Flav*Z2Flav==pow(11, 4))*1 + (Z1Flav*Z2Flav==pow(11*13, 2))*2;
    addByLowest(index[ic], mzz, ev);
    if (ic<2) addByLowest(index[1-ic], mzz, ev);
    ev_acc++;
  }

  for (unsigned int ic=0; ic<3; ic++){
    float firstVal=index[ic].at(0).first;
    float lastVal=index[ic].at(index[ic].size()-1).first;
    float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

    float divisor=15000;
    if (prod==TVar::Had_ZH_S || prod==TVar::Had_WH_S) divisor=12000;
    int nbins = index[ic].size()/divisor;
    const int nbins_th=10/*50*/;
    while (nbins<nbins_th){
      if (divisor>1000) divisor -= 1000;
      else if (divisor>100) divisor -= 100;
      else break;
      nbins=index[ic].size()/divisor;
    }
    cout << "nbins=" << nbins << endl;
    if (nbins<3) cerr << "Not enough bins!" << endl;
    vector<ExtBin> binList;
    float* binning = new float[nbins+1];
    binning[0]=infimum;
    binning[nbins]=supremum;
    int ev_stepsize = index[ic].size()/nbins;
    cout << "Event step size: " << ev_stepsize << endl;
    cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
    for (int ix=1; ix<nbins; ix++){
      binning[ix]=(index[ic][ix*ev_stepsize-1].first+index[ic][ix*ev_stepsize].first)*0.5;
      ExtBin tmpbin;
      tmpbin.binlow = binning[ix-1];
      tmpbin.binhigh = binning[ix];
      for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[ic][(ix-1)*ev_stepsize+bin].second);
      binList.push_back(tmpbin);
      cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ic][ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
    }
    ExtBin tmpbin;
    tmpbin.binlow = binning[nbins-1];
    tmpbin.binhigh = binning[nbins];
    for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index[ic].size(); bin++) tmpbin.events.push_back(index[ic][bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
    cout << "Bin list has the following bins:" << endl;
    for (unsigned int ib=0; ib<binList.size(); ib++){
      cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    }

    TProfile* hvar = new TProfile(Form("candMass_%s", strchannel[ic].Data()), "", nbins, binning); hvar->Sumw2();
    TProfile* hmesq_conserveDifermMass = new TProfile(Form("P_ConserveDifermionMass_%s", strchannel[ic].Data()), "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

    TTree* newtree=0;
    if (writeFinalTree){
      newtree = new TTree(Form("FinalTree_%s", strchannel[ic].Data()), "");
      newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
      newtree->Branch("ZZMass", &mzz);
      newtree->Branch("weight", &wgt);
    }

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

    double aL1=0, aR1=0, aL2=0, aR2=0;
    double mh=0, gah=0;
    double mv=0, gav=0;
    unsigned int ev_acc=0;
    for (unsigned int bin=0; bin<binList.size(); bin++){
      cout << "Bin " << bin << " is now being scrutinized..." << endl;
      for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
        int getEv = binList.at(bin).events.at(ev);
        getEntry(treeList, getEv);
        if (ev%1000==0) cout << "Doing event " << getEv << endl;
        if (JetPt->size()<2 || JetEta->size()<2 || JetPhi->size()<2 || JetMass->size()<2){
          cerr << "Jet array sizes are less than Njets!" << endl;
          continue;
        }
        TTree* tree = findTree(treeList, getEv);
        wgt = fabs(genxsec*genBR*genHEPMCweight*nGenMap[tree].first*rewgt);

        TLorentzVector jet[2], higgs;
        for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(JetPt->at(ij), JetEta->at(ij), JetPhi->at(ij), JetMass->at(ij));
        higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
        TVector3 boostH = higgs.BoostVector();

        SimpleParticleCollection_t associated;
        associated.push_back(SimpleParticle_t(0, jet[0]));
        associated.push_back(SimpleParticle_t(0, jet[1]));

        SimpleParticleCollection_t daughters;
        for (int id=0; id<4; id++){
          double mass=0;
          if (abs(LepLepId->at(id))==13) mass = 0.105658;
          else if (abs(LepLepId->at(id))==11) mass = 0.000511;
          TLorentzVector pDaughter;
          pDaughter.SetPtEtaPhiM(LepPt->at(id), LepEta->at(id), LepPhi->at(id), mass);
          daughters.push_back(SimpleParticle_t(LepLepId->at(id), pDaughter));
        }
        mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

        mela.setProcess(proc, me, prod);

        double propagator=1;
        double propagatorV=1;
        if (recalculate || ev_acc==0){
          TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
          mela.computeProdP(mesq_conserveDifermMass, false);
          if (!recalculate){
            mesq_calc /= cconst_calc;
            mesq_conserveDifermMass=mesq_calc;

            mh = mPOLE;
            gah = mela.getHiggsWidthAtPoleMass(mh);
          }
          else{
            mela.getIORecord()->getHiggsMassWidth(mh, gah, 0);
          }
          mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
          mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
          if (prod==TVar::Had_ZH_S){
            mv = mela.getPrimaryMass(23);
            gav = mela.getPrimaryWidth(23);
          }
          else if (prod==TVar::Had_WH_S){
            mv = mela.getPrimaryMass(24);
            gav = mela.getPrimaryWidth(24);
          }
        }
        else{
          mesq_calc /= cconst_calc;
          mesq_conserveDifermMass=mesq_calc;
        }
        propagator = 1./(pow(pow(mzz, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
        mesq_conserveDifermMass /= propagator;
        if (fabs(aL1)>0. || fabs(aR1)>0.) mesq_conserveDifermMass /= pow(aL1, 2)+pow(aR1, 2);
        if (fabs(aL2)>0. || fabs(aR2)>0.) mesq_conserveDifermMass /= pow(aL2, 2)+pow(aR2, 2);
        if (prod==TVar::Had_ZH_S || prod==TVar::Had_WH_S){
          double mjj = (jet[0] + jet[1]).M();
          propagatorV = 1./(pow(pow(mjj, 2)-pow(mv, 2), 2) + pow(mv*gav, 2));
          mesq_conserveDifermMass /= propagatorV;
          if (recomputePmavjj){
            mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, (prod==TVar::Had_ZH_S ? TVar::Had_ZH : TVar::Had_WH));
            mela.computeDijetConvBW(pmavjj, false);
          }
          mesq_conserveDifermMass *= pmavjj;
        }

        bool doFill = !(
          isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
          );

        if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0, wgt);

        mela.resetInputEvent();
        ev_acc++;
      }

      binList.at(bin).sift(); binList.at(bin).adjustWeights(10);

      for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
        mzz = binList.at(bin).masses.at(ev);
        mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
        wgt = binList.at(bin).weights.at(ev);
        hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
        hvar->Fill(mzz, mzz);
        if (writeFinalTree) newtree->Fill();
      }
    }
    cout << "Total events accumulated: " << ev_acc << endl;

    double* xexyey[4];
    for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[0][bin] = hvar->GetBinContent(bin+1);
      xexyey[1][bin] = hvar->GetBinError(bin+1);

      cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;

      xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
      xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
      xexyey[2][bin] = log10(xexyey[2][bin]);
    }

    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
    foutput->WriteTObject(tg);
    delete tg;

    for (int ix=0; ix<4; ix++) delete[] xexyey[ix];

    foutput->WriteTObject(hmesq_conserveDifermMass);
    foutput->WriteTObject(hvar);
    if (writeFinalTree) foutput->WriteTObject(newtree);
    if (writeFinalTree) delete newtree;
    delete hmesq_conserveDifermMass;
    delete hvar;
    delete[] binning;
  }
  foutput->Close();
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}

/*
SPECIFIC COMMENT: OUTPUT ME DIVIDED BY
- (aL1**2+aR1**2)*(aL2**2+aR2**2) TO REMAIN INDEPENDENT OF CHANNEL
- PROP(Z)/LINE_PROP(Z) AROUND M4L~=MZ
- Ideal mJJ propagator is taken out in WZZ or ZZZ.
- Reco mJJ propagator is multiplied.
*/
void get_PAvgProfile_MCFM_JJPROD_bkgZZ_13TeV(TString strprod, int sqrts=13, bool recalculate = true){
  if (!(strprod == "Had_ZH" || strprod == "Had_WH" || strprod == "JJVBF")) return;
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";
  bool writeFinalTree=true;

  TVar::Process proc = TVar::bkgZZ;
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod;
  for (int iprod=(int) TVar::ZZGG; iprod<(int) TVar::nProductions; iprod++){
    prod = (TVar::Production)iprod;
    if (TVar::ProductionName(prod)==strprod) break;
  }
  TString strsamples;
  if (strprod == "Had_ZH") strsamples="ZH";
  else if (strprod == "Had_WH") strsamples="WH";
  else if (strprod == "JJVBF") strsamples="JJVBF";

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  short NJets30;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float> myJetPt;
  std::vector<float> myJetEta;
  std::vector<float> myJetPhi;
  std::vector<float> myJetMass;

  float GenHMass=-1;
  std::vector<float>* LHEMotherPz=0;
  std::vector<float>* LHEMotherE=0;
  std::vector<short>* LHEMotherId=0;
  std::vector<float>* LHEDaughterPt=0;
  std::vector<float>* LHEDaughterEta=0;
  std::vector<float>* LHEDaughterPhi=0;
  std::vector<float>* LHEDaughterMass=0;
  std::vector<short>* LHEDaughterId=0;
  std::vector<float>* LHEAssociatedParticlePt=0;
  std::vector<float>* LHEAssociatedParticleEta=0;
  std::vector<float>* LHEAssociatedParticlePhi=0;
  std::vector<float>* LHEAssociatedParticleMass=0;
  std::vector<short>* LHEAssociatedParticleId=0;

  float mesq_calc=0., cconst_calc=1., pmavjj=1.;
  float mesq_conserveDifermMass=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  short Z1Flav, Z2Flav;
  float genxsec=0, genBR=0, ngen=0, genHEPMCweight=0;
  float wgt=1;
  float rewgt=1;
  bool doRewgt=false;
  bool recomputePmavjj=false;
  TString strrewgtbranch="";
  TString strrecalcbranch="";
  TString strrecalcconstbranch="";
  TString strpmavjjbranch="";

  TString cinput_main;
  if (sqrts==13) cinput_main = inputdir_13TeV;
  else return;

  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  vector<TString> strSamples = constructSamplesList(strsamples, sqrts);

  unordered_map<TTree*, pair<float, float>> nGenMap;
  unordered_map<TTree*, float> mass_map;
  unordered_map<float, vector<TTree*>> mass_sample_map;

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  for (unsigned int is=0; is<strSamples.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          if (!recalculate){
            if (prod==TVar::Had_ZH && tree->GetBranchStatus("pConst_HadZH_BKG_MCFM_JECNominal")!=0){
              strrecalcbranch = "p_HadZH_BKG_MCFM_JECNominal";
              strrecalcconstbranch = "pConst_HadZH_BKG_MCFM_JECNominal";
            }
            else if (prod==TVar::Had_WH && tree->GetBranchStatus("pConst_HadWH_BKG_MCFM_JECNominal")!=0){
              strrecalcbranch = "p_HadWH_BKG_MCFM_JECNominal";
              strrecalcconstbranch = "pConst_HadWH_BKG_MCFM_JECNominal";
            }
            else if (prod==TVar::JJVBF && tree->GetBranchStatus("pConst_JJVBF_BKG_MCFM_JECNominal")!=0){
              strrecalcbranch = "p_JJVBF_BKG_MCFM_JECNominal";
              strrecalcconstbranch = "pConst_JJVBF_BKG_MCFM_JECNominal";
            }
            else recalculate=true;
            if (prod==TVar::Had_ZH && tree->GetBranchStatus("p_HadZH_mavjj_JECNominal")!=0) strpmavjjbranch="p_HadZH_mavjj_JECNominal";
            else if (prod==TVar::Had_WH && tree->GetBranchStatus("p_HadWH_mavjj_JECNominal")!=0) strpmavjjbranch="p_HadWH_mavjj_JECNominal";
            else if (prod!=TVar::JJVBF) recomputePmavjj=true;
          }
          if (tree->GetBranchStatus("p_Gen_JJEW_BKG_MCFM")!=0) strrewgtbranch = "p_Gen_JJEW_BKG_MCFM";
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("genxsec", 1); tree->SetBranchAddress("genxsec", &genxsec); tree->SetBranchStatus("genBR", 1); tree->SetBranchAddress("genBR", &genBR);
          tree->SetBranchStatus("genHEPMCweight", 1); tree->SetBranchAddress("genHEPMCweight", &genHEPMCweight);
          tree->SetBranchStatus("GenHMass", 1); tree->SetBranchAddress("GenHMass", &GenHMass);
          tree->SetBranchStatus("nCleanedJetsPt30", 1); tree->SetBranchAddress("nCleanedJetsPt30", &NJets30);
          tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
          tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
          tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
          tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
          tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
          tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("LepLepId", 1); tree->SetBranchAddress("LepLepId", &LepLepId);
          tree->SetBranchStatus("LepPt", 1); tree->SetBranchAddress("LepPt", &LepPt);
          tree->SetBranchStatus("LepEta", 1); tree->SetBranchAddress("LepEta", &LepEta);
          tree->SetBranchStatus("LepPhi", 1); tree->SetBranchAddress("LepPhi", &LepPhi);
          tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
          tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
          // Do reweight signal to bkg.
          if (strrewgtbranch!=""){
            tree->SetBranchStatus(strrewgtbranch, 1); tree->SetBranchAddress(strrewgtbranch, &rewgt);
            cout << "Reweighting is read from branch " << strrewgtbranch << endl;
          }
          else{
            doRewgt=true;
            cout << "Reweighting needs to be computed." << endl;
            cerr << "NOT YET IMPLEMENTED!" << endl;
            assert(0);
          }
          if (!recalculate){
            tree->SetBranchStatus(strrecalcbranch, 1); tree->SetBranchAddress(strrecalcbranch, &mesq_calc);
            tree->SetBranchStatus(strrecalcconstbranch, 1); tree->SetBranchAddress(strrecalcconstbranch, &cconst_calc);
            cout << "Extracting ME from " << strrecalcbranch << " and const from " << strrecalcconstbranch << endl;
          }
          else cout << "Recalculating the ME" << endl;
          if (!recomputePmavjj && strpmavjjbranch!=""){
            tree->SetBranchStatus(strpmavjjbranch, 1);
            tree->SetBranchAddress(strpmavjjbranch, &pmavjj);
            cout << "Extracting P_V(mJJ) from " << strpmavjjbranch << endl;
          }
          else if (recomputePmavjj) cout << "Recomputing P_V(mJJ)" << endl;
          else cout << "No valid P_V(mJJ)" << endl;
          if (doRewgt){
            tree->SetBranchStatus("LHEMotherPz", 1); tree->SetBranchAddress("LHEMotherPz", &LHEMotherPz);
            tree->SetBranchStatus("LHEMotherE", 1); tree->SetBranchAddress("LHEMotherE", &LHEMotherE);
            tree->SetBranchStatus("LHEMotherId", 1); tree->SetBranchAddress("LHEMotherId", &LHEMotherId);

            tree->SetBranchStatus("LHEDaughterPt", 1); tree->SetBranchAddress("LHEDaughterPt", &LHEDaughterPt);
            tree->SetBranchStatus("LHEDaughterEta", 1); tree->SetBranchAddress("LHEDaughterEta", &LHEDaughterEta);
            tree->SetBranchStatus("LHEDaughterPhi", 1); tree->SetBranchAddress("LHEDaughterPhi", &LHEDaughterPhi);
            tree->SetBranchStatus("LHEDaughterMass", 1); tree->SetBranchAddress("LHEDaughterMass", &LHEDaughterMass);
            tree->SetBranchStatus("LHEDaughterId", 1); tree->SetBranchAddress("LHEDaughterId", &LHEDaughterId);

            tree->SetBranchStatus("LHEAssociatedParticlePt", 1); tree->SetBranchAddress("LHEAssociatedParticlePt", &LHEAssociatedParticlePt);
            tree->SetBranchStatus("LHEAssociatedParticleEta", 1); tree->SetBranchAddress("LHEAssociatedParticleEta", &LHEAssociatedParticleEta);
            tree->SetBranchStatus("LHEAssociatedParticlePhi", 1); tree->SetBranchAddress("LHEAssociatedParticlePhi", &LHEAssociatedParticlePhi);
            tree->SetBranchStatus("LHEAssociatedParticleMass", 1); tree->SetBranchAddress("LHEAssociatedParticleMass", &LHEAssociatedParticleMass);
            tree->SetBranchStatus("LHEAssociatedParticleId", 1); tree->SetBranchAddress("LHEAssociatedParticleId", &LHEAssociatedParticleId);
          }
          tree->GetEntry(0);
          cout << "Cross section = " << genxsec*genBR << endl;

          TH1F* htmp = (TH1F*)finput->Get(COUNTERS_NAME);
          pair<float, float> nsum(htmp->GetBinContent(1), htmp->GetBinContent(41)); // No PU reweighting
          nGenMap[tree]=nsum;

          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);

          float polemass = findPoleMass(strSamples[is]);
          cout << "Pole mass = " << polemass << endl;
          mass_map[tree]=polemass;
          if (mass_sample_map.find(polemass)==mass_sample_map.end()){
            cout << "Inserting new pole mass sample array" << endl;
            vector<TTree*> dumarr;
            mass_sample_map[polemass] = dumarr;
          }
          mass_sample_map[polemass].push_back(tree);
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }
  cout << "NEntries = " << nEntries << endl;

  for (auto it = mass_sample_map.begin(); it != mass_sample_map.end(); ++it){
    float sum_ngen=0;
    float sum_xsec=0;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      it->second.at(ix)->GetEntry(0);
      sum_ngen += nGenMap[it->second.at(ix)].first;
      sum_xsec += genxsec*genBR;
    }
    float overallWeight = sum_ngen/sum_xsec;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      cout << "Sum Hep MC weights in tree " << ix << " / " << it->second.size() << " was " << nGenMap[it->second.at(ix)].second << " over " << nGenMap[it->second.at(ix)].first << " total gen. events." << endl;
      nGenMap[it->second.at(ix)].first = overallWeight/nGenMap[it->second.at(ix)].second;
      cout << "Event scale for tree " << ix << " / " << it->second.size() << " at pole mass " << it->first << " = " << nGenMap[it->second.at(ix)].first << endl;
    }
  }

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  vector<pair<float, int>> index[3];
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=
      NJets30>=2 && (
      Z1Flav*Z2Flav==pow(13, 4)
      ||
      Z1Flav*Z2Flav==pow(11, 4)
      ||
      Z1Flav*Z2Flav==pow(11*13, 2)
      )
      ;
    if (!doProcess) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    unsigned int ic = (Z1Flav*Z2Flav==pow(13, 4))*0 + (Z1Flav*Z2Flav==pow(11, 4))*1 + (Z1Flav*Z2Flav==pow(11*13, 2))*2;
    addByLowest(index[ic], mzz, ev);
    if (ic<2) addByLowest(index[1-ic], mzz, ev);
    ev_acc++;
  }

  for (unsigned int ic=0; ic<3; ic++){
    float firstVal=index[ic].at(0).first;
    float lastVal=index[ic].at(index[ic].size()-1).first;
    float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

    float divisor=25000*(ic==2)+40000*(ic<2);
    if (prod==TVar::Had_ZH) divisor=15000;
    if (prod==TVar::Had_WH) divisor=15000;
    int nbins = index[ic].size()/divisor;
    int nbins_th=10/*50*/;
    if (prod==TVar::Had_WH) nbins_th=8;
    while (nbins<nbins_th){
      if (divisor>1000) divisor -= 1000;
      else if (divisor>100) divisor -= 100;
      else break;
      nbins=index[ic].size()/divisor;
    }
    cout << "nbins=" << nbins << endl;
    if (nbins<3) cerr << "Not enough bins!" << endl;
    vector<ExtBin> binList;
    float* binning = new float[nbins+1];
    binning[0]=infimum;
    binning[nbins]=supremum;
    int ev_stepsize = index[ic].size()/nbins;
    cout << "Event step size: " << ev_stepsize << endl;
    cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
    for (int ix=1; ix<nbins; ix++){
      binning[ix]=(index[ic][ix*ev_stepsize-1].first+index[ic][ix*ev_stepsize].first)*0.5;
      ExtBin tmpbin;
      tmpbin.binlow = binning[ix-1];
      tmpbin.binhigh = binning[ix];
      for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[ic][(ix-1)*ev_stepsize+bin].second);
      binList.push_back(tmpbin);
      cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ic][ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
    }
    ExtBin tmpbin;
    tmpbin.binlow = binning[nbins-1];
    tmpbin.binhigh = binning[nbins];
    for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index[ic].size(); bin++) tmpbin.events.push_back(index[ic][bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
    cout << "Bin list has the following bins:" << endl;
    for (unsigned int ib=0; ib<binList.size(); ib++){
      cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    }

    TProfile* hvar = new TProfile(Form("candMass_%s", strchannel[ic].Data()), "", nbins, binning); hvar->Sumw2();
    TProfile* hmesq_conserveDifermMass = new TProfile(Form("P_ConserveDifermionMass_%s", strchannel[ic].Data()), "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

    TTree* newtree=0;
    if (writeFinalTree){
      newtree = new TTree(Form("FinalTree_%s", strchannel[ic].Data()), "");
      newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
      newtree->Branch("ZZMass", &mzz);
      newtree->Branch("weight", &wgt);
    }

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

    double aL1=0, aR1=0, aL2=0, aR2=0;
    double mv=0, gav=0;
    unsigned int ev_acc=0;
    for (unsigned int bin=0; bin<binList.size(); bin++){
      cout << "Bin " << bin << " is now being scrutinized..." << endl;
      for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
        int getEv = binList.at(bin).events.at(ev);
        getEntry(treeList, getEv);
        if (ev%1000==0) cout << "Doing event " << getEv << endl;
        if (JetPt->size()<2 || JetEta->size()<2 || JetPhi->size()<2 || JetMass->size()<2){
          cerr << "Jet array sizes are less than Njets!" << endl;
          continue;
        }
        TTree* tree = findTree(treeList, getEv);

        // Divide out the signal H propagator in the reweighting denominator = Multiply with the H BW
        // Why? So that reweighting is not affected by the reweighting out of the BW
        {
          float& samplePoleMass = mass_map[tree];
          float samplePoleWidth = 0;
          if (samplePoleMass>0.){
            samplePoleWidth = mela.getHiggsWidthAtPoleMass(samplePoleMass);
            rewgt /= pow(pow(GenHMass, 2)-pow(samplePoleMass, 2), 2) + pow(samplePoleMass*samplePoleWidth, 2);
          }
        }
        wgt = fabs(genxsec*genBR*genHEPMCweight*nGenMap[tree].first*rewgt);

        TLorentzVector jet[2], higgs;
        for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(JetPt->at(ij), JetEta->at(ij), JetPhi->at(ij), JetMass->at(ij));
        higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
        TVector3 boostH = higgs.BoostVector();

        SimpleParticleCollection_t associated;
        associated.push_back(SimpleParticle_t(0, jet[0]));
        associated.push_back(SimpleParticle_t(0, jet[1]));

        SimpleParticleCollection_t daughters;
        for (int id=0; id<4; id++){
          double mass=0;
          if (abs(LepLepId->at(id))==13) mass = 0.105658;
          else if (abs(LepLepId->at(id))==11) mass = 0.000511;
          TLorentzVector pDaughter;
          pDaughter.SetPtEtaPhiM(LepPt->at(id), LepEta->at(id), LepPhi->at(id), mass);
          daughters.push_back(SimpleParticle_t(LepLepId->at(id), pDaughter));
        }
        mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

        mela.setProcess(proc, me, prod);

        double propagatorV=1;
        if (recalculate || ev_acc==0){
          TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
          mela.computeProdP(mesq_conserveDifermMass, false);
          if (!recalculate){
            mesq_calc /= cconst_calc;
            mesq_conserveDifermMass=mesq_calc;
          }
          mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
          mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
          if (prod==TVar::Had_ZH){
            mv = mela.getPrimaryMass(23);
            gav = mela.getPrimaryWidth(23);
          }
          else if (prod==TVar::Had_WH){
            mv = mela.getPrimaryMass(24);
            gav = mela.getPrimaryWidth(24);
          }
        }
        else{
          mesq_calc /= cconst_calc;
          mesq_conserveDifermMass=mesq_calc;
        }
        if (fabs(aL1)>0. || fabs(aR1)>0.) mesq_conserveDifermMass /= pow(aL1, 2)+pow(aR1, 2);
        if (fabs(aL2)>0. || fabs(aR2)>0.) mesq_conserveDifermMass /= pow(aL2, 2)+pow(aR2, 2);
        if (prod==TVar::Had_ZH || prod==TVar::Had_WH){
          double mjj = (jet[0] + jet[1]).M();
          propagatorV = 1./(pow(pow(mjj, 2)-pow(mv, 2), 2) + pow(mv*gav, 2));
          mesq_conserveDifermMass /= propagatorV;
          if (recomputePmavjj){
            mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, ((prod==TVar::Had_ZH_S || prod==TVar::Had_ZH) ? TVar::Had_ZH : TVar::Had_WH));
            mela.computeDijetConvBW(pmavjj, false);
          }
          mesq_conserveDifermMass *= pmavjj;
        }

        double mz, gaz, propagator;
        mz = mela.getPrimaryMass(23);
        gaz = mela.getPrimaryWidth(23);
        if (fabs(mzz-mz)<=4.*gaz){
          double sh = pow(mzz, 2);
          double shdn = pow(mz-4.*gaz, 2);
          double shup = pow(mz+4.*gaz, 2);
          double prop_sh = 1./(pow(sh-pow(mz, 2), 2) + pow(mz*gaz, 2));
          double prop_shdn = 1./(pow(shdn-pow(mz, 2), 2) + pow(mz*gaz, 2));
          double prop_shup = 1./(pow(shup-pow(mz, 2), 2) + pow(mz*gaz, 2));
          double fsh = (sh-shdn)/(shup-shdn);
          propagator = prop_sh / (prop_shdn*(1.-fsh) + prop_shup*fsh);
        }
        else propagator=1.;
        mesq_conserveDifermMass /= propagator;

        mela.resetInputEvent();
        if (doRewgt){
          // Reweighting of NLO is too complex to implement here
          exit(1);
        }

        bool doFill = !(
          isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
          );

        if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0, wgt);

        mela.resetInputEvent();
        ev_acc++;
      }

      binList.at(bin).sift(); binList.at(bin).adjustWeights(10);

      for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
        mzz = binList.at(bin).masses.at(ev);
        mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
        wgt = binList.at(bin).weights.at(ev);
        hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
        hvar->Fill(mzz, mzz);
        if (writeFinalTree) newtree->Fill();
      }
    }
    cout << "Total events accumulated: " << ev_acc << endl;

    double* xexyey[4];
    for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[0][bin] = hvar->GetBinContent(bin+1);
      xexyey[1][bin] = hvar->GetBinError(bin+1);

      cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;

      xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
      xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
      xexyey[2][bin] = log10(xexyey[2][bin]);
    }

    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
    foutput->WriteTObject(tg);
    delete tg;

    for (int ix=0; ix<4; ix++) delete[] xexyey[ix];

    foutput->WriteTObject(hmesq_conserveDifermMass);
    foutput->WriteTObject(hvar);
    if (writeFinalTree) foutput->WriteTObject(newtree);
    if (writeFinalTree) delete newtree;
    delete hmesq_conserveDifermMass;
    delete hvar;
    delete[] binning;
  }
  foutput->Close();
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}

/*
SPECIFIC COMMENT: OUTPUT ME DIVIDED BY
- ALPHAS(MZ)**2 TO REMAIN INDEPENDENT OF PDF CHOICE TO FIRST APPROXIMATION
- (aL1**2+aR1**2)*(aL2**2+aR2**2) TO REMAIN INDEPENDENT OF CHANNEL
- PROP(Z)/LINE_PROP(Z) AROUND M4L~=MZ
*/
void get_PAvgProfile_MCFM_JJQCD_bkgZZ_13TeV(int sqrts=13, bool recalculate = true){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";
  bool writeFinalTree=true;

  TVar::Process proc = TVar::bkgZZ;
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod = TVar::JJQCD;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  short NJets30;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float> myJetPt;
  std::vector<float> myJetEta;
  std::vector<float> myJetPhi;
  std::vector<float> myJetMass;

  float mesq_calc=0., cconst_calc=1.;
  float mesq_conserveDifermMass=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  short Z1Flav, Z2Flav;
  float genxsec=0, genBR=0, ngen=0, genHEPMCweight=0;
  float wgt=1;

  TString cinput_main;
  if (sqrts==13) cinput_main = inputdir_13TeV;
  else return;

  TString strchannel[3]={ "4mu", "4e", "2mu2e" };

  vector<TString> strSamples = constructSamplesList("qq_Bkg", sqrts);

  unordered_map<TTree*, pair<float, float>> nGenMap;
  unordered_map<float, vector<TTree*>> mass_sample_map;

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  for (unsigned int is=0; is<strSamples.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          if (!recalculate && tree->GetBranchStatus("pConst_JJQCD_BKG_MCFM_JECNominal")==0) recalculate=true;
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("genxsec", 1); tree->SetBranchAddress("genxsec", &genxsec); tree->SetBranchStatus("genBR", 1); tree->SetBranchAddress("genBR", &genBR);
          tree->SetBranchStatus("genHEPMCweight", 1); tree->SetBranchAddress("genHEPMCweight", &genHEPMCweight);
          tree->SetBranchStatus("nCleanedJetsPt30", 1); tree->SetBranchAddress("nCleanedJetsPt30", &NJets30);
          tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
          tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
          tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
          tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("ZZPt", 1); tree->SetBranchAddress("ZZPt", &ZZPt);
          tree->SetBranchStatus("ZZEta", 1); tree->SetBranchAddress("ZZEta", &ZZEta);
          tree->SetBranchStatus("ZZPhi", 1); tree->SetBranchAddress("ZZPhi", &ZZPhi);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("LepLepId", 1); tree->SetBranchAddress("LepLepId", &LepLepId);
          tree->SetBranchStatus("LepPt", 1); tree->SetBranchAddress("LepPt", &LepPt);
          tree->SetBranchStatus("LepEta", 1); tree->SetBranchAddress("LepEta", &LepEta);
          tree->SetBranchStatus("LepPhi", 1); tree->SetBranchAddress("LepPhi", &LepPhi);
          tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
          tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
          if (!recalculate){
            tree->SetBranchStatus("p_JJQCD_BKG_MCFM_JECNominal", 1);
            tree->SetBranchAddress("p_JJQCD_BKG_MCFM_JECNominal", &mesq_calc);
            tree->SetBranchStatus("pConst_JJQCD_BKG_MCFM_JECNominal", 1);
            tree->SetBranchAddress("pConst_JJQCD_BKG_MCFM_JECNominal", &cconst_calc);
          }
          else cout << "ME needs to be recalculated!" << endl;
          tree->GetEntry(0);
          cout << "Cross section = " << genxsec*genBR << endl;

          TH1F* htmp = (TH1F*)finput->Get(COUNTERS_NAME);
          pair<float, float> nsum(htmp->GetBinContent(1), htmp->GetBinContent(41)); // No PU reweighting
          nGenMap[tree]=nsum;

          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);

          float polemass = findPoleMass(strSamples[is]);
          cout << "Pole mass = " << polemass << endl;
          if (mass_sample_map.find(polemass)==mass_sample_map.end()){
            cout << "Inserting new pole mass sample array" << endl;
            vector<TTree*> dumarr;
            mass_sample_map[polemass] = dumarr;
          }
          mass_sample_map[polemass].push_back(tree);
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }
  cout << "NEntries = " << nEntries << endl;

  for (auto it = mass_sample_map.begin(); it != mass_sample_map.end(); ++it){
    float sum_ngen=0;
    float sum_xsec=0;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      it->second.at(ix)->GetEntry(0);
      sum_ngen += nGenMap[it->second.at(ix)].first;
      sum_xsec += genxsec*genBR;
    }
    float overallWeight = sum_ngen/sum_xsec;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      cout << "Sum Hep MC weights in tree " << ix << " / " << it->second.size() << " was " << nGenMap[it->second.at(ix)].second << " over " << nGenMap[it->second.at(ix)].first << " total gen. events." << endl;
      nGenMap[it->second.at(ix)].first = overallWeight/nGenMap[it->second.at(ix)].second;
      cout << "Event scale for tree " << ix << " / " << it->second.size() << " at pole mass " << it->first << " = " << nGenMap[it->second.at(ix)].first << endl;
    }
  }

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  vector<pair<float, int>> index[3];
  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=
      NJets30>=2 && (
      Z1Flav*Z2Flav==pow(13, 4)
      ||
      Z1Flav*Z2Flav==pow(11, 4)
      ||
      Z1Flav*Z2Flav==pow(11*13, 2)
      )
      ;
    if (!doProcess) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    unsigned int ic = (Z1Flav*Z2Flav==pow(13, 4))*0 + (Z1Flav*Z2Flav==pow(11, 4))*1 + (Z1Flav*Z2Flav==pow(11*13, 2))*2;
    addByLowest(index[ic], mzz, ev);
    if (ic<2) addByLowest(index[1-ic], mzz, ev);
    ev_acc++;
  }

  for (unsigned int ic=0; ic<3; ic++){
    float firstVal=index[ic].at(0).first;
    float lastVal=index[ic].at(index[ic].size()-1).first;
    float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

    float divisor=10000;
    int nbins = index[ic].size()/divisor;
    const int nbins_th=8/*50*/;
    while (nbins<nbins_th){
      if (divisor>1000) divisor -= 1000;
      else if (divisor>100) divisor -= 100;
      else break;
      nbins=index[ic].size()/divisor;
    }
    cout << "nbins=" << nbins << endl;
    if (nbins<3) cerr << "Not enough bins!" << endl;
    vector<ExtBin> binList;
    float* binning = new float[nbins+1];
    binning[0]=infimum;
    binning[nbins]=supremum;
    int ev_stepsize = index[ic].size()/nbins;
    cout << "Event step size: " << ev_stepsize << endl;
    cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
    for (int ix=1; ix<nbins; ix++){
      binning[ix]=(index[ic][ix*ev_stepsize-1].first+index[ic][ix*ev_stepsize].first)*0.5;
      ExtBin tmpbin;
      tmpbin.binlow = binning[ix-1];
      tmpbin.binhigh = binning[ix];
      for (int bin=0; bin<ev_stepsize; bin++) tmpbin.events.push_back(index[ic][(ix-1)*ev_stepsize+bin].second);
      binList.push_back(tmpbin);
      cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ic][ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
    }
    ExtBin tmpbin;
    tmpbin.binlow = binning[nbins-1];
    tmpbin.binhigh = binning[nbins];
    for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index[ic].size(); bin++) tmpbin.events.push_back(index[ic][bin].second);
    binList.push_back(tmpbin);
    cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
    cout << "Bin list has the following bins:" << endl;
    for (unsigned int ib=0; ib<binList.size(); ib++){
      cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    }

    TProfile* hvar = new TProfile(Form("candMass_%s", strchannel[ic].Data()), "", nbins, binning); hvar->Sumw2();
    TProfile* hmesq_conserveDifermMass = new TProfile(Form("P_ConserveDifermionMass_%s", strchannel[ic].Data()), "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

    TTree* newtree=0;
    if (writeFinalTree){
      newtree = new TTree(Form("FinalTree_%s", strchannel[ic].Data()), "");
      newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
      newtree->Branch("ZZMass", &mzz);
      newtree->Branch("weight", &wgt);
    }

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

    double aL1=0, aR1=0, aL2=0, aR2=0;
    double alphasVal=0;
    unsigned int ev_acc=0;
    for (unsigned int bin=0; bin<binList.size(); bin++){
      cout << "Bin " << bin << " is now being scrutinized..." << endl;
      for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
        int getEv = binList.at(bin).events.at(ev);
        getEntry(treeList, getEv);
        if (ev%1000==0) cout << "Doing event " << getEv << endl;
        if (JetPt->size()<2 || JetEta->size()<2 || JetPhi->size()<2 || JetMass->size()<2){
          cerr << "Jet array sizes are less than Njets!" << endl;
          continue;
        }
        TTree* tree = findTree(treeList, getEv);
        wgt = fabs(genxsec*genBR*genHEPMCweight*nGenMap[tree].first);

        TLorentzVector jet[2], higgs;
        for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(JetPt->at(ij), JetEta->at(ij), JetPhi->at(ij), JetMass->at(ij));
        higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
        TVector3 boostH = higgs.BoostVector();

        SimpleParticleCollection_t associated;
        associated.push_back(SimpleParticle_t(0, jet[0]));
        associated.push_back(SimpleParticle_t(0, jet[1]));

        SimpleParticleCollection_t daughters;
        for (int id=0; id<4; id++){
          double mass=0;
          if (abs(LepLepId->at(id))==13) mass = 0.105658;
          else if (abs(LepLepId->at(id))==11) mass = 0.000511;
          TLorentzVector pDaughter;
          pDaughter.SetPtEtaPhiM(LepPt->at(id), LepEta->at(id), LepPhi->at(id), mass);
          daughters.push_back(SimpleParticle_t(LepLepId->at(id), pDaughter));
        }
        mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

        mela.setProcess(proc, me, prod);

        if (recalculate || ev_acc==0){
          TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
          mela.computeProdP(mesq_conserveDifermMass, false);
          if (!recalculate){
            mesq_calc /= cconst_calc;
            mesq_conserveDifermMass=mesq_calc;
          }
          mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
          mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
          alphasVal = mela.getIORecord()->getAlphaSatMZ();
          if (ev_acc==0 && !recalculate){
            cout << "\taL1, aR1, aL2, aR2 = " << aL1 << " , " << aR1 << " , " << aL2 << " , " << aR2 << endl;
            cout << "\talphaS(mZ) = " << alphasVal << endl;
          }
        }
        else{
          mesq_calc /= cconst_calc;
          mesq_conserveDifermMass=mesq_calc;
        }
        if (fabs(aL1)>0. || fabs(aR1)>0.) mesq_conserveDifermMass /= pow(aL1, 2)+pow(aR1, 2);
        if (fabs(aL2)>0. || fabs(aR2)>0.) mesq_conserveDifermMass /= pow(aL2, 2)+pow(aR2, 2);
        mesq_conserveDifermMass /= pow(alphasVal, 2);

        double mz, gaz, propagator;
        mz = mela.getPrimaryMass(23);
        gaz = mela.getPrimaryWidth(23);
        if (fabs(mzz-mz)<=4.*gaz){
          double sh = pow(mzz, 2);
          double shdn = pow(mz-4.*gaz, 2);
          double shup = pow(mz+4.*gaz, 2);
          double prop_sh = 1./(pow(sh-pow(mz, 2), 2) + pow(mz*gaz, 2));
          double prop_shdn = 1./(pow(shdn-pow(mz, 2), 2) + pow(mz*gaz, 2));
          double prop_shup = 1./(pow(shup-pow(mz, 2), 2) + pow(mz*gaz, 2));
          double fsh = (sh-shdn)/(shup-shdn);
          propagator = prop_sh / (prop_shdn*(1.-fsh) + prop_shup*fsh);
        }
        else propagator=1.;
        mesq_conserveDifermMass /= propagator;

        bool doFill = !(
          isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
          );

        if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0, wgt);

        mela.resetInputEvent();
        ev_acc++;
      }

      binList.at(bin).sift(); binList.at(bin).adjustWeights(-1);

      for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
        mzz = binList.at(bin).masses.at(ev);
        mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
        wgt = binList.at(bin).weights.at(ev);
        hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
        hvar->Fill(mzz, mzz);
        if (writeFinalTree) newtree->Fill();
      }
    }

    double* xexyey[4];
    for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[0][bin] = hvar->GetBinContent(bin+1);
      xexyey[1][bin] = hvar->GetBinError(bin+1);

      cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;

      xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
      xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
      xexyey[2][bin] = log10(xexyey[2][bin]);
    }

    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
    foutput->WriteTObject(tg);
    delete tg;

    for (int ix=0; ix<4; ix++) delete[] xexyey[ix];

    foutput->WriteTObject(hmesq_conserveDifermMass);
    foutput->WriteTObject(hvar);
    if (writeFinalTree) foutput->WriteTObject(newtree);
    if (writeFinalTree) delete newtree;
    delete hmesq_conserveDifermMass;
    delete hvar;
    delete[] binning;
  }
  foutput->Close();
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
}


/* SPECIFIC COMMENT: OUTPUT ME DIVIDED BY (aL1**2+aR1**2)*(aL2**2+aR2**2) TO REMAIN INDEPENDENT OF CHANNEL */
void get_PAvgProfile_JHUGen_ZZINDEPENDENT_HSMHiggs(bool recalculate=false){
  int erg_tev=8;
  float mPOLE=125.;
  TString TREE_NAME;
  const bool writeFinalTree=false;
  const TString strMEBranchname = "p_ZZINDEPENDENT_SIG_ghz1_1_JHUGen";
  const TString strConstBranchname = "pConst_ZZINDEPENDENT_SIG_ghz1_1_JHUGen";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::ZZINDEPENDENT;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  float mesq_calc=0., cconst_calc=1.;
  float mesq_conserveDifermMass=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  short Z1Flav, Z2Flav;
  int LepID[4];

  vector<pair<float, int>> index[3];
  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s.root", strme.Data(), strprod.Data(), strproc.Data()), "recreate");

  vector<TString> dumappend;
  vector<TString> strSamples_13TeV = constructSamplesList("JJQCD", 13.);
  dumappend = constructSamplesList("JJVBF", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);

  vector<TString> strSamples_8TeV = constructSamplesList("JJQCD", 8.);
  dumappend = constructSamplesList("JJVBF", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);

  vector<TString> strSamples_7TeV = constructSamplesList("JJQCD", 7.);
  dumappend = constructSamplesList("JJVBF", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  TString cinput_main;

  TREE_NAME = "ZZTree/candTree";
  cinput_main = inputdir_13TeV;
  //for (int is=0; is<2; is++){
  for (int is=0; is<(int)strSamples_13TeV.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples_13TeV[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;

          bool doRecalculate = recalculate;
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0) doRecalculate = true;
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
          tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0){
            tree->SetBranchStatus(strMEBranchname, 1); tree->SetBranchAddress(strMEBranchname, &mesq_calc);
            tree->SetBranchStatus(strConstBranchname, 1); tree->SetBranchAddress(strConstBranchname, &cconst_calc);
          }
          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }

  for (int ic=0; ic<3; ic++){
    TREE_NAME = "SelectedTree";
    cinput_main = inputdir_8TeV;
    //for (int is=0; is<2; is++){
    for (int is=0; is<(int)strSamples_8TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_8TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
    cinput_main = inputdir_7TeV;
    //for (int is=0; is<2; is++){
    for (int is=0; is<(int)strSamples_7TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_7TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }

  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;

  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=
      (
      Z1Flav*Z2Flav==pow(13, 4)
      ||
      Z1Flav*Z2Flav==pow(11, 4)
      ||
      Z1Flav*Z2Flav==pow(11*13, 2)
      )
      ;
    if (!doProcess) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    unsigned int ic = (Z1Flav*Z2Flav==pow(13, 4))*0 + (Z1Flav*Z2Flav==pow(11, 4))*1 + (Z1Flav*Z2Flav==pow(11*13, 2))*2;
    addByLowest(index[ic], mzz, ev);
    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;

  for (unsigned int ic=0; ic<3; ic++){
    float firstVal=index[ic].at(0).first;
    float lastVal=index[ic].at(index[ic].size()-1).first;
    float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

    float divisor=95000;
    int nbins = index[ic].size()/divisor;
    const int nbins_th=10/*50*/;
    while (nbins<nbins_th){
      if (divisor>1000) divisor -= 1000;
      else if (divisor>100) divisor -= 100;
      else break;
      nbins=index[ic].size()/divisor;
    }
    cout << "nbins=" << nbins << endl;
    if (nbins<3) cerr << "Not enough bins!" << endl;
    vector<ExtBin> binList;
    float* binning = new float[nbins+1];
    binning[0]=infimum;
    binning[nbins]=supremum;
    int ev_stepsize = index[ic].size()/nbins;
    cout << "Event step size: " << ev_stepsize << endl;
    cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
    for (int ix=1; ix<nbins; ix++){
      binning[ix]=(index[ic][ix*ev_stepsize-1].first+index[ic][ix*ev_stepsize].first)*0.5;
      ExtBin tmpbin;
      tmpbin.binlow = binning[ix-1]; tmpbin.binhigh = binning[ix];
      for (int bin=0; bin<ev_stepsize; bin++){
        int evid = index[ic][(ix-1)*ev_stepsize+bin].second;
        tmpbin.events.push_back(evid);
      }
      binList.push_back(tmpbin);
      cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ic][ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
    }
    ExtBin tmpbin;
    tmpbin.binlow = binning[nbins-1]; tmpbin.binhigh = binning[nbins];
    for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index[ic].size(); bin++){
      int evid = index[ic][bin].second;
      tmpbin.events.push_back(evid);
    }
    binList.push_back(tmpbin);
    cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
    cout << "Bin list has the following bins:" << endl;
    for (unsigned int ib=0; ib<binList.size(); ib++){
      cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    }

    foutput->cd();

    TProfile* hvar = new TProfile(Form("candMass_%s", strchannel[ic].Data()), "", nbins, binning); hvar->Sumw2();
    TProfile* hmesq_conserveDifermMass = new TProfile(Form("P_ConserveDifermionMass_%s", strchannel[ic].Data()), "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

    TTree* newtree=0;
    if (writeFinalTree){
      newtree = new TTree(Form("FinalTree_%s", strchannel[ic].Data()), "");
      newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
      newtree->Branch("ZZMass", &mzz);
    }

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

    if (ic==1){
      LepID[0]=11;
      LepID[1]=-11;
    }
    else{
      LepID[0]=13;
      LepID[1]=-13;
    }
    if (ic==0){
      LepID[2]=13;
      LepID[3]=-13;
    }
    else{
      LepID[2]=11;
      LepID[3]=-11;
    }

    double aL1=0, aR1=0, aL2=0, aR2=0;
    unsigned int ctr=0;
    for (unsigned int bin=0; bin<binList.size(); bin++){
      cout << "Bin " << bin << " is now being scrutinized..." << endl;
      for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
        int getEv = binList.at(bin).events.at(ev);
        getEntry(treeList, getEv);
        TTree* tree = findTree(treeList, getEv);

        if (ev%1000==0) cout << "Doing event " << getEv << endl;

        TLorentzVector pDaughters[4];
        std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
        for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); }
        SimpleParticleCollection_t daughters;
        for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
        mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

        mela.setProcess(proc, me, prod);
        TUtil::setLeptonMassScheme(TVar::ConserveDifermionMass);

        bool hasConst = (tree->GetBranchStatus(strConstBranchname)==1);
        bool doCalc=(!hasConst || ctr==0);
        if (doCalc){
          mela.computeP(mesq_conserveDifermMass, false);
          if (hasConst) mesq_conserveDifermMass = mesq_calc / cconst_calc;

          mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
          mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
        }
        else mesq_conserveDifermMass = mesq_calc / cconst_calc;

        if (fabs(aL1)>0. || fabs(aR1)>0.) mesq_conserveDifermMass /= pow(aL1, 2)+pow(aR1, 2);
        if (fabs(aL2)>0. || fabs(aR2)>0.) mesq_conserveDifermMass /= pow(aL2, 2)+pow(aR2, 2);

        bool doFill = !(
          isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
          );

        if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0);

        mela.resetInputEvent();
        ctr++;
      }

      binList.at(bin).sift(); binList.at(bin).adjustWeights();

      for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
        mzz = binList.at(bin).masses.at(ev);
        mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
        hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
        hvar->Fill(mzz, mzz);
        if (writeFinalTree) newtree->Fill();
      }
    }

    double* xexyey[4];
    for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[0][bin] = hvar->GetBinContent(bin+1);
      xexyey[1][bin] = hvar->GetBinError(bin+1);

      cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;
      xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
      xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
      xexyey[2][bin] = log10(xexyey[2][bin]);
    }

    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
    foutput->WriteTObject(tg);
    delete tg;

    for (int ix=0; ix<4; ix++) delete[] xexyey[ix];
    foutput->WriteTObject(hmesq_conserveDifermMass);
    foutput->WriteTObject(hvar);
    if (writeFinalTree) foutput->WriteTObject(newtree);
    if (writeFinalTree) delete newtree;
    delete hmesq_conserveDifermMass;
    delete hvar;
    delete[] binning;
  }
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
  foutput->Close();
}

/*
SPECIFIC COMMENT: OUTPUT ME DIVIDED BY
- ALPHAS(MZ)**2 TO REMAIN INDEPENDENT OF PDF CHOICE TO FIRST APPROXIMATION
- H(1) PROPAGATOR
- (aL1**2+aR1**2)*(aL2**2+aR2**2) TO REMAIN INDEPENDENT OF CHANNEL
*/
void get_PAvgProfile_JHUGen_ZZGG_HSMHiggs(bool recalculate=false){
  int erg_tev=8;
  float mPOLE=125.;
  TString TREE_NAME;
  const bool writeFinalTree=false;
  const TString strMEBranchname = "p_GG_SIG_ghg2_1_ghz1_1_JHUGen";
  const TString strConstBranchname = "pConst_GG_SIG_ghg2_1_ghz1_1_JHUGen";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::ZZGG;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  float mesq_calc=0., cconst_calc=1.;
  float mesq_conserveDifermMass=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  short Z1Flav, Z2Flav;
  int LepID[4];

  vector<pair<float, int>> index[3];
  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s.root", strme.Data(), strprod.Data(), strproc.Data()), "recreate");

  vector<TString> dumappend;
  vector<TString> strSamples_13TeV = constructSamplesList("JJQCD", 13.);
  dumappend = constructSamplesList("JJVBF", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);

  vector<TString> strSamples_8TeV = constructSamplesList("JJQCD", 8.);
  dumappend = constructSamplesList("JJVBF", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);

  vector<TString> strSamples_7TeV = constructSamplesList("JJQCD", 7.);
  dumappend = constructSamplesList("JJVBF", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  TString cinput_main;

  TREE_NAME = "ZZTree/candTree";
  cinput_main = inputdir_13TeV;
  //for (int is=0; is<2; is++){
  for (int is=0; is<(int)strSamples_13TeV.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples_13TeV[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;

          bool doRecalculate = recalculate;
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0) doRecalculate = true;
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
          tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0){
            tree->SetBranchStatus(strMEBranchname, 1); tree->SetBranchAddress(strMEBranchname, &mesq_calc);
            tree->SetBranchStatus(strConstBranchname, 1); tree->SetBranchAddress(strConstBranchname, &cconst_calc);
          }
          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }

  for (int ic=0; ic<3; ic++){
    TREE_NAME = "SelectedTree";
    cinput_main = inputdir_8TeV;
    //for (int is=0; is<2; is++){
    for (int is=0; is<(int)strSamples_8TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_8TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
    cinput_main = inputdir_7TeV;
    //for (int is=0; is<2; is++){
    for (int is=0; is<(int)strSamples_7TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_7TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }

  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;

  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=
      (
      Z1Flav*Z2Flav==pow(13, 4)
      ||
      Z1Flav*Z2Flav==pow(11, 4)
      ||
      Z1Flav*Z2Flav==pow(11*13, 2)
      )
      ;
    if (!doProcess) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    unsigned int ic = (Z1Flav*Z2Flav==pow(13, 4))*0 + (Z1Flav*Z2Flav==pow(11, 4))*1 + (Z1Flav*Z2Flav==pow(11*13, 2))*2;
    addByLowest(index[ic], mzz, ev);
    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;

  for (unsigned int ic=0; ic<3; ic++){
    float firstVal=index[ic].at(0).first;
    float lastVal=index[ic].at(index[ic].size()-1).first;
    float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

    float divisor=95000;
    int nbins = index[ic].size()/divisor;
    const int nbins_th=10/*50*/;
    while (nbins<nbins_th){
      if (divisor>1000) divisor -= 1000;
      else if (divisor>100) divisor -= 100;
      else break;
      nbins=index[ic].size()/divisor;
    }
    cout << "nbins=" << nbins << endl;
    if (nbins<3) cerr << "Not enough bins!" << endl;
    vector<ExtBin> binList;
    float* binning = new float[nbins+1];
    binning[0]=infimum;
    binning[nbins]=supremum;
    int ev_stepsize = index[ic].size()/nbins;
    cout << "Event step size: " << ev_stepsize << endl;
    cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
    for (int ix=1; ix<nbins; ix++){
      binning[ix]=(index[ic][ix*ev_stepsize-1].first+index[ic][ix*ev_stepsize].first)*0.5;
      ExtBin tmpbin;
      tmpbin.binlow = binning[ix-1]; tmpbin.binhigh = binning[ix];
      for (int bin=0; bin<ev_stepsize; bin++){
        int evid = index[ic][(ix-1)*ev_stepsize+bin].second;
        tmpbin.events.push_back(evid);
      }
      binList.push_back(tmpbin);
      cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ic][ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
    }
    ExtBin tmpbin;
    tmpbin.binlow = binning[nbins-1]; tmpbin.binhigh = binning[nbins];
    for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index[ic].size(); bin++){
      int evid = index[ic][bin].second;
      tmpbin.events.push_back(evid);
    }
    binList.push_back(tmpbin);
    cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
    cout << "Bin list has the following bins:" << endl;
    for (unsigned int ib=0; ib<binList.size(); ib++){
      cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    }

    foutput->cd();

    TProfile* hvar = new TProfile(Form("candMass_%s", strchannel[ic].Data()), "", nbins, binning); hvar->Sumw2();
    TProfile* hmesq_conserveDifermMass = new TProfile(Form("P_ConserveDifermionMass_%s", strchannel[ic].Data()), "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

    TTree* newtree=0;
    if (writeFinalTree){
      newtree = new TTree(Form("FinalTree_%s", strchannel[ic].Data()), "");
      newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
      newtree->Branch("ZZMass", &mzz);
    }

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

    if (ic==1){
      LepID[0]=11;
      LepID[1]=-11;
    }
    else{
      LepID[0]=13;
      LepID[1]=-13;
    }
    if (ic==0){
      LepID[2]=13;
      LepID[3]=-13;
    }
    else{
      LepID[2]=11;
      LepID[3]=-11;
    }

    double alphasVal=0, mh=0, gah=0, aL1=0, aR1=0, aL2=0, aR2=0;
    unsigned int ctr=0;
    for (unsigned int bin=0; bin<binList.size(); bin++){
      cout << "Bin " << bin << " is now being scrutinized..." << endl;
      for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
        int getEv = binList.at(bin).events.at(ev);
        getEntry(treeList, getEv);
        TTree* tree = findTree(treeList, getEv);

        if (ev%1000==0) cout << "Doing event " << getEv << endl;

        TLorentzVector pDaughters[4];
        std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
        for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); }
        SimpleParticleCollection_t daughters;
        for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
        mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

        mela.setProcess(proc, me, prod);
        TUtil::setLeptonMassScheme(TVar::ConserveDifermionMass);

        bool hasConst = (tree->GetBranchStatus(strConstBranchname)==1);
        bool doCalc=(!hasConst || ctr==0);
        if (doCalc){
          mela.computeP(mesq_conserveDifermMass, false);
          if (hasConst) mesq_conserveDifermMass = mesq_calc / cconst_calc;

          alphasVal = mela.getIORecord()->getAlphaSatMZ();
          mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
          mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
          mela.getIORecord()->getHiggsMassWidth(mh, gah, 0);
        }
        else{
          mesq_conserveDifermMass = mesq_calc / cconst_calc;
          mh=mzz;
          gah = mela.getHiggsWidthAtPoleMass(mh);
        }

        double propagator = 1./(pow(pow(mzz, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
        mesq_conserveDifermMass /= propagator;
        mesq_conserveDifermMass /= pow(alphasVal, 2);
        if (fabs(aL1)>0. || fabs(aR1)>0.) mesq_conserveDifermMass /= pow(aL1, 2)+pow(aR1, 2);
        if (fabs(aL2)>0. || fabs(aR2)>0.) mesq_conserveDifermMass /= pow(aL2, 2)+pow(aR2, 2);

        bool doFill = !(
          isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
          );

        if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0);

        mela.resetInputEvent();
        ctr++;
      }

      binList.at(bin).sift(); binList.at(bin).adjustWeights();

      for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
        mzz = binList.at(bin).masses.at(ev);
        mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
        hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
        hvar->Fill(mzz, mzz);
        if (writeFinalTree) newtree->Fill();
      }
    }

    double* xexyey[4];
    for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[0][bin] = hvar->GetBinContent(bin+1);
      xexyey[1][bin] = hvar->GetBinError(bin+1);

      cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;
      xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
      xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
      xexyey[2][bin] = log10(xexyey[2][bin]);
    }


    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
    foutput->WriteTObject(tg);
    delete tg;

    for (int ix=0; ix<4; ix++) delete[] xexyey[ix];
    foutput->WriteTObject(hmesq_conserveDifermMass);
    foutput->WriteTObject(hvar);
    if (writeFinalTree) foutput->WriteTObject(newtree);
    if (writeFinalTree) delete newtree;
    delete hmesq_conserveDifermMass;
    delete hvar;
    delete[] binning;
  }
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
  foutput->Close();
}

/*
SPECIFIC COMMENT: OUTPUT ME DIVIDED BY
- ALPHAS(MZ)**2 TO REMAIN INDEPENDENT OF PDF CHOICE TO FIRST APPROXIMATION
- H(1) PROPAGATOR
- (aL1**2+aR1**2)*(aL2**2+aR2**2) TO REMAIN INDEPENDENT OF CHANNEL
*/
void get_PAvgProfile_MCFM_ZZGG_HSMHiggs(bool recalculate=false){
  int erg_tev=8;
  float mPOLE=125.;
  TString TREE_NAME;
  const bool writeFinalTree=false;
  const TString strMEBranchname = "p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM";
  const TString strConstBranchname = "pConst_GG_SIG_kappaTopBot_1_ghz1_1_MCFM";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod = TVar::ZZGG;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  float mesq_calc=0., cconst_calc=1.;
  float mesq_conserveDifermMass=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  short Z1Flav, Z2Flav;
  int LepID[4];

  vector<pair<float, int>> index[3];
  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s.root", strme.Data(), strprod.Data(), strproc.Data()), "recreate");

  vector<TString> dumappend;
  vector<TString> strSamples_13TeV = constructSamplesList("JJQCD", 13.);
  dumappend = constructSamplesList("JJVBF", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);

  vector<TString> strSamples_8TeV = constructSamplesList("JJQCD", 8.);
  dumappend = constructSamplesList("JJVBF", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);

  vector<TString> strSamples_7TeV = constructSamplesList("JJQCD", 7.);
  dumappend = constructSamplesList("JJVBF", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  TString cinput_main;

  TREE_NAME = "ZZTree/candTree";
  cinput_main = inputdir_13TeV;
  //for (int is=0; is<2; is++){
  for (int is=0; is<(int)strSamples_13TeV.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples_13TeV[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;

          bool doRecalculate = recalculate;
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0) doRecalculate = true;
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
          tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0){
            tree->SetBranchStatus(strMEBranchname, 1); tree->SetBranchAddress(strMEBranchname, &mesq_calc);
            tree->SetBranchStatus(strConstBranchname, 1); tree->SetBranchAddress(strConstBranchname, &cconst_calc);
          }
          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }

  for (int ic=0; ic<3; ic++){
    TREE_NAME = "SelectedTree";
    cinput_main = inputdir_8TeV;
    //for (int is=0; is<2; is++){
    for (int is=0; is<(int)strSamples_8TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_8TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
    cinput_main = inputdir_7TeV;
    //for (int is=0; is<2; is++){
    for (int is=0; is<(int)strSamples_7TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_7TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }

  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;

  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=
      (
      Z1Flav*Z2Flav==pow(13, 4)
      ||
      Z1Flav*Z2Flav==pow(11, 4)
      ||
      Z1Flav*Z2Flav==pow(11*13, 2)
      )
      ;
    if (!doProcess) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    unsigned int ic = (Z1Flav*Z2Flav==pow(13, 4))*0 + (Z1Flav*Z2Flav==pow(11, 4))*1 + (Z1Flav*Z2Flav==pow(11*13, 2))*2;
    addByLowest(index[ic], mzz, ev);
    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;

  for (unsigned int ic=0; ic<3; ic++){
    float firstVal=index[ic].at(0).first;
    float lastVal=index[ic].at(index[ic].size()-1).first;
    float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

    float divisor=95000;
    int nbins = index[ic].size()/divisor;
    const int nbins_th=10/*50*/;
    while (nbins<nbins_th){
      if (divisor>1000) divisor -= 1000;
      else if (divisor>100) divisor -= 100;
      else break;
      nbins=index[ic].size()/divisor;
    }
    cout << "nbins=" << nbins << endl;
    if (nbins<3) cerr << "Not enough bins!" << endl;
    vector<ExtBin> binList;
    float* binning = new float[nbins+1];
    binning[0]=infimum;
    binning[nbins]=supremum;
    int ev_stepsize = index[ic].size()/nbins;
    cout << "Event step size: " << ev_stepsize << endl;
    cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
    for (int ix=1; ix<nbins; ix++){
      binning[ix]=(index[ic][ix*ev_stepsize-1].first+index[ic][ix*ev_stepsize].first)*0.5;
      ExtBin tmpbin;
      tmpbin.binlow = binning[ix-1]; tmpbin.binhigh = binning[ix];
      for (int bin=0; bin<ev_stepsize; bin++){
        int evid = index[ic][(ix-1)*ev_stepsize+bin].second;
        tmpbin.events.push_back(evid);
      }
      binList.push_back(tmpbin);
      cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ic][ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
    }
    ExtBin tmpbin;
    tmpbin.binlow = binning[nbins-1]; tmpbin.binhigh = binning[nbins];
    for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index[ic].size(); bin++){
      int evid = index[ic][bin].second;
      tmpbin.events.push_back(evid);
    }
    binList.push_back(tmpbin);
    cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
    cout << "Bin list has the following bins:" << endl;
    for (unsigned int ib=0; ib<binList.size(); ib++){
      cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    }

    foutput->cd();

    TProfile* hvar = new TProfile(Form("candMass_%s", strchannel[ic].Data()), "", nbins, binning); hvar->Sumw2();
    TProfile* hmesq_conserveDifermMass = new TProfile(Form("P_ConserveDifermionMass_%s", strchannel[ic].Data()), "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

    TTree* newtree=0;
    if (writeFinalTree){
      newtree = new TTree(Form("FinalTree_%s", strchannel[ic].Data()), "");
      newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
      newtree->Branch("ZZMass", &mzz);
    }

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

    if (ic==1){
      LepID[0]=11;
      LepID[1]=-11;
    }
    else{
      LepID[0]=13;
      LepID[1]=-13;
    }
    if (ic==0){
      LepID[2]=13;
      LepID[3]=-13;
    }
    else{
      LepID[2]=11;
      LepID[3]=-11;
    }

    double alphasVal=0, mh=0, gah=0, aL1=0, aR1=0, aL2=0, aR2=0;
    unsigned int ctr=0;
    for (unsigned int bin=0; bin<binList.size(); bin++){
      cout << "Bin " << bin << " is now being scrutinized..." << endl;
      for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
        int getEv = binList.at(bin).events.at(ev);
        getEntry(treeList, getEv);
        TTree* tree = findTree(treeList, getEv);

        if (ev%1000==0) cout << "Doing event " << getEv << endl;

        TLorentzVector pDaughters[4];
        std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
        for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); }
        SimpleParticleCollection_t daughters;
        for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
        mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

        mela.setProcess(proc, me, prod);
        TUtil::setLeptonMassScheme(TVar::ConserveDifermionMass);

        bool hasConst = (tree->GetBranchStatus(strConstBranchname)==1);
        bool doCalc=(!hasConst || ctr==0);
        if (doCalc){
          mela.computeP(mesq_conserveDifermMass, false);
          if (hasConst) mesq_conserveDifermMass = mesq_calc / cconst_calc;

          alphasVal = mela.getIORecord()->getAlphaSatMZ();
          mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
          mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
          mela.getIORecord()->getHiggsMassWidth(mh, gah, 0);
        }
        else mesq_conserveDifermMass = mesq_calc / cconst_calc;

        double propagator = 1./(pow(pow(mzz, 2)-pow(mh, 2), 2) + pow(mh*gah, 2));
        mesq_conserveDifermMass /= propagator;
        mesq_conserveDifermMass /= pow(alphasVal, 2);
        if (fabs(aL1)>0. || fabs(aR1)>0.) mesq_conserveDifermMass /= pow(aL1, 2)+pow(aR1, 2);
        if (fabs(aL2)>0. || fabs(aR2)>0.) mesq_conserveDifermMass /= pow(aL2, 2)+pow(aR2, 2);

        bool doFill = !(
          isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
          );

        if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0);

        mela.resetInputEvent();
        ctr++;
      }

      binList.at(bin).sift(); binList.at(bin).adjustWeights();

      for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
        mzz = binList.at(bin).masses.at(ev);
        mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
        hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
        hvar->Fill(mzz, mzz);
        if (writeFinalTree) newtree->Fill();
      }
    }

    double* xexyey[4];
    for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[0][bin] = hvar->GetBinContent(bin+1);
      xexyey[1][bin] = hvar->GetBinError(bin+1);

      cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;
      xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
      xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
      xexyey[2][bin] = log10(xexyey[2][bin]);
    }


    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
    foutput->WriteTObject(tg);
    delete tg;

    for (int ix=0; ix<4; ix++) delete[] xexyey[ix];
    foutput->WriteTObject(hmesq_conserveDifermMass);
    foutput->WriteTObject(hvar);
    if (writeFinalTree) foutput->WriteTObject(newtree);
    if (writeFinalTree) delete newtree;
    delete hmesq_conserveDifermMass;
    delete hvar;
    delete[] binning;
  }
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
  foutput->Close();
}

/*
SPECIFIC COMMENT: OUTPUT ME DIVIDED BY
- ALPHAS(MZ)**2 TO REMAIN INDEPENDENT OF PDF CHOICE TO FIRST APPROXIMATION
- (aL1**2+aR1**2)*(aL2**2+aR2**2) TO REMAIN INDEPENDENT OF CHANNEL
*/
void get_PAvgProfile_MCFM_ZZGG_bkgZZ(bool recalculate=false){
  int erg_tev=8;
  float mPOLE=125.;
  TString TREE_NAME;
  const bool writeFinalTree=true;
  const TString strMEBranchname = "p_GG_BKG_MCFM";
  const TString strConstBranchname = "pConst_GG_BKG_MCFM";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::bkgZZ;
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod = TVar::ZZGG;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  float mesq_calc=0., cconst_calc=1.;
  float mesq_conserveDifermMass=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  float wgt=1.;
  float GenHMass;
  short Z1Flav, Z2Flav;
  int LepID[4];
  short useWeighted;

  vector<pair<float, int>> index[3];
  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s.root", strme.Data(), strprod.Data(), strproc.Data()), "recreate");

  vector<TString> dumappend;
  vector<TString> strSamples_13TeV = constructSamplesList("gg_Bkg_MCFM", 13.);
  dumappend = constructSamplesList("gg_Bkg_ggVV", 13.);
  appendVector<TString>(strSamples_13TeV, dumappend);

  vector<TString> strSamples_8TeV = constructSamplesList("gg_Bkg_MCFM", 8.);
  dumappend = constructSamplesList("gg_Bkg_ggVV", 8.);
  appendVector<TString>(strSamples_8TeV, dumappend);

  vector<TString> strSamples_7TeV = constructSamplesList("gg_Bkg_MCFM", 7.);
  dumappend = constructSamplesList("gg_Bkg_ggVV", 7.);
  appendVector<TString>(strSamples_7TeV, dumappend);

  // Consider if weights exist
  unordered_map<TTree*, pair<float, float>> samplePoleMasses;

  vector<TString> strSamples_weighted_13TeV = constructSamplesList("JJQCD", 13.);
  dumappend = constructSamplesList("gg_Sig_JHUGen", 13.);
  appendVector<TString>(strSamples_weighted_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_MCFM", 13.);
  appendVector<TString>(strSamples_weighted_13TeV, dumappend);
  dumappend = constructSamplesList("gg_Sig_ggVV", 13.);
  appendVector<TString>(strSamples_weighted_13TeV, dumappend);

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  TString cinput_main;

  TREE_NAME = "ZZTree/candTree";
  cinput_main = inputdir_13TeV;
  //for (int is=0; is<1; is++){
  for (int is=0; is<(int)strSamples_weighted_13TeV.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples_weighted_13TeV[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          if (tree->GetBranchStatus("p_Gen_GG_BKG_MCFM")==1){
            bool doRecalculate = recalculate;
            if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0) doRecalculate = true;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
            tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
            tree->SetBranchStatus("p_Gen_GG_BKG_MCFM", 1); tree->SetBranchAddress("p_Gen_GG_BKG_MCFM", &wgt);
            tree->SetBranchStatus("GenHMass", 1); tree->SetBranchAddress("GenHMass", &GenHMass);
            if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0){
              tree->SetBranchStatus(strMEBranchname, 1); tree->SetBranchAddress(strMEBranchname, &mesq_calc);
              tree->SetBranchStatus(strConstBranchname, 1); tree->SetBranchAddress(strConstBranchname, &cconst_calc);
            }
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);

            float polemass = findPoleMass(strSamples_weighted_13TeV[is]);
            float polewidth = mela.getHiggsWidthAtPoleMass(polemass);
            cout << strSamples_weighted_13TeV[is] << " pole mass = " << polemass << ", pole width = " << polewidth << endl;
            samplePoleMasses[tree] = pair<float, float>(polemass, polewidth);
          }
          else{
            cerr << TREE_NAME << " in " << cinput << " does not have weight p_Gen_GG_BKG_MCFM." << endl;
            finput->Close();
          }
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }

  TREE_NAME = "ZZTree/candTree";
  cinput_main = inputdir_13TeV;
  //for (int is=0; is<1; is++){
  for (int is=0; is<(int)strSamples_13TeV.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples_13TeV[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          bool doRecalculate = recalculate;
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0) doRecalculate = true;
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
          tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0){
            tree->SetBranchStatus(strMEBranchname, 1); tree->SetBranchAddress(strMEBranchname, &mesq_calc);
            tree->SetBranchStatus(strConstBranchname, 1); tree->SetBranchAddress(strConstBranchname, &cconst_calc);
          }
          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }

  TREE_NAME = "SelectedTree";
  for (int ic=0; ic<3; ic++){
    cinput_main = inputdir_8TeV;
    //for (int is=0; is<0; is++){
    for (int is=0; is<(int)strSamples_8TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_8TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
    cinput_main = inputdir_7TeV;
    //for (int is=0; is<0; is++){
    for (int is=0; is<(int)strSamples_7TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_7TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }

  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;

  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=
      (
      Z1Flav*Z2Flav==pow(13, 4)
      ||
      Z1Flav*Z2Flav==pow(11, 4)
      ||
      Z1Flav*Z2Flav==pow(11*13, 2)
      )
      ;
    if (!doProcess) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    unsigned int ic = (Z1Flav*Z2Flav==pow(13, 4))*0 + (Z1Flav*Z2Flav==pow(11, 4))*1 + (Z1Flav*Z2Flav==pow(11*13, 2))*2;
    addByLowest(index[ic], mzz, ev);
    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;

  for (unsigned int ic=0; ic<3; ic++){
    float firstVal=index[ic].at(0).first;
    float lastVal=index[ic].at(index[ic].size()-1).first;
    float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

    float divisor=100000;
    if (ic!=2) divisor=130000;
    int nbins = index[ic].size()/divisor;
    const int nbins_th=10/*50*/;
    while (nbins<nbins_th){
      if (divisor>1000) divisor -= 1000;
      else if (divisor>100) divisor -= 100;
      else break;
      nbins=index[ic].size()/divisor;
    }
    cout << "nbins=" << nbins << endl;
    if (nbins<3) cerr << "Not enough bins!" << endl;
    vector<ExtBin> binList;
    vector<ExtBin> weightedBinList;
    float* binning = new float[nbins+1];
    binning[0]=infimum;
    binning[nbins]=supremum;
    int ev_stepsize = index[ic].size()/nbins;
    cout << "Event step size: " << ev_stepsize << endl;
    cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
    for (int ix=1; ix<nbins; ix++){
      binning[ix]=(index[ic][ix*ev_stepsize-1].first+index[ic][ix*ev_stepsize].first)*0.5;
      ExtBin tmpbin, tmpbin_weighted;
      tmpbin.binlow = binning[ix-1]; tmpbin.binhigh = binning[ix];
      tmpbin_weighted.binlow = binning[ix-1]; tmpbin_weighted.binhigh = binning[ix];
      for (int bin=0; bin<ev_stepsize; bin++){
        int evid = index[ic][(ix-1)*ev_stepsize+bin].second;
        TTree* tree = findTree(treeList, evid);
        if (tree->GetBranchStatus("p_Gen_GG_BKG_MCFM")==1) tmpbin_weighted.events.push_back(evid);
        else tmpbin.events.push_back(evid);
      }
      binList.push_back(tmpbin);
      weightedBinList.push_back(tmpbin_weighted);
      cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ic][ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
    }
    ExtBin tmpbin, tmpbin_weighted;
    tmpbin.binlow = binning[nbins-1]; tmpbin.binhigh = binning[nbins];
    tmpbin_weighted.binlow = binning[nbins-1]; tmpbin_weighted.binhigh = binning[nbins];
    for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index[ic].size(); bin++){
      int evid = index[ic][bin].second;
      TTree* tree = findTree(treeList, evid);
      if (tree->GetBranchStatus("p_Gen_GG_BKG_MCFM")==1) tmpbin_weighted.events.push_back(evid);
      else tmpbin.events.push_back(evid);
    }
    binList.push_back(tmpbin);
    weightedBinList.push_back(tmpbin_weighted);
    cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
    cout << "Bin list has the following bins:" << endl;
    for (unsigned int ib=0; ib<binList.size(); ib++){
      cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    }

    foutput->cd();
    TProfile* hvar = new TProfile(Form("candMass_%s", strchannel[ic].Data()), "", nbins, binning); hvar->Sumw2();
    TProfile* hmesq_conserveDifermMass = new TProfile(Form("P_ConserveDifermionMass_%s", strchannel[ic].Data()), "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

    TTree* newtree=0;
    if (writeFinalTree){
      newtree = new TTree(Form("FinalTree_%s", strchannel[ic].Data()), "");
      newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
      newtree->Branch("ZZMass", &mzz);
      newtree->Branch("weight", &wgt);
      newtree->Branch("isWeighted", &useWeighted);
    }

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

    if (ic==1){
      LepID[0]=11;
      LepID[1]=-11;
    }
    else{
      LepID[0]=13;
      LepID[1]=-13;
    }
    if (ic==0){
      LepID[2]=13;
      LepID[3]=-13;
    }
    else{
      LepID[2]=11;
      LepID[3]=-11;
    }

    double alphasVal=0, aL1=0, aR1=0, aL2=0, aR2=0;
    unsigned int ctr=0;
    for (unsigned int bin=0; bin<binList.size(); bin++){
      cout << "Bin " << bin << " is now being scrutinized..." << endl;
      unsigned int nweighted = weightedBinList.at(bin).events.size();
      unsigned int nunweighted = binList.at(bin).events.size();
      unsigned int ntotal = nweighted + nunweighted;
      for (unsigned int ev = 0; ev < ntotal; ev++){
        int getEv;
        bool isweighted = (ev<nweighted);
        if (isweighted) getEv = weightedBinList.at(bin).events.at(ev);
        else getEv = binList.at(bin).events.at(ev-nweighted);
        wgt=1;
        getEntry(treeList, getEv);
        TTree* tree = findTree(treeList, getEv);

        if (ev%1000==0) cout << "Doing event " << getEv << endl;

        TLorentzVector pDaughters[4];
        std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
        for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); }
        SimpleParticleCollection_t daughters;
        for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
        mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

        mela.setProcess(proc, me, prod);
        TUtil::setLeptonMassScheme(TVar::ConserveDifermionMass);

        bool hasConst = (tree->GetBranchStatus(strConstBranchname)==1);
        bool doCalc=(!hasConst || ctr==0);
        if (doCalc){
          mela.computeP(mesq_conserveDifermMass, false);
          if (hasConst) mesq_conserveDifermMass = mesq_calc / cconst_calc;

          alphasVal = mela.getIORecord()->getAlphaSatMZ();
          mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
          mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
        }
        else mesq_conserveDifermMass = mesq_calc / cconst_calc;

        mesq_conserveDifermMass /= pow(alphasVal, 2);
        if (fabs(aL1)>0. || fabs(aR1)>0.) mesq_conserveDifermMass /= pow(aL1, 2)+pow(aR1, 2);
        if (fabs(aL2)>0. || fabs(aR2)>0.) mesq_conserveDifermMass /= pow(aL2, 2)+pow(aR2, 2);

        bool doFill = !(
          isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
          );

        if (doFill){
          if (isweighted){
            float mh = samplePoleMasses[tree].first;
            float gh = samplePoleMasses[tree].second;
            float prop = 1./(pow(pow(GenHMass, 2) - pow(mh, 2), 2) + pow(mh*gh, 2));
            wgt *= prop;
            weightedBinList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0, wgt);
          }
          else binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0, 1.);
        }

        mela.resetInputEvent();
        ctr++;
      }

      cout << "Adjusting weights in weighted samples" << endl;
      weightedBinList.at(bin).adjustWeights();
      cout << "Merging event info. from weighted samples" << endl;
      binList.at(bin).mergeBin(weightedBinList.at(bin));
      cout << "Sifting" << endl;
      binList.at(bin).sift();

      for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
        mzz = binList.at(bin).masses.at(ev);
        mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
        wgt = binList.at(bin).weights.at(ev);
        useWeighted = (wgt!=1.);
        hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass, wgt);
        hvar->Fill(mzz, mzz, wgt);
        if (writeFinalTree) newtree->Fill();
      }
    }

    double* xexyey[4];
    for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[0][bin] = hvar->GetBinContent(bin+1);
      xexyey[1][bin] = hvar->GetBinError(bin+1);

      cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;
      xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
      xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
      xexyey[2][bin] = log10(xexyey[2][bin]);
    }


    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
    foutput->WriteTObject(tg);
    delete tg;

    for (int ix=0; ix<4; ix++) delete[] xexyey[ix];
    foutput->WriteTObject(hmesq_conserveDifermMass);
    foutput->WriteTObject(hvar);
    if (writeFinalTree) foutput->WriteTObject(newtree);
    if (writeFinalTree) delete newtree;
    delete hmesq_conserveDifermMass;
    delete hvar;
    delete[] binning;
  }
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
  foutput->Close();
}

/*
SPECIFIC COMMENT: OUTPUT ME DIVIDED BY
- (aL1**2+aR1**2)*(aL2**2+aR2**2) TO REMAIN INDEPENDENT OF CHANNEL
*/
void get_PAvgProfile_MCFM_ZZQQB_bkgZZ(bool recalculate=false){
  int erg_tev=8;
  float mPOLE=125.;
  TString TREE_NAME;
  const bool writeFinalTree=true;
  const TString strMEBranchname = "p_QQB_BKG_MCFM";
  const TString strConstBranchname = "pConst_QQB_BKG_MCFM";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::bkgZZ;
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod = TVar::ZZQQB;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  std::vector<short>* LepLepId=0;
  std::vector<float>* LepPt=0;
  std::vector<float>* LepEta=0;
  std::vector<float>* LepPhi=0;

  float mesq_calc=0., cconst_calc=1.;
  float mesq_conserveDifermMass=0;
  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  short Z1Flav, Z2Flav;
  int LepID[4];

  vector<pair<float, int>> index[3];
  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s.root", strme.Data(), strprod.Data(), strproc.Data()), "recreate");

  vector<TString> strSamples_13TeV = constructSamplesList("qq_Bkg", 13.);
  vector<TString> strSamples_8TeV = constructSamplesList("qq_Bkg", 8.);
  vector<TString> strSamples_7TeV = constructSamplesList("qq_Bkg", 7.);

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  TString cinput_main;

  TREE_NAME = "ZZTree/candTree";
  cinput_main = inputdir_13TeV;
  //for (int is=0; is<1; is++){
  for (int is=0; is<(int)strSamples_13TeV.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples_13TeV[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found with " << tree->GetEntries() << " events." << endl;
          bool doRecalculate = recalculate;
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0) doRecalculate = true;
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
          tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
          tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
          tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
          tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
          tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
          tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
          tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
          tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
          tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
          if (!doRecalculate && tree->GetBranchStatus(strConstBranchname)==0){
            tree->SetBranchStatus(strMEBranchname, 1); tree->SetBranchAddress(strMEBranchname, &mesq_calc);
            tree->SetBranchStatus(strConstBranchname, 1); tree->SetBranchAddress(strConstBranchname, &cconst_calc);
          }
          nEntries += tree->GetEntries();
          treeList.push_back(tree);
          finputList.push_back(finput);
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }

  for (int ic=0; ic<3; ic++){
    TREE_NAME = "SelectedTree";
    cinput_main = inputdir_8TeV;
    //for (int is=0; is<2; is++){
    for (int is=0; is<(int)strSamples_8TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_8TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found with " << tree->GetEntries() << " events." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
    cinput_main = inputdir_7TeV;
    //for (int is=0; is<0; is++){
    for (int is=0; is<(int)strSamples_7TeV.size(); is++){
      TString cinput = Form("%s/%s/%s", cinput_main.Data(), strchannel[ic].Data(), (strSamples_7TeV[is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found with " << tree->GetEntries() << " events." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &mzz);
            tree->SetBranchStatus("Z1Mass", 1); tree->SetBranchAddress("Z1Mass", &m1);
            tree->SetBranchStatus("Z2Mass", 1); tree->SetBranchAddress("Z2Mass", &m2);
            tree->SetBranchStatus("helcosthetaZ1", 1); tree->SetBranchAddress("helcosthetaZ1", &h1);
            tree->SetBranchStatus("helcosthetaZ2", 1); tree->SetBranchAddress("helcosthetaZ2", &h2);
            tree->SetBranchStatus("helphi", 1); tree->SetBranchAddress("helphi", &phi);
            tree->SetBranchStatus("costhetastar", 1); tree->SetBranchAddress("costhetastar", &hs);
            tree->SetBranchStatus("phistarZ1", 1); tree->SetBranchAddress("phistarZ1", &phi1);
            tree->SetBranchStatus("Z1ids", 1); tree->SetBranchAddress("Z1ids", &Z1Flav);
            tree->SetBranchStatus("Z2ids", 1); tree->SetBranchAddress("Z2ids", &Z2Flav);
            nEntries += tree->GetEntries();
            treeList.push_back(tree);
            finputList.push_back(finput);
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
  }

  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;

  unsigned int ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=
      (
      Z1Flav*Z2Flav==pow(13, 4)
      ||
      Z1Flav*Z2Flav==pow(11, 4)
      ||
      Z1Flav*Z2Flav==pow(11*13, 2)
      )
      ;
    if (!doProcess) continue;
    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    unsigned int ic = (Z1Flav*Z2Flav==pow(13, 4))*0 + (Z1Flav*Z2Flav==pow(11, 4))*1 + (Z1Flav*Z2Flav==pow(11*13, 2))*2;
    addByLowest(index[ic], mzz, ev);
    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;

  for (unsigned int ic=0; ic<3; ic++){
    float firstVal=index[ic].at(0).first;
    float lastVal=index[ic].at(index[ic].size()-1).first;
    float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

    float divisor=55000;
    int nbins = index[ic].size()/divisor;
    const int nbins_th=10/*50*/;
    while (nbins<nbins_th){
      if (divisor>1000) divisor -= 1000;
      else if (divisor>100) divisor -= 100;
      else break;
      nbins=index[ic].size()/divisor;
    }
    cout << "nbins=" << nbins << endl;
    if (nbins<3) cerr << "Not enough bins!" << endl;
    vector<ExtBin> binList;
    float* binning = new float[nbins+1];
    binning[0]=infimum;
    binning[nbins]=supremum;
    int ev_stepsize = index[ic].size()/nbins;
    cout << "Event step size: " << ev_stepsize << endl;
    cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
    for (int ix=1; ix<nbins; ix++){
      binning[ix]=(index[ic][ix*ev_stepsize-1].first+index[ic][ix*ev_stepsize].first)*0.5;
      ExtBin tmpbin;
      tmpbin.binlow = binning[ix-1]; tmpbin.binhigh = binning[ix];
      for (int bin=0; bin<ev_stepsize; bin++){
        int evid = index[ic][(ix-1)*ev_stepsize+bin].second;
        tmpbin.events.push_back(evid);
      }
      binList.push_back(tmpbin);
      cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ic][ix*ev_stepsize].second << ", step " << ix*ev_stepsize << "]" << endl;
    }
    ExtBin tmpbin;
    tmpbin.binlow = binning[nbins-1]; tmpbin.binhigh = binning[nbins];
    for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index[ic].size(); bin++){
      int evid = index[ic][bin].second;
      tmpbin.events.push_back(evid);
    }
    binList.push_back(tmpbin);
    cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
    cout << "Bin list has the following bins:" << endl;
    for (unsigned int ib=0; ib<binList.size(); ib++){
      cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    }

    foutput->cd();

    TProfile* hvar = new TProfile(Form("candMass_%s", strchannel[ic].Data()), "", nbins, binning); hvar->Sumw2();
    TProfile* hmesq_conserveDifermMass = new TProfile(Form("P_ConserveDifermionMass_%s", strchannel[ic].Data()), "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

    TTree* newtree=0;
    if (writeFinalTree){
      newtree = new TTree(Form("FinalTree_%s", strchannel[ic].Data()), "");
      newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
      newtree->Branch("ZZMass", &mzz);
    }

    mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

    if (ic==1){
      LepID[0]=11;
      LepID[1]=-11;
    }
    else{
      LepID[0]=13;
      LepID[1]=-13;
    }
    if (ic==0){
      LepID[2]=13;
      LepID[3]=-13;
    }
    else{
      LepID[2]=11;
      LepID[3]=-11;
    }

    double aL1=0, aR1=0, aL2=0, aR2=0, mz=0, gaz=0;
    unsigned int ctr=0;
    for (unsigned int bin=0; bin<binList.size(); bin++){
      cout << "Bin " << bin << " is now being scrutinized..." << endl;
      for (unsigned int ev = 0; ev < binList.at(bin).events.size(); ev++){
        int getEv = binList.at(bin).events.at(ev);
        getEntry(treeList, getEv);
        TTree* tree = findTree(treeList, getEv);

        if (ev%1000==0) cout << "Doing event " << getEv << endl;

        TLorentzVector pDaughters[4];
        std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
        for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); }
        SimpleParticleCollection_t daughters;
        for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
        mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

        mela.setProcess(proc, me, prod);
        TUtil::setLeptonMassScheme(TVar::ConserveDifermionMass);

        bool hasConst = (tree->GetBranchStatus(strConstBranchname)==1);
        bool doCalc=(!hasConst || ctr==0);
        if (doCalc){
          mela.computeP(mesq_conserveDifermMass, false);
          if (hasConst) mesq_conserveDifermMass = mesq_calc / cconst_calc;

          mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
          mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
          mz = mela.getPrimaryMass(23);
          gaz = mela.getPrimaryWidth(23);
        }
        else mesq_conserveDifermMass = mesq_calc / cconst_calc;

        double propagator;
        if (fabs(aL1)>0. || fabs(aR1)>0.) mesq_conserveDifermMass /= pow(aL1, 2)+pow(aR1, 2);
        if (fabs(aL2)>0. || fabs(aR2)>0.) mesq_conserveDifermMass /= pow(aL2, 2)+pow(aR2, 2);
        if (fabs(mzz-mz)<=4.*gaz){
          double sh = pow(mzz, 2);
          double shdn = pow(mz-4.*gaz, 2);
          double shup = pow(mz+4.*gaz, 2);
          double prop_sh = 1./(pow(sh-pow(mz, 2), 2) + pow(mz*gaz, 2));
          double prop_shdn = 1./(pow(shdn-pow(mz, 2), 2) + pow(mz*gaz, 2));
          double prop_shup = 1./(pow(shup-pow(mz, 2), 2) + pow(mz*gaz, 2));
          double fsh = (sh-shdn)/(shup-shdn);
          propagator = prop_sh / (prop_shdn*(1.-fsh) + prop_shup*fsh);
        }
        else propagator=1.;
        mesq_conserveDifermMass /= propagator;

        bool doFill = !(
          isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)
          );

        if (doFill) binList.at(bin).addEvent(mzz, mesq_conserveDifermMass, 0);

        mela.resetInputEvent();
        ctr++;
      }

      binList.at(bin).sift(); binList.at(bin).adjustWeights();

      for (unsigned int ev=0; ev<binList.at(bin).masses.size(); ev++){
        mzz = binList.at(bin).masses.at(ev);
        mesq_conserveDifermMass = binList.at(bin).mevals.at(ev);
        hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
        hvar->Fill(mzz, mzz);
        if (writeFinalTree) newtree->Fill();
      }
    }

    double* xexyey[4];
    for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
    for (int bin=0; bin<nbins; bin++){
      xexyey[0][bin] = hvar->GetBinContent(bin+1);
      xexyey[1][bin] = hvar->GetBinError(bin+1);

      cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;
      xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
      xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
      xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
      xexyey[2][bin] = log10(xexyey[2][bin]);
    }


    TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
    tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
    foutput->WriteTObject(tg);
    delete tg;

    for (int ix=0; ix<4; ix++) delete[] xexyey[ix];
    foutput->WriteTObject(hmesq_conserveDifermMass);
    foutput->WriteTObject(hvar);
    if (writeFinalTree) foutput->WriteTObject(newtree);
    if (writeFinalTree) delete newtree;
    delete hmesq_conserveDifermMass;
    delete hvar;
    delete[] binning;
  }
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
  foutput->Close();
}

void get_PAvgProfile_ANALYTICAL_ZZQQB_bkgZZ(){
  const int erg_tev=13;
  const float mPOLE=125.;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::bkgZZ;
  TVar::MatrixElement me = TVar::ANALYTICAL;
  TVar::Production prod = TVar::ZZQQB;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  double m1_low=40;
  double m2_low=12;

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s.root", strme.Data(), strprod.Data(), strproc.Data()), "recreate");

  cout << "Acquiring Mela measurables and pdf..." << endl;
  RooSpin::modelMeasurables measurables = mela.getMeasurablesRRV();
  const double xrange[2]={ max(57., (double)((int)(m1_low+m2_low+0.5))), 15000. };
  vector<double> masses;
  double massmin=xrange[0];
  double mass=massmin;
  while (mass<=xrange[1]){
    masses.push_back(mass);
    double massinc;
    if (mass<90.) massinc=1;
    else if (mass<122.) massinc=4;
    else if (mass<194.) massinc=2;
    else if (mass<200.) massinc=3;
    else if (mass<600.) massinc=20.;
    else if (mass<1500.) massinc=100.;
    else if (mass<3000.) massinc=250.;
    else if (mass<10000.) massinc=500.;
    else massinc=1000.;
    mass += massinc;
  }
  const unsigned int npoints=masses.size();
  ((RooRealVar*)measurables.m12)->setRange(xrange[0], xrange[1]);
  ((RooRealVar*)measurables.m1)->setRange(m1_low, 120); ((RooRealVar*)measurables.m1)->setVal(m1_low);
  ((RooRealVar*)measurables.m2)->setRange(m2_low, 120); ((RooRealVar*)measurables.m2)->setVal(m2_low);
  double* xy[2];
  for (int i=0; i<2; i++) xy[i] = new double[npoints];

  mela.upFrac_rrv->setVal(1.);
  mela.upFrac_rrv->setConstant(true);

  RooAbsPdf* pdf = mela.qqZZmodel;

  pdf->defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
  pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 100);;
  pdf->defaultIntegratorConfig()->method2D().setLabel("RooAdaptiveGaussKronrodIntegrator2D");
  pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator2D").setRealValue("maxSeg", 100);;
  pdf->defaultIntegratorConfig()->methodND().setLabel("RooAdaptiveGaussKronrodIntegratorND");
  pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegratorND").setRealValue("maxSeg", 100);;
  pdf->defaultIntegratorConfig()->setEpsAbs(1e-5);
  pdf->defaultIntegratorConfig()->setEpsRel(1e-5);

  cout << "Computing pdf integral as a function of " << measurables.m12->GetName() << endl;
  RooRealIntegral* pdfInt = new RooRealIntegral(
    "pdfInt", "", *pdf,
    RooArgSet(
    *measurables.h1, *measurables.h2, *measurables.Phi,
    *measurables.hs, *measurables.Phi1,
    *measurables.m1, *measurables.m2/*,
    *measurables.Y*/
    )
    );
  pdfInt->Print("v");
  for (unsigned int ix=0; ix<npoints; ix++){
    xy[0][ix]=masses.at(ix);
    ((RooRealVar*)measurables.m12)->setVal(xy[0][ix]);
    xy[1][ix]=pdfInt->getVal();
    cout << "pdfInt(" << xy[0][ix] << ") = " << xy[1][ix] << endl;
    xy[1][ix] = log10(xy[1][ix]);
  }
  TGraph* tg = new TGraph(npoints, xy[0], xy[1]);
  tg->SetName("tg_anaPdfInt");

  foutput->WriteTObject(tg);
  delete tg;
  delete pdfInt;
  for (int i=0; i<2; i++) delete[] xy[i];

  foutput->Close();
}

void get_PAvgProfile_ANALYTICAL_ZZGG_HSMHiggs(){
  const int erg_tev=13;
  const float mPOLE=125.;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::ANALYTICAL;
  TVar::Production prod = TVar::ZZGG;

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);

  double m1_low=40;
  double m2_low=12;

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s.root", strme.Data(), strprod.Data(), strproc.Data()), "recreate");

  cout << "Acquiring Mela measurables and pdf..." << endl;
  RooSpin::modelMeasurables measurables = mela.getMeasurablesRRV();
  const double xrange[2]={ max(57., (double)((int)(m1_low+m2_low+0.5))), 15000. };
  vector<double> masses;
  double massmin=xrange[0];
  double mass=massmin;
  while (mass<=xrange[1]){
    masses.push_back(mass);
    double massinc;
    if (mass<90.) massinc=1;
    else if (mass<194.) massinc=2;
    else if (mass<200.) massinc=3;
    else if (mass<600.) massinc=20.;
    else if (mass<1500.) massinc=100.;
    else if (mass<3000.) massinc=250.;
    else if (mass<10000.) massinc=500.;
    else massinc=1000.;
    mass += massinc;
  }
  const unsigned int npoints=masses.size();
  ((RooRealVar*)measurables.m12)->setRange(xrange[0], xrange[1]);
  ((RooRealVar*)measurables.m1)->setRange(m1_low, 120); ((RooRealVar*)measurables.m1)->setVal(m1_low);
  ((RooRealVar*)measurables.m2)->setRange(m2_low, 120); ((RooRealVar*)measurables.m2)->setVal(m2_low);
  double* xy[2];
  for (int i=0; i<2; i++) xy[i] = new double[npoints];

  mela.upFrac_rrv->setVal(1.);
  mela.upFrac_rrv->setConstant(true);

  RooAbsPdf* pdf = mela.ggSpin0Model->getPDF();
  pdf->defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
  pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 100);;
  pdf->defaultIntegratorConfig()->method2D().setLabel("RooAdaptiveGaussKronrodIntegrator2D");
  pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator2D").setRealValue("maxSeg", 100);;
  pdf->defaultIntegratorConfig()->methodND().setLabel("RooAdaptiveGaussKronrodIntegratorND");
  pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegratorND").setRealValue("maxSeg", 100);;
  pdf->defaultIntegratorConfig()->setEpsAbs(1e-5);
  pdf->defaultIntegratorConfig()->setEpsRel(1e-5);

  cout << "Computing pdf integral as a function of " << measurables.m12->GetName() << endl;
  RooRealIntegral* pdfInt = new RooRealIntegral(
    "pdfInt", "", *pdf,
    RooArgSet(
    *measurables.h1, *measurables.h2, *measurables.Phi,
    *measurables.hs, *measurables.Phi1,
    *measurables.m1, *measurables.m2/*,
    *measurables.Y*/
    )
    );
  pdfInt->Print("v");
  for (unsigned int ix=0; ix<npoints; ix++){
    xy[0][ix]=masses.at(ix);
    ((RooRealVar*)measurables.m12)->setVal(xy[0][ix]);
    xy[1][ix]=pdfInt->getVal();
    cout << "pdfInt(" << xy[0][ix] << ") = " << xy[1][ix] << endl;
    xy[1][ix] = log10(xy[1][ix]);
  }
  TGraph* tg = new TGraph(npoints, xy[0], xy[1]);
  tg->SetName("tg_anaPdfInt");

  foutput->WriteTObject(tg);
  delete tg;
  delete pdfInt;
  for (int i=0; i<2; i++) delete[] xy[i];

  foutput->Close();
}

void get_PAvgProfile_ANALYTICAL_HadVH_HSMHiggs(TString strprod, int sqrts){
  if (!(strprod == "Had_ZH" || strprod == "Had_WH")) return;
  int erg_tev=sqrts;
  float mPOLE=125.;
  float mVPOLE;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  TVar::Process proc = TVar::HSMHiggs;
  TVar::MatrixElement me = TVar::ANALYTICAL;
  TVar::Production prod;
  for (int iprod=(int)TVar::JJVBF; iprod<(int)TVar::nProductions; iprod++){
    prod = (TVar::Production)iprod;
    if (TVar::ProductionName(prod)==strprod) break;
  }
  RooSpin::VdecayType pmode, dmode;
  if (prod==TVar::Had_ZH){
    pmode=RooSpin::kVdecayType_Zud;
    dmode=pmode;
    mVPOLE=mela.getPrimaryMass(23);
  }
  else{
    pmode=RooSpin::kVdecayType_Wany;
    dmode=pmode;
    mVPOLE=mela.getPrimaryMass(24);
  }

  TString strproc = ProcessName(proc);
  TString strme = MatrixElementName(me);


  cout << "Acquiring measurables and pdf..." << endl;
  const double xrange[2]={ 57, sqrts*1000.-1. };
  RooRealVar* m12 = new RooRealVar("m12", "m_{H} (GeV)", mPOLE, xrange[0], xrange[1]);
  RooRealVar* m1 = new RooRealVar("m1", "m1", 460, (mVPOLE)-10., sqrts*1000.);
  RooRealVar* m2 = new RooRealVar("m2", "m2", mVPOLE, mVPOLE-0.05, mVPOLE+0.05);
  RooRealVar* h1 = new RooRealVar("h1", "cos#theta_{V*}", -1, 1);
  RooRealVar* h2 = new RooRealVar("h2", "cos#theta_{V}", -1, 1);
  RooRealVar* Phi1 = new RooRealVar("Phi1", "#Phi_{V*}", -TMath::Pi(), TMath::Pi());
  RooRealVar* hs = new RooRealVar("hs", "cos#theta^{*}", -1, 1);
  RooRealVar* Phi = new RooRealVar("Phi", "#Phi", -TMath::Pi(), TMath::Pi());
  RooRealVar* Y = new RooRealVar("Y", "Y", 0, -4, 4);

  RooSpin::modelMeasurables measurables;
  measurables.h1 = h1;
  measurables.h2 = h2;
  measurables.Phi = Phi;
  measurables.m1 = m1;
  measurables.m2 = m2;
  measurables.m12 = m12;
  measurables.hs = hs;
  measurables.Phi1 = Phi1;
  measurables.Y = Y;

  ScalarPdfFactory_VH fac(measurables, sqrts, pmode, dmode, false);
  cout << "m2 is fixed to " << measurables.m2->getVal() << endl;
  fac.makeParamsConst(false);
  RooRealVar* g1List[8][2];
  RooRealVar* g2List[8][2];
  //RooRealVar* g3List[8][2];
  RooRealVar* g4List[8][2];
  for (int gg=0; gg<8; gg++){
    for (int im=0; im<2; im++){
      g1List[gg][im] = (RooRealVar*)fac.couplings.g1List[gg][im];
      g2List[gg][im] = (RooRealVar*)fac.couplings.g2List[gg][im];
      //g3List[gg][im] = (RooRealVar*)fac.couplings.g3List[gg][im];
      g4List[gg][im] = (RooRealVar*)fac.couplings.g4List[gg][im];
    }
  }
  g1List[0][0]->setVal(1);

  TFile* foutput = TFile::Open(Form("pAvgLinToLog_%s_%s_%s_%iTeV.root", strme.Data(), strprod.Data(), strproc.Data(), sqrts), "recreate");

  m12->setConstant(false);
  vector<double> masses;
  double massmin=xrange[0];
  double mass=massmin;
  while (mass<=xrange[1]){
    masses.push_back(mass);
    double massinc;
    if (mass<90.) massinc=1;
    else if (mass<194.) massinc=2;
    else if (mass<200.) massinc=3;
    else if (mass<600.) massinc=20.;
    else if (mass<1500.) massinc=100.;
    else if (mass<3000.) massinc=250.;
    else if (mass<10000.) massinc=500.;
    else massinc=1000.;
    mass += massinc;
  }
  const unsigned int npoints=masses.size();
  m12->setRange(xrange[0], xrange[1]);

  double* xy[2];
  for (int i=0; i<2; i++) xy[i] = new double[npoints];

  RooAbsPdf* pdf = fac.getPDF();
  /*
  pdf->defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
  pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 100);;
  pdf->defaultIntegratorConfig()->method2D().setLabel("RooAdaptiveGaussKronrodIntegrator2D");
  pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator2D").setRealValue("maxSeg", 100);;
  pdf->defaultIntegratorConfig()->methodND().setLabel("RooAdaptiveGaussKronrodIntegratorND");
  pdf->defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegratorND").setRealValue("maxSeg", 100);;
  pdf->defaultIntegratorConfig()->setEpsAbs(1e-5);
  pdf->defaultIntegratorConfig()->setEpsRel(1e-5);
  */

  cout << "Computing pdf integral as a function of " << measurables.m12->GetName() << endl;
  for (unsigned int ix=0; ix<npoints; ix++){
    xy[0][ix]=masses.at(ix);
    m12->setVal(xy[0][ix]);
    ((RooRealVar*)fac.parameters.mX)->setVal(m12->getVal());
    ((RooRealVar*)fac.parameters.gamX)->setVal(1./m12->getVal());

    float m1_low = xy[0][ix] + measurables.m2->getVal();
    m1->setRange(m1_low, 1000.*sqrts);
    m1->setVal((m1->getMin()+m1->getMax())/2.);
    cout << "m1 min, max, val = " << m1->getMin() << " " << m1->getMax() << " " << m1->getVal() << endl;

    RooRealIntegral* pdfInt = new RooRealIntegral(
      "pdfInt", "", *pdf,
      RooArgSet(
      *measurables.h1, *measurables.h2, *measurables.Phi,
      *measurables.hs, *measurables.Phi1//,
      //*measurables.m1//,
      /**measurables.m2,*/
      //*measurables.Y
      )
      );
    if (ix==0) pdfInt->Print("v");
    xy[1][ix]=pdfInt->getVal();
    delete pdfInt;

    cout << "pdfInt(" << xy[0][ix] << ") = " << xy[1][ix] << endl;
    xy[1][ix] = log10(xy[1][ix]);
  }
  TGraph* tg = new TGraph(npoints, xy[0], xy[1]);
  tg->SetName("tg_anaPdfInt");

  foutput->WriteTObject(tg);
  delete tg;
  for (int i=0; i<2; i++) delete[] xy[i];

  foutput->Close();

  delete m12;
  delete m1;
  delete m2;
  delete hs;
  delete h1;
  delete h2;
  delete Phi;
  delete Phi1;
  delete Y;
}


/*
SPECIFIC COMMENT: OUTPUT ME DIVIDED BY
- (aL1**2+aR1**2)*(aL2**2+aR2**2) TO REMAIN INDEPENDENT OF CHANNEL
*/
void get_PAvgProfile_MCFM_JJQCD_bkgZJets_13TeV_2l2q(){
  int erg_tev=13;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";
  const bool writeFinalTree=false;

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  float mesq_conserveDifermMass=0;
  float mzz;
  float m1;
  float m2;
  float h1;
  float h2;
  float phi;
  float hs;
  float phi1;
  vector<float>* mzz_array=0;
  vector<float>* m1_array=0;
  vector<float>* m2_array=0;
  vector<float>* h1_array=0;
  vector<float>* h2_array=0;
  vector<float>* phi_array=0;
  vector<float>* hs_array=0;
  vector<float>* phi1_array=0;
  vector<short>* ZZsel=0;
  vector<short>* ZZCandType=0;
  int LepID[4]={ 0, 0, 11, -11 };

  TString cinput_main = "/scratch0/hep/usarical/CJLST/LHC_13TeV/2l2q/";
  const int nSamples = 2;
  TString strSamples[nSamples]={
    "ZZ2l2qAnalysis_DY2JetsToLL.root",
    "ZZ2l2qAnalysis_DYJetsToLL.root"
  };

  TFile* foutput = new TFile("pAvgLinToLog_MCFM_JJQCD_bkgZJets_13TeV_2l2q.root", "recreate");

  gROOT->cd();
  TChain* tree = new TChain(TREE_NAME, "");
  for (int is=0; is<nSamples; is++) tree->Add(Form("%s/%s", cinput_main.Data(), (strSamples[is]).Data()));
  tree->SetBranchAddress("ZZCandType", &ZZCandType);
  tree->SetBranchAddress("ZZsel", &ZZsel);
  tree->SetBranchAddress("ZZMass", &mzz_array);
  tree->SetBranchAddress("Z1Mass", &m1_array);
  tree->SetBranchAddress("Z2Mass", &m2_array);
  tree->SetBranchAddress("helcosthetaZ1", &h1_array);
  tree->SetBranchAddress("helcosthetaZ2", &h2_array);
  tree->SetBranchAddress("helphi", &phi_array);
  tree->SetBranchAddress("costhetastar", &hs_array);
  tree->SetBranchAddress("phistarZ1", &phi1_array);

  int nTotalEntries = tree->GetEntries();

  TTree* tmptree = new TTree("IntermediateTree", "");
  tmptree->Branch("ZZMass", &mzz);
  tmptree->Branch("Z1Mass", &m1);
  tmptree->Branch("Z2Mass", &m2);
  tmptree->Branch("helcosthetaZ1", &h1);
  tmptree->Branch("helcosthetaZ2", &h2);
  tmptree->Branch("helphi", &phi);
  tmptree->Branch("costhetastar", &hs);
  tmptree->Branch("phistarZ1", &phi1);

  unsigned int ctr=0;
  cout << "Nrawentries = " << nTotalEntries << endl;
  for (int ev = 0; ev < nTotalEntries; ev++){
    tree->GetEntry(ev);
    for (unsigned int j=0; j < ZZsel->size(); j++){
      if (ZZsel->at(j)>=90 && (ZZCandType->at(j)==2 || ZZCandType->at(j)==1)){
        if (ctr%100==0) cout << "Event " << ctr << " being recorded at actual event " << ev << "." << endl;
        mzz=mzz_array->at(j);
        m1=m1_array->at(j);
        h1=h1_array->at(j);
        m2=m2_array->at(j);
        h2=h2_array->at(j);
        hs=hs_array->at(j);
        phi=phi_array->at(j);
        phi1=phi1_array->at(j);
        tmptree->Fill();
        ctr++;
      }
    }
  }

  nTotalEntries=tmptree->GetEntries();
  cout << "Ntotalentries = " << nTotalEntries << endl;

  TTree* tmptree2 = new TTree("IntermediateTree2", "");
  tmptree2->Branch("ZZMass", &mzz);
  tmptree2->Branch("Z1Mass", &m1);
  tmptree2->Branch("Z2Mass", &m2);
  tmptree2->Branch("helcosthetaZ1", &h1);
  tmptree2->Branch("helcosthetaZ2", &h2);
  tmptree2->Branch("helphi", &phi);
  tmptree2->Branch("costhetastar", &hs);
  tmptree2->Branch("phistarZ1", &phi1);

  TRandom3 randomthrow(1234567);
  double portion_to_keep = 1;
  if (nTotalEntries>1000000) portion_to_keep = 9.95e5/nTotalEntries;
  for (int ev = 0; ev < nTotalEntries; ev++){
    tmptree->GetEntry(ev);
    double rndnum = randomthrow.Uniform();
    if (rndnum<=portion_to_keep) tmptree2->Fill();
  }
  delete tmptree;
  tmptree=tmptree2;

  const int nEntries = tmptree->GetEntries();
  if (nEntries>=1000000){
    cerr << "TMath::Sort will experience problems. Aborting!" << endl;
    delete tmptree;
    delete tree;
    assert(0);
  }
  int* index = new int[nEntries];
  tmptree->Draw("ZZMass", "", "goff");
  TMath::Sort(nEntries, tmptree->GetV1(), index, false);

  tmptree->GetEntry(index[0]);
  float firstVal=mzz;
  tmptree->GetEntry(index[nEntries-1]);
  float lastVal=mzz;
  float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
  float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
  cout << "Nentries = " << nEntries << " | mzz = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

  int nbins=0;
  int divisor=11000;
  while (nbins<50){
    if (divisor>1000) divisor -= 1000;
    else if (divisor>100) divisor -= 100;
    else break;
    nbins=nEntries/divisor+1;
  }
  cout << "nbins=" << nbins << endl;
  if (nbins<3) cerr << "Not enough bins!" << endl;
  float* binning = new float[nbins+1];
  binning[0]=infimum;
  binning[nbins]=supremum;
  int ev_stepsize = nEntries/nbins;
  cout << "Event step size: " << ev_stepsize << endl;
  cout << "Boundary (" << 0 << ") = " << binning[0] << endl;
  for (int ix=1; ix<nbins; ix++){
    int ev = index[ix*ev_stepsize];
    tmptree->GetEntry(ev);
    float bhigh = mzz;
    ev = index[ix*ev_stepsize-1];
    float blow = mzz;
    binning[ix]=(bhigh+blow)*0.5;
    cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << ev << ", step " << ix*ev_stepsize << "]" << endl;
  }
  cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
  delete[] index;

  foutput->cd();
  TProfile* hvar = new TProfile("candMass", "", nbins, binning); hvar->Sumw2();
  TProfile* hmesq_conserveDifermMass = new TProfile("P_ConserveDifermionMass", "", nbins, binning); hmesq_conserveDifermMass->Sumw2();

  TTree* newtree=0;
  if (writeFinalTree){
    newtree = new TTree("FinalTree", "");
    newtree->Branch("mesq_conserveDifermMass", &mesq_conserveDifermMass);
    newtree->Branch("ZZMass", &mzz);
  }

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  for (int ev = 0; ev < nEntries; ev++){
    tmptree->GetEntry(ev); // No need for ordering anymore
    if (ev%10000==0) cout << "Doing event " << ev << endl;

    TLorentzVector pDaughters[4];
    std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
    for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); }
    SimpleParticleCollection_t daughters;
    for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
    mela.setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);

    //double alphasVal;
    bool doFill=true;
    mela.setProcess(TVar::bkgZJets, TVar::MCFM, TVar::JJQCD);

    TUtil::setLeptonMassScheme(TVar::ConserveDifermionMass);
    mela.computeP(mesq_conserveDifermMass, false);
    double aL1, aR1, aL2, aR2;
    mela.getIORecord()->getVDaughterCouplings(aL1, aR1, 0);
    mela.getIORecord()->getVDaughterCouplings(aL2, aR2, 1);
    if (fabs(aL1)>0. || fabs(aR1)>0.) mesq_conserveDifermMass /= pow(aL1, 2)+pow(aR1, 2);
    if (fabs(aL2)>0. || fabs(aR2)>0.) mesq_conserveDifermMass /= pow(aL2, 2)+pow(aR2, 2);
    cout << "aL1: " << aL1 << '\t';
    cout << "aR1: " << aR1 << '\t';
    cout << "aL2: " << aL2 << '\t';
    cout << "aR2: " << aR2 << endl;
    //mesq_conserveDifermMass = log10(mesq_conserveDifermMass);
    if (isnan(mesq_conserveDifermMass) || isinf(mesq_conserveDifermMass)) doFill=false;

    if (doFill){
      hmesq_conserveDifermMass->Fill(mzz, mesq_conserveDifermMass);
      hvar->Fill(mzz, mzz);
    }

    if (writeFinalTree) newtree->Fill();
    mela.resetInputEvent();
  }

  double* xexyey[4];
  for (int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
  for (int bin=0; bin<nbins; bin++){
    xexyey[0][bin] = hvar->GetBinContent(bin+1);
    xexyey[1][bin] = hvar->GetBinError(bin+1);

    cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;
    xexyey[2][bin] = hmesq_conserveDifermMass->GetBinContent(bin+1);
    xexyey[3][bin] = hmesq_conserveDifermMass->GetBinError(bin+1);
    xexyey[3][bin] = log10(xexyey[3][bin])/xexyey[2][bin];
    xexyey[2][bin] = log10(xexyey[2][bin]);
  }


  TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
  tg->SetName(Form("tg_%s", hmesq_conserveDifermMass->GetName()));
  foutput->WriteTObject(tg);
  delete tg;

  for (int ix=0; ix<4; ix++) delete[] xexyey[ix];
  foutput->WriteTObject(hmesq_conserveDifermMass);
  foutput->WriteTObject(hvar);
  if (writeFinalTree) foutput->WriteTObject(newtree);
  if (writeFinalTree) delete newtree;
  delete hmesq_conserveDifermMass;
  delete hvar;
  delete[] binning;
  delete tmptree;
  delete tree;

  foutput->Close();
}

/* SPECIFIC COMMENT: Convert a TGraph to a TSpline3 */
TSpline3* convertGraphToSpline3(TGraph* tg, double* dfirst=0, double* dlast=0){
  unsigned int nbins = tg->GetN();
  double* xy[2]={
    tg->GetX(),
    tg->GetY()
  };
  TSpline3* spline = new TSpline3("spline", tg, "b2e2", 0, 0);
  spline->SetName(Form("sp_%s", tg->GetName()));
  if (dfirst!=0) *dfirst = spline->Derivative(xy[0][0]);
  if (dlast!=0) *dlast = spline->Derivative(xy[0][nbins-1]);
  return spline;
}

/* SPECIFIC COMMENT: Convert a TGraph to a TSpline3 */
TSpline3* convertGraphToSpline3_FaithfulSlopes(TGraph* tg, double* dfirst=0, double* dlast=0, bool faithfulLow=true, bool faithfulHigh=true){
  if (!faithfulLow&&!faithfulHigh) return convertGraphToSpline3(tg, dfirst, dlast);
  else{
    unsigned int nbins = tg->GetN();
    double* xy[2]={
      tg->GetX(),
      tg->GetY()
    };
    double derivative_first = (xy[1][1]-xy[1][0])/(xy[0][1]-xy[0][0]);
    double derivative_last = (xy[1][nbins-1]-xy[1][nbins-2])/(xy[0][nbins-1]-xy[0][nbins-2]);
    TSpline3* spline;
    if (faithfulLow && faithfulHigh) spline = new TSpline3("spline", tg, "b1e1", derivative_first, derivative_last);
    else if (!faithfulLow && faithfulHigh) spline = new TSpline3("spline", tg, "b2e1", 0, derivative_last);
    else if (faithfulLow && !faithfulHigh) spline = new TSpline3("spline", tg, "b1e2", derivative_first, 0);
    else return nullptr; // Should never reach here, just to avoid warnings
    spline->SetName(Form("sp_%s", tg->GetName()));
    if (dfirst!=0) *dfirst = spline->Derivative(xy[0][0]);
    if (dlast!=0) *dlast = spline->Derivative(xy[0][nbins-1]);
    return spline;
  }
}
/* SPECIFIC COMMENT: Convert a TGraph to a TSpline5 */
TSpline5* convertGraphToSpline5(TGraph* tg, double* dfirst=0, double* dlast=0){
  unsigned int nbins = tg->GetN();
  double* xy[2]={
    tg->GetX(),
    tg->GetY()
  };
  double derivative_first = (xy[1][1]-xy[1][0])/(xy[0][1]-xy[0][0]);
  double derivative_last = (xy[1][nbins-1]-xy[1][nbins-2])/(xy[0][nbins-1]-xy[0][nbins-2]);
  TSpline5* spline = new TSpline5("spline", tg, "b1e1", derivative_first, derivative_last);
  spline->SetName(Form("sp_%s", tg->GetName()));
  if (dfirst!=0) *dfirst=derivative_first;
  if (dlast!=0) *dlast=derivative_last;
  return spline;
}
TSpline3* convertTSpline5ToTspline3(TSpline5* sp){
  double xmin = sp->GetXmin();
  double xmax = sp->GetXmax();
  const int nbins=500;
  double xyval[2][nbins+1];
  double interval = (xmax-xmin)/nbins;
  for (int ix=0; ix<=nbins; ix++){
    xyval[0][ix] = xmin + ix*interval;
    xyval[1][ix] = sp->Eval(xyval[0][ix]);
  }
  double dfirst = sp->Derivative(xmin);
  double dlast = sp->Derivative(xmax);
  TSpline3* sp_new = new TSpline3("spline", xyval[0], xyval[1], nbins+1, "b1e1", dfirst, dlast);
  sp_new->SetName(sp->GetName());
  sp_new->SetTitle(sp->GetTitle());
  return sp_new;
}

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*exp(x) */
TF1* getFcn_a0(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);

  double a0;
  a0 = y;

  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TF1* fcn = new TF1(fcnName, "[0]", xmin, xmax);
  fcn->SetParameter(0, a0);

  return fcn;
}

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*exp(x) */
TF1* getFcn_a0plusa1expX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);

  double a0, a1;
  a0 = y-s;
  a1 = s*exp(-x);

  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TF1* fcn = new TF1(fcnName, "[0]+[1]*exp(x)", xmin, xmax);
  fcn->SetParameter(0, a0);
  fcn->SetParameter(1, a1);

  return fcn;
}

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1/x */
TF1* getFcn_a0plusa1overX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);

  double a0, a1;
  a0 = y+s*x;
  a1 = -s*pow(x, 2);

  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TF1* fcn = new TF1(fcnName, "[0]+[1]/x", xmin, xmax);
  fcn->SetParameter(0, a0);
  fcn->SetParameter(1, a1);

  return fcn;
}

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0/x**2-a1/x */
TF1* getFcn_a0overX2minusa1overX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);

  double a0, a1;
  a0 = -y*pow(x, 2)-s*pow(x, 3);
  a1 = -2.*y*x-s*pow(x, 2);

  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TF1* fcn = new TF1(fcnName, "[0]/x/x-[1]/x", xmin, xmax);
  fcn->SetParameter(0, a0);
  fcn->SetParameter(1, a1);

  return fcn;
}

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1/x-(a1/x)**2 */
TF1* getFcn_a0plusXPinvminusXpsqinv(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);

  double a0, a1;
  double disc = 1.+8.*x*s;
  if (disc>0.){
    a1 = x/4.*(1.+sqrt(disc));
    a0 = y-a1/x+pow(a1/x, 2);
    cout << x << '\t' << a0 << '\t' << a1 << '\t' << s << '\t' << y << endl;

    TString fcnName;
    if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
    else fcnName = Form("highFcn_%s", sp->GetName());
    TF1* fcn = new TF1(fcnName, "[0]+[1]/x-pow([1]/x, 2)", xmin, xmax);
    fcn->SetParameter(0, a0);
    fcn->SetParameter(1, a1);
    return fcn;
  }
  else return getFcn_a0plusa1overX(sp, xmin, xmax, useLowBound);
}

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*x */
TF1* getFcn_a0plusa1timesX(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);

  double a0, a1;
  a0 = y-s*x;
  a1 = s;

  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TF1* fcn = new TF1(fcnName, "[0]+[1]*x", xmin, xmax);
  fcn->SetParameter(0, a0);
  fcn->SetParameter(1, a1);

  return fcn;
}

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*x+a2*x**2 */
TF1* getFcn_a0plusa1timesXplusa2overX2(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, y, s, d;
  double xh, sh;
  if (useLowBound){
    x = sp->GetXmin();
    xh = x+0.001;
  }
  else{
    x = sp->GetXmax();
    xh = x-0.001;
  }
  y = sp->Eval(x);
  s = sp->Derivative(x);
  sh = sp->Derivative(xh);
  d = (sh-s)/(xh-x);

  double a0, a1, a2;
  a2 = d*pow(x, 4)/6.;
  a1 = s+2.*a2/pow(x, 3);
  a0 = y-a1*x-a2/pow(x, 2);

  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TF1* fcn = new TF1(fcnName, "[0]+[1]*x+[2]/pow(x, 2)", xmin, xmax);
  fcn->SetParameter(0, a0);
  fcn->SetParameter(1, a1);
  fcn->SetParameter(2, a2);

  return fcn;
}

/* SPECIFIC COMMENT: Get a1 and a2 as well as a TF1 object for the formula a0+a1*x+a2/x**2 */
TF1* getFcn_a0plusa1timesXplusa2timesX2(TSpline3* sp, double xmin, double xmax, bool useLowBound){
  double x, y, s, d;
  double xh, sh;
  if (useLowBound){
    x = sp->GetXmin();
    xh = x+0.001;
  }
  else{
    x = sp->GetXmax();
    xh = x-0.001;
  }
  y = sp->Eval(x);
  s = sp->Derivative(x);
  sh = sp->Derivative(xh);
  d = (sh-s)/(xh-x);

  double a0, a1, a2;
  a2 = d/2.;
  a1 = s-2.*a2*x;
  a0 = y-a1*x-a2*pow(x, 2);

  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TF1* fcn = new TF1(fcnName, "[0]+[1]*x+[2]*pow(x, 2)", xmin, xmax);
  fcn->SetParameter(0, a0);
  fcn->SetParameter(1, a1);
  fcn->SetParameter(2, a2);

  return fcn;
}

/* SPECIFIC COMMENT: Get external input added by a0 and multiplied by a1 */
TGraph* getPatch_a0a1(TSpline3* sp, TSpline3* sppatch, double xmin, double xmax, bool useLowBound, bool forceOutput=true){
  double x, y, s;
  double y_patch, s_patch;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);
  y_patch = sppatch->Eval(x);
  s_patch = sppatch->Derivative(x);

  double a0, a1;
  a1 = s/s_patch;
  a0 = y-a1*y_patch;

  int n = (xmax-xmin+0.5);
  cout << "Patcher n: " << n << endl;
  double xwidth = (xmax-xmin)/n;
  double* xy[2];
  for (unsigned int i=0; i<2; i++) xy[i]=new double[n+1];
  for (int ip=0; ip<=n; ip++){
    double xval = xmin+xwidth*ip;
    double yval = sppatch->Eval(xval);
    yval = a0+a1*yval;
    xy[0][ip]=xval;
    xy[1][ip]=yval;
  }
  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TGraph* fcn = new TGraph(n+1, xy[0], xy[1]);
  fcn->SetName(fcnName);

  for (unsigned int i=0; i<2; i++) delete[] xy[i];
  return fcn;
}

/* SPECIFIC COMMENT: Get external input added by a0+a1/x */
TGraph* getPatch_a0plusa1overX(TSpline3* sp, TSpline3* sppatch, double xmin, double xmax, bool useLowBound, bool forceOutput=true){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);
  y -= sppatch->Eval(x);
  s -= sppatch->Derivative(x);

  double a0, a1;
  a0 = y+s*x;
  a1 = -s*pow(x, 2);

  cout << "getPatch_a0plusa1overX: Checking if patch will flip: " << endl;
  double xext;
  if (useLowBound) xext = x;
  else xext = x;
  if (1.-a1/pow(xext, 2)/sppatch->Derivative(xext)<0.){
    cout << "Sign flips at " << xext << endl;
    if (!forceOutput) return 0;
  }
  else cout << "Sign does not flip at " << xext << endl;

  int n = (xmax-xmin+0.5);
  cout << "Patcher n: " << n << endl;
  double xwidth = (xmax-xmin)/n;
  double* xy[2];
  for (unsigned int i=0; i<2; i++) xy[i]=new double[n+1];
  for (int ip=0; ip<=n; ip++){
    double xval = xmin+xwidth*ip;
    if (xval==0.) xval += xwidth*0.5;
    double yval = sppatch->Eval(xval);
    yval += (a0+a1/xval);
    xy[0][ip]=xval;
    xy[1][ip]=yval;
  }
  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TGraph* fcn = new TGraph(n+1, xy[0], xy[1]);
  fcn->SetName(fcnName);

  for (unsigned int i=0; i<2; i++) delete[] xy[i];
  return fcn;
}

/* SPECIFIC COMMENT: Get external input added by a0+a1*x */
TGraph* getPatch_a0plusa1timesX(TSpline3* sp, TSpline3* sppatch, double xmin, double xmax, bool useLowBound, bool forceOutput=true){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);
  y -= sppatch->Eval(x);
  s -= sppatch->Derivative(x);

  double a0, a1;
  a1 = s;
  a0 = y-a1*x;

  cout << "getPatch_a0plusa1timesX: Checking if patch will flip: " << endl;
  double xext;
  if (useLowBound) xext = x;
  else xext = x;
  if (1.+a1/sppatch->Derivative(xext)<0.){
    cout << "Sign flips at " << xext << endl;
    if (!forceOutput) return 0;
  }
  else cout << "Sign does not flip at " << xext << endl;

  int n = (xmax-xmin+0.5);
  cout << "Patcher n: " << n << endl;
  double xwidth = (xmax-xmin)/n;
  double* xy[2];
  for (unsigned int i=0; i<2; i++) xy[i]=new double[n+1];
  for (int ip=0; ip<=n; ip++){
    double xval = xmin+xwidth*ip;
    double yval = sppatch->Eval(xval);
    yval += (a0+a1*xval);
    xy[0][ip]=xval;
    xy[1][ip]=yval;
  }
  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TGraph* fcn = new TGraph(n+1, xy[0], xy[1]);
  fcn->SetName(fcnName);

  for (unsigned int i=0; i<2; i++) delete[] xy[i];
  return fcn;
}

/* SPECIFIC COMMENT: Get external input added by a0*exp(a1*x) */
TGraph* getPatch_a0expa1timesX(TSpline3* sp, TSpline3* sppatch, double xmin, double xmax, bool useLowBound){
  double x, y, s;
  if (useLowBound) x = sp->GetXmin();
  else x = sp->GetXmax();
  y = sp->Eval(x);
  s = sp->Derivative(x);
  y -= sppatch->Eval(x);
  s -= sppatch->Derivative(x);

  double a0, a1;
  a1 = s/y;
  a0 = y/exp(a1*x);

  int n = (xmax-xmin+0.5);
  cout << "Patcher n: " << n << endl;
  double xwidth = (xmax-xmin)/n;
  double* xy[2];
  for (unsigned int i=0; i<2; i++) xy[i]=new double[n+1];
  for (int ip=0; ip<=n; ip++){
    double xval = xmin+xwidth*ip;
    double yval = sppatch->Eval(xval);
    yval += (a0*exp(a1*xval));
    xy[0][ip]=xval;
    xy[1][ip]=yval;
  }
  TString fcnName;
  if (useLowBound) fcnName = Form("lowFcn_%s", sp->GetName());
  else fcnName = Form("highFcn_%s", sp->GetName());
  TGraph* fcn = new TGraph(n+1, xy[0], xy[1]);
  fcn->SetName(fcnName);

  for (unsigned int i=0; i<2; i++) delete[] xy[i];
  return fcn;
}

/* SPECIFIC COMMENT: THIS FUNCTION TAKES A TGRAPHERRORS AND MODIFIES IT DIRECTLY. THE OPTIONAL STD::VECTOR IS FOR FIXING <P> FOR CERTAIN X-VALUES. */
void regularizeSlice(TGraphErrors* tgSlice, std::vector<double>* fixedX=0, double omitbelow=0., int nIter_=-1, double threshold_ = -1){
  unsigned int nbins_slice = tgSlice->GetN();
  double* xexyey_slice[4]={
    tgSlice->GetX(),
    tgSlice->GetEX(),
    tgSlice->GetY(),
    tgSlice->GetEY()
  };

  double* xexyey_linear[4];
  for (unsigned int ix=0; ix<4; ix++){
    xexyey_linear[ix] = new double[nbins_slice];
    for (unsigned int iy=0; iy<nbins_slice; iy++){
      if (ix<2) xexyey_linear[ix][iy] = xexyey_slice[ix][iy];
      else if (ix==3) xexyey_linear[ix][iy] = exp(xexyey_slice[ix][iy]);
      else xexyey_linear[ix][iy] = xexyey_slice[ix][iy]*xexyey_linear[ix-1][iy];
    }
  }
  TGraphErrors* tgSlice_linear = new TGraphErrors(nbins_slice, xexyey_linear[0], xexyey_linear[2], xexyey_linear[1], xexyey_linear[3]);
  double integral_in=tgSlice_linear->Integral();
  delete tgSlice_linear;
  for (unsigned int ix=0; ix<4; ix++) delete[] xexyey_linear[ix];


  double* xexyey_mod[4];
  for (unsigned int ix=0; ix<4; ix++){
    xexyey_mod[ix] = new double[nbins_slice];
    for (unsigned int iy=0; iy<nbins_slice; iy++) xexyey_mod[ix][iy] = xexyey_slice[ix][iy];
  }
  unsigned int bin_first = 2, bin_last = nbins_slice-1;

  std::vector<int> fixedBins;
  if (fixedX!=0){
    for (unsigned int ifx=0; ifx<fixedX->size(); ifx++){
      double requestedVal = fixedX->at(ifx);
      double distance=1e15;
      int bin_to_fix=-1;
      for (unsigned int bin=0; bin<nbins_slice; bin++){ if (distance>fabs(xexyey_mod[0][bin]-requestedVal)){ bin_to_fix = bin; distance = fabs(xexyey_mod[0][bin]-requestedVal); } }
      if (bin_to_fix>=0) fixedBins.push_back(bin_to_fix);
      cout << "Requested to fix bin " << bin_to_fix << endl;
    }
  }
  if (omitbelow>0.){
    for (unsigned int bin=0; bin<nbins_slice; bin++){
      if (xexyey_mod[0][bin]<omitbelow) fixedBins.push_back(bin);
      cout << "Requested to fix bin " << bin << endl;
    }
  }

  double* xx_second;
  double* yy_second;

  int nIter = (nIter_<0 ? 1000 : nIter_);
  for (int it=0; it<nIter; it++){
    double threshold = (threshold_<0. ? 0.01 : threshold_);
    for (unsigned int binIt = bin_first; binIt<=bin_last; binIt++){
      bool doFix=false;
      for (unsigned int ifx=0; ifx<fixedBins.size(); ifx++){
        if ((int)(binIt-1)==fixedBins.at(ifx)){ doFix=true; /*cout << "Iteration " << it << " is fixing bin " << (binIt-1) << endl; */break; }
      }
      if (doFix) continue;

      int ctr = 0;
      int nbins_second = nbins_slice-1;
      xx_second = new double[nbins_second];
      yy_second = new double[nbins_second];
      for (unsigned int bin = 1; bin<=nbins_slice; bin++){
        if (bin==binIt) continue;
        xx_second[ctr] = xexyey_mod[0][bin-1];
        yy_second[ctr] = xexyey_mod[2][bin-1];
        ctr++;
      }

      TGraph* interpolator = new TGraph(nbins_second, xx_second, yy_second);
      double derivative_first = (yy_second[1]-yy_second[0])/(xx_second[1]-xx_second[0]);
      double derivative_last = (yy_second[nbins_second-1]-yy_second[nbins_second-2])/(xx_second[nbins_second-1]-xx_second[nbins_second-2]);
      TSpline3* spline = new TSpline3("spline", interpolator, "b1e1", derivative_first, derivative_last);

      double center = xexyey_mod[0][binIt-1];
      double val = spline->Eval(center);
      if (fabs(xexyey_mod[2][binIt-1]-val)>threshold*val) xexyey_mod[2][binIt-1]=val;

      delete spline;
      delete interpolator;

      delete[] yy_second;
      delete[] xx_second;
    }
  }

  for (unsigned int ix=0; ix<4; ix++){
    xexyey_linear[ix] = new double[nbins_slice];
    for (unsigned int iy=0; iy<nbins_slice; iy++){
      if (ix<2) xexyey_linear[ix][iy] = xexyey_mod[ix][iy];
      else if (ix==3) xexyey_linear[ix][iy] = exp(xexyey_mod[ix][iy]);
      else xexyey_linear[ix][iy] = xexyey_mod[ix][iy]*xexyey_linear[ix-1][iy];
    }
  }
  tgSlice_linear = new TGraphErrors(nbins_slice, xexyey_linear[0], xexyey_linear[2], xexyey_linear[1], xexyey_linear[3]);
  double integral_out=tgSlice_linear->Integral();
  delete tgSlice_linear;
  for (unsigned int ix=0; ix<4; ix++) delete[] xexyey_linear[ix];

  double scale = integral_out / integral_in;
  for (unsigned int iy=0; iy<nbins_slice; iy++){
    xexyey_slice[2][iy] = xexyey_mod[2][iy]*scale;
    xexyey_slice[3][iy] *= xexyey_mod[2][iy]/xexyey_slice[2][iy];
  }

  for (unsigned int ix=0; ix<4; ix++) delete[] xexyey_mod[ix];
}

TGraphErrors* removePointsBetween(TGraphErrors* tgSlice, double xmin, double xmax){
  const unsigned int nbins_slice = tgSlice->GetN();
  double* xexyey_slice[4]={
    tgSlice->GetX(),
    tgSlice->GetEX(),
    tgSlice->GetY(),
    tgSlice->GetEY()
  };

  double xexyey[4][nbins_slice];
  unsigned int ctr=0;
  for (unsigned int iy=0; iy<nbins_slice; iy++){
    if (xexyey_slice[0][iy]<=xmax && xexyey_slice[0][iy]>=xmin){
      cout << "Point to remove with X: " << xexyey_slice[0][iy] << endl;
      continue;
    }
    for (unsigned int ix=0; ix<4; ix++) xexyey[ix][ctr] = xexyey_slice[ix][iy];
    cout << "Point " << ctr << " X: " << xexyey[0][ctr] << endl;
    ctr++;
  }
  cout << "removePointsBetween: " << "Starting number of points " << nbins_slice << " final number " << ctr << endl;
  TGraphErrors* tgSlice_new = new TGraphErrors(ctr, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
  tgSlice_new->SetName(tgSlice->GetName());
  tgSlice_new->SetTitle(tgSlice->GetTitle());
  return tgSlice_new;
}

TGraphErrors* replacePointsBetween(TGraphErrors* tgSlice, double xmin, double xmax){
  const unsigned int nbins_slice = tgSlice->GetN();
  double* xexyey_slice[4]={
    tgSlice->GetX(),
    tgSlice->GetEX(),
    tgSlice->GetY(),
    tgSlice->GetEY()
  };

  double xexyey[4][nbins_slice];
  unsigned int lowbin=0, highbin=0;
  for (unsigned int iy=0; iy<nbins_slice; iy++){
    if (xexyey_slice[0][iy]<xmin) lowbin=iy;
    if (xexyey_slice[0][iy]>xmax){ highbin=iy; break; }
  }
  cout << "Low bin " << lowbin << " at " << xexyey_slice[0][lowbin] << endl;
  cout << "High bin " << highbin << " at " << xexyey_slice[0][highbin] << endl;

  for (unsigned int iy=0; iy<nbins_slice; iy++){
    if (xexyey_slice[0][iy]<=xmax && xexyey_slice[0][iy]>=xmin){
      for (unsigned int ix=0; ix<4; ix++){
        double xlow = xexyey_slice[0][lowbin];
        double xhigh = xexyey_slice[0][highbin];
        double ylow = xexyey_slice[ix][lowbin];
        double yhigh = xexyey_slice[ix][highbin];
        double val;
        if (ix==0) val = xexyey_slice[0][iy];
        else if (ix==2) val = ylow + (yhigh-ylow)/(xhigh-xlow)*(xexyey_slice[0][iy]-xlow);
        else val = xexyey_slice[ix][iy]*(xexyey[ix-1][iy]/xexyey_slice[ix-1][iy]);
        xexyey[ix][iy] = val;
      }
    }
    else{
      for (unsigned int ix=0; ix<4; ix++) xexyey[ix][iy] = xexyey_slice[ix][iy];
    }
    cout << "Point " << iy << " X: " << xexyey[0][iy] << endl;
  }
  TGraphErrors* tgSlice_new = new TGraphErrors(nbins_slice, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
  tgSlice_new->SetName(tgSlice->GetName());
  tgSlice_new->SetTitle(tgSlice->GetTitle());
  return tgSlice_new;
}

void addPoint(TGraphErrors*& tg, double x, double y, double ex, double ey){
  TString strname = tg->GetName();
  TString strtitle = tg->GetTitle();
  TString strxtitle = tg->GetXaxis()->GetTitle();
  TString strytitle = tg->GetYaxis()->GetTitle();

  vector<double> xarray;
  vector<double> yarray;
  vector<double> exarray;
  vector<double> eyarray;
  xarray.push_back(x);
  yarray.push_back(y);
  exarray.push_back(ex);
  eyarray.push_back(ey);
  for (int ip=0; ip<tg->GetN(); ip++){
    if (tg->GetX()[ip]!=x){
      xarray.push_back(tg->GetX()[ip]);
      yarray.push_back(tg->GetY()[ip]);
      exarray.push_back(tg->GetEX()[ip]);
      eyarray.push_back(tg->GetEY()[ip]);
    }
  }
  vector<pair<double, int>> xorder;
  for (unsigned int ip=0; ip<xarray.size(); ip++) addByLowest<double, int>(xorder, xarray.at(ip), ip);

  double* xynew[4];
  for (unsigned int i=0; i<4; i++) xynew[i] = new double[xorder.size()];
  for (unsigned int ip=0; ip<xarray.size(); ip++){
    unsigned int pos = xorder[ip].second;
    xynew[0][ip] = xarray[pos];
    xynew[1][ip] = yarray[pos];
    xynew[2][ip] = exarray[pos];
    xynew[3][ip] = eyarray[pos];
  }

  delete tg;

  tg = new TGraphErrors(xorder.size(), xynew[0], xynew[1], xynew[2], xynew[3]);
  tg->SetName(strname);
  tg->SetTitle(strtitle);
  tg->GetXaxis()->SetTitle(strxtitle);
  tg->GetYaxis()->SetTitle(strytitle);
  for (unsigned int i=0; i<4; i++) delete[] xynew[i];
}

TGraphErrors* addPoint(TGraphErrors* tgSlice, double x){
  const unsigned int nbins_slice = tgSlice->GetN();
  double* xexyey_slice[4]={
    tgSlice->GetX(),
    tgSlice->GetEX(),
    tgSlice->GetY(),
    tgSlice->GetEY()
  };

  double xexyey[4][nbins_slice+1];
  unsigned int lowbin=0, highbin=0;
  for (unsigned int iy=0; iy<nbins_slice; iy++){
    if (xexyey_slice[0][iy]<x) lowbin=iy;
    if (xexyey_slice[0][iy]>x){ highbin=iy; break; }
  }
  cout << "Low bin " << lowbin << " at " << xexyey_slice[0][lowbin] << endl;
  cout << "High bin " << highbin << " at " << xexyey_slice[0][highbin] << endl;

  int ctr=0;
  for (unsigned int iy=0; iy<nbins_slice; iy++){
    for (unsigned int ix=0; ix<4; ix++) xexyey[ix][ctr] = xexyey_slice[ix][iy];
    ctr++;
    if (iy==lowbin){
      for (unsigned int ix=0; ix<4; ix++){
        double ylow = xexyey_slice[ix][lowbin];
        double yhigh = xexyey_slice[ix][highbin];
        double val = (yhigh+ylow)*0.5;
        xexyey[ix][ctr] = val;
      }
      ctr++;
    }
  }
  TGraphErrors* tgSlice_new = new TGraphErrors(ctr, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
  tgSlice_new->SetName(tgSlice->GetName());
  tgSlice_new->SetTitle(tgSlice->GetTitle());
  return tgSlice_new;
}

TGraphErrors* addPointAfterBin(TGraphErrors* tgSlice, int abin){
  const unsigned int nbins_slice = tgSlice->GetN();
  double* xexyey_slice[4]={
    tgSlice->GetX(),
    tgSlice->GetEX(),
    tgSlice->GetY(),
    tgSlice->GetEY()
  };

  double xexyey[4][nbins_slice+1];
  unsigned int lowbin=abin, highbin=abin+1;
  cout << "Low bin " << lowbin << " at " << xexyey_slice[0][lowbin] << endl;
  cout << "High bin " << highbin << " at " << xexyey_slice[0][highbin] << endl;

  int ctr=0;
  for (unsigned int iy=0; iy<nbins_slice; iy++){
    for (unsigned int ix=0; ix<4; ix++) xexyey[ix][ctr] = xexyey_slice[ix][iy];
    ctr++;
    if (iy==lowbin){
      for (unsigned int ix=0; ix<4; ix++){
        double ylow = xexyey_slice[ix][lowbin];
        double yhigh = xexyey_slice[ix][highbin];
        double val = (yhigh+ylow)*0.5;
        xexyey[ix][ctr] = val;
      }
      ctr++;
      cout << "Adding additional point at " << xexyey[0][ctr-1] << '\t' << xexyey[2][ctr-1] << endl;
    }
  }
  TGraphErrors* tgSlice_new = new TGraphErrors(ctr, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
  tgSlice_new->SetName(tgSlice->GetName());
  tgSlice_new->SetTitle(tgSlice->GetTitle());
  return tgSlice_new;
}

struct PointRedivision{
  unsigned int npoints;
  double xlow;
  double xhigh;
  PointRedivision(unsigned int n, double xl, double xh) : npoints(n), xlow(xl), xhigh(xh){}
};

void rebinAverageME(
  TGraphErrors*& tg, TTree*& tree, TProfile*& pmass,
  PointRedivision& rediv,
  float& ZZMass, float& mesq, float& weight
){
  cout << "Old set of boundaries: " << endl;
  for (int ibin=1; ibin<=pmass->GetNbinsX()+1; ibin++) cout << pmass->GetXaxis()->GetBinLowEdge(ibin) << " ";
  cout << endl;

  int binlow=pmass->GetXaxis()->FindBin(rediv.xlow);
  int binhigh=pmass->GetXaxis()->FindBin(rediv.xhigh);
  float redivlow, redivhigh;
  if (binlow==0) binlow=1; redivlow=pmass->GetXaxis()->GetBinLowEdge(binlow);
  if (binhigh==0) binhigh=1; redivhigh=pmass->GetXaxis()->GetBinUpEdge(binhigh);
  cout << "Rebinning between " << redivlow << " and " << redivhigh << endl;

  vector<float> massvals;
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if (ZZMass<redivlow || ZZMass>=redivhigh) continue;
    massvals.push_back(ZZMass);
  }
  std::sort(massvals.begin(), massvals.end());
  cout << massvals.size() << " events in this range" << endl;
  vector<float> xboundaries;
  unsigned int nevts_collected=massvals.size();
  addByLowest(xboundaries, redivlow, true);
  addByLowest(xboundaries, redivhigh, true);
  for (unsigned int ip=1; ip<rediv.npoints; ip++) addByLowest(xboundaries, float((massvals.at(nevts_collected/rediv.npoints*ip) + massvals.at(nevts_collected/rediv.npoints*ip+1))*0.5), true);
  cout << "Additional x boundaries: ";
  for (unsigned int ip=0; ip<xboundaries.size(); ip++) cout << xboundaries.at(ip) << " ";
  cout << endl;

  TProfile* hvar = new TProfile("mass", "", rediv.npoints, xboundaries.data()); hvar->Sumw2();
  TProfile* hmesq = new TProfile("mesq", "", rediv.npoints, xboundaries.data()); hmesq->Sumw2();
  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);
    if (ZZMass<redivlow || ZZMass>=redivhigh) continue;
    hvar->Fill(ZZMass, ZZMass, weight);
    hmesq->Fill(ZZMass, mesq, weight);
  }

  { TGraphErrors* tgnew = removePointsBetween(tg, redivlow, redivhigh); delete tg; tg=tgnew; }
  for (int ibin=1; ibin<=hvar->GetNbinsX(); ibin++){
    double xadd=hvar->GetBinContent(ibin);
    double exadd=hvar->GetBinError(ibin);
    double yadd=hmesq->GetBinContent(ibin);
    double eyadd=hmesq->GetBinError(ibin);
    eyadd=log10(eyadd)/yadd;
    yadd=log10(yadd);
    addPoint(tg, xadd, yadd, exadd, eyadd);
    cout << "Added point " << xadd << ", " << yadd << " for range [ " << hvar->GetXaxis()->GetBinLowEdge(ibin) << ", " << hvar->GetXaxis()->GetBinUpEdge(ibin) << " ]" << endl;
  }
  delete hmesq;
  delete hvar;

  TString pmassname=pmass->GetName();
  TString pmasstitle=pmass->GetTitle();
  TString pmassxtitle=pmass->GetXaxis()->GetTitle();
  TString pmassytitle=pmass->GetYaxis()->GetTitle();
  vector<float> pmassboundaries;
  for (int ibin=1; ibin<=pmass->GetNbinsX()+1; ibin++){
    if (ibin>=binlow && ibin<=binhigh+1) continue;
    addByLowest(pmassboundaries, float(pmass->GetXaxis()->GetBinLowEdge(ibin)), true);
  }
  for (float& bb:xboundaries) addByLowest(pmassboundaries, bb, true);
  assert((int) pmassboundaries.size()==tg->GetN()+1);
  delete pmass;
  pmass = new TProfile(pmassname, pmasstitle, pmassboundaries.size()-1, pmassboundaries.data());
  for (int ibin=1; ibin<=tg->GetN(); ibin++){
    pmass->SetBinContent(ibin, tg->GetX()[ibin-1]);
    pmass->SetBinError(ibin, tg->GetEX()[ibin-1]);
  }

  cout << "New set of boundaries: " << endl;
  for (int ibin=1; ibin<=pmass->GetNbinsX()+1; ibin++) cout << pmass->GetXaxis()->GetBinLowEdge(ibin) << " ";
  cout << endl;
}


/* SPECIFIC COMMENT: NONE */
void generic_PAvgSmoothProducer(
  TVar::MatrixElement me, TVar::Production prod, TVar::Process proc,
  TF1* (*lowf)(TSpline3*, double, double, bool),
  TF1* (*highf)(TSpline3*, double, double, bool),
  int sqrts=-1,
  bool useFaithfulLowSlopes=false,
  bool useFaithfulHighSlopes=false,
  vector<PointRedivision>* redivision=nullptr
){
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);
  TString strproc = ProcessName(proc);
  TString strsqrts="";
  if (sqrts>0) strsqrts = Form("_%iTeV", sqrts);

  TString strappend = Form("%s_%s_%s%s", strme.Data(), strprod.Data(), strproc.Data(), strsqrts.Data());

  TFile* finput = TFile::Open(Form("pAvgLinToLog_%s%s", strappend.Data(), ".root"), "read");
  TFile* foutput = TFile::Open(Form("pAvgSmooth_%s%s", strappend.Data(), ".root"), "recreate");
  const unsigned int ngraphs=2;
  TString strtg[ngraphs]={
    "tg_P_ConserveDifermionMass",
    "tg_P_MomentumToEnergy"
  };
  const double xmin=0;
  const double xmax=(sqrts>0 ? (double) sqrts*1000. : 15000.);

  for (unsigned int ig=0; ig<ngraphs; ig++){
    finput->cd();
    TGraphErrors* tgin = 0;
    tgin = (TGraphErrors*) finput->Get(strtg[ig]);
    if (tgin==0){ cerr << strtg[ig] << " does not exist." << endl; continue; }
    TGraphErrors* tg = new TGraphErrors(*tgin);
    if (redivision){
      TTree* tree = (TTree*) finput->Get("FinalTree");
      if (!tree) cerr << "Final tree does not exist, so cannot apply redivision" << endl;
      else{
        float ZZMass, weight, mesq;
        tree->SetBranchAddress("ZZMass", &ZZMass);
        tree->SetBranchAddress("weight", &weight);
        if (ig==0) tree->SetBranchAddress("mesq_conserveDifermMass", &mesq);
        else tree->SetBranchAddress("mesq_jetPtoEScale", &mesq);
        TProfile* pmass = (TProfile*) finput->Get("candMass");
        for (PointRedivision& rediv:*redivision){
          rebinAverageME(
            tg, tree, pmass,
            rediv,
            ZZMass, mesq, weight
          );
        }
      }
    }

    foutput->WriteTObject(tg);

    int n = tg->GetN();
    double* xx = tg->GetX();
    double* ex = tg->GetEX();
    double* yy = tg->GetY();
    double* ey = tg->GetEY();

    TSpline3* sp = convertGraphToSpline3_FaithfulSlopes(tg, (double*)nullptr, (double*)nullptr, useFaithfulLowSlopes, useFaithfulHighSlopes);
    double tglow = xx[0];
    double tghigh = xx[tg->GetN()-1];
    TF1* lowFcn = lowf(sp, xmin, tglow, true);
    TF1* highFcn = highf(sp, tghigh, xmax, false);
    lowFcn->SetNpx((int) (tglow-xmin)*5);
    highFcn->SetNpx((int) (xmax-tghigh)*5);

    vector<pair<double, double>> points;
    for (double xval=xmin; xval<tglow; xval+=1){
      double yval = lowFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }
    for (int ix=0; ix<n; ix++){
      addByLowest<double, double>(points, xx[ix], yy[ix]);
    }
    int tghigh_int = ((int) ((tghigh+1.)/100.+0.5))*100;
    if (tghigh>=(double) tghigh_int) tghigh_int+=100;
    for (double xval=tghigh_int; xval<=xmax; xval+=100){
      double yval = highFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }

    int nn_new = points.size();
    cout << "Number of new points: " << nn_new-n << endl;
    double* xy_new[2];
    for (unsigned int i=0; i<2; i++) xy_new[i] = new double[nn_new];
    for (int ix=0; ix<nn_new; ix++){
      xy_new[0][ix] = points.at(ix).first;
      xy_new[1][ix] = points.at(ix).second;
    }

    delete highFcn;
    delete lowFcn;
    delete sp;

    TGraph* tg_updated = new TGraph(nn_new, xy_new[0], xy_new[1]);
    tg_updated->SetName(Form("%s_Smooth", tg->GetName()));
    foutput->WriteTObject(tg_updated);

    sp = convertGraphToSpline3(tg_updated);
    foutput->WriteTObject(sp);

    TCanvas* ctest = new TCanvas(Form("test_%s", strtg[ig].Data()), "", 8, 30, 800, 800);
    ctest->cd();
    tg->GetXaxis()->SetRangeUser(xmin, xmax);
    tg->Draw("ae1p");
    sp->Draw("csame");
    ctest->RedrawAxis();
    ctest->Modified();
    ctest->Update();
    foutput->WriteTObject(ctest);
    ctest->SaveAs(Form("%s%s", ctest->GetName(), ".pdf"));
    ctest->Close();

    delete sp;
    delete tg_updated;
    for (unsigned int i=0; i<2; i++) delete[] xy_new[i];
    delete tg;
  }
  foutput->Close();
  finput->Close();
}

/* SPECIFIC COMMENT: NONE */
void generic_PAvgSmoothProducer_withDecay(
  TVar::MatrixElement me, TVar::Production prod, TVar::Process proc,
  TF1* (*lowf)(TSpline3*, double, double, bool),
  TF1* (*highf)(TSpline3*, double, double, bool),
  int sqrts=-1,
  bool useFaithfulLowSlopes=false,
  bool useFaithfulHighSlopes=false,
  vector<PointRedivision>* redivision=nullptr
  ){
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);
  TString strproc = ProcessName(proc);
  TString strsqrts="";
  if (sqrts>0) strsqrts = Form("_%iTeV", sqrts);

  TString strappend = Form("%s_%s_%s%s", strme.Data(), strprod.Data(), strproc.Data(), strsqrts.Data());

  TFile* finput = TFile::Open(Form("pAvgLinToLog_%s%s", strappend.Data(), ".root"), "read");
  TFile* foutput = TFile::Open(Form("pAvgSmooth_%s%s", strappend.Data(), ".root"), "recreate");
  const unsigned int ngraphs=3;
  TString strtg[ngraphs]={
    "tg_P_ConserveDifermionMass_4mu",
    "tg_P_ConserveDifermionMass_4e",
    "tg_P_ConserveDifermionMass_2mu2e"
  };
  const double xmin=0;
  const double xmax=(sqrts>0 ? (double)sqrts*1000. : 15000.);

  for (unsigned int ig=0; ig<ngraphs; ig++){
    finput->cd();
    TGraphErrors* tgin = 0;
    tgin = (TGraphErrors*) finput->Get(strtg[ig]);
    if (tgin==0){ cerr << strtg[ig] << " does not exist." << endl; continue; }
    TGraphErrors* tg = new TGraphErrors(*tgin);
    if (redivision){
      TTree* tree;
      if (ig==0) tree = (TTree*) finput->Get("FinalTree_4mu");
      else if (ig==1) tree = (TTree*) finput->Get("FinalTree_4e");
      else tree = (TTree*) finput->Get("FinalTree_2mu2e");
      if (!tree) cerr << "Final tree does not exist, so cannot apply redivision" << endl;
      else{
        float ZZMass, weight, mesq;
        tree->SetBranchAddress("ZZMass", &ZZMass);
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("mesq_conserveDifermMass", &mesq);
        TProfile* pmass;
        if (ig==0) pmass = (TProfile*) finput->Get("candMass_4mu");
        else if (ig==1) pmass = (TProfile*) finput->Get("candMass_4e");
        else pmass = (TProfile*) finput->Get("candMass_2mu2e");
        for (PointRedivision& rediv:*redivision){
          rebinAverageME(
            tg, tree, pmass,
            rediv,
            ZZMass, mesq, weight
          );
        }
      }
    }

    foutput->cd();
    foutput->WriteTObject(tg);

    int n = tg->GetN();
    double* xx = tg->GetX();
    double* ex = tg->GetEX();
    double* yy = tg->GetY();
    double* ey = tg->GetEY();

    TSpline3* sp = convertGraphToSpline3_FaithfulSlopes(tg, (double*)nullptr, (double*)nullptr, useFaithfulLowSlopes, useFaithfulHighSlopes);
    double tglow = xx[0];
    double tghigh = xx[tg->GetN()-1];
    TF1* lowFcn = lowf(sp, xmin, tglow, true);
    TF1* highFcn = highf(sp, tghigh, xmax, false);
    lowFcn->SetNpx((int)(tglow-xmin)*5);
    highFcn->SetNpx((int)(xmax-tghigh)*5);

    vector<pair<double, double>> points;
    for (double xval=xmin; xval<tglow; xval+=1){
      double yval = lowFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }
    for (int ix=0; ix<n; ix++){
      addByLowest<double, double>(points, xx[ix], yy[ix]);
    }
    int tghigh_int = ((int)((tghigh+1.)/100.+0.5))*100;
    if (tghigh>=(double)tghigh_int) tghigh_int+=100;
    for (double xval=tghigh_int; xval<=xmax; xval+=100){
      double yval = highFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }

    int nn_new = points.size();
    cout << "Number of new points: " << nn_new-n << endl;
    double* xy_new[2];
    for (unsigned int i=0; i<2; i++) xy_new[i] = new double[nn_new];
    for (int ix=0; ix<nn_new; ix++){
      xy_new[0][ix] = points.at(ix).first;
      xy_new[1][ix] = points.at(ix).second;
    }

    delete highFcn;
    delete lowFcn;
    delete sp;

    TGraph* tg_updated = new TGraph(nn_new, xy_new[0], xy_new[1]);
    tg_updated->SetName(Form("%s_Smooth", tg->GetName()));
    foutput->WriteTObject(tg_updated);

    sp = convertGraphToSpline3(tg_updated);
    foutput->WriteTObject(sp);

    TCanvas* ctest = new TCanvas(Form("test_%s", strtg[ig].Data()), "", 8, 30, 800, 800);
    ctest->cd();
    tg->GetXaxis()->SetRangeUser(xmin, xmax);
    tg->Draw("ae1p");
    sp->Draw("csame");
    ctest->RedrawAxis();
    ctest->Modified();
    ctest->Update();
    foutput->WriteTObject(ctest);
    ctest->SaveAs(Form("%s%s", ctest->GetName(), ".pdf"));
    ctest->Close();

    delete sp;
    delete tg_updated;
    for (unsigned int i=0; i<2; i++) delete[] xy_new[i];
    delete tg;
  }
  foutput->Close();
  finput->Close();
}

  /* SPECIFIC COMMENT: NONE */
void generic_PAvgSmoothProducer_withDecay(
  TVar::MatrixElement me, TVar::Production prod, TVar::Process proc,
  TString strpatchpath, TString strpatchname,
  TGraph* (*lowpatcher)(TSpline3*, TSpline3*, double, double, bool, bool),
  TGraph* (*highpatcher)(TSpline3*, TSpline3*, double, double, bool, bool),
  int sqrts=-1,
  vector<PointRedivision>* redivision=nullptr
  ){
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);
  TString strproc = ProcessName(proc);
  TString strsqrts="";
  if (sqrts>0) strsqrts = Form("_%iTeV", sqrts);

  TString strappend = Form("%s_%s_%s%s", strme.Data(), strprod.Data(), strproc.Data(), strsqrts.Data());

  TFile* finput = TFile::Open(Form("pAvgLinToLog_%s%s", strappend.Data(), ".root"), "read");
  TFile* finput_patch = TFile::Open(strpatchpath, "read");
  TFile* foutput = TFile::Open(Form("pAvgSmooth_%s%s", strappend.Data(), ".root"), "recreate");
  const unsigned int ngraphs=3;
  TString strtg[ngraphs]={
    "tg_P_ConserveDifermionMass_4mu",
    "tg_P_ConserveDifermionMass_4e",
    "tg_P_ConserveDifermionMass_2mu2e"
  };

  TGraph* tgpatch = (TGraph*)finput_patch->Get(strpatchname);
  TSpline3* sppatch = convertGraphToSpline3(tgpatch);

  const double xmin=max(0., (tgpatch->GetX())[0]);
  const double xmax=min((sqrts>0 ? (double)sqrts*1000. : 15000.), (tgpatch->GetX())[tgpatch->GetN()-1]);

  for (unsigned int ig=0; ig<ngraphs; ig++){
    finput->cd();
    TGraphErrors* tgin = 0;
    tgin = (TGraphErrors*) finput->Get(strtg[ig]);
    if (tgin==0){ cerr << strtg[ig] << " does not exist." << endl; continue; }
    TGraphErrors* tg = new TGraphErrors(*tgin);
    if (redivision){
      TTree* tree;
      if (ig==0) tree = (TTree*) finput->Get("FinalTree_4mu");
      else if (ig==1) tree = (TTree*) finput->Get("FinalTree_4e");
      else tree = (TTree*) finput->Get("FinalTree_2mu2e");
      if (!tree) cerr << "Final tree does not exist, so cannot apply redivision" << endl;
      else{
        float ZZMass, weight, mesq;
        tree->SetBranchAddress("ZZMass", &ZZMass);
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("mesq_conserveDifermMass", &mesq);
        TProfile* pmass;
        if (ig==0) pmass = (TProfile*) finput->Get("candMass_4mu");
        else if (ig==1) pmass = (TProfile*) finput->Get("candMass_4e");
        else pmass = (TProfile*) finput->Get("candMass_2mu2e");
        for (PointRedivision& rediv:*redivision){
          rebinAverageME(
            tg, tree, pmass,
            rediv,
            ZZMass, mesq, weight
          );
        }
      }
    }

    foutput->cd();
    foutput->WriteTObject(tg);

    int n = tg->GetN();
    double* xx = tg->GetX();
    double* ex = tg->GetEX();
    double* yy = tg->GetY();
    double* ey = tg->GetEY();

    TSpline3* sp = convertGraphToSpline3(tg);
    sp->SetName(Form("%s_initial", sp->GetName()));
    foutput->WriteTObject(sp);
    double tglow = xx[0];
    double tghigh = xx[tg->GetN()-1];
    cout << "Extracting low patch..." << endl;
    TGraph* lowGraph = lowpatcher(sp, sppatch, xmin, tglow, true, true);
    cout << "Extracting high patch..." << endl;
    TGraph* highGraph = highpatcher(sp, sppatch, tghigh, xmax, false, true);
    cout << "Converting patches to splines..." << endl;
    TSpline3* lowFcn = convertGraphToSpline3(lowGraph);
    TSpline3* highFcn = convertGraphToSpline3(highGraph);

    vector<pair<double, double>> points;
    for (double xval=xmin; xval<tglow; xval+=1){
      double yval = lowFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }
    for (int ix=0; ix<n; ix++){
      addByLowest<double, double>(points, xx[ix], yy[ix]);
    }
    int tghigh_int = ((int)((tghigh+1.)/100.+0.5))*100;
    if (tghigh>=(double)tghigh_int) tghigh_int+=100;
    for (double xval=tghigh_int; xval<=xmax; xval+=100){
      double yval = highFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }

    int nn_new = points.size();
    cout << "Number of new points: " << nn_new-n << endl;
    double* xy_new[2];
    for (unsigned int i=0; i<2; i++) xy_new[i] = new double[nn_new];
    for (int ix=0; ix<nn_new; ix++){
      xy_new[0][ix] = points.at(ix).first;
      xy_new[1][ix] = points.at(ix).second;
    }

    delete highFcn;
    delete lowFcn;
    delete highGraph;
    delete lowGraph;
    delete sp;

    TGraph* tg_updated = new TGraph(nn_new, xy_new[0], xy_new[1]);
    tg_updated->SetName(Form("%s_Smooth", tg->GetName()));
    foutput->WriteTObject(tg_updated);

    sp = convertGraphToSpline3(tg_updated);
    foutput->WriteTObject(sp);

    TCanvas* ctest = new TCanvas(Form("test_%s", strtg[ig].Data()), "", 8, 30, 800, 800);
    ctest->cd();
    tg->GetXaxis()->SetRangeUser(xmin, xmax);
    tg->Draw("ae1p");
    sp->Draw("csame");
    ctest->RedrawAxis();
    ctest->Modified();
    ctest->Update();
    foutput->WriteTObject(ctest);
    ctest->SaveAs(Form("%s%s", ctest->GetName(), ".pdf"));
    ctest->Close();

    delete sp;
    delete tg_updated;
    for (unsigned int i=0; i<2; i++) delete[] xy_new[i];
    delete tg;
  }

  delete sppatch;
  foutput->Close();
  finput_patch->Close();
  finput->Close();
}

void generic_PAvgSmoothProducer_withDecay(
  TVar::MatrixElement me, TVar::Production prod, TVar::Process proc,
  TString strpatchpath, TString strpatchname,
  TString strpatchpath2, TString strpatchname2,
  TGraph* (*lowpatcher)(TSpline3*, TSpline3*, double, double, bool, bool),
  TGraph* (*highpatcher)(TSpline3*, TSpline3*, double, double, bool, bool),
  int sqrts=-1,
  vector<PointRedivision>* redivision=nullptr
  ){
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);
  TString strproc = ProcessName(proc);
  TString strsqrts="";
  if (sqrts>0) strsqrts = Form("_%iTeV", sqrts);

  TString strappend = Form("%s_%s_%s%s", strme.Data(), strprod.Data(), strproc.Data(), strsqrts.Data());

  TFile* finput = TFile::Open(Form("pAvgLinToLog_%s%s", strappend.Data(), ".root"), "read");
  TFile* finput_patch = TFile::Open(strpatchpath, "read");
  TFile* finput_patch2 = TFile::Open(strpatchpath2, "read");
  TFile* foutput = TFile::Open(Form("pAvgSmooth_%s%s", strappend.Data(), ".root"), "recreate");
  const unsigned int ngraphs=3;
  TString strtg[ngraphs]={
    "tg_P_ConserveDifermionMass_4mu",
    "tg_P_ConserveDifermionMass_4e",
    "tg_P_ConserveDifermionMass_2mu2e"
  };

  TGraph* tgpatch = (TGraph*)finput_patch->Get(strpatchname);
  TSpline3* sppatch = convertGraphToSpline3(tgpatch);

  TGraph* tgpatch2 = (TGraph*)finput_patch2->Get(strpatchname2);
  TSpline3* sppatch2 = convertGraphToSpline3(tgpatch2);

  const double xmin=max(0., (tgpatch->GetX())[0]);
  const double xmax=min((sqrts>0 ? (double)sqrts*1000. : 15000.), (tgpatch->GetX())[tgpatch->GetN()-1]);

  for (unsigned int ig=0; ig<ngraphs; ig++){
    finput->cd();
    TGraphErrors* tgin = 0;
    tgin = (TGraphErrors*) finput->Get(strtg[ig]);
    if (tgin==0){ cerr << strtg[ig] << " does not exist." << endl; continue; }
    TGraphErrors* tg = new TGraphErrors(*tgin);
    if (redivision){
      TTree* tree;
      if (ig==0) tree = (TTree*) finput->Get("FinalTree_4mu");
      else if (ig==1) tree = (TTree*) finput->Get("FinalTree_4e");
      else tree = (TTree*) finput->Get("FinalTree_2mu2e");
      if (!tree) cerr << "Final tree does not exist, so cannot apply redivision" << endl;
      else{
        float ZZMass, weight, mesq;
        tree->SetBranchAddress("ZZMass", &ZZMass);
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("mesq_conserveDifermMass", &mesq);
        TProfile* pmass;
        if (ig==0) pmass = (TProfile*) finput->Get("candMass_4mu");
        else if (ig==1) pmass = (TProfile*) finput->Get("candMass_4e");
        else pmass = (TProfile*) finput->Get("candMass_2mu2e");
        for (PointRedivision& rediv:*redivision){
          rebinAverageME(
            tg, tree, pmass,
            rediv,
            ZZMass, mesq, weight
          );
        }
      }
    }

    foutput->cd();
    foutput->WriteTObject(tg);

    int n = tg->GetN();
    double* xx = tg->GetX();
    double* ex = tg->GetEX();
    double* yy = tg->GetY();
    double* ey = tg->GetEY();

    TSpline3* sp = convertGraphToSpline3(tg);
    sp->SetName(Form("%s_initial", sp->GetName()));
    foutput->WriteTObject(sp);
    double tglow = xx[0];
    double tghigh = xx[tg->GetN()-1];
    cout << "Extracting low patch..." << endl;
    TGraph* lowGraph = lowpatcher(sp, sppatch, xmin, tglow, true, false);
    if (lowGraph==0) lowGraph = lowpatcher(sp, sppatch2, xmin, tglow, true, true);
    cout << "Extracting high patch..." << endl;
    TGraph* highGraph = highpatcher(sp, sppatch, tghigh, xmax, false, false);
    if (highGraph==0) highGraph = highpatcher(sp, sppatch2, tghigh, xmax, false, true);
    cout << "Converting patches to splines..." << endl;
    TSpline3* lowFcn = convertGraphToSpline3(lowGraph);
    TSpline3* highFcn = convertGraphToSpline3(highGraph);

    vector<pair<double, double>> points;
    for (double xval=xmin; xval<tglow; xval+=1){
      double yval = lowFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }
    for (int ix=0; ix<n; ix++){
      addByLowest<double, double>(points, xx[ix], yy[ix]);
    }
    int tghigh_int = ((int)((tghigh+1.)/100.+0.5))*100;
    if (tghigh>=(double)tghigh_int) tghigh_int+=100;
    for (double xval=tghigh_int; xval<=xmax; xval+=100){
      double yval = highFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }

    int nn_new = points.size();
    cout << "Number of new points: " << nn_new-n << endl;
    double* xy_new[2];
    for (unsigned int i=0; i<2; i++) xy_new[i] = new double[nn_new];
    for (int ix=0; ix<nn_new; ix++){
      xy_new[0][ix] = points.at(ix).first;
      xy_new[1][ix] = points.at(ix).second;
    }

    delete highFcn;
    delete lowFcn;
    delete highGraph;
    delete lowGraph;
    delete sp;

    TGraph* tg_updated = new TGraph(nn_new, xy_new[0], xy_new[1]);
    tg_updated->SetName(Form("%s_Smooth", tg->GetName()));
    foutput->WriteTObject(tg_updated);

    sp = convertGraphToSpline3(tg_updated);
    foutput->WriteTObject(sp);

    TCanvas* ctest = new TCanvas(Form("test_%s", strtg[ig].Data()), "", 8, 30, 800, 800);
    ctest->cd();
    tg->GetXaxis()->SetRangeUser(xmin, xmax);
    tg->Draw("ae1p");
    sp->Draw("csame");
    ctest->RedrawAxis();
    ctest->Modified();
    ctest->Update();
    foutput->WriteTObject(ctest);
    ctest->SaveAs(Form("%s%s", ctest->GetName(), ".pdf"));
    ctest->Close();

    delete sp;
    delete tg_updated;
    for (unsigned int i=0; i<2; i++) delete[] xy_new[i];
    delete tg;
  }

  delete sppatch2;
  delete sppatch;
  foutput->Close();
  finput_patch2->Close();
  finput_patch->Close();
  finput->Close();
}

/* SPECIFIC COMMENT: NONE */
void generic_PAvgSmoothProducer_withDecay(
  TVar::MatrixElement me, TVar::Production prod, TVar::Process proc,
  TString strpatchpath, TString strpatchname,
  TF1* (*lowf)(TSpline3*, double, double, bool),
  TGraph* (*highpatcher)(TSpline3*, TSpline3*, double, double, bool, bool),
  int sqrts=-1,
  vector<PointRedivision>* redivision=nullptr
  ){
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);
  TString strproc = ProcessName(proc);
  TString strsqrts="";
  if (sqrts>0) strsqrts = Form("_%iTeV", sqrts);

  TString strappend = Form("%s_%s_%s%s", strme.Data(), strprod.Data(), strproc.Data(), strsqrts.Data());

  TFile* finput = TFile::Open(Form("pAvgLinToLog_%s%s", strappend.Data(), ".root"), "read");
  TFile* finput_patch = TFile::Open(strpatchpath, "read");
  TFile* foutput = TFile::Open(Form("pAvgSmooth_%s%s", strappend.Data(), ".root"), "recreate");
  const unsigned int ngraphs=3;
  TString strtg[ngraphs]={
    "tg_P_ConserveDifermionMass_4mu",
    "tg_P_ConserveDifermionMass_4e",
    "tg_P_ConserveDifermionMass_2mu2e"
  };

  TGraph* tgpatch = (TGraph*)finput_patch->Get(strpatchname);
  TSpline3* sppatch = convertGraphToSpline3(tgpatch);

  const double xmin=max(0., (tgpatch->GetX())[0]);
  const double xmax=min((sqrts>0 ? (double)sqrts*1000. : 15000.), (tgpatch->GetX())[tgpatch->GetN()-1]);

  for (unsigned int ig=0; ig<ngraphs; ig++){
    finput->cd();
    TGraphErrors* tgin = 0;
    tgin = (TGraphErrors*) finput->Get(strtg[ig]);
    if (tgin==0){ cerr << strtg[ig] << " does not exist." << endl; continue; }
    TGraphErrors* tg = new TGraphErrors(*tgin);
    if (redivision){
      TTree* tree;
      if (ig==0) tree = (TTree*) finput->Get("FinalTree_4mu");
      else if (ig==1) tree = (TTree*) finput->Get("FinalTree_4e");
      else tree = (TTree*) finput->Get("FinalTree_2mu2e");
      if (!tree) cerr << "Final tree does not exist, so cannot apply redivision" << endl;
      else{
        float ZZMass, weight, mesq;
        tree->SetBranchAddress("ZZMass", &ZZMass);
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("mesq_conserveDifermMass", &mesq);
        TProfile* pmass;
        if (ig==0) pmass = (TProfile*) finput->Get("candMass_4mu");
        else if (ig==1) pmass = (TProfile*) finput->Get("candMass_4e");
        else pmass = (TProfile*) finput->Get("candMass_2mu2e");
        for (PointRedivision& rediv:*redivision){
          rebinAverageME(
            tg, tree, pmass,
            rediv,
            ZZMass, mesq, weight
          );
        }
      }
    }

    foutput->cd();
    foutput->WriteTObject(tg);

    int n = tg->GetN();
    double* xx = tg->GetX();
    double* ex = tg->GetEX();
    double* yy = tg->GetY();
    double* ey = tg->GetEY();

    TSpline3* sp = convertGraphToSpline3(tg);
    sp->SetName(Form("%s_initial", sp->GetName()));
    foutput->WriteTObject(sp);
    double tglow = xx[0];
    double tghigh = xx[tg->GetN()-1];
    cout << "Extracting high patch..." << endl;
    TGraph* highGraph = highpatcher(sp, sppatch, tghigh, xmax, false, true);
    cout << "Converting patches to splines..." << endl;
    TSpline3* highFcn = convertGraphToSpline3(highGraph);

    TF1* lowFcn = lowf(sp, xmin, tglow, true);
    lowFcn->SetNpx((int)(tglow-xmin)*5);

    vector<pair<double, double>> points;
    for (double xval=xmin; xval<tglow; xval+=1){
      double yval = lowFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }
    for (int ix=0; ix<n; ix++){
      addByLowest<double, double>(points, xx[ix], yy[ix]);
    }
    int tghigh_int = ((int)((tghigh+1.)/100.+0.5))*100;
    if (tghigh>=(double)tghigh_int) tghigh_int+=100;
    for (double xval=tghigh_int; xval<=xmax; xval+=100){
      double yval = highFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }

    int nn_new = points.size();
    cout << "Number of new points: " << nn_new-n << endl;
    double* xy_new[2];
    for (unsigned int i=0; i<2; i++) xy_new[i] = new double[nn_new];
    for (int ix=0; ix<nn_new; ix++){
      xy_new[0][ix] = points.at(ix).first;
      xy_new[1][ix] = points.at(ix).second;
    }

    delete lowFcn;
    delete highFcn;
    delete highGraph;
    delete sp;

    TGraph* tg_updated = new TGraph(nn_new, xy_new[0], xy_new[1]);
    tg_updated->SetName(Form("%s_Smooth", tg->GetName()));
    foutput->WriteTObject(tg_updated);

    sp = convertGraphToSpline3(tg_updated);
    foutput->WriteTObject(sp);

    TCanvas* ctest = new TCanvas(Form("test_%s", strtg[ig].Data()), "", 8, 30, 800, 800);
    ctest->cd();
    tg->GetXaxis()->SetRangeUser(xmin, xmax);
    tg->Draw("ae1p");
    sp->Draw("csame");
    ctest->RedrawAxis();
    ctest->Modified();
    ctest->Update();
    foutput->WriteTObject(ctest);
    ctest->SaveAs(Form("%s%s", ctest->GetName(), ".pdf"));
    ctest->Close();

    delete sp;
    delete tg_updated;
    for (unsigned int i=0; i<2; i++) delete[] xy_new[i];
    delete tg;
  }

  delete sppatch;
  foutput->Close();
  finput_patch->Close();
  finput->Close();
}

void generic_PAvgSmoothProducer_withDecay(
  TVar::MatrixElement me, TVar::Production prod, TVar::Process proc,
  TString strpatchpath, TString strpatchname,
  TGraph* (*lowpatcher)(TSpline3*, TSpline3*, double, double, bool, bool),
  TF1* (*highf)(TSpline3*, double, double, bool),
  int sqrts=-1,
  vector<PointRedivision>* redivision=nullptr
  ){
  TString strme = MatrixElementName(me);
  TString strprod = ProductionName(prod);
  TString strproc = ProcessName(proc);
  TString strsqrts="";
  if (sqrts>0) strsqrts = Form("_%iTeV", sqrts);

  TString strappend = Form("%s_%s_%s%s", strme.Data(), strprod.Data(), strproc.Data(), strsqrts.Data());

  TFile* finput = TFile::Open(Form("pAvgLinToLog_%s%s", strappend.Data(), ".root"), "read");
  TFile* finput_patch = TFile::Open(strpatchpath, "read");
  TFile* foutput = TFile::Open(Form("pAvgSmooth_%s%s", strappend.Data(), ".root"), "recreate");
  const unsigned int ngraphs=3;
  TString strtg[ngraphs]={
    "tg_P_ConserveDifermionMass_4mu",
    "tg_P_ConserveDifermionMass_4e",
    "tg_P_ConserveDifermionMass_2mu2e"
  };

  TGraph* tgpatch = (TGraph*)finput_patch->Get(strpatchname);
  TSpline3* sppatch = convertGraphToSpline3(tgpatch);

  const double xmin=max(0., (tgpatch->GetX())[0]);
  const double xmax=min((sqrts>0 ? (double)sqrts*1000. : 15000.), (tgpatch->GetX())[tgpatch->GetN()-1]);

  for (unsigned int ig=0; ig<ngraphs; ig++){
    finput->cd();
    TGraphErrors* tgin = 0;
    tgin = (TGraphErrors*) finput->Get(strtg[ig]);
    if (tgin==0){ cerr << strtg[ig] << " does not exist." << endl; continue; }
    TGraphErrors* tg = new TGraphErrors(*tgin);
    if (redivision){
      TTree* tree;
      if (ig==0) tree = (TTree*) finput->Get("FinalTree_4mu");
      else if (ig==1) tree = (TTree*) finput->Get("FinalTree_4e");
      else tree = (TTree*) finput->Get("FinalTree_2mu2e");
      if (!tree) cerr << "Final tree does not exist, so cannot apply redivision" << endl;
      else{
        float ZZMass, weight, mesq;
        tree->SetBranchAddress("ZZMass", &ZZMass);
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("mesq_conserveDifermMass", &mesq);
        TProfile* pmass;
        if (ig==0) pmass = (TProfile*) finput->Get("candMass_4mu");
        else if (ig==1) pmass = (TProfile*) finput->Get("candMass_4e");
        else pmass = (TProfile*) finput->Get("candMass_2mu2e");
        for (PointRedivision& rediv:*redivision){
          rebinAverageME(
            tg, tree, pmass,
            rediv,
            ZZMass, mesq, weight
          );
        }
      }
    }

    foutput->cd();
    foutput->WriteTObject(tg);

    int n = tg->GetN();
    double* xx = tg->GetX();
    double* ex = tg->GetEX();
    double* yy = tg->GetY();
    double* ey = tg->GetEY();

    TSpline3* sp = convertGraphToSpline3(tg);
    sp->SetName(Form("%s_initial", sp->GetName()));
    foutput->WriteTObject(sp);
    double tglow = xx[0];
    double tghigh = xx[tg->GetN()-1];
    cout << "Extracting low patch..." << endl;
    TGraph* lowGraph = lowpatcher(sp, sppatch, tglow, xmax, true, true);
    cout << "Converting patches to splines..." << endl;
    TSpline3* lowFcn = convertGraphToSpline3(lowGraph);

    TF1* highFcn = highf(sp, tghigh, xmax, false);
    highFcn->SetNpx((int)(tglow-xmin)*5);

    vector<pair<double, double>> points;
    for (double xval=xmin; xval<tglow; xval+=1){
      double yval = lowFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }
    for (int ix=0; ix<n; ix++){
      addByLowest<double, double>(points, xx[ix], yy[ix]);
    }
    int tghigh_int = ((int)((tghigh+1.)/100.+0.5))*100;
    if (tghigh>=(double)tghigh_int) tghigh_int+=100;
    for (double xval=tghigh_int; xval<=xmax; xval+=100){
      double yval = highFcn->Eval(xval);
      addByLowest<double, double>(points, xval, yval);
    }

    int nn_new = points.size();
    cout << "Number of new points: " << nn_new-n << endl;
    double* xy_new[2];
    for (unsigned int i=0; i<2; i++) xy_new[i] = new double[nn_new];
    for (int ix=0; ix<nn_new; ix++){
      xy_new[0][ix] = points.at(ix).first;
      xy_new[1][ix] = points.at(ix).second;
    }

    delete lowFcn;
    delete highFcn;
    delete lowGraph;
    delete sp;

    TGraph* tg_updated = new TGraph(nn_new, xy_new[0], xy_new[1]);
    tg_updated->SetName(Form("%s_Smooth", tg->GetName()));
    foutput->WriteTObject(tg_updated);

    sp = convertGraphToSpline3(tg_updated);
    foutput->WriteTObject(sp);

    TCanvas* ctest = new TCanvas(Form("test_%s", strtg[ig].Data()), "", 8, 30, 800, 800);
    ctest->cd();
    tg->GetXaxis()->SetRangeUser(xmin, xmax);
    tg->Draw("ae1p");
    sp->Draw("csame");
    ctest->RedrawAxis();
    ctest->Modified();
    ctest->Update();
    foutput->WriteTObject(ctest);
    ctest->SaveAs(Form("%s%s", ctest->GetName(), ".pdf"));
    ctest->Close();

    delete sp;
    delete tg_updated;
    for (unsigned int i=0; i<2; i++) delete[] xy_new[i];
    delete tg;
  }

  delete sppatch;
  foutput->Close();
  finput_patch->Close();
  finput->Close();
}

/* SPECIFIC COMMENT: NONE */
void produce_PAvgSmooth_MCFM_JJQCD_bkgZJets_2l2q(int sqrts=13){
  TFile* finput = new TFile(Form("pAvgLinToLog_MCFM_JJQCD_bkgZJets_%iTeV_2l2q.root", sqrts), "read");
  TFile* foutput = new TFile(Form("pAvgSmooth_MCFM_JJQCD_bkgZJets_%iTeV_2l2q.root", sqrts), "recreate");
  const unsigned int ngraphs=1;
  TString strtg[ngraphs]={
    "tg_P_ConserveDifermionMass"
  };

  for (unsigned int ig=0; ig<ngraphs; ig++){
    finput->cd();
    TGraphErrors* tg = 0;
    tg = (TGraphErrors*)finput->Get(strtg[ig].Data());
    if (tg==0){ cerr << strtg[ig] << " does not exist." << endl; continue; }
    tg->SetName(strtg[ig]);
    foutput->cd();
    foutput->WriteTObject(tg);
    tg->SetName(Form("%s_Smooth", tg->GetName()));

    TGraphErrors* tg_new = replacePointsBetween(tg, 510, 520); tg=tg_new;
    tg_new = replacePointsBetween(tg, 329.5, 330.5); delete tg; tg=tg_new;
    tg_new = replacePointsBetween(tg, 338, 340); delete tg; tg=tg_new;

    regularizeSlice(tg);
    foutput->WriteTObject(tg);

    TSpline3* sp = convertGraphToSpline3(tg);
    TF1* lowFcn = getFcn_a0plusa1overX(sp, 0, (tg->GetX())[0], true);
    TF1* highFcn = getFcn_a0plusa1timesX(sp, (tg->GetX())[tg->GetN()-1], 20000., false);
    lowFcn->SetNpx(1000);
    highFcn->SetNpx(10000);

    foutput->WriteTObject(sp);
    foutput->WriteTObject(lowFcn);
    foutput->WriteTObject(highFcn);

    TCanvas* ctest = new TCanvas("test", "", 8, 30, 800, 800);
    ctest->cd();
    tg->GetXaxis()->SetRangeUser(0, 20000);
    tg->Draw("ae1p");
    sp->Draw("csame");
    lowFcn->Draw("csame");
    highFcn->Draw("csame");
    ctest->RedrawAxis();
    ctest->Modified();
    ctest->Update();
    foutput->WriteTObject(ctest);
    ctest->Close();

    delete highFcn;
    delete lowFcn;
    delete sp;
    delete tg;
  }
  foutput->Close();
  finput->Close();
}

/* SPECIFIC COMMENT: NONE */
void produce_PAvgSmooth_JHUGen_JJVBF_HSMHiggs(int sqrts=8){
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::JJVBF;
  TVar::Process proc = TVar::HSMHiggs;
  generic_PAvgSmoothProducer(
    me, prod, proc,
    &getFcn_a0plusa1timesX,
    &getFcn_a0plusa1timesX,
    sqrts
    );
}

/* SPECIFIC COMMENT: NONE */
void produce_PAvgSmooth_JHUGen_JJQCD_HSMHiggs(int sqrts=8){
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::JJQCD;
  TVar::Process proc = TVar::HSMHiggs;
  generic_PAvgSmoothProducer(
    me, prod, proc,
    &getFcn_a0plusa1timesX,
    &getFcn_a0plusa1overX,
    sqrts
    );
}

/* SPECIFIC COMMENT: NONE */
void produce_PAvgSmooth_JHUGen_JQCD_HSMHiggs(int sqrts=8){
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::JQCD;
  TVar::Process proc = TVar::HSMHiggs;
  generic_PAvgSmoothProducer(
    me, prod, proc,
    &getFcn_a0plusa1timesX,
    &getFcn_a0plusa1overX,
    sqrts
    );
}

void produce_PAvgSmooth_JHUGen_HadVH_HSMHiggs(TString strprod, int sqrts=13){
  if (!(strprod == "Had_ZH" || strprod == "Had_WH")) return;

  TVar::Production prod;
  for (int iprod=(int)TVar::JJVBF; iprod<(int)TVar::nProductions; iprod++){
    prod = (TVar::Production)iprod;
    if (TVar::ProductionName(prod)==strprod) break;
  }
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Process proc = TVar::HSMHiggs;
  vector<PointRedivision> redivs;
  redivs.emplace_back(2, 90, 90);
  redivs.emplace_back(2, 90, 90);
  redivs.emplace_back(2, 90, 90);
  redivs.emplace_back(4, 1600, 3000);
  generic_PAvgSmoothProducer(
    me, prod, proc,
    &getFcn_a0plusa1expX,
    &getFcn_a0plusa1overX,
    sqrts,
    false, false,
    &redivs
    );
}

/* SPECIFIC COMMENT: PATCHING DONE USING ANALYTICAL ggH PDF */
void produce_PAvgSmooth_JHUGen_ZZGG_HSMHiggs(){
  TVar::MatrixElement me = TVar::JHUGen;
  TVar::Production prod = TVar::ZZGG;
  TVar::Process proc = TVar::HSMHiggs;
  generic_PAvgSmoothProducer_withDecay(
    me, prod, proc,
    "../data/pAvgLinToLog_ANALYTICAL_ZZGG_HSMHiggs.root", "tg_anaPdfInt",
    &getPatch_a0plusa1timesX,
    &getPatch_a0plusa1overX,
    -1
    );
}

/* SPECIFIC COMMENT: PATCHING DONE USING ANALYTICAL ggH PDF */
void produce_PAvgSmooth_MCFM_ZZGG_HSMHiggs(){
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod = TVar::ZZGG;
  TVar::Process proc = TVar::HSMHiggs;
  generic_PAvgSmoothProducer_withDecay(
    me, prod, proc,
    "../data/pAvgLinToLog_ANALYTICAL_ZZGG_HSMHiggs.root", "tg_anaPdfInt",
    &getPatch_a0plusa1timesX,
    &getPatch_a0plusa1overX,
    -1
    );
}

/* SPECIFIC COMMENT: PATCHING DONE USING ANALYTICAL QQZZ PDF */
void produce_PAvgSmooth_MCFM_ZZGG_bkgZZ(){
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod = TVar::ZZGG;
  TVar::Process proc = TVar::bkgZZ;
  generic_PAvgSmoothProducer_withDecay(
    me, prod, proc,
    "../data/pAvgLinToLog_ANALYTICAL_ZZGG_HSMHiggs.root", "tg_anaPdfInt",
    "../data/pAvgLinToLog_ANALYTICAL_ZZQQB_bkgZZ.root", "tg_anaPdfInt",
    &getPatch_a0plusa1timesX,
    &getPatch_a0plusa1overX,
    -1
    );
}

/* SPECIFIC COMMENT: PATCHING DONE USING ANALYTICAL QQZZ PDF */
void produce_PAvgSmooth_MCFM_ZZQQB_bkgZZ(){
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod = TVar::ZZQQB;
  TVar::Process proc = TVar::bkgZZ;
  generic_PAvgSmoothProducer_withDecay(
    me, prod, proc,
    "../data/pAvgLinToLog_ANALYTICAL_ZZQQB_bkgZZ.root","tg_anaPdfInt",
    &getPatch_a0plusa1timesX,
    &getPatch_a0plusa1overX,
    -1
    );
}

/* SPECIFIC COMMENT: PATCHING DONE USING ANALYTICAL ggH PDF */
void produce_PAvgSmooth_MCFM_JJPROD_S_HSMHiggs(TString strprod, int sqrts=13){
  if (!(strprod == "Had_ZH_S" || strprod == "Had_WH_S" || strprod == "JJVBF_S")) return;

  TVar::Production prod;
  for (int iprod=(int)TVar::JJVBF; iprod<(int)TVar::nProductions; iprod++){
    prod = (TVar::Production)iprod;
    if (TVar::ProductionName(prod)==strprod) break;
  }
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Process proc = TVar::HSMHiggs;
  //generic_PAvgSmoothProducer_withDecay(
  //  me, prod, proc,
  //  "../data/pAvgLinToLog_ANALYTICAL_ZZGG_HSMHiggs.root", "tg_anaPdfInt",
  //  &getPatch_a0plusa1timesX,
  //  &getPatch_a0plusa1overX,
  //  sqrts
  //  );
  vector<PointRedivision> redivs;
  if (prod==TVar::Had_ZH_S || prod==TVar::Had_WH_S){
    redivs.emplace_back(2, 90, 90);
    redivs.emplace_back(2, 90, 90);
    redivs.emplace_back(2, 90, 90);

    redivs.emplace_back(4, 1600, 3000);

    vector<float> redivmultiplesA{ 192.072, 210.622, 243.901, 272.591, 331.305, 392.596, 449.076, 504.783, 565.813, 640.583, 729.157, 822.262, 945.274 };
    vector<float> redivmultiplesB{ 195.798, 215.901, 249.163, 286.182, 346.839, 406.022, 464.092, 522.738, 588.661, 672.415, 764.014, 862.423, 1015.89 };
    for (unsigned int irm=0; irm<redivmultiplesA.size()-1; irm++){
      if (prod==TVar::Had_ZH_S && (irm==0 || irm==2)) continue;
      float xlow=std::max(redivmultiplesA.at(irm), redivmultiplesB.at(irm))+1e-2;
      float xhigh=std::min(redivmultiplesA.at(irm+1), redivmultiplesB.at(irm+1))-1e-2;
      redivs.emplace_back(1, xlow, xhigh);
    }
  }
  generic_PAvgSmoothProducer_withDecay(
    me, prod, proc,
    &getFcn_a0plusa1timesX,
    &getFcn_a0plusa1overX,
    sqrts,
    true, false,
    &redivs
  );
}

/* SPECIFIC COMMENT: PATCHING DONE USING FUNCTIONS */
void produce_PAvgSmooth_MCFM_JJPROD_bkgZZ(TString strprod, int sqrts=13){
  if (!(strprod == "Had_ZH" || strprod == "Had_WH" || strprod == "JJVBF")) return;

  TVar::Production prod;
  for (int iprod=(int) TVar::JJVBF; iprod<(int) TVar::nProductions; iprod++){
    prod = (TVar::Production)iprod;
    if (TVar::ProductionName(prod)==strprod) break;
  }
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Process proc = TVar::bkgZZ;
  if (prod!=TVar::JJVBF){
    vector<PointRedivision> redivs;
    redivs.emplace_back(2, 90, 90);
    redivs.emplace_back(2, 90, 90);
    redivs.emplace_back(2, 90, 90);
    if (prod==TVar::Had_ZH) redivs.emplace_back(5, 125.3, 170);
    else{
      redivs.emplace_back(1, 105, 118);
      redivs.emplace_back(1, 120, 125);
    }
    generic_PAvgSmoothProducer_withDecay(
      me, prod, proc,
      //"../data/pAvgLinToLog_ANALYTICAL_ZZQQB_bkgZZ.root", "tg_anaPdfInt",
      &getFcn_a0plusa1timesX,
      &getFcn_a0plusa1overX,
      sqrts,
      true, false,
      &redivs
    );
  }
  else{
    generic_PAvgSmoothProducer_withDecay(
      me, prod, proc,
      &getFcn_a0plusa1timesX,
      &getFcn_a0plusa1timesX,
      sqrts,
      false, false
    );
  }
}

/* SPECIFIC COMMENT: PATCHING DONE USING ANALYTICAL QQZZ PDF */
void produce_PAvgSmooth_MCFM_JJQCD_bkgZZ(int sqrts=13){
  TVar::MatrixElement me = TVar::MCFM;
  TVar::Production prod = TVar::JJQCD;
  TVar::Process proc = TVar::bkgZZ;
  generic_PAvgSmoothProducer_withDecay(
    me, prod, proc,
    "../data/pAvgLinToLog_ANALYTICAL_ZZQQB_bkgZZ.root", "tg_anaPdfInt",
    &getFcn_a0plusa1overX,
    &getPatch_a0plusa1overX,
    sqrts
    );
}


void check_JJVBF_vs_JJQCD_7or8TeV(int sqrts=8){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

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
  float jetptetaphimass[2][4];

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  int LepID[4]={ 13, -13, 11, -11 };

  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TString cinput_main;
  if (sqrts==8) cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  else if (sqrts==7) cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR";
  else return;
  const int nSamples_JJVBF = 42;
  TString strSamples_JJVBF[nSamples_JJVBF]={
    "HZZ4lTree_VBFH116.root",
    "HZZ4lTree_VBFH117.root",
    "HZZ4lTree_VBFH118.root",
    "HZZ4lTree_VBFH119.root",
    "HZZ4lTree_VBFH120.root",
    "HZZ4lTree_VBFH121.root",
    "HZZ4lTree_VBFH122.root",
    "HZZ4lTree_VBFH123.root",
    "HZZ4lTree_VBFH124.root",
    "HZZ4lTree_VBFH125.root",
    "HZZ4lTree_VBFH126.root",
    "HZZ4lTree_VBFH127.root",
    "HZZ4lTree_VBFH128.root",
    "HZZ4lTree_VBFH129.root",
    "HZZ4lTree_VBFH130.root",
    "HZZ4lTree_VBFH135.root",
    "HZZ4lTree_VBFH140.root",
    "HZZ4lTree_VBFH145.root",
    "HZZ4lTree_VBFH150.root",
    "HZZ4lTree_VBFH160.root",
    "HZZ4lTree_VBFH170.root",
    "HZZ4lTree_VBFH180.root",
    "HZZ4lTree_VBFH190.root",
    "HZZ4lTree_powheg15VBFH200.root",
    "HZZ4lTree_powheg15VBFH225.root",
    "HZZ4lTree_powheg15VBFH250.root",
    "HZZ4lTree_powheg15VBFH275.root",
    "HZZ4lTree_powheg15VBFH300.root",
    "HZZ4lTree_powheg15VBFH350.root",
    "HZZ4lTree_powheg15VBFH400.root",
    "HZZ4lTree_powheg15VBFH450.root",
    "HZZ4lTree_powheg15VBFH500.root",
    "HZZ4lTree_powheg15VBFH550.root",
    "HZZ4lTree_powheg15VBFH600.root",
    "HZZ4lTree_powheg15VBFH650.root",
    "HZZ4lTree_powheg15VBFH700.root",
    "HZZ4lTree_powheg15VBFH750.root",
    "HZZ4lTree_powheg15VBFH800.root",
    "HZZ4lTree_powheg15VBFH850.root",
    "HZZ4lTree_powheg15VBFH900.root",
    "HZZ4lTree_powheg15VBFH950.root",
    "HZZ4lTree_powheg15VBFH1000.root"
  };
  const int nSamples_JJQCD = 37;
  TString strSamples_JJQCD[nSamples_JJQCD]={
    "HZZ4lTree_minloH90.root",
    "HZZ4lTree_minloH95.root",
    "HZZ4lTree_minloH100.root",
    "HZZ4lTree_minloH105.root",
    "HZZ4lTree_minloH110.root",
    "HZZ4lTree_minloH115.root",
    "HZZ4lTree_minloH120.root",
    "HZZ4lTree_minloH124.root",
    "HZZ4lTree_minloH125.root",
    "HZZ4lTree_minloH126.root",
    "HZZ4lTree_minloH130.root",
    "HZZ4lTree_minloH135.root",
    "HZZ4lTree_minloH140.root",
    "HZZ4lTree_minloH145.root",
    "HZZ4lTree_minloH150.root",
    "HZZ4lTree_minloH155.root",
    "HZZ4lTree_minloH160.root",
    "HZZ4lTree_minloH170.root",
    "HZZ4lTree_minloH180.root",
    "HZZ4lTree_minloH190.root",
    "HZZ4lTree_minloH200.root",
    "HZZ4lTree_minloH250.root",
    "HZZ4lTree_minloH300.root",
    "HZZ4lTree_minloH350.root",
    "HZZ4lTree_minloH400.root",
    "HZZ4lTree_minloH450.root",
    "HZZ4lTree_minloH500.root",
    "HZZ4lTree_minloH550.root",
    "HZZ4lTree_minloH600.root",
    "HZZ4lTree_minloH650.root",
    "HZZ4lTree_minloH700.root",
    "HZZ4lTree_minloH750.root",
    "HZZ4lTree_minloH800.root",
    "HZZ4lTree_minloH850.root",
    "HZZ4lTree_minloH900.root",
    "HZZ4lTree_minloH950.root",
    "HZZ4lTree_minloH1000.root"
  };

  TChain* tree[2] ={
    new TChain(TREE_NAME, ""),
    new TChain(TREE_NAME, "")
  };
  for (int ic=0; ic<3; ic++){
    for (int is=0; is<nSamples_JJVBF; is++) tree[0]->Add(Form("%s/%s/%s", cinput_main.Data(), (strchannel[ic]).Data(), (strSamples_JJVBF[is]).Data()));
    for (int is=0; is<nSamples_JJQCD; is++) tree[1]->Add(Form("%s/%s/%s", cinput_main.Data(), (strchannel[ic]).Data(), (strSamples_JJQCD[is]).Data()));
  }
  for (int it=0; it<2; it++){
    tree[it]->SetBranchAddress("NJets30", &NJets30);
    tree[it]->SetBranchAddress("JetPt", &JetPt, &bJetPt);
    tree[it]->SetBranchAddress("JetEta", &JetEta, &bJetEta);
    tree[it]->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
    tree[it]->SetBranchAddress("JetMass", &JetMass, &bJetMass);
    tree[it]->SetBranchAddress("ZZMass", &mzz);
    tree[it]->SetBranchAddress("ZZPt", &ZZPt);
    tree[it]->SetBranchAddress("ZZEta", &ZZEta);
    tree[it]->SetBranchAddress("ZZPhi", &ZZPhi);
    tree[it]->SetBranchAddress("Z1Mass", &m1);
    tree[it]->SetBranchAddress("Z2Mass", &m2);
    tree[it]->SetBranchAddress("helcosthetaZ1", &h1);
    tree[it]->SetBranchAddress("helcosthetaZ2", &h2);
    tree[it]->SetBranchAddress("helphi", &phi);
    tree[it]->SetBranchAddress("costhetastar", &hs);
    tree[it]->SetBranchAddress("phistarZ1", &phi1);
  }

  TFile* foutput = new TFile(Form("pJHUGen_JJVBF_JJQCD_HSMHiggs_Comparison_%iTeV.root", sqrts), "recreate");

  TH2F* hJJVBF = new TH2F("JJVBF", "", 50, 70, 1070, 10, 0, 1);
  TH2F* hJJQCD = new TH2F("JJQCD", "", 50, 70, 1070, 10, 0, 1);

  TProfile* prJJVBF = new TProfile("prJJVBF", "", 50, 70, 1070); prJJVBF->Sumw2();
  TProfile* prJJQCD = new TProfile("prJJQCD", "", 50, 70, 1070); prJJQCD->Sumw2();

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  int nTotalEntries = tree[0]->GetEntries();
  for (int ev = 0; ev < nTotalEntries; ev++){
    tree[0]->GetEntry(ev);
    if (ev%10000==0) cout << "Doing event " << ev << endl;
    if (NJets30>=2){
      for (int ij=0; ij<2; ij++){
        jetptetaphimass[ij][0]=JetPt->at(ij);
        jetptetaphimass[ij][1]=JetEta->at(ij);
        jetptetaphimass[ij][2]=JetPhi->at(ij);
        jetptetaphimass[ij][3]=JetMass->at(ij);
      }

      TLorentzVector jet[2], higgs;
      for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(jetptetaphimass[ij][0], jetptetaphimass[ij][1], jetptetaphimass[ij][2], jetptetaphimass[ij][3]);
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));
      associated.push_back(SimpleParticle_t(0, jet[1]));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      float mesq_jjvbf=0, mesq_jjqcd=0;
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_jjvbf, true);
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_jjqcd, true);
      float kd = mesq_jjvbf/(mesq_jjvbf+0.06*mesq_jjqcd);

      hJJVBF->Fill(mzz, kd);
      prJJVBF->Fill(mzz, mesq_jjvbf);

      mela.resetInputEvent();
    }
  }
  nTotalEntries = tree[1]->GetEntries();
  for (int ev = 0; ev < nTotalEntries; ev++){
    tree[1]->GetEntry(ev);
    if (ev%10000==0) cout << "Doing event " << ev << endl;
    if (NJets30>=2){
      for (int ij=0; ij<2; ij++){
        jetptetaphimass[ij][0]=JetPt->at(ij);
        jetptetaphimass[ij][1]=JetEta->at(ij);
        jetptetaphimass[ij][2]=JetPhi->at(ij);
        jetptetaphimass[ij][3]=JetMass->at(ij);
      }

      TLorentzVector jet[2], higgs;
      for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(jetptetaphimass[ij][0], jetptetaphimass[ij][1], jetptetaphimass[ij][2], jetptetaphimass[ij][3]);
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));
      associated.push_back(SimpleParticle_t(0, jet[1]));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      float mesq_jjvbf=0, mesq_jjqcd=0;
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_jjvbf, true);
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_jjqcd, true);
      float kd = mesq_jjvbf/(mesq_jjvbf+0.06*mesq_jjqcd);

      hJJQCD->Fill(mzz, kd);
      prJJQCD->Fill(mzz, mesq_jjqcd);

      mela.resetInputEvent();
    }
  }

  for (int ix=1; ix<=hJJVBF->GetNbinsX(); ix++){
    double integral = hJJVBF->Integral(ix, ix, 0, hJJVBF->GetNbinsY()+1);
    for (int iy=0; iy<=hJJVBF->GetNbinsY(); iy++){ if (integral!=0.) hJJVBF->SetBinContent(ix, iy, hJJVBF->GetBinContent(ix, iy)/integral); }
  }
  for (int ix=1; ix<=hJJQCD->GetNbinsX(); ix++){
    double integral = hJJQCD->Integral(ix, ix, 0, hJJQCD->GetNbinsY()+1);
    for (int iy=0; iy<=hJJQCD->GetNbinsY(); iy++){ if (integral!=0.) hJJQCD->SetBinContent(ix, iy, hJJQCD->GetBinContent(ix, iy)/integral); }
  }

  foutput->WriteTObject(hJJVBF);
  foutput->WriteTObject(hJJQCD);
  foutput->WriteTObject(prJJVBF);
  foutput->WriteTObject(prJJQCD);
  foutput->Close();
  for (int it=0; it<2; it++) delete tree[it];
}

void check_JJVBF_vs_JJQCD_13TeV(int sqrts=13){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  short NJets30;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float> myJetPt;
  std::vector<float> myJetEta;
  std::vector<float> myJetPhi;
  std::vector<float> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;
  float jetptetaphimass[2][4];

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  int LepID[4]={ 13, -13, 11, -11 };

  TString cinput_main;
  if (sqrts==13) cinput_main = "/scratch0/hep/usarical/CJLST/LHC_13TeV/4l/160225";
  else return;
  const int nSamples_JJVBF = 33;
  TString strSamples_JJVBF[nSamples_JJVBF]={
    "VBFH115/ZZ4lAnalysis.root",
    "VBFH120/ZZ4lAnalysis.root",
    "VBFH124/ZZ4lAnalysis.root",
    "VBFH125/ZZ4lAnalysis.root",
    "VBFH126/ZZ4lAnalysis.root",
    "VBFH130/ZZ4lAnalysis.root",
    "VBFH135/ZZ4lAnalysis.root",
    "VBFH140/ZZ4lAnalysis.root",
    "VBFH145/ZZ4lAnalysis.root",
    "VBFH150/ZZ4lAnalysis.root",
    "VBFH155/ZZ4lAnalysis.root",
    "VBFH160/ZZ4lAnalysis.root",
    "VBFH165/ZZ4lAnalysis.root",
    "VBFH170/ZZ4lAnalysis.root",
    "VBFH175/ZZ4lAnalysis.root",
    "VBFH180/ZZ4lAnalysis.root",
    "VBFH190/ZZ4lAnalysis.root",
    "VBFH200/ZZ4lAnalysis.root",
    "VBFH210/ZZ4lAnalysis.root",
    "VBFH230/ZZ4lAnalysis.root",
    "VBFH250/ZZ4lAnalysis.root",
    "VBFH270/ZZ4lAnalysis.root",
    "VBFH300/ZZ4lAnalysis.root",
    "VBFH350/ZZ4lAnalysis.root",
    "VBFH400/ZZ4lAnalysis.root",
    "VBFH450/ZZ4lAnalysis.root",
    "VBFH500/ZZ4lAnalysis.root",
    "VBFH550/ZZ4lAnalysis.root",
    "VBFH600/ZZ4lAnalysis.root",
    "VBFH700/ZZ4lAnalysis.root",
    "VBFH750/ZZ4lAnalysis.root",
    "VBFH800/ZZ4lAnalysis.root",
    "VBFH900/ZZ4lAnalysis.root"
  };
  const int nSamples_JJQCD = 35;
  TString strSamples_JJQCD[nSamples_JJQCD]={
    "ggH91_GaZ/ZZ4lAnalysis.root",
    "ggH115/ZZ4lAnalysis.root",
    "ggH120/ZZ4lAnalysis.root",
    "ggH124/ZZ4lAnalysis.root",
    "ggH125/ZZ4lAnalysis.root",
    "ggH126/ZZ4lAnalysis.root",
    "ggH130/ZZ4lAnalysis.root",
    "ggH135/ZZ4lAnalysis.root",
    "ggH140/ZZ4lAnalysis.root",
    "ggH145/ZZ4lAnalysis.root",
    "ggH150/ZZ4lAnalysis.root",
    "ggH155/ZZ4lAnalysis.root",
    "ggH160/ZZ4lAnalysis.root",
    "ggH165/ZZ4lAnalysis.root",
    "ggH170/ZZ4lAnalysis.root",
    "ggH175/ZZ4lAnalysis.root",
    "ggH180/ZZ4lAnalysis.root",
    "ggH190/ZZ4lAnalysis.root",
    "ggH200/ZZ4lAnalysis.root",
    "ggH210/ZZ4lAnalysis.root",
    "ggH230/ZZ4lAnalysis.root",
    "ggH250/ZZ4lAnalysis.root",
    "ggH270/ZZ4lAnalysis.root",
    "ggH300/ZZ4lAnalysis.root",
    "ggH350/ZZ4lAnalysis.root",
    "ggH400/ZZ4lAnalysis.root",
    "ggH450/ZZ4lAnalysis.root",
    "ggH500/ZZ4lAnalysis.root",
    "ggH550/ZZ4lAnalysis.root",
    "ggH600/ZZ4lAnalysis.root",
    "ggH700/ZZ4lAnalysis.root",
    "ggH750/ZZ4lAnalysis.root",
    "ggH800/ZZ4lAnalysis.root",
    "ggH900/ZZ4lAnalysis.root",
    "ggH1000/ZZ4lAnalysis.root"
  };

  TChain* tree[2] ={
    new TChain(TREE_NAME, ""),
    new TChain(TREE_NAME, "")
  };
  for (int is=0; is<nSamples_JJVBF; is++) tree[0]->Add(Form("%s/%s", cinput_main.Data(), (strSamples_JJVBF[is]).Data()));
  for (int is=0; is<nSamples_JJQCD; is++) tree[1]->Add(Form("%s/%s", cinput_main.Data(), (strSamples_JJQCD[is]).Data()));

  for (int it=0; it<2; it++){
    tree[it]->SetBranchAddress("nCleanedJetsPt30", &NJets30);
    tree[it]->SetBranchAddress("JetPt", &JetPt, &bJetPt);
    tree[it]->SetBranchAddress("JetEta", &JetEta, &bJetEta);
    tree[it]->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
    tree[it]->SetBranchAddress("JetMass", &JetMass, &bJetMass);
    tree[it]->SetBranchAddress("ZZMass", &mzz);
    tree[it]->SetBranchAddress("ZZPt", &ZZPt);
    tree[it]->SetBranchAddress("ZZEta", &ZZEta);
    tree[it]->SetBranchAddress("ZZPhi", &ZZPhi);
    tree[it]->SetBranchAddress("Z1Mass", &m1);
    tree[it]->SetBranchAddress("Z2Mass", &m2);
    tree[it]->SetBranchAddress("helcosthetaZ1", &h1);
    tree[it]->SetBranchAddress("helcosthetaZ2", &h2);
    tree[it]->SetBranchAddress("helphi", &phi);
    tree[it]->SetBranchAddress("costhetastar", &hs);
    tree[it]->SetBranchAddress("phistarZ1", &phi1);
  }

  TFile* foutput = new TFile(Form("pJHUGen_JJVBF_JJQCD_HSMHiggs_Comparison_%iTeV.root", sqrts), "recreate");

  TH2F* hJJVBF = new TH2F("JJVBF", "", 75, 70, 1570, 10, 0, 1);
  TH2F* hJJQCD = new TH2F("JJQCD", "", 75, 70, 1570, 10, 0, 1);

  TProfile* prJJVBF = new TProfile("prJJVBF", "", 75, 70, 1570); prJJVBF->Sumw2();
  TProfile* prJJQCD = new TProfile("prJJQCD", "", 75, 70, 1570); prJJQCD->Sumw2();

  TProfile* prRatioForJJVBF = new TProfile("prRatioForJJVBF", "", 75, 70, 1570); prRatioForJJVBF->Sumw2();
  TProfile* prRatioForJJQCD = new TProfile("prRatioForJJQCD", "", 75, 70, 1570); prRatioForJJQCD->Sumw2();

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  int nTotalEntries = tree[0]->GetEntries();
  for (int ev = 0; ev < nTotalEntries; ev++){
    tree[0]->GetEntry(ev);
    if (ev%10000==0) cout << "Doing event " << ev << endl;
    if (NJets30>=2){
      for (int ij=0; ij<2; ij++){
        jetptetaphimass[ij][0]=JetPt->at(ij);
        jetptetaphimass[ij][1]=JetEta->at(ij);
        jetptetaphimass[ij][2]=JetPhi->at(ij);
        jetptetaphimass[ij][3]=JetMass->at(ij);
      }

      TLorentzVector jet[2], higgs;
      for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(jetptetaphimass[ij][0], jetptetaphimass[ij][1], jetptetaphimass[ij][2], jetptetaphimass[ij][3]);
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));
      associated.push_back(SimpleParticle_t(0, jet[1]));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      float mesq_jjvbf=0, mesq_jjqcd=0;
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_jjvbf, true);
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_jjqcd, true);
      float kd = mesq_jjvbf/(mesq_jjvbf+0.06*mesq_jjqcd);

      hJJVBF->Fill(mzz, kd);
      prJJVBF->Fill(mzz, mesq_jjvbf);
      if (mesq_jjvbf>0.) prRatioForJJVBF->Fill(mzz, mesq_jjqcd/mesq_jjvbf);

      mela.resetInputEvent();
    }
  }
  nTotalEntries = tree[1]->GetEntries();
  for (int ev = 0; ev < nTotalEntries; ev++){
    tree[1]->GetEntry(ev);
    if (ev%10000==0) cout << "Doing event " << ev << endl;
    if (NJets30>=2){
      for (int ij=0; ij<2; ij++){
        jetptetaphimass[ij][0]=JetPt->at(ij);
        jetptetaphimass[ij][1]=JetEta->at(ij);
        jetptetaphimass[ij][2]=JetPhi->at(ij);
        jetptetaphimass[ij][3]=JetMass->at(ij);
      }

      TLorentzVector jet[2], higgs;
      for (int ij=0; ij<2; ij++) jet[ij].SetPtEtaPhiM(jetptetaphimass[ij][0], jetptetaphimass[ij][1], jetptetaphimass[ij][2], jetptetaphimass[ij][3]);
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));
      associated.push_back(SimpleParticle_t(0, jet[1]));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      float mesq_jjvbf=0, mesq_jjqcd=0;
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_jjvbf, true);
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
      TUtil::setJetMassScheme(TVar::ConserveDifermionMass);
      mela.computeProdP(mesq_jjqcd, true);
      float kd = mesq_jjvbf/(mesq_jjvbf+0.06*mesq_jjqcd);

      hJJQCD->Fill(mzz, kd);
      prJJQCD->Fill(mzz, mesq_jjqcd);
      if (mesq_jjvbf>0.) prRatioForJJQCD->Fill(mzz, mesq_jjqcd/mesq_jjvbf);

      mela.resetInputEvent();
    }
  }

  for (int ix=1; ix<=hJJVBF->GetNbinsX(); ix++){
    double integral = hJJVBF->Integral(ix, ix, 0, hJJVBF->GetNbinsY()+1);
    for (int iy=0; iy<=hJJVBF->GetNbinsY(); iy++){ if (integral!=0.) hJJVBF->SetBinContent(ix, iy, hJJVBF->GetBinContent(ix, iy)/integral); }
  }
  for (int ix=1; ix<=hJJQCD->GetNbinsX(); ix++){
    double integral = hJJQCD->Integral(ix, ix, 0, hJJQCD->GetNbinsY()+1);
    for (int iy=0; iy<=hJJQCD->GetNbinsY(); iy++){ if (integral!=0.) hJJQCD->SetBinContent(ix, iy, hJJQCD->GetBinContent(ix, iy)/integral); }
  }

  foutput->WriteTObject(hJJVBF);
  foutput->WriteTObject(hJJQCD);
  foutput->WriteTObject(prJJVBF);
  foutput->WriteTObject(prJJQCD);
  foutput->WriteTObject(prRatioForJJVBF);
  foutput->WriteTObject(prRatioForJJQCD);
  foutput->Close();
  for (int it=0; it<2; it++) delete tree[it];
}


void check_JQCD_7or8TeV(int sqrts=8){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "SelectedTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

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
  float jetptetaphimass[2][4];

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  int LepID[4]={ 13, -13, 11, -11 };

  const int nMEs=1;
  TChain* tree[nMEs] ={
    new TChain(TREE_NAME, "")
  };

  TString strchannel[3]={ "4mu", "4e", "2mu2e" };
  TString cinput_main;
  if (sqrts==8) cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  else if (sqrts==7) cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR";
  else return;
  const int nSamples_JQCD = 37;
  TString strSamples_JQCD[nSamples_JQCD]={
    "HZZ4lTree_minloH90.root",
    "HZZ4lTree_minloH95.root",
    "HZZ4lTree_minloH100.root",
    "HZZ4lTree_minloH105.root",
    "HZZ4lTree_minloH110.root",
    "HZZ4lTree_minloH115.root",
    "HZZ4lTree_minloH120.root",
    "HZZ4lTree_minloH124.root",
    "HZZ4lTree_minloH125.root",
    "HZZ4lTree_minloH126.root",
    "HZZ4lTree_minloH130.root",
    "HZZ4lTree_minloH135.root",
    "HZZ4lTree_minloH140.root",
    "HZZ4lTree_minloH145.root",
    "HZZ4lTree_minloH150.root",
    "HZZ4lTree_minloH155.root",
    "HZZ4lTree_minloH160.root",
    "HZZ4lTree_minloH170.root",
    "HZZ4lTree_minloH180.root",
    "HZZ4lTree_minloH190.root",
    "HZZ4lTree_minloH200.root",
    "HZZ4lTree_minloH250.root",
    "HZZ4lTree_minloH300.root",
    "HZZ4lTree_minloH350.root",
    "HZZ4lTree_minloH400.root",
    "HZZ4lTree_minloH450.root",
    "HZZ4lTree_minloH500.root",
    "HZZ4lTree_minloH550.root",
    "HZZ4lTree_minloH600.root",
    "HZZ4lTree_minloH650.root",
    "HZZ4lTree_minloH700.root",
    "HZZ4lTree_minloH750.root",
    "HZZ4lTree_minloH800.root",
    "HZZ4lTree_minloH850.root",
    "HZZ4lTree_minloH900.root",
    "HZZ4lTree_minloH950.root",
    "HZZ4lTree_minloH1000.root"
  };
  for (int ic=0; ic<3; ic++){
    for (int is=0; is<nSamples_JQCD; is++) tree[0]->Add(Form("%s/%s/%s", cinput_main.Data(), (strchannel[ic]).Data(), (strSamples_JQCD[is]).Data()));
  }
  for (int it=0; it<nMEs; it++){
    tree[it]->SetBranchAddress("NJets30", &NJets30);
    tree[it]->SetBranchAddress("JetPt", &JetPt, &bJetPt);
    tree[it]->SetBranchAddress("JetEta", &JetEta, &bJetEta);
    tree[it]->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
    tree[it]->SetBranchAddress("JetMass", &JetMass, &bJetMass);
    tree[it]->SetBranchAddress("ZZMass", &mzz);
    tree[it]->SetBranchAddress("ZZPt", &ZZPt);
    tree[it]->SetBranchAddress("ZZEta", &ZZEta);
    tree[it]->SetBranchAddress("ZZPhi", &ZZPhi);
    tree[it]->SetBranchAddress("Z1Mass", &m1);
    tree[it]->SetBranchAddress("Z2Mass", &m2);
    tree[it]->SetBranchAddress("helcosthetaZ1", &h1);
    tree[it]->SetBranchAddress("helcosthetaZ2", &h2);
    tree[it]->SetBranchAddress("helphi", &phi);
    tree[it]->SetBranchAddress("costhetastar", &hs);
    tree[it]->SetBranchAddress("phistarZ1", &phi1);
  }

  TFile* foutput = new TFile(Form("pJHUGen_JQCD_HSMHiggs_Comparison_%iTeV.root", sqrts), "recreate");

  TProfile* prJQCD = new TProfile("prJQCD", "", 75, 70, 1570); prJQCD->Sumw2();

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  int nTotalEntries = tree[0]->GetEntries();
  for (int ev = 0; ev < nTotalEntries; ev++){
    tree[0]->GetEntry(ev);
    if (ev%10000==0) cout << "Doing event " << ev << endl;
    if (NJets30==1){
      for (int ij=0; ij<1; ij++){
        jetptetaphimass[ij][0]=JetPt->at(ij);
        jetptetaphimass[ij][1]=JetEta->at(ij);
        jetptetaphimass[ij][2]=JetPhi->at(ij);
        jetptetaphimass[ij][3]=JetMass->at(ij);
      }

      TLorentzVector jet[1], higgs;
      for (int ij=0; ij<1; ij++) jet[ij].SetPtEtaPhiM(jetptetaphimass[ij][0], jetptetaphimass[ij][1], jetptetaphimass[ij][2], jetptetaphimass[ij][3]);
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      float mesq_jqcd=0;
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JQCD);
      mela.computeProdP(mesq_jqcd, true);

      prJQCD->Fill(mzz, mesq_jqcd);

      mela.resetInputEvent();
    }
  }

  foutput->WriteTObject(prJQCD);
  foutput->Close();
  for (int it=0; it<nMEs; it++) delete tree[it];
}

void check_JQCD_13TeV(int sqrts=13){
  int erg_tev=sqrts;
  float mPOLE=125.;
  TString TREE_NAME = "ZZTree/candTree";

  TVar::VerbosityLevel verbosity = TVar::ERROR;
  Mela mela(erg_tev, mPOLE, verbosity);

  short NJets30;
  std::vector<float>* JetPt=0;
  std::vector<float>* JetEta=0;
  std::vector<float>* JetPhi=0;
  std::vector<float>* JetMass=0;
  std::vector<float> myJetPt;
  std::vector<float> myJetEta;
  std::vector<float> myJetPhi;
  std::vector<float> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;
  float jetptetaphimass[2][4];

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float ZZPt, ZZPhi, ZZEta;
  int LepID[4]={ 13, -13, 11, -11 };

  const int nMEs=1;
  TChain* tree[nMEs] ={
    new TChain(TREE_NAME, "")
  };

  TString cinput_main;
  if (sqrts==13) cinput_main = "/scratch0/hep/usarical/CJLST/LHC_13TeV/4l/160225";
  else return;
  const int nSamples_JQCD = 35;
  TString strSamples_JQCD[nSamples_JQCD]={
    "ggH91_GaZ/ZZ4lAnalysis.root",
    "ggH115/ZZ4lAnalysis.root",
    "ggH120/ZZ4lAnalysis.root",
    "ggH124/ZZ4lAnalysis.root",
    "ggH125/ZZ4lAnalysis.root",
    "ggH126/ZZ4lAnalysis.root",
    "ggH130/ZZ4lAnalysis.root",
    "ggH135/ZZ4lAnalysis.root",
    "ggH140/ZZ4lAnalysis.root",
    "ggH145/ZZ4lAnalysis.root",
    "ggH150/ZZ4lAnalysis.root",
    "ggH155/ZZ4lAnalysis.root",
    "ggH160/ZZ4lAnalysis.root",
    "ggH165/ZZ4lAnalysis.root",
    "ggH170/ZZ4lAnalysis.root",
    "ggH175/ZZ4lAnalysis.root",
    "ggH180/ZZ4lAnalysis.root",
    "ggH190/ZZ4lAnalysis.root",
    "ggH200/ZZ4lAnalysis.root",
    "ggH210/ZZ4lAnalysis.root",
    "ggH230/ZZ4lAnalysis.root",
    "ggH250/ZZ4lAnalysis.root",
    "ggH270/ZZ4lAnalysis.root",
    "ggH300/ZZ4lAnalysis.root",
    "ggH350/ZZ4lAnalysis.root",
    "ggH400/ZZ4lAnalysis.root",
    "ggH450/ZZ4lAnalysis.root",
    "ggH500/ZZ4lAnalysis.root",
    "ggH550/ZZ4lAnalysis.root",
    "ggH600/ZZ4lAnalysis.root",
    "ggH700/ZZ4lAnalysis.root",
    "ggH750/ZZ4lAnalysis.root",
    "ggH800/ZZ4lAnalysis.root",
    "ggH900/ZZ4lAnalysis.root",
    "ggH1000/ZZ4lAnalysis.root"
  };
  for (int is=0; is<nSamples_JQCD; is++) tree[0]->Add(Form("%s/%s", cinput_main.Data(), (strSamples_JQCD[is]).Data()));
  for (int it=0; it<nMEs; it++){
    tree[it]->SetBranchAddress("nCleanedJetsPt30", &NJets30);
    tree[it]->SetBranchAddress("JetPt", &JetPt, &bJetPt);
    tree[it]->SetBranchAddress("JetEta", &JetEta, &bJetEta);
    tree[it]->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
    tree[it]->SetBranchAddress("JetMass", &JetMass, &bJetMass);
    tree[it]->SetBranchAddress("ZZMass", &mzz);
    tree[it]->SetBranchAddress("ZZPt", &ZZPt);
    tree[it]->SetBranchAddress("ZZEta", &ZZEta);
    tree[it]->SetBranchAddress("ZZPhi", &ZZPhi);
    tree[it]->SetBranchAddress("Z1Mass", &m1);
    tree[it]->SetBranchAddress("Z2Mass", &m2);
    tree[it]->SetBranchAddress("helcosthetaZ1", &h1);
    tree[it]->SetBranchAddress("helcosthetaZ2", &h2);
    tree[it]->SetBranchAddress("helphi", &phi);
    tree[it]->SetBranchAddress("costhetastar", &hs);
    tree[it]->SetBranchAddress("phistarZ1", &phi1);
  }

  TFile* foutput = new TFile(Form("pJHUGen_JQCD_HSMHiggs_Comparison_%iTeV.root", sqrts), "recreate");

  TProfile* prJQCD = new TProfile("prJQCD", "", 35, 70, 1575); prJQCD->Sumw2();

  mela.setCandidateDecayMode(TVar::CandidateDecay_ZZ);

  int nTotalEntries = tree[0]->GetEntries();
  for (int ev = 0; ev < nTotalEntries; ev++){
    tree[0]->GetEntry(ev);
    if (ev%10000==0) cout << "Doing event " << ev << endl;
    if (NJets30==1){
      for (int ij=0; ij<1; ij++){
        jetptetaphimass[ij][0]=JetPt->at(ij);
        jetptetaphimass[ij][1]=JetEta->at(ij);
        jetptetaphimass[ij][2]=JetPhi->at(ij);
        jetptetaphimass[ij][3]=JetMass->at(ij);
      }

      TLorentzVector jet[1], higgs;
      for (int ij=0; ij<1; ij++) jet[ij].SetPtEtaPhiM(jetptetaphimass[ij][0], jetptetaphimass[ij][1], jetptetaphimass[ij][2], jetptetaphimass[ij][3]);
      higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
      TVector3 boostH = higgs.BoostVector();

      SimpleParticleCollection_t associated;
      associated.push_back(SimpleParticle_t(0, jet[0]));

      TLorentzVector pDaughters[4];
      std::vector<TLorentzVector> daus = mela.calculate4Momentum(mzz, m1, m2, acos(hs), acos(h1), acos(h2), phi1, phi);
      for (int ip=0; ip<min(4, (int)daus.size()); ip++){ pDaughters[ip]=daus.at(ip); pDaughters[ip].Boost(boostH); }
      SimpleParticleCollection_t daughters;
      for (unsigned int idau=0; idau<4; idau++) daughters.push_back(SimpleParticle_t(LepID[idau], pDaughters[idau]));
      mela.setInputEvent(&daughters, &associated, (SimpleParticleCollection_t*)0, false);

      float mesq_jqcd=0;
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JQCD);
      mela.computeProdP(mesq_jqcd, true);

      prJQCD->Fill(mzz, mesq_jqcd);

      mela.resetInputEvent();
    }
  }

  foutput->WriteTObject(prJQCD);
  foutput->Close();
  for (int it=0; it<nMEs; it++) delete tree[it];
}
