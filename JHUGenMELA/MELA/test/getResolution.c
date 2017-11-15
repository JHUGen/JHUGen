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
#include "RooGenericPdf.h"
#include "RooProdPdf.h"
#include "RooNumIntConfig.h"
#include "RooRealIntegral.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooRandom.h"
#include "Mela.h"
#include "MELALinearInterpFunc.h"
#include "MELANCSplineFactory_1D.h"
#include "MELANCSplineFactory_2D.h"


using namespace std;
using namespace RooFit;


TString inputdir_7TeV = "/work-zfs/lhc/ianderso/hep/CJLST/140519/PRODFSR";
TString inputdir_8TeV = "/work-zfs/lhc/ianderso/hep/CJLST/140519b/PRODFSR_8TeV";
TString inputdir_13TeV = "root://lxcms03//data3/Higgs/170222";

struct SimpleEntry{
  int id;
  float trueval;
  float recoval;
  float recotrackval;
  float weight;

  SimpleEntry() : id(0), trueval(0), recoval(0), recotrackval(0), weight(0) {}
  SimpleEntry(int id_, float trueval_, float recoval_, float recotrackval_=0, float weight_=1) : id(id_), trueval(trueval_), recoval(recoval_), recotrackval(recotrackval_), weight(weight_) {}

  bool operator != (const SimpleEntry& other)const{ return trueval!=other.trueval; }
  bool operator == (const SimpleEntry& other)const{ return trueval==other.trueval; }
  bool operator > (const SimpleEntry& other)const{ return trueval>other.trueval; }
  bool operator >= (const SimpleEntry& other)const{ return trueval>=other.trueval; }
  bool operator < (const SimpleEntry& other)const{ return trueval<other.trueval; }
  bool operator <= (const SimpleEntry& other)const{ return trueval<=other.trueval; }
};

template<typename T> void appendVector(std::vector<T>& a, std::vector<T>& b){ a.insert(a.end(), b.begin(), b.end()); }

template<typename T> void addByLowest(std::vector<T>& valArray, T val, bool unique){
  bool inserted = false;
  for (typename std::vector<T>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if (*it>val || (!unique && *it==val)){
      inserted=true;
      valArray.insert(it, val);
      break;
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

template<typename T, typename U> void addByLowest(std::vector<std::pair<T, U>>& valArray, std::vector<std::pair<T, U>>& inArray, bool consecutive=false, bool inputordered=false){
  if (consecutive){
    bool inserted = false;
    typename std::vector<std::pair<T, U>>::iterator inbegin = inArray.begin();
    typename std::vector<std::pair<T, U>>::iterator inend = inArray.end();
    for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
      if ((*it).first>=(*inbegin).first){
        inserted=true;
        if ((*it).second!=(*inbegin).second) valArray.insert(it, inbegin, inend);
        break;
      }
    }
    if (!inserted) appendVector<pair<T, U>>(valArray, inArray);
  }
  else if (!inputordered){
    for (typename std::vector<std::pair<T, U>>::iterator init = inArray.begin(); init<inArray.end(); init++){
      bool inserted = false;
      for (typename std::vector<std::pair<T, U>>::iterator it = valArray.begin(); it<valArray.end(); it++){
        if ((*it).first>=(*init).first){
          inserted=true;
          if ((*it).second!=(*init).second) valArray.insert(it, *init);
          break;
        }
      }
      if (!inserted) valArray.push_back(*init);
    }
  }
  else if (inArray.size()>0){
    typename std::vector<std::pair<T, U>>::iterator infirst = inArray.begin();
    typename std::vector<std::pair<T, U>>::iterator inlast = inArray.end()-1;
    typename std::vector<std::pair<T, U>>::iterator valfirst = valArray.begin();
    typename std::vector<std::pair<T, U>>::iterator vallast = valArray.end()-1;
    while ((*valfirst).first<(*infirst).first) valfirst++;
    while ((*vallast).first>=(*inlast).first) vallast--;
    vallast++;
    inlast++;

    for (typename std::vector<std::pair<T, U>>::iterator init = infirst; init<inlast; init++){
      bool inserted = false;
      for (typename std::vector<std::pair<T, U>>::iterator it = valfirst; it<vallast; it++){
        if ((*it).first>=(*init).first){
          inserted=true;
          if ((*it).second!=(*init).second) valArray.insert(it, *init);
          break;
        }
      }
      if (!inserted) valArray.insert(vallast, *init);
    }
  }
}

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

  vector<SimpleEntry> collection;

  void addEvent(SimpleEntry evt){
    collection.push_back(evt);
  }
  void addEvent(int ev, float trueval, float recoval, float recotrackval=0, float weight=1){
    SimpleEntry tmp(ev, trueval, recoval, recotrackval, weight);
    collection.push_back(tmp);
  }
  void sortByRecoMinusTrueVals(vector<pair<float, int>>& res){
    res.clear();
    for (unsigned int ev=0; ev<collection.size(); ev++){
      float val = collection.at(ev).recoval/collection.at(ev).trueval-1.;
      addByLowest<float, int>(res, val, ev);
    }
  }
  void sortByRecoVals(vector<pair<float, int>>& res){
    res.clear();
    for (unsigned int ev=0; ev<collection.size(); ev++) addByLowest<float, int>(res, collection.at(ev).recoval, ev);
  }
  void sift(double threshold=0.998, bool use_recominustrue=false){
    if (collection.size()==0) return;

    vector<int> take_out;
    vector<pair<float, int>> theEntries;
    if (!use_recominustrue) sortByRecoVals(theEntries);
    else sortByRecoMinusTrueVals(theEntries);
    for (unsigned int im=0; im<2; im++){
      int at99p8ev = std::floor((float(theEntries.size()))*threshold);
      int bin=theEntries.size()-1;
      while (bin>at99p8ev){
        if (
          theEntries.at(at99p8ev).first*2.<theEntries.at(bin).first
          ) addByLowest<int>(take_out, theEntries.at(bin).second, true);
        bin--;
      }
    }
    for (int bin=take_out.size()-1; bin>=0; bin--){
      int t_ev = take_out.at(bin);
      cout << "Attempting to discard event at position " << t_ev << endl;
      if (t_ev>=(int)collection.size()) cerr << "Collection size " << collection.size() << " >= " << t_ev << endl;
      else{
        cout << "Discarding event " << collection.at(t_ev).id << " at position " << t_ev << endl;
        vector<SimpleEntry>::iterator it;
        it=collection.begin()+t_ev; collection.erase(it);
      }
    }
  }
  void adjustWeights(float refth=-1){
    if (collection.size()==0) return;

    vector<pair<float, int>> weight_entry;
    for (unsigned int ev=0; ev<collection.size(); ev++) addByLowest<float, int>(weight_entry, collection.at(ev).weight, ev);
    int bin=weight_entry.size()-1;

    int at99p0ev = (float(weight_entry.size()))*0.99;
    if (at99p0ev==(int)weight_entry.size()) at99p0ev--;
    float threshold = weight_entry.at(at99p0ev).first;
    while (bin>at99p0ev){
      if (
        threshold*2.<weight_entry.at(bin).first
        ){
        if (weight_entry.at(bin).first==0.) cerr << "Weight=0!" << endl;
        collection.at(weight_entry.at(bin).second).weight = pow(threshold, 2)/weight_entry.at(bin).first;
        cout << "Adjusting weight for event " << bin << endl;
      }
      bin--;
    }
    float sum[2]={ 0 };
    for (unsigned int ev=0; ev<collection.size(); ev++){
      sum[0] += collection.at(ev).weight;
      sum[1] += 1;
    }
    if (sum[0]==0. || sum[1]==0.) cerr << "sum[0]=" << sum[0] << " or sum[1]=" << sum[1] << " is invalid." << endl;
    for (unsigned int ev=0; ev<collection.size(); ev++) collection.at(ev).weight *= sum[1]/sum[0];
    if (refth>0.){
      for (unsigned int ev=0; ev<collection.size(); ev++){
        if (collection.at(ev).weight<=refth) continue;
        collection.at(ev).weight = pow(refth, 2)/collection.at(ev).weight;
      }
    }
  }
  void mergeBin(const ExtBin& other){ for (unsigned int ev=0; ev<other.collection.size(); ev++) collection.push_back(other.collection.at(ev)); }
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
      samples.push_back("VBFH150");
      samples.push_back("VBFH155");
      samples.push_back("VBFH160");
      samples.push_back("VBFH165");
      samples.push_back("VBFH170");
      samples.push_back("VBFH175");
      samples.push_back("VBFH180");
      samples.push_back("VBFH190");
      samples.push_back("VBFH210");
      samples.push_back("VBFH230");
      samples.push_back("VBFH250");
      samples.push_back("VBFH270");
      samples.push_back("VBFH300");
      samples.push_back("VBFH350");
      samples.push_back("VBFH450");
      samples.push_back("VBFH500");
      samples.push_back("VBFH550");
      samples.push_back("VBFH600");
      samples.push_back("VBFH700");
      samples.push_back("VBFH750");
      samples.push_back("VBFH800");
      samples.push_back("VBFH900");
      samples.push_back("VBFH1000");
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
      //samples.push_back("ZH210");
      //samples.push_back("ZH230");
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

TGraphErrors* makeGraphFromTH1(TH1* hx, TH1* hy, TString name){
  if (hx->GetNbinsX()!=hy->GetNbinsX()){
    cerr << "Number of bins for x coordinate != those for y" << endl;
    assert(0);
  }
  unsigned int nbins = hx->GetNbinsX();
  double* xexyey[4];
  for (unsigned int ix=0; ix<4; ix++) xexyey[ix] = new double[nbins];
  for (unsigned int bin=0; bin<nbins; bin++){
    xexyey[0][bin] = hx->GetBinContent(bin+1);
    xexyey[1][bin] = hx->GetBinError(bin+1);

    cout << "Bin " << bin << " x-center: " << xexyey[0][bin] << " +- " << xexyey[1][bin] << endl;
    xexyey[2][bin] = hy->GetBinContent(bin+1);
    xexyey[3][bin] = hy->GetBinError(bin+1);
  }
  TGraphErrors* tg = new TGraphErrors(nbins, xexyey[0], xexyey[2], xexyey[1], xexyey[3]);
  tg->SetName(name);
  for (unsigned int ix=0; ix<4; ix++) delete[] xexyey[ix];
  return tg;
}

template<typename T> TGraph* makeGraphFromPair(vector<pair<T,T>> points, TString name){
  if (points.size()==0) return 0;
  unsigned int nbins = points.size();
  double* xy[2];
  for (unsigned int ix=0; ix<2; ix++) xy[ix] = new double[nbins];
  for (unsigned int bin=0; bin<nbins; bin++){
    xy[0][bin] = points[bin].first;
    xy[1][bin] = points[bin].second;
  }
  TGraph* tg = new TGraph(nbins, xy[0], xy[1]);
  tg->SetName(name);
  for (unsigned int ix=0; ix<2; ix++) delete[] xy[ix];
  return tg;
}

// Return A, B, C for A + B*x + C*x**2
vector<double> solveQuadratic(vector<pair<double, double>> xy){
  vector<double> res;
  if (xy.size()==3){
    double& x1 = xy[0].first;
    double& x2 = xy[1].first;
    double& x3 = xy[2].first;
    double det = pow(x1, 2)*(x3-x2) + pow(x2, 2)*(x1-x3) + pow(x3, 2)*(x2-x1);
    double Ainv[3][3]={
      { x2*x3*(x3-x2), x1*x3*(x1-x3), x1*x2*(x2-x1) },
      { (x2-x3)*(x2+x3), (x3-x1)*(x3+x1), (x1-x2)*(x1+x2) },
      { -(x2-x3), -(x3-x1), -(x1-x2) }
    };
    for (unsigned int i=0; i<3; i++){
      double sum=0;
      for (unsigned int j=0; j<3; j++) sum += Ainv[i][j]/det*xy[j].second;
      res.push_back(sum);
    }
  }
  return res;
}
// Return A, B, C for A*exp(-(x-B)**2/(2*C))
vector<double> solveGaussian(vector<pair<double, double>> xy){
  vector<double> res;
  if (xy.size()==3){
    vector<pair<double, double>> xlny;
    for(auto ipair : xy) xlny.push_back(pair<double, double>(ipair.first, log(ipair.second)));
    vector<double> qc = solveQuadratic(xlny);

    double A = exp(qc[0]-pow(qc[1], 2)/qc[2]);
    double B = -qc[1]*0.5/qc[2];
    double C = -0.5/qc[2];
    res.push_back(A);
    res.push_back(B);
    res.push_back(C);
  }
  return res;
}

void addPoint(TGraph*& tg, double x, double y){
  TString strname = tg->GetName();
  TString strtitle = tg->GetTitle();
  TString strxtitle = tg->GetXaxis()->GetTitle();
  TString strytitle = tg->GetYaxis()->GetTitle();

  bool hasErrors = false;

  vector<double> xarray;
  vector<double> yarray;
  xarray.push_back(x);
  yarray.push_back(y);
  for (int ip=0; ip<tg->GetN(); ip++){
    if (tg->GetX()[ip]!=x){
      xarray.push_back(tg->GetX()[ip]);
      yarray.push_back(tg->GetY()[ip]);
    }
  }
  vector<pair<double, int>> xorder;
  for (unsigned int ip=0; ip<xarray.size(); ip++) addByLowest<double, int>(xorder, xarray.at(ip), ip);

  double* xynew[2];
  for (unsigned int i=0; i<2; i++) xynew[i] = new double[xorder.size()];
  for (unsigned int ip=0; ip<xarray.size(); ip++){
    unsigned int pos = xorder[ip].second;
    xynew[0][ip] = xarray[pos];
    xynew[1][ip] = yarray[pos];
  }

  delete tg;

  tg = new TGraph(xorder.size(), xynew[0], xynew[1]);
  tg->SetName(strname);
  tg->SetTitle(strtitle);
  tg->GetXaxis()->SetTitle(strxtitle);
  tg->GetYaxis()->SetTitle(strytitle);
  for (unsigned int i=0; i<2; i++) delete[] xynew[i];
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


/* SPECIFIC COMMENT: NONE */
void getResolution_mH_signal(float sqrts, TString strprod){
  TString TREE_NAME;
  const bool writeFinalTree=false;
  const TString strTrueBranch = "GenHMass";
  const TString strRecoBranch = "ZZMass";

  short Z1Flav, Z2Flav;
  float varTrue, varReco;

  vector<SimpleEntry> index[3];
  TString strchannel[3]={ "4mu", "4e", "2mu2e" };

  vector<TString> dumappend;
  vector<TString> strSamples;
  if (strprod=="gg"){
    dumappend = constructSamplesList("JJQCD", sqrts);
    appendVector<TString>(strSamples, dumappend);
    appendVector<TString>(strSamples, dumappend);
    dumappend = constructSamplesList("gg_Sig_JHUGen", sqrts);
  }
  else if (strprod=="VBF"){
    dumappend = constructSamplesList("JJVBF", sqrts);
    appendVector<TString>(strSamples, dumappend);
  }
  else if (strprod=="ZH"){
    dumappend = constructSamplesList("ZH", sqrts);
    appendVector<TString>(strSamples, dumappend);
  }
  else if (strprod=="WH"){
    dumappend = constructSamplesList("WH", sqrts);
    appendVector<TString>(strSamples, dumappend);
  }
  else{
    cerr << "Production " << strprod << " is unknown." << endl;
    assert(0);
  }

  TFile* foutput = TFile::Open(Form("resolution_m4l_truem4l_%s.root", strprod.Data()), "recreate");

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  TString cinput_main;

  if (sqrts<13){
    TREE_NAME = "SelectedTree";
    if (sqrts==8.)cinput_main = inputdir_8TeV;
    else cinput_main = inputdir_7TeV;

    for (int is=0; is<(int)strSamples.size(); is++){
      for (unsigned int ic=0; ic<3; ic++){
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
              tree->SetBranchStatus(strTrueBranch, 1); tree->SetBranchAddress(strTrueBranch, &varTrue);
              tree->SetBranchStatus(strRecoBranch, 1); tree->SetBranchAddress(strRecoBranch, &varReco);
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

  }
  else /*if (sqrts==13.)*/{
    TREE_NAME = "ZZTree/candTree";
    cinput_main = inputdir_13TeV;
    //for (int is=0; is<2; is++){
    for (int is=0; is<(int)strSamples.size(); is++){
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
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus(strTrueBranch, 1); tree->SetBranchAddress(strTrueBranch, &varTrue);
            tree->SetBranchStatus(strRecoBranch, 1); tree->SetBranchAddress(strRecoBranch, &varReco);
            tree->SetBranchStatus("Z1Flav", 1); tree->SetBranchAddress("Z1Flav", &Z1Flav);
            tree->SetBranchStatus("Z2Flav", 1); tree->SetBranchAddress("Z2Flav", &Z2Flav);
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

  unsigned ev_acc=0;
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
    SimpleEntry theEntry(ev, varTrue, varReco, 0, 1);
    addByLowest(index[ic], theEntry, false);
    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;

  for (unsigned int ic=0; ic<3; ic++){
    float firstVal=index[ic].at(0).trueval;
    float lastVal=index[ic].at(index[ic].size()-1).trueval;
    float infimum = (float)((int)firstVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nEntries << " | truth = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

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
      binning[ix]=(index[ic][ix*ev_stepsize-1].trueval+index[ic][ix*ev_stepsize].trueval)*0.5;
      ExtBin tmpbin;
      tmpbin.binlow = binning[ix-1]; tmpbin.binhigh = binning[ix];
      for (int bin=0; bin<ev_stepsize; bin++){ tmpbin.addEvent(index[ic][(ix-1)*ev_stepsize+bin]); }
      binList.push_back(tmpbin);
      cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ic][ix*ev_stepsize].id << ", step " << ix*ev_stepsize << "]" << endl;
    }
    ExtBin tmpbin;
    tmpbin.binlow = binning[nbins-1]; tmpbin.binhigh = binning[nbins];
    for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index[ic].size(); bin++){ tmpbin.addEvent(index[ic][bin]); }
    binList.push_back(tmpbin);
    cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
    cout << "Bin list has the following bins:" << endl;
    for (unsigned int ib=0; ib<binList.size(); ib++){
      cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
    }

    foutput->cd();

    TProfile* hvar = new TProfile(Form("varTrue_%s", strchannel[ic].Data()), "", nbins, binning); hvar->Sumw2();
    TProfile* havg_varReco = new TProfile(Form("avg_varReco_%s", strchannel[ic].Data()), "", nbins, binning); havg_varReco->Sumw2();
    TH2F* hvarTrue_varReco = new TH2F(Form("varTrue_varReco_%s", strchannel[ic].Data()), "", nbins, binning, nbins, binning); hvarTrue_varReco->Sumw2();

    TTree* newtree=0;
    if (writeFinalTree){
      newtree = new TTree(Form("FinalTree_%s", strchannel[ic].Data()), "");
      newtree->Branch("varTrue", &varTrue);
      newtree->Branch("varReco", &varReco);
    }

    unsigned int ctr=0;
    for (unsigned int bin=0; bin<binList.size(); bin++){
      cout << "Bin " << bin << " is now being scrutinized..." << endl;

      //binList.at(bin).sift();
      binList.at(bin).adjustWeights();

      vector<pair<float, int>> recominustruevals_entries;
      binList.at(bin).sortByRecoMinusTrueVals(recominustruevals_entries);
      vector<float> recominustruevalbounds;
      unsigned int nrecototal = recominustruevals_entries.size();
      unsigned int nrecobins = 10;
      unsigned int recoinc = nrecototal / nrecobins;
      float minval = std::floor(recominustruevals_entries.at(0).first);
      float maxval = std::ceil(recominustruevals_entries.at(nrecototal-1).first);
      recominustruevalbounds.push_back(minval);
      for (unsigned int i=1; i<nrecobins; i++){
        float a = recominustruevals_entries.at(i*recoinc).first;
        float b = recominustruevals_entries.at(i*recoinc-1).first;
        recominustruevalbounds.push_back((a+b)*0.5);
      }
      recominustruevalbounds.push_back(maxval);
      cout << "Bin " << bin << " has the following reco-true boundaries:";
      for (unsigned int j=0; j<recominustruevalbounds.size(); j++) cout << " " << recominustruevalbounds.at(j);
      cout << endl;

      foutput->cd();
      TH1F* h_RecoMinusTrueVals = new TH1F(Form("RecoMinusTrueVal_%s_Bin%i", strchannel[ic].Data(), bin), "", nrecobins, recominustruevalbounds.data()); h_RecoMinusTrueVals->Sumw2();
      TProfile* p_RecoMinusTrueVals = new TProfile(Form("avg_RecoMinusTrueVal_%s_Bin%i", strchannel[ic].Data(), bin), "", nrecobins, recominustruevalbounds.data()); p_RecoMinusTrueVals->Sumw2();
      for (unsigned int i=0; i<nrecototal; i++){
        float& val = recominustruevals_entries.at(i).first;
        int& ev = recominustruevals_entries.at(i).second;
        float& weight = binList.at(bin).collection.at(ev).weight;
        if (weight==0.) cerr << "Weight is 0!" << endl;
        else if (isnan(weight) || isinf(weight)) cerr << "Weight " << weight << " is nan or inf" << endl;
        if (isnan(val) || isinf(val)) cerr << "Value " << val << " is nan or inf" << endl;
        h_RecoMinusTrueVals->Fill(val, weight);
        p_RecoMinusTrueVals->Fill(val, val, weight);
      }
      const TAxis* axis = h_RecoMinusTrueVals->GetXaxis();
      for (int bin=1; bin<=h_RecoMinusTrueVals->GetNbinsX(); bin++){
        float width = axis->GetBinWidth(bin);
        h_RecoMinusTrueVals->SetBinContent(bin, h_RecoMinusTrueVals->GetBinContent(bin)/width);
        h_RecoMinusTrueVals->SetBinError(bin, h_RecoMinusTrueVals->GetBinError(bin)/width);
      }
      foutput->WriteTObject(h_RecoMinusTrueVals);
      foutput->WriteTObject(p_RecoMinusTrueVals);
      delete p_RecoMinusTrueVals;
      delete h_RecoMinusTrueVals;

      for (unsigned int ev=0; ev<binList.at(bin).collection.size(); ev++){
        varTrue = binList.at(bin).collection.at(ev).trueval;
        varReco = binList.at(bin).collection.at(ev).recoval;
        havg_varReco->Fill(varTrue, varReco);
        hvar->Fill(varTrue, varTrue);
        hvarTrue_varReco->Fill(varTrue, varReco);
        if (writeFinalTree) newtree->Fill();
      }
    }

    TGraphErrors* tg = makeGraphFromTH1(hvar, havg_varReco, Form("tg_%s", havg_varReco->GetName()));
    foutput->WriteTObject(tg);
    delete tg;

    for (int ix=0; ix<=hvarTrue_varReco->GetNbinsX()+1; ix++){
      double integral = hvarTrue_varReco->Integral(ix, ix, 0, hvarTrue_varReco->GetNbinsY()+1);
      if (integral!=0.){
        for (int iy=0; iy<=hvarTrue_varReco->GetNbinsY()+1; iy++){
          double bincontent = hvarTrue_varReco->GetBinContent(ix, iy);
          double binerror = hvarTrue_varReco->GetBinError(ix, iy);
          bincontent /= integral;
          binerror /= integral; // I know this is wrong.
          hvarTrue_varReco->SetBinContent(ix, iy, bincontent);
          hvarTrue_varReco->SetBinError(ix, iy, binerror);
        }
      }
    }
    foutput->WriteTObject(hvarTrue_varReco);
    foutput->WriteTObject(havg_varReco);
    foutput->WriteTObject(hvar);
    if (writeFinalTree){
      foutput->WriteTObject(newtree);
      delete newtree;
    }
    delete hvarTrue_varReco;
    delete havg_varReco;
    delete hvar;
    delete[] binning;
  }
  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
  foutput->Close();
}

/*
SPECIFIC COMMENT:
- Final product is reco. BW(mJJ), so if division/multiplication is needed, divide by 2*mJJ.
- BW is already convoluted. For future use, added also the transfer functions to output.
*/
void getResolution_mJJ(float sqrts, TString strprod, bool debug=false){
  if (sqrts<13) assert(0);

  Mela mela(sqrts, 125);

  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";
  const bool writeFinalTree = false && debug;

  short Z1Flav, Z2Flav, njets;
  float varTrue, varReco, ZZMass, weight;
  typedef float T;

  std::vector<T>* LepPt=0;
  std::vector<T>* LepEta=0;
  std::vector<T>* LepPhi=0;
  std::vector<short>* LepLepId=0;
  std::vector<T>* JetPt=0;
  std::vector<T>* JetEta=0;
  std::vector<T>* JetPhi=0;
  std::vector<T>* JetMass=0;

  std::vector<T>* LHEMotherPz=0;
  std::vector<T>* LHEMotherE=0;
  std::vector<short>* LHEMotherId=0;
  std::vector<T>* LHEDaughterPt=0;
  std::vector<T>* LHEDaughterEta=0;
  std::vector<T>* LHEDaughterPhi=0;
  std::vector<T>* LHEDaughterMass=0;
  std::vector<short>* LHEDaughterId=0;
  std::vector<T>* LHEAssociatedParticlePt=0;
  std::vector<T>* LHEAssociatedParticleEta=0;
  std::vector<T>* LHEAssociatedParticlePhi=0;
  std::vector<T>* LHEAssociatedParticleMass=0;
  std::vector<short>* LHEAssociatedParticleId=0;

  float xsec=1, overallEventWeight=1;

  vector<SimpleEntry> index;
  TString strchannel[3]={ "4mu", "4e", "2mu2e" };

  vector<TString> dumappend;
  vector<TString> strSamples;
  if (strprod=="gg"){
    dumappend = constructSamplesList("JJQCD", sqrts);
    appendVector<TString>(strSamples, dumappend);
    dumappend = constructSamplesList("gg_Sig_JHUGen", sqrts);
  }
  else if (strprod=="VBF"){
    dumappend = constructSamplesList("JJVBF", sqrts);
    appendVector<TString>(strSamples, dumappend);
  }
  else if (strprod=="ZH"){
    dumappend = constructSamplesList("ZH", sqrts);
    appendVector<TString>(strSamples, dumappend);
  }
  else if (strprod=="WH"){
    dumappend = constructSamplesList("WH", sqrts);
    appendVector<TString>(strSamples, dumappend);
  }
  else if (strprod=="JJEW"){
    dumappend = constructSamplesList("JJVBF", sqrts);
    appendVector<TString>(strSamples, dumappend);
    dumappend = constructSamplesList("ZH", sqrts);
    appendVector<TString>(strSamples, dumappend);
    dumappend = constructSamplesList("WH", sqrts);
    appendVector<TString>(strSamples, dumappend);
  }
  else{
    cerr << "Production " << strprod << " is unknown." << endl;
    assert(0);
  }

  TFile* foutput = TFile::Open(Form("resolution_mJJ_recoVStrue_%s_%.0fTeV.root", strprod.Data(), sqrts), "recreate"); // Y vs X

  vector<TFile*> finputList;
  vector<TTree*> treeList;
  int nEntries=0;
  TString cinput_main;

  unordered_map<TTree*, pair<float, float>> nGenMap;
  unordered_map<TTree*, float> mass_map;
  unordered_map<float, vector<TTree*>> mass_sample_map;

  if (sqrts<13){
    TREE_NAME = "SelectedTree";
    if (sqrts==8.)cinput_main = inputdir_8TeV;
    else cinput_main = inputdir_7TeV;

    for (int is=0; is<(int)strSamples.size(); is++){
      for (unsigned int ic=0; ic<3; ic++){
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
              tree->SetBranchStatus("MC_weight", 1); tree->SetBranchAddress("MC_weight", &overallEventWeight);
              tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &ZZMass);
              tree->SetBranchStatus("NJets30", 1); tree->SetBranchAddress("nCleanedJetsPt30", &njets);
              tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
              tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
              tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
              tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);

              TH1F* htmp = (TH1F*)finput->Get(COUNTERS_NAME);
              pair<float, float> nsum(1, 1); // MC_weight already divides by number of generated events
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
    }

  }
  else /*if (sqrts==13.)*/{
    TREE_NAME = "ZZTree/candTree";
    cinput_main = inputdir_13TeV;
    //for (int is=0; is<2; is++){
    for (int is=0; is<(int)strSamples.size(); is++){
      TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[is]).Data());
      cout << "Opening file " << cinput << "..." << endl;
      TFile* finput = TFile::Open(cinput, "read");
      cout << "File open attempt passed." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("xsec", 1); tree->SetBranchAddress("xsec", &xsec);
            tree->SetBranchStatus("overallEventWeight", 1); tree->SetBranchAddress("overallEventWeight", &overallEventWeight);
            tree->SetBranchStatus("ZZMass", 1); tree->SetBranchAddress("ZZMass", &ZZMass);

            tree->SetBranchStatus("LepPt", 1); tree->SetBranchAddress("LepPt", &LepPt);
            tree->SetBranchStatus("LepEta", 1); tree->SetBranchAddress("LepEta", &LepEta);
            tree->SetBranchStatus("LepPhi", 1); tree->SetBranchAddress("LepPhi", &LepPhi);
            tree->SetBranchStatus("LepLepId", 1); tree->SetBranchAddress("LepLepId", &LepLepId);

            tree->SetBranchStatus("nCleanedJetsPt30", 1); tree->SetBranchAddress("nCleanedJetsPt30", &njets);
            tree->SetBranchStatus("JetPt", 1); tree->SetBranchAddress("JetPt", &JetPt);
            tree->SetBranchStatus("JetEta", 1); tree->SetBranchAddress("JetEta", &JetEta);
            tree->SetBranchStatus("JetPhi", 1); tree->SetBranchAddress("JetPhi", &JetPhi);
            tree->SetBranchStatus("JetMass", 1); tree->SetBranchAddress("JetMass", &JetMass);

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

            tree->GetEntry(0);
            cout << "Cross section = " << xsec << endl;

            TH1F* htmp = (TH1F*)finput->Get(COUNTERS_NAME);
            pair<float, float> nsum(htmp->GetBinContent(1), htmp->GetBinContent(40)); // Add PU reweighting
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
  }

  for (auto it = mass_sample_map.begin(); it != mass_sample_map.end(); ++it){
    float sum_ngen=0;
    float sum_xsec=0;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      it->second.at(ix)->GetEntry(0);
      sum_ngen += nGenMap[it->second.at(ix)].first;
      sum_xsec += xsec;
    }
    float overallWeight = sum_ngen/sum_xsec;
    for (unsigned int ix=0; ix<it->second.size(); ix++){
      cout << "Sum Hep MC weights in tree " << ix << " / " << it->second.size() << " was " << nGenMap[it->second.at(ix)].second << " over " << nGenMap[it->second.at(ix)].first << " total gen. events." << endl;
      nGenMap[it->second.at(ix)].first = overallWeight/nGenMap[it->second.at(ix)].second;
      cout << "Event scale for tree " << ix << " / " << it->second.size() << " at pole mass " << it->first << " = " << nGenMap[it->second.at(ix)].first << endl;
    }
  }

  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;

  foutput->cd();
  TTree* newtree=0;
  if (writeFinalTree){
    newtree = new TTree("FinalTree", "");
    newtree->Branch("varTrue", &varTrue);
    newtree->Branch("varReco", &varReco);
    newtree->Branch("ZZMass", &ZZMass);
    newtree->Branch("weight", &weight);
  }

  unsigned ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    getEntry(treeList, ev);
    bool doProcess=
      njets>=2 && LHEAssociatedParticlePt->size()>=2 && LepPt->size()==4
      ;
    if (!doProcess) continue;
    TTree* tree = findTree(treeList, ev);
    weight = nGenMap[tree].first*xsec*overallEventWeight;

    vector<MELAParticle*> partList;

    // Reconstruct true associated jet information
    vector<MELAParticle*> ttot, tg, tq, tothers;
    for (unsigned int i=0; i<LHEAssociatedParticlePt->size(); i++){
      int id = LHEAssociatedParticleId->at(i);
      TLorentzVector tmp(0, 0, 0, 0);
      tmp.SetPtEtaPhiM(
        LHEAssociatedParticlePt->at(i),
        LHEAssociatedParticleEta->at(i),
        LHEAssociatedParticlePhi->at(i),
        LHEAssociatedParticleMass->at(i)
        );

      MELAParticle* part = new MELAParticle(id, tmp);
      if (PDGHelpers::isAQuark(id)) tq.push_back(part);
      else if (PDGHelpers::isAGluon(id)) tg.push_back(part);
      else tothers.push_back(part);
      partList.push_back(part);
    }

    // Reconstruct reco daughter information
    vector<MELAParticle*> daus;
    TLorentzVector pSumDaus(0, 0, 0, 0);
    int Z1Z2id=1;
    for (unsigned int i=0; i<LepPt->size(); i++){
      int id = LepLepId->at(i);
      float mass=0;
      if (std::abs(id)==11) mass = PDGHelpers::m_el;
      else if (std::abs(id)==13) mass = PDGHelpers::m_mu;
      else if (std::abs(id)==15) mass = PDGHelpers::m_tau;
      TLorentzVector tmp(0, 0, 0, 0);
      tmp.SetPtEtaPhiM(
        LepPt->at(i),
        LepEta->at(i),
        LepPhi->at(i),
        mass
        );
      MELAParticle* part = new MELAParticle(id, tmp);
      daus.push_back(part);
      pSumDaus = pSumDaus + tmp;
      Z1Z2id *= id;
      partList.push_back(part);
    }

    // Reconstruct true daughter information, swap true daughters with associated true quarks and leptons if needed
    vector<MELACandidate*> candList;
    MELACandidate* theCand=0; float minDM=9e9;

    vector<MELAParticle*> tdaus;
    TLorentzVector pSumTrueDaus(0, 0, 0, 0);
    int tZ1Z2id=1;
    vector<MELAParticle*> tdaus_split[2][2];
    for (unsigned int i=0; i<LHEDaughterPt->size(); i++){
      int id = LHEDaughterId->at(i);
      TLorentzVector tmp(0, 0, 0, 0);
      tmp.SetPtEtaPhiM(
        LHEDaughterPt->at(i),
        LHEDaughterEta->at(i),
        LHEDaughterPhi->at(i),
        LHEDaughterMass->at(i)
        );
      MELAParticle* part = new MELAParticle(id, tmp);
      tdaus.push_back(part);
      int ipart=-1;
      int jdau=-1;
      if (part->id>0) ipart=0;
      else ipart=1;
      if (std::abs(id)==11) jdau=0;
      else if (std::abs(id)==13) jdau=1;
      if (ipart>=0 && jdau>=0) tdaus_split[jdau][ipart].push_back(part);
      if (PDGHelpers::isAQuark(std::abs(id))) tq.push_back(part);
      partList.push_back(part);
    }

    for (auto it=tothers.begin(); it<tothers.end(); it++){
      MELAParticle* part = *it;
      int ipart=-1;
      int jdau=-1;
      if (part->id>0) ipart=0;
      else ipart=1;
      if (std::abs(part->id)==11) jdau=0;
      else if (std::abs(part->id)==13) jdau=1;
      if (ipart>=0 && jdau>=0) tdaus_split[jdau][ipart].push_back(part);
    }

    //cout << "Begin tmpV construction" << endl;
    std::vector<MELAParticle*> tmpVhandle;
    for (int c=0; c<2; c++){
      for (unsigned int i=0; i<tdaus_split[c][0].size(); i++){
        for (unsigned int j=0; j<tdaus_split[c][1].size(); j++){
          TLorentzVector pV = tdaus_split[c][0].at(i)->p4+tdaus_split[c][1].at(j)->p4;
          MELAParticle* V = new MELAParticle(23, pV);
          V->addDaughter(tdaus_split[c][0].at(i));
          V->addDaughter(tdaus_split[c][1].at(j));
          tmpVhandle.push_back(V);
        }
      }
    }
    //cout << "Begin candList construction" << endl;
    for (unsigned int i=0; i<tmpVhandle.size(); i++){
      for (unsigned int j=i+1; j<tmpVhandle.size(); j++){
        MELAParticle* Vi1 = tmpVhandle.at(i)->getDaughter(0);
        MELAParticle* Vi2 = tmpVhandle.at(i)->getDaughter(1);
        MELAParticle* Vj1 = tmpVhandle.at(j)->getDaughter(0);
        MELAParticle* Vj2 = tmpVhandle.at(j)->getDaughter(1);
        if (Vi1==Vj1 || (Vi2==Vj2 && Vi2 != 0)) continue;

        TLorentzVector pH(0, 0, 0, 0);
        if (Vi1!=0) pH = pH + Vi1->p4;
        if (Vi2!=0) pH = pH + Vi2->p4;
        if (Vj1!=0) pH = pH + Vj1->p4;
        if (Vj2!=0) pH = pH + Vj2->p4;
        MELACandidate* cand = new MELACandidate(25, pH, true);

        if (Vi1!=0) cand->addDaughter(Vi1);
        if (Vi2!=0) cand->addDaughter(Vi2);
        if (Vj1!=0) cand->addDaughter(Vj1);
        if (Vj2!=0) cand->addDaughter(Vj2);

        TVar::CandidateDecayMode defaultHDecayMode = PDGHelpers::HDecayMode;
        PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ZZ);
        cand->sortDaughters();
        PDGHelpers::setCandidateDecayMode(defaultHDecayMode);
        candList.push_back(cand);
      }
    }
    //cout << "Delete tmpVs" << endl;
    for (unsigned int i=0; i<tmpVhandle.size(); i++) delete tmpVhandle.at(i);
    //cout << "Choose best cand out of " << candList.size() << " cands" << endl;
    for (unsigned int icand=0; icand<candList.size(); icand++){
      MELACandidate* cand = candList.at(icand);
      float DM = fabs(cand->m()-ZZMass);
      if (minDM>DM){
        minDM=DM;
        theCand=cand;
      }
    }
    /*
    if (theCand==0){
    cout << "Best cand dne. True daughters are " << endl;
    for (unsigned int i=0; i<tdaus.size(); i++) cout << " - " << tdaus.at(i)->id << " " << tdaus.at(i)->p4.X() << " " << tdaus.at(i)->p4.Y() << " " << tdaus.at(i)->p4.Z() << " " << tdaus.at(i)->p4.T() << endl;
    cout << "Others:" << endl;
    for (unsigned int i=0; i<tothers.size(); i++) cout << " - " << tothers.at(i)->id << " " << tothers.at(i)->p4.X() << " " << tothers.at(i)->p4.Y() << " " << tothers.at(i)->p4.Z() << " " << tothers.at(i)->p4.T() << endl;
    }
    */
    if (theCand!=0){ for (int idau=0; idau<theCand->getNDaughters(); idau++) tZ1Z2id *= theCand->getDaughter(idau)->id; }
    if (tq.size()+tg.size()<2 || Z1Z2id!=tZ1Z2id || theCand==0){
      for (unsigned int icand=0; icand<candList.size(); icand++) delete candList.at(icand);
      for (unsigned int i=0; i<partList.size(); i++) delete partList.at(i);
      continue;
    }
    for (int idau=0; idau<theCand->getNDaughters(); idau++) pSumTrueDaus += theCand->getDaughter(idau)->p4;

    TVector3 boostTrueDaus = pSumTrueDaus.BoostVector();
    TVector3 boostDaus = pSumDaus.BoostVector();

    // True mothers
    vector<MELAParticle*> tmothers;
    for (unsigned int i=0; i<LHEMotherPz->size(); i++){
      int id = LHEMotherId->at(i);
      TLorentzVector tmp(0, 0, LHEMotherPz->at(i), LHEMotherE->at(i));
      MELAParticle* part = new MELAParticle(id, tmp);
      tmothers.push_back(part);
      partList.push_back(part);
    }

    vector<MELAParticle*> jets;
    for (unsigned int i=0; i<(unsigned int)min(2, (int)JetPt->size()); i++){
      int id = 0;
      TLorentzVector tmp(0, 0, 0, 0);
      tmp.SetPtEtaPhiM(
        JetPt->at(i),
        JetEta->at(i),
        JetPhi->at(i),
        JetMass->at(i)
        );
      MELAParticle* part = new MELAParticle(id, tmp);
      jets.push_back(part);
      partList.push_back(part);
    }

    for (unsigned int i=0; i<daus.size(); i++){ daus.at(i)->p4.Boost(-boostDaus); daus.at(i)->p4.Boost(boostTrueDaus); }
    for (unsigned int i=0; i<jets.size(); i++){ jets.at(i)->p4.Boost(-boostDaus); jets.at(i)->p4.Boost(boostTrueDaus); }
    varReco = (jets.at(0)->p4 + jets.at(1)->p4).M();

    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;

    for (unsigned int i=0; i<tg.size(); i++){
      int jfound=-1; float dRmin=9e10;
      MELAParticle* gluon = tg.at(i);

      for (unsigned int j=0; j<tq.size(); j++){
        MELAParticle* quark = tq.at(j);
        float dR = 1./fabs(quark->dot(gluon));
        if (dRmin>dR){
          dRmin=dR;
          jfound=j;
        }
      }
      for (unsigned int j=0; j<tmothers.size(); j++){
        MELAParticle* quark = tmothers.at(j);
        if (!PDGHelpers::isAQuark(quark->id)) continue;
        float dR = 1./fabs(quark->dot(gluon));
        if (dRmin>dR){
          dRmin=dR;
          jfound=-1;
        }
      }

      if (jfound>=0){
        int& id = tq.at(jfound)->id;
        TLorentzVector mom = tq.at(jfound)->p4 + gluon->p4;
        MELAParticle* composite = new MELAParticle(id, mom);
        composite->addDaughter(tq.at(jfound));
        composite->addDaughter(gluon);
        ttot.push_back(composite);
        partList.push_back(composite);
      }
    }

    TLorentzVector sum(0, 0, 0, 0);
    vector<MELAParticle*> matchedpartons;
    for (unsigned int ijet=0; ijet<2; ijet++){
      MELAParticle* jet = jets.at(ijet);
      MELAParticle* matched=0; float dRmin=9e9;

      // First check Q+G
      for (unsigned int i=0; i<ttot.size(); i++){
        MELAParticle* part = ttot.at(i);
        bool unmatched=true;
        for (unsigned int im=0; im<matchedpartons.size(); im++){
          if (
            matchedpartons.at(im)==part
            ||
            matchedpartons.at(im)==part->getDaughter(0)
            ||
            matchedpartons.at(im)==part->getDaughter(1)
            ){
            unmatched=false; break;
          }
        }
        if (!unmatched) continue;
        float dR = part->deltaR(jet);
        if (dRmin>dR){
          dRmin=dR;
          matched = part;
        }
      }
      // Check individual Q
      for (unsigned int i=0; i<tq.size(); i++){
        MELAParticle* part = tq.at(i);
        bool unmatched=true;
        for (unsigned int im=0; im<matchedpartons.size(); im++){
          if (matchedpartons.at(im)==part) unmatched=false;
          for (int idau=0; idau<matchedpartons.at(im)->getNDaughters(); idau++){
            if (matchedpartons.at(im)->getDaughter(idau)==part){ unmatched=false; break; }
          }
          if (!unmatched) break;
        }
        if (!unmatched) continue;
        float dR = part->deltaR(jet);
        if (dRmin>dR){
          dRmin=dR;
          matched = part;
        }
      }
      // Check individual G
      for (unsigned int i=0; i<tg.size(); i++){
        MELAParticle* part = tg.at(i);
        bool unmatched=true;
        for (unsigned int im=0; im<matchedpartons.size(); im++){
          if (matchedpartons.at(im)==part) unmatched=false;
          for (int idau=0; idau<matchedpartons.at(im)->getNDaughters(); idau++){
            if (matchedpartons.at(im)->getDaughter(idau)==part){ unmatched=false; break; }
          }
          if (!unmatched) break;
        }
        if (!unmatched) continue;
        float dR = part->deltaR(jet);
        if (dRmin>dR){
          dRmin=dR;
          matched = part;
        }
      }

      if (matched!=0){
        sum = sum + matched->p4;
        matchedpartons.push_back(matched);
      }
    }

    varTrue = sum.M();

    if (varTrue<0. || ev_acc%100000==0/* || varReco/varTrue<0.4 || varReco/varTrue>1./0.4*/){
      cout << "True mJJ = " << varTrue << ", reco mJJ = " << varReco << ". Event info:" << endl;
      cout << "- Best cand mass = " << pSumTrueDaus.M() << " vs. reco mass = " << ZZMass << endl;
      cout << " - Reco jets:" << endl;
      for (unsigned int i=0; i<jets.size(); i++) cout << " - " << jets.at(i)->p4.X() << " " << jets.at(i)->p4.Y() << " " << jets.at(i)->p4.Z() << " " << jets.at(i)->p4.T() << endl;
      cout << " - " << " - Matched truth:" << endl;
      for (unsigned int i=0; i<matchedpartons.size(); i++) cout << " - " << matchedpartons.at(i)->id << " " << matchedpartons.at(i)->p4.X() << " " << matchedpartons.at(i)->p4.Y() << " " << matchedpartons.at(i)->p4.Z() << " " << matchedpartons.at(i)->p4.T() << endl;
      cout << " - " << " - True partons:" << endl;
      for (unsigned int i=0; i<tg.size(); i++) cout << " - " << tg.at(i)->id << " " << tg.at(i)->p4.X() << " " << tg.at(i)->p4.Y() << " " << tg.at(i)->p4.Z() << " " << tg.at(i)->p4.T() << endl;
      for (unsigned int i=0; i<tq.size(); i++) cout << " - " << tq.at(i)->id << " " << tq.at(i)->p4.X() << " " << tq.at(i)->p4.Y() << " " << tq.at(i)->p4.Z() << " " << tq.at(i)->p4.T() << endl;
      cout << " - " << " - Merged partons:" << endl;
      for (unsigned int i=0; i<ttot.size(); i++){
        cout << " - " << ttot.at(i)->id << " " << ttot.at(i)->p4.X() << " " << ttot.at(i)->p4.Y() << " " << ttot.at(i)->p4.Z() << " " << ttot.at(i)->p4.T() << endl;
        for (int j=0; j<ttot.at(i)->getNDaughters(); j++)cout << " ---Dau" << j+1 << ": " << ttot.at(i)->getDaughter(j)->id << " " << ttot.at(i)->getDaughter(j)->p4.X() << " " << ttot.at(i)->getDaughter(j)->p4.Y() << " " << ttot.at(i)->getDaughter(j)->p4.Z() << " " << ttot.at(i)->getDaughter(j)->p4.T() << endl;
      }
      cout << " - " << " - Other true asssociated:" << endl;
      for (unsigned int i=0; i<tothers.size(); i++) cout << " - " << tothers.at(i)->id << " " << tothers.at(i)->p4.X() << " " << tothers.at(i)->p4.Y() << " " << tothers.at(i)->p4.Z() << " " << tothers.at(i)->p4.T() << endl;
      cout << "True candidate summary:" << endl;
      TUtil::PrintCandidateSummary(theCand);
    }

    for (unsigned int icand=0; icand<candList.size(); icand++) delete candList.at(icand);
    for (unsigned int i=0; i<partList.size(); i++) delete partList.at(i);

    SimpleEntry theEntry(ev, varTrue, varReco, ZZMass, weight);
    addByLowest(index, theEntry, false);

    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;

  float firstVal=index.at(0).trueval;
  float lastVal=index.at(index.size()-1).trueval;
  float infimum = std::floor(firstVal);
  float supremum = std::ceil(lastVal);
  cout << "Nentries = " << nEntries << " | truth = " << firstVal << " - " << lastVal << "(" << infimum << ", " << supremum << ")" << endl;

  float divisor=20000;
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
    binning[ix]=(index[ix*ev_stepsize-1].trueval+index[ix*ev_stepsize].trueval)*0.5;
    ExtBin tmpbin;
    tmpbin.binlow = binning[ix-1]; tmpbin.binhigh = binning[ix];
    for (int bin=0; bin<ev_stepsize; bin++){ tmpbin.addEvent(index[(ix-1)*ev_stepsize+bin]); }
    binList.push_back(tmpbin);
    cout << "Boundary (" << ix << ")= " << binning[ix] << " [event " << index[ix*ev_stepsize].id << ", step " << ix*ev_stepsize << "]" << endl;
  }
  ExtBin tmpbin;
  tmpbin.binlow = binning[nbins-1]; tmpbin.binhigh = binning[nbins];
  for (unsigned int bin=(nbins-1)*ev_stepsize; bin<index.size(); bin++){ tmpbin.addEvent(index[bin]); }
  binList.push_back(tmpbin);
  cout << "Boundary (" << nbins << ") = " << binning[nbins] << endl;
  cout << "Bin list has the following bins:" << endl;
  for (unsigned int ib=0; ib<binList.size(); ib++){
    cout << ib << " / " << binList.size() << ": [" << binList.at(ib).binlow << "," << binList.at(ib).binhigh << "]" << endl;
  }

  foutput->cd();

  TProfile* hvar = new TProfile("varTrue", "", nbins, binning); hvar->Sumw2();
  TProfile* havg_varReco = new TProfile("avg_varReco", "", nbins, binning); havg_varReco->Sumw2();
  TH2F* hvarTrue_varReco = new TH2F("varTrue_varReco", "", nbins, binning, nbins, binning); hvarTrue_varReco->Sumw2();

  vector<TGraphErrors*> tglist;

  unsigned int ctr=0;
  for (unsigned int bin=0; bin<binList.size(); bin++){
    cout << "Bin " << bin << " is now being scrutinized..." << endl;

    binList.at(bin).sift(0.99, true);
    binList.at(bin).adjustWeights();

    for (unsigned int ev=0; ev<binList.at(bin).collection.size(); ev++){
      varTrue = binList.at(bin).collection.at(ev).trueval;
      varReco = binList.at(bin).collection.at(ev).recoval;
      ZZMass = binList.at(bin).collection.at(ev).recotrackval;
      weight = binList.at(bin).collection.at(ev).weight;
      havg_varReco->Fill(varTrue, varReco, weight);
      hvar->Fill(varTrue, varTrue, weight);
      hvarTrue_varReco->Fill(varTrue, varReco, weight);
      if (writeFinalTree) newtree->Fill();
    }

    float avgbin_varTrue = hvar->GetBinContent(bin+1);

    vector<pair<float, int>> recominustruevals_entries;
    binList.at(bin).sortByRecoMinusTrueVals(recominustruevals_entries);
    vector<float> recominustruevalbounds;
    unsigned int nrecototal = recominustruevals_entries.size();
    for (unsigned int ev=0; ev<nrecototal; ev++) recominustruevals_entries.at(ev).first = (recominustruevals_entries.at(ev).first+1.)*avgbin_varTrue;

    float recodivisor=1600;
    unsigned int nrecobins = nrecototal/recodivisor;
    const unsigned int nrecobins_th=12;
    while (nrecobins<nrecobins_th){
      if (recodivisor>100) recodivisor -= 100;
      else if (recodivisor>10) recodivisor -= 10;
      else break;
      nrecobins = nrecototal/recodivisor;
    }
    unsigned int recoinc = nrecototal / nrecobins;
    float minval = min(float(0), std::floor(recominustruevals_entries.at(0).first));
    float maxval = max(std::ceil(recominustruevals_entries.at(nrecototal-1).first), float(sqrts*1000.));
    recominustruevalbounds.push_back(minval);
    for (unsigned int i=1; i<nrecobins; i++){
      float a = recominustruevals_entries.at(i*recoinc).first;
      float b = recominustruevals_entries.at(i*recoinc-1).first;
      recominustruevalbounds.push_back((a+b)*0.5);
    }
    recominustruevalbounds.push_back(maxval);
    cout << "Bin " << bin << " has the following reco-true boundaries:";
    for (unsigned int j=0; j<recominustruevalbounds.size(); j++) cout << " " << recominustruevalbounds.at(j);
    cout << endl;

    foutput->cd();

    TH1F* h_RecoMinusTrueVals = new TH1F(Form("RecoMinusTrueVal_Bin%i", bin), "", nrecobins, recominustruevalbounds.data()); h_RecoMinusTrueVals->Sumw2();
    TProfile* p_RecoMinusTrueVals = new TProfile(Form("avg_RecoMinusTrueVal_Bin%i", bin), "", nrecobins, recominustruevalbounds.data()); p_RecoMinusTrueVals->Sumw2();
    for (unsigned int i=0; i<nrecototal; i++){
      float& val = recominustruevals_entries.at(i).first;
      int& ev = recominustruevals_entries.at(i).second;
      float& weight = binList.at(bin).collection.at(ev).weight;
      if (weight==0.) cerr << "Weight is 0!" << endl;
      else if (isnan(weight) || isinf(weight)) cerr << "Weight " << weight << " is nan or inf" << endl;
      if (isnan(val) || isinf(val)) cerr << "Value " << val << " is nan or inf" << endl;
      h_RecoMinusTrueVals->Fill(val, weight);
      p_RecoMinusTrueVals->Fill(val, val, weight);
    }
    const TAxis* axis = h_RecoMinusTrueVals->GetXaxis();
    for (int bin=1; bin<=h_RecoMinusTrueVals->GetNbinsX(); bin++){
      float width = axis->GetBinWidth(bin);
      h_RecoMinusTrueVals->SetBinContent(bin, h_RecoMinusTrueVals->GetBinContent(bin)/width);
      h_RecoMinusTrueVals->SetBinError(bin, h_RecoMinusTrueVals->GetBinError(bin)/width);
    }
    if (debug){
      foutput->WriteTObject(h_RecoMinusTrueVals);
      foutput->WriteTObject(p_RecoMinusTrueVals);
    }

    TGraphErrors* tgbin = makeGraphFromTH1(p_RecoMinusTrueVals, h_RecoMinusTrueVals, Form("tg_%s", h_RecoMinusTrueVals->GetName()));
    tgbin->SetTitle(Form("<m^{true}_{JJ}> = %.3f #pm %.3f GeV, <m^{reco}_{JJ}> = %.3f #pm %.3f GeV", hvar->GetBinContent(bin+1), hvar->GetBinError(bin+1), havg_varReco->GetBinContent(bin+1), havg_varReco->GetBinError(bin+1)));
    tgbin->GetYaxis()->SetTitle("Sum of weights");
    tgbin->GetXaxis()->SetTitle("m^{reco}_{JJ} / m^{true}_{JJ} - 1");
    cout << "Adding point (0,1)" << endl;
    addPoint(tgbin, 0, 1, 0, 0);
    for (int ip=0; ip<tgbin->GetN(); ip++) tgbin->GetY()[ip] = log(tgbin->GetY()[ip]);

    double dfirst, dlast;
    dfirst = (tgbin->GetY()[1]-tgbin->GetY()[0])/(tgbin->GetX()[1]-tgbin->GetX()[0]);
    dlast = (0-tgbin->GetY()[tgbin->GetN()-1])/(maxval-tgbin->GetX()[tgbin->GetN()-1]);
    TSpline3* sp_tgbin = new TSpline3("sp_tgbin", tgbin, "b1e1", dfirst, dlast);
    unsigned int ninsert;
    vector<double> xcoord;
    for (int ip=0; ip<tgbin->GetN(); ip++) xcoord.push_back(tgbin->GetX()[ip]);
    /*
    for (unsigned int ip=0; ip<xcoord.size()-1; ip++){
    const double xwidth = (xcoord[ip+1] - xcoord[ip])/((double)(ninsert+1));
    for (unsigned int j=1; j<=ninsert; j++){
    double xtmp = xwidth*j+xcoord[ip];
    double ytmp = sp_tgbin->Eval(xtmp);
    addPoint(tgbin, xtmp, ytmp, 0, 0);
    }
    }
    */
    unsigned int binzero=0;
    for (unsigned int ip=0; ip<xcoord.size()-1; ip++){
      if (xcoord[ip]<avgbin_varTrue && xcoord[ip+1]>=avgbin_varTrue){
        binzero = ip;
        break;
      }
    }
    for (unsigned int ip=0; ip<xcoord.size()-1; ip++){
      if (ip==binzero || ip==binzero+1 || ip==(xcoord.size()-2)) ninsert=9;
      else if (ip<binzero) ninsert=5;
      else ninsert=2;
      delete sp_tgbin;
      sp_tgbin = new TSpline3("sp_tgbin", tgbin, "b1e1", dfirst, dlast);
      const double xwidth = (xcoord[ip+1] - xcoord[ip])/((double)(ninsert+1));
      for (unsigned int j=1; j<=ninsert; j++){
        double xtmp = xwidth*j+xcoord[ip];
        double ytmp = sp_tgbin->Eval(xtmp);
        cout << "Inserting (x,y) = ( " << xtmp << " , " << ytmp << " )" << endl;
        addPoint(tgbin, xtmp, ytmp, 0, 0);
      }
    }

    for (int ip=0; ip<tgbin->GetN(); ip++) tgbin->GetY()[ip] = exp(tgbin->GetY()[ip]);
    delete sp_tgbin;
    //dfirst = (tgbin->GetY()[1]-tgbin->GetY()[0])/(tgbin->GetX()[1]-tgbin->GetX()[0]);
    dlast = (tgbin->GetY()[tgbin->GetN()-1]-tgbin->GetY()[tgbin->GetN()-2])/(tgbin->GetX()[tgbin->GetN()-1]-tgbin->GetX()[tgbin->GetN()-2]);
    float tgxmax = tgbin->GetX()[tgbin->GetN()-1];
    float tgyatxmax = tgbin->GetY()[tgbin->GetN()-1];
    float tgyatinf = tgyatxmax*0.001;
    float tgZeroPoint = (tgyatinf-tgyatxmax)/dlast + tgxmax;
    float tgLastPoint = min(tgZeroPoint, (float)(sqrts*1000.));
    ninsert=2;
    for (unsigned int j=1; j<=ninsert; j++){
      const double xwidth = (tgLastPoint-tgxmax)/((double)ninsert);
      double xtmp = tgxmax+xwidth*j;
      double ytmp = tgyatxmax + dlast*(xtmp-tgxmax);
      if (j==ninsert){
        xtmp += xwidth;
        xtmp = min(xtmp, (double)(sqrts*1000.));
        cout << "Inserting (x,y) = ( " << xtmp << " , " << ytmp << " )" << endl;
        addPoint(tgbin, xtmp, ytmp, 0, 0);
        for (unsigned int k=0; k<1; k++){
          if (xtmp+xwidth<=sqrts*1000.){
            xtmp += xwidth;
            cout << "Inserting (x,y) = ( " << xtmp << " , " << ytmp << " )" << endl;
            addPoint(tgbin, xtmp, ytmp, 0, 0);
          }
        }
        tgLastPoint=xtmp;
      }
      else{
        cout << "Inserting (x,y) = ( " << xtmp << " , " << ytmp << " )" << endl;
        addPoint(tgbin, xtmp, ytmp, 0, 0);
      }
    }
    ninsert=3;
    for (unsigned int j=1; j<=ninsert; j++){
      const double xwidth = (sqrts*1000.-tgLastPoint)/((double)ninsert);
      double xtmp = tgLastPoint+xwidth*j;
      double ytmp = tgyatinf;
      cout << "Inserting (x,y) = ( " << xtmp << " , " << ytmp << " )" << endl;
      addPoint(tgbin, xtmp, ytmp, 0, 0);
    }
    tgbin->SetMarkerStyle(20);
    tgbin->SetMarkerSize(1.2);
    tgbin->SetLineColor(kBlack);
    tgbin->SetLineWidth(2);

    RooRealVar* reco_minus_true_var = new RooRealVar("plotVar", "", minval, maxval);
    reco_minus_true_var->setBins((int)(maxval-minval));
    MELANCSplineFactory_1D* fac_sp = new MELANCSplineFactory_1D(*reco_minus_true_var, Form("Bin%i", bin));
    fac_sp->setPoints(tgbin);
    MELANCSplineCore* spPDF_bin = fac_sp->getFunc();
    double ymin=9e9, ymax=-9e9;
    //for (int ip=0; ip<tgbin->GetN(); ip++){
    //  reco_minus_true_var->setVal(tgbin->GetX()[ip]);
    //  cout << "Spline value at " << reco_minus_true_var->getVal() << " = " << spPDF_bin->getVal() << " =? " << tgbin->GetY()[ip] << endl;
    //}
    RooRealIntegral* pdf_int = new RooRealIntegral("pdf_int", "", *spPDF_bin, RooArgSet(*reco_minus_true_var));
    //double pdfintval = 1;
    double pdfintval = pdf_int->getVal();
    cout << "PDF integral: " << pdfintval << endl;
    cout << "Sum of weights: " << h_RecoMinusTrueVals->Integral("width") << endl;
    double riemannintval=0;
    for (int ip=0; ip<tgbin->GetN(); ip++){
      tgbin->GetY()[ip] /= pdfintval;
      tgbin->GetEY()[ip] /= pdfintval;
      ymin = min(ymin, tgbin->GetY()[ip]);
      ymax = max(ymax, tgbin->GetY()[ip]);

      if (ip!=tgbin->GetN()-1) riemannintval += (tgbin->GetY()[ip]+tgbin->GetY()[ip-1])*(tgbin->GetX()[ip]-tgbin->GetX()[ip-1])*0.5;
    }
    cout << "Post-normalization Riemann integral: " << riemannintval << endl;
    delete pdf_int;
    fac_sp->setPoints(tgbin);
    spPDF_bin = fac_sp->getFunc();
    //spPDF_bin->doFloor(false);
    //cout << "Analytical integral of spline: " << spPDF_bin->analyticalIntegral(2) << endl;
    pdf_int = new RooRealIntegral("pdf_int", "", *spPDF_bin, RooArgSet(*reco_minus_true_var));
    //cout << "Post-normalization PDF integral: " << pdf_int->getVal() << endl;
    spPDF_bin->forceNumInt(true);
    delete pdf_int;
    pdf_int = new RooRealIntegral("pdf_int", "", *spPDF_bin, RooArgSet(*reco_minus_true_var));
    //cout << "Post-normalization PDF numerical integral: " << pdf_int->getVal() << endl;
    delete pdf_int;

    RooPlot* plot_bin = reco_minus_true_var->frame();
    cout << "Frame created for [" << reco_minus_true_var->getMin() << " , " << reco_minus_true_var->getMax() << " , " << reco_minus_true_var->getBins() << "]" << endl;
    plot_bin->GetXaxis()->CenterTitle();
    plot_bin->GetYaxis()->SetTitleOffset(1.2);
    plot_bin->GetYaxis()->CenterTitle();
    plot_bin->GetXaxis()->SetTitle(tgbin->GetXaxis()->GetTitle());
    plot_bin->GetYaxis()->SetTitle(tgbin->GetYaxis()->GetTitle());
    plot_bin->SetTitle(tgbin->GetTitle()); tgbin->SetTitle("");
    //spPDF_bin->setVerbosity(MELANCSplineCore::kVerbose);
    spPDF_bin->plotOn(plot_bin, LineColor(kRed), LineWidth(2), LineStyle(1));
    plot_bin->GetYaxis()->SetRangeUser(ymin, ymax*1.2);
    plot_bin->GetXaxis()->SetNdivisions(-505);

    TCanvas* canvas = new TCanvas(Form("c_%s", h_RecoMinusTrueVals->GetName()), "", 800, 800);
    plot_bin->Draw();
    tgbin->Draw("e1psame");
    canvas->Modified();
    canvas->Update();
    if (debug) foutput->WriteTObject(canvas);
    canvas->Close();

    delete plot_bin;

    reco_minus_true_var->setRange(-1, reco_minus_true_var->getMax()/avgbin_varTrue-1.);
    for (int ip=0; ip<tgbin->GetN(); ip++){
      tgbin->GetY()[ip] *= avgbin_varTrue;
      tgbin->GetEY()[ip] *= avgbin_varTrue;
      tgbin->GetX()[ip] = tgbin->GetX()[ip]/avgbin_varTrue-1.;
      tgbin->GetEX()[ip] /= avgbin_varTrue;
    }
    tglist.push_back(tgbin);

    delete fac_sp;
    delete reco_minus_true_var;
    delete p_RecoMinusTrueVals;
    delete h_RecoMinusTrueVals;
  }

  RooArgList splineList;
  vector<MELALinearInterpFunc::T> truthList;
  vector<MELANCSplineFactory_1D*> facList;

  RooRealVar finalX("true", "", 0.1, sqrts*1000.);
  finalX.setBins((int)(finalX.getMax()-finalX.getMin())/50);
  RooRealVar finalY("reco", "", finalX.getMin(), finalX.getMax());
  finalY.setBins(finalX.getBins());
  RooFormulaVar finalYoverX("RecoOverTrueVar", "@1/max(0.1, @0)-1.", RooArgList(finalX, finalY));

  for (int is=-1; is<=int(tglist.size()); is++){
    TGraph* tg;
    MELALinearInterpFunc::T tval;
    if (is>=0 && is<int(tglist.size())){
      tval = hvar->GetBinContent(is+1);
      tg = (TGraphErrors*)tglist[is]->Clone(Form("tg_bin%i", is+1));
    }
    else if (is==-1){
      tval=finalX.getMin();
      tg = (TGraphErrors*)tglist[is+1]->Clone(Form("tg_bin%i", is+1));
    }
    else{
      tval = finalX.getMax();
      tg = (TGraphErrors*)tglist[is-1]->Clone(Form("tg_bin%i", is+1));

      TSpline3 sptmp("sptmp", tg, "b2e2", 0, 0);
      double* xy[2];
      int n = 1; // Point at exactly 0
      for (int ip=0; ip<tg->GetN(); ip++){ if (tg->GetX()[ip]<0.) n++; }
      for (unsigned int i=0; i<2; i++) xy[i]=new double[n];
      for (int ip=0; ip<n; ip++){
        if (ip==n-1){ xy[0][ip]=0; xy[1][ip]=sptmp.Eval(xy[0][ip]); }
        else{ xy[0][ip]=tg->GetX()[ip]; xy[1][ip]=tg->GetY()[ip]; }
      }
      delete tg;
      tg = new TGraph(n, xy[0], xy[1]);
      tg->SetName(Form("tg_bin%i", is+1));

      for (unsigned int i=0; i<2; i++) delete[] xy[i];
    }
    tg->SetName(Form("tg_bin%i", is+1));
    tg->SetTitle(Form("Truth at %.2f", tval));
    tg->GetXaxis()->SetTitle("Reco/Truth-1");
    tg->GetYaxis()->SetTitle("Probability");

    MELANCSplineFactory_1D* fac_tmp = new MELANCSplineFactory_1D(finalYoverX, Form("tmp1DFactory_bin%i", is+1));
    fac_tmp->setPoints(tg);
    MELANCSpline_1D_fast* functmp = fac_tmp->getFunc();
    functmp->setRangeValidity(tg->GetX()[0], tg->GetX()[tg->GetN()-1]);

    finalX.setVal(tval);
    RooRealIntegral sp_int = RooRealIntegral("pdftmp_int", "", *functmp, RooArgSet(finalY));
    double pdftmpintval = sp_int.getVal();
    cout << "Bin " << is+1 << " at truth = " << tval << " renormalized by 1/" << pdftmpintval << endl;
    for (int ip=0; ip<tg->GetN(); ip++){ tg->GetY()[ip] /= pdftmpintval; }
    fac_tmp->setPoints(tg);
    functmp = fac_tmp->getFunc();
    functmp->setRangeValidity(tg->GetX()[0], tg->GetX()[tg->GetN()-1]);
    MELAFuncPdf* pdftmp = fac_tmp->getPDF();

    splineList.add(*pdftmp);
    truthList.push_back(tval);
    facList.push_back(fac_tmp);

    foutput->WriteTObject(tg);
    delete tg;
  }

  MELALinearInterpFunc* theFunc = new MELALinearInterpFunc("linearFunc", "", finalX, truthList, splineList);
  theFunc->setRangeValidity(finalX.getMin(), finalX.getMax());
  MELAFuncPdf* thePdf = new MELAFuncPdf("linearPdf", "", *theFunc);
  TCanvas* canvas;
  {
    finalX.setRange(35, 535);
    finalY.setRange(35, 535);
    finalX.setBins(finalX.getMax()-finalX.getMin());
    finalY.setBins(finalY.getMax()-finalY.getMin());
    TH2F* h_frm = new TH2F("hFinalResolution", "", finalX.getBins(), finalX.getMin(), finalX.getMax(), finalY.getBins(), finalY.getMin(), finalY.getMax());
    for (int ix=1; ix<=h_frm->GetNbinsX(); ix++){
      double xval = h_frm->GetXaxis()->GetBinCenter(ix);
      finalX.setVal(xval);
      for (int iy=1; iy<=h_frm->GetNbinsY(); iy++){
        double yval = h_frm->GetYaxis()->GetBinCenter(iy);
        finalY.setVal(yval);
        h_frm->SetBinContent(ix, iy, thePdf->getVal());
      }
    }
    h_frm->SetStats(0);
    h_frm->GetXaxis()->SetTitle("True m_{JJ} (GeV)");
    h_frm->GetYaxis()->SetTitle("Reco. m_{JJ} (GeV)");
    h_frm->SetTitle(TString("m_{JJ} resolution for ")+strprod);
    canvas = new TCanvas(TString("cFinalResolution")+strprod, "", 800, 800);
    h_frm->Draw("colz");
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(Form("%s%s", canvas->GetName(), ".pdf"));
    canvas->Close();
    if (debug) foutput->WriteTObject(h_frm);
    delete h_frm;
  }

  if (strprod=="WH" || strprod=="ZH"){
    RooRealVar* mV = new RooRealVar("mV", "", 0, 1000.*sqrts);
    RooRealVar* GaV = new RooRealVar("GaV", "", 0, 1000.*sqrts);
    if (strprod=="ZH"){
      mV->setVal(mela.getPrimaryMass(23));
      GaV->setVal(mela.getPrimaryWidth(23));
    }
    else{
      mV->setVal(mela.getPrimaryMass(24));
      GaV->setVal(mela.getPrimaryWidth(24));
    }
    mV->setConstant(true);
    GaV->setConstant(true);
    cout << "Constructing a BW on top of the resolution. mV = " << mV->getVal() << ", GammaV = " << GaV->getVal() << endl;

    // First sample the true BW
    RooGenericPdf* trueBW = new RooGenericPdf("trueBW", "2.*@0/(pow(pow(@0,2)-pow(@1,2),2) + pow(@1*@2,2))", RooArgList(finalX, *mV, *GaV));
    vector<float> samplePointsList;

    vector<pair<double, double>> generateRanges;
    vector<pair<int, int>> generateReqStep;
    generateRanges.push_back(pair<double, double>(0.2, 64.9)); generateReqStep.push_back(pair<int, int>(100, 20));
    generateRanges.push_back(pair<double, double>(65., mV->getVal()-0.05)); generateReqStep.push_back(pair<int, int>(100, 10));
    generateRanges.push_back(pair<double, double>(mV->getVal()+0.05, 120.)); generateReqStep.push_back(pair<int, int>(100, 10));
    generateRanges.push_back(pair<double, double>(120.1, 200.)); generateReqStep.push_back(pair<int, int>(25, 40));
    generateRanges.push_back(pair<double, double>(200.1, 600.)); generateReqStep.push_back(pair<int, int>(25, 40));
    generateRanges.push_back(pair<double, double>(600.1, 1500.)); generateReqStep.push_back(pair<int, int>(25, 40));
    generateRanges.push_back(pair<double, double>(1500.1, 4000.)); generateReqStep.push_back(pair<int, int>(25, 40));
    generateRanges.push_back(pair<double, double>(4000.1, 1000.*sqrts-0.1)); generateReqStep.push_back(pair<int, int>(50, 20));
    for (unsigned int ig=0; ig<generateRanges.size(); ig++){
      int& nrequested = generateReqStep[ig].first;
      int& nstep = generateReqStep[ig].second;
      double& rmin = generateRanges[ig].first;
      double& rmax = generateRanges[ig].second;
      finalX.setRange(rmin, rmax);
      RooDataSet* dsSamplePoints = trueBW->generate(finalX, nrequested*nstep);
      vector<float> tmpList;
      for (int ev=0; ev<dsSamplePoints->numEntries(); ev++){
        float val = dynamic_cast<RooAbsReal*>(dsSamplePoints->get(ev)->find(finalX))->getVal();
        addByLowest(tmpList, val, true);
      }
      delete dsSamplePoints;
      for (int ev=0; ev<(int)tmpList.size(); ev+=nstep) samplePointsList.push_back(tmpList.at(ev));
    }
    finalX.setRange(0.1, 1000.*sqrts);
    finalY.setRange(0.1, 1000.*sqrts);
    addByLowest(samplePointsList, float(finalX.getMin()), true);
    addByLowest(samplePointsList, float(finalX.getMax()), true);

    vector<pair<double, double>> idealBWpoints;
    for (unsigned int ev=0; ev<samplePointsList.size(); ev++){
      finalX.setVal(samplePointsList[ev]);
      idealBWpoints.push_back(pair<double, double>(finalX.getVal(), trueBW->getVal()));
      cout << "Sampled point " << idealBWpoints.at(ev).first << " , " << idealBWpoints.at(ev).second << " from ideal BW." << endl;
    }
    vector<pair<double, double>> recoBWpoints;
    typedef MELANCSplineCore::T stype;
    for (unsigned int ev=0; ev<samplePointsList.size(); ev++){
      double recoval = samplePointsList[ev];
      finalY.setVal(recoval);

      vector<pair<stype, stype>> idealBWconvpoints;
      for (unsigned int ip=0; ip<idealBWpoints.size(); ip++){
        double& trueval = idealBWpoints.at(ip).first;
        finalX.setVal(trueval);
        double resoval = thePdf->getVal();
        idealBWconvpoints.push_back(pair<stype, stype>((stype)trueval, (stype)(idealBWpoints.at(ip).second*resoval)));
      }
      MELANCSplineFactory_1D truespline(finalX, "tmp");
      truespline.setPoints(idealBWconvpoints);
      double recobwval = truespline.getPDF()->analyticalIntegral(2);
      recoBWpoints.push_back(pair<double, double>(recoval, recobwval));

      cout << "Conv. ideal BW integral at reco = " << recoval << " : " << recobwval << endl;

      /*
      TGraph* tg_idealBWconv = makeGraphFromPair(idealBWconvpoints, Form("tg_idealBWconv_%i", ev));
      tg_idealBWconv->SetTitle(Form("Conv. ideal BW at reco = %.5f", recoval));
      tg_idealBWconv->SetDrawOption("pl");
      tg_idealBWconv->GetXaxis()->SetRangeUser(0.1, 250.);
      foutput->WriteTObject(tg_idealBWconv);
      delete tg_idealBWconv;
      */
    }
    TGraph* tg_recoBW = makeGraphFromPair(recoBWpoints, "tg_recoBW");
    TGraph* tg_idealBW = makeGraphFromPair(idealBWpoints, "tg_idealBW");

    tg_recoBW->SetMarkerStyle(20);
    tg_recoBW->SetMarkerSize(1.2);
    tg_recoBW->SetMarkerColor(kRed);
    tg_recoBW->SetLineColor(kRed);
    tg_recoBW->GetXaxis()->SetRangeUser(0.1, 250.);
    tg_recoBW->GetXaxis()->SetTitle("m_{JJ} (GeV)");
    tg_recoBW->GetYaxis()->SetTitle("BW (GeV^{-3})");
    tg_recoBW->SetTitle(TString("Reco. BW for ")+strprod);
    tg_idealBW->SetMarkerStyle(21);
    tg_idealBW->SetMarkerSize(0.8);
    tg_idealBW->SetMarkerColor(kBlack);
    tg_idealBW->SetLineColor(kBlack);
    tg_idealBW->GetXaxis()->SetRangeUser(0.1, 250.);
    tg_idealBW->GetXaxis()->SetTitle("m_{JJ} (GeV)");
    tg_idealBW->GetYaxis()->SetTitle("BW (GeV^{-3})");
    tg_idealBW->SetTitle(TString("True BW for ")+strprod);
    foutput->WriteTObject(tg_idealBW);
    foutput->WriteTObject(tg_recoBW);

    canvas = new TCanvas(TString("cFinalBWComparison_")+strprod, "", 800, 800);
    tg_idealBW->SetTitle(TString("BW for ")+strprod);
    tg_recoBW->SetTitle("");
    tg_idealBW->Draw("acp");
    tg_recoBW->Draw("cpsame");
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(Form("%s%s", canvas->GetName(), ".pdf"));
    canvas->Close();

    delete tg_idealBW;
    delete tg_recoBW;

    delete trueBW;
    delete GaV;
    delete mV;
  }
  delete thePdf;
  delete theFunc;

  for (unsigned int is=0; is<facList.size(); is++) delete facList.at(is);
  for (unsigned int is=0; is<tglist.size(); is++) delete tglist.at(is);

  TGraphErrors* tg_avg_varReco = makeGraphFromTH1(hvar, havg_varReco, Form("tg_%s", havg_varReco->GetName()));
  foutput->WriteTObject(tg_avg_varReco);
  delete tg_avg_varReco;

  for (int ix=0; ix<=hvarTrue_varReco->GetNbinsX()+1; ix++){
    double integral = hvarTrue_varReco->Integral(ix, ix, 0, hvarTrue_varReco->GetNbinsY()+1);
    if (integral!=0.){
      for (int iy=0; iy<=hvarTrue_varReco->GetNbinsY()+1; iy++){
        double bincontent = hvarTrue_varReco->GetBinContent(ix, iy);
        double binerror = hvarTrue_varReco->GetBinError(ix, iy);
        bincontent /= integral;
        binerror /= integral; // I know this is wrong, but I don't care about the errors that much.
        hvarTrue_varReco->SetBinContent(ix, iy, bincontent);
        hvarTrue_varReco->SetBinError(ix, iy, binerror);
      }
    }
  }
  if (debug){
    foutput->WriteTObject(hvar);
    foutput->WriteTObject(havg_varReco);
    foutput->WriteTObject(hvarTrue_varReco);
    if (writeFinalTree) foutput->WriteTObject(newtree);
  }
  if (writeFinalTree) delete newtree;
  delete hvarTrue_varReco;
  delete havg_varReco;
  delete hvar;
  delete[] binning;

  for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();
  foutput->Close();
}

void test(){
  vector<pair<MELANCSplineCore::T, MELANCSplineCore::T>> pList;
  MELANCSplineCore::T xrange[2]={ 0, 2 };
  const int n = 11;
  for (int i=0; i<n; i++){
    MELANCSplineCore::T x = (xrange[1]-xrange[0])/((MELANCSplineCore::T)(n-1))*((MELANCSplineCore::T)i)+xrange[0];
    MELANCSplineCore::T y = 1./(xrange[1]-xrange[0]);
    pList.push_back(pair<MELANCSplineCore::T, MELANCSplineCore::T>(x, y));
  }
  RooRealVar var("var", "", xrange[0], xrange[1]);
  var.setBins(2);
  MELANCSplineFactory_1D fac_sp(var, "tmp");
  fac_sp.setPoints(pList);
  MELANCSplineCore* sp = fac_sp.getFunc();

  RooRealIntegral* sp_int = new RooRealIntegral("sp_int", "", *sp, RooArgSet(var));
  double spintval = sp_int->getVal();
  cout << "PDF int: " << spintval << endl;
  delete sp_int;

  RooPlot* plot = var.frame();
  plot->GetXaxis()->CenterTitle();
  plot->GetYaxis()->SetTitleOffset(1.2);
  plot->GetYaxis()->CenterTitle();
  sp->plotOn(plot, LineColor(kRed), LineWidth(2), LineStyle(1));
  plot->GetXaxis()->SetNdivisions(510);
  TCanvas* canvas = new TCanvas("ctest", "", 800, 800);
  plot->Draw();
  canvas->Modified();
  canvas->Update();
  canvas->SaveAs("test.pdf");
  canvas->Close();
  delete plot;
}
