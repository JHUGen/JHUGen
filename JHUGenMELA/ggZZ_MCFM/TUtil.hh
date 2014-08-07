// March 28 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)


#ifndef ZZ_COMMON
#define ZZ_COMMON
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <TFile.h>
#include <TF1.h>
#include "TVar.hh"
#include "TMCFM.hh"
#include "TModHiggsMatEl.hh"
#include "TModGravitonMatEl.hh"
#include "TModZprimeMatEl.hh"
#include "TModHiggsJJMatEl.hh"

using namespace std;
TString DbnEventLepSelName(int i);
void My_choose(TVar::Process process, TVar::Production production, TVar::LeptonInterference leptonInterf, int flavor);
void SetEwkCouplingParameters();
void SetAlphaS(double Q, int mynloop, int mynflav, string mypartons);
void SetMCFMHiggsDecayCouplings(bool useBSM, double Hvvcoupl[SIZE_HVV][2]);
bool My_smalls(double s[][mxpart], int npart);
double SumMatrixElementPDF(TVar::Process procees, TVar::Production production,TVar::MatrixElement myME,mcfm_event_type* mcfm_event,double flavor_msq[][11],double* flux,double EBEAM, double coupling[SIZE_HVV_FREENORM]);
double JHUGenMatEl(TVar::Process process, TVar::Production production, mcfm_event_type* mcfm_event, double MReso, double GaReso, 
		   double Hggcoupl[SIZE_HGG][2], double Hvvcoupl[SIZE_HVV][2], double Zqqcoupl[SIZE_ZQQ][2], double Zvvcoupl[SIZE_ZVV][2],
		   double Gqqcoupl[SIZE_GQQ][2], double Gggcoupl[SIZE_GGG][2], double Gvvcoupl[SIZE_GVV][2]);
double HJJMatEl(TVar::Process process,TVar::Production production, const TLorentzVector p[5], double Hggcoupl[SIZE_HGG][2], double Hvvcoupl[SIZE_HVV_VBF][2], double Hwwcoupl[SIZE_HWW_VBF][2], TVar::VerbosityLevel verb, double EBEAM);
double SumMEPDF(const TLorentzVector p0, const TLorentzVector p1, double flavor_msq[nmsq][nmsq],  TVar::VerbosityLevel verbosity, double EBEAM);

#endif
