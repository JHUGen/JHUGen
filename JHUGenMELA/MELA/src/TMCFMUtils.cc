#include "MELAStreamHelpers.hh"
#include "TMCFMUtils.hh"
#include "TUtilHelpers.hh"
#include "TMath.h"


using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;
using namespace std;
using namespace PDGHelpers;
using namespace TNumericUtil;

namespace TMCFMUtils{
  const std::vector<intQuad_t> MCFMHash_QQVVQQAny = Hash_QQVVQQAny();
}


void TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(int iSel, int jSel, int rSel, int sSel, int order[2]){
  const std::vector<intQuad_t>& hash = MCFMHash_QQVVQQAny;
  bool outFound=false;
  for (intQuad_t const& hashel:hash){
    if (
      !(
      (PDGHelpers::isAnUnknownJet(iSel) || iSel==hashel[0])
      &&
      (PDGHelpers::isAnUnknownJet(jSel) || jSel==hashel[1])
      )
      ) continue;
    // Final particles are q
    if (PDGHelpers::isAJet(rSel) && PDGHelpers::isAJet(sSel)){
      if (
        (PDGHelpers::isAnUnknownJet(rSel) || rSel==hashel[2])
        &&
        (PDGHelpers::isAnUnknownJet(sSel) || sSel==hashel[3])
        ){
        order[0]=0;
        order[1]=1;
        outFound=true;
        //MELAout << "Hash requested outgoing "<< hashel[2] << " " << hashel[3] << ", unswapped r, s = " << rSel << " " << sSel << endl;
      }
      else if (
        (PDGHelpers::isAnUnknownJet(rSel) || rSel==hashel[3])
        &&
        (PDGHelpers::isAnUnknownJet(sSel) || sSel==hashel[2])
        ){
        order[0]=1;
        order[1]=0;
        outFound=true;
        //MELAout << "Hash requested outgoing "<< hashel[2] << " " << hashel[3] << ", swapped r, s = " << rSel << " " << sSel << endl;
      }
    }
    // Final particles l/nu
    else if ((PDGHelpers::isALepton(rSel) || PDGHelpers::isANeutrino(rSel)) && (PDGHelpers::isALepton(sSel) || PDGHelpers::isANeutrino(sSel))){
      if (abs(hashel[0])==abs(hashel[1]) && abs(hashel[0])==abs(hashel[2]) && abs(hashel[0])==abs(hashel[3])) continue; // Do not consider the ordering in uquq_uquq or dqdq_dqdq

      if (
        (
        TMath::Sign(1, rSel)==TMath::Sign(1, hashel[2]) &&
        ((PDGHelpers::isALepton(rSel) && PDGHelpers::isDownTypeQuark(hashel[2])) || (PDGHelpers::isANeutrino(rSel) && PDGHelpers::isUpTypeQuark(hashel[2])))
        )
        &&
        (
        TMath::Sign(1, sSel)==TMath::Sign(1, hashel[3]) &&
        ((PDGHelpers::isALepton(sSel) && PDGHelpers::isDownTypeQuark(hashel[3])) || (PDGHelpers::isANeutrino(sSel) && PDGHelpers::isUpTypeQuark(hashel[3])))
        )
        ){
        order[0]=0;
        order[1]=1;
        outFound=true;
        //MELAout << "Hash requested outgoing "<< hashel[2] << " " << hashel[3] << ", unswapped r, s = " << rSel << " " << sSel << endl;
      }
      else if (
        (
        TMath::Sign(1, rSel)==TMath::Sign(1, hashel[3]) &&
        ((PDGHelpers::isALepton(rSel) && PDGHelpers::isDownTypeQuark(hashel[3])) || (PDGHelpers::isANeutrino(rSel) && PDGHelpers::isUpTypeQuark(hashel[3])))
        )
        &&
        (
        TMath::Sign(1, sSel)==TMath::Sign(1, hashel[2]) &&
        ((PDGHelpers::isALepton(sSel) && PDGHelpers::isDownTypeQuark(hashel[2])) || (PDGHelpers::isANeutrino(sSel) && PDGHelpers::isUpTypeQuark(hashel[2])))
        )
        ){
        order[0]=1;
        order[1]=0;
        outFound=true;
        //MELAout << "Hash requested outgoing "<< hashel[2] << " " << hashel[3] << ", swapped r, s = " << rSel << " " << sSel << endl;
      }
      if (PDGHelpers::getCoupledVertex(rSel, sSel)!=PDGHelpers::getCoupledVertex(hashel[2], hashel[3])) outFound=false;
    }
    if (outFound) break;
  }
  if (!outFound){ for (unsigned int ip=0; ip<2; ip++) order[ip]=-1; }
}
std::vector<intQuad_t> TMCFMUtils::Hash_QQVVQQAny(){
  std::vector<intQuad_t> pcfg;
  std::vector<intQuad_t> hash_qqvvqq = TMCFMUtils::Hash_QQVVQQ();
  std::vector<intQuad_t> hash_qqvvqqstrong = TMCFMUtils::Hash_QQVVQQStrong();
  TUtilHelpers::copyVector(hash_qqvvqq, pcfg);
  TUtilHelpers::copyVector(hash_qqvvqqstrong, pcfg);
  return pcfg;
}
std::vector<intQuad_t> TMCFMUtils::Hash_QQVVQQ(){
  /*
  Based on the following cases in MCFM:
  parameter(
  & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
  & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6,
  & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9
  & jmax=12)
  integer,parameter:: j1(jmax)=(/1,2,8,8, 7,2,7,1, 1,7,2,7/)
  integer,parameter:: j2(jmax)=(/2,1,7,7, 2,7,1,7, 7,1,7,2/)
  integer,parameter:: j7(jmax)=(/7,7,2,1, 1,8,2,8, 2,8,1,8/)
  integer,parameter:: j8(jmax)=(/8,8,1,2, 8,1,8,2, 8,2,8,1/)
  */
  std::vector<intQuad_t> base_cfg; base_cfg.reserve(17);
  // uc_uc
  base_cfg.emplace_back(2, 4, 2, 4);
  // ds_ds
  base_cfg.emplace_back(1, 3, 1, 3);
  base_cfg.emplace_back(1, 5, 1, 5);
  base_cfg.emplace_back(3, 5, 3, 5);
  // ub_ub
  base_cfg.emplace_back(2, 3, 2, 3);
  base_cfg.emplace_back(2, 5, 2, 5);
  base_cfg.emplace_back(4, 5, 4, 5);
  // dc_dc
  base_cfg.emplace_back(1, 4, 1, 4);
  // du_du
  base_cfg.emplace_back(1, 2, 1, 2);
  base_cfg.emplace_back(3, 4, 3, 4);
  // dc_us
  base_cfg.emplace_back(1, 4, 2, 3);
  // us_dc
  base_cfg.emplace_back(2, 3, 1, 4);
  // uu_uu
  base_cfg.emplace_back(2);
  base_cfg.emplace_back(4);
  // dd_dd
  base_cfg.emplace_back(1);
  base_cfg.emplace_back(3);
  base_cfg.emplace_back(5);

  std::vector<intQuad_t> jcfg; jcfg.reserve(12);
  jcfg.emplace_back(0, 1, 2, 3);
  jcfg.emplace_back(1, 0, 2, 3);
  jcfg.emplace_back(3, 2, 1, 0);
  jcfg.emplace_back(3, 2, 0, 1);
  jcfg.emplace_back(2, 1, 0, 3);
  jcfg.emplace_back(1, 2, 3, 0);
  jcfg.emplace_back(2, 0, 1, 3);
  jcfg.emplace_back(0, 2, 3, 1);
  jcfg.emplace_back(0, 2, 1, 3);
  jcfg.emplace_back(2, 0, 3, 1);
  jcfg.emplace_back(1, 2, 0, 3);
  jcfg.emplace_back(2, 1, 3, 0);

  std::vector<intQuad_t> pcfg;
  for (unsigned int j=0; j<jcfg.size(); j++){
    for (unsigned int p=0; p<base_cfg.size(); p++){
      intQuad_t cfg;
      for (unsigned int ipos=0; ipos<4; ipos++){
        int idpos = jcfg.at(j)[ipos];
        int idAssigned = base_cfg.at(p)[ipos];
        if ((idpos<2 && ipos>=2) || (idpos>=2 && ipos<2)) idAssigned = -idAssigned;
        cfg[jcfg.at(j)[ipos]] = idAssigned;
      }
      if (!TUtilHelpers::checkElementExists(cfg, pcfg)) pcfg.push_back(cfg);
    }
  }
  // Uncommenting the lines below prints out the hash when the library is loaded.
  /*
  for (unsigned int ic=0; ic<pcfg.size(); ic++) MELAout
  << "TMCFMUtils::Hash_QQVVQQAny: Hash configuration " << ic << " requests ids=( "
  << pcfg.at(ic)[0] << ", " << pcfg.at(ic)[1] << ", " << pcfg.at(ic)[2] << ", " << pcfg.at(ic)[3]
  << ")" << endl;
  */
  return pcfg;
}
std::vector<intQuad_t> TMCFMUtils::Hash_QQVVQQStrong(){
  /*
  Based on the following cases in MCFM:
  call qq4lggampf(1,2,3,4,5,6,7,8,3,4,za,zb,msqgg)
  msq(1,-1)=msq(1,-1)+stat*aveqq*msqgg(1)
  call qq4lggampf(2,1,3,4,5,6,7,8,3,4,za,zb,msqgg)
  msq(-1,1)=msq(-1,1)+stat*aveqq*msqgg(1)
  call qq4lggampf(7,8,3,4,5,6,1,2,3,4,za,zb,msqgg)
  msq(0,0)=msq(0,0)+avegg*(3d0*msqgg(1)+2d0*msqgg(2))
  call qq4lggampf(7,2,3,4,5,6,1,8,3,4,za,zb,msqgg)
  msq(0,-1)=msq(0,-1)+aveqg*msqgg(1)
  call qq4lggampf(7,1,3,4,5,6,2,8,3,4,za,zb,msqgg)
  msq(-1,0)=msq(-1,0)+aveqg*msqgg(1)
  call qq4lggampf(2,8,3,4,5,6,1,7,3,4,za,zb,msqgg)
  msq(0,1)=msq(0,1)+aveqg*msqgg(1)
  call qq4lggampf(1,8,3,4,5,6,2,7,3,4,za,zb,msqgg)
  msq(1,0)=msq(1,0)+aveqg*msqgg(1)
  */
  std::vector<intQuad_t> base_cfg; base_cfg.reserve(5);
  // Start with qqb_gg
  for (int iq=1; iq<=5; iq++) base_cfg.emplace_back(iq, -iq, 21, 21);

  std::vector<intQuad_t> jcfg; jcfg.reserve(7);
  jcfg.emplace_back(0, 1, 2, 3); // qqb->gg
  jcfg.emplace_back(1, 0, 2, 3); // qbq->gg
  jcfg.emplace_back(2, 3, 0, 1); // gg->qbq
  jcfg.emplace_back(2, 1, 0, 3); // gqb->qbg
  jcfg.emplace_back(2, 0, 1, 3); // qbg->qbg
  jcfg.emplace_back(1, 3, 2, 0); // gq->gq
  jcfg.emplace_back(0, 3, 2, 1); // qg->gq

  std::vector<intQuad_t> pcfg;
  for (unsigned int j=0; j<jcfg.size(); j++){
    for (unsigned int p=0; p<base_cfg.size(); p++){
      intQuad_t cfg;
      for (unsigned int ipos=0; ipos<4; ipos++){
        int idpos = jcfg.at(j)[ipos];
        int idAssigned = base_cfg.at(p)[ipos];
        if (((idpos<2 && ipos>=2) || (idpos>=2 && ipos<2)) && !PDGHelpers::isAGluon(idAssigned)) idAssigned = -idAssigned;
        cfg[jcfg.at(j)[ipos]] = idAssigned;
      }
      if (!TUtilHelpers::checkElementExists(cfg, pcfg)) pcfg.push_back(cfg);
    }
  }
  // Uncommenting the lines below prints out the hash when the library is loaded.
  /*
  for (unsigned int ic=0; ic<pcfg.size(); ic++) MELAout
  << "TMCFMUtils::Hash_QQVVQQAny: Hash configuration " << ic << " requests ids=( "
  << pcfg.at(ic)[0] << ", " << pcfg.at(ic)[1] << ", " << pcfg.at(ic)[2] << ", " << pcfg.at(ic)[3]
  << ")" << endl;
  */
  return pcfg;
}
