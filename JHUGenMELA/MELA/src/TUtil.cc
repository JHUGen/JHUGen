#include <iostream>
#include <cstdio>
#include <cmath>
#include <iterator>
#include <utility>
#include <algorithm>
#include <cassert>
#include "MELAStreamHelpers.hh"
#include "TJHUGenUtils.hh"
#include "TUtilHelpers.hh"
#include "TUtil.hh"
#include "TMath.h"
#include "TLorentzRotation.h"


using namespace std;
using TVar::event_scales_type;
using TVar::simple_event_record;
using MELAStreamHelpers::MELAout;
using MELAStreamHelpers::MELAerr;
using namespace TJHUGenUtils;


namespace TUtil{
  bool forbidMassiveLeptons = true;
  bool forbidMassiveJets = true;
  TVar::FermionMassRemoval LeptonMassScheme = TVar::ConserveDifermionMass;
  TVar::FermionMassRemoval JetMassScheme = TVar::ConserveDifermionMass;
}

/***************************************************/
/***** Scripts for decay and production angles *****/
/***************************************************/

void TUtil::applyLeptonMassCorrection(bool flag){ TUtil::forbidMassiveLeptons = flag; }
void TUtil::applyJetMassCorrection(bool flag){ TUtil::forbidMassiveJets = flag; }
void TUtil::setLeptonMassScheme(TVar::FermionMassRemoval scheme){ LeptonMassScheme=scheme; }
void TUtil::setJetMassScheme(TVar::FermionMassRemoval scheme){ JetMassScheme=scheme; }
void TUtil::constrainedRemovePairMass(TLorentzVector& p1, TLorentzVector& p2, double m1, double m2){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  TLorentzVector p1hat, p2hat;
  if (p1==nullFourVector || p2==nullFourVector) return;

  /***** shiftMass in C++ *****/
  TLorentzVector p12=p1+p2;
  TLorentzVector diffp2p1=p2-p1;
  double p1sq = p1.M2();
  double p2sq = p2.M2();
  double p1p2 = p1.Dot(p2);
  double m1sq = m1*fabs(m1);
  double m2sq = m2*fabs(m2);
  double p12sq = p12.M2();

  TLorentzVector avec=(p1sq*p2 - p2sq*p1 + p1p2*diffp2p1);
  double a = avec.M2();
  double b = (p12sq + m2sq - m1sq) * (pow(p1p2, 2) - p1sq*p2sq);
  double c = pow((p12sq + m2sq - m1sq), 2)*p1sq/4. - pow((p1sq + p1p2), 2)*m2sq;
  double eta =  (-b - sqrt(fabs(b*b -4.*a*c)))/(2.*a);
  double xi = (p12sq + m2sq - m1sq - 2.*eta*(p2sq + p1p2))/(2.*(p1sq + p1p2));

  p1hat = (1.-xi)*p1 + (1.-eta)*p2;
  p2hat = xi*p1 + eta*p2;
  p1=p1hat;
  p2=p2hat;
}
void TUtil::scaleMomentumToEnergy(const TLorentzVector&  massiveJet, TLorentzVector& masslessJet, double mass){
  double energy, p3, newp3, ratio;
  energy = massiveJet.T();
  p3 = massiveJet.P();
  newp3 = sqrt(max(pow(energy, 2)-mass*fabs(mass), 0.));
  ratio = (p3>0. ? (newp3/p3) : 1.);
  masslessJet.SetXYZT(massiveJet.X()*ratio, massiveJet.Y()*ratio, massiveJet.Z()*ratio, energy);
}
pair<TLorentzVector, TLorentzVector> TUtil::removeMassFromPair(
  TLorentzVector jet1, int jet1Id,
  TLorentzVector jet2, int jet2Id,
  double m1, double m2
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  TLorentzVector jet1massless(0, 0, 0, 0), jet2massless(0, 0, 0, 0);

  if (TUtil::forbidMassiveJets && (PDGHelpers::isAJet(jet1Id) || PDGHelpers::isAJet(jet2Id))){
    if (JetMassScheme==TVar::NoRemoval){
      jet1massless=jet1;
      jet2massless=jet2;
    }
    else if (jet1==nullFourVector || jet2==nullFourVector || jet1==jet2 || JetMassScheme==TVar::MomentumToEnergy){
      TUtil::scaleMomentumToEnergy(jet1, jet1massless, m1);
      TUtil::scaleMomentumToEnergy(jet2, jet2massless, m2);
    }
    else if (JetMassScheme==TVar::ConserveDifermionMass){
      jet1massless=jet1;
      jet2massless=jet2;
      TUtil::constrainedRemovePairMass(jet1massless, jet2massless, m1, m2);
    }
    else{
      jet1massless=jet1;
      jet2massless=jet2;
    }
  }
  else if (TUtil::forbidMassiveLeptons && (PDGHelpers::isALepton(jet1Id) || PDGHelpers::isANeutrino(jet1Id) || PDGHelpers::isALepton(jet2Id) || PDGHelpers::isANeutrino(jet2Id))){
    if (LeptonMassScheme==TVar::NoRemoval){
      jet1massless=jet1;
      jet2massless=jet2;
    }
    else if (jet1==nullFourVector || jet2==nullFourVector || jet1==jet2 || LeptonMassScheme==TVar::MomentumToEnergy){
      TUtil::scaleMomentumToEnergy(jet1, jet1massless, m1);
      TUtil::scaleMomentumToEnergy(jet2, jet2massless, m2);
    }
    else if (LeptonMassScheme==TVar::ConserveDifermionMass){
      jet1massless=jet1;
      jet2massless=jet2;
      TUtil::constrainedRemovePairMass(jet1massless, jet2massless, m1, m2);
    }
    else{
      jet1massless=jet1;
      jet2massless=jet2;
    }
  }
  else{
    jet1massless=jet1;
    jet2massless=jet2;
  }

  pair<TLorentzVector, TLorentzVector> result(jet1massless, jet2massless);
  return result;
}
void TUtil::adjustTopDaughters(SimpleParticleCollection_t& daughters){ // Daughters are arranged as b, Wf, Wfb
  if (daughters.size()!=3) return; // Cannot work if the number of daughters is not exactly 3.
  pair<TLorentzVector, TLorentzVector> corrWpair = TUtil::removeMassFromPair(
    daughters.at(1).second, daughters.at(1).first,
    daughters.at(2).second, daughters.at(2).first
    );
  daughters.at(1).second = corrWpair.first;
  daughters.at(2).second = corrWpair.second;
  TLorentzVector pW = daughters.at(1).second+daughters.at(2).second;
  TVector3 pW_boost_old = -pW.BoostVector();
  pair<TLorentzVector, TLorentzVector> corrbW = TUtil::removeMassFromPair(
    daughters.at(0).second, daughters.at(0).first,
    pW, -daughters.at(0).first, // Trick the function
    0., pW.M() // Conserve W mass, ensures Wf+Wfb=W after re-boosting Wf and Wfb to te new W frame.
    );
  daughters.at(0).second=corrbW.first;
  pW=corrbW.second;
  TVector3 pW_boost_new = pW.BoostVector();
  for (unsigned int idau=1; idau<daughters.size(); idau++){
    daughters.at(idau).second.Boost(pW_boost_old);
    daughters.at(idau).second.Boost(pW_boost_new);
  }
}
// Compute a fake jet
void TUtil::computeFakeJet(TLorentzVector realJet, TLorentzVector others, TLorentzVector& fakeJet){
  TLorentzVector masslessRealJet(0, 0, 0, 0);
  if (TUtil::forbidMassiveJets) TUtil::scaleMomentumToEnergy(realJet, masslessRealJet);
  else masslessRealJet = realJet;
  fakeJet = others + masslessRealJet;
  fakeJet.SetVect(-fakeJet.Vect());
  fakeJet.SetE(fakeJet.P());
}

/***** Complex Boost *****/
std::pair<TLorentzVector, TLorentzVector> TUtil::ComplexBoost(TVector3 beta, TLorentzVector p4){
  double bx=beta.X();
  double by=beta.Y();
  double bz=beta.Z();

  double bp = bx*p4.X() + by*p4.Y() + bz*p4.Z();
  double b2 = bx*bx + by*by + bz*bz;
  double gammasqinv = 1.-b2;

  double gamma=0.;
  double gammap_real=0;
  double gammap_imag=0;
  TLorentzVector p4new_real(0, 0, 0, 0), p4new_imag(0, 0, 0, 0);
  if (gammasqinv>0.){
    gamma = 1./sqrt(gammasqinv);
    if (b2>0.) gammap_real = (gamma-1.)/b2;

    p4new_real.SetX(p4.X() + gammap_real*bp*bx + gamma*bx*p4.T());
    p4new_real.SetY(p4.Y() + gammap_real*bp*by + gamma*by*p4.T());
    p4new_real.SetZ(p4.Z() + gammap_real*bp*bz + gamma*bz*p4.T());
    p4new_real.SetT(gamma*(p4.T() + bp));
  }
  else if (gammasqinv<0.){
    gamma = -1./sqrt(-gammasqinv);
    if (b2>0.){
      gammap_real = -1./b2;
      gammap_imag = gamma/b2;
    }

    p4new_real.SetX(p4.X() + gammap_real*bp*bx);
    p4new_real.SetY(p4.Y() + gammap_real*bp*by);
    p4new_real.SetZ(p4.Z() + gammap_real*bp*bz);
    p4new_real.SetT(0.);
    p4new_imag.SetX(gammap_imag*bp*bx + gamma*bx*p4.T());
    p4new_imag.SetY(gammap_imag*bp*by + gamma*by*p4.T());
    p4new_imag.SetZ(gammap_imag*bp*bz + gamma*bz*p4.T());
    p4new_imag.SetT(gamma*(p4.T() + bp));
  }

  return (pair<TLorentzVector, TLorentzVector>(p4new_real, p4new_imag));
}

/***** Decay angles *****/
void TUtil::computeAngles(
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (p4M12==nullFourVector || p4M22==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f13Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M21, Z2_lept1Id);
    p4M11 = f13Pair.first;
    p4M21 = f13Pair.second;
  }
  else if (p4M11==nullFourVector || p4M21==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f24Pair = TUtil::removeMassFromPair(p4M12, Z1_lept2Id, p4M22, Z2_lept2Id);
    p4M12 = f24Pair.first;
    p4M22 = f24Pair.second;
  }
  else{
    pair<TLorentzVector, TLorentzVector> f12Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M12, Z1_lept2Id);
    pair<TLorentzVector, TLorentzVector> f34Pair = TUtil::removeMassFromPair(p4M21, Z2_lept1Id, p4M22, Z2_lept2Id);
    p4M11 = f12Pair.first;
    p4M12 = f12Pair.second;
    p4M21 = f34Pair.first;
    p4M22 = f34Pair.second;
  }

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;

  // Sort Z1 leptons so that:
  if (
    Z1_lept2Id!=-9000
    &&
    (
    (Z1_lept1Id*Z1_lept2Id<0 && Z1_lept1Id<0) // for OS pairs: lep1 must be the negative one
    ||
    ((Z1_lept1Id*Z1_lept2Id>0 || (Z1_lept1Id==0 && Z1_lept2Id==0)) && p4M11.Phi()<=p4M12.Phi()) //for SS pairs: use random deterministic convention
    )
    ) swap(p4M11, p4M12);
  // Same for Z2 leptons
  if (
    Z2_lept2Id!=-9000
    &&
    (
    (Z2_lept1Id*Z2_lept2Id<0 && Z2_lept1Id<0)
    ||
    ((Z2_lept1Id*Z2_lept2Id>0 || (Z2_lept1Id==0 && Z2_lept2Id==0)) && p4M21.Phi()<=p4M22.Phi())
    )
    ) swap(p4M21, p4M22);

  // BEGIN THE CALCULATION

  // build H 4-vectors
  TLorentzVector p4H = p4Z1 + p4Z2;

  // -----------------------------------

  //// costhetastar
  TVector3 boostX = -(p4H.BoostVector());
  TLorentzVector thep4Z1inXFrame(p4Z1);
  TLorentzVector thep4Z2inXFrame(p4Z2);
  thep4Z1inXFrame.Boost(boostX);
  thep4Z2inXFrame.Boost(boostX);
  TVector3 theZ1X_p3 = TVector3(thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z());
  TVector3 theZ2X_p3 = TVector3(thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z());
  costhetastar = theZ1X_p3.CosTheta();

  TVector3 boostV1(0, 0, 0);
  TVector3 boostV2(0, 0, 0);
  //// --------------------------- costheta1
  if (!(fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==22 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==22)){
    boostV1 = -(p4Z1.BoostVector());
    if (boostV1.Mag()>=1.) {
      MELAout << "Warning: Mela::computeAngles: Z1 boost with beta=1, scaling down" << endl;
      boostV1*=0.9999/boostV1.Mag();
    }
    TLorentzVector p4M11_BV1(p4M11);
    TLorentzVector p4M12_BV1(p4M12);
    TLorentzVector p4M21_BV1(p4M21);
    TLorentzVector p4M22_BV1(p4M22);
    p4M11_BV1.Boost(boostV1);
    p4M12_BV1.Boost(boostV1);
    p4M21_BV1.Boost(boostV1);
    p4M22_BV1.Boost(boostV1);

    TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
    //// costheta1
    costheta1 = -p4V2_BV1.Vect().Unit().Dot(p4M11_BV1.Vect().Unit());
  }
  else costheta1 = 0;

  //// --------------------------- costheta2
  if (!(fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==22 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==22)){
    boostV2 = -(p4Z2.BoostVector());
    if (boostV2.Mag()>=1.) {
      MELAout << "Warning: Mela::computeAngles: Z2 boost with beta=1, scaling down" << endl;
      boostV2*=0.9999/boostV2.Mag();
    }
    TLorentzVector p4M11_BV2(p4M11);
    TLorentzVector p4M12_BV2(p4M12);
    TLorentzVector p4M21_BV2(p4M21);
    TLorentzVector p4M22_BV2(p4M22);
    p4M11_BV2.Boost(boostV2);
    p4M12_BV2.Boost(boostV2);
    p4M21_BV2.Boost(boostV2);
    p4M22_BV2.Boost(boostV2);

    TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
    //// costheta2
    costheta2 = -p4V1_BV2.Vect().Unit().Dot(p4M21_BV2.Vect().Unit());
  }
  else costheta2 = 0;

  //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
  TLorentzVector p4M11_BX(p4M11);
  TLorentzVector p4M12_BX(p4M12);
  TLorentzVector p4M21_BX(p4M21);
  TLorentzVector p4M22_BX(p4M22);

  p4M11_BX.Boost(boostX);
  p4M12_BX.Boost(boostX);
  p4M21_BX.Boost(boostX);
  p4M22_BX.Boost(boostX);
  TLorentzVector p4V1_BX = p4M11_BX + p4M12_BX;

  TVector3 beamAxis(0, 0, 1);
  TVector3 p3V1_BX = p4V1_BX.Vect().Unit();
  TVector3 normal1_BX = (p4M11_BX.Vect().Cross(p4M12_BX.Vect())).Unit();
  TVector3 normal2_BX = (p4M21_BX.Vect().Cross(p4M22_BX.Vect())).Unit();
  TVector3 normalSC_BX = (beamAxis.Cross(p3V1_BX)).Unit();


  //// Phi
  float tmpSgnPhi = p3V1_BX.Dot(normal1_BX.Cross(normal2_BX));
  float sgnPhi = 0;
  if (fabs(tmpSgnPhi)>0.) sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  float dot_BX12 = normal1_BX.Dot(normal2_BX);
  if (fabs(dot_BX12)>=1.) dot_BX12 *= 1./fabs(dot_BX12);
  Phi = sgnPhi * acos(-1.*dot_BX12);


  //// Phi1
  float tmpSgnPhi1 = p3V1_BX.Dot(normal1_BX.Cross(normalSC_BX));
  float sgnPhi1 = 0;
  if (fabs(tmpSgnPhi1)>0.) sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  float dot_BX1SC = normal1_BX.Dot(normalSC_BX);
  if (fabs(dot_BX1SC)>=1.) dot_BX1SC *= 1./fabs(dot_BX1SC);
  Phi1 = sgnPhi1 * acos(dot_BX1SC);

  if (isnan(costhetastar) || isnan(costheta1) || isnan(costheta2) || isnan(Phi) || isnan(Phi1)){
    MELAout << "WARNING: NaN in computeAngles: "
      << costhetastar << " "
      << costheta1  << " "
      << costheta2  << " "
      << Phi  << " "
      << Phi1  << " " << endl;
    MELAout << "   boostV1: " <<boostV1.Pt() << " " << boostV1.Eta() << " " << boostV1.Phi() << " " << boostV1.Mag() << endl;
    MELAout << "   boostV2: " <<boostV2.Pt() << " " << boostV2.Eta() << " " << boostV2.Phi() << " " << boostV2.Mag() << endl;
  }
}
void TUtil::computeAnglesCS(
  float pbeam,
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (p4M12==nullFourVector || p4M22==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f13Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M21, Z2_lept1Id);
    p4M11 = f13Pair.first;
    p4M21 = f13Pair.second;
  }
  else if (p4M11==nullFourVector || p4M21==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f24Pair = TUtil::removeMassFromPair(p4M12, Z1_lept2Id, p4M22, Z2_lept2Id);
    p4M12 = f24Pair.first;
    p4M22 = f24Pair.second;
  }
  else{
    pair<TLorentzVector, TLorentzVector> f12Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M12, Z1_lept2Id);
    pair<TLorentzVector, TLorentzVector> f34Pair = TUtil::removeMassFromPair(p4M21, Z2_lept1Id, p4M22, Z2_lept2Id);
    p4M11 = f12Pair.first;
    p4M12 = f12Pair.second;
    p4M21 = f34Pair.first;
    p4M22 = f34Pair.second;
  }

  TVector3 LabXaxis(1.0, 0.0, 0.0);
  TVector3 LabYaxis(0.0, 1.0, 0.0);
  TVector3 LabZaxis(0.0, 0.0, 1.0);

  float Mprot = 0.938;
  float Ebeam = sqrt(pbeam*pbeam + Mprot*Mprot);

  TLorentzVector targ(0., 0., -pbeam, Ebeam);
  TLorentzVector beam(0., 0., pbeam, Ebeam);

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;

  // Sort Z1 leptons so that:
  if (
    Z1_lept2Id!=9000
    &&
    (
    (Z1_lept1Id*Z1_lept2Id<0 && Z1_lept1Id<0) // for OS pairs: lep1 must be the negative one
    ||
    ((Z1_lept1Id*Z1_lept2Id>0 || (Z1_lept1Id==0 && Z1_lept2Id==0)) && p4M11.Phi()<=p4M12.Phi()) //for SS pairs: use random deterministic convention
    )
    ) swap(p4M11, p4M12);
  // Same for Z2 leptons
  if (
    Z2_lept2Id!=9000
    &&
    (
    (Z2_lept1Id*Z2_lept2Id<0 && Z2_lept1Id<0)
    ||
    ((Z2_lept1Id*Z2_lept2Id>0 || (Z2_lept1Id==0 && Z2_lept2Id==0)) && p4M21.Phi()<=p4M22.Phi())
    )
    ) swap(p4M21, p4M22);


  //build H 4-vectors
  TLorentzVector p4H = p4Z1 + p4Z2;
  TVector3 boostX = -(p4H.BoostVector());

  /////////////////////////////
  // Collin-Sopper calculation:
  // in the CS frame, the z-axis is along the bisectrice of one beam and the opposite of the other beam,
  // after their boost in X
  ///////////////////////////////
  // Rotation for the CS Frame

  TRotation rotationCS;

  TLorentzVector beaminX(beam);
  TLorentzVector targinX(targ);
  targinX.Boost(boostX);
  beaminX.Boost(boostX);

  //Bisectrice: sum of unit vectors (remember: you need to invert one beam vector)
  TVector3 beam_targ_bisecinX((beaminX.Vect().Unit() - targinX.Vect().Unit()).Unit());

  // Define a rotationCS Matrix, with Z along the bisectric, 
  TVector3 newZaxisCS(beam_targ_bisecinX.Unit());
  TVector3 newYaxisCS(beaminX.Vect().Unit().Cross(newZaxisCS).Unit());
  TVector3 newXaxisCS(newYaxisCS.Unit().Cross(newZaxisCS).Unit());
  rotationCS.RotateAxes(newXaxisCS, newYaxisCS, newZaxisCS);
  rotationCS.Invert();

  //// costhetastar
  TLorentzVector thep4Z1inXFrame_rotCS(p4Z1);
  TLorentzVector thep4Z2inXFrame_rotCS(p4Z2);
  thep4Z1inXFrame_rotCS.Transform(rotationCS);
  thep4Z2inXFrame_rotCS.Transform(rotationCS);
  thep4Z1inXFrame_rotCS.Boost(boostX);
  thep4Z2inXFrame_rotCS.Boost(boostX);
  TVector3 theZ1XrotCS_p3 = TVector3(thep4Z1inXFrame_rotCS.X(), thep4Z1inXFrame_rotCS.Y(), thep4Z1inXFrame_rotCS.Z());
  costhetastar = theZ1XrotCS_p3.CosTheta();

  //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
  //    TVector3 boostX = -(thep4H.BoostVector());
  TLorentzVector p4M11_BX_rotCS(p4M11);
  TLorentzVector p4M12_BX_rotCS(p4M12);
  TLorentzVector p4M21_BX_rotCS(p4M21);
  TLorentzVector p4M22_BX_rotCS(p4M22);
  p4M11_BX_rotCS.Transform(rotationCS);
  p4M12_BX_rotCS.Transform(rotationCS);
  p4M21_BX_rotCS.Transform(rotationCS);
  p4M22_BX_rotCS.Transform(rotationCS);
  p4M11_BX_rotCS.Boost(boostX);
  p4M12_BX_rotCS.Boost(boostX);
  p4M21_BX_rotCS.Boost(boostX);
  p4M22_BX_rotCS.Boost(boostX);

  TLorentzVector p4Z1_BX_rotCS = p4M11_BX_rotCS + p4M12_BX_rotCS;
  TVector3 p3V1_BX_rotCS = (p4Z1_BX_rotCS.Vect()).Unit();
  TVector3 beamAxis(0, 0, 1);
  TVector3 normal1_BX_rotCS = (p4M11_BX_rotCS.Vect().Cross(p4M12_BX_rotCS.Vect())).Unit();
  TVector3 normal2_BX_rotCS = (p4M21_BX_rotCS.Vect().Cross(p4M22_BX_rotCS.Vect())).Unit();
  TVector3 normalSC_BX_rotCS = (beamAxis.Cross(p3V1_BX_rotCS)).Unit();

  //// Phi
  float tmpSgnPhi = p3V1_BX_rotCS.Dot(normal1_BX_rotCS.Cross(normal2_BX_rotCS));
  float sgnPhi = 0;
  if (fabs(tmpSgnPhi)>0.) sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  float dot_BX12 = normal1_BX_rotCS.Dot(normal2_BX_rotCS);
  if (fabs(dot_BX12)>=1.) dot_BX12 *= 1./fabs(dot_BX12);
  Phi = sgnPhi * acos(-1.*dot_BX12);

  //// Phi1
  float tmpSgnPhi1 = p3V1_BX_rotCS.Dot(normal1_BX_rotCS.Cross(normalSC_BX_rotCS));
  float sgnPhi1 = 0;
  if (fabs(tmpSgnPhi1)>0.) sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  float dot_BX1SC = normal1_BX_rotCS.Dot(normalSC_BX_rotCS);
  if (fabs(dot_BX1SC)>=1.) dot_BX1SC *= 1./fabs(dot_BX1SC);
  Phi1 = sgnPhi1 * acos(dot_BX1SC);

  //// --------------------------- costheta1
  if (!(fabs(Z1_lept1Id)==21 || fabs(Z1_lept1Id)==22 || fabs(Z1_lept2Id)==21 || fabs(Z1_lept2Id)==22)){
    //define Z1 rotation 
    TRotation rotationZ1;
    TVector3 newZaxisZ1(thep4Z1inXFrame_rotCS.Vect().Unit());
    TVector3 newXaxisZ1(newYaxisCS.Cross(newZaxisZ1).Unit());
    TVector3 newYaxisZ1(newZaxisZ1.Cross(newXaxisZ1).Unit());
    rotationZ1.RotateAxes(newXaxisZ1, newYaxisZ1, newZaxisZ1);
    rotationZ1.Invert();

    TLorentzVector thep4Z1inXFrame_rotCS_rotZ1(thep4Z1inXFrame_rotCS);
    thep4Z1inXFrame_rotCS_rotZ1.Transform(rotationZ1);
    TVector3 boostZ1inX_rotCS_rotZ1= -(thep4Z1inXFrame_rotCS_rotZ1.BoostVector());

    TLorentzVector p4M11_BX_rotCS_rotZ1(p4M11_BX_rotCS);
    TLorentzVector p4M12_BX_rotCS_rotZ1(p4M12_BX_rotCS);
    TLorentzVector p4M21_BX_rotCS_rotZ1(p4M21_BX_rotCS);
    TLorentzVector p4M22_BX_rotCS_rotZ1(p4M22_BX_rotCS);
    p4M11_BX_rotCS_rotZ1.Transform(rotationZ1);
    p4M12_BX_rotCS_rotZ1.Transform(rotationZ1);
    p4M21_BX_rotCS_rotZ1.Transform(rotationZ1);
    p4M22_BX_rotCS_rotZ1.Transform(rotationZ1);
    p4M11_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
    p4M12_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
    p4M21_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);
    p4M22_BX_rotCS_rotZ1.Boost(boostZ1inX_rotCS_rotZ1);

    TLorentzVector p4V2_BX_rotCS_rotZ1 = p4M21_BX_rotCS_rotZ1 + p4M22_BX_rotCS_rotZ1;
    //// costheta1
    costheta1 = -p4V2_BX_rotCS_rotZ1.Vect().Dot(p4M11_BX_rotCS_rotZ1.Vect())/p4V2_BX_rotCS_rotZ1.Vect().Mag()/p4M11_BX_rotCS_rotZ1.Vect().Mag();
  }
  else costheta1=0;

  //// --------------------------- costheta2
  if (!(fabs(Z2_lept1Id)==21 || fabs(Z2_lept1Id)==22 || fabs(Z2_lept2Id)==21 || fabs(Z2_lept2Id)==22)){
    //define Z2 rotation 
    TRotation rotationZ2;
    TVector3 newZaxisZ2(thep4Z2inXFrame_rotCS.Vect().Unit());
    TVector3 newXaxisZ2(newYaxisCS.Cross(newZaxisZ2).Unit());
    TVector3 newYaxisZ2(newZaxisZ2.Cross(newXaxisZ2).Unit());
    rotationZ2.RotateAxes(newXaxisZ2, newYaxisZ2, newZaxisZ2);
    rotationZ2.Invert();

    TLorentzVector thep4Z2inXFrame_rotCS_rotZ2(thep4Z2inXFrame_rotCS);
    thep4Z2inXFrame_rotCS_rotZ2.Transform(rotationZ2);
    TVector3 boostZ2inX_rotCS_rotZ2= -(thep4Z2inXFrame_rotCS_rotZ2.BoostVector());

    TLorentzVector p4M11_BX_rotCS_rotZ2(p4M11_BX_rotCS);
    TLorentzVector p4M12_BX_rotCS_rotZ2(p4M12_BX_rotCS);
    TLorentzVector p4M21_BX_rotCS_rotZ2(p4M21_BX_rotCS);
    TLorentzVector p4M22_BX_rotCS_rotZ2(p4M22_BX_rotCS);
    p4M11_BX_rotCS_rotZ2.Transform(rotationZ2);
    p4M12_BX_rotCS_rotZ2.Transform(rotationZ2);
    p4M21_BX_rotCS_rotZ2.Transform(rotationZ2);
    p4M22_BX_rotCS_rotZ2.Transform(rotationZ2);
    p4M11_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
    p4M12_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
    p4M21_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);
    p4M22_BX_rotCS_rotZ2.Boost(boostZ2inX_rotCS_rotZ2);


    TLorentzVector p4V1_BX_rotCS_rotZ2= p4M11_BX_rotCS_rotZ2 + p4M12_BX_rotCS_rotZ2;
    //// costheta2
    costheta2 = -p4V1_BX_rotCS_rotZ2.Vect().Dot(p4M21_BX_rotCS_rotZ2.Vect())/p4V1_BX_rotCS_rotZ2.Vect().Mag()/p4M21_BX_rotCS_rotZ2.Vect().Mag();
  }
  else costheta2=0;

  if (isnan(costhetastar) || isnan(costheta1) || isnan(costheta2) || isnan(Phi) || isnan(Phi1)){
    MELAout << "WARNING: NaN in computeAngles: "
      << costhetastar << " "
      << costheta1  << " "
      << costheta2  << " "
      << Phi  << " "
      << Phi1  << " " << endl;
  }
}

/***** Associated production angles *****/
void TUtil::computeVBFAngles(
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  float& Q2V1,
  float& Q2V2,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  TLorentzVector jet1, int jet1Id,
  TLorentzVector jet2, int jet2Id,
  TLorentzVector* injet1, int injet1Id, // Gen. partons in lab frame
  TLorentzVector* injet2, int injet2Id
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (p4M12==nullFourVector || p4M22==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f13Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M21, Z2_lept1Id);
    p4M11 = f13Pair.first;
    p4M21 = f13Pair.second;
  }
  else if (p4M11==nullFourVector || p4M21==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f24Pair = TUtil::removeMassFromPair(p4M12, Z1_lept2Id, p4M22, Z2_lept2Id);
    p4M12 = f24Pair.first;
    p4M22 = f24Pair.second;
  }
  else{
    pair<TLorentzVector, TLorentzVector> f12Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M12, Z1_lept2Id);
    pair<TLorentzVector, TLorentzVector> f34Pair = TUtil::removeMassFromPair(p4M21, Z2_lept1Id, p4M22, Z2_lept2Id);
    p4M11 = f12Pair.first;
    p4M12 = f12Pair.second;
    p4M21 = f34Pair.first;
    p4M22 = f34Pair.second;
  }

  TLorentzVector jet1massless, jet2massless;
  pair<TLorentzVector, TLorentzVector> jetPair = TUtil::removeMassFromPair(jet1, jet1Id, jet2, jet2Id);
  jet1massless = jetPair.first;
  jet2massless = jetPair.second;

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;
  TLorentzVector pH = p4Z1+p4Z2;
  //jet1 is defined as going forwards (bigger pz), jet2 going backwards (smaller pz)
  if (jet1massless.Z() < jet2massless.Z()) { swap(jet1massless, jet2massless); swap(jet1Id, jet2Id); }

  //Find the incoming partons - first boost so the pT(HJJ) = 0, then boost away the pz.
  //This preserves the z direction.  Then assume that the partons come in this z direction.
  //This is exactly correct at LO (since pT=0 anyway).
  //Then associate the one going forwards with jet1 and the one going backwards with jet2
  TLorentzRotation movingframe;
  TLorentzVector pHJJ = pH+jet1massless+jet2massless;
  TLorentzVector pHJJ_perp(pHJJ.X(), pHJJ.Y(), 0, pHJJ.T());
  movingframe.Boost(-pHJJ_perp.BoostVector());
  pHJJ.Boost(-pHJJ_perp.BoostVector());
  movingframe.Boost(-pHJJ.BoostVector());
  pHJJ.Boost(-pHJJ.BoostVector());   //make sure to boost HJJ AFTER boosting movingframe

  TLorentzVector P1(0, 0, pHJJ.T()/2, pHJJ.T()/2);
  TLorentzVector P2(0, 0, -pHJJ.T()/2, pHJJ.T()/2);
  // Transform incoming partons back to original frame
  P1.Transform(movingframe.Inverse());
  P2.Transform(movingframe.Inverse());
  pHJJ.Transform(movingframe.Inverse());
  // movingframe and HJJ_T will not be used anymore
  // Handle gen. partons if they are available
  if (injet1 && injet2 && fabs((*injet1+*injet2).P()-pHJJ.P())<pHJJ.P()*1e-4){
    P1=*injet1;
    P2=*injet2;
    if (P1.Z() < P2.Z()){
      swap(P1, P2);
      swap(injet1Id, injet2Id);
    }
    // In the case of gen. partons, check if the intermediates are a Z or a W.
    int diff1Id = jet1Id-injet1Id;
    int diff2Id = jet2Id-injet2Id;
    if (
      !( // THIS IS A NOT-IF!
      (diff1Id==0 && diff2Id==0 && !(injet1Id==21 || injet2Id==21)) // Two Z bosons
        ||
        ((fabs(diff1Id)==1 || fabs(diff1Id)==3 || fabs(diff1Id)==5) && (fabs(diff2Id)==1 || fabs(diff2Id)==3 || fabs(diff2Id)==5)) // Two W bosons, do not check W+ vs W-
        )
      ){
      int diff12Id = jet1Id-injet2Id;
      int diff21Id = jet2Id-injet1Id;
      if (
        ((diff12Id==0 || diff21Id==0) && !(injet1Id==21 || injet2Id==21)) // At least one Z boson
        ||
        ((fabs(diff12Id)==1 || fabs(diff12Id)==3 || fabs(diff12Id)==5) || (fabs(diff21Id)==1 || fabs(diff21Id)==3 || fabs(diff21Id)==5)) // At least one W boson
        ){
        swap(P1, P2);
        swap(injet1Id, injet2Id);
      }
    }
  }

  TLorentzRotation ZZframe;
  ZZframe.Boost(-pH.BoostVector());
  P1.Transform(ZZframe);
  P2.Transform(ZZframe);
  p4Z1.Transform(ZZframe);
  p4Z2.Transform(ZZframe);
  jet1massless.Transform(ZZframe);
  jet2massless.Transform(ZZframe);

  TLorentzVector fermion1, fermion2, antifermion1, antifermion2;
  // Consider cases with gen. partons
  // By default, fermion1/2 are jet1/2 unless incoming partons are specifically anti-quarks!
  if (injet1!=0 && injet1Id<0){
    fermion1=-P1;
    antifermion1 = jet1massless;
  }
  else{
    fermion1 = jet1massless;
    antifermion1=-P1;
  }
  if (injet2!=0 && injet2Id<0){
    fermion2=-P2;
    antifermion2 = jet2massless;
  }
  else{
    fermion2 = jet2massless;
    antifermion2=-P2;
  }

  // Computations in the frame of X
  TLorentzVector V1 = fermion1 + antifermion1; // Outgoing V1
  TLorentzVector V2 = fermion2 + antifermion2; // Outgoing V2
  TVector3 normvec1 = fermion1.Vect().Cross(antifermion1.Vect()).Unit(); // p11 x p12
  TVector3 normvec2 = fermion2.Vect().Cross(antifermion2.Vect()).Unit(); // p21 x p22
  TVector3 normvec3 = p4Z2.Vect().Cross(V1.Vect()).Unit(); // z x V1
  double cosPhi = normvec1.Dot(normvec2);
  double sgnPhi = normvec1.Cross(normvec2).Dot(V1.Vect());
  if (fabs(sgnPhi)>0.) sgnPhi = sgnPhi/fabs(sgnPhi);
  double cosPhi1 = normvec1.Dot(normvec3);
  double sgnPhi1 = normvec1.Cross(normvec3).Dot(V1.Vect());
  if (fabs(sgnPhi1)>0.) sgnPhi1 = sgnPhi1/fabs(sgnPhi1);
  if (fabs(cosPhi)>1) cosPhi *= 1./fabs(cosPhi);
  if (fabs(cosPhi1)>1) cosPhi1 *= 1./fabs(cosPhi1);
  Phi = acos(-cosPhi)*sgnPhi;
  Phi1 = acos(cosPhi1)*sgnPhi1;
  costhetastar = V1.Vect().Unit().Dot(p4Z2.Vect().Unit()); // Note that p4Z2 is still in outgoing convention, so in incoming terms, this is equivalent to putting p4Z1.
  Q2V1 = -(V1.M2());
  Q2V2 = -(V2.M2());

  // Computations that would have been in the frame of X had V1 and V2 not been virtual
  costheta1 = V1.Vect().Unit().Dot(fermion1.Vect().Unit());
  costheta2 = V2.Vect().Unit().Dot(fermion2.Vect().Unit());
}
void TUtil::computeVBFAngles_ComplexBoost(
  float& costhetastar,
  float& costheta1_real, float& costheta1_imag,
  float& costheta2_real, float& costheta2_imag,
  float& Phi,
  float& Phi1,
  float& Q2V1,
  float& Q2V2,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  TLorentzVector jet1, int jet1Id,
  TLorentzVector jet2, int jet2Id,
  TLorentzVector* injet1, int injet1Id, // Gen. partons in lab frame
  TLorentzVector* injet2, int injet2Id
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (p4M12==nullFourVector || p4M22==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f13Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M21, Z2_lept1Id);
    p4M11 = f13Pair.first;
    p4M21 = f13Pair.second;
  }
  else if (p4M11==nullFourVector || p4M21==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f24Pair = TUtil::removeMassFromPair(p4M12, Z1_lept2Id, p4M22, Z2_lept2Id);
    p4M12 = f24Pair.first;
    p4M22 = f24Pair.second;
  }
  else{
    pair<TLorentzVector, TLorentzVector> f12Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M12, Z1_lept2Id);
    pair<TLorentzVector, TLorentzVector> f34Pair = TUtil::removeMassFromPair(p4M21, Z2_lept1Id, p4M22, Z2_lept2Id);
    p4M11 = f12Pair.first;
    p4M12 = f12Pair.second;
    p4M21 = f34Pair.first;
    p4M22 = f34Pair.second;
  }

  TLorentzVector jet1massless, jet2massless;
  pair<TLorentzVector, TLorentzVector> jetPair = TUtil::removeMassFromPair(jet1, jet1Id, jet2, jet2Id);
  jet1massless = jetPair.first;
  jet2massless = jetPair.second;

  //build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;
  TLorentzVector pH = p4Z1+p4Z2;
  //jet1 is defined as going forwards (bigger pz), jet2 going backwards (smaller pz)
  if (jet1massless.Z() < jet2massless.Z()) { swap(jet1massless, jet2massless); swap(jet1Id, jet2Id); }

  //Find the incoming partons - first boost so the pT(HJJ) = 0, then boost away the pz.
  //This preserves the z direction.  Then assume that the partons come in this z direction.
  //This is exactly correct at LO (since pT=0 anyway).
  //Then associate the one going forwards with jet1 and the one going backwards with jet2
  TLorentzRotation movingframe;
  TLorentzVector pHJJ = pH+jet1massless+jet2massless;
  TLorentzVector pHJJ_perp(pHJJ.X(), pHJJ.Y(), 0, pHJJ.T());
  movingframe.Boost(-pHJJ_perp.BoostVector());
  pHJJ.Boost(-pHJJ_perp.BoostVector());
  movingframe.Boost(-pHJJ.BoostVector());
  pHJJ.Boost(-pHJJ.BoostVector());   //make sure to boost HJJ AFTER boosting movingframe

  TLorentzVector P1(0, 0, pHJJ.T()/2, pHJJ.T()/2);
  TLorentzVector P2(0, 0, -pHJJ.T()/2, pHJJ.T()/2);
  // Transform incoming partons back to original frame
  P1.Transform(movingframe.Inverse());
  P2.Transform(movingframe.Inverse());
  pHJJ.Transform(movingframe.Inverse());
  // movingframe and HJJ_T will not be used anymore
  // Handle gen. partons if they are available
  if (injet1 && injet2 && fabs((*injet1+*injet2).P()-pHJJ.P())<pHJJ.P()*1e-4){
    P1=*injet1;
    P2=*injet2;
    if (P1.Z() < P2.Z()){
      swap(P1, P2);
      swap(injet1Id, injet2Id);
    }
    // In the case of gen. partons, check if the intermediates are a Z or a W.
    int diff1Id = jet1Id-injet1Id;
    int diff2Id = jet2Id-injet2Id;
    if (
      !( // THIS IS A NOT-IF!
      (diff1Id==0 && diff2Id==0 && !(injet1Id==21 || injet2Id==21)) // Two Z bosons
        ||
        ((fabs(diff1Id)==1 || fabs(diff1Id)==3 || fabs(diff1Id)==5) && (fabs(diff2Id)==1 || fabs(diff2Id)==3 || fabs(diff2Id)==5)) // Two W bosons, do not check W+ vs W-
        )
      ){
      int diff12Id = jet1Id-injet2Id;
      int diff21Id = jet2Id-injet1Id;
      if (
        ((diff12Id==0 || diff21Id==0) && !(injet1Id==21 || injet2Id==21)) // At least one Z boson
        ||
        ((fabs(diff12Id)==1 || fabs(diff12Id)==3 || fabs(diff12Id)==5) || (fabs(diff21Id)==1 || fabs(diff21Id)==3 || fabs(diff21Id)==5)) // At least one W boson
        ){
        swap(P1, P2);
        swap(injet1Id, injet2Id);
      }
    }
  }

  TLorentzRotation ZZframe;
  ZZframe.Boost(-pH.BoostVector());
  P1.Transform(ZZframe);
  P2.Transform(ZZframe);
  p4Z1.Transform(ZZframe);
  p4Z2.Transform(ZZframe);
  jet1massless.Transform(ZZframe);
  jet2massless.Transform(ZZframe);

  TLorentzVector fermion1, fermion2, antifermion1, antifermion2;
  // Consider cases with gen. partons
  // By default, fermion1/2 are jet1/2 unless incoming partons are specifically anti-quarks!
  if (injet1!=0 && injet1Id<0){
    fermion1=-P1;
    antifermion1 = jet1massless;
  }
  else{
    fermion1 = jet1massless;
    antifermion1=-P1;
  }
  if (injet2!=0 && injet2Id<0){
    fermion2=-P2;
    antifermion2 = jet2massless;
  }
  else{
    fermion2 = jet2massless;
    antifermion2=-P2;
  }

  // Computations in the frame of X
  TLorentzVector V1 = fermion1 + antifermion1; // Outgoing V1
  TLorentzVector V2 = fermion2 + antifermion2; // Outgoing V2
  TVector3 normvec1 = fermion1.Vect().Cross(antifermion1.Vect()).Unit(); // p11 x p12
  TVector3 normvec2 = fermion2.Vect().Cross(antifermion2.Vect()).Unit(); // p21 x p22
  TVector3 normvec3 = p4Z2.Vect().Cross(V1.Vect()).Unit(); // z x V1
  double cosPhi = normvec1.Dot(normvec2);
  double sgnPhi = normvec1.Cross(normvec2).Dot(V1.Vect());
  if (fabs(sgnPhi)>0.) sgnPhi = sgnPhi/fabs(sgnPhi);
  double cosPhi1 = normvec1.Dot(normvec3);
  double sgnPhi1 = normvec1.Cross(normvec3).Dot(V1.Vect());
  if (fabs(sgnPhi1)>0.) sgnPhi1 = sgnPhi1/fabs(sgnPhi1);
  if (fabs(cosPhi)>1) cosPhi *= 1./fabs(cosPhi);
  if (fabs(cosPhi1)>1) cosPhi1 *= 1./fabs(cosPhi1);
  Phi = acos(-cosPhi)*sgnPhi;
  Phi1 = acos(cosPhi1)*sgnPhi1;
  costhetastar = V1.Vect().Unit().Dot(p4Z2.Vect().Unit()); // Note that p4Z2 is still in outgoing convention, so in incoming terms, this is equivalent to putting p4Z1.
  Q2V1 = -(V1.M2());
  Q2V2 = -(V2.M2());

  // Up to here, everything has to be the same as TUtil::computeVBFAngles
  // Computations that would have been truly in the frame of X had V1 and V2 not been virtual:
  // Use TUtil::ComplexBoost to evade imaginary gamma problems when beta**2<0
  pair<TLorentzVector, TLorentzVector> V2_BV1 = TUtil::ComplexBoost(V1.BoostVector(), V2);
  pair<TLorentzVector, TLorentzVector> fermion1_BV1 = TUtil::ComplexBoost(V1.BoostVector(), fermion1);
  costheta1_real = -(V2_BV1.first.Vect().Unit().Dot(fermion1_BV1.first.Vect().Unit()) - V2_BV1.second.Vect().Unit().Dot(fermion1_BV1.second.Vect().Unit()));
  costheta1_imag = -(V2_BV1.first.Vect().Unit().Dot(fermion1_BV1.second.Vect().Unit()) + V2_BV1.second.Vect().Unit().Dot(fermion1_BV1.first.Vect().Unit()));

  pair<TLorentzVector, TLorentzVector> V1_BV2 = TUtil::ComplexBoost(V2.BoostVector(), V1);
  pair<TLorentzVector, TLorentzVector> fermion2_BV2 = TUtil::ComplexBoost(V2.BoostVector(), fermion2);
  costheta2_real = -(V1_BV2.first.Vect().Unit().Dot(fermion2_BV2.first.Vect().Unit()) - V1_BV2.second.Vect().Unit().Dot(fermion2_BV2.second.Vect().Unit()));
  costheta2_imag = -(V1_BV2.first.Vect().Unit().Dot(fermion2_BV2.second.Vect().Unit()) + V1_BV2.second.Vect().Unit().Dot(fermion2_BV2.first.Vect().Unit()));
}
void TUtil::computeVHAngles(
  float& costhetastar,
  float& costheta1,
  float& costheta2,
  float& Phi,
  float& Phi1,
  TLorentzVector p4M11, int Z1_lept1Id,
  TLorentzVector p4M12, int Z1_lept2Id,
  TLorentzVector p4M21, int Z2_lept1Id,
  TLorentzVector p4M22, int Z2_lept2Id,
  TLorentzVector jet1, int jet1Id,
  TLorentzVector jet2, int jet2Id,
  TLorentzVector* injet1, int injet1Id, // Gen. partons in lab frame
  TLorentzVector* injet2, int injet2Id
  ){
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (p4M12==nullFourVector || p4M22==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f13Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M21, Z2_lept1Id);
    p4M11 = f13Pair.first;
    p4M21 = f13Pair.second;
  }
  else if (p4M11==nullFourVector || p4M21==nullFourVector){
    pair<TLorentzVector, TLorentzVector> f24Pair = TUtil::removeMassFromPair(p4M12, Z1_lept2Id, p4M22, Z2_lept2Id);
    p4M12 = f24Pair.first;
    p4M22 = f24Pair.second;
  }
  else{
    pair<TLorentzVector, TLorentzVector> f12Pair = TUtil::removeMassFromPair(p4M11, Z1_lept1Id, p4M12, Z1_lept2Id);
    pair<TLorentzVector, TLorentzVector> f34Pair = TUtil::removeMassFromPair(p4M21, Z2_lept1Id, p4M22, Z2_lept2Id);
    p4M11 = f12Pair.first;
    p4M12 = f12Pair.second;
    p4M21 = f34Pair.first;
    p4M22 = f34Pair.second;
  }

  TLorentzVector jet1massless, jet2massless;
  pair<TLorentzVector, TLorentzVector> jetPair = TUtil::removeMassFromPair(jet1, jet1Id, jet2, jet2Id);
  jet1massless = jetPair.first;
  jet2massless = jetPair.second;

  // Build Z 4-vectors
  TLorentzVector p4Z1 = p4M11 + p4M12;
  TLorentzVector p4Z2 = p4M21 + p4M22;
  TLorentzVector pH = p4Z1 + p4Z2;

  // Apply convention for outgoing particles
  if (
    (jet1Id*jet2Id<0 && jet1Id<0) // for OS pairs: jet1 must be the particle
    ||
    (jet1Id*jet2Id>0 && jet1massless.Phi()<=jet2massless.Phi()) // for SS pairs: use random deterministic convention
    ){
    swap(jet1massless, jet2massless);
    swap(jet1Id, jet2Id);
  }

  //Find the incoming partons - first boost so the pT(HJJ) = 0, then boost away the pz.
  //This preserves the z direction.  Then assume that the partons come in this z direction.
  //This is exactly correct at LO (since pT=0 anyway).
  //Then associate the one going forwards with jet1 and the one going backwards with jet2
  TLorentzRotation movingframe;
  TLorentzVector pHJJ = pH+jet1massless+jet2massless;
  TLorentzVector pHJJ_perp(pHJJ.X(), pHJJ.Y(), 0, pHJJ.T());
  movingframe.Boost(-pHJJ_perp.BoostVector());
  pHJJ.Boost(-pHJJ_perp.BoostVector());
  movingframe.Boost(-pHJJ.BoostVector());
  pHJJ.Boost(-pHJJ.BoostVector());   //make sure to boost HJJ AFTER boosting movingframe

  TLorentzVector P1(0, 0, -pHJJ.T()/2, pHJJ.T()/2);
  TLorentzVector P2(0, 0, pHJJ.T()/2, pHJJ.T()/2);
  // Transform incoming partons back to the original frame
  P1.Transform(movingframe.Inverse());
  P2.Transform(movingframe.Inverse());
  pHJJ.Transform(movingframe.Inverse());
  // movingframe and HJJ_T will not be used anymore
  // Handle gen. partons if they are available
  if (injet1 && injet2 && fabs((*injet1+*injet2).P()-pHJJ.P())<pHJJ.P()*1e-4){
    P1=*injet1;
    P2=*injet2;
    // Apply convention for incoming (!) particles
    if (
      (injet1Id*injet2Id<0 && injet1Id>0) // for OS pairs: parton 2 must be the particle
      ||
      (injet1Id*injet2Id>0 && P1.Z()>=P2.Z()) //for SS pairs: use random deterministic convention
      ){
      swap(P1, P2);
      swap(injet1Id, injet2Id);
    }
  }

  // Rotate every vector such that Z1 - Z2 axis is the "beam axis" analogue of decay
  TLorentzRotation ZZframe;
  TVector3 beamAxis(0, 0, 1);
  if (p4Z1==nullFourVector || p4Z2==nullFourVector){
    TVector3 pNewAxis = (p4Z2-p4Z1).Vect().Unit(); // Let Z2 be in the z direction so that once the direction of H is reversed, Z1 is in the z direction
    if (pNewAxis != nullFourVector.Vect()){
      TVector3 pNewAxisPerp = pNewAxis.Cross(beamAxis);
      ZZframe.Rotate(acos(pNewAxis.Dot(beamAxis)), pNewAxisPerp);

      P1.Transform(ZZframe);
      P2.Transform(ZZframe);
      jet1massless = -jet1massless;
      jet2massless = -jet2massless;
      jet1massless.Transform(ZZframe);
      jet2massless.Transform(ZZframe);
      jet1massless = -jet1massless;
      jet2massless = -jet2massless;
    }
  }
  else{
    ZZframe.Boost(-pH.BoostVector());
    p4Z1.Boost(-pH.BoostVector());
    p4Z2.Boost(-pH.BoostVector());
    TVector3 pNewAxis = (p4Z2-p4Z1).Vect().Unit(); // Let Z2 be in the z direction so that once the direction of H is reversed, Z1 is in the z direction
    TVector3 pNewAxisPerp = pNewAxis.Cross(beamAxis);
    ZZframe.Rotate(acos(pNewAxis.Dot(beamAxis)), pNewAxisPerp);
    P1.Transform(ZZframe);
    P2.Transform(ZZframe);
    jet1massless = -jet1massless;
    jet2massless = -jet2massless;
    jet1massless.Transform(ZZframe);
    jet2massless.Transform(ZZframe);
    jet1massless = -jet1massless;
    jet2massless = -jet2massless;
  }

  TUtil::computeAngles(
    costhetastar,
    costheta1,
    costheta2,
    Phi,
    Phi1,
    -P1, 23, // Id is 23 to avoid an attempt to remove quark mass
    -P2, 0, // Id is 0 to avoid swapping
    jet1massless, 23,
    jet2massless, 0
  );
}


/****************************************************/
/***** JHUGen- and MCFM-related ME computations *****/
/****************************************************/

void TUtil::SetEwkCouplingParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme){
  // Set JHUGen couplings
  const double GeV=1./100.;
  double ext_mZ_jhu = ext_mZ*GeV;
  double ext_mW_jhu = ext_mW*GeV;
  double ext_Gf_jhu = ext_Gf/pow(GeV, 2);
  __modjhugenmela_MOD_setewparameters(&ext_mZ_jhu, &ext_mW_jhu, &ext_Gf_jhu, &ext_aemmz, &ext_xW);

  // Set MCFM couplings
  if (ext_ewscheme<-1 || ext_ewscheme>3) ext_ewscheme=3;
  ewinput_.Gf_inp = ext_Gf;
  ewinput_.aemmz_inp = ext_aemmz;
  ewinput_.wmass_inp = ext_mW;
  ewinput_.zmass_inp = ext_mZ;
  ewinput_.xw_inp = ext_xW;
  ewscheme_.ewscheme = ext_ewscheme;
  coupling_();

  // SETTINGS TO MATCH JHUGen ME/generator:
  /*
  ewscheme_.ewscheme = 3; // Switch ewscheme to full control, default is 1
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.aemmz_inp=1./128.;
  ewinput_.wmass_inp=80.399;
  ewinput_.zmass_inp=91.1876;
  ewinput_.xw_inp=0.23119;
  */

  // Settings used in Run I MC, un-synchronized to JHUGen and Run II generation
  /*
  ewinput_.Gf_inp = 1.16639E-05;
  ewinput_.wmass_inp = 80.385;
  ewinput_.zmass_inp = 91.1876;
  //ewinput_.aemmz_inp = 7.81751e-3; // Not used in the compiled new default MCFM ewcheme=1
  //ewinput_.xw_inp = 0.23116864; // Not used in the compiled new default MCFM ewcheme=1
  ewinput_.aemmz_inp = 7.562468901984759e-3;
  ewinput_.xw_inp = 0.2228972225239183;
  */

  // INPUT SETTINGS in default MCFM generator:
  /*
  ewscheme_.ewscheme = 1;
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.wmass_inp=80.398;
  ewinput_.zmass_inp=91.1876;
  ewinput_.aemmz_inp=7.7585538055706e-3;
  ewinput_.xw_inp=0.2312;
  */

  // ACTUAL SETTINGS in default MCFM generator, gives the same values as above but is more explicit
  /*
  ewscheme_.ewscheme = 3;
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.wmass_inp=80.398;
  ewinput_.zmass_inp=91.1876;
  ewinput_.aemmz_inp=7.55638390728736e-3;
  ewinput_.xw_inp=0.22264585341299625;
  */
}
void TUtil::SetMass(double inmass, int ipart){
  const int ipartabs = abs(ipart);
  bool runcoupling_mcfm=false;
  bool runcoupling_jhugen=false;

  // MCFM masses
  // Tprime and bprime masses are not defined in masses.f
  if (ipartabs==8) spinzerohiggs_anomcoupl_.mt_4gen = inmass;
  else if (ipartabs==7) spinzerohiggs_anomcoupl_.mb_4gen = inmass;
  else if (ipartabs==6){ masses_mcfm_.mt=inmass; runcoupling_mcfm=true; }
  else if (ipartabs==5){ masses_mcfm_.mb=inmass; masses_mcfm_.mbsq = pow(masses_mcfm_.mb, 2); runcoupling_mcfm=true; }
  else if (ipartabs==4){ masses_mcfm_.mc=inmass; masses_mcfm_.mcsq = pow(masses_mcfm_.mc, 2); runcoupling_mcfm=true; }
  else if (ipartabs==3) masses_mcfm_.ms=inmass;
  else if (ipartabs==2) masses_mcfm_.mu=inmass;
  else if (ipartabs==1) masses_mcfm_.md=inmass;
  else if (ipartabs==11) masses_mcfm_.mel=inmass;
  else if (ipartabs==13) masses_mcfm_.mmu=inmass;
  else if (ipartabs==15){ masses_mcfm_.mtau=inmass; masses_mcfm_.mtausq = pow(masses_mcfm_.mtau, 2); }
  else if (ipartabs==23){ masses_mcfm_.zmass=inmass; ewinput_.zmass_inp = inmass; runcoupling_mcfm=true; }
  else if (ipartabs==24){ masses_mcfm_.wmass=inmass; ewinput_.wmass_inp = inmass; runcoupling_mcfm=true; }
  else if (ipartabs==25) masses_mcfm_.hmass=inmass;

  // JHUGen masses
  if (
    ipartabs<=6
    ||
    (ipartabs>=11 && ipartabs<=16)
    ||
    ipartabs==23 || ipartabs==24 || ipartabs==25
    ||
    ipartabs==32 || ipartabs==34
    ){
    runcoupling_jhugen=(
      ipartabs==23
      ||
      ipartabs==24
      );
    const double GeV=1./100.;
    double jinmass = inmass*GeV;
    int jpart = convertLHEreverse(&ipart);
    __modparameters_MOD_setmass(&jinmass, &jpart);
  }

  // Recalculate couplings
  if (runcoupling_mcfm || runcoupling_jhugen) SetEwkCouplingParameters(ewcouple_.Gf, em_.aemmz, masses_mcfm_.wmass, masses_mcfm_.zmass, ewcouple_.xw, ewscheme_.ewscheme);
}
void TUtil::SetDecayWidth(double inwidth, int ipart){
  const int ipartabs = abs(ipart);
  // No need to recalculate couplings
  // MCFM masses
  if (ipartabs==6) masses_mcfm_.twidth=inwidth;
  else if (ipartabs==15) masses_mcfm_.tauwidth=inwidth;
  else if (ipartabs==23) masses_mcfm_.zwidth=inwidth;
  else if (ipartabs==24) masses_mcfm_.wwidth=inwidth;
  else if (ipartabs==25) masses_mcfm_.hwidth=inwidth;

  // JHUGen masses
  const double GeV=1./100.;
  double jinwidth = inwidth*GeV;
  int jpart = convertLHEreverse(&ipart);
  __modparameters_MOD_setdecaywidth(&jinwidth, &jpart);
}
void TUtil::SetCKMElements(double* invckm_ud, double* invckm_us, double* invckm_cd, double* invckm_cs, double* invckm_ts, double* invckm_tb, double* invckm_ub, double* invckm_cb, double* invckm_td){
  __modparameters_MOD_computeckmelements(invckm_ud, invckm_us, invckm_cd, invckm_cs, invckm_ts, invckm_tb, invckm_ub, invckm_cb, invckm_td);

  int i, j;
  i=2; j=1;
  cabib_.Vud = __modparameters_MOD_ckmbare(&i, &j);
  i=2; j=3;
  cabib_.Vus = __modparameters_MOD_ckmbare(&i, &j);
  i=2; j=5;
  cabib_.Vub = __modparameters_MOD_ckmbare(&i, &j);
  i=4; j=1;
  cabib_.Vcd = __modparameters_MOD_ckmbare(&i, &j);
  i=4; j=3;
  cabib_.Vcs = __modparameters_MOD_ckmbare(&i, &j);
  i=4; j=5;
  cabib_.Vcb = __modparameters_MOD_ckmbare(&i, &j);
  // Do not call ckmfill_(), it is called by MCFM_chooser!
}
double TUtil::GetCKMElement(int iquark, int jquark){
  return __modparameters_MOD_ckmbare(&iquark, &jquark);
}

double TUtil::GetMass(int ipart){
  const int ipartabs = abs(ipart);

  if (ipartabs==8) return spinzerohiggs_anomcoupl_.mt_4gen;
  else if (ipartabs==7) return spinzerohiggs_anomcoupl_.mb_4gen;
  else if (ipartabs==6) return masses_mcfm_.mt;
  else if (ipartabs==5) return masses_mcfm_.mb;
  else if (ipartabs==4) return masses_mcfm_.mc;
  else if (ipartabs==3) return masses_mcfm_.ms;
  else if (ipartabs==2) return masses_mcfm_.mu;
  else if (ipartabs==1) return masses_mcfm_.md;
  else if (ipartabs==11) return masses_mcfm_.mel;
  else if (ipartabs==13) return masses_mcfm_.mmu;
  else if (ipartabs==15) return masses_mcfm_.mtau;
  else if (ipartabs==23) return masses_mcfm_.zmass;
  else if (ipartabs==24) return masses_mcfm_.wmass;
  else if (ipartabs==25) return masses_mcfm_.hmass;
  else{
    // JHUGen masses
    if (
      ipartabs<=6 // (d, u, s, c, b, t)
      ||
      (ipartabs>=11 && ipartabs<=16) // (l, nu) x (e, mu, tau)
      ||
      ipartabs==23 || ipartabs==24 || ipartabs==25 // Z, W, H
      ||
      ipartabs==32 || ipartabs==34 // Z', W+-'
      ){
      const double GeV=1./100.;
      int jpart = convertLHEreverse(&ipart);
      double joutmass = __modparameters_MOD_getmass(&jpart);
      double outmass = joutmass/GeV;
      return outmass;
    }
    else return 0;
  }
}
double TUtil::GetDecayWidth(int ipart){
  const int ipartabs = abs(ipart);

  // MCFM widths
  if (ipartabs==6) return masses_mcfm_.twidth;
  else if (ipartabs==15) return masses_mcfm_.tauwidth;
  else if (ipartabs==23) return masses_mcfm_.zwidth;
  else if (ipartabs==24) return masses_mcfm_.wwidth;
  else if (ipartabs==25) return masses_mcfm_.hwidth;
  else{
    // JHUGen widths
    const double GeV=1./100.;
    int jpart = convertLHEreverse(&ipart);
    double joutwidth = __modparameters_MOD_getdecaywidth(&jpart);
    double outwidth = joutwidth/GeV;
    return outwidth;
  }
}
double TUtil::GetMass(const MELAParticle* part){
  if (part==nullptr) return GetMass(-9000);
  return GetMass(part->id);
}
double TUtil::GetDecayWidth(const MELAParticle* part){
  if (part==nullptr) return GetDecayWidth(-9000);
  return GetDecayWidth(part->id);
}
void TUtil::GetMassWidth(int ipart, double& m, double& ga){
  m=TUtil::GetMass(ipart);
  ga=TUtil::GetDecayWidth(ipart);
}
void TUtil::GetMassWidth(const MELAParticle* part, double& m, double& ga){
  if (part==nullptr) return GetMassWidth(-9000, m, ga);
  return GetMassWidth(part->id, m, ga);
}

double TUtil::InterpretScaleScheme(const TVar::Production& production, const TVar::MatrixElement& matrixElement, const TVar::EventScaleScheme& scheme, TLorentzVector p[mxpart]){
  double Q=0;
  TLorentzVector nullFourVector(0, 0, 0, 0);
  if (scheme == TVar::Fixed_mH) Q = masses_mcfm_.hmass;
  else if (scheme == TVar::Fixed_mW) Q = masses_mcfm_.wmass;
  else if (scheme == TVar::Fixed_mZ) Q = masses_mcfm_.zmass;
  else if (scheme == TVar::Fixed_mWPlusmH) Q = (masses_mcfm_.wmass+masses_mcfm_.hmass);
  else if (scheme == TVar::Fixed_mZPlusmH) Q = (masses_mcfm_.zmass+masses_mcfm_.hmass);
  else if (scheme == TVar::Fixed_TwomtPlusmH) Q = (2.*masses_mcfm_.mt+masses_mcfm_.hmass);
  else if (scheme == TVar::Fixed_mtPlusmH) Q = masses_mcfm_.mt+masses_mcfm_.hmass;
  else if (scheme == TVar::Dynamic_qH){
    TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5];
    Q = fabs(pTotal.M());
  }
  else if (scheme == TVar::Dynamic_qJJH){
    TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5]+p[6]+p[7];
    Q = fabs(pTotal.M());
  }
  else if (scheme == TVar::Dynamic_qJJ_qH){
    TLorentzVector qH = p[2]+p[3]+p[4]+p[5];
    TLorentzVector qJJ = p[6]+p[7];
    Q = fabs(qH.M()+qJJ.M());
  }
  else if (scheme == TVar::Dynamic_qJ_qJ_qH){
    TLorentzVector qH = p[2]+p[3]+p[4]+p[5];
    Q = fabs(qH.M()+p[6].M()+p[7].M());
  }
  else if (scheme == TVar::Dynamic_HT){
    for (int c=2; c<mxpart; c++) Q += p[c].Pt(); // Scalar sum of all pTs
  }
  else if (scheme == TVar::DefaultScaleScheme){
    // Defaults are dynamic scales except in ttH and bbH.
    if (matrixElement==TVar::JHUGen){
      if (
        production==TVar::JJQCD
        || production==TVar::JJVBF
        || production==TVar::JQCD
        || production==TVar::ZZGG
        || production==TVar::ZZQQB
        || production==TVar::ZZQQB_STU
        || production==TVar::ZZQQB_S
        || production==TVar::ZZQQB_TU
        || production==TVar::ZZINDEPENDENT
        || production==TVar::Lep_WH || production==TVar::Had_WH
        || production==TVar::Lep_ZH || production==TVar::Had_ZH
        || production==TVar::GammaH
        ){
        TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5];
        Q = fabs(pTotal.M());
      }
      else if (
        production==TVar::ttH
        ||
        production==TVar::bbH
        ) Q = (2.*masses_mcfm_.mt+masses_mcfm_.hmass);
    }
    else if (matrixElement==TVar::MCFM){
      if (
        production==TVar::ZZGG
        || production==TVar::ZZINDEPENDENT
        || production==TVar::ZZQQB
        || production==TVar::ZZQQB_STU

        || production==TVar::JQCD

        || production==TVar::JJQCD
        || production==TVar::JJVBF
        || production==TVar::JJEW
        || production==TVar::JJEWQCD
        || production==TVar::Lep_WH || production==TVar::Had_WH
        || production==TVar::Lep_ZH || production==TVar::Had_ZH

        || production==TVar::ZZQQB_S
        || production==TVar::JJQCD_S
        || production==TVar::JJVBF_S
        || production==TVar::JJEW_S
        || production==TVar::JJEWQCD_S
        || production==TVar::Lep_WH_S || production==TVar::Had_WH_S
        || production==TVar::Lep_ZH_S || production==TVar::Had_ZH_S

        || production==TVar::ZZQQB_TU
        || production==TVar::JJQCD_TU
        || production==TVar::JJVBF_TU
        || production==TVar::JJEW_TU
        || production==TVar::JJEWQCD_TU
        || production==TVar::Lep_WH_TU || production==TVar::Had_WH_TU
        || production==TVar::Lep_ZH_TU || production==TVar::Had_ZH_TU

        || production==TVar::GammaH
        ){
        TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5];
        Q = fabs(pTotal.M());
      }
      else if (
        production==TVar::ttH
        || production==TVar::bbH
        ){ // ttH and bbH are not implemented through MCFM, will need revision if they are implemented.
        TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5]+p[6]+p[7];
        Q = fabs(pTotal.M());
      }
    }
  }

  if (Q<=0.){
    MELAerr << "Scaling fails for production " << production << ", defaulting to dynamic scheme m3456 " << endl;
    TLorentzVector pTotal = p[2]+p[3]+p[4]+p[5];
    Q = fabs(pTotal.M());
  }
  return Q;
}
void TUtil::SetAlphaS(double& Q_ren, double& Q_fac, double multiplier_ren, double multiplier_fac, int mynloop, int mynflav, string mypartons){
  bool hasReset=false;
  if (multiplier_ren<=0. || multiplier_fac<=0.){
    MELAerr << "TUtil::SetAlphaS: Invalid scale multipliers" << endl;
    return;
  }
  if (Q_ren<=1. || Q_fac<=1. || mynloop<=0 || mypartons.compare("Default")==0){
    if (Q_ren<0.) MELAout << "TUtil::SetAlphaS: Invalid QCD scale for alpha_s, setting to mH/2..." << endl;
    if (Q_fac<0.) MELAout << "TUtil::SetAlphaS: Invalid factorization scale, setting to mH/2..." << endl;
    Q_ren = (masses_mcfm_.hmass)*0.5;
    Q_fac = Q_ren;
    mynloop = 1;
    hasReset=true;
  }
  if (!hasReset){
    Q_ren *= multiplier_ren;
    Q_fac *= multiplier_fac;
  }

  /***** MCFM Alpha_S *****/
  bool nflav_is_same = (nflav_.nflav == mynflav);
  if (!nflav_is_same) MELAout << "TUtil::SetAlphaS: nflav=" << nflav_.nflav << " is the only one supported." << endl;
  scale_.scale = Q_ren;
  scale_.musq = Q_ren*Q_ren;
  facscale_.facscale = Q_fac;
  if (mynloop!=1){
    MELAout << "TUtil::SetAlphaS: Only nloop=1 is supported!" << endl;
    mynloop=1;
  }
  nlooprun_.nlooprun = mynloop;

  /***** JHUGen Alpha_S *****/
  /*
  ///// Disabling alpha_s computation from MCFM to replace with the JHUGen implementation, allows LHAPDF interface readily /////
  // For proper pdfwrapper_linux.f execution (alpha_s computation does not use pdf but examines the pdf name to initialize amz.)
  if (mypartons.compare("Default")!=0 && mypartons.compare("cteq6_l")!=0 && mypartons.compare("cteq6l1")!=0){
  MELAout << "Only default=cteq6l1 or cteq6_l are supported. Modify mela.cc symlinks, put the pdf table into data/Pdfdata and retry. Setting mypartons to Default..." << endl;
  mypartons = "Default";
  }
  // From pdfwrapper_linux.f:
  if (mypartons.compare("cteq6_l")==0) couple_.amz = 0.118;
  else if (mypartons.compare("cteq6l1")==0 || mypartons.compare("Default")==0) couple_.amz = 0.130;
  else couple_.amz = 0.118; // Add pdf as appropriate

  if (!nflav_is_same){
    nflav_.nflav = mynflav;

    if (mypartons.compare("Default")!=0) sprintf(pdlabel_.pdlabel, "%s", mypartons.c_str());
    else sprintf(pdlabel_.pdlabel, "%s", "cteq6l1"); // Default pdf is cteq6l1
    coupling2_();
  }
  else qcdcouple_.as = alphas_(&(scale_.scale), &(couple_.amz), &(nlooprun_.nlooprun));
  */

  const double GeV=1./100.;
  double muren_jhu = scale_.scale*GeV;
  double mufac_jhu = facscale_.facscale*GeV;
  __modjhugenmela_MOD_setmurenfac(&muren_jhu, &mufac_jhu);
  __modkinematics_MOD_evalalphas();
  TUtil::GetAlphaS(&(qcdcouple_.as), &(couple_.amz));

  qcdcouple_.gsq = 4.0*TMath::Pi()*qcdcouple_.as;
  qcdcouple_.ason2pi = qcdcouple_.as/(2.0*TMath::Pi());
  qcdcouple_.ason4pi = qcdcouple_.as/(4.0*TMath::Pi());

  // TEST RUNNING SCALE PER EVENT:
  /*
  if(verbosity >= TVar::DEBUG){
  MELAout << "My pdf is: " << pdlabel_.pdlabel << endl;
  MELAout << "My Q_ren: " << Q_ren << " | My alpha_s: " << qcdcouple_.as << " at order " << nlooprun_.nlooprun << " with a(m_Z): " << couple_.amz << '\t'
  << "Nflav: " << nflav_.nflav << endl;
  */
}
void TUtil::GetAlphaS(double* alphas_, double* alphasmz_){
  double alphasVal=0, alphasmzVal=0;
  __modjhugenmela_MOD_getalphasalphasmz(&alphasVal, &alphasmzVal);
  if (alphas_!=0) *alphas_ = alphasVal;
  if (alphasmz_!=0) *alphasmz_ = alphasmzVal;
}

// chooser.f split into 2 different functions
bool TUtil::MCFM_chooser(
  const TVar::Process& process, const TVar::Production& production, const TVar::LeptonInterference& leptonInterf,
  const TVar::VerbosityLevel& verbosity,
  const simple_event_record& mela_event
  ){
  bool result = true;

  unsigned int ndau = mela_event.pDaughters.size();
  int* pId = new int[ndau];
  for (unsigned int ip=0; ip<ndau; ip++){
    pId[ip]=mela_event.pDaughters.at(ip).first;
  }
  bool isWW = (ndau>=4 && PDGHelpers::isAWBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAWBoson(mela_event.intermediateVid.at(1)));
  bool isZZ = (ndau>=4 && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(1)));
  bool hasZZ4fInterf = isZZ && abs(pId[0])==abs(pId[2]) && abs(pId[1])==abs(pId[3]) && !PDGHelpers::isAnUnknownJet(pId[0]) && !PDGHelpers::isAnUnknownJet(pId[3]);
  //bool isZJJ = false;
  //if (ndau>=4) isZJJ = (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAJet(pId[2]) && PDGHelpers::isAJet(pId[3])); // Notice both isZZ and isZJJ could be true
  //bool isZG = (ndau>=3 && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(1)));
  //bool isGG = (ndau>=2 && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(0)) && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(1)));
  if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: isWW=" << (int)isWW << ", isZZ=" << (int)isZZ << ", hasZZ4fInterf=" << (int)hasZZ4fInterf << endl;

  npart_.npart=4; // Default number of particles is just 4, indicating no associated particles
  sprintf(runstring_.runstring, "test");

  if (
    ndau>=4
    &&
    (
    ((production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU) && process == TVar::bkgZZ)
    ||
    ((production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB) && process == TVar::bkgZZ)
    )
    ){
    //81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
    //86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'
    //90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))' 'L'

    nqcdjets_.nqcdjets=0;
    srdiags_.srdiags=true;
    nwz_.nwz=0; ckmfill_(&(nwz_.nwz));
    bveg1_mcfm_.ndim=10;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    vsymfact_.vsymfact=1.0;
    interference_.interference=false;
    if (hasZZ4fInterf && (leptonInterf==TVar::DefaultLeptonInterf || leptonInterf==TVar::InterfOn)){
      //90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))' 'L'
      vsymfact_.vsymfact=0.5;
      interference_.interference=true;
    }

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
  else if (
    ndau>=4
    &&
    ((production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB) && process == TVar::bkgWW)
    ){
    // Processes 61 (4l), 62 (2l2q), 64 (2q2l)

    nqcdjets_.nqcdjets=0;
    nwz_.nwz=1; ckmfill_(&(nwz_.nwz));
    bveg1_mcfm_.ndim=10;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2=masses_mcfm_.wmass;
    breit_.width2=masses_mcfm_.wwidth;
    breit_.mass3=masses_mcfm_.wmass;
    breit_.width3=masses_mcfm_.wwidth;
    srdiags_.srdiags=true;

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
  else if (
    ndau>=4
    &&
    production == TVar::JJQCD
    &&
    process == TVar::bkgZJets
    ){
    // -- 44 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)+f(p6)'
    // these settings are identical to use the chooser_() function

    flags_.Gflag=true;
    flags_.Qflag=true;
    flags_.QandGflag=true; // This is in case NLO is implemented.
    bveg1_mcfm_.ndim=10;
    breit_.n2=0;
    breit_.n3=1;
    nqcdjets_.nqcdjets=2;
    nwz_.nwz=0; ckmfill_(&(nwz_.nwz));
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
  else if (
    ndau==3
    &&
    (production == TVar::ZZQQB || production == TVar::ZZINDEPENDENT)
    &&
    process == TVar::bkgZGamma
    ){
    // -- 300 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+gamma(p5)'
    // -- 305 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4)))-(sum over 3 nu)+gamma(p5)'

    nqcdjets_.nqcdjets=0;
    bveg1_mcfm_.ndim=7;
    breit_.n2=0;
    breit_.n3=1;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;
    nwz_.nwz=0; ckmfill_(&(nwz_.nwz));

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
  else if (
    ndau>=4
    &&
    production == TVar::ZZGG
    &&
    (isZZ && (process==TVar::bkgZZ || process==TVar::HSMHiggs || process == TVar::bkgZZ_SMHiggs))
    ){
    /*
    nprocs:
    c--- 128 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [top, bottom loops, exact]' 'L'
    c--- 129 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [only H, gg->ZZ intf.]' 'L' -> NOT IMPLEMENTED
    c--- 130 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [H squared and H, gg->ZZ intf.]' 'L'
    c--- 131 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6) [gg only, (H + gg->ZZ) squared]' 'L'
    c--- 132 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6) [(gg->ZZ) squared]' 'L'
    */

    nqcdjets_.nqcdjets=0;
    nwz_.nwz=0; ckmfill_(&(nwz_.nwz));
    bveg1_mcfm_.ndim=10;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2 =masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3 =masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    vsymfact_.vsymfact=1.0;
    interference_.interference=false;
    if (hasZZ4fInterf && (leptonInterf==TVar::DefaultLeptonInterf || leptonInterf==TVar::InterfOn)){
      vsymfact_.vsymfact=0.5;
      interference_.interference=true;
    }

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
  else if (
    ndau>=4
    &&
    production == TVar::ZZGG
    &&
    (
    (
    (isWW && (process==TVar::bkgWW || process==TVar::HSMHiggs || process == TVar::bkgWW_SMHiggs))
    )
    ||
    (
    ((isWW || isZZ) && (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs))
    )
    )
    ){ // gg->VV
    // Processes 1281, 1311, 1321

    nwz_.nwz=0; ckmfill_(&(nwz_.nwz));
    nqcdjets_.nqcdjets=0;
    bveg1_mcfm_.ndim=10;
    breit_.n2=1;
    breit_.n3=1;
    nuflav_.nuflav=1; // Keep this at 1. Mela controls how many flavors in a more exact way.

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
  // JJ + ZZ->4f
  else if ( // Check for support in qq'H+2J
    ndau>=4
    &&
    (
    production==TVar::Had_WH || production==TVar::Had_ZH
    || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
    || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
    || production==TVar::Lep_WH || production==TVar::Lep_ZH
    || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
    || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
    || production==TVar::JJVBF || production==TVar::JJEW
    || production==TVar::JJVBF_S || production==TVar::JJEW_S
    || production==TVar::JJVBF_TU || production==TVar::JJEW_TU
    )
    &&
    (isZZ && (process==TVar::bkgZZ || process==TVar::HSMHiggs || process == TVar::bkgZZ_SMHiggs))
    ){
    // 220 '  f(p1)+f(p2) --> Z(e-(p3),e^+(p4))Z(mu-(p5),mu+(p6)))+f(p7)+f(p8) [weak]' 'L'

    if (process==TVar::bkgZZ) sprintf(runstring_.runstring, "wbfBO");
    else if (process==TVar::HSMHiggs) sprintf(runstring_.runstring, "wbfHO");

    npart_.npart=6;
    nwz_.nwz=2; ckmfill_(&(nwz_.nwz));
    bveg1_mcfm_.ndim=16;
    nqcdjets_.nqcdjets=2;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2 =masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3 =masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    vsymfact_.vsymfact=1.0;
    interference_.interference=false;
    if (hasZZ4fInterf && (leptonInterf==TVar::DefaultLeptonInterf || leptonInterf==TVar::InterfOn)){
      vsymfact_.vsymfact=0.5;
      interference_.interference=true;
    }

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
  // JJ + WW->4f
  else if ( // Check for support in qq'H+2J
    ndau>=4
    &&
    (
    production==TVar::Had_WH || production==TVar::Had_ZH
    || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
    || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
    || production==TVar::Lep_WH || production==TVar::Lep_ZH
    || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
    || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
    || production==TVar::JJVBF || production==TVar::JJEW
    || production==TVar::JJVBF_S || production==TVar::JJEW_S
    || production==TVar::JJVBF_TU || production==TVar::JJEW_TU
    )
    &&
    (isWW && (process==TVar::bkgWW || process==TVar::HSMHiggs || process == TVar::bkgWW_SMHiggs))
    ){
    // Process 224

    if (process==TVar::bkgWW) sprintf(runstring_.runstring, "wbfBO");
    else if (process==TVar::HSMHiggs) sprintf(runstring_.runstring, "wbfHO");

    npart_.npart=6;
    nwz_.nwz=2; ckmfill_(&(nwz_.nwz));
    bveg1_mcfm_.ndim=16;
    nqcdjets_.nqcdjets=2;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2 =masses_mcfm_.wmass;
    breit_.width2=masses_mcfm_.wwidth;
    breit_.mass3 =masses_mcfm_.wmass;
    breit_.width3=masses_mcfm_.wwidth;

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
  // JJ + VV->4f
  else if ( // Check for support in qq'H+2J
    ndau>=4
    &&
    (
    production==TVar::Had_WH || production==TVar::Had_ZH
    || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
    || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
    || production==TVar::Lep_WH || production==TVar::Lep_ZH
    || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
    || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
    || production==TVar::JJVBF || production==TVar::JJEW
    || production==TVar::JJVBF_S || production==TVar::JJEW_S
    || production==TVar::JJVBF_TU || production==TVar::JJEW_TU
    )
    &&
    ((isZZ || isWW) && (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs))
    ){
    // Procsss 226

    if (process==TVar::bkgWWZZ) sprintf(runstring_.runstring, "wbfBO");
    else if (process==TVar::HSMHiggs_WWZZ) sprintf(runstring_.runstring, "wbfHO");

    npart_.npart=6;
    nwz_.nwz=2; ckmfill_(&(nwz_.nwz));
    bveg1_mcfm_.ndim=16;
    nqcdjets_.nqcdjets=2;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2 =masses_mcfm_.wmass;
    breit_.width2=masses_mcfm_.wwidth;
    breit_.mass3 =masses_mcfm_.wmass;
    breit_.width3=masses_mcfm_.wwidth;

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
  // JJ + ZZ->4f (QCD)
  else if (
    ndau>=4
    &&
    (
    production == TVar::JJQCD/* || production == TVar::JJQCD_S || production == TVar::JJQCD_TU*/
    || production == TVar::JJEWQCD || production == TVar::JJEWQCD_S || production == TVar::JJEWQCD_TU
    )
    &&
    (isZZ && process==TVar::bkgZZ) // ME is bkg-only
    ){
    // 2201 '  f(p1)+f(p2) --> Z(e-(p3),e^+(p4))Z(mu-(p5),mu+(p6)))+f(p7)+f(p8) [strong]' 'L'
    /*
    NOTE TO DEVELOPER
    MCFM also sets common/VVstrong/VVstrong (logical) here, but it only controls which ww_ZZqq function is called.
    It is irrelevant as long as eventgeneration through lowint is not used directly through MELA.
    */
    npart_.npart=6;
    nwz_.nwz=2; ckmfill_(&(nwz_.nwz));
    bveg1_mcfm_.ndim=16;
    nqcdjets_.nqcdjets=2;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2 =masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3 =masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    vsymfact_.vsymfact=1.0;
    interference_.interference=false;
    if (hasZZ4fInterf && (leptonInterf==TVar::DefaultLeptonInterf || leptonInterf==TVar::InterfOn)){
      vsymfact_.vsymfact=0.5;
      interference_.interference=true;
    }

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
    // JJ + WW->4f (QCD)
  else if (
    ndau>=4
    &&
    (
    production == TVar::JJQCD/* || production == TVar::JJQCD_S || production == TVar::JJQCD_TU*/
    || production == TVar::JJEWQCD || production == TVar::JJEWQCD_S || production == TVar::JJEWQCD_TU
    )
    &&
    (isWW && process==TVar::bkgWW) // ME is bkg-only
    ){
    // Process 2241
    /*
    NOTE TO DEVELOPER
    MCFM also sets common/VVstrong/VVstrong (logical) here, but it only controls which ww_ZZqq function is called.
    It is irrelevant as long as eventgeneration through lowint is not used directly through MELA.
    */
    npart_.npart=6;
    nwz_.nwz=2; ckmfill_(&(nwz_.nwz));
    bveg1_mcfm_.ndim=16;
    nqcdjets_.nqcdjets=2;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2 =masses_mcfm_.wmass;
    breit_.width2=masses_mcfm_.wwidth;
    breit_.mass3 =masses_mcfm_.wmass;
    breit_.width3=masses_mcfm_.wwidth;

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_chooser: Setup is (production, process)=(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << ")" << endl;

  }
    // JJ + VV->4f (QCD)
  else if (
    ndau>=4
    &&
    (
    production == TVar::JJQCD/* || production == TVar::JJQCD_S || production == TVar::JJQCD_TU*/
    || production == TVar::JJEWQCD || production == TVar::JJEWQCD_S || production == TVar::JJEWQCD_TU
    )
    &&
    ((isZZ || isWW) && process==TVar::bkgWWZZ)
    ){
    // Procsss 2261 (not supported yet)
    MELAerr << "TUtil::MCFM_chooser: MCFM does not support QCD JJ+VV->4f interaction yet. Please contact the MELA authors if you need more information." << endl;
    result = false;

  }
  else{
    MELAerr << "TUtil::MCFM_chooser: Can't identify (process, production) = (" << process << ", " << production << ")" << endl;
    MELAerr << "TUtil::MCFM_chooser: ndau: " << ndau << '\t';
    MELAerr << "TUtil::MCFM_chooser: isZZ: " << isZZ << '\t';
    MELAerr << "TUtil::MCFM_chooser: isWW: " << isWW << '\t';
    MELAerr << endl;
    result = false;
  }

  delete[] pId;
  return result;
}
bool TUtil::MCFM_SetupParticleCouplings(
  const TVar::Process& process, const TVar::Production& production,
  const TVar::VerbosityLevel& verbosity,
  const simple_event_record& mela_event,
  vector<int>* partOrder, vector<int>* apartOrder
  ){
  bool result=true;

  // Initialize Z couplings
  zcouple_.q1=0;
  zcouple_.l1=0;
  zcouple_.r1=0;
  zcouple_.q2=0;
  zcouple_.l2=0;
  zcouple_.r2=0;

  // Initialize plabels
  TString strplabel[mxpart];
  for (int ip=0; ip<mxpart; ip++) strplabel[ip]="  ";

  // Channel checks
  unsigned int ndau = mela_event.pDaughters.size();
  if (ndau<1) return false;
  unsigned int napart = mela_event.pAssociated.size();
  bool isWW = (ndau>=4 && PDGHelpers::isAWBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAWBoson(mela_event.intermediateVid.at(1)));
  bool isZZ = (ndau>=4 && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(1)));
  //bool hasZZ4fInterf = isZZ && abs(pId[0])==abs(pId[2]) && abs(pId[1])==abs(pId[3]) && !PDGHelpers::isAnUnknownJet(pId[0]) && !PDGHelpers::isAnUnknownJet(pId[3]);
  bool isZJJ = false;
  if (ndau>=4) isZJJ = PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)); // No check on whether daughters 2 and 3 are jets. Notice both isZZ and isZJJ could be true
  bool isZG = (ndau>=3 && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(1)));
  bool isGG = (ndau>=2 && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(0)) && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(1)));
  bool hasZ1 = (isZZ || isZG || isZJJ);
  bool hasZ2 = isZZ;
  bool hasW1 = isWW;
  bool hasW2 = isWW;
  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::MCFM_SetupParticleCouplings(" << TVar::ProductionName(production) << ", " << TVar::ProcessName(process) << "):\nInitial configuration is ";
    vector<TString> strcfgs;
    if (isGG) strcfgs.push_back(TString("GG"));
    if (isZG) strcfgs.push_back(TString("ZG"));
    if (isZZ) strcfgs.push_back(TString("ZZ"));
    if (isZJJ) strcfgs.push_back(TString("ZJJ"));
    if (isWW) strcfgs.push_back(TString("WW"));
    if (!strcfgs.empty()){
      for (unsigned int istr=0; istr<strcfgs.size(); istr++){
        if (istr==0) MELAout << strcfgs.at(istr);
        else MELAout << ", " << strcfgs.at(istr);
      }
      MELAout << ".";
    }
    else MELAout << "has no valid decay!";
    MELAout << endl;
  }

  // Setup the ordering arrays
  int* pApartOrder = 0;
  int* pApartId = 0;
  if (napart>0){
    pApartOrder = new int[napart];
    pApartId = new int[napart];
    for (unsigned int ip=0; ip<napart; ip++){
      pApartOrder[ip]=ip; // This order USUALLY does not change!
      pApartId[ip]=mela_event.pAssociated.at(ip).first; // This order does not change.
    }
  }
  int* pOrder = new int[ndau];
  int* pZOrder = new int[ndau];
  int* pWOrder = new int[ndau];
  int* pId = new int[ndau];
  for (unsigned int ip=0; ip<ndau; ip++){
    pOrder[ip]=ip; // This order does not change
    pZOrder[ip]=ip;
    pWOrder[ip]=ip;
    pId[ip]=mela_event.pDaughters.at(ip).first; // This order does not change.
  }

  // Special case checks
  bool useQQBZGAM = (isZG && process == TVar::bkgZGamma && (production == TVar::ZZQQB || production == TVar::ZZINDEPENDENT));
  bool useQQVVQQ =
    (
    ((isZZ && (process==TVar::bkgZZ || process==TVar::HSMHiggs || process == TVar::bkgZZ_SMHiggs))
    || (isWW && (process==TVar::bkgWW || process==TVar::HSMHiggs || process == TVar::bkgWW_SMHiggs))
    || ((isZZ || isWW) && (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs)))
    &&
    (production==TVar::Had_WH || production==TVar::Had_ZH
    || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
    || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
    || production==TVar::Lep_WH || production==TVar::Lep_ZH
    || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
    || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
    || production==TVar::JJVBF || production==TVar::JJEW
    || production==TVar::JJVBF_S || production==TVar::JJEW_S
    || production==TVar::JJVBF_TU || production==TVar::JJEW_TU)
    );
  bool useQQVVQQstrong =
    (
    ((isZZ && process==TVar::bkgZZ) || (isWW && process==TVar::bkgWW) || ((isZZ || isWW) && process==TVar::bkgWWZZ))
    &&
    (production == TVar::JJQCD || production == TVar::JJQCD_S || production == TVar::JJQCD_TU)
    );
  bool useQQVVQQboth =
    (
    ((isZZ && process==TVar::bkgZZ) || (isWW && process==TVar::bkgWW) || ((isZZ || isWW) && process==TVar::bkgWWZZ))
    &&
    (production == TVar::JJEWQCD || production == TVar::JJEWQCD_S || production == TVar::JJEWQCD_TU)
    );
  bool useQQVVQQany = useQQVVQQ || useQQVVQQstrong || useQQVVQQboth;

  /**************************/
  /* Begin the setup checks */
  /**************************/
  // Check first if the decay mode is valid and the number of associated particles is consistent with the requested number
  if (
    !(isWW || isZZ || isZJJ || isZG || isGG) // Only ZZ, WW, ZG or GG supported in MCFM
    ||
    ((int)napart<(mela_event.nRequested_AssociatedJets+mela_event.nRequested_AssociatedLeptons)) // Associated particle not found
    ) result=false;
  else{

    if (isWW && !hasZ1 && !hasZ2){
      // Default swap is 3-5
      // For ggVV amplitudes, MCFM generates W+W- and does 3-5 swap inside the ME, so by definition, this is ok.
      // For VBFWW amplitudes, MCFM generates W-W+ in the phase space and does a 4-6 swap before passing to the ME. The ordering here is W+W-, so this is equivalent to doing a 3-5 swap to get the corresponding ZZ-like combinations.
      swap(pZOrder[0], pZOrder[2]);

      int V1id, V2id;
      if (!PDGHelpers::isAnUnknownJet(pId[pZOrder[0]]) && !PDGHelpers::isAnUnknownJet(pId[pZOrder[1]])) V1id = PDGHelpers::getCoupledVertex(pId[pZOrder[0]], pId[pZOrder[1]]);
      else V1id = 23;
      if (!PDGHelpers::isAnUnknownJet(pId[pZOrder[2]]) && !PDGHelpers::isAnUnknownJet(pId[pZOrder[3]])) V2id = PDGHelpers::getCoupledVertex(pId[pZOrder[2]], pId[pZOrder[3]]);
      else V2id = 23;
      hasZ1=hasZ1 || PDGHelpers::isAZBoson(V1id);
      hasZ2=hasZ2 || PDGHelpers::isAZBoson(V2id);
      if (!hasZ1 && hasZ2){
        swap(pZOrder[0], pZOrder[2]);
        swap(pZOrder[1], pZOrder[3]);
        swap(hasZ1, hasZ2);
      }
      if (verbosity>=TVar::DEBUG){
        if (hasZ1) MELAout << "TUtil::MCFM_SetupParticleCouplings: Found a Z1";
        if (hasZ2) MELAout << " and a Z2";
        if (hasZ1 || hasZ2) MELAout << " in WW." << endl;
      }
    }
    if (isZZ && !hasW1 && !hasW2){
      swap(pWOrder[0], pWOrder[2]);

      int V1id, V2id;
      if (!PDGHelpers::isAnUnknownJet(pId[pWOrder[0]]) && !PDGHelpers::isAnUnknownJet(pId[pWOrder[1]])) V1id = PDGHelpers::getCoupledVertex(pId[pWOrder[0]], pId[pWOrder[1]]);
      else V1id = 24;
      if (!PDGHelpers::isAnUnknownJet(pId[pWOrder[2]]) && !PDGHelpers::isAnUnknownJet(pId[pWOrder[3]])) V2id = PDGHelpers::getCoupledVertex(pId[pWOrder[2]], pId[pWOrder[3]]);
      else V2id = -24;
      if (PDGHelpers::isAWBoson(V1id) && PDGHelpers::isAWBoson(V2id)){
        if (V1id<0 || V2id>0){
          swap(pWOrder[0], pWOrder[2]);
          swap(pWOrder[1], pWOrder[3]);
        }
        hasW1=true;
        hasW2=true;
      }
      if (verbosity>=TVar::DEBUG){
        if (hasW1) MELAout << "TUtil::MCFM_SetupParticleCouplings: Found a W1(" << V1id << ")";
        if (hasW2) MELAout << " and a W2(" << V2id << ")";
        if (hasW1 || hasW2) MELAout << " in ZZ." << endl;
      }
    }

    /*******************/
    /* Particle labels */
    /*******************/
    // Mother particles
    // Default labels for most processes
    strplabel[0]="pp";
    strplabel[1]="pp";
    if (useQQVVQQany){ // Special case if using qqZZ/VVqq*
      if (verbosity>=TVar::DEBUG) MELAout << "TUtil::MCFM_SetupParticleCouplings: Setting up mother labels for MCFM:";
      for (int ip=0; ip<min(2, (int)mela_event.pMothers.size()); ip++){
        const int& idmot = mela_event.pMothers.at(ip).first;
        if (!PDGHelpers::isAnUnknownJet(idmot)) strplabel[ip]=TUtil::GetMCFMParticleLabel(idmot, false, useQQVVQQany);
        if (verbosity>=TVar::DEBUG) MELAout << " " << idmot << "=" << strplabel[ip];
        // No need to check unknown parton case, already "pp"
      }
      if (verbosity>=TVar::DEBUG) MELAout << endl;
    }

    // Decay and associated particles
    // 0-jet processes, set only decay labels (and plabel[5/6]=pp due to NLO stuff)
    if (isZJJ && production == TVar::JJQCD && process == TVar::bkgZJets){
      strplabel[2]="el";
      strplabel[3]="ea";
      strplabel[4]="pp";
      strplabel[5]="pp";
      strplabel[6]="pp";
    }
    else if (production == TVar::ZZGG && (
      (isWW && (process==TVar::bkgWW || process==TVar::HSMHiggs || process == TVar::bkgWW_SMHiggs))
      ||
      ((isWW || isZZ) && (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs))
      )
      ){
      strplabel[2]="el";
      strplabel[3]="ea";
      strplabel[4]="nl";
      strplabel[5]="na";
      strplabel[6]="pp";

      string strrun = runstring_.runstring;
      if (isWW && ((!hasZ1 || !hasZ2) || (process==TVar::bkgWW || process==TVar::HSMHiggs || process == TVar::bkgWW_SMHiggs))) strrun += "_ww";
      else if (isZZ && (!hasW1 || !hasW2)) strrun += "_zz";
      sprintf(runstring_.runstring, strrun.c_str());
    }
    else if (isZG && (production == TVar::ZZQQB || production == TVar::ZZINDEPENDENT) && process == TVar::bkgZGamma){
      lastphot_.lastphot=5;
      strplabel[4]="ga";
      strplabel[5]="pp";
    }
    else if (isWW && (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB) && process == TVar::bkgWW){
      strplabel[6]="pp";
      zcouple_.l1=1.;

      if (PDGHelpers::isANeutrino(pId[pWOrder[0]])){
        strplabel[2]="nl";
        strplabel[3]="ea";
      }
      else if (PDGHelpers::isAJet(pId[pWOrder[1]])){
        strplabel[2]="qj";
        strplabel[3]="qj";
        zcouple_.l1 *= sqrt(6.);
        nqcdjets_.nqcdjets += 2;
      }
      else result = false;

      if (PDGHelpers::isANeutrino(pId[pWOrder[3]])){
        strplabel[4]="el";
        strplabel[5]="na";
      }
      else if (PDGHelpers::isAJet(pId[pWOrder[3]])){
        strplabel[4]="qj";
        strplabel[5]="qj";
        zcouple_.l1 *= sqrt(6.);
        nqcdjets_.nqcdjets += 2;
      }
      else result = false;
      if (PDGHelpers::isAJet(pId[pWOrder[0]]) && PDGHelpers::isAJet(pId[pWOrder[3]])) result = false; // MCFM does not support WW->4q
    }
    // 2-jet processes
    // Most of the associated particle labels are set below
    // VH productions loop over all possible associated particles to assign V daughters as the first two particles in particle-antiparticle order
    else if ((isWW || isZZ) && napart>=2 && (production==TVar::Lep_ZH || production==TVar::Lep_ZH_S || production==TVar::Lep_ZH_TU)){
      spinzerohiggs_anomcoupl_.channeltoggle_stu = int(production==TVar::Lep_ZH)*2 + int(production==TVar::Lep_ZH_TU)*1 + int(production==TVar::Lep_ZH_S)*0;
      spinzerohiggs_anomcoupl_.vvhvvtoggle_vbfvh = 1;

      if (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs){
        string strrun = runstring_.runstring;
        if (isWW && (!hasZ1 || !hasZ2)) strrun += "_ww";
        else if (isZZ && (!hasW1 || !hasW2)) strrun += "_zz";
        sprintf(runstring_.runstring, strrun.c_str());
      }

      bool hasZll=false;
      bool hasZnn=false;
      for (unsigned int ix=0; ix<napart; ix++){
        if (PDGHelpers::isAJet(pApartId[ix])) continue;
        for (unsigned int iy=ix+1; iy<napart; iy++){
          if (PDGHelpers::isAJet(pApartId[iy])) continue;
          int Vid = PDGHelpers::getCoupledVertex(pApartId[ix], pApartId[iy]);
          if (PDGHelpers::isAZBoson(Vid) && PDGHelpers::isALepton(pApartId[ix])){
            int i1=ix; int i2=iy;
            if (pApartId[ix]<0){ swap(i1, i2); }
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasZll=true;
            break;
          }
          else if (PDGHelpers::isAZBoson(Vid) && PDGHelpers::isANeutrino(pApartId[ix])){
            int i1=ix; int i2=iy;
            if (pApartId[ix]<0){ swap(i1, i2); }
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasZnn=true;
            break;
          }
          if (hasZll || hasZnn) break;
        }
      }
      unsigned int iout=6;
      unsigned int jout=7;
      if (useQQVVQQany){
        int qqvvqq_apartordering[2]={ -1, -1 };
        TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(
          mela_event.pMothers.at(0).first, mela_event.pMothers.at(1).first,
          pApartId[pApartOrder[0]], pApartId[pApartOrder[1]],
          qqvvqq_apartordering
          );
        if (qqvvqq_apartordering[0]==-1 || qqvvqq_apartordering[1]==-1) result=false;
        else if (qqvvqq_apartordering[0]>qqvvqq_apartordering[1]){
          swap(iout, jout);
          swap(pApartOrder[0], pApartOrder[1]);
        }
      }
      if (hasZll){
        strplabel[iout]="el";
        strplabel[jout]="ea";
      }
      else if (hasZnn){
        strplabel[iout]="nl";
        strplabel[jout]="na";
      }
      else result = false;
    }
    else if ((isWW || isZZ) && napart>=2 && (production==TVar::Lep_WH || production==TVar::Lep_WH_S || production==TVar::Lep_WH_TU)){
      spinzerohiggs_anomcoupl_.channeltoggle_stu = int(production==TVar::Lep_WH)*2 + int(production==TVar::Lep_WH_TU)*1 + int(production==TVar::Lep_WH_S)*0;
      spinzerohiggs_anomcoupl_.vvhvvtoggle_vbfvh = 1;

      if (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs){
        string strrun = runstring_.runstring;
        if (isWW && (!hasZ1 || !hasZ2)) strrun += "_ww";
        else if (isZZ && (!hasW1 || !hasW2)) strrun += "_zz";
        sprintf(runstring_.runstring, strrun.c_str());
      }

      bool hasWplus=false;
      bool hasWminus=false;
      for (unsigned int ix=0; ix<napart; ix++){
        if (PDGHelpers::isAJet(pApartId[ix])) continue;
        for (unsigned int iy=ix+1; iy<napart; iy++){
          if (PDGHelpers::isAJet(pApartId[iy])) continue;
          int Vid = PDGHelpers::getCoupledVertex(pApartId[ix], pApartId[iy]);
          if (
            PDGHelpers::isAWBoson(Vid)
            &&
            (
            (PDGHelpers::isALepton(pApartId[ix]) && pApartId[ix]<0)
            ||
            (PDGHelpers::isALepton(pApartId[iy]) && pApartId[iy]<0)
            )
            ){
            int i1=ix; int i2=iy;
            if (pApartId[ix]<0){ swap(i1, i2); }
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasWplus=true;
            break;
          }
          else if (
            PDGHelpers::isAWBoson(Vid)
            &&
            (
            (PDGHelpers::isALepton(pApartId[ix]) && pApartId[ix]>0)
            ||
            (PDGHelpers::isALepton(pApartId[iy]) && pApartId[iy]>0)
            )
            ){
            int i1=ix; int i2=iy;
            if (pApartId[ix]<0){ swap(i1, i2); }
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasWminus=true;
            break;
          }
        }
        if (hasWplus || hasWminus) break;
      }
      unsigned int iout=6;
      unsigned int jout=7;
      if (useQQVVQQany){
        int qqvvqq_apartordering[2]={ -1, -1 };
        TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(
          mela_event.pMothers.at(0).first, mela_event.pMothers.at(1).first,
          pApartId[pApartOrder[0]], pApartId[pApartOrder[1]],
          qqvvqq_apartordering
          );
        if (qqvvqq_apartordering[0]==-1 || qqvvqq_apartordering[1]==-1) result=false;
        else if (qqvvqq_apartordering[0]>qqvvqq_apartordering[1]){
          swap(iout, jout);
          swap(pApartOrder[0], pApartOrder[1]);
        }
      }
      if (hasWplus){ // W+
        strplabel[iout]="nl";
        strplabel[jout]="ea";
      }
      else if (hasWminus){ // W-
        strplabel[iout]="el";
        strplabel[jout]="na";
      }
      else result = false;
    }
    else if ((isWW || isZZ) && napart>=2 && (production==TVar::Had_ZH || production==TVar::Had_ZH_S || production==TVar::Had_ZH_TU)){
      spinzerohiggs_anomcoupl_.channeltoggle_stu = int(production==TVar::Had_ZH)*2 + int(production==TVar::Had_ZH_TU)*1 + int(production==TVar::Had_ZH_S)*0;
      spinzerohiggs_anomcoupl_.vvhvvtoggle_vbfvh = 1;

      if (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs){
        string strrun = runstring_.runstring;
        if (isWW && (!hasZ1 || !hasZ2)) strrun += "_ww";
        else if (isZZ && (!hasW1 || !hasW2)) strrun += "_zz";
        sprintf(runstring_.runstring, strrun.c_str());
      }

      bool hasZuu=false;
      bool hasZdd=false;
      bool hasZjj=false;
      for (unsigned int ix=0; ix<napart; ix++){
        if (!PDGHelpers::isAJet(pApartId[ix])) continue;
        for (unsigned int iy=ix+1; iy<napart; iy++){
          if (!PDGHelpers::isAJet(pApartId[iy])) continue;
          int Vid;
          if (!PDGHelpers::isAnUnknownJet(pApartId[ix]) && !PDGHelpers::isAnUnknownJet(pApartId[iy])) Vid = PDGHelpers::getCoupledVertex(pApartId[ix], pApartId[iy]);
          else Vid = 23;
          if (PDGHelpers::isAZBoson(Vid) && (PDGHelpers::isUpTypeQuark(pApartId[ix]) || PDGHelpers::isUpTypeQuark(pApartId[iy]))){
            int i1=ix; int i2=iy;
            if (pApartId[ix]<0 || pApartId[iy]>0){ swap(i1, i2); }
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasZuu=true;
            break;
          }
          else if (PDGHelpers::isAZBoson(Vid) && (PDGHelpers::isDownTypeQuark(pApartId[ix]) || PDGHelpers::isDownTypeQuark(pApartId[iy]))){
            int i1=ix; int i2=iy;
            if (pApartId[ix]<0 || pApartId[iy]>0){ swap(i1, i2); }
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasZdd=true;
            break;
          }
          else if (PDGHelpers::isAZBoson(Vid)){
            int i1=ix; int i2=iy;
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasZjj=true;
            break;
          }
          if (hasZuu || hasZdd || hasZjj) break;
        }
      }
      if (useQQVVQQany){
        int qqvvqq_apartordering[2]={ -1, -1 };
        TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(
          mela_event.pMothers.at(0).first, mela_event.pMothers.at(1).first,
          pApartId[pApartOrder[0]], pApartId[pApartOrder[1]],
          qqvvqq_apartordering
          );
        if (qqvvqq_apartordering[0]==-1 || qqvvqq_apartordering[1]==-1) result=false;
        else if (qqvvqq_apartordering[0]>qqvvqq_apartordering[1]) swap(pApartOrder[0], pApartOrder[1]);
      }
      if (hasZuu || hasZdd){
        strplabel[6]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[0]], true, useQQVVQQany);
        strplabel[7]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[1]], true, useQQVVQQany);
      }
      else if (hasZjj){
        strplabel[6]="qj";
        strplabel[7]="qj";
      }
      else result = false;
    }
    else if ((isWW || isZZ) && napart>=2 && (production==TVar::Had_WH || production==TVar::Had_WH_S || production==TVar::Had_WH_TU)){
      spinzerohiggs_anomcoupl_.channeltoggle_stu = int(production==TVar::Had_WH)*2 + int(production==TVar::Had_WH_TU)*1 + int(production==TVar::Had_WH_S)*0;
      spinzerohiggs_anomcoupl_.vvhvvtoggle_vbfvh = 1;

      if (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs){
        string strrun = runstring_.runstring;
        if (isWW && (!hasZ1 || !hasZ2)) strrun += "_ww";
        else if (isZZ && (!hasW1 || !hasW2)) strrun += "_zz";
        sprintf(runstring_.runstring, strrun.c_str());
      }

      bool hasWplus=false;
      bool hasWminus=false;
      bool hasWjj=false;
      for (unsigned int ix=0; ix<napart; ix++){
        if (!PDGHelpers::isAJet(pApartId[ix])) continue;
        for (unsigned int iy=ix+1; iy<napart; iy++){
          if (!PDGHelpers::isAJet(pApartId[iy])) continue;
          int Vid;
          if (!PDGHelpers::isAnUnknownJet(pApartId[ix]) && !PDGHelpers::isAnUnknownJet(pApartId[iy])) Vid = PDGHelpers::getCoupledVertex(pApartId[ix], pApartId[iy]);
          else if (
            (PDGHelpers::isUpTypeQuark(pApartId[ix]) && pApartId[ix]<0) || (PDGHelpers::isUpTypeQuark(pApartId[iy]) && pApartId[iy]<0) ||
            (PDGHelpers::isDownTypeQuark(pApartId[ix]) && pApartId[ix]>0) || (PDGHelpers::isDownTypeQuark(pApartId[iy]) && pApartId[iy]>0)
            ) Vid=-24;
          else Vid = 24;
          if (
            PDGHelpers::isAWBoson(Vid)
            &&
            (
            ((PDGHelpers::isDownTypeQuark(pApartId[ix]) && pApartId[ix]<0) || (PDGHelpers::isUpTypeQuark(pApartId[iy]) && pApartId[iy]>0))
            ||
            ((PDGHelpers::isDownTypeQuark(pApartId[iy]) && pApartId[iy]<0) || (PDGHelpers::isUpTypeQuark(pApartId[ix]) && pApartId[ix]>0))
            )
            ){
            int i1=ix; int i2=iy;
            if (pApartId[ix]<0 || pApartId[iy]>0){ swap(i1, i2); }
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasWplus=true;
            break;
          }
          else if (
            PDGHelpers::isAWBoson(Vid)
            &&
            (
            ((PDGHelpers::isDownTypeQuark(pApartId[ix]) && pApartId[ix]>0) || (PDGHelpers::isUpTypeQuark(pApartId[iy]) && pApartId[iy]<0))
            ||
            ((PDGHelpers::isDownTypeQuark(pApartId[iy]) && pApartId[iy]>0) || (PDGHelpers::isUpTypeQuark(pApartId[ix]) && pApartId[ix]<0))
            )
            ){
            int i1=ix; int i2=iy;
            if (pApartId[ix]<0 || pApartId[iy]>0){ swap(i1, i2); }
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasWminus=true;
            break;
          }
          else if (PDGHelpers::isAWBoson(Vid)){
            int i1=ix; int i2=iy;
            int* tmpOrder = new int[napart];
            tmpOrder[0] = i1; tmpOrder[1] = i2;
            unsigned int ctr=2;
            for (unsigned int ipart=0; ipart<napart; ipart++){ if (ctr<napart && pApartOrder[ipart]!=i1 && pApartOrder[ipart]!=i2){ tmpOrder[ctr] = pApartOrder[ipart]; ctr++; } }
            delete[] pApartOrder; pApartOrder = tmpOrder;
            hasWjj=true;
            break;
          }
        }
        if (hasWplus || hasWminus || hasWjj) break;
      }
      if (useQQVVQQany){
        int qqvvqq_apartordering[2]={ -1, -1 };
        TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(
          mela_event.pMothers.at(0).first, mela_event.pMothers.at(1).first,
          pApartId[pApartOrder[0]], pApartId[pApartOrder[1]],
          qqvvqq_apartordering
          );
        if (qqvvqq_apartordering[0]==-1 || qqvvqq_apartordering[1]==-1) result=false;
        else if (qqvvqq_apartordering[0]>qqvvqq_apartordering[1]) swap(pApartOrder[0], pApartOrder[1]);
      }
      if (hasWplus || hasWminus){ // W+ or W-
        strplabel[6]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[0]], true, useQQVVQQany);
        strplabel[7]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[1]], true, useQQVVQQany);
      }
      else if (hasWjj){ // W+/-
        strplabel[6]="qj";
        strplabel[7]="qj";
      }
      else result = false;
    }
    else if ((isWW || isZZ) && napart>=2 && (production==TVar::JJQCD || production==TVar::JJQCD_S || production==TVar::JJQCD_TU)){
      spinzerohiggs_anomcoupl_.channeltoggle_stu = int(production==TVar::JJQCD)*2 + int(production==TVar::JJQCD_TU)*1 + int(production==TVar::JJQCD_S)*0;;
      spinzerohiggs_anomcoupl_.vvhvvtoggle_vbfvh = 2; // Doesn't matter, already JJQCD

      if (process==TVar::bkgWWZZ/* || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs*/){
        string strrun = runstring_.runstring;
        if (isWW && (!hasZ1 || !hasZ2)) strrun += "_ww";
        else if (isZZ && (!hasW1 || !hasW2)) strrun += "_zz";
        sprintf(runstring_.runstring, strrun.c_str());
      }

      if (useQQVVQQany){
        int qqvvqq_apartordering[2]={ -1, -1 };
        TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(
          mela_event.pMothers.at(0).first, mela_event.pMothers.at(1).first,
          pApartId[pApartOrder[0]], pApartId[pApartOrder[1]],
          qqvvqq_apartordering
          );
        if (qqvvqq_apartordering[0]==-1 || qqvvqq_apartordering[1]==-1) result=false;
        else if (qqvvqq_apartordering[0]>qqvvqq_apartordering[1]) swap(pApartOrder[0], pApartOrder[1]);
      }
      strplabel[6]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[0]], false, useQQVVQQany);
      strplabel[7]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[1]], false, useQQVVQQany);
    }
    else if ((isWW || isZZ) && napart>=2 && (production==TVar::JJVBF || production==TVar::JJVBF_S || production==TVar::JJVBF_TU)){
      spinzerohiggs_anomcoupl_.channeltoggle_stu = int(production==TVar::JJVBF)*2 + int(production==TVar::JJVBF_TU)*1 + int(production==TVar::JJVBF_S)*0;;
      spinzerohiggs_anomcoupl_.vvhvvtoggle_vbfvh = 0; // VBF-only process

      if (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs){
        string strrun = runstring_.runstring;
        if (isWW && (!hasZ1 || !hasZ2)) strrun += "_ww";
        else if (isZZ && (!hasW1 || !hasW2)) strrun += "_zz";
        sprintf(runstring_.runstring, strrun.c_str());
      }

      if (useQQVVQQany){
        int qqvvqq_apartordering[2]={ -1, -1 };
        TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(
          mela_event.pMothers.at(0).first, mela_event.pMothers.at(1).first,
          pApartId[pApartOrder[0]], pApartId[pApartOrder[1]],
          qqvvqq_apartordering
          );
        if (qqvvqq_apartordering[0]==-1 || qqvvqq_apartordering[1]==-1) result=false;
        else if (qqvvqq_apartordering[0]>qqvvqq_apartordering[1]) swap(pApartOrder[0], pApartOrder[1]);
      }
      strplabel[6]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[0]], true, useQQVVQQany);
      strplabel[7]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[1]], true, useQQVVQQany);
    }
    else if ((isWW || isZZ) && napart>=2 && (production==TVar::JJEW || production==TVar::JJEWQCD || production==TVar::JJEW_S || production==TVar::JJEWQCD_S || production==TVar::JJEW_TU || production==TVar::JJEWQCD_TU)){
      spinzerohiggs_anomcoupl_.channeltoggle_stu = int(production==TVar::JJEWQCD || production==TVar::JJEW)*2 + int(production==TVar::JJEWQCD_TU || production==TVar::JJEW_TU)*1 + int(production==TVar::JJEWQCD_S || production==TVar::JJEW_S)*0;;
      spinzerohiggs_anomcoupl_.vvhvvtoggle_vbfvh = 2; // VBF+VH process

      if (process==TVar::bkgWWZZ || process==TVar::HSMHiggs_WWZZ || process == TVar::bkgWWZZ_SMHiggs){
        string strrun = runstring_.runstring;
        if (isWW && (!hasZ1 || !hasZ2)) strrun += "_ww";
        else if (isZZ && (!hasW1 || !hasW2)) strrun += "_zz";
        sprintf(runstring_.runstring, strrun.c_str());
      }

      if (useQQVVQQany){
        int qqvvqq_apartordering[2]={ -1, -1 };
        TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(
          mela_event.pMothers.at(0).first, mela_event.pMothers.at(1).first,
          pApartId[pApartOrder[0]], pApartId[pApartOrder[1]],
          qqvvqq_apartordering
          );
        if (qqvvqq_apartordering[0]==-1 || qqvvqq_apartordering[1]==-1) result=false;
        else if (qqvvqq_apartordering[0]>qqvvqq_apartordering[1]) swap(pApartOrder[0], pApartOrder[1]);
      }
      strplabel[6]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[0]], false, useQQVVQQany);
      strplabel[7]=TUtil::GetMCFMParticleLabel(pApartId[pApartOrder[1]], false, useQQVVQQany);
    }

    // Couplings for Z1
    if (isWW && (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB) && process == TVar::bkgWW){} // Skip this one, already handled above
    else if (hasZ1){
      if (PDGHelpers::isALepton(pId[pZOrder[0]]) && PDGHelpers::isALepton(pId[pZOrder[1]])){
        if (verbosity>=TVar::DEBUG) MELAout << "- Setting Z1->ll couplings." << endl;

        zcouple_.q1=-1.0;
        zcouple_.l1=zcouple_.le;
        zcouple_.r1=zcouple_.re;

        // Special Z1->ll cases
        if (useQQBZGAM){
          strplabel[2]="el";
          strplabel[3]="ea";
        }
        else if (useQQVVQQany){
          strplabel[2]=TUtil::GetMCFMParticleLabel(pId[pZOrder[0]], false, useQQVVQQany);
          strplabel[3]=TUtil::GetMCFMParticleLabel(pId[pZOrder[1]], false, useQQVVQQany);
        }
        // End special Z1->ll cases
      }
      else if (PDGHelpers::isANeutrino(pId[pZOrder[0]]) && PDGHelpers::isANeutrino(pId[pZOrder[1]])){
        if (verbosity>=TVar::DEBUG) MELAout << "- Setting Z1->nn couplings." << endl;

        zcouple_.q1=0;
        zcouple_.l1=zcouple_.ln;
        zcouple_.r1=zcouple_.rn;

        // Special Z1->nn cases
        if (useQQBZGAM){
          strplabel[2]="nl";
          strplabel[3]="na";
        }
        else if (useQQVVQQany){
          strplabel[2]=TUtil::GetMCFMParticleLabel(pId[pZOrder[0]], false, useQQVVQQany);
          strplabel[3]=TUtil::GetMCFMParticleLabel(pId[pZOrder[1]], false, useQQVVQQany);
        }
        // End special Z1->nn cases
      }
      else if (PDGHelpers::isAJet(pId[pZOrder[0]]) && PDGHelpers::isAJet(pId[pZOrder[1]])){
        if (verbosity>=TVar::DEBUG) MELAout << "- Setting Z1->jj couplings." << endl;

        nqcdjets_.nqcdjets += 2;
        int jetid=(PDGHelpers::isAnUnknownJet(pId[pZOrder[0]]) ? abs(pId[pZOrder[1]]) : abs(pId[pZOrder[0]]));
        if (!PDGHelpers::isAnUnknownJet(jetid)){
          if (jetid==6) jetid = 2;
          zcouple_.q1=ewcharge_.Q[5+jetid];
          zcouple_.l1=zcouple_.l[-1+jetid];
          zcouple_.r1=zcouple_.r[-1+jetid];

          // Special Z1->qq cases
          if (useQQBZGAM){
            strplabel[2]="el";// Trick MCFM, not implemented properly
            strplabel[3]="ea";
          }
          else if (useQQVVQQany){
            strplabel[2]=TUtil::GetMCFMParticleLabel(pId[pZOrder[0]], false, useQQVVQQany);
            strplabel[3]=TUtil::GetMCFMParticleLabel(pId[pZOrder[1]], false, useQQVVQQany);
          }
        }
        else{
          double gV_up = (zcouple_.l[-1+2]+zcouple_.r[-1+2])/2.;
          double gV_dn = (zcouple_.l[-1+1]+zcouple_.r[-1+1])/2.;
          double gA_up = (zcouple_.l[-1+2]-zcouple_.r[-1+2])/2.;
          double gA_dn = (zcouple_.l[-1+1]-zcouple_.r[-1+1])/2.;
          double yy_up = pow(gV_up, 2) + pow(gA_up, 2);
          double yy_dn = pow(gV_dn, 2) + pow(gA_dn, 2);
          double xx_up = gV_up*gA_up;
          double xx_dn = gV_dn*gA_dn;
          double yy = (2.*yy_up+3.*yy_dn);
          double xx = (2.*xx_up+3.*xx_dn);
          double discriminant = pow(yy, 2)-4.*pow(xx, 2);
          double gVsq = (yy+sqrt(fabs(discriminant)))/2.;
          double gAsq = pow(xx, 2)/gVsq;
          double gV=-sqrt(gVsq);
          double gA=-sqrt(gAsq);
          zcouple_.l1 = (gV+gA);
          zcouple_.r1 = (gV-gA);
          zcouple_.q1=sqrt(pow(ewcharge_.Q[5+1], 2)*3.+pow(ewcharge_.Q[5+2], 2)*2.);

          // Special Z1->qq cases
          if (useQQBZGAM){
            strplabel[2]="el";// Trick MCFM, not implemented properly
            strplabel[3]="ea";
          }
          else if (useQQVVQQany){
            strplabel[2]="qj";
            strplabel[3]="qj";
          }
        }
        if (!useQQVVQQany){ // colfac34_56 handles all couplings instead of this simple scaling (Reason: WW final states)
          zcouple_.l1 *= sqrt(3.);
          zcouple_.r1 *= sqrt(3.);
          zcouple_.q1 *= sqrt(3.);
        }
      } // End Z1 daughter id tests
    } // End ZZ/ZG/ZJJ Z1 couplings
    else if (useQQVVQQany){
      strplabel[2]=TUtil::GetMCFMParticleLabel(pId[pZOrder[0]], false, useQQVVQQany);
      strplabel[3]=TUtil::GetMCFMParticleLabel(pId[pZOrder[1]], false, useQQVVQQany);
    }

    // Couplings for Z2
    if (
      (isWW && (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB) && process == TVar::bkgWW) // Skip this one, already handled above
      ||
      (isZJJ && production == TVar::JJQCD && process == TVar::bkgZJets) // No need to handle
      ){}
    else if (hasZ2){
      if (PDGHelpers::isALepton(pId[pZOrder[2]]) && PDGHelpers::isALepton(pId[pZOrder[3]])){
        if (verbosity>=TVar::DEBUG) MELAout << "- Setting Z2->ll couplings." << endl;
        zcouple_.q2=-1.0;
        zcouple_.l2=zcouple_.le;
        zcouple_.r2=zcouple_.re;

        // Special Z2->ll cases
        if (useQQVVQQany){
          strplabel[4]=TUtil::GetMCFMParticleLabel(pId[pZOrder[2]], false, useQQVVQQany);
          strplabel[5]=TUtil::GetMCFMParticleLabel(pId[pZOrder[3]], false, useQQVVQQany);
        }
        // End special Z2->ll cases
      }
      else if (PDGHelpers::isANeutrino(pId[pZOrder[2]]) && PDGHelpers::isANeutrino(pId[pZOrder[3]])){
        if (verbosity>=TVar::DEBUG) MELAout << "- Setting Z2->nn couplings." << endl;
        zcouple_.q2=0;
        zcouple_.l2=zcouple_.ln;
        zcouple_.r2=zcouple_.rn;

        // Special Z2->nn cases
        if (useQQVVQQany){
          strplabel[4]=TUtil::GetMCFMParticleLabel(pId[pZOrder[2]], false, useQQVVQQany);
          strplabel[5]=TUtil::GetMCFMParticleLabel(pId[pZOrder[3]], false, useQQVVQQany);
        }
        // End special Z2->nn cases
      }
      else if (PDGHelpers::isAJet(pId[pZOrder[2]]) && PDGHelpers::isAJet(pId[pZOrder[3]])){
        if (verbosity>=TVar::DEBUG) MELAout << "- Setting Z2->jj couplings." << endl;

        nqcdjets_.nqcdjets += 2;
        int jetid=(PDGHelpers::isAnUnknownJet(pId[pZOrder[2]]) ? abs(pId[pZOrder[3]]) : abs(pId[pZOrder[2]]));
        if (!PDGHelpers::isAnUnknownJet(jetid)){
          if (jetid==6) jetid = 2;
          zcouple_.q2=ewcharge_.Q[5+jetid];
          zcouple_.l2=zcouple_.l[-1+jetid];
          zcouple_.r2=zcouple_.r[-1+jetid];

          // Special Z2->qq cases
          if (useQQVVQQany){
            strplabel[4]=TUtil::GetMCFMParticleLabel(pId[pZOrder[2]], false, useQQVVQQany);
            strplabel[5]=TUtil::GetMCFMParticleLabel(pId[pZOrder[3]], false, useQQVVQQany);
          }
        }
        else{
          double gV_up = (zcouple_.l[-1+2]+zcouple_.r[-1+2])/2.;
          double gV_dn = (zcouple_.l[-1+1]+zcouple_.r[-1+1])/2.;
          double gA_up = (zcouple_.l[-1+2]-zcouple_.r[-1+2])/2.;
          double gA_dn = (zcouple_.l[-1+1]-zcouple_.r[-1+1])/2.;
          double yy_up = pow(gV_up, 2) + pow(gA_up, 2);
          double yy_dn = pow(gV_dn, 2) + pow(gA_dn, 2);
          double xx_up = gV_up*gA_up;
          double xx_dn = gV_dn*gA_dn;
          double yy = (2.*yy_up+3.*yy_dn);
          double xx = (2.*xx_up+3.*xx_dn);
          double discriminant = pow(yy, 2)-4.*pow(xx, 2);
          double gVsq = (yy+sqrt(fabs(discriminant)))/2.;
          double gAsq = pow(xx, 2)/gVsq;
          double gV=-sqrt(gVsq);
          double gA=-sqrt(gAsq);
          zcouple_.l2 = (gV+gA);
          zcouple_.r2 = (gV-gA);
          zcouple_.q2=sqrt(pow(ewcharge_.Q[5+1], 2)*3.+pow(ewcharge_.Q[5+2], 2)*2.);

          // Special Z2->qq cases
          if (useQQVVQQany){
            strplabel[4]="qj";
            strplabel[5]="qj";
          }
        }
        if (!useQQVVQQany){ // colfac34_56 handles all couplings instead of this simple scaling (Reason: WW Final states)
          zcouple_.l2 *= sqrt(3.);
          zcouple_.r2 *= sqrt(3.);
          zcouple_.q2 *= sqrt(3.);
        }
      } // End Z2 daughter id tests
    } // End ZZ Z2 couplings
    else if (useQQVVQQany){
      strplabel[4]=TUtil::GetMCFMParticleLabel(pId[pZOrder[2]], false, useQQVVQQany);
      strplabel[5]=TUtil::GetMCFMParticleLabel(pId[pZOrder[3]], false, useQQVVQQany);
    }

  } // End check WW, ZZ, ZG etc.

  int* ordering;
  if (ndau<=3){
    if (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB){
      if (process == TVar::bkgZGamma){ ordering = pZOrder; if (!hasZ1 || !isZG) result = false; }
    }
  }
  else{

    if (
      (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB)
      ||
      ((production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU) && process == TVar::bkgZZ)
      ){
      // Just get the default ordering
      if (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB){
        if (process == TVar::bkgZGamma){ ordering = pZOrder; if (!hasZ1) result = false; }
        else if (process == TVar::bkgZZ){ ordering = pZOrder; if (!hasZ1 || !hasZ2) result = false; }
        else if (process == TVar::bkgWW){ ordering = pWOrder; if (!hasW1 || !hasW2) result = false; }
      }
      // bkgZZ process with S/T+U/STU separation needs ZZ ordering
      else{ ordering = pZOrder; if (!hasZ1 || !hasZ2) result = false; }
    }
    else if (production == TVar::ZZGG){
      // ggZZ/WW/VV
      // All of these require ZZ ordering since they use either ZZ ME or VV ME
      if (process==TVar::HSMHiggs || process==TVar::bkgZZ_SMHiggs || process==TVar::bkgZZ){
        ordering = pZOrder;
        if (
          !(
          (hasZ1 && hasZ2)
          ||
          (process==TVar::HSMHiggs && hasW1 && hasW2)
          )
          ) result = false;
      }
      else if (process==TVar::HSMHiggs_WWZZ || process==TVar::bkgWWZZ_SMHiggs || process==TVar::bkgWWZZ){
        ordering = pZOrder;
        if (!hasZ1 || !hasZ2 || !hasW1 || !hasW2) result = false;
      }
      else if (process==TVar::bkgWW_SMHiggs || process==TVar::bkgWW){
        ordering = pZOrder;
        if (!hasW1 || !hasW2) result = false;
      }
    }
    else if (
      production==TVar::Had_WH || production==TVar::Had_ZH
      || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
      || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
      || production==TVar::Lep_WH || production==TVar::Lep_ZH
      || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
      || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
      || production==TVar::JJVBF || production==TVar::JJEW || production==TVar::JJEWQCD
      || production==TVar::JJVBF_S || production==TVar::JJEW_S || production==TVar::JJEWQCD_S
      || production==TVar::JJVBF_TU || production==TVar::JJEW_TU || production==TVar::JJEWQCD_TU
      || production==TVar::JJQCD || production==TVar::JJQCD_S || production==TVar::JJQCD_TU
      ){
      // Z+2 jets
      // Just get the default ordering
      if (process == TVar::bkgZJets){ ordering = pZOrder; if (!hasZ1 || !isZJJ) result = false; }
      // VBF or QCD MCFM SBI, S or B
      // All of these require ZZ ordering since they use either ZZ ME or VV ME
      else if (process==TVar::HSMHiggs || process==TVar::bkgZZ_SMHiggs || process==TVar::bkgZZ){
        ordering = pZOrder;
        if (
          !(
          (hasZ1 && hasZ2)
          ||
          (process==TVar::HSMHiggs && hasW1 && hasW2)
          )
          ) result = false;
      }
      else if (process==TVar::HSMHiggs_WWZZ || process==TVar::bkgWWZZ_SMHiggs || process==TVar::bkgWWZZ){
        ordering = pZOrder;
        if (!hasZ1 || !hasZ2 || !hasW1 || !hasW2) result = false;
      }
      else if (process==TVar::bkgWW_SMHiggs || process==TVar::bkgWW){
        ordering = pZOrder;
        if (!hasW1 || !hasW2) result = false;
      }
    }

  }

  if (partOrder!=0) { for (unsigned int ip=0; ip<ndau; ip++) partOrder->push_back(ordering[ip]); }
  if (apartOrder!=0 && pApartOrder!=0) { for (unsigned int ip=0; ip<napart; ip++) apartOrder->push_back(pApartOrder[ip]); }
  for (int ip=0; ip<mxpart; ip++) sprintf((plabel_.plabel)[ip], strplabel[ip].Data());

  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::MCFM_SetupParticleCouplings: Summary (result=" << (int)result << "):\n";

    MELAout << "\trunstring=" << runstring_.runstring << endl;

    if (hasZ1) MELAout << "\tProcess found a Z1." << endl;
    if (hasZ2) MELAout << "\tProcess found a Z2." << endl;
    if (hasW1) MELAout << "\tProcess found a W1." << endl;
    if (hasW2) MELAout << "\tProcess found a W2." << endl;
    MELAout << "\t(l1, l2) = (" << zcouple_.l1 << ", " << zcouple_.l2 << ")" << endl;
    MELAout << "\t(r1, r2) = (" << zcouple_.r1 << ", " << zcouple_.r2 << ")" << endl;
    MELAout << "\t(q1, q2) = (" << zcouple_.q1 << ", " << zcouple_.q2 << ")" << endl;

    MELAout << "\tplabels:\n";
    for (int ip=0; ip<mxpart; ip++) MELAout << "\t[" << ip << "]=" << (plabel_.plabel)[ip] << endl;
    if (partOrder!=0){
      if (partOrder->size()>0) MELAout << "\tpartOrder[" << partOrder->size() << "]:\n";
      for (unsigned int ip=0; ip<partOrder->size(); ip++) MELAout << "\t[" << ip << "] -> " << partOrder->at(ip) << endl;
    }
    if (apartOrder!=0){
      if (apartOrder->size()>0) MELAout << "\tapartOrder[" << apartOrder->size() << "]:\n";
      for (unsigned int ip=0; ip<apartOrder->size(); ip++) MELAout << "\t[" << ip << "] -> " << apartOrder->at(ip) << endl;
    }
  }

  delete[] pId;
  delete[] pWOrder;
  delete[] pZOrder;
  delete[] pOrder;
  if (pApartOrder!=0) delete[] pApartOrder;
  if (pApartId!=0) delete[] pApartId;
  return result;
}
TString TUtil::GetMCFMParticleLabel(const int& pid, bool useQJ, bool useExtendedConventions){
  if (PDGHelpers::isAnUnknownJet(pid)){
    if (useQJ) return TString("qj");
    else return TString("pp");
  }
  else if (PDGHelpers::isAGluon(pid)) return TString("ig");
  else if (pid==1) return TString("dq");
  else if (pid==2) return TString("uq");
  else if (pid==3) return TString("sq");
  else if (pid==4) return TString("cq");
  else if (pid==5) return TString("bq");
  else if (pid==-1) return TString("da");
  else if (pid==-2) return TString("ua");
  else if (pid==-3) return TString("sa");
  else if (pid==-4) return TString("ca");
  else if (pid==-5) return TString("ba");
  else if (abs(pid)==6) return TString("qj"); // No tops!
  else if (pid==11) return TString("el");
  else if (pid==13) return TString("ml");
  else if (pid==15) return TString("tl");
  else if (pid==-11) return TString("ea");
  else if (pid==-13) return TString("ma");
  else if (pid==-15) return TString("ta");
  else if (std::abs(pid)>=12 && std::abs(pid)<=16){
    if (!useExtendedConventions){
      if (pid>0) return TString("nl");
      else return TString("na");
    }
    else{
      if (pid==12) return TString("ne");
      else if (pid==14) return TString("nm");
      else if (pid==16) return TString("nt");
      else if (pid==-12) return TString("ke");
      else if (pid==-14) return TString("km");
      else/* if (pid==-16)*/ return TString("kt");
    }
  }
  else return TString("  ");
}

bool TUtil::MCFM_masscuts(double s[][mxpart], const TVar::Process& process){
  double minZmassSqr=10*10;
  if (
    process==TVar::bkgZZ
    &&
    (s[2][3]<minZmassSqr || s[4][5]<minZmassSqr)
    ) return true;
  return false;
}
bool TUtil::MCFM_smalls(double s[][mxpart], int npart){

  // Reject event if any s(i,j) is too small
  // cutoff is defined in technical.Dat

  if (
    npart == 3 &&
    (
    (-s[5-1][1-1]< cutoff_.cutoff)  //gamma p1
    || (-s[5-1][2-1]< cutoff_.cutoff)  //gamma p2
    || (-s[4-1][1-1]< cutoff_.cutoff)  //e+    p1
    || (-s[4-1][2-1]< cutoff_.cutoff)  //e-    p2
    || (-s[3-1][1-1]< cutoff_.cutoff)  //nu    p1
    || (-s[3-1][2-1]< cutoff_.cutoff)  //nu    p2
    || (+s[5-1][4-1]< cutoff_.cutoff)  //gamma e+
    || (+s[5-1][3-1]< cutoff_.cutoff)  //gamma nu
    || (+s[4-1][3-1]< cutoff_.cutoff)  //e+    nu
    )
    )
    return true;

  else if (
    npart == 4 &&
    (
    (-s[5-1][1-1]< cutoff_.cutoff)  //e-    p1
    || (-s[5-1][2-1]< cutoff_.cutoff)  //e-    p2
    || (-s[6-1][1-1]< cutoff_.cutoff)  //nb    p1
    || (-s[6-1][2-1]< cutoff_.cutoff)  //nb    p2
    || (+s[6-1][5-1]< cutoff_.cutoff)  //e-    nb
    )

    )

    return true;

  return false;
}


void TUtil::InitJHUGenMELA(const char* pathtoPDFSet, int PDFMember){
  char path_pdf_c[200];
  sprintf(path_pdf_c, "%s", pathtoPDFSet);
  int pathpdfLength = strlen(path_pdf_c);
  __modjhugen_MOD_initfirsttime(path_pdf_c, &pathpdfLength, &PDFMember);
}
void TUtil::SetJHUGenHiggsMassWidth(double MReso, double GaReso){
  const double GeV = 1./100.;
  MReso *= GeV; // GeV units in JHUGen
  GaReso *= GeV; // GeV units in JHUGen
  __modjhugenmela_MOD_sethiggsmasswidth(&MReso, &GaReso);
}
void TUtil::SetJHUGenDistinguishWWCouplings(bool doAllow){
  int iAllow = (doAllow ? 1 : 0);
  __modjhugenmela_MOD_setdistinguishwwcouplingsflag(&iAllow);
}
void TUtil::ResetAmplitudeIncludes(){
  __modjhugenmela_MOD_resetamplitudeincludes();
}
void TUtil::SetMCFMSpinZeroCouplings(bool useBSM, SpinZeroCouplings* Hcouplings, bool forceZZ){
  if (!useBSM){
    spinzerohiggs_anomcoupl_.AllowAnomalousCouplings = 0;
    spinzerohiggs_anomcoupl_.distinguish_HWWcouplings = 0;

    /***** REGULAR RESONANCE *****/
    //
    spinzerohiggs_anomcoupl_.cz_q1sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_z11 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z21 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z31 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z41 = 100;
    spinzerohiggs_anomcoupl_.cz_q2sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_z12 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z22 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z32 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z42 = 100;
    spinzerohiggs_anomcoupl_.cz_q12sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_z10 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z20 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z30 = 100;
    spinzerohiggs_anomcoupl_.Lambda_z40 = 100;
    //
    spinzerohiggs_anomcoupl_.cw_q1sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_w11 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w21 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w31 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w41 = 100;
    spinzerohiggs_anomcoupl_.cw_q2sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_w12 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w22 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w32 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w42 = 100;
    spinzerohiggs_anomcoupl_.cw_q12sq = 0;
    spinzerohiggs_anomcoupl_.Lambda_w10 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w20 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w30 = 100;
    spinzerohiggs_anomcoupl_.Lambda_w40 = 100;
    //
    // Just in case, not really needed to set to the default value
    spinzerohiggs_anomcoupl_.kappa_top[0] = 1; spinzerohiggs_anomcoupl_.kappa_top[1] = 0;
    spinzerohiggs_anomcoupl_.kappa_bot[0] = 1; spinzerohiggs_anomcoupl_.kappa_bot[1] = 0;
    spinzerohiggs_anomcoupl_.ghz1[0] = 1; spinzerohiggs_anomcoupl_.ghz1[1] = 0;
    spinzerohiggs_anomcoupl_.ghw1[0] = 1; spinzerohiggs_anomcoupl_.ghw1[1] = 0;
    for (int im=0; im<2; im++){
      spinzerohiggs_anomcoupl_.kappa_tilde_top[im] = 0;
      spinzerohiggs_anomcoupl_.kappa_tilde_bot[im] = 0;
      spinzerohiggs_anomcoupl_.ghg2[im] = 0;
      spinzerohiggs_anomcoupl_.ghg3[im] = 0;
      spinzerohiggs_anomcoupl_.ghg4[im] = 0;
      spinzerohiggs_anomcoupl_.kappa_4gen_top[im] = 0;
      spinzerohiggs_anomcoupl_.kappa_4gen_bot[im] = 0;
      spinzerohiggs_anomcoupl_.kappa_tilde_4gen_top[im] = 0;
      spinzerohiggs_anomcoupl_.kappa_tilde_4gen_bot[im] = 0;
      spinzerohiggs_anomcoupl_.ghg2_4gen[im] = 0;
      spinzerohiggs_anomcoupl_.ghg3_4gen[im] = 0;
      spinzerohiggs_anomcoupl_.ghg4_4gen[im] = 0;
      //
      spinzerohiggs_anomcoupl_.ghz2[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghz1_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghz2_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghz3_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghz4_prime7[im] = 0;
      //
      spinzerohiggs_anomcoupl_.ghzgs1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghzgs2[im] = 0;
      spinzerohiggs_anomcoupl_.ghzgs3[im] = 0;
      spinzerohiggs_anomcoupl_.ghzgs4[im] = 0;
      spinzerohiggs_anomcoupl_.ghgsgs2[im] = 0;
      spinzerohiggs_anomcoupl_.ghgsgs3[im] = 0;
      spinzerohiggs_anomcoupl_.ghgsgs4[im] = 0;
      //
      spinzerohiggs_anomcoupl_.ghw2[im] =  0;
      spinzerohiggs_anomcoupl_.ghw3[im] =  0;
      spinzerohiggs_anomcoupl_.ghw4[im] =  0;
      spinzerohiggs_anomcoupl_.ghw1_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghw1_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghw2_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghw3_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.ghw4_prime7[im] = 0;

    }
    /***** END REGULAR RESONANCE *****/
    //
    /***** SECOND RESONANCE *****/
    //
    spinzerohiggs_anomcoupl_.c2z_q1sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_z11 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z21 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z31 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z41 = 100;
    spinzerohiggs_anomcoupl_.c2z_q2sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_z12 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z22 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z32 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z42 = 100;
    spinzerohiggs_anomcoupl_.c2z_q12sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_z10 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z20 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z30 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_z40 = 100;
    //
    spinzerohiggs_anomcoupl_.c2w_q1sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_w11 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w21 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w31 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w41 = 100;
    spinzerohiggs_anomcoupl_.c2w_q2sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_w12 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w22 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w32 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w42 = 100;
    spinzerohiggs_anomcoupl_.c2w_q12sq = 0;
    spinzerohiggs_anomcoupl_.Lambda2_w10 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w20 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w30 = 100;
    spinzerohiggs_anomcoupl_.Lambda2_w40 = 100;
    //
    // AllowAnomalousCouplings==0, so these are still treated as 1 when h2mass>=0
    spinzerohiggs_anomcoupl_.kappa2_top[0] = 0; spinzerohiggs_anomcoupl_.kappa2_top[1] = 0;
    spinzerohiggs_anomcoupl_.kappa2_bot[0] = 0; spinzerohiggs_anomcoupl_.kappa2_bot[1] = 0;
    spinzerohiggs_anomcoupl_.gh2z1[0] = 0; spinzerohiggs_anomcoupl_.gh2z1[1] = 0;
    spinzerohiggs_anomcoupl_.gh2w1[0] = 0; spinzerohiggs_anomcoupl_.gh2w1[1] = 0;
    for (int im=0; im<2; im++){
      spinzerohiggs_anomcoupl_.kappa2_tilde_top[im] = 0;
      spinzerohiggs_anomcoupl_.kappa2_tilde_bot[im] = 0;
      spinzerohiggs_anomcoupl_.gh2g2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2g3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2g4[im] = 0;
      spinzerohiggs_anomcoupl_.kappa2_4gen_top[im] = 0;
      spinzerohiggs_anomcoupl_.kappa2_4gen_bot[im] = 0;
      spinzerohiggs_anomcoupl_.kappa2_tilde_4gen_top[im] = 0;
      spinzerohiggs_anomcoupl_.kappa2_tilde_4gen_bot[im] = 0;
      spinzerohiggs_anomcoupl_.gh2g2_4gen[im] = 0;
      spinzerohiggs_anomcoupl_.gh2g3_4gen[im] = 0;
      spinzerohiggs_anomcoupl_.gh2g4_4gen[im] = 0;
      //
      spinzerohiggs_anomcoupl_.gh2z2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z1_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z2_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z3_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2z4_prime7[im] = 0;
      //
      spinzerohiggs_anomcoupl_.gh2zgs1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2zgs2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2zgs3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2zgs4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2gsgs2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2gsgs3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2gsgs4[im] = 0;
      //
      spinzerohiggs_anomcoupl_.gh2w2[im] =  0;
      spinzerohiggs_anomcoupl_.gh2w3[im] =  0;
      spinzerohiggs_anomcoupl_.gh2w4[im] =  0;
      spinzerohiggs_anomcoupl_.gh2w1_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w1_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w2_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w3_prime7[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime2[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime3[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime4[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime5[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime6[im] = 0;
      spinzerohiggs_anomcoupl_.gh2w4_prime7[im] = 0;
    }
    /***** END SECOND RESONANCE *****/
  }
  else{
    spinzerohiggs_anomcoupl_.AllowAnomalousCouplings = 1;
    spinzerohiggs_anomcoupl_.distinguish_HWWcouplings = ((Hcouplings->separateWWZZcouplings && !forceZZ) ? 1 : 0);

    /***** REGULAR RESONANCE *****/
    //
    spinzerohiggs_anomcoupl_.cz_q1sq = (Hcouplings->HzzCLambda_qsq)[cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.Lambda_z11 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.Lambda_z21 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.Lambda_z31 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.Lambda_z41 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.cz_q2sq = (Hcouplings->HzzCLambda_qsq)[cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.Lambda_z12 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.Lambda_z22 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.Lambda_z32 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.Lambda_z42 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.cz_q12sq = (Hcouplings->HzzCLambda_qsq)[cLambdaHIGGS_VV_QSQ12];
    spinzerohiggs_anomcoupl_.Lambda_z10 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12];
    spinzerohiggs_anomcoupl_.Lambda_z20 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12];
    spinzerohiggs_anomcoupl_.Lambda_z30 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12];
    spinzerohiggs_anomcoupl_.Lambda_z40 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12];
    //
    for (int im=0; im<2; im++){
      spinzerohiggs_anomcoupl_.kappa_top[im] = (Hcouplings->Httcoupl)[gHIGGS_KAPPA][im];
      spinzerohiggs_anomcoupl_.kappa_bot[im] = (Hcouplings->Hbbcoupl)[gHIGGS_KAPPA][im];
      spinzerohiggs_anomcoupl_.kappa_tilde_top[im] = (Hcouplings->Httcoupl)[gHIGGS_KAPPA_TILDE][im];
      spinzerohiggs_anomcoupl_.kappa_tilde_bot[im] = (Hcouplings->Hbbcoupl)[gHIGGS_KAPPA_TILDE][im];
      spinzerohiggs_anomcoupl_.ghg2[im] = (Hcouplings->Hggcoupl)[gHIGGS_GG_2][im];
      spinzerohiggs_anomcoupl_.ghg3[im] = (Hcouplings->Hggcoupl)[gHIGGS_GG_3][im];
      spinzerohiggs_anomcoupl_.ghg4[im] = (Hcouplings->Hggcoupl)[gHIGGS_GG_4][im];
      spinzerohiggs_anomcoupl_.kappa_4gen_top[im] = (Hcouplings->Ht4t4coupl)[gHIGGS_KAPPA][im];
      spinzerohiggs_anomcoupl_.kappa_4gen_bot[im] = (Hcouplings->Hb4b4coupl)[gHIGGS_KAPPA][im];
      spinzerohiggs_anomcoupl_.kappa_tilde_4gen_top[im] = (Hcouplings->Ht4t4coupl)[gHIGGS_KAPPA_TILDE][im];
      spinzerohiggs_anomcoupl_.kappa_tilde_4gen_bot[im] = (Hcouplings->Hb4b4coupl)[gHIGGS_KAPPA_TILDE][im];
      spinzerohiggs_anomcoupl_.ghg2_4gen[im] = (Hcouplings->Hg4g4coupl)[gHIGGS_GG_2][im];
      spinzerohiggs_anomcoupl_.ghg3_4gen[im] = (Hcouplings->Hg4g4coupl)[gHIGGS_GG_3][im];
      spinzerohiggs_anomcoupl_.ghg4_4gen[im] = (Hcouplings->Hg4g4coupl)[gHIGGS_GG_4][im];
      //
      spinzerohiggs_anomcoupl_.ghz1[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1][im];
      spinzerohiggs_anomcoupl_.ghz2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2][im];
      spinzerohiggs_anomcoupl_.ghz3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3][im];
      spinzerohiggs_anomcoupl_.ghz4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4][im];
      spinzerohiggs_anomcoupl_.ghz1_prime[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME][im];
      spinzerohiggs_anomcoupl_.ghz1_prime2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME2][im];
      spinzerohiggs_anomcoupl_.ghz1_prime3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME3][im];
      spinzerohiggs_anomcoupl_.ghz1_prime4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME4][im];
      spinzerohiggs_anomcoupl_.ghz1_prime5[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME5][im];
      spinzerohiggs_anomcoupl_.ghz1_prime6[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME6][im];
      spinzerohiggs_anomcoupl_.ghz1_prime7[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME7][im];
      spinzerohiggs_anomcoupl_.ghz2_prime[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME][im];
      spinzerohiggs_anomcoupl_.ghz2_prime2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME2][im];
      spinzerohiggs_anomcoupl_.ghz2_prime3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME3][im];
      spinzerohiggs_anomcoupl_.ghz2_prime4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME4][im];
      spinzerohiggs_anomcoupl_.ghz2_prime5[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME5][im];
      spinzerohiggs_anomcoupl_.ghz2_prime6[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME6][im];
      spinzerohiggs_anomcoupl_.ghz2_prime7[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME7][im];
      spinzerohiggs_anomcoupl_.ghz3_prime[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME][im];
      spinzerohiggs_anomcoupl_.ghz3_prime2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME2][im];
      spinzerohiggs_anomcoupl_.ghz3_prime3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME3][im];
      spinzerohiggs_anomcoupl_.ghz3_prime4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME4][im];
      spinzerohiggs_anomcoupl_.ghz3_prime5[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME5][im];
      spinzerohiggs_anomcoupl_.ghz3_prime6[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME6][im];
      spinzerohiggs_anomcoupl_.ghz3_prime7[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME7][im];
      spinzerohiggs_anomcoupl_.ghz4_prime[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME][im];
      spinzerohiggs_anomcoupl_.ghz4_prime2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME2][im];
      spinzerohiggs_anomcoupl_.ghz4_prime3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME3][im];
      spinzerohiggs_anomcoupl_.ghz4_prime4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME4][im];
      spinzerohiggs_anomcoupl_.ghz4_prime5[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME5][im];
      spinzerohiggs_anomcoupl_.ghz4_prime6[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME6][im];
      spinzerohiggs_anomcoupl_.ghz4_prime7[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME7][im];
      //
      spinzerohiggs_anomcoupl_.ghzgs1_prime2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_ZA_1_PRIME2][im];
      spinzerohiggs_anomcoupl_.ghzgs2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_ZA_2][im];
      spinzerohiggs_anomcoupl_.ghzgs3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_ZA_3][im];
      spinzerohiggs_anomcoupl_.ghzgs4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_ZA_4][im];
      spinzerohiggs_anomcoupl_.ghgsgs2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_AA_2][im];
      spinzerohiggs_anomcoupl_.ghgsgs3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_AA_3][im];
      spinzerohiggs_anomcoupl_.ghgsgs4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_AA_4][im];
    }
    //
    if (spinzerohiggs_anomcoupl_.distinguish_HWWcouplings==1){
      //
      spinzerohiggs_anomcoupl_.cw_q1sq = (Hcouplings->HwwCLambda_qsq)[cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda_w11 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda_w21 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda_w31 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda_w41 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.cw_q2sq = (Hcouplings->HwwCLambda_qsq)[cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda_w12 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda_w22 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda_w32 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda_w42 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.cw_q12sq = (Hcouplings->HwwCLambda_qsq)[cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda_w10 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda_w20 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda_w30 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda_w40 = (Hcouplings->HwwLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12];
      //
      for (int im=0; im<2; im++){
        spinzerohiggs_anomcoupl_.ghw1[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_1][im];
        spinzerohiggs_anomcoupl_.ghw2[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_2][im];
        spinzerohiggs_anomcoupl_.ghw3[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_3][im];
        spinzerohiggs_anomcoupl_.ghw4[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_4][im];
        spinzerohiggs_anomcoupl_.ghw1_prime[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_1_PRIME][im];
        spinzerohiggs_anomcoupl_.ghw1_prime2[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_1_PRIME2][im];
        spinzerohiggs_anomcoupl_.ghw1_prime3[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_1_PRIME3][im];
        spinzerohiggs_anomcoupl_.ghw1_prime4[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_1_PRIME4][im];
        spinzerohiggs_anomcoupl_.ghw1_prime5[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_1_PRIME5][im];
        spinzerohiggs_anomcoupl_.ghw1_prime6[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_1_PRIME6][im];
        spinzerohiggs_anomcoupl_.ghw1_prime7[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_1_PRIME7][im];
        spinzerohiggs_anomcoupl_.ghw2_prime[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_2_PRIME][im];
        spinzerohiggs_anomcoupl_.ghw2_prime2[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_2_PRIME2][im];
        spinzerohiggs_anomcoupl_.ghw2_prime3[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_2_PRIME3][im];
        spinzerohiggs_anomcoupl_.ghw2_prime4[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_2_PRIME4][im];
        spinzerohiggs_anomcoupl_.ghw2_prime5[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_2_PRIME5][im];
        spinzerohiggs_anomcoupl_.ghw2_prime6[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_2_PRIME6][im];
        spinzerohiggs_anomcoupl_.ghw2_prime7[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_2_PRIME7][im];
        spinzerohiggs_anomcoupl_.ghw3_prime[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_3_PRIME][im];
        spinzerohiggs_anomcoupl_.ghw3_prime2[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_3_PRIME2][im];
        spinzerohiggs_anomcoupl_.ghw3_prime3[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_3_PRIME3][im];
        spinzerohiggs_anomcoupl_.ghw3_prime4[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_3_PRIME4][im];
        spinzerohiggs_anomcoupl_.ghw3_prime5[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_3_PRIME5][im];
        spinzerohiggs_anomcoupl_.ghw3_prime6[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_3_PRIME6][im];
        spinzerohiggs_anomcoupl_.ghw3_prime7[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_3_PRIME7][im];
        spinzerohiggs_anomcoupl_.ghw4_prime[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_4_PRIME][im];
        spinzerohiggs_anomcoupl_.ghw4_prime2[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_4_PRIME2][im];
        spinzerohiggs_anomcoupl_.ghw4_prime3[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_4_PRIME3][im];
        spinzerohiggs_anomcoupl_.ghw4_prime4[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_4_PRIME4][im];
        spinzerohiggs_anomcoupl_.ghw4_prime5[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_4_PRIME5][im];
        spinzerohiggs_anomcoupl_.ghw4_prime6[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_4_PRIME6][im];
        spinzerohiggs_anomcoupl_.ghw4_prime7[im] = (Hcouplings->Hwwcoupl)[gHIGGS_VV_4_PRIME7][im];
      }
    }
    else{
      //
      spinzerohiggs_anomcoupl_.cw_q1sq = (Hcouplings->HzzCLambda_qsq)[cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda_w11 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda_w21 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda_w31 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda_w41 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.cw_q2sq = (Hcouplings->HzzCLambda_qsq)[cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda_w12 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda_w22 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda_w32 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda_w42 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.cw_q12sq = (Hcouplings->HzzCLambda_qsq)[cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda_w10 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda_w20 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda_w30 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda_w40 = (Hcouplings->HzzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12];
      //
      for (int im=0; im<2; im++){
        spinzerohiggs_anomcoupl_.ghw1[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1][im];
        spinzerohiggs_anomcoupl_.ghw2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2][im];
        spinzerohiggs_anomcoupl_.ghw3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3][im];
        spinzerohiggs_anomcoupl_.ghw4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4][im];
        spinzerohiggs_anomcoupl_.ghw1_prime[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME][im];
        spinzerohiggs_anomcoupl_.ghw1_prime2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME2][im];
        spinzerohiggs_anomcoupl_.ghw1_prime3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME3][im];
        spinzerohiggs_anomcoupl_.ghw1_prime4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME4][im];
        spinzerohiggs_anomcoupl_.ghw1_prime5[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME5][im];
        spinzerohiggs_anomcoupl_.ghw1_prime6[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME6][im];
        spinzerohiggs_anomcoupl_.ghw1_prime7[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_1_PRIME7][im];
        spinzerohiggs_anomcoupl_.ghw2_prime[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME][im];
        spinzerohiggs_anomcoupl_.ghw2_prime2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME2][im];
        spinzerohiggs_anomcoupl_.ghw2_prime3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME3][im];
        spinzerohiggs_anomcoupl_.ghw2_prime4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME4][im];
        spinzerohiggs_anomcoupl_.ghw2_prime5[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME5][im];
        spinzerohiggs_anomcoupl_.ghw2_prime6[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME6][im];
        spinzerohiggs_anomcoupl_.ghw2_prime7[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_2_PRIME7][im];
        spinzerohiggs_anomcoupl_.ghw3_prime[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME][im];
        spinzerohiggs_anomcoupl_.ghw3_prime2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME2][im];
        spinzerohiggs_anomcoupl_.ghw3_prime3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME3][im];
        spinzerohiggs_anomcoupl_.ghw3_prime4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME4][im];
        spinzerohiggs_anomcoupl_.ghw3_prime5[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME5][im];
        spinzerohiggs_anomcoupl_.ghw3_prime6[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME6][im];
        spinzerohiggs_anomcoupl_.ghw3_prime7[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_3_PRIME7][im];
        spinzerohiggs_anomcoupl_.ghw4_prime[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME][im];
        spinzerohiggs_anomcoupl_.ghw4_prime2[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME2][im];
        spinzerohiggs_anomcoupl_.ghw4_prime3[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME3][im];
        spinzerohiggs_anomcoupl_.ghw4_prime4[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME4][im];
        spinzerohiggs_anomcoupl_.ghw4_prime5[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME5][im];
        spinzerohiggs_anomcoupl_.ghw4_prime6[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME6][im];
        spinzerohiggs_anomcoupl_.ghw4_prime7[im] = (Hcouplings->Hzzcoupl)[gHIGGS_VV_4_PRIME7][im];
      }
    }
    /***** END REGULAR RESONANCE *****/
    //
    /***** SECOND RESONANCE *****/
    //
    spinzerohiggs_anomcoupl_.c2z_q1sq = (Hcouplings->H2zzCLambda_qsq)[cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.Lambda2_z11 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.Lambda2_z21 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.Lambda2_z31 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.Lambda2_z41 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1];
    spinzerohiggs_anomcoupl_.c2z_q2sq = (Hcouplings->H2zzCLambda_qsq)[cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.Lambda2_z12 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.Lambda2_z22 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.Lambda2_z32 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.Lambda2_z42 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2];
    spinzerohiggs_anomcoupl_.c2z_q12sq = (Hcouplings->H2zzCLambda_qsq)[cLambdaHIGGS_VV_QSQ12];
    spinzerohiggs_anomcoupl_.Lambda2_z10 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12];
    spinzerohiggs_anomcoupl_.Lambda2_z20 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12];
    spinzerohiggs_anomcoupl_.Lambda2_z30 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12];
    spinzerohiggs_anomcoupl_.Lambda2_z40 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12];
    //
    for (int im=0; im<2; im++){
      spinzerohiggs_anomcoupl_.kappa2_top[im] = (Hcouplings->H2ttcoupl)[gHIGGS_KAPPA][im];
      spinzerohiggs_anomcoupl_.kappa2_bot[im] = (Hcouplings->H2bbcoupl)[gHIGGS_KAPPA][im];
      spinzerohiggs_anomcoupl_.kappa2_tilde_top[im] = (Hcouplings->H2ttcoupl)[gHIGGS_KAPPA_TILDE][im];
      spinzerohiggs_anomcoupl_.kappa2_tilde_bot[im] = (Hcouplings->H2bbcoupl)[gHIGGS_KAPPA_TILDE][im];
      spinzerohiggs_anomcoupl_.gh2g2[im] = (Hcouplings->H2ggcoupl)[gHIGGS_GG_2][im];
      spinzerohiggs_anomcoupl_.gh2g3[im] = (Hcouplings->H2ggcoupl)[gHIGGS_GG_3][im];
      spinzerohiggs_anomcoupl_.gh2g4[im] = (Hcouplings->H2ggcoupl)[gHIGGS_GG_4][im];
      spinzerohiggs_anomcoupl_.kappa2_4gen_top[im] = (Hcouplings->H2t4t4coupl)[gHIGGS_KAPPA][im];
      spinzerohiggs_anomcoupl_.kappa2_4gen_bot[im] = (Hcouplings->H2b4b4coupl)[gHIGGS_KAPPA][im];
      spinzerohiggs_anomcoupl_.kappa2_tilde_4gen_top[im] = (Hcouplings->H2t4t4coupl)[gHIGGS_KAPPA_TILDE][im];
      spinzerohiggs_anomcoupl_.kappa2_tilde_4gen_bot[im] = (Hcouplings->H2b4b4coupl)[gHIGGS_KAPPA_TILDE][im];
      spinzerohiggs_anomcoupl_.gh2g2_4gen[im] = (Hcouplings->H2g4g4coupl)[gHIGGS_GG_2][im];
      spinzerohiggs_anomcoupl_.gh2g3_4gen[im] = (Hcouplings->H2g4g4coupl)[gHIGGS_GG_3][im];
      spinzerohiggs_anomcoupl_.gh2g4_4gen[im] = (Hcouplings->H2g4g4coupl)[gHIGGS_GG_4][im];
      //
      spinzerohiggs_anomcoupl_.gh2z1[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1][im];
      spinzerohiggs_anomcoupl_.gh2z2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2][im];
      spinzerohiggs_anomcoupl_.gh2z3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3][im];
      spinzerohiggs_anomcoupl_.gh2z4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME2][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME3][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME4][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime5[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME5][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime6[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME6][im];
      spinzerohiggs_anomcoupl_.gh2z1_prime7[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME7][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME2][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME3][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME4][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime5[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME5][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime6[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME6][im];
      spinzerohiggs_anomcoupl_.gh2z2_prime7[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME7][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME2][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME3][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME4][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime5[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME5][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime6[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME6][im];
      spinzerohiggs_anomcoupl_.gh2z3_prime7[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME7][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME2][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME3][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME4][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime5[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME5][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime6[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME6][im];
      spinzerohiggs_anomcoupl_.gh2z4_prime7[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME7][im];
      //
      spinzerohiggs_anomcoupl_.gh2zgs1_prime2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_ZA_1_PRIME2][im];
      spinzerohiggs_anomcoupl_.gh2zgs2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_ZA_2][im];
      spinzerohiggs_anomcoupl_.gh2zgs3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_ZA_3][im];
      spinzerohiggs_anomcoupl_.gh2zgs4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_ZA_4][im];
      spinzerohiggs_anomcoupl_.gh2gsgs2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_AA_2][im];
      spinzerohiggs_anomcoupl_.gh2gsgs3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_AA_3][im];
      spinzerohiggs_anomcoupl_.gh2gsgs4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_AA_4][im];
    }
    //
    if (spinzerohiggs_anomcoupl_.distinguish_HWWcouplings==1){
      //
      spinzerohiggs_anomcoupl_.c2w_q1sq = (Hcouplings->H2wwCLambda_qsq)[cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda2_w11 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda2_w21 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda2_w31 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda2_w41 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.c2w_q2sq = (Hcouplings->H2wwCLambda_qsq)[cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda2_w12 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda2_w22 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda2_w32 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda2_w42 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.c2w_q12sq = (Hcouplings->H2wwCLambda_qsq)[cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda2_w10 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda2_w20 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda2_w30 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda2_w40 = (Hcouplings->H2wwLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12];
      //
      for (int im=0; im<2; im++){
        spinzerohiggs_anomcoupl_.gh2w1[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_1][im];
        spinzerohiggs_anomcoupl_.gh2w2[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_2][im];
        spinzerohiggs_anomcoupl_.gh2w3[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_3][im];
        spinzerohiggs_anomcoupl_.gh2w4[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_4][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_1_PRIME][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime2[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_1_PRIME2][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime3[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_1_PRIME3][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime4[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_1_PRIME4][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime5[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_1_PRIME5][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime6[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_1_PRIME6][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime7[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_1_PRIME7][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_2_PRIME][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime2[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_2_PRIME2][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime3[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_2_PRIME3][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime4[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_2_PRIME4][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime5[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_2_PRIME5][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime6[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_2_PRIME6][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime7[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_2_PRIME7][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_3_PRIME][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime2[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_3_PRIME2][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime3[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_3_PRIME3][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime4[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_3_PRIME4][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime5[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_3_PRIME5][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime6[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_3_PRIME6][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime7[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_3_PRIME7][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_4_PRIME][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime2[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_4_PRIME2][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime3[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_4_PRIME3][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime4[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_4_PRIME4][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime5[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_4_PRIME5][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime6[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_4_PRIME6][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime7[im] = (Hcouplings->H2wwcoupl)[gHIGGS_VV_4_PRIME7][im];
      }
    }
    else{
      //
      spinzerohiggs_anomcoupl_.c2w_q1sq = (Hcouplings->H2zzCLambda_qsq)[cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda2_w11 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda2_w21 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda2_w31 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.Lambda2_w41 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ1];
      spinzerohiggs_anomcoupl_.c2w_q2sq = (Hcouplings->H2zzCLambda_qsq)[cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda2_w12 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda2_w22 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda2_w32 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.Lambda2_w42 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ2];
      spinzerohiggs_anomcoupl_.c2w_q12sq = (Hcouplings->H2zzCLambda_qsq)[cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda2_w10 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_1][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda2_w20 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_2][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda2_w30 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_3][cLambdaHIGGS_VV_QSQ12];
      spinzerohiggs_anomcoupl_.Lambda2_w40 = (Hcouplings->H2zzLambda_qsq)[LambdaHIGGS_QSQ_VV_4][cLambdaHIGGS_VV_QSQ12];
      //
      for (int im=0; im<2; im++){
        spinzerohiggs_anomcoupl_.gh2w1[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1][im];
        spinzerohiggs_anomcoupl_.gh2w2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2][im];
        spinzerohiggs_anomcoupl_.gh2w3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3][im];
        spinzerohiggs_anomcoupl_.gh2w4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME2][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME3][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME4][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime5[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME5][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime6[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME6][im];
        spinzerohiggs_anomcoupl_.gh2w1_prime7[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_1_PRIME7][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME2][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME3][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME4][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime5[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME5][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime6[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME6][im];
        spinzerohiggs_anomcoupl_.gh2w2_prime7[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_2_PRIME7][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME2][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME3][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME4][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime5[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME5][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime6[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME6][im];
        spinzerohiggs_anomcoupl_.gh2w3_prime7[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_3_PRIME7][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime2[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME2][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime3[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME3][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime4[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME4][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime5[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME5][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime6[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME6][im];
        spinzerohiggs_anomcoupl_.gh2w4_prime7[im] = (Hcouplings->H2zzcoupl)[gHIGGS_VV_4_PRIME7][im];
      }
    }
    /***** END SECOND RESONANCE *****/
  }
}
void TUtil::SetJHUGenSpinZeroVVCouplings(double Hvvcoupl[SIZE_HVV][2], double Hvvpcoupl[SIZE_HVV][2], double Hvpvpcoupl[SIZE_HVV][2], int Hvvcoupl_cqsq[SIZE_HVV_CQSQ], double HvvLambda_qsq[SIZE_HVV_LAMBDAQSQ][SIZE_HVV_CQSQ], bool useWWcoupl){
  const double GeV = 1./100.;
  int iWWcoupl = (useWWcoupl ? 1 : 0);
  for (int c=0; c<SIZE_HVV_LAMBDAQSQ; c++){ for (int k=0; k<SIZE_HVV_CQSQ; k++) HvvLambda_qsq[c][k] *= GeV; } // GeV units in JHUGen
  __modjhugenmela_MOD_setspinzerovvcouplings(Hvvcoupl, Hvvpcoupl, Hvpvpcoupl, Hvvcoupl_cqsq, HvvLambda_qsq, &iWWcoupl);
}
void TUtil::SetJHUGenSpinZeroGGCouplings(double Hggcoupl[SIZE_HGG][2]){ __modjhugenmela_MOD_setspinzeroggcouplings(Hggcoupl); }
void TUtil::SetJHUGenSpinZeroQQCouplings(double Hqqcoupl[SIZE_HQQ][2]){ __modjhugenmela_MOD_setspinzeroqqcouplings(Hqqcoupl); }
void TUtil::SetJHUGenSpinOneCouplings(double Zqqcoupl[SIZE_ZQQ][2], double Zvvcoupl[SIZE_ZVV][2]){ __modjhugenmela_MOD_setspinonecouplings(Zqqcoupl, Zvvcoupl); }
void TUtil::SetJHUGenSpinTwoCouplings(double Gacoupl[SIZE_GGG][2], double Gvvcoupl[SIZE_GVV][2], double Gvvpcoupl[SIZE_GVV][2], double Gvpvpcoupl[SIZE_GVV][2], double qLeftRightcoupl[SIZE_GQQ][2]){
  __modjhugenmela_MOD_setspintwocouplings(Gacoupl, Gvvcoupl, Gvvpcoupl, Gvpvpcoupl, qLeftRightcoupl);
}
void TUtil::SetJHUGenVprimeContactCouplings(double Zpffcoupl[SIZE_Vpff][2], double Wpffcoupl[SIZE_Vpff][2]){
  __modjhugenmela_MOD_setvprimecontactcouplings(Zpffcoupl, Wpffcoupl);
}

//Make sure
// 1. tot Energy Sum < 2EBEAM
// 2. PartonEnergy Fraction minimum<x0,x1<1
// 3. number of final state particle is defined
//
double TUtil::SumMatrixElementPDF(
  const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement, const TVar::LeptonInterference& leptonInterf,
  event_scales_type* event_scales, MelaIO* RcdME,
  const double& EBEAM,
  TVar::VerbosityLevel verbosity
  ){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin SumMatrixElementPDF" << endl;
  double msqjk=0;

  int nRequested_AssociatedJets=0;
  int nRequested_AssociatedLeptons=0;
  int AssociationVCompatibility=0;
  int partIncCode=TVar::kNoAssociated; // Do not use associated particles in the pT=0 frame boost
  if (production==TVar::JQCD){ // Use asociated jets in the pT=0 frame boost
    partIncCode=TVar::kUseAssociated_Jets;
    nRequested_AssociatedJets = 1;
    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::SumMatrixElementPDF: Requesting " << nRequested_AssociatedJets << " jets" << endl;
  }
  else if (
    production==TVar::Had_WH || production==TVar::Had_ZH
    || production==TVar::JJVBF || production==TVar::JJEW || production==TVar::JJEWQCD
    || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
    || production==TVar::JJVBF_S || production==TVar::JJEW_S || production==TVar::JJEWQCD_S
    || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
    || production==TVar::JJVBF_TU || production==TVar::JJEW_TU || production==TVar::JJEWQCD_TU
    || ((process==TVar::bkgZZ || process==TVar::bkgWW || process==TVar::bkgWWZZ) && (production==TVar::JJQCD || production==TVar::JJQCD_S || production==TVar::JJQCD_TU))
    ){ // Use asociated jets in the pT=0 frame boost
    partIncCode=TVar::kUseAssociated_Jets;
    nRequested_AssociatedJets = 2;
    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::SumMatrixElementPDF: Requesting " << nRequested_AssociatedJets << " jets" << endl;
  }
  else if (
    production==TVar::Lep_ZH || production==TVar::Lep_WH
    || production==TVar::Lep_ZH_S || production==TVar::Lep_WH_S
    || production==TVar::Lep_ZH_TU || production==TVar::Lep_WH_TU
    ){ // Only use associated leptons(+)neutrinos
    partIncCode=TVar::kUseAssociated_Leptons;
    nRequested_AssociatedLeptons = 2;
  }
  if (
    production==TVar::Had_WH || production==TVar::Lep_WH
    || production==TVar::Had_WH_S || production==TVar::Lep_WH_S
    || production==TVar::Had_WH_TU || production==TVar::Lep_WH_TU
    ) AssociationVCompatibility=24;
  else if (
    production==TVar::Had_ZH || production==TVar::Lep_ZH
    || production==TVar::Had_ZH_S || production==TVar::Lep_ZH_S
    || production==TVar::Had_ZH_TU || production==TVar::Lep_ZH_TU
    ) AssociationVCompatibility=23;
  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.AssociationVCompatibility=AssociationVCompatibility;
  mela_event.nRequested_AssociatedJets=nRequested_AssociatedJets;
  mela_event.nRequested_AssociatedLeptons=nRequested_AssociatedLeptons;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event,
    verbosity
    );

  double xx[2]={ 0 };
  vector<int> partOrder;
  vector<int> apartOrder;
  bool doProceed =
    CheckPartonMomFraction(mela_event.pMothers.at(0).second, mela_event.pMothers.at(1).second, xx, EBEAM, verbosity) // Check momentum transfers
    &&
    TUtil::MCFM_chooser(process, production, leptonInterf, verbosity, mela_event); // Set some of the specifics of the process through this function
  if (doProceed) doProceed = TUtil::MCFM_SetupParticleCouplings(process, production, verbosity, mela_event, &partOrder, &apartOrder); // Set the specifics of the daughter or associated particle couplings through this function
  if (doProceed){
    if (partOrder.size()!=mela_event.pDaughters.size()){
      if (verbosity >= TVar::ERROR){
        MELAerr << "TUtil::SumMatrixElementPDF: Ordering size " << partOrder.size() << " and number of daughter particles " << mela_event.pDaughters.size() << " are not the same!" << endl;
        TUtil::PrintCandidateSummary(&mela_event);
      }
      doProceed=false;
    }
    if (apartOrder.size()!=mela_event.pAssociated.size()){
      if (verbosity >= TVar::ERROR){
        MELAerr << "TUtil::SumMatrixElementPDF: Ordering size " << apartOrder.size() << " and number of associated particles " << mela_event.pAssociated.size() << " are not the same!" << endl;
        TUtil::PrintCandidateSummary(&mela_event);
      }
      doProceed=false;
    }
  }
  if (doProceed){
    int NPart=npart_.npart+2; // +2 for mothers
    double p4[4][mxpart]={ { 0 } };
    double p4_tmp[4][mxpart]={ { 0 } };
    int id[mxpart]; for (int ipar=0; ipar<mxpart; ipar++) id[ipar]=-9000;
    double msq[nmsq][nmsq]={ { 0 } };
    double msq_tmp[nmsq][nmsq]={ { 0 } };
    int channeltoggle=0;

    TLorentzVector MomStore[mxpart];
    for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);

    //Convert TLorentzVector into 4xNPart Matrix
    //reverse sign of incident partons
    for (int ipar=0; ipar<2; ipar++){
      if (mela_event.pMothers.at(ipar).second.T()>0.){
        p4[0][ipar] = -mela_event.pMothers.at(ipar).second.X();
        p4[1][ipar] = -mela_event.pMothers.at(ipar).second.Y();
        p4[2][ipar] = -mela_event.pMothers.at(ipar).second.Z();
        p4[3][ipar] = -mela_event.pMothers.at(ipar).second.T();
        MomStore[ipar] = mela_event.pMothers.at(ipar).second;
        id[ipar] = mela_event.pMothers.at(ipar).first;
      }
      else{
        p4[0][ipar] = mela_event.pMothers.at(ipar).second.X();
        p4[1][ipar] = mela_event.pMothers.at(ipar).second.Y();
        p4[2][ipar] = mela_event.pMothers.at(ipar).second.Z();
        p4[3][ipar] = mela_event.pMothers.at(ipar).second.T();
        MomStore[ipar] = -mela_event.pMothers.at(ipar).second;
        id[ipar] = mela_event.pMothers.at(ipar).first;
      }
    }

    // Determine if the decay mode involves WW or ZZ, to be used for ZZ or WW-specific signal MEs
    //bool isZG = (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(1)));
    bool isWW = (PDGHelpers::isAWBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAWBoson(mela_event.intermediateVid.at(1)));
    bool isZZ = (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(1)));
    //initialize decayed particles
    for (int ix=0; ix<(int)partOrder.size(); ix++){
      int ipar = min((int)mxpart, min(NPart, 2))+ix;
      if (ipar>=mxpart) break;
      TLorentzVector* momTmp = &(mela_event.pDaughters.at(partOrder.at(ix)).second);
      p4[0][ipar] = momTmp->X();
      p4[1][ipar] = momTmp->Y();
      p4[2][ipar] = momTmp->Z();
      p4[3][ipar] = momTmp->T();
      MomStore[ipar]=*momTmp;
      id[ipar] = mela_event.pDaughters.at(partOrder.at(ix)).first;
    }
    for (int ix=0; ix<(int)apartOrder.size(); ix++){
      int ipar = min((int)mxpart, min(NPart, (int)(partOrder.size()+2)))+ix;
      if (ipar>=mxpart) break;
      TLorentzVector* momTmp = &(mela_event.pAssociated.at(apartOrder.at(ix)).second);
      p4[0][ipar] = momTmp->X();
      p4[1][ipar] = momTmp->Y();
      p4[2][ipar] = momTmp->Z();
      p4[3][ipar] = momTmp->T();
      MomStore[ipar]=*momTmp;
      id[ipar] = mela_event.pAssociated.at(apartOrder.at(ix)).first;
    }

    if (verbosity >= TVar::DEBUG){ for (int i=0; i<NPart; i++) MELAout << "p["<<i<<"] (Px, Py, Pz, E):\t" << p4[0][i] << '\t' << p4[1][i] << '\t' << p4[2][i] << '\t' << p4[3][i] << endl; }

    double defaultRenScale = scale_.scale;
    double defaultFacScale = facscale_.facscale;
    int defaultNloop = nlooprun_.nlooprun;
    int defaultNflav = nflav_.nflav;
    string defaultPdflabel = pdlabel_.pdlabel;
    double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
    double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
    SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l");
    double alphasVal, alphasmzVal;
    GetAlphaS(&alphasVal, &alphasmzVal);
    RcdME->setRenormalizationScale(renQ);
    RcdME->setFactorizationScale(facQ);
    RcdME->setAlphaS(alphasVal);
    RcdME->setAlphaSatMZ(alphasmzVal);
    RcdME->setHiggsMassWidth(masses_mcfm_.hmass, masses_mcfm_.hwidth, 0);
    RcdME->setHiggsMassWidth(spinzerohiggs_anomcoupl_.h2mass, spinzerohiggs_anomcoupl_.h2width, 1);
    if (verbosity>=TVar::DEBUG){
      MELAout
        << "TUtil::SumMatrixElementPDF: Set AlphaS:\n"
        << "\tBefore set, alphas scale: " << defaultRenScale << ", PDF scale: " << defaultFacScale << '\n'
        << "\trenQ: " << renQ << " ( x " << event_scales->ren_scale_factor << "), facQ: " << facQ << " ( x " << event_scales->fac_scale_factor << ")\n"
        << "\tAfter set, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
    }

    int nInstances=0;
    if (
      (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB)
      ||
      ((production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU) && process == TVar::bkgZZ)
      ){
      if (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB){
        if (process == TVar::bkgZGamma) qqb_zgam_(p4[0], msq[0]);
        else if (process == TVar::bkgZZ) qqb_zz_(p4[0], msq[0]);
        else if (process == TVar::bkgWW) qqb_ww_(p4[0], msq[0]);
      }
      else{
        if (production == TVar::ZZQQB_STU) channeltoggle=0;
        else if (production == TVar::ZZQQB_S) channeltoggle=1;
        else/* if (production == TVar::ZZQQB_TU)*/ channeltoggle=2;
        qqb_zz_stu_(p4[0], msq[0], &channeltoggle);
      }
      nInstances=WipeMEArray(process, production, id, msq, verbosity);

      // Sum over valid MEs without PDF weights
      // By far the dominant contribution is uub initial state.
      for (int iquark=-5; iquark<=5; iquark++){
        for (int jquark=-5; jquark<=5; jquark++){
          if (
            (PDGHelpers::isAnUnknownJet(id[0]) || (PDGHelpers::isAGluon(id[0]) && iquark==0) || iquark==id[0])
            &&
            (PDGHelpers::isAnUnknownJet(id[1]) || (PDGHelpers::isAGluon(id[1]) && jquark==0) || jquark==id[1])
            ){
            if (
              PDGHelpers::isAnUnknownJet(id[0])
              &&
              PDGHelpers::isAnUnknownJet(id[1])
              ) msqjk = msq[3][7]+msq[7][3]; // Use only uub initial state
            else if (
              PDGHelpers::isAnUnknownJet(id[0])
              ) msqjk = msq[jquark+5][-jquark+5];
            else if (
              PDGHelpers::isAnUnknownJet(id[1])
              ) msqjk = msq[-iquark+5][iquark+5];
            else msqjk += msq[jquark+5][iquark+5];
          }
        }
      }
      SumMEPDF(MomStore[0], MomStore[1], msq, RcdME, EBEAM, verbosity);
    }
    else if (production == TVar::ZZGG){
      if (isZZ){
        if (process==TVar::HSMHiggs) gg_hzz_tb_(p4[0], msq[0]); // |ggHZZ|**2
        else if (process==TVar::bkgZZ_SMHiggs) gg_zz_all_(p4[0], msq[0]); // |ggZZ + ggHZZ|**2
        else if (process==TVar::bkgZZ) gg_zz_(p4[0], &(msq[5][5])); // |ggZZ|**2
        // The following WWZZ processes need no swapping because ZZ->4f (same type) is ordered such that 3/5 (C++: 2/4) swap gives W+W-.
        else if (process==TVar::HSMHiggs_WWZZ) gg_hvv_tb_(p4[0], msq[0]); // |ggHZZ+WW|**2
        else if (process==TVar::bkgWWZZ_SMHiggs) gg_vv_all_(p4[0], msq[0]); // |ggZZ + ggHZZ (+WW)|**2
        else if (process==TVar::bkgWWZZ) gg_vv_(p4[0], &(msq[5][5])); // |ggZZ+WW|**2

        if (verbosity>=TVar::DEBUG){
          MELAout << "\tTUtil::SumMatrixElementPDF: ZZGG && ZZ/WW/ZZ+WW MEs using ZZ (runstring: " << runstring_.runstring << ")" << endl;
          for (int i=0; i<NPart; i++) MELAout << "\tp["<<i<<"] (Px, Py, Pz, E):\t" << p4[0][i] << '\t' << p4[1][i] << '\t' << p4[2][i] << '\t' << p4[3][i] << endl;
        }
      }
      else if (isWW){
        if (process==TVar::HSMHiggs || process==TVar::HSMHiggs_WWZZ) gg_hvv_tb_(p4[0], msq[0]); // |ggHZZ+WW|**2
        else if (process==TVar::bkgWW_SMHiggs || process==TVar::bkgWWZZ_SMHiggs) gg_vv_all_(p4[0], msq[0]); // |ggZZ + ggHZZ (+WW)|**2
        else if (process==TVar::bkgWW || process==TVar::bkgWWZZ) gg_vv_(p4[0], &(msq[5][5])); // |ggZZ+WW|**2
        if (verbosity>=TVar::DEBUG){
          MELAout << "\tTUtil::SumMatrixElementPDF: ZZGG && ZZ/WW/ZZ+WW MEs using WW (runstring: " << runstring_.runstring << ")" << endl;
          for (int i=0; i<NPart; i++) MELAout << "\tp["<<i<<"] (Px, Py, Pz, E):\t" << p4[0][i] << '\t' << p4[1][i] << '\t' << p4[2][i] << '\t' << p4[3][i] << endl;
        }
      }

      nInstances=WipeMEArray(process, production, id, msq, verbosity);
      if (
        (PDGHelpers::isAnUnknownJet(id[0]) || PDGHelpers::isAGluon(id[0]))
        &&
        (PDGHelpers::isAnUnknownJet(id[1]) || PDGHelpers::isAGluon(id[1]))
        ){
        msqjk = msq[5][5];
      }
      SumMEPDF(MomStore[0], MomStore[1], msq, RcdME, EBEAM, verbosity);
    }
    else if (
      production==TVar::Had_WH || production==TVar::Had_ZH
      || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
      || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
      || production==TVar::Lep_WH || production==TVar::Lep_ZH
      || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
      || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
      || production==TVar::JJVBF || production==TVar::JJEW || production==TVar::JJEWQCD
      || production==TVar::JJVBF_S || production==TVar::JJEW_S || production==TVar::JJEWQCD_S
      || production==TVar::JJVBF_TU || production==TVar::JJEW_TU || production==TVar::JJEWQCD_TU
      || production==TVar::JJQCD || production==TVar::JJQCD_S || production==TVar::JJQCD_TU
      ){
      // Z+2 jets
      if (process == TVar::bkgZJets) qqb_z2jet_(p4[0], msq[0]);
      // VBF or QCD MCFM SBI, S or B
      else if (isZZ && (process==TVar::bkgZZ_SMHiggs || process==TVar::HSMHiggs || process==TVar::bkgZZ)){
        if (
          production==TVar::Had_WH || production==TVar::Had_ZH
          || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
          || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
          || production==TVar::Lep_WH || production==TVar::Lep_ZH
          || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
          || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
          || production==TVar::JJVBF || production==TVar::JJEW
          || production==TVar::JJVBF_S || production==TVar::JJEW_S
          || production==TVar::JJVBF_TU || production==TVar::JJEW_TU
          ){
          /*
          if (process==TVar::bkgZZ){//
            qq_zzqq_bkg_(p4[0], msq[0]);
            if (PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) && PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])){
              for (unsigned int ix=0; ix<4; ix++){
                for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
                swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
              }
              qq_zzqq_bkg_(p4_tmp[0], msq_tmp[0]);
              for (int iquark=-5; iquark<=5; iquark++){ for (int jquark=-5; jquark<=5; jquark++){ msq[jquark+5][iquark+5] = (msq[jquark+5][iquark+5] + msq_tmp[jquark+5][iquark+5]); if (iquark==jquark) msq[jquark+5][iquark+5]*=0.5; } }
            }
          }
          else{//
            */
            qq_zzqq_(p4[0], msq[0]);
            if (PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) && PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])){
              for (unsigned int ix=0; ix<4; ix++){
                for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
                swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
              }
              qq_zzqq_(p4_tmp[0], msq_tmp[0]);
              for (int iquark=-5; iquark<=5; iquark++){ for (int jquark=-5; jquark<=5; jquark++){ msq[jquark+5][iquark+5] = (msq[jquark+5][iquark+5] + msq_tmp[jquark+5][iquark+5]); if (iquark==jquark) msq[jquark+5][iquark+5]*=0.5; } }
            }
            else if (PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) || PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])){
              for (unsigned int ix=0; ix<4; ix++){
                for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
                swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
              }
              TString strlabel[2] ={ (plabel_.plabel)[7], (plabel_.plabel)[6] }; for (unsigned int il=0; il<2; il++) strlabel[il].Resize(2);
              for (int ip=0; ip<mxpart; ip++){
                if (ip!=6 && ip!=7) sprintf((plabel_.plabel)[ip], (plabel_.plabel)[ip]);
                else sprintf((plabel_.plabel)[ip], strlabel[ip-6].Data());
              }
              qq_zzqq_(p4_tmp[0], msq_tmp[0]);
              if (verbosity>=TVar::DEBUG){
                MELAout << "TUtil::SumMatrixElementPDF: Adding missing contributions:\n";
                MELAout << "\tplabels:\n";
                for (int ip=0; ip<mxpart; ip++) MELAout << "\t[" << ip << "]=" << (plabel_.plabel)[ip] << endl;
                MELAout << "\tMEsq initial:" << endl;
                for (int iquark=-5; iquark<=5; iquark++){
                  for (int jquark=-5; jquark<=5; jquark++) MELAout << msq[jquark+5][iquark+5] << '\t';
                  MELAout << endl;
                }
                MELAout << "\tMEsq added:" << endl;
                for (int iquark=-5; iquark<=5; iquark++){
                  for (int jquark=-5; jquark<=5; jquark++) MELAout << msq_tmp[jquark+5][iquark+5] << '\t';
                  MELAout << endl;
                }
              }
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++){
                  if (msq[jquark+5][iquark+5]==0.) msq[jquark+5][iquark+5] = msq_tmp[jquark+5][iquark+5];
                }
              }
            }
          //}//
        }
        else if (
          production==TVar::JJQCD || production==TVar::JJQCD_S || production==TVar::JJQCD_TU
          ){
          qq_zzqqstrong_(p4[0], msq[0]);
          if (PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) && PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])){
            for (unsigned int ix=0; ix<4; ix++){
              for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
              swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
            }
            qq_zzqqstrong_(p4_tmp[0], msq_tmp[0]);
            for (int iquark=-5; iquark<=5; iquark++){ for (int jquark=-5; jquark<=5; jquark++){ msq[jquark+5][iquark+5] = (msq[jquark+5][iquark+5] + msq_tmp[jquark+5][iquark+5]); if (iquark==jquark && iquark!=0) msq[jquark+5][iquark+5]*=0.5; } }
            // Subtract qqb/qbq->gg that was counted twice.
            TString gglabel = TUtil::GetMCFMParticleLabel(21, false, true);
            for (int ip=0; ip<mxpart; ip++){
              if (ip!=6 && ip!=7) sprintf((plabel_.plabel)[ip], (plabel_.plabel)[ip]);
              else sprintf((plabel_.plabel)[ip], gglabel.Data());
            }
            qq_zzqqstrong_(p4_tmp[0], msq_tmp[0]);
            for (int iquark=-5; iquark<=5; iquark++){ int jquark=-iquark; msq[jquark+5][iquark+5] = (msq[jquark+5][iquark+5] - msq_tmp[jquark+5][iquark+5]); }
            if (verbosity>=TVar::DEBUG){
              MELAout << "TUtil::SumMatrixElementPDF: Subtracting qqb/qbq->gg double-counted contribution:\n";
              MELAout << "\tplabels:\n";
              for (int ip=0; ip<mxpart; ip++) MELAout << "\t[" << ip << "]=" << (plabel_.plabel)[ip] << endl;
              MELAout << "\tMEsq:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq_tmp[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
            }
          }
          else if ((PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) || PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])) && !(PDGHelpers::isAGluon(id[partOrder.size()+2]) || PDGHelpers::isAGluon(id[partOrder.size()+3]))){
            for (unsigned int ix=0; ix<4; ix++){
              for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
              swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
            }
            TString strlabel[2] ={ (plabel_.plabel)[7], (plabel_.plabel)[6] }; for (unsigned int il=0; il<2; il++) strlabel[il].Resize(2);
            for (int ip=0; ip<mxpart; ip++){
              if (ip!=6 && ip!=7) sprintf((plabel_.plabel)[ip], (plabel_.plabel)[ip]);
              else sprintf((plabel_.plabel)[ip], strlabel[ip-6].Data());
            }
            qq_zzqqstrong_(p4_tmp[0], msq_tmp[0]);
            if (verbosity>=TVar::DEBUG){
              MELAout << "TUtil::SumMatrixElementPDF: Adding missing contributions:\n";
              MELAout << "\tplabels:\n";
              for (int ip=0; ip<mxpart; ip++) MELAout << "\t[" << ip << "]=" << (plabel_.plabel)[ip] << endl;
              MELAout << "\tMEsq initial:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
              MELAout << "\tMEsq added:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq_tmp[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
            }
            for (int iquark=-5; iquark<=5; iquark++){ 
              for (int jquark=-5; jquark<=5; jquark++){
                if (msq[jquark+5][iquark+5]==0.) msq[jquark+5][iquark+5] = msq_tmp[jquark+5][iquark+5];
              }
            }
          }
        }
      }
      else if (isWW && (process==TVar::bkgWW_SMHiggs || process==TVar::HSMHiggs || process==TVar::bkgWW)){
        if (
          production==TVar::Had_WH || production==TVar::Had_ZH
          || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
          || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
          || production==TVar::Lep_WH || production==TVar::Lep_ZH
          || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
          || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
          || production==TVar::JJVBF || production==TVar::JJEW
          || production==TVar::JJVBF_S || production==TVar::JJEW_S
          || production==TVar::JJVBF_TU || production==TVar::JJEW_TU
          ){
          qq_wwqq_(p4[0], msq[0]);
          if (PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) && PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])){
            for (unsigned int ix=0; ix<4; ix++){
              for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
              swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
            }
            qq_wwqq_(p4_tmp[0], msq_tmp[0]);
            for (int iquark=-5; iquark<=5; iquark++){ for (int jquark=-5; jquark<=5; jquark++){ msq[jquark+5][iquark+5] = (msq[jquark+5][iquark+5] + msq_tmp[jquark+5][iquark+5]); if (iquark==jquark) msq[jquark+5][iquark+5]*=0.5; } }
          }
          else if (PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) || PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])){
            for (unsigned int ix=0; ix<4; ix++){
              for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
              swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
            }
            TString strlabel[2] ={ (plabel_.plabel)[7], (plabel_.plabel)[6] }; for (unsigned int il=0; il<2; il++) strlabel[il].Resize(2);
            for (int ip=0; ip<mxpart; ip++){
              if (ip!=6 && ip!=7) sprintf((plabel_.plabel)[ip], (plabel_.plabel)[ip]);
              else sprintf((plabel_.plabel)[ip], strlabel[ip-6].Data());
            }
            qq_wwqq_(p4_tmp[0], msq_tmp[0]);
            if (verbosity>=TVar::DEBUG){
              MELAout << "TUtil::SumMatrixElementPDF: Adding missing contributions:\n";
              MELAout << "\tplabels:\n";
              for (int ip=0; ip<mxpart; ip++) MELAout << "\t[" << ip << "]=" << (plabel_.plabel)[ip] << endl;
              MELAout << "\tMEsq initial:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
              MELAout << "\tMEsq added:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq_tmp[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
            }
            for (int iquark=-5; iquark<=5; iquark++){
              for (int jquark=-5; jquark<=5; jquark++){
                if (msq[jquark+5][iquark+5]==0.) msq[jquark+5][iquark+5] = msq_tmp[jquark+5][iquark+5];
              }
            }
          }
        }
        else if (
          production==TVar::JJQCD || production==TVar::JJQCD_S || production==TVar::JJQCD_TU
          ){
          qq_wwqqstrong_(p4[0], msq[0]);
          if (PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) && PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])){
            for (unsigned int ix=0; ix<4; ix++){
              for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
              swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
            }
            qq_wwqqstrong_(p4_tmp[0], msq_tmp[0]);
            for (int iquark=-5; iquark<=5; iquark++){ for (int jquark=-5; jquark<=5; jquark++){ msq[jquark+5][iquark+5] = (msq[jquark+5][iquark+5] + msq_tmp[jquark+5][iquark+5]); if (iquark==jquark && iquark!=0) msq[jquark+5][iquark+5]*=0.5; } }
            // Subtract qqb/qbq->gg that was counted twice.
            TString gglabel = TUtil::GetMCFMParticleLabel(21, false, true);
            for (int ip=0; ip<mxpart; ip++){
              if (ip!=6 && ip!=7) sprintf((plabel_.plabel)[ip], (plabel_.plabel)[ip]);
              else sprintf((plabel_.plabel)[ip], gglabel.Data());
            }
            qq_wwqqstrong_(p4_tmp[0], msq_tmp[0]);
            for (int iquark=-5; iquark<=5; iquark++){ int jquark=-iquark; msq[jquark+5][iquark+5] = (msq[jquark+5][iquark+5] - msq_tmp[jquark+5][iquark+5]); }
            if (verbosity>=TVar::DEBUG){
              MELAout << "TUtil::SumMatrixElementPDF: Subtracting qqb/qbq->gg double-counted contribution:\n";
              MELAout << "\tplabels:\n";
              for (int ip=0; ip<mxpart; ip++) MELAout << "\t[" << ip << "]=" << (plabel_.plabel)[ip] << endl;
              MELAout << "\tMEsq:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq_tmp[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
            }
          }
          else if ((PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) || PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])) && !(PDGHelpers::isAGluon(id[partOrder.size()+2]) || PDGHelpers::isAGluon(id[partOrder.size()+3]))){
            for (unsigned int ix=0; ix<4; ix++){
              for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
              swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
            }
            TString strlabel[2] ={ (plabel_.plabel)[7], (plabel_.plabel)[6] }; for (unsigned int il=0; il<2; il++) strlabel[il].Resize(2);
            for (int ip=0; ip<mxpart; ip++){
              if (ip!=6 && ip!=7) sprintf((plabel_.plabel)[ip], (plabel_.plabel)[ip]);
              else sprintf((plabel_.plabel)[ip], strlabel[ip-6].Data());
            }
            qq_wwqqstrong_(p4_tmp[0], msq_tmp[0]);
            if (verbosity>=TVar::DEBUG){
              MELAout << "TUtil::SumMatrixElementPDF: Adding missing contributions:\n";
              MELAout << "\tplabels:\n";
              for (int ip=0; ip<mxpart; ip++) MELAout << "\t[" << ip << "]=" << (plabel_.plabel)[ip] << endl;
              MELAout << "\tMEsq initial:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
              MELAout << "\tMEsq added:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq_tmp[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
            }
            for (int iquark=-5; iquark<=5; iquark++){
              for (int jquark=-5; jquark<=5; jquark++){
                if (msq[jquark+5][iquark+5]==0.) msq[jquark+5][iquark+5] = msq_tmp[jquark+5][iquark+5];
              }
            }
          }
        }
      }
      else if ((isWW || isZZ) && (process==TVar::bkgWWZZ_SMHiggs || process==TVar::HSMHiggs_WWZZ || process==TVar::bkgWWZZ)){
        if (
          production==TVar::Had_WH || production==TVar::Had_ZH
          || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
          || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
          || production==TVar::Lep_WH || production==TVar::Lep_ZH
          || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
          || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
          || production==TVar::JJVBF || production==TVar::JJEW
          || production==TVar::JJVBF_S || production==TVar::JJEW_S
          || production==TVar::JJVBF_TU || production==TVar::JJEW_TU
          ){
          qq_vvqq_(p4[0], msq[0]);
          if (PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) && PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])){
            for (unsigned int ix=0; ix<4; ix++){
              for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
              swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
            }
            qq_vvqq_(p4_tmp[0], msq_tmp[0]);
            for (int iquark=-5; iquark<=5; iquark++){ for (int jquark=-5; jquark<=5; jquark++){ msq[jquark+5][iquark+5] = (msq[jquark+5][iquark+5] + msq_tmp[jquark+5][iquark+5]); if (iquark==jquark) msq[jquark+5][iquark+5]*=0.5; } }
          }
          else if (PDGHelpers::isAnUnknownJet(id[partOrder.size()+2]) || PDGHelpers::isAnUnknownJet(id[partOrder.size()+3])){
            for (unsigned int ix=0; ix<4; ix++){
              for (int ip=0; ip<mxpart; ip++) p4_tmp[ix][ip]=p4[ix][ip];
              swap(p4_tmp[ix][partOrder.size()+2], p4_tmp[ix][partOrder.size()+3]);
            }
            TString strlabel[2] ={ (plabel_.plabel)[7], (plabel_.plabel)[6] }; for (unsigned int il=0; il<2; il++) strlabel[il].Resize(2);
            for (int ip=0; ip<mxpart; ip++){
              if (ip!=6 && ip!=7) sprintf((plabel_.plabel)[ip], (plabel_.plabel)[ip]);
              else sprintf((plabel_.plabel)[ip], strlabel[ip-6].Data());
            }
            qq_vvqq_(p4_tmp[0], msq_tmp[0]);
            if (verbosity>=TVar::DEBUG){
              MELAout << "TUtil::SumMatrixElementPDF: Adding missing contributions:\n";
              MELAout << "\tplabels:\n";
              for (int ip=0; ip<mxpart; ip++) MELAout << "\t[" << ip << "]=" << (plabel_.plabel)[ip] << endl;
              MELAout << "\tMEsq initial:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
              MELAout << "\tMEsq added:" << endl;
              for (int iquark=-5; iquark<=5; iquark++){
                for (int jquark=-5; jquark<=5; jquark++) MELAout << msq_tmp[jquark+5][iquark+5] << '\t';
                MELAout << endl;
              }
            }
            for (int iquark=-5; iquark<=5; iquark++){
              for (int jquark=-5; jquark<=5; jquark++){
                if (msq[jquark+5][iquark+5]==0.) msq[jquark+5][iquark+5] = msq_tmp[jquark+5][iquark+5];
              }
            }
          }
        }
        //else if (
        //  production==TVar::JJQCD || production==TVar::JJQCD_S || production==TVar::JJQCD_TU
        //  ); // No JJQCD-VV ME in MCFM
      }

      nInstances=WipeMEArray(process, production, id, msq, verbosity);
      msqjk = SumMEPDF(MomStore[0], MomStore[1], msq, RcdME, EBEAM, verbosity);
    }

    // Set aL/R 1,2 into RcdME
    if (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0))) RcdME->setVDaughterCouplings(zcouple_.l1, zcouple_.r1, 0);
    else if (PDGHelpers::isAWBoson(mela_event.intermediateVid.at(0))) RcdME->setVDaughterCouplings(1., 0., 0);
    else RcdME->setVDaughterCouplings(0., 0., 0);
    if (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(1))) RcdME->setVDaughterCouplings(zcouple_.l2, zcouple_.r2, 1);
    else if (PDGHelpers::isAWBoson(mela_event.intermediateVid.at(1))) RcdME->setVDaughterCouplings(1., 0., 1);
    else RcdME->setVDaughterCouplings(0., 0., 1);
    // Reset some of these couplings if the decay is not actually ZZ or WW.
    if (
      process==TVar::bkgZJets
      ||
      process==TVar::bkgZGamma
      ) RcdME->setVDaughterCouplings(0., 0., 1);

    if (msqjk != msqjk){
      if (verbosity>=TVar::ERROR) MELAout << "TUtil::SumMatrixElementPDF: "<< TVar::ProcessName(process) << " msqjk="  << msqjk << endl;
      for (int i=0; i<NPart; i++) MELAout << "p["<<i<<"] (Px, Py, Pz, E):\t" << p4[0][i] << '\t' << p4[1][i] << '\t' << p4[2][i] << '\t' << p4[3][i] << endl;
      MELAout << "TUtil::SumMatrixElementPDF: The number of ME instances: " << nInstances << ". MEsq[ip][jp] = " << endl;
      for (int iquark=-5; iquark<=5; iquark++){
        for (int jquark=-5; jquark<=5; jquark++) MELAout << msq[jquark+5][iquark+5] << '\t';
        MELAout << endl;
      }
      msqjk=0;
    }
    else if (verbosity>=TVar::DEBUG){
      MELAout << "TUtil::SumMatrixElementPDF: The number of ME instances: " << nInstances << ". MEsq[ip][jp] = " << endl;
      for (int iquark=-5; iquark<=5; iquark++){
        for (int jquark=-5; jquark<=5; jquark++) MELAout << msq[jquark+5][iquark+5] << '\t';
        MELAout << endl;
      }
    }

    if (verbosity>=TVar::DEBUG) MELAout
        << "TUtil::SumMatrixElementPDF: Reset AlphaS:\n"
        << "\tBefore reset, alphas scale: " << scale_.scale
        << ", PDF scale: " << facscale_.facscale
        << endl;
    SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel);
    if (verbosity>=TVar::DEBUG) MELAout
        << "TUtil::SumMatrixElementPDF: Reset AlphaS result:\n"
        << "\tAfter reset, alphas scale: " << scale_.scale
        << ", PDF scale: " << facscale_.facscale
        << endl;
  } // End if doProceed
  if (verbosity>=TVar::DEBUG) MELAout << "End TUtil::SumMatrixElementPDF(" << msqjk << ")" << endl;
  return msqjk;
}

// ggHdec+0J or Hdec+0J
double TUtil::JHUGenMatEl(
  const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  const double& EBEAM,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double MatElSq=0; // Return value

  if (matrixElement!=TVar::JHUGen){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::JHUGenMatEl: Non-JHUGen MEs are not supported" << endl; return MatElSq; }
  bool isSpinZero = (
    process == TVar::HSMHiggs
    || process == TVar::H0minus
    || process == TVar::H0hplus
    || process == TVar::H0_g1prime2
    || process== TVar::H0_Zgs
    || process ==TVar::H0_gsgs
    || process ==TVar::H0_Zgs_PS
    || process ==TVar::H0_gsgs_PS
    || process ==TVar::H0_Zgsg1prime2
    || process == TVar::SelfDefine_spin0
    );
  bool isSpinOne = (
    process == TVar::H1minus
    || process == TVar::H1plus
    || process == TVar::SelfDefine_spin1
    );
  bool isSpinTwo = (
    process == TVar::H2_g1g5
    || process == TVar::H2_g1
    || process == TVar::H2_g8
    || process == TVar::H2_g4
    || process == TVar::H2_g5
    || process == TVar::H2_g2
    || process == TVar::H2_g3
    || process == TVar::H2_g6
    || process == TVar::H2_g7
    || process == TVar::H2_g9
    || process == TVar::H2_g10
    || process == TVar::SelfDefine_spin2
    );
  if (!(isSpinZero || isSpinOne || isSpinTwo)){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::JHUGenMatEl: Process " << TVar::ProcessName(process) << " (" << process << ")" << " not supported." << endl; return MatElSq; }
  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::JHUGenMatEl: Process " << TVar::ProcessName(process) << " is computing a spin-";
    if (isSpinZero) MELAout << "0";
    else if (isSpinOne) MELAout << "1";
    else if (isSpinTwo) MELAout << "2";
    else MELAout << "?";
    MELAout << " ME with production " << TVar::ProductionName(production) << "." << endl;
  }

  double msq[nmsq][nmsq]={ { 0 } }; // ME**2[parton2][parton1] for each incoming parton 1 and 2, used in RcdME
  int MYIDUP_tmp[4]={ 0 }; // Initial assignment array, unconverted. 0==Unassigned
  vector<pair<int, int>> idarray[2]; // All possible ids for each daughter based on the value of MYIDUP_tmp[0:3] and the desired V ids taken from mela_event.intermediateVid.at(0:1)
  int MYIDUP[4]={ 0 }; // "Working" assignment, converted
  int idfirst[2]={ 0 }; // Used to set DecayMode1, = MYIDUP[0:1]
  int idsecond[2]={ 0 }; // Used to set DecayMode2, = MYIDUP[2:3]
  double p4[6][4]={ { 0 } }; // Mom (*GeV) to pass into JHUGen
  TLorentzVector MomStore[mxpart]; // Mom (in natural units) to compute alphaS
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);

  // Notice that partIncCode is specific for this subroutine
  int partIncCode=TVar::kNoAssociated; // Do not use associated particles in the pT=0 frame boost
  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event,
    verbosity
    );
  if (mela_event.pDaughters.size()<2 || mela_event.intermediateVid.size()!=2){
    if (verbosity>=TVar::ERROR) MELAerr << "TUtil::JHUGenMatEl: Number of daughters " << mela_event.pDaughters.size() << " or number of intermediate Vs " << mela_event.intermediateVid << " not supported!" << endl;
    return MatElSq;
  }

  // p(i,0:3) = (E(i),px(i),py(i),pz(i))
  // i=0,1: g1,g2 or q1, qb2 (outgoing convention)
  // i=2,3: correspond to MY_IDUP(0),MY_IDUP(1)
  // i=4,5: correspond to MY_IDUP(2),MY_IDUP(3)
  for (int ipar=0; ipar<2; ipar++){
    if (mela_event.pMothers.at(ipar).second.T()>0.){
      p4[ipar][0] = -mela_event.pMothers.at(ipar).second.T()*GeV;
      p4[ipar][1] = -mela_event.pMothers.at(ipar).second.X()*GeV;
      p4[ipar][2] = -mela_event.pMothers.at(ipar).second.Y()*GeV;
      p4[ipar][3] = -mela_event.pMothers.at(ipar).second.Z()*GeV;
      MomStore[ipar] = mela_event.pMothers.at(ipar).second;
    }
    else{
      p4[ipar][0] = mela_event.pMothers.at(ipar).second.T()*GeV;
      p4[ipar][1] = mela_event.pMothers.at(ipar).second.X()*GeV;
      p4[ipar][2] = mela_event.pMothers.at(ipar).second.Y()*GeV;
      p4[ipar][3] = mela_event.pMothers.at(ipar).second.Z()*GeV;
      MomStore[ipar] = -mela_event.pMothers.at(ipar).second;
    }
    // From Markus:
    // Note that the momentum no.2, p(1:4, 2), is a dummy which is not used in case production==TVar::ZZINDEPENDENT.
    if (ipar==1 && production==TVar::ZZINDEPENDENT){ for (int ix=0; ix<4; ix++){ p4[0][ix] += p4[ipar][ix]; p4[ipar][ix]=0.; } }
  }
  //initialize decayed particles
  if (mela_event.pDaughters.size()==2){
    for (unsigned int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
      TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
      int* idtmp = &(mela_event.pDaughters.at(ipar).first);

      int arrindex = 2*ipar;
      if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_tmp[arrindex] = *idtmp;
      else MYIDUP_tmp[arrindex] = 0;
      MYIDUP_tmp[arrindex+1] = -9000;
      p4[arrindex+2][0] = momTmp->T()*GeV;
      p4[arrindex+2][1] = momTmp->X()*GeV;
      p4[arrindex+2][2] = momTmp->Y()*GeV;
      p4[arrindex+2][3] = momTmp->Z()*GeV;
      MomStore[arrindex+2] = *momTmp;
      if (verbosity >= TVar::DEBUG) MELAout << "MYIDUP_tmp[" << arrindex << "(" << ipar << ")" << "]=" << MYIDUP_tmp[arrindex] << endl;
    }
  }
  else{
    for (unsigned int ipar=0; ipar<4; ipar++){
      if (ipar<mela_event.pDaughters.size()){
        TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
        int* idtmp = &(mela_event.pDaughters.at(ipar).first);

        if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_tmp[ipar] = *idtmp;
        else MYIDUP_tmp[ipar] = 0;
        p4[ipar+2][0] = momTmp->T()*GeV;
        p4[ipar+2][1] = momTmp->X()*GeV;
        p4[ipar+2][2] = momTmp->Y()*GeV;
        p4[ipar+2][3] = momTmp->Z()*GeV;
        MomStore[ipar+2] = *momTmp;
      }
      else MYIDUP_tmp[ipar] = -9000; // No need to set p4, which is already 0 by initialization
      // __modparameters_MOD_not_a_particle__?
      if (verbosity >= TVar::DEBUG) MELAout << "MYIDUP_tmp[" << ipar << "]=" << MYIDUP_tmp[ipar] << endl;
    }
  }

  // Set alphas
  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l"); // Set AlphaS(|Q|/2, mynloop, mynflav, mypartonPDF)
  double alphasVal, alphasmzVal;
  GetAlphaS(&alphasVal, &alphasmzVal);
  RcdME->setRenormalizationScale(renQ);
  RcdME->setFactorizationScale(facQ);
  RcdME->setAlphaS(alphasVal);
  RcdME->setAlphaSatMZ(alphasmzVal);
  RcdME->setHiggsMassWidth(masses_mcfm_.hmass, masses_mcfm_.hwidth, 0);
  RcdME->setHiggsMassWidth(spinzerohiggs_anomcoupl_.h2mass, spinzerohiggs_anomcoupl_.h2width, 1);
  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::JHUGenMatEl: Set AlphaS:\n"
      << "\tBefore set, alphas scale: " << defaultRenScale << ", PDF scale: " << defaultFacScale << '\n'
      << "\trenQ: " << renQ << " ( x " << event_scales->ren_scale_factor << "), facQ: " << facQ << " ( x " << event_scales->fac_scale_factor << ")\n"
      << "\tAfter set, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }

  // Determine te actual ids to compute the ME. Assign ids if any are unknown.
  for (int iv=0; iv<2; iv++){ // Max. 2 vector bosons
    if (MYIDUP_tmp[2*iv+0]!=0 && MYIDUP_tmp[2*iv+1]!=0){ // Z->2l,2nu,2q, W->lnu,qq', G
      // OSSF pairs or just one "V-daughter"
      if (
        TMath::Sign(1, MYIDUP_tmp[2*iv+0])==TMath::Sign(1, -MYIDUP_tmp[2*iv+1])
        ||
        (MYIDUP_tmp[2*iv+0]==-9000 || MYIDUP_tmp[2*iv+1]==-9000)
        ) idarray[iv].push_back(pair<int, int>(MYIDUP_tmp[2*iv+0], MYIDUP_tmp[2*iv+1]));
      // SSSF pairs, ordered already by phi, so avoid the re-ordering inside the ME
      else if (MYIDUP_tmp[2*iv+0]<0) idarray[iv].push_back(pair<int, int>(-MYIDUP_tmp[2*iv+0], MYIDUP_tmp[2*iv+1])); // Reverse sign of first daughter if SS(--)SF pair
      else idarray[iv].push_back(pair<int, int>(MYIDUP_tmp[2*iv+0], -MYIDUP_tmp[2*iv+1])); // Reverse sign of daughter if SS(++)SF pair
    }
    else if (MYIDUP_tmp[2*iv+0]!=0){ // ->f?,  one quark is known
      if (PDGHelpers::isAWBoson(mela_event.intermediateVid.at(iv))){ // (W+/-)->f?
        if (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2*iv+0])){ // (W+)->u?
          int id_dn = TMath::Sign(1, -MYIDUP_tmp[2*iv+0]);
          int id_st = TMath::Sign(3, -MYIDUP_tmp[2*iv+0]);
          int id_bt = TMath::Sign(5, -MYIDUP_tmp[2*iv+0]);
          idarray[iv].push_back(pair<int, int>(MYIDUP_tmp[2*iv+0], id_dn));
          idarray[iv].push_back(pair<int, int>(MYIDUP_tmp[2*iv+0], id_st));
          idarray[iv].push_back(pair<int, int>(MYIDUP_tmp[2*iv+0], id_bt));
        }
        else if (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2*iv+0])){ // (W-)->d?
          int id_up = TMath::Sign(2, -MYIDUP_tmp[2*iv+0]);
          int id_ch = TMath::Sign(4, -MYIDUP_tmp[2*iv+0]);
          idarray[iv].push_back(pair<int, int>(MYIDUP_tmp[2*iv+0], id_up));
          idarray[iv].push_back(pair<int, int>(MYIDUP_tmp[2*iv+0], id_ch));
        }
      }
      else if (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(iv))){ // Z->f?
        idarray[iv].push_back(pair<int, int>(MYIDUP_tmp[2*iv+0], -MYIDUP_tmp[2*iv+0]));
      }
    }
    else if (MYIDUP_tmp[2*iv+1]!=0){ // ->?fb,  one quark is known
      if (PDGHelpers::isAWBoson(mela_event.intermediateVid.at(iv))){ // (W+/-)->?fb
        if (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2*iv+1])){ // (W-)->?ub
          int id_dn = TMath::Sign(1, -MYIDUP_tmp[2*iv+1]);
          int id_st = TMath::Sign(3, -MYIDUP_tmp[2*iv+1]);
          int id_bt = TMath::Sign(5, -MYIDUP_tmp[2*iv+1]);
          idarray[iv].push_back(pair<int, int>(id_dn, MYIDUP_tmp[2*iv+1]));
          idarray[iv].push_back(pair<int, int>(id_st, MYIDUP_tmp[2*iv+1]));
          idarray[iv].push_back(pair<int, int>(id_bt, MYIDUP_tmp[2*iv+1]));
        }
        else if (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2*iv+1])){ // (W+)->?db
          int id_up = TMath::Sign(2, -MYIDUP_tmp[2*iv+1]);
          int id_ch = TMath::Sign(4, -MYIDUP_tmp[2*iv+1]);
          idarray[iv].push_back(pair<int, int>(id_up, MYIDUP_tmp[2*iv+1]));
          idarray[iv].push_back(pair<int, int>(id_ch, MYIDUP_tmp[2*iv+1]));
        }
      }
      else if (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(iv))){ // Z->?fb
        idarray[iv].push_back(pair<int, int>(-MYIDUP_tmp[2*iv+1], MYIDUP_tmp[2*iv+1]));
      }
    }
    else{ // Both fermions unknown
      if (PDGHelpers::isAWBoson(mela_event.intermediateVid.at(iv))){ // (W+/-)->??
        for (int iqu=1; iqu<=2; iqu++){
          int id_up = iqu*2;
          for (int iqd=1; iqd<=3; iqd++){
            int id_dn = iqd*2-1;
            if (iv==0){ idarray[iv].push_back(pair<int, int>(id_up, -id_dn)); idarray[iv].push_back(pair<int, int>(-id_dn, id_up)); }
            else{ idarray[iv].push_back(pair<int, int>(id_dn, -id_up)); idarray[iv].push_back(pair<int, int>(-id_up, id_dn)); }
          }
        }
      }
      else if (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(iv))){ // Z->??
        for (int iq=1; iq<=5; iq++){ idarray[iv].push_back(pair<int, int>(iq, -iq)); idarray[iv].push_back(pair<int, int>(-iq, iq)); }
      }
    }
  }

  // Variables to set L1/2 and R1/2 couplings into RcdME
  const int InvalidMode=-1, WWMode=0, ZZMode=1, ZgMode=5, ggMode=7;
  int VVMode=InvalidMode;
  if (PDGHelpers::isAWBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAWBoson(mela_event.intermediateVid.at(1))) VVMode=WWMode;
  else if (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAZBoson(mela_event.intermediateVid.at(1))) VVMode=ZZMode;
  else if (PDGHelpers::isAPhoton(mela_event.intermediateVid.at(0)) && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(1))) VVMode=ggMode;
  else if (
    (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(0)) && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(1)))
    ||
    (PDGHelpers::isAZBoson(mela_event.intermediateVid.at(1)) && PDGHelpers::isAPhoton(mela_event.intermediateVid.at(0)))
    ) VVMode=ZgMode;
  double aL1=0, aL2=0, aR1=0, aR2=0;

  int nNonZero=0;
  for (unsigned int v1=0; v1<idarray[0].size(); v1++){
    for (unsigned int v2=0; v2<idarray[1].size(); v2++){
      // Convert the particle ids to JHU convention
      if (idarray[0].at(v1).first!=-9000) MYIDUP[0] = convertLHEreverse(&(idarray[0].at(v1).first));
      if (idarray[0].at(v1).second!=-9000) MYIDUP[1] = convertLHEreverse(&(idarray[0].at(v1).second));
      if (idarray[1].at(v2).first!=-9000) MYIDUP[2] = convertLHEreverse(&(idarray[1].at(v2).first));
      if (idarray[1].at(v2).second!=-9000) MYIDUP[3] = convertLHEreverse(&(idarray[1].at(v2).second));

      // Check working ids
      if (verbosity>=TVar::DEBUG){ for (unsigned int idau=0; idau<4; idau++) MELAout << "MYIDUP[" << idau << "]=" << MYIDUP[idau] << endl; }

      // Determine M_V and Ga_V in JHUGen, needed for g1 vs everything else.
      for (int ip=0; ip<2; ip++){ idfirst[ip]=MYIDUP[ip]; idsecond[ip]=MYIDUP[ip+2]; }
      __modjhugenmela_MOD_setdecaymodes(idfirst, idsecond); // Set M_V and Ga_V in JHUGen
      if (verbosity>=TVar::DEBUG){
        double mv, gv;
        __modjhugenmela_MOD_getmvgv(&mv, &gv);
        MELAout << "TUtil::JHUGenMatEl: M_V=" << mv/GeV << ", Ga_V=" << gv/GeV << endl;
        __modjhugenmela_MOD_getmvprimegvprime(&mv, &gv);
        MELAout << "TUtil::JHUGenMatEl: M_Vprime=" << mv/GeV << ", Ga_Vprime=" << gv/GeV << endl;
      }

      // Sum over possible left/right couplings of the Vs
      double aLRtmp[4]={ 0 };
      __modjhugenmela_MOD_getdecaycouplings(&VVMode, MYIDUP, &(aLRtmp[0]), &(aLRtmp[1]), &(aLRtmp[2]), &(aLRtmp[3]));
      if (idarray[0].size()>1){
        aL1 = sqrt(pow(aL1, 2)+pow(aLRtmp[0], 2));
        aR1 = sqrt(pow(aR1, 2)+pow(aLRtmp[1], 2));
      }
      else{
        aL1 = aLRtmp[0];
        aR1 = aLRtmp[1];
      }
      if (idarray[1].size()>1){
        aL2 = sqrt(pow(aL2, 2)+pow(aLRtmp[2], 2));
        aR2 = sqrt(pow(aR2, 2)+pow(aLRtmp[3], 2));
      }
      else{
        aL2 = aLRtmp[2];
        aR2 = aLRtmp[3];
      }

      double MatElTmp=0.;
      if (production == TVar::ZZGG){
        if (isSpinZero) __modhiggs_MOD_evalamp_gg_h_vv(p4, MYIDUP, &MatElTmp);
        else if (isSpinTwo) __modgraviton_MOD_evalamp_gg_g_vv(p4, MYIDUP, &MatElTmp);
      }
      else if (production == TVar::ZZQQB){
        if (isSpinOne) __modzprime_MOD_evalamp_qqb_zprime_vv(p4, MYIDUP, &MatElTmp);
        else if (isSpinTwo) __modgraviton_MOD_evalamp_qqb_g_vv(p4, MYIDUP, &MatElTmp);
      }
      else if (production == TVar::ZZINDEPENDENT){
        if (isSpinZero) __modhiggs_MOD_evalamp_h_vv(p4, MYIDUP, &MatElTmp);
        else if (isSpinOne) __modzprime_MOD_evalamp_zprime_vv(p4, MYIDUP, &MatElTmp);
        else if (isSpinTwo) __modgraviton_MOD_evalamp_g_vv(p4, MYIDUP, &MatElTmp);
      }
      // Add CKM elements since they are not included
      if (PDGHelpers::isAWBoson(mela_event.intermediateVid.at(0))) MatElTmp *= pow(__modparameters_MOD_ckmbare(&(idarray[0].at(v1).first), &(idarray[0].at(v1).second)), 2);
      if (PDGHelpers::isAWBoson(mela_event.intermediateVid.at(1))) MatElTmp *= pow(__modparameters_MOD_ckmbare(&(idarray[1].at(v2).first), &(idarray[1].at(v2).second)), 2);

      if (verbosity >= TVar::DEBUG) MELAout << "=====\nTUtil::JHUGenMatEl: Instance MatElTmp = " << MatElTmp << "\n=====" << endl;
      MatElSq += MatElTmp;
      if (MatElTmp>0.) nNonZero++;
    }
  }
  if (nNonZero>0) MatElSq /= ((double)nNonZero);
  if (verbosity >= TVar::DEBUG){
    MELAout << "TUtil::JHUGenMatEl: Number of matrix element instances computed: " << nNonZero << endl;
    MELAout << "TUtil::JHUGenMatEl: MatElSq after division = " << MatElSq << endl;
  }

  // Set aL/R 1,2 into RcdME
  RcdME->setVDaughterCouplings(aL1, aR1, 0);
  RcdME->setVDaughterCouplings(aL2, aR2, 1);
  if (verbosity >= TVar::DEBUG_VERBOSE) MELAout
    << "TUtil::JHUGenMatEl: aL1, aR1, aL2, aR2: "
    << aL1 << ", " << aR1 << ", " << aL2 << ", " << aR2
    << endl;

  // This constant is needed to account for the different units used in
  // JHUGen compared to the MCFM
  int GeVexponent_MEsq = 4-((int)mela_event.pDaughters.size())*2;
  if (production==TVar::ZZINDEPENDENT) GeVexponent_MEsq += 2; // Amplitude missing m from production and 1/m**2 from propagator == Amplitude missing 1/m == MEsq missing 1/m**2
  double constant = pow(GeV, -GeVexponent_MEsq);
  MatElSq *= constant;

  // Set RcdME information for ME and parton distributions, taking into account the mothers if id!=0 (i.e. if not unknown).
  if (mela_event.pMothers.at(0).first!=0 && mela_event.pMothers.at(1).first!=0){
    int ix=5, iy=5;
    if (abs(mela_event.pMothers.at(0).first)<6) ix=mela_event.pMothers.at(0).first + 5;
    if (abs(mela_event.pMothers.at(1).first)<6) iy=mela_event.pMothers.at(1).first + 5;
    msq[iy][ix]=MatElSq; // Note that SumMEPdf receives a transposed msq
  }
  else{
    if (production == TVar::ZZGG || production == TVar::ZZINDEPENDENT) msq[5][5]=MatElSq;
    else if (production == TVar::ZZQQB){
      for (int ix=0; ix<5; ix++){ msq[ix][10-ix]=MatElSq; msq[10-ix][ix]=MatElSq; }
    }
  }
  if (production!=TVar::ZZINDEPENDENT) SumMEPDF(MomStore[0], MomStore[1], msq, RcdME, EBEAM, verbosity);
  else{ // If production is ZZINDEPENDENT, only set gg index with fx1,2[g,g]=1.
    double fx_dummy[nmsq]={ 0 }; fx_dummy[5]=1.;
    RcdME->setPartonWeights(fx_dummy, fx_dummy);
    RcdME->setMEArray(msq, true);
    RcdME->computeWeightedMEArray();
  }

  if (verbosity >= TVar::DEBUG) MELAout << "TUtil::JHUGenMatEl: Final MatElSq = " << MatElSq << endl;

  // Reset alphas
  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::JHUGenMatEl: Reset AlphaS:\n"
      << "\tBefore reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << endl;
  }
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel);
  if (verbosity>=TVar::DEBUG){
    GetAlphaS(&alphasVal, &alphasmzVal);
    MELAout
      << "TUtil::JHUGenMatEl: Reset AlphaS result:\n"
      << "\tAfter reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }
  return MatElSq;
}

// VBF and HJJ (QCD), H undecayed, includes H BW
double TUtil::HJJMatEl(
  const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  const double& EBEAM,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double sum_msqjk = 0;
  // by default assume only gg productions
  // FOTRAN convention -5    -4   -3  -2    -1  0 1 2 3 4 5
  //     parton flavor bbar cbar sbar ubar dbar g d u s c b
  // C++ convention     0     1   2    3    4   5 6 7 8 9 10
  //2-D matrix is reversed in fortran
  // msq[ parton2 ] [ parton1 ]
  //      flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];
  double MatElsq[nmsq][nmsq]={ { 0 } };
  double MatElsq_tmp[nmsq][nmsq]={ { 0 } }; // For H+J
  double msq_tmp=0; // For "*_exact"

  if (matrixElement!=TVar::JHUGen){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::HJJMatEl: Non-JHUGen MEs are not supported" << endl; return sum_msqjk; }
  if (!(production==TVar::JJQCD || production==TVar::JJVBF || production==TVar::JQCD)){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::HJJMatEl: Production is not supported!" << endl; return sum_msqjk; }

  // Notice that partIncCode is specific for this subroutine
  int nRequested_AssociatedJets=2;
  if (production == TVar::JQCD) nRequested_AssociatedJets=1;
  int partIncCode=TVar::kUseAssociated_Jets; // Only use associated partons in the pT=0 frame boost
  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.nRequested_AssociatedJets=nRequested_AssociatedJets;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event,
    verbosity
    );
  if (mela_event.pAssociated.size()==0){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::HJJMatEl: Number of associated particles is 0!" << endl; return sum_msqjk; }

  int MYIDUP_tmp[4]={ 0 }; // "Incoming" partons 1, 2, "outgoing" partons 3, 4
  double p4[5][4]={ { 0 } };
  double pOneJet[4][4] ={ { 0 } }; // For HJ
  TLorentzVector MomStore[mxpart];
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);

  // p4(0:4,i) = (E(i),px(i),py(i),pz(i))
  // i=0,1: g1,g2 or q1, qb2 (outgoing convention)
  // i=2,3: J1, J2 in outgoing convention when possible
  // i=4: H
  for (int ipar=0; ipar<2; ipar++){
    TLorentzVector* momTmp = &(mela_event.pMothers.at(ipar).second);
    int* idtmp = &(mela_event.pMothers.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_tmp[ipar] = *idtmp;
    else MYIDUP_tmp[ipar] = 0;
    if (momTmp->T()>0.){
      p4[ipar][0] = momTmp->T()*GeV;
      p4[ipar][1] = momTmp->X()*GeV;
      p4[ipar][2] = momTmp->Y()*GeV;
      p4[ipar][3] = momTmp->Z()*GeV;
      MomStore[ipar] = (*momTmp);
    }
    else{
      p4[ipar][0] = -momTmp->T()*GeV;
      p4[ipar][1] = -momTmp->X()*GeV;
      p4[ipar][2] = -momTmp->Y()*GeV;
      p4[ipar][3] = -momTmp->Z()*GeV;
      MomStore[ipar] = -(*momTmp);
      MYIDUP_tmp[ipar] = -MYIDUP_tmp[ipar];
    }
  }
  for (unsigned int ipar=0; ipar<2; ipar++){
    if (ipar<mela_event.pAssociated.size()){
      TLorentzVector* momTmp = &(mela_event.pAssociated.at(ipar).second);
      int* idtmp = &(mela_event.pAssociated.at(ipar).first);
      if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_tmp[ipar+2] = *idtmp;
      else MYIDUP_tmp[ipar+2] = 0;
      p4[ipar+2][0] = momTmp->T()*GeV;
      p4[ipar+2][1] = momTmp->X()*GeV;
      p4[ipar+2][2] = momTmp->Y()*GeV;
      p4[ipar+2][3] = momTmp->Z()*GeV;
      MomStore[ipar+6] = (*momTmp); // i==(2, 3, 4) is (J1, J2, H), recorded as MomStore (I1, I2, 0, 0, 0, H, J1, J2)
    }
    else MYIDUP_tmp[ipar+2] = -9000; // No need to set p4, which is already 0 by initialization
  }
  for (unsigned int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
    p4[4][0] += momTmp->T()*GeV;
    p4[4][1] += momTmp->X()*GeV;
    p4[4][2] += momTmp->Y()*GeV;
    p4[4][3] += momTmp->Z()*GeV;
    MomStore[5] = MomStore[5] + (*momTmp); // i==(2, 3, 4) is (J1, J2, H), recorded as MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }
  // Momenta for HJ
  for (unsigned int i = 0; i < 4; i++){
    if (i<2){ for (unsigned int j = 0; j < 4; j++) pOneJet[i][j] = p4[i][j]; } // p1 p2
    else if (i==2){ for (unsigned int j = 0; j < 4; j++) pOneJet[i][j] = p4[4][j]; } // H
    else{ for (unsigned int j = 0; j < 4; j++) pOneJet[i][j] = p4[2][j]; } // J1
  }
  if (verbosity >= TVar::DEBUG){ 
    for (unsigned int i=0; i<5; i++) MELAout << "p["<<i<<"] (Px, Py, Pz, E, M):\t" << p4[i][1]/GeV << '\t' << p4[i][2]/GeV << '\t' << p4[i][3]/GeV << '\t' << p4[i][0]/GeV << '\t' << sqrt(fabs(pow(p4[i][0], 2)-pow(p4[i][1], 2)-pow(p4[i][2], 2)-pow(p4[i][3], 2)))/GeV << endl;
    MELAout << "One-jet momenta: " << endl;
    for (unsigned int i=0; i<4; i++) MELAout << "pOneJet["<<i<<"] (Px, Py, Pz, E, M):\t" << pOneJet[i][1]/GeV << '\t' << pOneJet[i][2]/GeV << '\t' << pOneJet[i][3]/GeV << '\t' << pOneJet[i][0]/GeV << '\t' << sqrt(fabs(pow(pOneJet[i][0], 2)-pow(pOneJet[i][1], 2)-pow(pOneJet[i][2], 2)-pow(pOneJet[i][3], 2)))/GeV << endl;
  }

  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l");
  double alphasVal, alphasmzVal;
  GetAlphaS(&alphasVal, &alphasmzVal);
  RcdME->setRenormalizationScale(renQ);
  RcdME->setFactorizationScale(facQ);
  RcdME->setAlphaS(alphasVal);
  RcdME->setAlphaSatMZ(alphasmzVal);
  RcdME->setHiggsMassWidth(masses_mcfm_.hmass, masses_mcfm_.hwidth, 0);
  RcdME->setHiggsMassWidth(spinzerohiggs_anomcoupl_.h2mass, spinzerohiggs_anomcoupl_.h2width, 1);
  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::HJJMatEl: Set AlphaS:\n"
      << "\tBefore set, alphas scale: " << defaultRenScale << ", PDF scale: " << defaultFacScale << '\n'
      << "\trenQ: " << renQ << " ( x " << event_scales->ren_scale_factor << "), facQ: " << facQ << " ( x " << event_scales->fac_scale_factor << ")\n"
      << "\tAfter set, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }

  // Since we have a lot of these checks, do them here.
  bool partonIsUnknown[4];
  for (unsigned int ip=0; ip<4; ip++) partonIsUnknown[ip] = (MYIDUP_tmp[ip]==0);

  // NOTE ON HJ AND HJJ CHANNEL HASHES:
  // THEY ONLY RETURN ISEL>=JSEL CASES. ISEL<JSEL NEEDS TO BE DONE MANUALLY.
  if (production == TVar::JQCD){ // Computation is already for all possible qqb/qg/qbg/gg, and incoming q, qb and g flavor have 1-1 correspondence to the outgoing jet flavor.
    __modhiggsj_MOD_evalamp_hj(pOneJet, MatElsq_tmp);
    for (int isel=-5; isel<=5; isel++){
      if (!partonIsUnknown[0] && !((PDGHelpers::isAGluon(MYIDUP_tmp[0]) && isel==0) || MYIDUP_tmp[0]==isel)) continue;
      for (int jsel=-5; jsel<=5; jsel++){
        if (!partonIsUnknown[1] && !((PDGHelpers::isAGluon(MYIDUP_tmp[1]) && jsel==0) || MYIDUP_tmp[1]==jsel)) continue;
        int rsel;
        if (isel!=0 && jsel!=0) rsel=0; // Covers qqb->Hg
        else if (isel==0) rsel=jsel; // Covers gg->Hg, gq->Hq, gqb->Hqb
        else rsel=isel; // Covers qg->Hq, qbg->Hqb
        if (!partonIsUnknown[2] && !((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && rsel==0) || MYIDUP_tmp[2]==rsel)) continue;
        MatElsq[jsel+5][isel+5] = MatElsq_tmp[jsel+5][isel+5]; // Assign only those that match gen. info, if present at all.
        if (verbosity >= TVar::DEBUG) MELAout << "Channel (isel, jsel)=" << isel << ", " << jsel << endl;
      }
    }
  }
  else if (production==TVar::JJQCD){
    const std::vector<TNumericUtil::intTriplet_t>& ijsel = Get_JHUGenHash_OnshellHJJHash();
    int nijchannels = ijsel.size();
    for (int ic=0; ic<nijchannels; ic++){
      // Emulate EvalWeighted_HJJ_test
      int isel = ijsel[ic][0];
      int jsel = ijsel[ic][1];
      int code = ijsel[ic][2];

      if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "HJJ channel " << ic << " code " << code << endl;

      // Default assignments
      if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "HJJ mother unswapped case" << endl;
      int rsel=isel;
      int ssel=jsel;

      if (
        (partonIsUnknown[0] || ((PDGHelpers::isAGluon(MYIDUP_tmp[0]) && isel==0) || MYIDUP_tmp[0]==isel))
        &&
        (partonIsUnknown[1] || ((PDGHelpers::isAGluon(MYIDUP_tmp[1]) && jsel==0) || MYIDUP_tmp[1]==jsel))
        ){ // Do it this way to be able to swap isel and jsel later

        if (isel==0 && jsel==0){ // gg->?
          if (code==2){ // gg->qqb
            // Only compute u-ub. The amplitude is multiplied by nf=5
            rsel=1;
            ssel=-1;
            if (
              (partonIsUnknown[2] || (PDGHelpers::isAQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0))
              &&
              (partonIsUnknown[3] || (PDGHelpers::isAQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
            }
            if (
              (partonIsUnknown[2] || (PDGHelpers::isAQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0))
              &&
              (partonIsUnknown[3] || (PDGHelpers::isAQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
            }
          }
          else{ // gg->gg
            // rsel=ssel=g already
            if (
              (partonIsUnknown[2] || (PDGHelpers::isAGluon(MYIDUP_tmp[2]) && rsel==0))
              &&
              (partonIsUnknown[3] || (PDGHelpers::isAGluon(MYIDUP_tmp[3]) && ssel==0))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
            }
          }
        }
        else if (isel==0 || jsel==0){ // qg/qbg/gq/gqb->qg/qbg/gq/gqb
          if (
            (partonIsUnknown[2] || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && rsel==0) || MYIDUP_tmp[2]==rsel))
            &&
            (partonIsUnknown[3] || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && ssel==0) || MYIDUP_tmp[3]==ssel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
          }
          if (
            (partonIsUnknown[2] || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && ssel==0) || MYIDUP_tmp[2]==ssel))
            &&
            (partonIsUnknown[3] || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && rsel==0) || MYIDUP_tmp[3]==rsel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
          }
        }
        else if ((isel>0 && jsel<0) || (isel<0 && jsel>0)){ // qQb/qbQ->?
          if (code==1 && isel==-jsel){ // qqb/qbq->gg
            rsel=0; ssel=0;
            if (
              (partonIsUnknown[2] || PDGHelpers::isAGluon(MYIDUP_tmp[2]))
              &&
              (partonIsUnknown[3] || PDGHelpers::isAGluon(MYIDUP_tmp[3]))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
            }
          }
          else if (code==3 && isel==-jsel){ // qqb->QQb
            if (abs(isel)!=1){ rsel=1; ssel=-1; } // Make sure rsel, ssel are not of same flavor as isel, jsel
            else{ rsel=2; ssel=-2; }
            // The amplitude is aready multiplied by nf-1, so no need to calculate everything (nf-1) times.
            if (
              (partonIsUnknown[2] && partonIsUnknown[3])
              ||
              (partonIsUnknown[2] && std::abs(MYIDUP_tmp[3])!=std::abs(isel) && MYIDUP_tmp[3]<0)
              ||
              (std::abs(MYIDUP_tmp[2])!=std::abs(isel) && MYIDUP_tmp[2]>0 && partonIsUnknown[3])
              ||
              (std::abs(MYIDUP_tmp[2])!=std::abs(isel) && MYIDUP_tmp[2]==-MYIDUP_tmp[3] && MYIDUP_tmp[2]>0)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
            }
            if (
              (partonIsUnknown[2] && partonIsUnknown[3])
              ||
              (partonIsUnknown[2] && std::abs(MYIDUP_tmp[3])!=std::abs(isel) && MYIDUP_tmp[3]>0)
              ||
              (std::abs(MYIDUP_tmp[2])!=std::abs(isel) && MYIDUP_tmp[2]<0 && partonIsUnknown[3])
              ||
              (std::abs(MYIDUP_tmp[2])!=std::abs(isel) && MYIDUP_tmp[2]==-MYIDUP_tmp[3] && MYIDUP_tmp[2]<0)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
            }
          }
          else{ // qQb/qbQ->qQb/qbQ
            if (
              (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
              &&
              (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
            }
            if (
              (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
              &&
              (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
            }
          }
        }
        else{
          if (
            (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
            &&
            (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
          }
          if (
            rsel!=ssel
            &&
            (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
            &&
            (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
          }
        }
      } // End unswapped isel>=jsel cases

      if (isel==jsel) continue;
      isel = ijsel[ic][1];
      jsel = ijsel[ic][0];

      if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "HJJ mother swapped case" << endl;
      // Reset to default assignments
      rsel=isel;
      ssel=jsel;

      if (
        (partonIsUnknown[0] || ((PDGHelpers::isAGluon(MYIDUP_tmp[0]) && isel==0) || MYIDUP_tmp[0]==isel))
        &&
        (partonIsUnknown[1] || ((PDGHelpers::isAGluon(MYIDUP_tmp[1]) && jsel==0) || MYIDUP_tmp[1]==jsel))
        ){
        // isel==jsel==0 is already eliminated by isel!=jsel condition
        if (isel==0 || jsel==0){ // qg/qbg/gq/gqb->qg/qbg/gq/gqb
          if (
            (partonIsUnknown[2] || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && rsel==0) || MYIDUP_tmp[2]==rsel))
            &&
            (partonIsUnknown[3] || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && ssel==0) || MYIDUP_tmp[3]==ssel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
          }
          if (
            (partonIsUnknown[2] || ((PDGHelpers::isAGluon(MYIDUP_tmp[2]) && ssel==0) || MYIDUP_tmp[2]==ssel))
            &&
            (partonIsUnknown[3] || ((PDGHelpers::isAGluon(MYIDUP_tmp[3]) && rsel==0) || MYIDUP_tmp[3]==rsel))
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
          }
        }
        else if ((isel>0 && jsel<0) || (isel<0 && jsel>0)){ // qQb/qbQ->?
          if (code==1 && isel==-jsel){ // qqb/qbq->gg
            rsel=0; ssel=0;
            if (
              (partonIsUnknown[2] || PDGHelpers::isAGluon(MYIDUP_tmp[2]))
              &&
              (partonIsUnknown[3] || PDGHelpers::isAGluon(MYIDUP_tmp[3]))
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
            }
          }
          else if(code==3 && isel==-jsel){ // qqb->QQb
            if (abs(isel)!=1){ rsel=1; ssel=-1; } // Make sure rsel, ssel are not of same flavor as isel, jsel
            else{ rsel=2; ssel=-2; }
            // The amplitude is aready multiplied by nf-1, so no need to calculate everything (nf-1) times.
            if (
              (partonIsUnknown[2] && partonIsUnknown[3])
              ||
              (partonIsUnknown[2] && std::abs(MYIDUP_tmp[3])!=std::abs(isel) && MYIDUP_tmp[3]<0)
              ||
              (std::abs(MYIDUP_tmp[2])!=std::abs(isel) && MYIDUP_tmp[2]>0 && partonIsUnknown[3])
              ||
              (std::abs(MYIDUP_tmp[2])!=std::abs(isel) && MYIDUP_tmp[2]==-MYIDUP_tmp[3] && MYIDUP_tmp[2]>0)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
            }
            if (
              (partonIsUnknown[2] && partonIsUnknown[3])
              ||
              (partonIsUnknown[2] && std::abs(MYIDUP_tmp[3])!=std::abs(isel) && MYIDUP_tmp[3]>0)
              ||
              (std::abs(MYIDUP_tmp[2])!=std::abs(isel) && MYIDUP_tmp[2]<0 && partonIsUnknown[3])
              ||
              (std::abs(MYIDUP_tmp[2])!=std::abs(isel) && MYIDUP_tmp[2]==-MYIDUP_tmp[3] && MYIDUP_tmp[2]<0)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
            }
          }
          else{ // qQb/qbQ->qQb/qbQ
            if (
              (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
              &&
              (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
            }
            if (
              (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
              &&
              (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
              ){
              __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
              MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
              if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
            }
          }
        }
        else{
          if (
            (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
            &&
            (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
          }
          if (
            rsel!=ssel
            &&
            (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
            &&
            (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
            ){
            __modhiggsjj_MOD_evalamp_sbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
          }
        }
      } // End swapped isel<jsel cases
    } // End loop over ic<nijchannels
  } // End production==TVar::JJQCD
  else if (production==TVar::JJVBF){
    int isel, jsel, rsel, ssel;
    /*
    1234/1243 [0] vs [1] correspond to transposing 12 and 34 at the same time.
    Eg. 12->34 == udbar->udbar
    1234[0] : 1=u, 2=dbar, 3=u, 4=dbar
    1243[0] : 1=u, 2=dbar, 4=u, 3=dbar
    1234[1] : 2=u, 1=dbar, 4=u, 3=dbar
    1243[1] : 2=u, 1=dbar, 3=u, 4=dbar
    */
    double ckm_takeout=1;
    double msq_uu_zz_ijrs1234[2]={ 0 };
    double msq_uu_zz_ijrs1243[2]={ 0 };
    double msq_dd_zz_ijrs1234[2]={ 0 };
    double msq_dd_zz_ijrs1243[2]={ 0 };
    double msq_ubarubar_zz_ijrs1234[2]={ 0 };
    double msq_ubarubar_zz_ijrs1243[2]={ 0 };
    double msq_dbardbar_zz_ijrs1234[2]={ 0 };
    double msq_dbardbar_zz_ijrs1243[2]={ 0 };
    double msq_uu_zzid_ijrs1234[2]={ 0 };
    double msq_uu_zzid_ijrs1243[2]={ 0 };
    double msq_dd_zzid_ijrs1234[2]={ 0 };
    double msq_dd_zzid_ijrs1243[2]={ 0 };
    double msq_ubarubar_zzid_ijrs1234[2]={ 0 };
    double msq_ubarubar_zzid_ijrs1243[2]={ 0 };
    double msq_dbardbar_zzid_ijrs1234[2]={ 0 };
    double msq_dbardbar_zzid_ijrs1243[2]={ 0 };
    double msq_udbar_zz_ijrs1234[2]={ 0 };
    double msq_udbar_zz_ijrs1243[2]={ 0 };
    double msq_dubar_zz_ijrs1234[2]={ 0 };
    double msq_dubar_zz_ijrs1243[2]={ 0 };

    double msq_uubar_zz_ijrs1234[2]={ 0 };
    double msq_uubar_zz_ijrs1243[2]={ 0 };
    double msq_ddbar_zz_ijrs1234[2]={ 0 };
    double msq_ddbar_zz_ijrs1243[2]={ 0 };
    double msq_uubar_ww_ijrs1234[2]={ 0 };
    double msq_uubar_ww_ijrs1243[2]={ 0 };
    double msq_ddbar_ww_ijrs1234[2]={ 0 };
    double msq_ddbar_ww_ijrs1243[2]={ 0 };

    double msq_ud_wwonly_ijrs1234[2]={ 0 };
    double msq_ud_wwonly_ijrs1243[2]={ 0 };
    double msq_ubardbar_wwonly_ijrs1234[2]={ 0 };
    double msq_ubardbar_wwonly_ijrs1243[2]={ 0 };
    // ud, ubardbar ZZ(+)WW case to be added separately inside a loop.
    
    // NOTE: The order should be u>c>d>s>b; particle>antiparticle
    if (
      (partonIsUnknown[0] || (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0))
      &&
      (partonIsUnknown[1] || (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0))
      &&
      (partonIsUnknown[2] || (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0))
      &&
      (partonIsUnknown[3] || (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0))
      ){
      isel=2; jsel=4; // uu'->(ZZ)->uu'
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_uu_zz_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_uu_zz_ijrs1243[0]));
      // Swapping i/j vs r/s makes no physical difference here.
      msq_uu_zz_ijrs1234[1] = msq_uu_zz_ijrs1234[0];
      msq_uu_zz_ijrs1243[1] = msq_uu_zz_ijrs1243[0];

      isel=2; jsel=2; // uu->(ZZ)->uu
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_uu_zzid_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_uu_zzid_ijrs1243[0]));
      // Swapping i/j vs r/s makes no physical difference here.
      msq_uu_zzid_ijrs1234[1] = msq_uu_zzid_ijrs1234[0];
      msq_uu_zzid_ijrs1243[1] = msq_uu_zzid_ijrs1243[0];
    }

    if (
      (partonIsUnknown[0] || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0))
      &&
      (partonIsUnknown[1] || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0))
      &&
      (partonIsUnknown[2] || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0))
      &&
      (partonIsUnknown[3] || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0))
      ){
      isel=1; jsel=3; // dd'->(ZZ)->dd'
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_dd_zz_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_dd_zz_ijrs1243[0]));
      // Swapping i/j vs r/s makes no physical difference here.
      msq_dd_zz_ijrs1234[1] = msq_dd_zz_ijrs1234[0];
      msq_dd_zz_ijrs1243[1] = msq_dd_zz_ijrs1243[0];

      isel=1; jsel=1; // dd->(ZZ)->dd
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_dd_zzid_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_dd_zzid_ijrs1243[0]));
      // Swapping i/j vs r/s makes no physical difference here.
      msq_dd_zzid_ijrs1234[1] = msq_dd_zzid_ijrs1234[0];
      msq_dd_zzid_ijrs1243[1] = msq_dd_zzid_ijrs1243[0];
    }

    if (
      (partonIsUnknown[0] || (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0))
      &&
      (partonIsUnknown[1] || (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0))
      &&
      (partonIsUnknown[2] || (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0))
      &&
      (partonIsUnknown[3] || (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0))
      ){
      isel=-2; jsel=-4; // ubarubar'->(ZZ)->ubarubar'
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_ubarubar_zz_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_ubarubar_zz_ijrs1243[0]));
      // Swapping i/j vs r/s makes no physical difference here.
      msq_ubarubar_zz_ijrs1234[1] = msq_ubarubar_zz_ijrs1234[0];
      msq_ubarubar_zz_ijrs1243[1] = msq_ubarubar_zz_ijrs1243[0];

      isel=-2; jsel=-2; // ubarubar->(ZZ)->ubarubar
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_ubarubar_zzid_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_ubarubar_zzid_ijrs1243[0]));
      // Swapping i/j vs r/s makes no physical difference here.
      msq_ubarubar_zzid_ijrs1234[1] = msq_ubarubar_zzid_ijrs1234[0];
      msq_ubarubar_zzid_ijrs1243[1] = msq_ubarubar_zzid_ijrs1243[0];
    }

    if (
      (partonIsUnknown[0] || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0))
      &&
      (partonIsUnknown[1] || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0))
      &&
      (partonIsUnknown[2] || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0))
      &&
      (partonIsUnknown[3] || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0))
      ){
      isel=-1; jsel=-3; // dbardbar'->(ZZ)->dbardbar'
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_dbardbar_zz_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_dbardbar_zz_ijrs1243[0]));
      // Swapping i/j vs r/s makes no physical difference here.
      msq_dbardbar_zz_ijrs1234[1] = msq_dbardbar_zz_ijrs1234[0];
      msq_dbardbar_zz_ijrs1243[1] = msq_dbardbar_zz_ijrs1243[0];

      isel=-1; jsel=-1; // dbardbar->(ZZ)->dbardbar
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_dbardbar_zzid_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_dbardbar_zzid_ijrs1243[0]));
      // Swapping i/j vs r/s makes no physical difference here.
      msq_dbardbar_zzid_ijrs1234[1] = msq_dbardbar_zzid_ijrs1234[0];
      msq_dbardbar_zzid_ijrs1243[1] = msq_dbardbar_zzid_ijrs1243[0];
    }

    if (
      (
      (partonIsUnknown[0] && partonIsUnknown[1])
      ||
      (
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0)
      ||
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0)
      )
      ||
      (partonIsUnknown[0] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0)))
      ||
      (partonIsUnknown[1] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0)))
      )
      &&
      (
      (partonIsUnknown[2] && partonIsUnknown[3])
      ||
      (
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0)
      ||
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0)
      )
      ||
      (partonIsUnknown[2] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0)))
      ||
      (partonIsUnknown[3] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0)))
      )
      ){
      isel=2; jsel=-1; // udbar->(ZZ)->udbar
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_udbar_zz_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_udbar_zz_ijrs1243[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &ssel, &rsel, &(msq_udbar_zz_ijrs1234[1]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &rsel, &ssel, &(msq_udbar_zz_ijrs1243[1]));
    }

    if (
      (
      (partonIsUnknown[0] && partonIsUnknown[1])
      ||
      (
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0)
      ||
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0)
      )
      ||
      (partonIsUnknown[0] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0)))
      ||
      (partonIsUnknown[1] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0)))
      )
      &&
      (
      (partonIsUnknown[2] && partonIsUnknown[3])
      ||
      (
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0)
      ||
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0)
      )
      ||
      (partonIsUnknown[2] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0)))
      ||
      (partonIsUnknown[3] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0)))
      )
      ){
      isel=1; jsel=-2; // dubar->(ZZ)->dubar
      rsel=isel; ssel=jsel;
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_dubar_zz_ijrs1234[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_dubar_zz_ijrs1243[0]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &ssel, &rsel, &(msq_dubar_zz_ijrs1234[1]));
      __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &rsel, &ssel, &(msq_dubar_zz_ijrs1243[1]));
    }

    if (
      (partonIsUnknown[0] && partonIsUnknown[1])
      ||
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[0]*MYIDUP_tmp[1]<0)
      ||
      (partonIsUnknown[0] && PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]))
      ||
      (partonIsUnknown[1] && PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]))
      ){
      isel=2; jsel=-2; // uubar->(ZZ)->uubar
      if (
        (partonIsUnknown[2] && partonIsUnknown[3])
        ||
        (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[2]*MYIDUP_tmp[3]<0)
        ||
        (partonIsUnknown[2] && PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]))
        ||
        (partonIsUnknown[3] && PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]))
        ){
        rsel=isel; ssel=jsel; // uubar->(ZZ)->uubar
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_uubar_zz_ijrs1234[0]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_uubar_zz_ijrs1243[0]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &ssel, &rsel, &(msq_uubar_zz_ijrs1234[1]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &rsel, &ssel, &(msq_uubar_zz_ijrs1243[1]));
      }
      if (
        (partonIsUnknown[2] && partonIsUnknown[3])
        ||
        (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[2]*MYIDUP_tmp[3]<0)
        ||
        (partonIsUnknown[2] && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]))
        ||
        (partonIsUnknown[3] && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]))
        ){
        rsel=-1; ssel=1; // uubar->(WW)->ddbar
        double ckm_ir = 0.;
        double ckm_js = 0.;
        while (ckm_ir==0.){
          rsel += 2;
          ssel -= 2;
          if (abs(rsel)>=6){
            rsel=-1; ssel=1;
            isel+=2; jsel-=2;
            continue;
          }
          if (abs(isel)>=6) break;
          ckm_ir = __modparameters_MOD_ckmbare(&isel, &rsel);
          ckm_js = __modparameters_MOD_ckmbare(&jsel, &ssel);
        }
        if (ckm_ir!=0.){
          __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_uubar_ww_ijrs1234[0]));
          __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_uubar_ww_ijrs1243[0]));
          __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &ssel, &rsel, &(msq_uubar_ww_ijrs1234[1]));
          __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &rsel, &ssel, &(msq_uubar_ww_ijrs1243[1]));
          ckm_takeout = pow(ckm_ir*ckm_js, 2);
          for (unsigned int iswap=0; iswap<2; iswap++){
            msq_uubar_ww_ijrs1234[iswap] /= ckm_takeout;
            msq_uubar_ww_ijrs1243[iswap] /= ckm_takeout;
          }
        }
      }
    }

    if (
      (partonIsUnknown[0] && partonIsUnknown[1])
      ||
      (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[0]*MYIDUP_tmp[1]<0)
      ||
      (partonIsUnknown[0] && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]))
      ||
      (partonIsUnknown[1] && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]))
      ){
      isel=1; jsel=-1;
      if (
        (partonIsUnknown[2] && partonIsUnknown[3])
        ||
        (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[2]*MYIDUP_tmp[3]<0)
        ||
        (partonIsUnknown[2] && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]))
        ||
        (partonIsUnknown[3] && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]))
        ){
        rsel=isel; ssel=jsel; // ddbar->(ZZ)->ddbar
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_ddbar_zz_ijrs1234[0]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_ddbar_zz_ijrs1243[0]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &ssel, &rsel, &(msq_ddbar_zz_ijrs1234[1]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &rsel, &ssel, &(msq_ddbar_zz_ijrs1243[1]));
      }
      if (
        (partonIsUnknown[2] && partonIsUnknown[3])
        ||
        (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[2]*MYIDUP_tmp[3]<0)
        ||
        (partonIsUnknown[2] && PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]))
        ||
        (partonIsUnknown[3] && PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]))
        ){
        rsel=0; ssel=0; // ddbar->(WW)->uubar
        double ckm_ir = 0.;
        double ckm_js = 0.;
        while (ckm_ir==0.){
          rsel += 2;
          ssel -= 2;
          if (abs(rsel)>=6){
            rsel=0; ssel=0;
            isel+=2; jsel-=2;
            continue;
          }
          if (abs(isel)>=6) break;
          ckm_ir = __modparameters_MOD_ckmbare(&isel, &rsel);
          ckm_js = __modparameters_MOD_ckmbare(&jsel, &ssel);
        }
        if (ckm_ir!=0.){
          __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_ddbar_ww_ijrs1234[0]));
          __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_ddbar_ww_ijrs1243[0]));
          __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &ssel, &rsel, &(msq_ddbar_ww_ijrs1234[1]));
          __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &rsel, &ssel, &(msq_ddbar_ww_ijrs1243[1]));
          ckm_takeout = pow(ckm_ir*ckm_js, 2);
          for (unsigned int iswap=0; iswap<2; iswap++){
            msq_ddbar_ww_ijrs1234[iswap] /= ckm_takeout;
            msq_ddbar_ww_ijrs1243[iswap] /= ckm_takeout;
          }
        }
      }
    }

    if (
      (
      (partonIsUnknown[0] && partonIsUnknown[1])
      ||
      (
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0)
      ||
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0)
      )
      ||
      (partonIsUnknown[0] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]>0)))
      ||
      (partonIsUnknown[1] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]>0)))
      )
      &&
      (
      (partonIsUnknown[2] && partonIsUnknown[3])
      ||
      (
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0)
      ||
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0)
      )
      ||
      (partonIsUnknown[2] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]>0)))
      ||
      (partonIsUnknown[3] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]>0)))
      )
      ){
      // ud->(WW)->d'u'
      double ckm_ir = 0.;
      double ckm_js = 0.;
      const unsigned int nCombinations = 6;
      int ijrs_arr[nCombinations][4] ={
        { 2, 3, 1, 4 },
        { 2, 3, 5, 4 },
        { 2, 5, 1, 4 },
        { 2, 5, 3, 4 },
        { 2, 1, 3, 4 },
        { 2, 1, 5, 4 }
      };
      for (unsigned int ijrs=0; ijrs<nCombinations; ijrs++){
        isel = ijrs_arr[ijrs][0];
        jsel = ijrs_arr[ijrs][1];
        rsel = ijrs_arr[ijrs][2];
        ssel = ijrs_arr[ijrs][3];
        ckm_ir = __modparameters_MOD_ckmbare(&isel, &rsel);
        ckm_js = __modparameters_MOD_ckmbare(&jsel, &ssel);
        if (ckm_ir!=0. && ckm_js!=0.) break;
      }
      if (ckm_ir!=0. && ckm_js!=0.){
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_ud_wwonly_ijrs1234[0]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_ud_wwonly_ijrs1243[0]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &ssel, &rsel, &(msq_ud_wwonly_ijrs1234[1]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &rsel, &ssel, &(msq_ud_wwonly_ijrs1243[1]));
        ckm_takeout = pow(ckm_ir*ckm_js, 2);
        for (unsigned int iswap=0; iswap<2; iswap++){
          msq_ud_wwonly_ijrs1234[iswap] /= ckm_takeout;
          msq_ud_wwonly_ijrs1243[iswap] /= ckm_takeout;
        }
      }
    }

    if (
      (
      (partonIsUnknown[0] && partonIsUnknown[1])
      ||
      (
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0)
      ||
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0)
      )
      ||
      (partonIsUnknown[0] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[1]) && MYIDUP_tmp[1]<0)))
      ||
      (partonIsUnknown[1] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[0]) && MYIDUP_tmp[0]<0)))
      )
      &&
      (
      (partonIsUnknown[2] && partonIsUnknown[3])
      ||
      (
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0)
      ||
      (PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0 && PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0)
      )
      ||
      (partonIsUnknown[2] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[3]) && MYIDUP_tmp[3]<0)))
      ||
      (partonIsUnknown[3] && ((PDGHelpers::isUpTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0) || (PDGHelpers::isDownTypeQuark(MYIDUP_tmp[2]) && MYIDUP_tmp[2]<0)))
      )
      ){
      // ubardbar->(WW)->dbar'ubar'
      double ckm_ir = 0.;
      double ckm_js = 0.;
      const unsigned int nCombinations = 6;
      int ijrs_arr[nCombinations][4] ={
        { -2, -3, -1, -4 },
        { -2, -3, -5, -4 },
        { -2, -5, -1, -4 },
        { -2, -5, -3, -4 },
        { -2, -1, -3, -4 },
        { -2, -1, -5, -4 }
      };
      for (unsigned int ijrs=0; ijrs<nCombinations; ijrs++){
        isel = ijrs_arr[ijrs][0];
        jsel = ijrs_arr[ijrs][1];
        rsel = ijrs_arr[ijrs][2];
        ssel = ijrs_arr[ijrs][3];
        ckm_ir = __modparameters_MOD_ckmbare(&isel, &rsel);
        ckm_js = __modparameters_MOD_ckmbare(&jsel, &ssel);
        if (ckm_ir!=0. && ckm_js!=0.) break;
      }
      if (ckm_ir!=0. && ckm_js!=0.){
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &(msq_ubardbar_wwonly_ijrs1234[0]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &(msq_ubardbar_wwonly_ijrs1243[0]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &ssel, &rsel, &(msq_ubardbar_wwonly_ijrs1234[1]));
        __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &jsel, &isel, &rsel, &ssel, &(msq_ubardbar_wwonly_ijrs1243[1]));
        ckm_takeout = pow(ckm_ir*ckm_js, 2);
        for (unsigned int iswap=0; iswap<2; iswap++){
          msq_ubardbar_wwonly_ijrs1234[iswap] /= ckm_takeout;
          msq_ubardbar_wwonly_ijrs1243[iswap] /= ckm_takeout;
        }
      }
    }

    if (verbosity>=TVar::DEBUG_VERBOSE){
      MELAout << "TUtil::HJJMatEl: The pre-computed MEs:" << endl;
      for (unsigned int iswap=0; iswap<2; iswap++){
        MELAout
          << "\tmsq_uu_zz_ijrs1234[" << iswap << "] = " << msq_uu_zz_ijrs1234[iswap] << '\n'
          << "\tmsq_uu_zz_ijrs1243[" << iswap << "] = " << msq_uu_zz_ijrs1243[iswap] << '\n'
          << "\tmsq_dd_zz_ijrs1234[" << iswap << "] = " << msq_dd_zz_ijrs1234[iswap] << '\n'
          << "\tmsq_dd_zz_ijrs1243[" << iswap << "] = " << msq_dd_zz_ijrs1243[iswap] << '\n'
          << "\tmsq_ubarubar_zz_ijrs1234[" << iswap << "] = " << msq_ubarubar_zz_ijrs1234[iswap] << '\n'
          << "\tmsq_ubarubar_zz_ijrs1243[" << iswap << "] = " << msq_ubarubar_zz_ijrs1243[iswap] << '\n'
          << "\tmsq_dbardbar_zz_ijrs1234[" << iswap << "] = " << msq_dbardbar_zz_ijrs1234[iswap] << '\n'
          << "\tmsq_dbardbar_zz_ijrs1243[" << iswap << "] = " << msq_dbardbar_zz_ijrs1243[iswap] << '\n'
          << "\tmsq_uu_zzid_ijrs1234[" << iswap << "] = " << msq_uu_zzid_ijrs1234[iswap] << '\n'
          << "\tmsq_uu_zzid_ijrs1243[" << iswap << "] = " << msq_uu_zzid_ijrs1243[iswap] << '\n'
          << "\tmsq_dd_zzid_ijrs1234[" << iswap << "] = " << msq_dd_zzid_ijrs1234[iswap] << '\n'
          << "\tmsq_dd_zzid_ijrs1243[" << iswap << "] = " << msq_dd_zzid_ijrs1243[iswap] << '\n'
          << "\tmsq_ubarubar_zzid_ijrs1234[" << iswap << "] = " << msq_ubarubar_zzid_ijrs1234[iswap] << '\n'
          << "\tmsq_ubarubar_zzid_ijrs1243[" << iswap << "] = " << msq_ubarubar_zzid_ijrs1243[iswap] << '\n'
          << "\tmsq_dbardbar_zzid_ijrs1234[" << iswap << "] = " << msq_dbardbar_zzid_ijrs1234[iswap] << '\n'
          << "\tmsq_dbardbar_zzid_ijrs1243[" << iswap << "] = " << msq_dbardbar_zzid_ijrs1243[iswap] << '\n'
          << "\tmsq_udbar_zz_ijrs1234[" << iswap << "] = " << msq_udbar_zz_ijrs1234[iswap] << '\n'
          << "\tmsq_udbar_zz_ijrs1243[" << iswap << "] = " << msq_udbar_zz_ijrs1243[iswap] << '\n'
          << "\tmsq_dubar_zz_ijrs1234[" << iswap << "] = " << msq_dubar_zz_ijrs1234[iswap] << '\n'
          << "\tmsq_dubar_zz_ijrs1243[" << iswap << "] = " << msq_dubar_zz_ijrs1243[iswap] << '\n'
          << "\tmsq_uubar_zz_ijrs1234[" << iswap << "] = " << msq_uubar_zz_ijrs1234[iswap] << '\n'
          << "\tmsq_uubar_zz_ijrs1243[" << iswap << "] = " << msq_uubar_zz_ijrs1243[iswap] << '\n'
          << "\tmsq_ddbar_zz_ijrs1234[" << iswap << "] = " << msq_ddbar_zz_ijrs1234[iswap] << '\n'
          << "\tmsq_ddbar_zz_ijrs1243[" << iswap << "] = " << msq_ddbar_zz_ijrs1243[iswap] << '\n'
          << "\tmsq_uubar_ww_ijrs1234[" << iswap << "] = " << msq_uubar_ww_ijrs1234[iswap] << '\n'
          << "\tmsq_uubar_ww_ijrs1243[" << iswap << "] = " << msq_uubar_ww_ijrs1243[iswap] << '\n'
          << "\tmsq_ddbar_ww_ijrs1234[" << iswap << "] = " << msq_ddbar_ww_ijrs1234[iswap] << '\n'
          << "\tmsq_ddbar_ww_ijrs1243[" << iswap << "] = " << msq_ddbar_ww_ijrs1243[iswap] << '\n'
          << "\tmsq_ud_wwonly_ijrs1234[" << iswap << "] = " << msq_ud_wwonly_ijrs1234[iswap] << '\n'
          << "\tmsq_ud_wwonly_ijrs1243[" << iswap << "] = " << msq_ud_wwonly_ijrs1243[iswap] << '\n'
          << "\tmsq_ubardbar_wwonly_ijrs1234[" << iswap << "] = " << msq_ubardbar_wwonly_ijrs1234[iswap] << '\n'
          << "\tmsq_ubardbar_wwonly_ijrs1243[" << iswap << "] = " << msq_ubardbar_wwonly_ijrs1243[iswap]
          << endl;
      }
    }

    const std::vector<TNumericUtil::intTriplet_t>& ijsel = Get_JHUGenHash_OnshellVBFHash();
    int nijchannels = ijsel.size();
    // BEGIN COMPUTATION
    for (int ic=0; ic<nijchannels; ic++){
      // Emulate EvalWeighted_HJJ_test
      isel = ijsel[ic][0];
      jsel = ijsel[ic][1];
      int code = ijsel[ic][2];
      bool ijselIsUpType[2];
      bool ijselIsDownType[2];
      bool ijselIsParticle[2];
      bool ijselIsAntiparticle[2];

      if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "VBF channel " << ic << " code " << code << endl;

      // Default assignments
      if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "VBF mother unswapped case" << endl;
      // Set r=i, s=j default values
      rsel=isel;
      ssel=jsel;
      // Set i-j bools
      ijselIsUpType[0] = (PDGHelpers::isUpTypeQuark(isel));
      ijselIsUpType[1] = (PDGHelpers::isUpTypeQuark(jsel));
      ijselIsDownType[0] = (PDGHelpers::isDownTypeQuark(isel));
      ijselIsDownType[1] = (PDGHelpers::isDownTypeQuark(jsel));
      ijselIsParticle[0] = (isel>0); ijselIsAntiparticle[0] = (isel<0);
      ijselIsParticle[1] = (jsel>0); ijselIsAntiparticle[1] = (jsel<0);

      if (
        (partonIsUnknown[0] || MYIDUP_tmp[0]==isel)
        &&
        (partonIsUnknown[1] || MYIDUP_tmp[1]==jsel)
        ){ // Do it this way to be able to swap isel and jsel later
        if (code==1){ // Only ZZ->H possible
          // rsel=isel and ssel=jsel already
          if (
            (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
            &&
            (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
            ){
            msq_tmp=0;
            if (ijselIsUpType[0] && ijselIsUpType[1]){
              if (ijselIsParticle[0] && ijselIsParticle[1]){
                if (isel!=jsel) msq_tmp = msq_uu_zz_ijrs1234[0];
                else msq_tmp = msq_uu_zzid_ijrs1234[0];
              }
              else if (ijselIsAntiparticle[0] && ijselIsAntiparticle[1]){
                if (isel!=jsel) msq_tmp = msq_ubarubar_zz_ijrs1234[0];
                else msq_tmp = msq_ubarubar_zzid_ijrs1234[0];
              }
              else if (ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_uubar_zz_ijrs1234[0];
            }
            else if (ijselIsDownType[0] && ijselIsDownType[1]){
              if (ijselIsParticle[0] && ijselIsParticle[1]){
                if (isel!=jsel) msq_tmp = msq_dd_zz_ijrs1234[0];
                else msq_tmp = msq_dd_zzid_ijrs1234[0];
              }
              else if (ijselIsAntiparticle[0] && ijselIsAntiparticle[1]){
                if (isel!=jsel) msq_tmp = msq_dbardbar_zz_ijrs1234[0];
                else msq_tmp = msq_dbardbar_zzid_ijrs1234[0];
              }
              else if (ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_ddbar_zz_ijrs1234[0];
            }
            else if (ijselIsUpType[0] && ijselIsDownType[1] && ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_udbar_zz_ijrs1234[0];
            else if (ijselIsDownType[0] && ijselIsUpType[1] && ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_dubar_zz_ijrs1234[0];
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
          }
          if (
            rsel!=ssel
            &&
            (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
            &&
            (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
            ){
            msq_tmp=0;
            if (ijselIsUpType[0] && ijselIsUpType[1]){
              if (ijselIsParticle[0] && ijselIsParticle[1]){
                if (isel!=jsel) msq_tmp = msq_uu_zz_ijrs1243[0];
                else msq_tmp = msq_uu_zzid_ijrs1243[0];
              }
              else if (ijselIsAntiparticle[0] && ijselIsAntiparticle[1]){
                if (isel!=jsel) msq_tmp = msq_ubarubar_zz_ijrs1243[0];
                else msq_tmp = msq_ubarubar_zzid_ijrs1243[0];
              }
              else if (ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_uubar_zz_ijrs1243[0];
            }
            else if (ijselIsDownType[0] && ijselIsDownType[1]){
              if (ijselIsParticle[0] && ijselIsParticle[1]){
                if (isel!=jsel) msq_tmp = msq_dd_zz_ijrs1243[0];
                else msq_tmp = msq_dd_zzid_ijrs1243[0];
              }
              else if (ijselIsAntiparticle[0] && ijselIsAntiparticle[1]){
                if (isel!=jsel) msq_tmp = msq_dbardbar_zz_ijrs1243[0];
                else msq_tmp = msq_dbardbar_zzid_ijrs1243[0];
              }
              else if (ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_ddbar_zz_ijrs1243[0];
            }
            else if (ijselIsUpType[0] && ijselIsDownType[1] && ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_udbar_zz_ijrs1243[0];
            else if (ijselIsDownType[0] && ijselIsUpType[1] && ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_dubar_zz_ijrs1243[0];
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
          }
        }
        else if (code==0){ // code==0 means WW->H is also possible with no interference to ZZ->H, for example u ub -> d db.
          vector<int> possible_rsel;
          vector<int> possible_ssel;
          vector<double> possible_Vsqir;
          vector<double> possible_Vsqjs;
          if (ijselIsUpType[0]){ possible_rsel.push_back(1); possible_rsel.push_back(3); possible_rsel.push_back(5); }
          else if (ijselIsDownType[0]){ possible_rsel.push_back(2); possible_rsel.push_back(4); }
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            rsel=possible_rsel.at(ix);
            double ckmval = pow(__modparameters_MOD_ckmbare(&isel, &rsel), 2);
            possible_Vsqir.push_back(ckmval);
          }
          if (ijselIsUpType[1]){ possible_ssel.push_back(1); possible_ssel.push_back(3); possible_ssel.push_back(5); }
          else if (ijselIsDownType[1]){ possible_ssel.push_back(2); possible_ssel.push_back(4); }
          for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
            ssel=possible_ssel.at(iy);
            double ckmval = pow(__modparameters_MOD_ckmbare(&jsel, &ssel), 2);
            possible_Vsqjs.push_back(ckmval);
          }

          // Combine the ME and ME_swap based on actual ids
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            rsel=possible_rsel.at(ix)*TMath::Sign(1, isel);
            for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
              ssel=possible_ssel.at(iy)*TMath::Sign(1, jsel);
              double ckmval = possible_Vsqir.at(ix)*possible_Vsqjs.at(iy);
              if (
                (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
                &&
                (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
                ){
                msq_tmp=0;
                if (ijselIsUpType[0] && ijselIsUpType[1] && ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_uubar_ww_ijrs1234[0];
                else if (ijselIsDownType[0] && ijselIsDownType[1] && ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_ddbar_ww_ijrs1234[0];
                msq_tmp *= ckmval;
                MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
              }
              if (
                rsel!=ssel
                &&
                (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
                &&
                (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
                ){
                msq_tmp=0;
                if (ijselIsUpType[0] && ijselIsUpType[1] && ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_uubar_ww_ijrs1243[0];
                else if (ijselIsDownType[0] && ijselIsDownType[1] && ijselIsParticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_ddbar_ww_ijrs1243[0];
                msq_tmp *= ckmval;
                MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
              }
            }
          }
        }
        else{ // code==2 means states with WW/ZZ interference allowed, for example u1 d2 -> u3 d4 + d4 u3.
          vector<int> possible_rsel;
          vector<int> possible_ssel;
          vector<double> possible_Vsqir;
          vector<double> possible_Vsqjs;
          if (ijselIsUpType[0]){ possible_rsel.push_back(1); possible_rsel.push_back(3); possible_rsel.push_back(5); }
          else if (ijselIsDownType[0]){ possible_rsel.push_back(2); possible_rsel.push_back(4); }
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            rsel=possible_rsel.at(ix);
            double ckmval = pow(__modparameters_MOD_ckmbare(&isel, &rsel), 2);
            possible_Vsqir.push_back(ckmval);
          }
          if (ijselIsUpType[1]){ possible_ssel.push_back(1); possible_ssel.push_back(3); possible_ssel.push_back(5); }
          else if (ijselIsDownType[1]){ possible_ssel.push_back(2); possible_ssel.push_back(4); }
          for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
            ssel=possible_ssel.at(iy);
            double ckmval = pow(__modparameters_MOD_ckmbare(&jsel, &ssel), 2);
            possible_Vsqjs.push_back(ckmval);
          }

          // Loop over all possible combinations to get interference correct
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            rsel=possible_rsel.at(ix)*TMath::Sign(1, isel);
            for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
              ssel=possible_ssel.at(iy)*TMath::Sign(1, jsel);
              double ckmval = possible_Vsqir.at(ix)*possible_Vsqjs.at(iy);
              if (
                (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
                &&
                (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
                ){
                if (rsel!=jsel && ssel!=isel){
                  msq_tmp=0;
                  if (ijselIsUpType[0] && ijselIsDownType[1]){
                    if (ijselIsParticle[0] && ijselIsParticle[1]) msq_tmp = msq_ud_wwonly_ijrs1234[0];
                    else if (ijselIsAntiparticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_ubardbar_wwonly_ijrs1234[0];
                  }
                  msq_tmp *= ckmval;
                  MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                  if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
                }
                else{
                  __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
                  MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                  if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
                }
              }
              if (
                rsel!=ssel
                &&
                (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
                &&
                (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
                ){
                if (rsel!=jsel && ssel!=isel){
                  msq_tmp=0;
                  if (ijselIsUpType[0] && ijselIsDownType[1]){
                    if (ijselIsParticle[0] && ijselIsParticle[1]) msq_tmp = msq_ud_wwonly_ijrs1243[0];
                    else if (ijselIsAntiparticle[0] && ijselIsAntiparticle[1]) msq_tmp = msq_ubardbar_wwonly_ijrs1243[0];
                  }
                  msq_tmp *= ckmval;
                  MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                  if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
                }
                else{
                  __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
                  MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                  if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
                }
              }
            }
          }
        }
      } // End unswapped isel>=jsel cases

      if (isel==jsel) continue;
      isel = ijsel[ic][1];
      jsel = ijsel[ic][0];

      // Reset to default assignments
      if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "VBF mother swapped case" << endl;
      // Set r=i, s=j default values
      rsel=isel;
      ssel=jsel;
      // Set i-j bools
      ijselIsUpType[0] = (PDGHelpers::isUpTypeQuark(isel));
      ijselIsUpType[1] = (PDGHelpers::isUpTypeQuark(jsel));
      ijselIsDownType[0] = (PDGHelpers::isDownTypeQuark(isel));
      ijselIsDownType[1] = (PDGHelpers::isDownTypeQuark(jsel));
      ijselIsParticle[0] = (isel>0); ijselIsAntiparticle[0] = (isel<0);
      ijselIsParticle[1] = (jsel>0); ijselIsAntiparticle[1] = (jsel<0);

      if (
        (partonIsUnknown[0] || ((PDGHelpers::isAGluon(MYIDUP_tmp[0]) && isel==0) || MYIDUP_tmp[0]==isel))
        &&
        (partonIsUnknown[1] || ((PDGHelpers::isAGluon(MYIDUP_tmp[1]) && jsel==0) || MYIDUP_tmp[1]==jsel))
        ){
        // isel==jsel==0 is already eliminated by isel!=jsel condition
        if (code==1){ // Only ZZ->H possible
          // rsel=isel and ssel=jsel already
          if (
            (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
            &&
            (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
            ){
            msq_tmp=0;
            if (ijselIsUpType[1] && ijselIsUpType[0]){
              if (ijselIsParticle[1] && ijselIsParticle[0]){
                if (isel!=jsel) msq_tmp = msq_uu_zz_ijrs1234[1];
                else msq_tmp = msq_uu_zzid_ijrs1234[1];
              }
              else if (ijselIsAntiparticle[1] && ijselIsAntiparticle[0]){
                if (isel!=jsel) msq_tmp = msq_ubarubar_zz_ijrs1234[1];
                else msq_tmp = msq_ubarubar_zzid_ijrs1234[1];
              }
              else if (ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_uubar_zz_ijrs1234[1];
            }
            else if (ijselIsDownType[1] && ijselIsDownType[0]){
              if (ijselIsParticle[1] && ijselIsParticle[0]){
                if (isel!=jsel) msq_tmp = msq_dd_zz_ijrs1234[1];
                else msq_tmp = msq_dd_zzid_ijrs1234[1];
              }
              else if (ijselIsAntiparticle[1] && ijselIsAntiparticle[0]){
                if (isel!=jsel) msq_tmp = msq_dbardbar_zz_ijrs1234[1];
                else msq_tmp = msq_dbardbar_zzid_ijrs1234[1];
              }
              else if (ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_ddbar_zz_ijrs1234[1];
            }
            else if (ijselIsUpType[1] && ijselIsDownType[0] && ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_udbar_zz_ijrs1234[1];
            else if (ijselIsDownType[1] && ijselIsUpType[0] && ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_dubar_zz_ijrs1234[1];
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
          }
          if (
            rsel!=ssel
            &&
            (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
            &&
            (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
            ){
            msq_tmp=0;
            if (ijselIsUpType[1] && ijselIsUpType[0]){
              if (ijselIsParticle[1] && ijselIsParticle[0]){
                if (isel!=jsel) msq_tmp = msq_uu_zz_ijrs1243[1];
                else msq_tmp = msq_uu_zzid_ijrs1243[1];
              }
              else if (ijselIsAntiparticle[1] && ijselIsAntiparticle[0]){
                if (isel!=jsel) msq_tmp = msq_ubarubar_zz_ijrs1243[1];
                else msq_tmp = msq_ubarubar_zzid_ijrs1243[1];
              }
              else if (ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_uubar_zz_ijrs1243[1];
            }
            else if (ijselIsDownType[1] && ijselIsDownType[0]){
              if (ijselIsParticle[1] && ijselIsParticle[0]){
                if (isel!=jsel) msq_tmp = msq_dd_zz_ijrs1243[1];
                else msq_tmp = msq_dd_zzid_ijrs1243[1];
              }
              else if (ijselIsAntiparticle[1] && ijselIsAntiparticle[0]){
                if (isel!=jsel) msq_tmp = msq_dbardbar_zz_ijrs1243[1];
                else msq_tmp = msq_dbardbar_zzid_ijrs1243[1];
              }
              else if (ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_ddbar_zz_ijrs1243[1];
            }
            else if (ijselIsUpType[1] && ijselIsDownType[0] && ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_udbar_zz_ijrs1243[1];
            else if (ijselIsDownType[1] && ijselIsUpType[0] && ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_dubar_zz_ijrs1243[1];
            MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
            if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
          }
        }
        else if (code==0){ // code==0 means WW->H is also possible with no interference to ZZ->H, for example u ub -> d db.
          vector<int> possible_rsel;
          vector<int> possible_ssel;
          vector<double> possible_Vsqir;
          vector<double> possible_Vsqjs;
          if (ijselIsUpType[0]){ possible_rsel.push_back(1); possible_rsel.push_back(3); possible_rsel.push_back(5); }
          else if (ijselIsDownType[0]){ possible_rsel.push_back(2); possible_rsel.push_back(4); }
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            rsel=possible_rsel.at(ix);
            double ckmval = pow(__modparameters_MOD_ckmbare(&isel, &rsel), 2);
            possible_Vsqir.push_back(ckmval);
          }
          if (ijselIsUpType[1]){ possible_ssel.push_back(1); possible_ssel.push_back(3); possible_ssel.push_back(5); }
          else if (ijselIsDownType[1]){ possible_ssel.push_back(2); possible_ssel.push_back(4); }
          for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
            ssel=possible_ssel.at(iy);
            double ckmval = pow(__modparameters_MOD_ckmbare(&jsel, &ssel), 2);
            possible_Vsqjs.push_back(ckmval);
          }

          // Combine the ME and ME_swap based on actual ids
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            rsel=possible_rsel.at(ix)*TMath::Sign(1, isel);
            for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
              ssel=possible_ssel.at(iy)*TMath::Sign(1, jsel);
              double ckmval = possible_Vsqir.at(ix)*possible_Vsqjs.at(iy);
              if (
                (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
                &&
                (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
                ){
                msq_tmp=0;
                if (ijselIsUpType[1] && ijselIsUpType[0] && ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_uubar_ww_ijrs1234[1];
                else if (ijselIsDownType[1] && ijselIsDownType[0] && ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_ddbar_ww_ijrs1234[1];
                msq_tmp *= ckmval;
                MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
              }
              if (
                rsel!=ssel
                &&
                (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
                &&
                (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
                ){
                msq_tmp=0;
                if (ijselIsUpType[1] && ijselIsUpType[0] && ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_uubar_ww_ijrs1243[1];
                else if (ijselIsDownType[1] && ijselIsDownType[0] && ijselIsParticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_ddbar_ww_ijrs1243[1];
                msq_tmp *= ckmval;
                MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
              }
            }
          }
        }
        else{ // code==2 means states with WW/ZZ interference allowed, for example u1 d2 -> u3 d4 + d4 u3.
          vector<int> possible_rsel;
          vector<int> possible_ssel;
          vector<double> possible_Vsqir;
          vector<double> possible_Vsqjs;
          if (ijselIsUpType[0]){ possible_rsel.push_back(1); possible_rsel.push_back(3); possible_rsel.push_back(5); }
          else if (ijselIsDownType[0]){ possible_rsel.push_back(2); possible_rsel.push_back(4); }
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            rsel=possible_rsel.at(ix);
            double ckmval = pow(__modparameters_MOD_ckmbare(&isel, &rsel), 2);
            possible_Vsqir.push_back(ckmval);
          }
          if (ijselIsUpType[1]){ possible_ssel.push_back(1); possible_ssel.push_back(3); possible_ssel.push_back(5); }
          else if (ijselIsDownType[1]){ possible_ssel.push_back(2); possible_ssel.push_back(4); }
          for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
            ssel=possible_ssel.at(iy);
            double ckmval = pow(__modparameters_MOD_ckmbare(&jsel, &ssel), 2);
            possible_Vsqjs.push_back(ckmval);
          }

          // Loop over all possible combinations to get interference correct
          for (unsigned int ix=0; ix<possible_rsel.size(); ix++){
            rsel=possible_rsel.at(ix)*TMath::Sign(1, isel);
            for (unsigned int iy=0; iy<possible_ssel.size(); iy++){
              ssel=possible_ssel.at(iy)*TMath::Sign(1, jsel);
              double ckmval = possible_Vsqir.at(ix)*possible_Vsqjs.at(iy);
              if (
                (partonIsUnknown[2] || MYIDUP_tmp[2]==rsel)
                &&
                (partonIsUnknown[3] || MYIDUP_tmp[3]==ssel)
                ){
                if (rsel!=jsel && ssel!=isel){
                  msq_tmp=0;
                  if (ijselIsUpType[1] && ijselIsDownType[0]){
                    if (ijselIsParticle[1] && ijselIsParticle[0]) msq_tmp = msq_ud_wwonly_ijrs1234[1];
                    else if (ijselIsAntiparticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_ubardbar_wwonly_ijrs1234[1];
                  }
                  msq_tmp *= ckmval;
                  MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                  if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
                }
                else{
                  __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &rsel, &ssel, &msq_tmp);
                  MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                  if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << rsel << ", " << ssel << '\t' <<  msq_tmp << endl;
                }
              }
              if (
                rsel!=ssel
                &&
                (partonIsUnknown[2] || MYIDUP_tmp[2]==ssel)
                &&
                (partonIsUnknown[3] || MYIDUP_tmp[3]==rsel)
                ){
                if (rsel!=jsel && ssel!=isel){
                  msq_tmp=0;
                  if (ijselIsUpType[1] && ijselIsDownType[0]){
                    if (ijselIsParticle[1] && ijselIsParticle[0]) msq_tmp = msq_ud_wwonly_ijrs1243[1];
                    else if (ijselIsAntiparticle[1] && ijselIsAntiparticle[0]) msq_tmp = msq_ubardbar_wwonly_ijrs1243[1];
                  }
                  msq_tmp *= ckmval;
                  MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                  if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
                }
                else{
                  __modhiggsjj_MOD_evalamp_wbfh_unsymm_sa_select_exact(p4, &isel, &jsel, &ssel, &rsel, &msq_tmp);
                  MatElsq[jsel+5][isel+5] += msq_tmp; // Assign only those that match gen. info, if present at all.
                  if (verbosity >= TVar::DEBUG_VERBOSE) MELAout << "Channel (isel, jsel, rsel, ssel)=" << isel << ", " << jsel << ", " << ssel << ", " << rsel << '\t' <<  msq_tmp << endl;
                }
              }
            }
          }
        }
      } // End swapped isel<jsel cases
    } // End loop over ic<nijchannels
    // END COMPUTATION
  } // End production==TVar::JJVBF

  int GeVexponent_MEsq = 4-(1+nRequested_AssociatedJets)*2;
  double constant = pow(GeV, -GeVexponent_MEsq);
  for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) MatElsq[jj][ii] *= constant; }
  //    FOTRAN convention    -5    -4   -3   -2   -1    0   1   2   3  4  5
  //     parton flavor      bbar  cbar  sbar ubar dbar  g   d   u   s  c  b
  //      C++ convention     0      1    2    3    4    5   6   7   8  9  10
  if (verbosity >= TVar::DEBUG){
    MELAout << "MatElsq:\n";
    for (int ii = 0; ii < nmsq; ii++){ for (int jj = 0; jj < nmsq; jj++) MELAout << MatElsq[jj][ii] << '\t'; MELAout << endl; }
  }
  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, EBEAM, verbosity);
  /*
  if (verbosity >= TVar::ERROR && (std::isnan(sum_msqjk) || std::isinf(sum_msqjk))){
    MELAout << "TUtil::HJJMatEl: FAILURE!" << endl;
    MELAout << "MatElsq:\n";
    for (int ii = 0; ii < nmsq; ii++){ for (int jj = 0; jj < nmsq; jj++) MELAout << MatElsq[jj][ii] << '\t'; MELAout << endl; }
    double fx[2][nmsq];
    RcdME->getPartonWeights(fx[0], fx[1]);
    for (int ii = 0; ii < nmsq; ii++){ for (int jj = 0; jj < 2; jj++) MELAout << fx[jj][ii] << '\t'; MELAout << endl; }
    for (int i=0; i<5; i++) MELAout << "p["<<i<<"] (Px, Py, Pz, E, M):\t" << p4[i][1]/GeV << '\t' << p4[i][2]/GeV << '\t' << p4[i][3]/GeV << '\t' << p4[i][0]/GeV << '\t' << sqrt(fabs(pow(p4[i][0], 2)-pow(p4[i][1], 2)-pow(p4[i][2], 2)-pow(p4[i][3], 2)))/GeV << endl;
    TUtil::PrintCandidateSummary(RcdME->melaCand);
    MELAout << endl;
  }
  */

  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::HJJMatEl: Reset AlphaS:\n"
      << "\tBefore reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << endl;
  }
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel);
  if (verbosity>=TVar::DEBUG){
    GetAlphaS(&alphasVal, &alphasmzVal);
    MELAout
      << "TUtil::HJJMatEl: Reset AlphaS result:\n"
      << "\tAfter reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }
  return sum_msqjk;
}

// VH, H undecayed or Hbb
double TUtil::VHiggsMatEl(
  const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  const double& EBEAM,
  bool includeHiggsDecay,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double sum_msqjk = 0;
  // by default assume only gg productions
  // FOTRAN convention -5    -4   -3  -2    -1  0 1 2 3 4 5
  //     parton flavor bbar cbar sbar ubar dbar g d u s c b
  // C++ convention     0     1   2    3    4   5 6 7 8 9 10
  //2-D matrix is reversed in fortran
  // msq[ parton2 ] [ parton1 ]
  //      flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];
  double MatElsq[nmsq][nmsq]={ { 0 } };

  if (matrixElement!=TVar::JHUGen){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::VHiggsMatEl: Non-JHUGen MEs are not supported" << endl; return sum_msqjk; }
  if (!(production == TVar::Lep_ZH || production == TVar::Lep_WH || production == TVar::Had_ZH || production == TVar::Had_WH || production == TVar::GammaH)){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::VHiggsMatEl: Production is not supported!" << endl; return sum_msqjk; }

  int nRequested_AssociatedJets=0;
  int nRequested_AssociatedLeptons=0;
  int nRequested_AssociatedPhotons=0;
  int AssociationVCompatibility=0;
  int partIncCode=TVar::kNoAssociated; // Just to avoid warnings
  if (production == TVar::Had_ZH || production == TVar::Had_WH){ // Only use associated partons
    partIncCode=TVar::kUseAssociated_Jets;
    nRequested_AssociatedJets=2;
  }
  else if (production == TVar::Lep_ZH || production == TVar::Lep_WH){ // Only use associated leptons(+)neutrinos
    partIncCode=TVar::kUseAssociated_Leptons;
    nRequested_AssociatedLeptons=2;
  }
  else if (production == TVar::GammaH){ // Only use associated photon
    partIncCode=TVar::kUseAssociated_Photons;
    nRequested_AssociatedPhotons=1;
  }
  if (production==TVar::Lep_WH || production==TVar::Had_WH) AssociationVCompatibility=24;
  else if (production==TVar::Lep_ZH || production==TVar::Had_ZH) AssociationVCompatibility=23;
  else if (production==TVar::GammaH) AssociationVCompatibility=22;
  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.AssociationVCompatibility=AssociationVCompatibility;
  mela_event.nRequested_AssociatedJets=nRequested_AssociatedJets;
  mela_event.nRequested_AssociatedLeptons=nRequested_AssociatedLeptons;
  mela_event.nRequested_AssociatedPhotons=nRequested_AssociatedPhotons;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event,
    verbosity
    );
  if ((mela_event.pAssociated.size()<(unsigned int)(nRequested_AssociatedJets+nRequested_AssociatedLeptons) && production!=TVar::GammaH) || (mela_event.pAssociated.size()<(unsigned int)nRequested_AssociatedPhotons && production==TVar::GammaH)){
    if (verbosity>=TVar::ERROR){
      MELAerr << "TUtil::VHiggsMatEl: Number of associated particles (" << mela_event.pAssociated.size() << ") is less than ";
      if (production!=TVar::GammaH) MELAerr << (nRequested_AssociatedJets+nRequested_AssociatedLeptons);
      else MELAerr << nRequested_AssociatedPhotons;
      MELAerr << endl;
    }
    return sum_msqjk;
  }

  int MYIDUP_prod[4]={ 0 }; // "Incoming" partons 1, 2, "outgoing" partons 3, 4
  int MYIDUP_dec[2]={ -9000, -9000 }; // "Outgoing" partons 1, 2 from the Higgs (->bb)
  double p4[9][4] ={ { 0 } };
  double helicities[9] ={ 0 };
  int vh_ids[9] ={ 0 };
  TLorentzVector MomStore[mxpart];
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);

  // p4(0:8,i) = (E(i),px(i),py(i),pz(i))
  // i=0,1: q1, qb2 (outgoing convention)
  // i=2,3: V*, V
  // i=4: H
  // i=5,6: f, fb from V
  // i=7,8: b, bb from H
  for (int ipar=0; ipar<2; ipar++){
    TLorentzVector* momTmp = &(mela_event.pMothers.at(ipar).second);
    int* idtmp = &(mela_event.pMothers.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_prod[ipar] = *idtmp;
    else MYIDUP_prod[ipar] = 0;
    if (momTmp->T()>0.){
      p4[ipar][0] = momTmp->T()*GeV;
      p4[ipar][1] = momTmp->X()*GeV;
      p4[ipar][2] = momTmp->Y()*GeV;
      p4[ipar][3] = momTmp->Z()*GeV;
      MomStore[ipar] = (*momTmp);
    }
    else{
      p4[ipar][0] = -momTmp->T()*GeV;
      p4[ipar][1] = -momTmp->X()*GeV;
      p4[ipar][2] = -momTmp->Y()*GeV;
      p4[ipar][3] = -momTmp->Z()*GeV;
      MomStore[ipar] = -(*momTmp);
      MYIDUP_prod[ipar] = -MYIDUP_prod[ipar];
    }
  }
  // Associated particles
  for (int ipar=0; ipar<(production!=TVar::GammaH ? 2 : 1); ipar++){
    TLorentzVector* momTmp = &(mela_event.pAssociated.at(ipar).second);
    int* idtmp = &(mela_event.pAssociated.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_prod[ipar+2] = *idtmp;
    else MYIDUP_prod[ipar+2] = 0;
    p4[ipar+5][0] = momTmp->T()*GeV;
    p4[ipar+5][1] = momTmp->X()*GeV;
    p4[ipar+5][2] = momTmp->Y()*GeV;
    p4[ipar+5][3] = momTmp->Z()*GeV;
    MomStore[ipar+6] = (*momTmp);
  }
  if (production==TVar::GammaH) MYIDUP_prod[3]=-9000;

  if (PDGHelpers::isAGluon(MYIDUP_prod[0]) || PDGHelpers::isAGluon(MYIDUP_prod[1])){ if (verbosity>=TVar::INFO) MELAerr << "TUtil::VHiggsMatEl: Initial state gluons are not permitted!" << endl; return sum_msqjk; }
  if (PDGHelpers::isAGluon(MYIDUP_prod[2]) || PDGHelpers::isAGluon(MYIDUP_prod[3])){ if (verbosity>=TVar::INFO) MELAerr << "TUtil::VHiggsMatEl: Final state gluons are not permitted!" << endl; return sum_msqjk; }
  if (production==TVar::GammaH && !PDGHelpers::isAPhoton(MYIDUP_prod[2])){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::VHiggsMatEl: GammaH associated photon (id=" << MYIDUP_prod[2] << ") is not a photon! Please fix its id." << endl; return sum_msqjk; }

  // Decay V/f ids
  // MYIDUP_dec as size=2 because JHUGen supports b-bbar decay
  for (int iv=0; iv<2; iv++){
    int idtmp = mela_event.intermediateVid.at(iv);
    if (!PDGHelpers::isAnUnknownJet(idtmp)) MYIDUP_dec[iv] = idtmp;
    else MYIDUP_dec[iv] = 0;
  }
  // Decay daughters
  for (unsigned int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
    if (mela_event.pDaughters.size()==1){
      p4[ipar+7][0] = momTmp->T()*GeV;
      p4[ipar+7][1] = momTmp->X()*GeV;
      p4[ipar+7][2] = momTmp->Y()*GeV;
      p4[ipar+7][3] = momTmp->Z()*GeV;
      MomStore[5] = (*momTmp); // 5
    }
    else if (mela_event.pDaughters.size()==2){
      p4[ipar+7][0] = momTmp->T()*GeV;
      p4[ipar+7][1] = momTmp->X()*GeV;
      p4[ipar+7][2] = momTmp->Y()*GeV;
      p4[ipar+7][3] = momTmp->Z()*GeV;
      MomStore[2*ipar+2] = (*momTmp); // 2,4
      if (PDGHelpers::isAQuark(mela_event.pDaughters.at(ipar).first)) MYIDUP_dec[ipar]=mela_event.pDaughters.at(ipar).first;
    }
    else if (mela_event.pDaughters.size()==3){
      if (ipar<2){
        p4[7][0] += momTmp->T()*GeV;
        p4[7][1] += momTmp->X()*GeV;
        p4[7][2] += momTmp->Y()*GeV;
        p4[7][3] += momTmp->Z()*GeV;
      }
      else{
        p4[8][0] = momTmp->T()*GeV;
        p4[8][1] = momTmp->X()*GeV;
        p4[8][2] = momTmp->Y()*GeV;
        p4[8][3] = momTmp->Z()*GeV;
      }
      MomStore[ipar+2] = (*momTmp); // 2,3,4
    }
    else if (mela_event.pDaughters.size()==4){
      if (ipar<2){
        p4[7][0] += momTmp->T()*GeV;
        p4[7][1] += momTmp->X()*GeV;
        p4[7][2] += momTmp->Y()*GeV;
        p4[7][3] += momTmp->Z()*GeV;
      }
      else{
        p4[8][0] += momTmp->T()*GeV;
        p4[8][1] += momTmp->X()*GeV;
        p4[8][2] += momTmp->Y()*GeV;
        p4[8][3] += momTmp->Z()*GeV;
      }
      MomStore[ipar+2] = (*momTmp); // 2,3,4,5
    }
    else{ // Should never happen
      p4[7][0] += momTmp->T()*GeV;
      p4[7][1] += momTmp->X()*GeV;
      p4[7][2] += momTmp->Y()*GeV;
      p4[7][3] += momTmp->Z()*GeV;
      MomStore[5] = MomStore[5] + (*momTmp);
    }
  }
  for (int ix=0; ix<4; ix++){
    p4[3][ix] = p4[5][ix] + p4[6][ix];
    p4[4][ix] = p4[7][ix] + p4[8][ix];
    p4[2][ix] = p4[3][ix] + p4[4][ix];
  }

  vh_ids[4] = 25;
  if (production==TVar::Lep_ZH || production==TVar::Had_ZH || production==TVar::GammaH) vh_ids[2] = 23;
  else if (production==TVar::Lep_WH || production==TVar::Had_WH) vh_ids[2] = 24; // To be changed later
  if (production!=TVar::GammaH) vh_ids[3] = vh_ids[2]; // To be changed later for WH
  else vh_ids[3] = 22;

  // H->ffb decay is turned off, so no need to loop over helicities[7]=helicities[8]=+-1
  vh_ids[7] = 5; helicities[7] = 1;
  vh_ids[8] = -5; helicities[8] = 1;
  int HDKon = 0;
  if (includeHiggsDecay && MYIDUP_dec[0]!=-9000 && MYIDUP_dec[1]!=-9000 && MYIDUP_dec[0]==-MYIDUP_dec[1]){ // H->ffb
    HDKon=1;
    __modjhugenmela_MOD_sethdk(&HDKon);
    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::VHiggsMatEl: HDKon" << endl;
  }
  else if (verbosity>=TVar::INFO && includeHiggsDecay) MELAerr << "TUtil::VHiggsMatEl: includeHiggsDecay=true is not supported for the present decay mode." << endl;

  if (verbosity>=TVar::DEBUG){
    for (int i=0; i<9; i++) MELAout << "p4(" << vh_ids[i] << ") = "  << p4[i][0] << ", " <<  p4[i][1] << ", "  <<  p4[i][2] << ", "  <<  p4[i][3] << "\n";
  }

  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l");
  double alphasVal, alphasmzVal;
  GetAlphaS(&alphasVal, &alphasmzVal);
  RcdME->setRenormalizationScale(renQ);
  RcdME->setFactorizationScale(facQ);
  RcdME->setAlphaS(alphasVal);
  RcdME->setAlphaSatMZ(alphasmzVal);
  RcdME->setHiggsMassWidth(masses_mcfm_.hmass, masses_mcfm_.hwidth, 0);
  RcdME->setHiggsMassWidth(spinzerohiggs_anomcoupl_.h2mass, spinzerohiggs_anomcoupl_.h2width, 1);
  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::VHiggsMatEl: Set AlphaS:\n"
      << "\tBefore set, alphas scale: " << defaultRenScale << ", PDF scale: " << defaultFacScale << '\n'
      << "\trenQ: " << renQ << " ( x " << event_scales->ren_scale_factor << "), facQ: " << facQ << " ( x " << event_scales->fac_scale_factor << ")\n"
      << "\tAfter set, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }

  // Since we have a lot of these checks, do them here.
  bool partonIsKnown[4];
  for (unsigned int ip=0; ip<4; ip++) partonIsKnown[ip] = (MYIDUP_prod[ip]!=0);
  if ((production==TVar::Lep_WH || production==TVar::Lep_ZH || production==TVar::GammaH) && !(partonIsKnown[2] && partonIsKnown[3])){ if (verbosity>=TVar::INFO) MELAerr << "TUtil::VHiggsMatEl: Final state particles in leptonic/photonic VH have to have a definite id!" << endl; return sum_msqjk; }

  const double allowed_helicities[2] ={ -1, 1 }; // L,R
  // Setup outgoing H decay products (H->f fbar), templated with H->b bar if both fermions are unknown.
  vector<pair<int, int>> Hffparticles;
  double Hffscale=1;
  if (HDKon!=0){
    if (!PDGHelpers::isAnUnknownJet(MYIDUP_dec[0]) || !PDGHelpers::isAnUnknownJet(MYIDUP_dec[1])){ // If one particle is known, pick that line
      if (!PDGHelpers::isAnUnknownJet(MYIDUP_dec[0])) Hffparticles.push_back(pair<int, int>(MYIDUP_dec[0], -MYIDUP_dec[0]));
      else Hffparticles.push_back(pair<int, int>(-MYIDUP_dec[1], MYIDUP_dec[1]));
    }
    else{ // Else loop over possible quark lines
      double Hffscalesum=0;
      for (int hquark=1; hquark<=5; hquark++){
        double Hffmass = __modparameters_MOD_getmass(&hquark);
        Hffscalesum += Hffmass;
        if (hquark==5){
          Hffparticles.push_back(pair<int, int>(hquark, -hquark));
          Hffparticles.push_back(pair<int, int>(-hquark, hquark));
          Hffscale /= Hffmass;
        }
      }
      Hffscale *= Hffscalesum;
    }
  }
  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::VHiggsMatEl: Outgoing H-> f fbar particles to compute for the ME template:" << endl;
    for (unsigned int ihf=0; ihf<Hffparticles.size(); ihf++) MELAout << "\t - (id8, id9) = (" << Hffparticles.at(ihf).first << ", " << Hffparticles.at(ihf).second << ")" << endl;
    MELAout << "TUtil::VHiggsMatEl: ME scale for the H-> f fbar particles: " << Hffscale << endl;
  }

  if (production==TVar::Lep_WH || production==TVar::Had_WH){
    // Setup incoming partons
    vector<pair<int, int>> incomingPartons;
    if (partonIsKnown[0] && partonIsKnown[1]) incomingPartons.push_back(pair<int, int>(MYIDUP_prod[0], MYIDUP_prod[1])); // Parton 0 and 1 are both known
    // FIXME: The following 2/1 assignments assume the CKM element for this pair is non-zero (there are divisions by Vckmsq_in/out later on).
    else if (!partonIsKnown[0] && !partonIsKnown[1]){ // Parton 0 and 1 are unknown
      // Consider all 4 incoming cases: d au, au d, u ad, ad u
      incomingPartons.push_back(pair<int, int>(1, -2)); // du~ -> W-
      incomingPartons.push_back(pair<int, int>(-2, 1)); // u~d -> W-
      incomingPartons.push_back(pair<int, int>(2, -1)); // ud~ -> W+
      incomingPartons.push_back(pair<int, int>(-1, 2)); // d~u -> W+
    }
    else if (!partonIsKnown[1] && partonIsKnown[0]){ // Parton 0 is known
      // Consider the only possible general cases
      if (PDGHelpers::isUpTypeQuark(MYIDUP_prod[0])){ // ud~ or u~d
        for (int iqf=1; iqf<=5; iqf++){
          if (iqf%2==0) continue;
          int jqf = -TMath::Sign(iqf, MYIDUP_prod[0]);
          double ckm_test = __modparameters_MOD_ckmbare(&(MYIDUP_prod[0]), &jqf);
          if (ckm_test!=0.){ incomingPartons.push_back(pair<int, int>(MYIDUP_prod[0], jqf)); break; }
        }
      }
      else if (PDGHelpers::isDownTypeQuark(MYIDUP_prod[0])){ // du~ or d~u
        for (int iqf=2; iqf<=4; iqf++){
          if (iqf%2==1) continue;
          int jqf = -TMath::Sign(iqf, MYIDUP_prod[0]);
          double ckm_test = __modparameters_MOD_ckmbare(&(MYIDUP_prod[0]), &jqf);
          if (ckm_test!=0.){ incomingPartons.push_back(pair<int, int>(MYIDUP_prod[0], jqf)); break; }
        }
      }
    }
    else/* if (!partonIsKnown[0] && partonIsKnown[1])*/{ // Parton 1 is known
      // Consider the only possible general cases
      if (PDGHelpers::isUpTypeQuark(MYIDUP_prod[1])){ // ud~ or u~d
        for (int iqf=1; iqf<=5; iqf++){
          if (iqf%2==0) continue;
          int jqf = -TMath::Sign(iqf, MYIDUP_prod[1]);
          double ckm_test = __modparameters_MOD_ckmbare(&jqf, &(MYIDUP_prod[1]));
          if (ckm_test!=0.){ incomingPartons.push_back(pair<int, int>(jqf, MYIDUP_prod[1])); break; }
        }
      }
      else if (PDGHelpers::isDownTypeQuark(MYIDUP_prod[1])){ // du~ or d~u
        for (int iqf=2; iqf<=4; iqf++){
          if (iqf%2==1) continue;
          int jqf = -TMath::Sign(iqf, MYIDUP_prod[1]);
          double ckm_test = __modparameters_MOD_ckmbare(&jqf, &(MYIDUP_prod[1]));
          if (ckm_test!=0.){ incomingPartons.push_back(pair<int, int>(jqf, MYIDUP_prod[1])); break; }
        }
      }
    }
    if (verbosity>=TVar::DEBUG){
      MELAout << "TUtil::VHiggsMatEl: Incoming partons to compute for the ME template:" << endl;
      for (auto& inpart:incomingPartons) MELAout << "\t - (id1, id2) = (" << inpart.first << ", " << inpart.second << ")" << endl;
    }

    // Setup outgoing partons
    vector<pair<int, int>> outgoingPartons;
    if (partonIsKnown[2] && partonIsKnown[3]) outgoingPartons.push_back(pair<int, int>(MYIDUP_prod[2], MYIDUP_prod[3])); // Parton 0 and 1 are both known or Lep_WH
    // FIXME: The following 2/1 assignments assume the CKM element for this pair is non-zero (there are divisions by Vckmsq_in/out later on).
    else if (!partonIsKnown[2] && !partonIsKnown[3]){ // Parton 0 and 1 are unknown
      // Consider all 4 outgoing cases: d au, au d, u ad, ad u
      outgoingPartons.push_back(pair<int, int>(1, -2)); // W- -> du~
      outgoingPartons.push_back(pair<int, int>(-2, 1)); // W- -> u~d
      outgoingPartons.push_back(pair<int, int>(2, -1)); // W+ -> ud~
      outgoingPartons.push_back(pair<int, int>(-1, 2)); // W+ -> d~u
    }
    else if (!partonIsKnown[3] && partonIsKnown[2]){ // Parton 0 is known
      // Consider the only possible general cases
      if (PDGHelpers::isUpTypeQuark(MYIDUP_prod[2])){ // ud~ or u~d
        for (int iqf=1; iqf<=5; iqf++){
          if (iqf%2==0) continue;
          int jqf = -TMath::Sign(iqf, MYIDUP_prod[2]);
          double ckm_test = __modparameters_MOD_ckmbare(&(MYIDUP_prod[2]), &jqf);
          if (ckm_test!=0.){ outgoingPartons.push_back(pair<int, int>(MYIDUP_prod[2], jqf)); break; }
        }
      }
      else if (PDGHelpers::isDownTypeQuark(MYIDUP_prod[2])){ // du~ or d~u
        for (int iqf=2; iqf<=4; iqf++){
          if (iqf%2==1) continue;
          int jqf = -TMath::Sign(iqf, MYIDUP_prod[2]);
          double ckm_test = __modparameters_MOD_ckmbare(&(MYIDUP_prod[2]), &jqf);
          if (ckm_test!=0.){ outgoingPartons.push_back(pair<int, int>(MYIDUP_prod[2], jqf)); break; }
        }
      }
    }
    else/* if (!partonIsKnown[2] && partonIsKnown[3])*/{ // Parton 1 is known
      // Consider the only possible general cases
      if (PDGHelpers::isUpTypeQuark(MYIDUP_prod[3])){ // ud~ or u~d
        for (int iqf=1; iqf<=5; iqf++){
          if (iqf%2==0) continue;
          int jqf = -TMath::Sign(iqf, MYIDUP_prod[3]);
          double ckm_test = __modparameters_MOD_ckmbare(&jqf, &(MYIDUP_prod[3]));
          if (ckm_test!=0.){ outgoingPartons.push_back(pair<int, int>(jqf, MYIDUP_prod[3])); break; }
        }
      }
      else if (PDGHelpers::isDownTypeQuark(MYIDUP_prod[3])){ // du~ or d~u
        for (int iqf=2; iqf<=4; iqf++){
          if (iqf%2==1) continue;
          int jqf = -TMath::Sign(iqf, MYIDUP_prod[3]);
          double ckm_test = __modparameters_MOD_ckmbare(&jqf, &(MYIDUP_prod[3]));
          if (ckm_test!=0.){ outgoingPartons.push_back(pair<int, int>(jqf, MYIDUP_prod[3])); break; }
        }
      }
    }
    if (verbosity>=TVar::DEBUG){
      MELAout << "TUtil::VHiggsMatEl: Outgoing particles to compute for the ME template:" << endl;
      for (auto& outpart:outgoingPartons) MELAout << "\t - (id6, id7) = (" << outpart.first << ", " << outpart.second << ")" << endl;
    }

    for (auto& inpart:incomingPartons){
      vh_ids[0] = inpart.first;
      vh_ids[1] = inpart.second;
      vh_ids[2] = PDGHelpers::getCoupledVertex(vh_ids[0], vh_ids[1]);
      if (!PDGHelpers::isAWBoson(vh_ids[2])) continue;
      if (verbosity>=TVar::DEBUG) MELAout << "\tIncoming " << vh_ids[0] << "," << vh_ids[1] << " -> " << vh_ids[2] << endl;

      double Vckmsq_in = pow(__modparameters_MOD_ckmbare(&(vh_ids[0]), &(vh_ids[1])), 2);
      if (verbosity>=TVar::DEBUG) MELAout << "\tNeed to divide the ME by |VCKM_incoming|**2 = " << Vckmsq_in << endl;

      for (auto& outpart:outgoingPartons){
        vh_ids[5] = outpart.first;
        vh_ids[6] = outpart.second;
        vh_ids[3] = PDGHelpers::getCoupledVertex(vh_ids[5], vh_ids[6]);
        if (vh_ids[2]!=vh_ids[3]) continue;
        if (verbosity>=TVar::DEBUG) MELAout << "\t\tOutgoing " << vh_ids[3] << " -> " << vh_ids[5] << "," << vh_ids[6] << endl;

        // Compute a raw ME
        double msq=0;
        for (int h01 = 0; h01 < 2; h01++){
          helicities[0] = allowed_helicities[h01];
          helicities[1] = -helicities[0];
          for (int h56 = 0; h56 < 2; h56++){
            helicities[5] = allowed_helicities[h56];
            helicities[6] = -helicities[5];

            double msq_inst=0;
            if (HDKon==0) __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, &msq_inst);
            else{
              for (int h78=0; h78<2; h78++){
                helicities[7]=allowed_helicities[h78];
                helicities[8]=allowed_helicities[h78];
                for (unsigned int ihf=0; ihf<Hffparticles.size(); ihf++){
                  vh_ids[7]=Hffparticles.at(ihf).first;
                  vh_ids[8]=Hffparticles.at(ihf).second;
                  double msq_inst_LR=0;
                  __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, &msq_inst_LR);
                  msq_inst += msq_inst_LR;
                } // End loop over template H->f fbar products
              } // End loop over the spin of H->f fbar line
              msq_inst *= Hffscale;
            } // End HDKon!=0
            msq += msq_inst;
          } // End loop over h56
        } // End loop over h01

        // Determine the outgoing scale
        double scalesum_out=0;
        if (!(partonIsKnown[2] && partonIsKnown[3])){
          double Vckmsq_out = pow(__modparameters_MOD_ckm(&(vh_ids[5]), &(vh_ids[6])), 2);
          if (verbosity>=TVar::DEBUG) MELAout << "\t\tDividing ME by |VCKM_outgoing|**2 = " << Vckmsq_out << endl;
          msq /= Vckmsq_out;
          for (int outgoing1=1; outgoing1<=nf; outgoing1++){
            if (partonIsKnown[2] && outgoing1!=abs(vh_ids[5])) continue;
            if (outgoing1%2!=abs(vh_ids[5])%2 || outgoing1==6) continue;
            int iout = outgoing1 * TMath::Sign(1, vh_ids[5]);
            for (int outgoing2=1; outgoing2<=nf; outgoing2++){
              if (partonIsKnown[3] && outgoing2!=abs(vh_ids[6])) continue;
              if (outgoing2%2!=abs(vh_ids[6])%2 || outgoing2==6) continue;
              int jout = outgoing2 * TMath::Sign(1, vh_ids[6]);
              scalesum_out += pow(__modparameters_MOD_ckm(&(iout), &(jout)), 2);
            }
          }
        }
        else scalesum_out = 1;
        if (verbosity>=TVar::DEBUG) MELAout << "\t\tScale for outgoing particles: " << scalesum_out << endl;

        // Divide ME by the incoming scale factor (will be multiplied again inside the loop)
        if (!partonIsKnown[0] || !partonIsKnown[1]){
          if (verbosity>=TVar::DEBUG) MELAout << "\t\tDividing ME by |VCKM_incoming|**2 = " << Vckmsq_in << endl;
          msq /= Vckmsq_in;
        }

        // Sum all possible combinations
        for (int incoming1=1; incoming1<=nf; incoming1++){
          if (partonIsKnown[0] && incoming1!=abs(vh_ids[0])) continue;
          if (incoming1%2!=abs(vh_ids[0])%2 || incoming1==6) continue;
          int iin = incoming1 * TMath::Sign(1, vh_ids[0]);
          for (int incoming2=1; incoming2<=nf; incoming2++){
            if (partonIsKnown[1] && incoming2!=abs(vh_ids[1])) continue;
            if (incoming2%2!=abs(vh_ids[1])%2 || incoming2==6) continue;
            int jin = incoming2 * TMath::Sign(1, vh_ids[1]);
            if (verbosity>=TVar::DEBUG) MELAout << "\t\t\tiin,jin = " << iin << "," << jin << endl;
            double scale_in=1;
            if (!partonIsKnown[0] || !partonIsKnown[1]){
              scale_in = pow(__modparameters_MOD_ckmbare(&(iin), &(jin)), 2);
              if (verbosity>=TVar::DEBUG) MELAout << "\t\tScale for incoming particles: " << scale_in << endl;
            }
            MatElsq[jin+5][iin+5] += msq * 0.25 * scale_in *scalesum_out;
          }
        }
        // msq is now added to MatElSq.
      } // End loop over outgoing particle templates
    } // End loop over incoming parton templates
  } // End WH
  else{ // ZH, GammaH
    // Setup incoming partons
    vector<pair<int, int>> incomingPartons;
    if (partonIsKnown[0] && partonIsKnown[1]) incomingPartons.push_back(pair<int, int>(MYIDUP_prod[0], MYIDUP_prod[1])); // Parton 0 and 1 are both known
    else if (!partonIsKnown[0] && !partonIsKnown[1]){ // Parton 0 and 1 are unknown
      // Consider all 4 incoming cases: d ad, ad d, u au, au u
      incomingPartons.push_back(pair<int, int>(1, -1)); // dd~ -> Z
      incomingPartons.push_back(pair<int, int>(-1, 1)); // d~d -> Z
      incomingPartons.push_back(pair<int, int>(2, -2)); // uu~ -> Z
      incomingPartons.push_back(pair<int, int>(-2, 2)); // u~u -> Z
    }
    else if (!partonIsKnown[1] && partonIsKnown[0]) // Parton 0 is known
      incomingPartons.push_back(pair<int, int>(MYIDUP_prod[0], -MYIDUP_prod[0])); // id1, -id1
    else/* if (!partonIsKnown[0] && partonIsKnown[1])*/ // Parton 1 is known
      incomingPartons.push_back(pair<int, int>(-MYIDUP_prod[1], MYIDUP_prod[1])); // -id2, id2
    if (verbosity>=TVar::DEBUG){
      MELAout << "TUtil::VHiggsMatEl: Incoming partons to compute for the ME template:" << endl;
      for (auto& inpart:incomingPartons) MELAout << "\t - (id1, id2) = (" << inpart.first << ", " << inpart.second << ")" << endl;
    }

    // Setup outgoing partons
    vector<pair<int, int>> outgoingPartons;
    if ((partonIsKnown[2] && partonIsKnown[3]) || production==TVar::GammaH) outgoingPartons.push_back(pair<int, int>(MYIDUP_prod[2], MYIDUP_prod[3])); // Parton 0 and 1 are both known or Lep_ZH
    else if (!partonIsKnown[2] && !partonIsKnown[3]){ // Parton 0 and 1 are unknown
      // Consider all 4 outgoing cases: d au, au d, u ad, ad u
      outgoingPartons.push_back(pair<int, int>(1, -1)); // Z -> dd~
      outgoingPartons.push_back(pair<int, int>(-1, 1)); // Z -> d~d
      outgoingPartons.push_back(pair<int, int>(2, -2)); // Z -> uu~
      outgoingPartons.push_back(pair<int, int>(-2, 2)); // Z -> u~u
    }
    else if (!partonIsKnown[3] && partonIsKnown[2]) // Parton 0 is known
      outgoingPartons.push_back(pair<int, int>(MYIDUP_prod[2], -MYIDUP_prod[2])); // id1, -id1
    else/* if (!partonIsKnown[2] && partonIsKnown[3])*/ // Parton 1 is known
      outgoingPartons.push_back(pair<int, int>(-MYIDUP_prod[3], MYIDUP_prod[3])); // -id2, id2
    if (verbosity>=TVar::DEBUG){
      MELAout << "TUtil::VHiggsMatEl: Outgoing particles to compute for the ME template:" << endl;
      for (unsigned int op=0; op<outgoingPartons.size(); op++) MELAout << "\t - (id6, id7) = (" << outgoingPartons.at(op).first << ", " << outgoingPartons.at(op).second << ")" << endl;
    }

    for (auto& inpart:incomingPartons){
      vh_ids[0] = inpart.first;
      vh_ids[1] = inpart.second;
      vh_ids[2] = PDGHelpers::getCoupledVertex(vh_ids[0], vh_ids[1]);
      if (!PDGHelpers::isAZBoson(vh_ids[2])) continue; // Notice, Z-> Gamma + H should also have the id of the Z!
      if (verbosity>=TVar::DEBUG) MELAout << "\tIncoming " << vh_ids[0] << "," << vh_ids[1] << " -> " << vh_ids[2] << endl;

      // Scale for the incoming state is 1

      for (auto& outpart:outgoingPartons){
        vh_ids[5] = outpart.first;
        vh_ids[6] = outpart.second;
        if (production==TVar::GammaH) vh_ids[3] = 22;
        else{
          vh_ids[3] = PDGHelpers::getCoupledVertex(vh_ids[5], vh_ids[6]);
          if (vh_ids[2]!=vh_ids[3]) continue;
        }
        if (verbosity>=TVar::DEBUG) MELAout << "\t\tOutgoing " << vh_ids[3] << " -> " << vh_ids[5] << "," << vh_ids[6] << endl;

        // Compute a raw ME
        double msq=0;
        for (int h01 = 0; h01 < 2; h01++){
          helicities[0] = allowed_helicities[h01];
          helicities[1] = -helicities[0];
          for (int h56 = 0; h56 < 2; h56++){
            helicities[5] = allowed_helicities[h56];
            helicities[6] = -helicities[5];

            double msq_inst=0;
            if (HDKon==0) __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, &msq_inst);
            else{
              for (int h78=0; h78<2; h78++){
                helicities[7]=allowed_helicities[h78];
                helicities[8]=allowed_helicities[h78];
                for (unsigned int ihf=0; ihf<Hffparticles.size(); ihf++){
                  vh_ids[7]=Hffparticles.at(ihf).first;
                  vh_ids[8]=Hffparticles.at(ihf).second;
                  double msq_inst_LR=0;
                  __modvhiggs_MOD_evalamp_vhiggs(vh_ids, helicities, p4, &msq_inst_LR);
                  msq_inst += msq_inst_LR;
                } // End loop over template H->f fbar products
              } // End loop over the spin of H->f fbar line
              msq_inst *= Hffscale;
            } // End HDKon!=0
            msq += msq_inst;
          } // End loop over h56
        } // End loop over h01

        // Determine the outgoing scale
        double scale_out=1;
        if (!partonIsKnown[2] && !partonIsKnown[3] && production!=TVar::GammaH){
          if (PDGHelpers::isDownTypeQuark(vh_ids[5])) scale_out=3;
          else if (PDGHelpers::isUpTypeQuark(vh_ids[5])) scale_out=2;
        }
        if (verbosity>=TVar::DEBUG) MELAout << "\t\tScale for outgoing particles: " << scale_out << endl;

        // Sum all possible combinations
        for (int incoming1=1; incoming1<=nf; incoming1++){
          if (partonIsKnown[0] && incoming1!=abs(vh_ids[0])) continue;
          if (incoming1%2!=abs(vh_ids[0])%2 || incoming1==6) continue;
          int iin = incoming1 * TMath::Sign(1, vh_ids[0]);
          int jin=-iin;
          if (verbosity>=TVar::DEBUG) MELAout << "\t\t\tiin,jin = " << iin << "," << jin << endl;
          MatElsq[jin+5][iin+5] += msq * 0.25 * scale_out; // No incoming state scale
        }
        // msq is now added to MatElSq.
      } // End loop over outgoing particle templates
    } // End loop over incoming parton templates
  } // End ZH, GammaH

  int GeVexponent_MEsq;
  if (HDKon==0) GeVexponent_MEsq = 4-(1+nRequested_AssociatedJets+nRequested_AssociatedLeptons+nRequested_AssociatedPhotons)*2-4; // Amplitude has additional 1/m**2 from propagator == MEsq has additional 1/m**4
  else GeVexponent_MEsq = 4-(2+nRequested_AssociatedJets+nRequested_AssociatedLeptons+nRequested_AssociatedPhotons)*2;
  double constant = pow(GeV, -GeVexponent_MEsq);
  for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) MatElsq[jj][ii] *= constant; }
  if (verbosity >= TVar::DEBUG){
    MELAout << "MatElsq:\n";
    for (int ii = 0; ii < nmsq; ii++){ for (int jj = 0; jj < nmsq; jj++) MELAout << MatElsq[jj][ii] << '\t'; MELAout << endl; }
  }
  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, EBEAM, verbosity);

  // Turn H_DK off
  if (HDKon!=0){
    HDKon=0;
    __modjhugenmela_MOD_sethdk(&HDKon);
  }

  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::VHiggsMatEl: Reset AlphaS:\n"
      << "\tBefore reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << endl;
  }
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel);
  if (verbosity>=TVar::DEBUG){
    GetAlphaS(&alphasVal, &alphasmzVal);
    MELAout
      << "TUtil::VHiggsMatEl: Reset AlphaS result:\n"
      << "\tAfter reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }
  return sum_msqjk;
}

// ttH, H undecayed
double TUtil::TTHiggsMatEl(
  const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  const double& EBEAM,
  int topDecay, int topProcess,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double sum_msqjk = 0;
  double MatElsq[nmsq][nmsq]={ { 0 } };

  if (matrixElement!=TVar::JHUGen){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::TTHiggsMatEl: Non-JHUGen MEs are not supported." << endl; return sum_msqjk; }
  if (production!=TVar::ttH){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::TTHiggsMatEl: Only ttH is supported." << endl; return sum_msqjk; }

  int partIncCode;
  int nRequested_Tops=1;
  int nRequested_Antitops=1;
  if (topDecay>0) partIncCode=TVar::kUseAssociated_UnstableTops; // Look for unstable tops
  else partIncCode=TVar::kUseAssociated_StableTops; // Look for unstable tops

  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.nRequested_Tops=nRequested_Tops;
  mela_event.nRequested_Antitops=nRequested_Antitops;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event,
    verbosity
    );


  if (topDecay==0 && (mela_event.pStableTops.size()<1 || mela_event.pStableAntitops.size()<1)){
    if (verbosity>=TVar::ERROR) MELAerr
      << "TUtil::TTHiggsMatEl: Number of stable tops (" << mela_event.pStableTops.size() << ")"
      << " and number of stable antitops (" << mela_event.pStableAntitops.size() << ")"
      << " in ttH process are not 1!" << endl;
    return sum_msqjk;
  }
  else if (topDecay>0 && (mela_event.pTopDaughters.size()<1 || mela_event.pAntitopDaughters.size()<1)){
    if (verbosity>=TVar::ERROR) MELAerr
      << "TUtil::TTHiggsMatEl: Number of set of top daughters (" << mela_event.pTopDaughters.size() << ")"
      << " and number of set of antitop daughters (" << mela_event.pAntitopDaughters.size() << ")"
      << " in ttH process are not 1!" << endl;
    return sum_msqjk;
  }
  else if (topDecay>0 && (mela_event.pTopDaughters.at(0).size()!=3 || mela_event.pAntitopDaughters.at(0).size()!=3)){
    if (verbosity>=TVar::ERROR) MELAerr
      << "TUtil::TTHiggsMatEl: Number of top daughters (" << mela_event.pTopDaughters.at(0).size() << ")"
      << " and number of antitop daughters (" << mela_event.pAntitopDaughters.at(0).size() << ")"
      << " in ttH process are not 3!" << endl;
    return sum_msqjk;
  }

  SimpleParticleCollection_t topDaughters;
  SimpleParticleCollection_t antitopDaughters;
  bool isUnknown[2]={ true, true };

  if (topDecay>0){
    // Daughters are assumed to have been ordered as b, Wf, Wfb already.
    for (unsigned int itd=0; itd<mela_event.pTopDaughters.at(0).size(); itd++) topDaughters.push_back(mela_event.pTopDaughters.at(0).at(itd));
    for (unsigned int itd=0; itd<mela_event.pAntitopDaughters.at(0).size(); itd++) antitopDaughters.push_back(mela_event.pAntitopDaughters.at(0).at(itd));
  }
  else{
    for (unsigned int itop=0; itop<mela_event.pStableTops.size(); itop++) topDaughters.push_back(mela_event.pStableTops.at(itop));
    for (unsigned int itop=0; itop<mela_event.pStableAntitops.size(); itop++) antitopDaughters.push_back(mela_event.pStableAntitops.at(itop));
  }
  // Check if either top is definitely identified
  for (unsigned int itd=0; itd<topDaughters.size(); itd++){ if (!PDGHelpers::isAnUnknownJet(topDaughters.at(itd).first)){ isUnknown[0]=false; break; } }
  for (unsigned int itd=0; itd<antitopDaughters.size(); itd++){ if (!PDGHelpers::isAnUnknownJet(antitopDaughters.at(itd).first)){ isUnknown[1]=false; break; } }

  // Start assigning the momenta
  // 0,1: p1 p2
  // 2-4: H,tb,t
  // 5-8: bb,W-,f,fb
  // 9-12: b,W+,fb,f
  double p4[13][4]={ { 0 } };
  double MYIDUP_prod[2]={ 0 };
  TLorentzVector MomStore[mxpart];
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);
  for (int ipar=0; ipar<2; ipar++){
    TLorentzVector* momTmp = &(mela_event.pMothers.at(ipar).second);
    int* idtmp = &(mela_event.pMothers.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_prod[ipar] = *idtmp;
    else MYIDUP_prod[ipar] = 0;
    if (momTmp->T()>0.){
      p4[ipar][0] = -momTmp->T()*GeV;
      p4[ipar][1] = -momTmp->X()*GeV;
      p4[ipar][2] = -momTmp->Y()*GeV;
      p4[ipar][3] = -momTmp->Z()*GeV;
      MomStore[ipar] = (*momTmp);
    }
    else{
      p4[ipar][0] = momTmp->T()*GeV;
      p4[ipar][1] = momTmp->X()*GeV;
      p4[ipar][2] = momTmp->Y()*GeV;
      p4[ipar][3] = momTmp->Z()*GeV;
      MomStore[ipar] = -(*momTmp);
      MYIDUP_prod[ipar] = -MYIDUP_prod[ipar];
    }
  }

  // Assign top momenta
  // t(4) -> b(9) W+(10) (-> f(12) fb(11))
  const unsigned int t_pos=4;
  const unsigned int b_pos=9;
  const unsigned int Wp_pos=10;
  const unsigned int Wpf_pos=12;
  const unsigned int Wpfb_pos=11;
  for (unsigned int ipar=0; ipar<topDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(topDaughters.at(ipar).second);
    if (topDaughters.size()==1){
      p4[t_pos][0] = momTmp->T()*GeV;
      p4[t_pos][1] = momTmp->X()*GeV;
      p4[t_pos][2] = momTmp->Y()*GeV;
      p4[t_pos][3] = momTmp->Z()*GeV;
    }
    // size==3
    else if (ipar==0){ // b
      p4[b_pos][0] = momTmp->T()*GeV;
      p4[b_pos][1] = momTmp->X()*GeV;
      p4[b_pos][2] = momTmp->Y()*GeV;
      p4[b_pos][3] = momTmp->Z()*GeV;
    }
    else if (ipar==2){ // Wfb
      p4[Wpfb_pos][0] = momTmp->T()*GeV;
      p4[Wpfb_pos][1] = momTmp->X()*GeV;
      p4[Wpfb_pos][2] = momTmp->Y()*GeV;
      p4[Wpfb_pos][3] = momTmp->Z()*GeV;
      p4[Wp_pos][0] += p4[Wpfb_pos][0];
      p4[Wp_pos][1] += p4[Wpfb_pos][1];
      p4[Wp_pos][2] += p4[Wpfb_pos][2];
      p4[Wp_pos][3] += p4[Wpfb_pos][3];
    }
    else/* if (ipar==1)*/{ // Wf
      p4[Wpf_pos][0] = momTmp->T()*GeV;
      p4[Wpf_pos][1] = momTmp->X()*GeV;
      p4[Wpf_pos][2] = momTmp->Y()*GeV;
      p4[Wpf_pos][3] = momTmp->Z()*GeV;
      p4[Wp_pos][0] += p4[Wpf_pos][0];
      p4[Wp_pos][1] += p4[Wpf_pos][1];
      p4[Wp_pos][2] += p4[Wpf_pos][2];
      p4[Wp_pos][3] += p4[Wpf_pos][3];
    }
    MomStore[6] = MomStore[6] + (*momTmp); // MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }
  if (topDaughters.size()!=1){ for (unsigned int ix=0; ix<4; ix++){ for (unsigned int ip=b_pos; ip<=Wp_pos; ip++) p4[t_pos][ix] = p4[ip][ix]; } }

  // Assign antitop momenta
  // tb(3) -> bb(5) W-(6) (-> f(7) fb(8))
  const unsigned int tb_pos=3;
  const unsigned int bb_pos=5;
  const unsigned int Wm_pos=6;
  const unsigned int Wmf_pos=7;
  const unsigned int Wmfb_pos=8;
  for (unsigned int ipar=0; ipar<antitopDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(antitopDaughters.at(ipar).second);
    if (antitopDaughters.size()==1){
      p4[tb_pos][0] = momTmp->T()*GeV;
      p4[tb_pos][1] = momTmp->X()*GeV;
      p4[tb_pos][2] = momTmp->Y()*GeV;
      p4[tb_pos][3] = momTmp->Z()*GeV;
    }
    // size==3
    else if (ipar==0){ // bb
      p4[bb_pos][0] = momTmp->T()*GeV;
      p4[bb_pos][1] = momTmp->X()*GeV;
      p4[bb_pos][2] = momTmp->Y()*GeV;
      p4[bb_pos][3] = momTmp->Z()*GeV;
    }
    else if (ipar==1){ // Wf
      p4[Wmf_pos][0] = momTmp->T()*GeV;
      p4[Wmf_pos][1] = momTmp->X()*GeV;
      p4[Wmf_pos][2] = momTmp->Y()*GeV;
      p4[Wmf_pos][3] = momTmp->Z()*GeV;
      p4[Wm_pos][0] += p4[Wmf_pos][0];
      p4[Wm_pos][1] += p4[Wmf_pos][1];
      p4[Wm_pos][2] += p4[Wmf_pos][2];
      p4[Wm_pos][3] += p4[Wmf_pos][3];
    }
    else/* if (ipar==1)*/{ // Wfb
      p4[Wmfb_pos][0] = momTmp->T()*GeV;
      p4[Wmfb_pos][1] = momTmp->X()*GeV;
      p4[Wmfb_pos][2] = momTmp->Y()*GeV;
      p4[Wmfb_pos][3] = momTmp->Z()*GeV;
      p4[Wm_pos][0] += p4[Wmfb_pos][0];
      p4[Wm_pos][1] += p4[Wmfb_pos][1];
      p4[Wm_pos][2] += p4[Wmfb_pos][2];
      p4[Wm_pos][3] += p4[Wmfb_pos][3];
    }
    MomStore[7] = MomStore[7] + (*momTmp); // MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }
  if (antitopDaughters.size()!=1){ for (unsigned int ix=0; ix<4; ix++){ for (unsigned int ip=5; ip<=6; ip++) p4[tb_pos][ix] = p4[ip][ix]; } }

  for (unsigned int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
    p4[2][0] += momTmp->T()*GeV;
    p4[2][1] += momTmp->X()*GeV;
    p4[2][2] += momTmp->Y()*GeV;
    p4[2][3] += momTmp->Z()*GeV;
    MomStore[5] = MomStore[5] + (*momTmp); // i==(2, 3, 4) is (J1, J2, H), recorded as MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }

  if (verbosity >= TVar::DEBUG){
    for (int ii=0; ii<13; ii++){ MELAout << "p4[" << ii << "] = "; for (int jj=0; jj<4; jj++) MELAout << p4[ii][jj]/GeV << '\t'; MELAout << endl; }
  }

  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l");
  double alphasVal, alphasmzVal;
  GetAlphaS(&alphasVal, &alphasmzVal);
  RcdME->setRenormalizationScale(renQ);
  RcdME->setFactorizationScale(facQ);
  RcdME->setAlphaS(alphasVal);
  RcdME->setAlphaSatMZ(alphasmzVal);
  RcdME->setHiggsMassWidth(masses_mcfm_.hmass, masses_mcfm_.hwidth, 0);
  RcdME->setHiggsMassWidth(spinzerohiggs_anomcoupl_.h2mass, spinzerohiggs_anomcoupl_.h2width, 1);
  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::TTHiggsMatEl: Set AlphaS:\n"
      << "\tBefore set, alphas scale: " << defaultRenScale << ", PDF scale: " << defaultFacScale << '\n'
      << "\trenQ: " << renQ << " ( x " << event_scales->ren_scale_factor << "), facQ: " << facQ << " ( x " << event_scales->fac_scale_factor << ")\n"
      << "\tAfter set, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }

  __modjhugenmela_MOD_settopdecays(&topDecay);
  /***** BEGIN TTH ME CALCULATION *****/
  for (unsigned int ib1=0; ib1<3; ib1++){
    if (topDaughters.size()==1 && ib1!=0) continue;
    else if (topDaughters.size()==3 && !PDGHelpers::isAnUnknownJet(topDaughters.at(0).first) && ib1!=0) continue;
    unsigned int b1index;
    if (ib1==0) b1index = b_pos;
    else if (ib1==1) b1index = Wpf_pos;
    else b1index = Wpfb_pos;
    for (unsigned int if1=0; if1<3; if1++){
      if (topDaughters.size()==1 && if1!=0) continue;
      else if (topDaughters.size()==3 && !PDGHelpers::isAnUnknownJet(topDaughters.at(1).first) && if1!=0) continue;
      unsigned int f1index;
      if (if1==0) f1index = Wpf_pos;
      else if (if1==2) f1index = b_pos;
      else f1index = Wpfb_pos;
      for (unsigned int ifb1=0; ifb1<3; ifb1++){
        if (topDaughters.size()==1 && ifb1!=0) continue;
        else if (topDaughters.size()==3 && !PDGHelpers::isAnUnknownJet(topDaughters.at(2).first) && ifb1!=0) continue;
        unsigned int fb1index;
        if (ifb1==2) fb1index = Wpf_pos;
        else if (ifb1==1) fb1index = b_pos;
        else fb1index = Wpfb_pos;

        if (b1index==f1index || b1index==fb1index || f1index==fb1index) continue;

        for (unsigned int ib2=0; ib2<3; ib2++){
          if (antitopDaughters.size()==1 && ib2!=0) continue;
          else if (antitopDaughters.size()==3 && !PDGHelpers::isAnUnknownJet(antitopDaughters.at(0).first) && ib2!=0) continue;
          unsigned int b2index;
          if (ib2==0) b2index = bb_pos;
          else if (ib2==1) b2index = Wmf_pos;
          else b2index = Wmfb_pos;
          for (unsigned int if2=0; if2<3; if2++){
            if (antitopDaughters.size()==1 && if2!=0) continue;
            else if (antitopDaughters.size()==3 && !PDGHelpers::isAnUnknownJet(antitopDaughters.at(1).first) && if2!=0) continue;
            unsigned int f2index;
            if (if2==0) f2index = Wmf_pos;
            else if (if2==2) f2index = bb_pos;
            else f2index = Wmfb_pos;
            for (unsigned int ifb2=0; ifb2<3; ifb2++){
              if (antitopDaughters.size()==1 && ifb2!=0) continue;
              else if (antitopDaughters.size()==3 && !PDGHelpers::isAnUnknownJet(antitopDaughters.at(2).first) && ifb2!=0) continue;
              unsigned int fb2index;
              if (ifb2==2) fb2index = Wmf_pos;
              else if (ifb2==1) fb2index = bb_pos;
              else fb2index = Wmfb_pos;

              if (b2index==f2index || b2index==fb2index || f2index==fb2index) continue;

              double p4_current[13][4]={ { 0 } };
              for (unsigned int ix=0; ix<4; ix++){
                for (unsigned int ip=0; ip<=2; ip++) p4_current[ip][ix] = p4[ip][ix]; // I1, I2, H do not change.
                p4_current[t_pos][ix] = p4[t_pos][ix]; // t does not change.
                p4_current[b1index][ix] = p4[b_pos][ix]; // Assign b to different position.
                p4_current[f1index][ix] = p4[Wpf_pos][ix]; // Assign Wp->f? to different position.
                p4_current[fb1index][ix] = p4[Wpfb_pos][ix]; // Assign Wp->?fb to different position.
                p4_current[Wp_pos][ix] = p4_current[Wpf_pos][ix] + p4_current[Wpfb_pos][ix]; // Re-sum W+ momentum.

                p4_current[tb_pos][ix] = p4[tb_pos][ix]; // tb does not change.
                p4_current[b2index][ix] = p4[bb_pos][ix]; // Assign bb to different position.
                p4_current[f2index][ix] = p4[Wmf_pos][ix]; // Assign Wm->f? to different position.
                p4_current[fb2index][ix] = p4[Wmfb_pos][ix]; // Assign Wm->?fb to different position.
                p4_current[Wm_pos][ix] = p4_current[Wmf_pos][ix] + p4_current[Wmfb_pos][ix]; // Re-sum W- momentum.
              }
              if (verbosity>=TVar::DEBUG){
                MELAout
                  << "TUtil::TTHiggsMatEl: Unswapped instance for "
                  << "b(" << b_pos << ") -> " << b1index << ", "
                  << "Wpf(" << Wpf_pos << ") -> " << f1index << ", "
                  << "Wpfb(" << Wpfb_pos << ") -> " << fb1index << ", "
                  << "bb(" << bb_pos << ") -> " << b2index << ", "
                  << "Wmf(" << Wmf_pos << ") -> " << f2index << ", "
                  << "Wmfb(" << Wmfb_pos << ") -> " << fb2index << endl;
                for (int ii=0; ii<13; ii++){ MELAout << "p4_instance[" << ii << "] = "; for (int jj=0; jj<4; jj++) MELAout << p4_current[ii][jj]/GeV << '\t'; MELAout << endl; }
                MELAout << endl;
              }

              double MatElsq_tmp[nmsq][nmsq]={ { 0 } };
              double MatElsq_tmp_swap[nmsq][nmsq]={ { 0 } };
              __modttbhiggs_MOD_evalxsec_pp_ttbh(p4_current, &topProcess, MatElsq_tmp);
              if (isUnknown[0] && isUnknown[1]){
                for (unsigned int ix=0; ix<4; ix++){
                  swap(p4_current[t_pos][ix], p4_current[tb_pos][ix]);
                  swap(p4_current[b_pos][ix], p4_current[bb_pos][ix]);
                  swap(p4_current[Wp_pos][ix], p4_current[Wm_pos][ix]);
                  swap(p4_current[Wpf_pos][ix], p4_current[Wmf_pos][ix]);
                  swap(p4_current[Wpfb_pos][ix], p4_current[Wmfb_pos][ix]);
                }
                __modttbhiggs_MOD_evalxsec_pp_ttbh(p4_current, &topProcess, MatElsq_tmp_swap);
                for (int ix=0; ix<11; ix++){ for (int iy=0; iy<11; iy++) MatElsq_tmp[iy][ix] = (MatElsq_tmp[iy][ix]+MatElsq_tmp_swap[iy][ix])/2.; }
              }
              for (int ix=0; ix<11; ix++){ for (int iy=0; iy<11; iy++) MatElsq[iy][ix] += MatElsq_tmp[iy][ix]; }
              if (verbosity>=TVar::DEBUG){
                MELAout
                  << "TUtil::TTHiggsMatEl: Swapped instance for "
                  << "b(" << b_pos << ") -> " << b1index << ", "
                  << "Wpf(" << Wpf_pos << ") -> " << f1index << ", "
                  << "Wpfb(" << Wpfb_pos << ") -> " << fb1index << ", "
                  << "bb(" << bb_pos << ") -> " << b2index << ", "
                  << "Wmf(" << Wmf_pos << ") -> " << f2index << ", "
                  << "Wmfb(" << Wmfb_pos << ") -> " << fb2index << endl;
                for (int ii=0; ii<13; ii++){ MELAout << "p4_instance[" << ii << "] = "; for (int jj=0; jj<4; jj++) MELAout << p4_current[ii][jj]/GeV << '\t'; MELAout << endl; }
                MELAout << endl;
              }
            } // End loop over ifb2
          } // End loop over if2
        } // End loop over ib2
      } // End loop over ifb1
    } // End loop over if1
  } // End loop over ib1
  /***** END TTH ME CALCULATION *****/
  int defaultTopDecay=-1;
  __modjhugenmela_MOD_settopdecays(&defaultTopDecay); // reset top decay

  int GeVexponent_MEsq;
  if (topDecay>0) GeVexponent_MEsq = 4-(1+3*(nRequested_Tops+nRequested_Antitops))*2;
  else GeVexponent_MEsq = 4-(1+nRequested_Tops+nRequested_Antitops)*2;
  double constant = pow(GeV, -GeVexponent_MEsq);
  for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) MatElsq[jj][ii] *= constant; }
  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::TTHiggsMatEl: MEsq[ip][jp] = " << endl;
    for (int iquark=-5; iquark<=5; iquark++){
      for (int jquark=-5; jquark<=5; jquark++) MELAout << MatElsq[jquark+5][iquark+5] << '\t';
      MELAout << endl;
    }
  }
  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, EBEAM, verbosity);

  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::TTHiggsMatEl: Reset AlphaS:\n"
      << "\tBefore reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << endl;
  }
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel);
  if (verbosity>=TVar::DEBUG){
    GetAlphaS(&alphasVal, &alphasmzVal);
    MELAout
      << "TUtil::TTHiggsMatEl: Reset AlphaS result:\n"
      << "\tAfter reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }
  return sum_msqjk;
}

// bbH, H undecayed
double TUtil::BBHiggsMatEl(
  const TVar::Process& process, const TVar::Production& production, const TVar::MatrixElement& matrixElement,
  event_scales_type* event_scales, MelaIO* RcdME,
  const double& EBEAM,
  int botProcess,
  TVar::VerbosityLevel verbosity
  ){
  const double GeV=1./100.; // JHUGen mom. scale factor
  double sum_msqjk = 0;
  double MatElsq[nmsq][nmsq]={ { 0 } };
  double MatElsq_tmp[nmsq][nmsq]={ { 0 } };

  if (matrixElement!=TVar::JHUGen){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::BBHiggsMatEl: Non-JHUGen MEs are not supported." << endl; return sum_msqjk; }
  if (production!=TVar::bbH){ if (verbosity>=TVar::ERROR) MELAerr << "TUtil::BBHiggsMatEl: Only bbH is supported." << endl; return sum_msqjk; }

  int partIncCode=TVar::kUseAssociated_Jets; // Look for jets
  int nRequested_AssociatedJets=2;

  simple_event_record mela_event;
  mela_event.AssociationCode=partIncCode;
  mela_event.nRequested_AssociatedJets=nRequested_AssociatedJets;
  GetBoostedParticleVectors(
    RcdME->melaCand,
    mela_event,
    verbosity
    );

  if (mela_event.pAssociated.size()<2){
    if (verbosity>=TVar::ERROR) MELAerr
      << "TUtil::BBHiggsMatEl: Number of stable bs (" << mela_event.pAssociated.size() << ")"
      <<" in bbH process is not 2!" << endl;
    return sum_msqjk;
  }

  SimpleParticleCollection_t topDaughters;
  SimpleParticleCollection_t antitopDaughters;
  bool isUnknown[2]; isUnknown[0]=false; isUnknown[1]=false;

  if (mela_event.pAssociated.at(0).first>=0){ topDaughters.push_back(mela_event.pAssociated.at(0)); isUnknown[0]=(PDGHelpers::isAnUnknownJet(mela_event.pAssociated.at(0).first)); }
  else antitopDaughters.push_back(mela_event.pAssociated.at(0));
  if (mela_event.pAssociated.at(1).first<=0){ antitopDaughters.push_back(mela_event.pAssociated.at(1)); isUnknown[1]=(PDGHelpers::isAnUnknownJet(mela_event.pAssociated.at(1).first)); }
  else topDaughters.push_back(mela_event.pAssociated.at(1));

  // Start assigning the momenta
  // 0,1: p1 p2
  // 2-4: H,tb,t
  // 5-8: bb,W-,f,fb
  // 9-12: b,W+,fb,f
  double p4[13][4]={ { 0 } };
  double MYIDUP_prod[2]={ 0 };
  TLorentzVector MomStore[mxpart];
  for (int i = 0; i < mxpart; i++) MomStore[i].SetXYZT(0, 0, 0, 0);
  for (int ipar=0; ipar<2; ipar++){
    TLorentzVector* momTmp = &(mela_event.pMothers.at(ipar).second);
    int* idtmp = &(mela_event.pMothers.at(ipar).first);
    if (!PDGHelpers::isAnUnknownJet(*idtmp)) MYIDUP_prod[ipar] = *idtmp;
    else MYIDUP_prod[ipar] = 0;
    if (momTmp->T()>0.){
      p4[ipar][0] = -momTmp->T()*GeV;
      p4[ipar][1] = -momTmp->X()*GeV;
      p4[ipar][2] = -momTmp->Y()*GeV;
      p4[ipar][3] = -momTmp->Z()*GeV;
      MomStore[ipar] = (*momTmp);
    }
    else{
      p4[ipar][0] = momTmp->T()*GeV;
      p4[ipar][1] = momTmp->X()*GeV;
      p4[ipar][2] = momTmp->Y()*GeV;
      p4[ipar][3] = momTmp->Z()*GeV;
      MomStore[ipar] = -(*momTmp);
      MYIDUP_prod[ipar] = -MYIDUP_prod[ipar];
    }
  }

  // Assign b momenta
  for (unsigned int ipar=0; ipar<topDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(topDaughters.at(ipar).second);
    p4[4][0] = momTmp->T()*GeV;
    p4[4][1] = momTmp->X()*GeV;
    p4[4][2] = momTmp->Y()*GeV;
    p4[4][3] = momTmp->Z()*GeV;
    MomStore[6] = MomStore[6] + (*momTmp); // MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }

  // Assign bb momenta
  for (unsigned int ipar=0; ipar<antitopDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(antitopDaughters.at(ipar).second);
    p4[3][0] = momTmp->T()*GeV;
    p4[3][1] = momTmp->X()*GeV;
    p4[3][2] = momTmp->Y()*GeV;
    p4[3][3] = momTmp->Z()*GeV;
    MomStore[7] = MomStore[7] + (*momTmp); // MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }

  for (unsigned int ipar=0; ipar<mela_event.pDaughters.size(); ipar++){
    TLorentzVector* momTmp = &(mela_event.pDaughters.at(ipar).second);
    p4[2][0] += momTmp->T()*GeV;
    p4[2][1] += momTmp->X()*GeV;
    p4[2][2] += momTmp->Y()*GeV;
    p4[2][3] += momTmp->Z()*GeV;
    MomStore[5] = MomStore[5] + (*momTmp); // i==(2, 3, 4) is (J1, J2, H), recorded as MomStore (I1, I2, 0, 0, 0, H, J1, J2)
  }

  if (verbosity >= TVar::DEBUG){
    for (int ii=0; ii<13; ii++){ MELAout << "p4[" << ii << "] = "; for (int jj=0; jj<4; jj++) MELAout << p4[ii][jj]/GeV << '\t'; MELAout << endl; }
  }

  double defaultRenScale = scale_.scale;
  double defaultFacScale = facscale_.facscale;
  //MELAout << "Default scales: " << defaultRenScale << '\t' << defaultFacScale << endl;
  int defaultNloop = nlooprun_.nlooprun;
  int defaultNflav = nflav_.nflav;
  string defaultPdflabel = pdlabel_.pdlabel;
  double renQ = InterpretScaleScheme(production, matrixElement, event_scales->renomalizationScheme, MomStore);
  double facQ = InterpretScaleScheme(production, matrixElement, event_scales->factorizationScheme, MomStore);
  SetAlphaS(renQ, facQ, event_scales->ren_scale_factor, event_scales->fac_scale_factor, 1, 5, "cteq6_l");
  double alphasVal, alphasmzVal;
  GetAlphaS(&alphasVal, &alphasmzVal);
  RcdME->setRenormalizationScale(renQ);
  RcdME->setFactorizationScale(facQ);
  RcdME->setAlphaS(alphasVal);
  RcdME->setAlphaSatMZ(alphasmzVal);
  RcdME->setHiggsMassWidth(masses_mcfm_.hmass, masses_mcfm_.hwidth, 0);
  RcdME->setHiggsMassWidth(spinzerohiggs_anomcoupl_.h2mass, spinzerohiggs_anomcoupl_.h2width, 1);
  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::BBHiggsMatEl: Set AlphaS:\n"
      << "\tBefore set, alphas scale: " << defaultRenScale << ", PDF scale: " << defaultFacScale << '\n'
      << "\trenQ: " << renQ << " ( x " << event_scales->ren_scale_factor << "), facQ: " << facQ << " ( x " << event_scales->fac_scale_factor << ")\n"
      << "\tAfter set, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }

  __modttbhiggs_MOD_evalxsec_pp_bbbh(p4, &botProcess, MatElsq);
  if (isUnknown[0] && isUnknown[1]){
    for (unsigned int ix=0; ix<4; ix++) swap(p4[3][ix], p4[4][ix]);
    __modttbhiggs_MOD_evalxsec_pp_bbbh(p4, &botProcess, MatElsq_tmp);
    for (int ix=0; ix<11; ix++){ for (int iy=0; iy<11; iy++) MatElsq[iy][ix] = (MatElsq[iy][ix]+MatElsq_tmp[iy][ix])/2.; }
  }

  int GeVexponent_MEsq = 4-(1+nRequested_AssociatedJets)*2;
  double constant = pow(GeV, -GeVexponent_MEsq);
  for (int ii=0; ii<nmsq; ii++){ for (int jj=0; jj<nmsq; jj++) MatElsq[jj][ii] *= constant; }
  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::BBHiggsMatEl: MEsq[ip][jp] = " << endl;
    for (int iquark=-5; iquark<=5; iquark++){
      for (int jquark=-5; jquark<=5; jquark++) MELAout << MatElsq[jquark+5][iquark+5] << '\t';
      MELAout << endl;
    }
  }
  sum_msqjk = SumMEPDF(MomStore[0], MomStore[1], MatElsq, RcdME, EBEAM, verbosity);

  if (verbosity>=TVar::DEBUG){
    MELAout
      << "TUtil::BBHiggsMatEl: Reset AlphaS:\n"
      << "\tBefore reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << endl;
  }
  SetAlphaS(defaultRenScale, defaultFacScale, 1., 1., defaultNloop, defaultNflav, defaultPdflabel);
  if (verbosity>=TVar::DEBUG){
    GetAlphaS(&alphasVal, &alphasmzVal);
    MELAout
      << "TUtil::BBHiggsMatEl: Reset AlphaS result:\n"
      << "\tAfter reset, alphas scale: " << scale_.scale << ", PDF scale: " << facscale_.facscale << ", alphas(Qren): " << alphasVal << ", alphas(MZ): " << alphasmzVal << endl;
  }
  return sum_msqjk;
}

// Wipe MEs
int TUtil::WipeMEArray(const TVar::Process& process, const TVar::Production& production, const int id[mxpart], double msq[nmsq][nmsq], const TVar::VerbosityLevel& verbosity){
  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::WipeMEArray: Initial MEsq:" << endl;
    for (int iquark=-5; iquark<=5; iquark++){
      for (int jquark=-5; jquark<=5; jquark++) MELAout << msq[jquark+5][iquark+5] << '\t';
      MELAout << endl;
    }
  }
  int nInstances=0;
  if (
    (production == TVar::ZZINDEPENDENT || production == TVar::ZZQQB)
    ||
    ((production == TVar::ZZQQB_STU || production == TVar::ZZQQB_S || production == TVar::ZZQQB_TU) && process == TVar::bkgZZ)
    ){
    // Sum over valid MEs without PDF weights
    // By far the dominant contribution is uub initial state.
    for (int iquark=-5; iquark<=5; iquark++){
      for (int jquark=-5; jquark<=5; jquark++){
        if (
          (PDGHelpers::isAnUnknownJet(id[0]) || (PDGHelpers::isAGluon(id[0]) && iquark==0) || iquark==id[0])
          &&
          (PDGHelpers::isAnUnknownJet(id[1]) || (PDGHelpers::isAGluon(id[1]) && jquark==0) || jquark==id[1])
          ){
          if (msq[jquark+5][iquark+5]>0.) nInstances++;
        }
        else msq[jquark+5][iquark+5]=0; // Kill the ME instance if mothers are known and their ids do not match the PDF indices.
      }
    }
  }
  else if (production == TVar::ZZGG){
    if (
      (PDGHelpers::isAnUnknownJet(id[0]) || PDGHelpers::isAGluon(id[0]))
      &&
      (PDGHelpers::isAnUnknownJet(id[1]) || PDGHelpers::isAGluon(id[1]))
      ){
      if (msq[5][5]>0.) nInstances=1;
    }
    else msq[5][5]=0; // Kill the ME instance if mothers are known and their ids do not match to gluons.
  }
  else if (
    production==TVar::Had_WH || production==TVar::Had_ZH
    || production==TVar::Had_WH_S || production==TVar::Had_ZH_S
    || production==TVar::Had_WH_TU || production==TVar::Had_ZH_TU
    || production==TVar::Lep_WH || production==TVar::Lep_ZH
    || production==TVar::Lep_WH_S || production==TVar::Lep_ZH_S
    || production==TVar::Lep_WH_TU || production==TVar::Lep_ZH_TU
    || production==TVar::JJVBF || production==TVar::JJEW || production==TVar::JJEWQCD
    || production==TVar::JJVBF_S || production==TVar::JJEW_S || production==TVar::JJEWQCD_S
    || production==TVar::JJVBF_TU || production==TVar::JJEW_TU || production==TVar::JJEWQCD_TU
    || production==TVar::JJQCD || production==TVar::JJQCD_S || production==TVar::JJQCD_TU
    ){
    // Z+2 jets
    if (process == TVar::bkgZJets){
      for (int iquark=-5; iquark<=5; iquark++){
        for (int jquark=-5; jquark<=5; jquark++){
          if (
            (PDGHelpers::isAnUnknownJet(id[0]) || (PDGHelpers::isAGluon(id[0]) && iquark==0) || iquark==id[0])
            &&
            (PDGHelpers::isAnUnknownJet(id[1]) || (PDGHelpers::isAGluon(id[1]) && jquark==0) || jquark==id[1])
            ){
            if (msq[jquark+5][iquark+5]>0.) nInstances++;
          }
          else msq[jquark+5][iquark+5]=0; // Kill the ME instance if mothers are known and their ids do not match the PDF indices.
        }
      }
    }
    // VBF or QCD MCFM SBI, S or B
    else{
      for (int iquark=-5; iquark<=5; iquark++){
        for (int jquark=-5; jquark<=5; jquark++){
          if (
            (
            PDGHelpers::isAnUnknownJet(id[0]) || (PDGHelpers::isAGluon(id[0]) && iquark==0) || iquark==id[0]
            )
            &&
            (
            PDGHelpers::isAnUnknownJet(id[1]) || (PDGHelpers::isAGluon(id[1]) && jquark==0) || jquark==id[1]
            )
            ){
            // Kill all non-VH initial parton states
            if (
              (
              (
              production==TVar::Had_WH || production==TVar::Lep_WH
              || production==TVar::Had_WH_S || production==TVar::Lep_WH_S
              || production==TVar::Had_WH_TU || production==TVar::Lep_WH_TU
              )
              &&
              (
              (TMath::Sign(1, iquark)==TMath::Sign(1, jquark))
              ||
              (PDGHelpers::isUpTypeQuark(iquark) && PDGHelpers::isUpTypeQuark(jquark))
              ||
              (PDGHelpers::isDownTypeQuark(iquark) && PDGHelpers::isDownTypeQuark(jquark))
              )
              )
              ||
              (
              (
              production==TVar::Had_ZH || production==TVar::Lep_ZH
              || production==TVar::Had_ZH_S || production==TVar::Lep_ZH_S
              || production==TVar::Had_ZH_TU || production==TVar::Lep_ZH_TU
              )
              &&
              (iquark!=-jquark)
              )
              ) msq[jquark+5][iquark+5]=0;

            // Check against a hash
            int order[2]={ -1, -1 };
            TMCFMUtils::AssociatedParticleOrdering_QQVVQQAny(iquark, jquark, id[6], id[7], order);
            if (order[0]<0 || order[1]<0) msq[jquark+5][iquark+5]=0;

            if (msq[jquark+5][iquark+5]>0.) nInstances++;
          }
          else msq[jquark+5][iquark+5]=0; // Kill the ME instance if mothers are known and their ids do not match the PDF indices.
        }
      }
    }
  }
  return nInstances;
}
// CheckPartonMomFraction computes xx[0:1] based on p0, p1
bool TUtil::CheckPartonMomFraction(const TLorentzVector& p0, const TLorentzVector& p1, double xx[2], const double& EBEAM, const TVar::VerbosityLevel& verbosity){
  //Make sure parton Level Energy fraction is [0,1]
  //phase space function already makes sure the parton energy fraction between [min,1]
  xx[0]=p0.P()/EBEAM;
  xx[1]=p1.P()/EBEAM;
  if (
    xx[0]>1. || xx[0]<=xmin_.xmin
    ||
    xx[1]>1. || xx[1]<=xmin_.xmin
    ||
    EBEAM<=0.
    ){
    if (verbosity>=TVar::ERROR){
      if (xx[0]>1. || xx[1]>1.) MELAerr << "TUtil::CheckPartonMomFraction: At least one of the parton momentum fractions is greater than 1." << endl;
      else if (xx[0]<=xmin_.xmin || xx[1]<=xmin_.xmin) MELAerr << "TUtil::CheckPartonMomFraction: At least one of the parton momentum fractions is less than or equal to " << xmin_.xmin << "." << endl;
      else MELAerr << "TUtil::CheckPartonMomFraction: EBEAM=" << EBEAM << "<=0." << endl;
    }
    return false;
  }
  else{
    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::CheckPartonMomFraction: xx[0]: " << xx[0] << ", xx[1] = " << xx[1] << ", xmin = " << xmin_.xmin << endl;
    return true;
  }
}
// ComputePDF does the PDF computation
void TUtil::ComputePDF(const TLorentzVector& p0, const TLorentzVector& p1, double fx1[nmsq], double fx2[nmsq], const double& EBEAM, const TVar::VerbosityLevel& verbosity){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin TUtil::ComputePDF"<< endl;
  double xx[2]={ 0 };
  if (CheckPartonMomFraction(p0, p1, xx, EBEAM, verbosity)){
    ///// USE JHUGEN SUBROUTINE (Accomodates LHAPDF) /////
    double fx1x2_jhu[2][13]={ { 0 } };
    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::ComputePDF: Calling setpdfs with xx[0]: " << xx[0] << ", xx[1] = " << xx[1] << endl;
    __modkinematics_MOD_setpdfs(&(xx[0]), &(xx[1]), fx1x2_jhu);
    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::ComputePDF: called"<< endl;
    for (int ip=-6; ip<=6; ip++){
      // JHUGen assignment is in JHU id coventions (up <-> down, g==g)
      int fac=0;
      if (ip!=0 && (abs(ip)%2==0)) fac=-1;
      else if (ip!=0) fac=+1;
      if (ip<0) fac=-fac;
      int jp=ip+fac;
      fx1[jp+5]=fx1x2_jhu[0][ip+6];
      fx2[jp+5]=fx1x2_jhu[1][ip+6];
    }
    /*
    ///// USE MCFM SUBROUTINE fdist_linux /////
    //Calculate Pdf
    //Parton Density Function is always evalualted at pT=0 frame
    //Always pass address through fortran function
    fdist_(&density_.ih1, &xx[0], &facscale_.facscale, fx1);
    fdist_(&density_.ih2, &xx[1], &facscale_.facscale, fx2);
    */
  }
  if (verbosity>=TVar::DEBUG){
    MELAout << "End TUtil::ComputePDF:"<< endl;
    for (int ip=-nf; ip<=nf; ip++) MELAout << "(fx1, fx2)[" << ip << "] = (" << fx1[ip+5] << " , " << fx2[ip+5] << ")" << endl;
  }
}
// SumMEPDF sums over all production parton flavors according to PDF and calls ComputePDF
double TUtil::SumMEPDF(const TLorentzVector& p0, const TLorentzVector& p1, double msq[nmsq][nmsq], MelaIO* RcdME, const double& EBEAM, const TVar::VerbosityLevel& verbosity){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin TUtil::SumMEPDF"<< endl;
  double fx1[nmsq]={ 0 };
  double fx2[nmsq]={ 0 };
  //double wgt_msq[nmsq][nmsq]={ { 0 } };

  ComputePDF(p0, p1, fx1, fx2, EBEAM, verbosity);
  if (verbosity>=TVar::DEBUG) MELAout << "TUtil::SumMEPDF: Setting RcdME"<< endl;
  RcdME->setPartonWeights(fx1, fx2);
  RcdME->setMEArray(msq, true);
  RcdME->computeWeightedMEArray();
  //RcdME->getWeightedMEArray(wgt_msq);
  if (verbosity>=TVar::DEBUG) MELAout << "End TUtil::SumMEPDF"<< endl;
  return RcdME->getSumME();
}


// Propagator reweighting
double TUtil::ResonancePropagator(double sqrts, TVar::ResonancePropagatorScheme scheme){
  __modjhugenmela_MOD_resetmubarhgabarh();

  const double GeV=1./100.; // JHUGen mom. scale factor
  int isch=(int)scheme;
  double shat_jhu = pow(sqrts*GeV, 2);
  double prop = __modkinematics_MOD_getbwpropagator(&shat_jhu, &isch);
  if (scheme!=TVar::NoPropagator) prop *= pow(GeV, 4);
  return prop;
}


// GetBoostedParticleVectors decomposes the MELACandidate object melaCand into mothers, daughters, associateds and tops
// and boosts them to the pT=0 frame for the particular mela_event.AssociationCode, which is a product of TVar::kUseAssociated* (prime numbers).
// If no mothers are present, it assigns mothers as Z>0, Z<0. If they are present, it orders them as "incoming" q-qbar / g-qbar / q-g / g-g
// Associated particles passed are different based on the code. If code==1, no associated particles are passed.
// mela_event.intermediateVids are needed to keep track of the decay mode. TVar::Process or TVar::Production do not keep track of the V/f decay modes.
// This is the major replacement functionality of the TVar lepton flavors.
void TUtil::GetBoostedParticleVectors(
  MELACandidate* melaCand,
  simple_event_record& mela_event,
  TVar::VerbosityLevel verbosity
  ){
  if (verbosity>=TVar::DEBUG) MELAout << "Begin GetBoostedParticleVectors" << endl;
  // This is the beginning of one long function.

  int code = mela_event.AssociationCode;
  if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: code=" << code << endl;
  int aVhypo = mela_event.AssociationVCompatibility;
  TLorentzVector nullFourVector(0, 0, 0, 0);

  SimpleParticleCollection_t daughters;
  vector<int> idVstar;
  if (melaCand->getNDaughters()==0){
    // Undecayed Higgs has V1=H, V2=empty, no sortedDaughters!
    daughters.push_back(SimpleParticle_t(melaCand->id, melaCand->p4));
    idVstar.push_back(melaCand->id);
    idVstar.push_back(-9000);
  }
  else{
    // H->ffb has V1=f->f, V2=fb->fb
    // H->GG has V1=G->G, V2=G->G
    // H->ZG has V1=Z->ffb, V2=G->G
    // Everything else is as expected.
    for (int iv=0; iv<2; iv++){ // 2 Vs are guaranteed in MELACandidate.
      MELAParticle* Vdau = melaCand->getSortedV(iv);
      if (Vdau!=0){
        int idtmp = Vdau->id;
        for (MELAParticle* Vdau_i:Vdau->getDaughters()){
          if (Vdau_i && Vdau_i->passSelection) daughters.push_back(SimpleParticle_t(Vdau_i->id, Vdau_i->p4));
        }
        if (idtmp!=0 || Vdau->getNDaughters()>0){ // Avoid "empty" intermediate Vs of the MELACandidate object
          if (Vdau->getNDaughters()>=2 && PDGHelpers::isAPhoton(idtmp)) idtmp=23; // Special case to avoid V->2f with massless decay mode (could happen by mistake)
          idVstar.push_back(idtmp);
        }
      }
      else idVstar.push_back(-9000);
    }
  }
  if (daughters.size()>=2){
    unsigned int nffs = daughters.size()/2;
    for (unsigned int iv=0; iv<nffs; iv++){
      pair<TLorentzVector, TLorentzVector> corrPair = TUtil::removeMassFromPair(
        daughters.at(2*iv+0).second, daughters.at(2*iv+0).first,
        daughters.at(2*iv+1).second, daughters.at(2*iv+1).first
        );
      daughters.at(2*iv+0).second = corrPair.first;
      daughters.at(2*iv+1).second = corrPair.second;
    }
    if (2*nffs<daughters.size()){
      TLorentzVector tmp = nullFourVector;
      pair<TLorentzVector, TLorentzVector> corrPair = TUtil::removeMassFromPair(
        daughters.at(daughters.size()-1).second, daughters.at(daughters.size()-1).first,
        tmp, -9000
        );
      daughters.at(daughters.size()-1).second = corrPair.first;
    }
  }

  /***** ASSOCIATED PARTICLES *****/
  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::GetBoostedParticleVectors: nRequestedLeps=" << mela_event.nRequested_AssociatedLeptons << endl;
    MELAout << "TUtil::GetBoostedParticleVectors: nRequestedJets=" << mela_event.nRequested_AssociatedJets << endl;
    MELAout << "TUtil::GetBoostedParticleVectors: nRequestedPhotons=" << mela_event.nRequested_AssociatedPhotons << endl;
  }
  int nsatisfied_jets=0;
  int nsatisfied_lnus=0;
  int nsatisfied_gammas=0;
  vector<MELAParticle*> candidateVs; // Used if aVhypo!=0
  SimpleParticleCollection_t associated;
  if (aVhypo!=0){
    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: aVhypo!=0 case start" << endl;

    vector<MELAParticle*> asortedVs = melaCand->getAssociatedSortedVs();
    for (MELAParticle* Vdau:asortedVs){ // Loop over associated Vs
      if (Vdau!=0){
        bool doAdd=false;
        int idV = Vdau->id;
        if ((abs(idV)==aVhypo || idV==0) && Vdau->getNDaughters()>0 && Vdau->passSelection){ // If the V is unknown or compatible with the requested hypothesis
          doAdd=true;
          for (MELAParticle* Vdau_i:Vdau->getDaughters()){ // Loop over the daughters of V
            if (Vdau_i==0){ doAdd=false; break; }
            else if (
              (mela_event.nRequested_AssociatedLeptons==0 && (PDGHelpers::isALepton(Vdau_i->id) || PDGHelpers::isANeutrino(Vdau_i->id)))
              ||
              (mela_event.nRequested_AssociatedJets==0 && PDGHelpers::isAJet(Vdau_i->id))
              ||
              !Vdau_i->passSelection // Protection against incorrect Vdau passSelection flag
              ){
              doAdd=false; break;
            }
          }
        }
        if (doAdd) candidateVs.push_back(Vdau);
      }
    } // End loop over associated Vs

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: candidateVs size = " << candidateVs.size() << endl;

    // Pick however many candidates necessary to fill up the requested number of jets or lepton(+)neutrinos
    for (MELAParticle* Vdau:candidateVs){
      SimpleParticleCollection_t associated_tmp;
      for (MELAParticle* part:Vdau->getDaughters()){ // Loop over the daughters of V
        if (
          part->passSelection
          &&
          (PDGHelpers::isALepton(part->id) || PDGHelpers::isANeutrino(part->id))
          &&
          nsatisfied_lnus<mela_event.nRequested_AssociatedLeptons
          ){
          nsatisfied_lnus++;
          associated_tmp.push_back(SimpleParticle_t(part->id, part->p4));
        }
        else if (
          part->passSelection
          &&
          PDGHelpers::isAJet(part->id)
          &&
          nsatisfied_jets<mela_event.nRequested_AssociatedJets
          ){
          nsatisfied_jets++;
          associated_tmp.push_back(SimpleParticle_t(part->id, part->p4));
        }
      }
      // Add the id of the intermediate V into the idVstar array
      if (associated_tmp.size()>=2 || (associated_tmp.size()==1 && PDGHelpers::isAPhoton(associated_tmp.at(0).first))) idVstar.push_back(Vdau->id); // Only add the V-id of pairs that passed
      // Adjust the kinematics of associated V-originating particles
      if (associated_tmp.size()>=2){ // ==1 means a photon, so omit it here.
        unsigned int nffs = associated_tmp.size()/2;
        for (unsigned int iv=0; iv<nffs; iv++){
          pair<TLorentzVector, TLorentzVector> corrPair = TUtil::removeMassFromPair(
            associated_tmp.at(2*iv+0).second, associated_tmp.at(2*iv+0).first,
            associated_tmp.at(2*iv+1).second, associated_tmp.at(2*iv+1).first
            );
          associated_tmp.at(2*iv+0).second = corrPair.first;
          associated_tmp.at(2*iv+1).second = corrPair.second;
        }
        if (2*nffs<associated_tmp.size()){
          TLorentzVector tmp = nullFourVector;
          pair<TLorentzVector, TLorentzVector> corrPair = TUtil::removeMassFromPair(
            associated_tmp.at(associated_tmp.size()-1).second, associated_tmp.at(associated_tmp.size()-1).first,
            tmp, -9000
            );
          associated_tmp.at(associated_tmp.size()-1).second = corrPair.first;
        }
      }
      for (SimpleParticle_t& sp:associated_tmp) associated.push_back(sp); // Fill associated at the last step
    }

    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: aVhypo!=0 case associated.size=" << associated.size() << endl;
  }
  else{ // Could be split to aVhypo==0 and aVhypo<0 if associated V+jets is needed
    // Associated leptons
    if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: aVhypo==0 case begin" << endl;

    if (code%TVar::kUseAssociated_Leptons==0){
      SimpleParticleCollection_t associated_tmp;
      for (MELAParticle* part:melaCand->getAssociatedLeptons()){
        if (part!=0 && part->passSelection && nsatisfied_lnus<mela_event.nRequested_AssociatedLeptons){
          nsatisfied_lnus++;
          associated_tmp.push_back(SimpleParticle_t(part->id, part->p4));
        }
      }

      if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: aVhypo==0 lep case associated_tmp.size=" << associated_tmp.size() << endl;

      // Adjust the kinematics of associated non-V-originating particles
      if (associated_tmp.size()>=1){
        unsigned int nffs = associated_tmp.size()/2;
        for (unsigned int iv=0; iv<nffs; iv++){
          if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: Removing mass from lepton pair " << 2*iv+0 << '\t' << 2*iv+1 << endl;
          pair<TLorentzVector, TLorentzVector> corrPair = TUtil::removeMassFromPair(
            associated_tmp.at(2*iv+0).second, associated_tmp.at(2*iv+0).first,
            associated_tmp.at(2*iv+1).second, associated_tmp.at(2*iv+1).first
            );
          associated_tmp.at(2*iv+0).second = corrPair.first;
          associated_tmp.at(2*iv+1).second = corrPair.second;
        }
        if (2*nffs<associated_tmp.size()){
          if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: Removing mass from last lepton  " << associated_tmp.size()-1 << endl;
          TLorentzVector tmp = nullFourVector;
          pair<TLorentzVector, TLorentzVector> corrPair = TUtil::removeMassFromPair(
            associated_tmp.at(associated_tmp.size()-1).second, associated_tmp.at(associated_tmp.size()-1).first,
            tmp, -9000
            );
          associated_tmp.at(associated_tmp.size()-1).second = corrPair.first;
        }
      }
      for (SimpleParticle_t& sp:associated_tmp) associated.push_back(sp); // Fill associated at the last step
    }
    // Associated jets
    if (code%TVar::kUseAssociated_Jets==0){
      SimpleParticleCollection_t associated_tmp;
      for (MELAParticle* part:melaCand->getAssociatedJets()){
        if (part!=0 && part->passSelection && nsatisfied_jets<mela_event.nRequested_AssociatedJets){
          nsatisfied_jets++;
          associated_tmp.push_back(SimpleParticle_t(part->id, part->p4));
        }
      }

      if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: aVhypo==0 jet case associated_tmp.size=" << associated_tmp.size() << endl;

      // Adjust the kinematics of associated non-V-originating particles
      if (associated_tmp.size()>=1){
        unsigned int nffs = associated_tmp.size()/2;
        for (unsigned int iv=0; iv<nffs; iv++){
          if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: Removing mass from jet pair " << 2*iv+0 << '\t' << 2*iv+1 << endl;
          pair<TLorentzVector, TLorentzVector> corrPair = TUtil::removeMassFromPair(
            associated_tmp.at(2*iv+0).second, associated_tmp.at(2*iv+0).first,
            associated_tmp.at(2*iv+1).second, associated_tmp.at(2*iv+1).first
            );
          associated_tmp.at(2*iv+0).second = corrPair.first;
          associated_tmp.at(2*iv+1).second = corrPair.second;
        }
        if (2*nffs<associated_tmp.size()){
          if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: Removing mass from last jet  " << associated_tmp.size()-1 << endl;
          TLorentzVector tmp = nullFourVector;
          pair<TLorentzVector, TLorentzVector> corrPair = TUtil::removeMassFromPair(
            associated_tmp.at(associated_tmp.size()-1).second, associated_tmp.at(associated_tmp.size()-1).first,
            tmp, -9000
            );
          associated_tmp.at(associated_tmp.size()-1).second = corrPair.first;
        }
      }
      for (SimpleParticle_t& sp:associated_tmp) associated.push_back(sp); // Fill associated at the last step
    }

  } // End if(aVhypo!=0)-else statement

  if (code%TVar::kUseAssociated_Photons==0){
    for (MELAParticle* part:melaCand->getAssociatedPhotons()){
      if (part!=0 && part->passSelection && nsatisfied_gammas<mela_event.nRequested_AssociatedPhotons){
        nsatisfied_gammas++;
        associated.push_back(SimpleParticle_t(part->id, part->p4));
      }
    }
  }
  /***** END ASSOCIATED PARTICLES *****/

  if (verbosity>=TVar::DEBUG) MELAout << "TUtil::GetBoostedParticleVectors: associated.size=" << associated.size() << endl;

  /***** ASSOCIATED TOP OBJECTS *****/
  int nsatisfied_tops=0;
  int nsatisfied_antitops=0;
  vector<SimpleParticleCollection_t> topDaughters;
  vector<SimpleParticleCollection_t> antitopDaughters;
  SimpleParticleCollection_t stableTops;
  SimpleParticleCollection_t stableAntitops;

  vector<MELATopCandidate*> tops;
  vector<MELATopCandidate*> topbars;
  vector<MELATopCandidate*> unknowntops;
  if (code%TVar::kUseAssociated_StableTops==0 && code%TVar::kUseAssociated_UnstableTops==0 && verbosity>=TVar::INFO) MELAerr << "TUtil::GetBoostedParticleVectors: Stable and unstable tops are not supported at the same time!"  << endl;
  else if (code%TVar::kUseAssociated_StableTops==0 || code%TVar::kUseAssociated_UnstableTops==0){

    for (MELATopCandidate* theTop:melaCand->getAssociatedTops()){
      if (theTop!=0 && theTop->passSelection){
        vector<MELATopCandidate*>* particleArray;
        if (theTop->id==6) particleArray = &tops;
        else if (theTop->id==-6) particleArray = &topbars;
        else particleArray = &unknowntops;
        if (
          (code%TVar::kUseAssociated_StableTops==0)
          ||
          (theTop->getNDaughters()==3 && code%TVar::kUseAssociated_UnstableTops==0)
          ) particleArray->push_back((MELATopCandidate*)theTop);
      }
    }
    if (verbosity>=TVar::DEBUG){ MELAout << "TUtil::GetBoostedParticleVectors: tops.size=" << tops.size() << ", topbars.size=" << topbars.size() << ", unknowntops.size=" << unknowntops.size() << endl; }

    // Fill the stable/unstable top arrays
    for (MELATopCandidate* theTop:tops){
      if (code%TVar::kUseAssociated_StableTops==0 && nsatisfied_tops<mela_event.nRequested_Tops){ // Case with no daughters needed
        nsatisfied_tops++;
        stableTops.push_back(SimpleParticle_t(theTop->id, theTop->p4));
      }
      else if (code%TVar::kUseAssociated_UnstableTops==0 && theTop->getNDaughters()==3 && nsatisfied_tops<mela_event.nRequested_Tops){ // Case with daughters needed
        nsatisfied_tops++;
        SimpleParticleCollection_t vdaughters;

        MELAParticle* bottom = theTop->getLightQuark();
        MELAParticle* Wf = theTop->getWFermion();
        MELAParticle* Wfb = theTop->getWAntifermion();
        if (bottom!=0) vdaughters.push_back(SimpleParticle_t(bottom->id, bottom->p4));
        if (Wf!=0) vdaughters.push_back(SimpleParticle_t(Wf->id, Wf->p4));
        if (Wfb!=0) vdaughters.push_back(SimpleParticle_t(Wfb->id, Wfb->p4));

        TUtil::adjustTopDaughters(vdaughters); // Adjust top daughter kinematics
        if (vdaughters.size()==3) topDaughters.push_back(vdaughters);
      }
    }
    for (MELATopCandidate* theTop:topbars){
      if (code%TVar::kUseAssociated_StableTops==0 && nsatisfied_antitops<mela_event.nRequested_Antitops){ // Case with no daughters needed
        nsatisfied_antitops++;
        stableAntitops.push_back(SimpleParticle_t(theTop->id, theTop->p4));
      }
      else if (code%TVar::kUseAssociated_UnstableTops==0 && nsatisfied_antitops<mela_event.nRequested_Antitops){ // Case with daughters needed
        nsatisfied_antitops++;
        SimpleParticleCollection_t vdaughters;

        MELAParticle* bottom = theTop->getLightQuark();
        MELAParticle* Wf = theTop->getWFermion();
        MELAParticle* Wfb = theTop->getWAntifermion();
        if (bottom!=0) vdaughters.push_back(SimpleParticle_t(bottom->id, bottom->p4));
        if (Wf!=0) vdaughters.push_back(SimpleParticle_t(Wf->id, Wf->p4));
        if (Wfb!=0) vdaughters.push_back(SimpleParticle_t(Wfb->id, Wfb->p4));

        TUtil::adjustTopDaughters(vdaughters); // Adjust top daughter kinematics
        if (vdaughters.size()==3) antitopDaughters.push_back(vdaughters);
      }
      else break;
    }
    // Loop over the unknown-id tops
    // Fill tops, then antitops from the unknown tops
    for (MELATopCandidate* theTop:unknowntops){
      // t, then tb cases with no daughters needed
      if (code%TVar::kUseAssociated_StableTops==0 && nsatisfied_tops<mela_event.nRequested_Tops){
        nsatisfied_tops++;
        stableTops.push_back(SimpleParticle_t(theTop->id, theTop->p4));
      }
      else if (code%TVar::kUseAssociated_StableTops==0 && nsatisfied_antitops<mela_event.nRequested_Antitops){
        nsatisfied_antitops++;
        stableAntitops.push_back(SimpleParticle_t(theTop->id, theTop->p4));
      }
      // t, then tb cases with daughters needed
      else if (code%TVar::kUseAssociated_UnstableTops==0 && nsatisfied_tops<mela_event.nRequested_Tops){
        nsatisfied_tops++;
        SimpleParticleCollection_t vdaughters;

        MELAParticle* bottom = theTop->getLightQuark();
        MELAParticle* Wf = theTop->getWFermion();
        MELAParticle* Wfb = theTop->getWAntifermion();
        if (bottom!=0) vdaughters.push_back(SimpleParticle_t(bottom->id, bottom->p4));
        if (Wf!=0) vdaughters.push_back(SimpleParticle_t(Wf->id, Wf->p4));
        if (Wfb!=0) vdaughters.push_back(SimpleParticle_t(Wfb->id, Wfb->p4));

        TUtil::adjustTopDaughters(vdaughters); // Adjust top daughter kinematics
        if (vdaughters.size()==3) topDaughters.push_back(vdaughters);
      }
      else if (code%TVar::kUseAssociated_UnstableTops==0 && nsatisfied_antitops<mela_event.nRequested_Antitops){
        nsatisfied_antitops++;
        SimpleParticleCollection_t vdaughters;

        MELAParticle* bottom = theTop->getLightQuark();
        MELAParticle* Wf = theTop->getWFermion();
        MELAParticle* Wfb = theTop->getWAntifermion();
        if (bottom!=0) vdaughters.push_back(SimpleParticle_t(bottom->id, bottom->p4));
        if (Wf!=0) vdaughters.push_back(SimpleParticle_t(Wf->id, Wf->p4));
        if (Wfb!=0) vdaughters.push_back(SimpleParticle_t(Wfb->id, Wfb->p4));

        TUtil::adjustTopDaughters(vdaughters); // Adjust top daughter kinematics
        if (vdaughters.size()==3) antitopDaughters.push_back(vdaughters);
      }
    }

  }
  /***** END ASSOCIATED TOP OBJECTS *****/
  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::GetBoostedParticleVectors: stableTops.size=" << stableTops.size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors: stableAntitops.size=" << stableAntitops.size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors: topDaughters.size=" << topDaughters.size() << endl;
    for (unsigned int itop=0; itop<topDaughters.size(); itop++) MELAout << "TUtil::GetBoostedParticleVectors: topDaughters.at(" << itop << ").size=" << topDaughters.at(itop).size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors: antitopDaughters.size=" << antitopDaughters.size() << endl;
    for (unsigned int itop=0; itop<antitopDaughters.size(); itop++) MELAout << "TUtil::GetBoostedParticleVectors: antitopDaughters.at(" << itop << ").size=" << antitopDaughters.at(itop).size() << endl;
  }

  /***** BOOSTS TO THE CORRECT PT=0 FRAME *****/
  // Gather all final state particles collected for this frame
  TLorentzVector pTotal(0, 0, 0, 0);
  for (SimpleParticle_t const& sp:daughters) pTotal = pTotal + sp.second;
  for (SimpleParticle_t const& sp:associated) pTotal = pTotal + sp.second;
  for (SimpleParticle_t const& sp:stableTops) pTotal = pTotal + sp.second;
  for (SimpleParticle_t const& sp:stableAntitops) pTotal = pTotal + sp.second;
  for (SimpleParticleCollection_t const& spc:topDaughters){ for (SimpleParticle_t const& sp:spc) pTotal = pTotal + sp.second; }
  for (SimpleParticleCollection_t const& spc:antitopDaughters){ for (SimpleParticle_t const& sp:spc) pTotal = pTotal + sp.second; }

  // Get the boost vector and boost all final state particles
  double qX = pTotal.X();
  double qY = pTotal.Y();
  double qE = pTotal.T();
  if ((qX*qX+qY*qY)>0.){
    TVector3 boostV(-qX/qE, -qY/qE, 0.);
    for (SimpleParticle_t& sp:daughters) sp.second.Boost(boostV);
    for (SimpleParticle_t& sp:associated) sp.second.Boost(boostV);
    for (SimpleParticle_t& sp:stableTops) sp.second.Boost(boostV);
    for (SimpleParticle_t& sp:stableAntitops) sp.second.Boost(boostV);
    for (SimpleParticleCollection_t& spc:topDaughters){ for (SimpleParticle_t& sp:spc) sp.second.Boost(boostV); }
    for (SimpleParticleCollection_t& spc:antitopDaughters){ for (SimpleParticle_t& sp:spc) sp.second.Boost(boostV); }
    pTotal.Boost(boostV);
  }

  // Mothers need special treatment:
  // In case they are undefined, mother id is unknown, and mothers are ordered as M1==(pz>0) and M2==(pz<0).
  // In case they are defined, mothers are first matched to the assumed momenta through their pz's.
  // They are ordered by M1==f, M2==fb afterward if this is possible.
  // Notice that the momenta of the mother objects are not actually used. If the event is truly a gen. particle, everything would be in the pT=0 frame to begin with, and everybody is happy in reweighting.
  double sysPz= pTotal.Z();
  double sysE = pTotal.T();
  double pz0 = (sysE+sysPz)/2.;
  double pz1 = -(sysE-sysPz)/2.;
  double E0 = pz0;
  double E1 = -pz1;
  int motherId[2]={ 0, 0 };
  if (melaCand->getNMothers()==2){
    for (int ip=0; ip<2; ip++) motherId[ip]=melaCand->getMother(ip)->id;
    // Match the "assumed" M'1==(larger signed pz) and M'2==(smaller signed pz) to the pz of the actual mothers:
    // Swap pZs to get the correct momentum matching default M1, M2.
    if (TMath::Sign(1., melaCand->getMother(0)->z()-melaCand->getMother(1)->z())!=TMath::Sign(1., pz0-pz1)){ swap(pz0, pz1); swap(E0, E1); }
    // Swap the ids of mothers and their pT=0-assumed momenta to achieve ordering as "incoming" q-qbar
    if ((motherId[0]<0 && motherId[1]>=0) || (motherId[1]>0 && motherId[0]<=0)){
      swap(pz0, pz1);
      swap(E0, E1);
      swap(motherId[0], motherId[1]);
    }
  }
  TLorentzVector pM[2];
  pM[0].SetXYZT(0., 0., pz0, E0);
  pM[1].SetXYZT(0., 0., pz1, E1);

  // Fill the ids of the V intermediates to the candidate daughters
  mela_event.intermediateVid.clear();
  TUtilHelpers::copyVector(idVstar, mela_event.intermediateVid);
  // Fill the mothers
  mela_event.pMothers.clear();
  for (unsigned int ip=0; ip<2; ip++){ mela_event.pMothers.push_back(SimpleParticle_t(motherId[ip], pM[ip])); }
  // Fill the daughters
  mela_event.pDaughters.clear();
  TUtilHelpers::copyVector(daughters, mela_event.pDaughters);
  // Fill the associated particles
  mela_event.pAssociated.clear();
  TUtilHelpers::copyVector(associated, mela_event.pAssociated);
  // Fill the stable tops and antitops
  mela_event.pStableTops.clear();
  TUtilHelpers::copyVector(stableTops, mela_event.pStableTops);
  mela_event.pStableAntitops.clear();
  TUtilHelpers::copyVector(stableAntitops, mela_event.pStableAntitops);
  // Fill the daughters of unstable tops and antitops
  mela_event.pTopDaughters.clear();
  TUtilHelpers::copyVector(topDaughters, mela_event.pTopDaughters);
  mela_event.pAntitopDaughters.clear();
  TUtilHelpers::copyVector(antitopDaughters, mela_event.pAntitopDaughters);

  // This is the end of one long function.
  if (verbosity>=TVar::DEBUG){
    MELAout << "TUtil::GetBoostedParticleVectors mela_event.intermediateVid.size=" << mela_event.intermediateVid.size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors mela_event.pMothers.size=" << mela_event.pMothers.size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors mela_event.pDaughters.size=" << mela_event.pDaughters.size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors mela_event.pAssociated.size=" << mela_event.pAssociated.size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors mela_event.pStableTops.size=" << mela_event.pStableTops.size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors mela_event.pStableAntitops.size=" << mela_event.pStableAntitops.size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors mela_event.pTopDaughters.size=" << mela_event.pTopDaughters.size() << endl;
    MELAout << "TUtil::GetBoostedParticleVectors mela_event.pAntitopDaughters.size=" << mela_event.pAntitopDaughters.size() << endl;
    MELAout << "End GetBoostedParticleVectors" << endl;
  }
}

// Convert vectors of simple particles to MELAParticles and create a MELACandidate
// The output lists could be members of TEvtProb directly.
MELACandidate* TUtil::ConvertVectorFormat(
  // Inputs
  SimpleParticleCollection_t* pDaughters,
  SimpleParticleCollection_t* pAssociated,
  SimpleParticleCollection_t* pMothers,
  bool isGen,
  // Outputs
  std::vector<MELAParticle*>* particleList,
  std::vector<MELACandidate*>* candList
  ){
  MELACandidate* cand=0;

  if (pDaughters==0){ MELAerr << "TUtil::ConvertVectorFormat: No daughters!" << endl; return cand; }
  else if (pDaughters->size()==0){ MELAerr << "TUtil::ConvertVectorFormat: Daughter size==0!" << endl; return cand; }
  else if (pDaughters->size()>4){ MELAerr << "TUtil::ConvertVectorFormat: Daughter size " << pDaughters->size() << ">4 is not supported!" << endl; return cand; }
  if (pMothers!=0 && pMothers->size()!=2){ MELAerr << "TUtil::ConvertVectorFormat: Mothers momentum size (" << pMothers->size() << ") has to have had been 2! Continuing by omitting mothers." << endl; /*return cand;*/ }

  // Create mother, daughter and associated particle MELAParticle objects
  std::vector<MELAParticle*> daughters;
  std::vector<MELAParticle*> aparticles;
  std::vector<MELAParticle*> mothers;
  for (auto& spart:(*pDaughters)){
    MELAParticle* onePart = new MELAParticle(spart.first, spart.second);
    onePart->setGenStatus(1); // Final state status
    if (particleList!=0) particleList->push_back(onePart);
    daughters.push_back(onePart);
  }
  if (pAssociated!=0){
    for (auto& spart:(*pAssociated)){
      MELAParticle* onePart = new MELAParticle(spart.first, spart.second);
      onePart->setGenStatus(1); // Final state status
      if (particleList!=0) particleList->push_back(onePart);
      aparticles.push_back(onePart);
    }
  }
  if (pMothers!=0 && pMothers->size()==2){
    for (auto& spart:(*pMothers)){
      MELAParticle* onePart = new MELAParticle(spart.first, spart.second);
      onePart->setGenStatus(-1); // Mother status
      if (particleList!=0) particleList->push_back(onePart);
      mothers.push_back(onePart);
    }
  }

  // Create the candidate
  /***** Adaptation of LHEAnalyzer::Event::constructVVCandidates *****/
  /*
  The assumption is that the daughters make sense for either ffb, gamgam, Zgam, ZZ or WW.
  No checking is done on whether particle-antiparticle pairing is correct when necessary.
  If not, you will get a seg. fault!
  */

  // Undecayed Higgs
  if (daughters.size()==1){
    cand = new MELACandidate(25, (daughters.at(0))->p4); // No sorting!
    cand->setDecayMode(TVar::CandidateDecay_Stable); // Explicitly set the decay mode in case the default of MELACandidate changes in the future
  }
  // GG / ff final states
  else if (daughters.size()==2){
    MELAParticle* F1 = daughters.at(0);
    MELAParticle* F2 = daughters.at(1);
    TLorentzVector pH = F1->p4+F2->p4;
    cand = new MELACandidate(25, pH);
    cand->addDaughter(F1);
    cand->addDaughter(F2);
    TVar::CandidateDecayMode defaultHDecayMode = PDGHelpers::HDecayMode;
    if (PDGHelpers::isAPhoton(F1->id) && PDGHelpers::isAPhoton(F2->id)) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_GG);
    else PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ff);
    cand->sortDaughters();
    PDGHelpers::setCandidateDecayMode(defaultHDecayMode);
  }
  // ZG / WG
  else if (daughters.size()==3){
    MELAParticle* F1 = daughters.at(0);
    MELAParticle* F2 = daughters.at(1);
    MELAParticle* gamma = daughters.at(2);
    if (PDGHelpers::isAPhoton(F1->id)){
      MELAParticle* tmp = F1;
      F1 = gamma;
      gamma = tmp;
    }
    else if (PDGHelpers::isAPhoton(F2->id)){
      MELAParticle* tmp = F2;
      F2 = gamma;
      gamma = tmp;
    }
    TLorentzVector pH = F1->p4+F2->p4+gamma->p4;
    double charge = F1->charge()+F2->charge()+gamma->charge();
    cand = new MELACandidate(25, pH);
    cand->addDaughter(F1);
    cand->addDaughter(F2);
    cand->addDaughter(gamma);
    TVar::CandidateDecayMode defaultHDecayMode = PDGHelpers::HDecayMode;
    if (fabs(charge)<0.01) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ZG); // ZG
    else PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_WG); // WG, un-tested
    cand->sortDaughters();
    PDGHelpers::setCandidateDecayMode(defaultHDecayMode);
  }
  // ZZ / WW / ZW (!)
  else/* if (daughters.size()==4)*/{
    TLorentzVector pH(0, 0, 0, 0);
    double charge = 0.;
    for (int ip=0; ip<4; ip++){ pH = pH + (daughters.at(ip))->p4; charge += (daughters.at(ip))->charge(); }
    cand = new MELACandidate(25, pH);
    for (int ip=0; ip<4; ip++) cand->addDaughter(daughters.at(ip));
    // FIXME/REIMPLEMENT: ZW trickier than I thought: Summing over charges over all 4f is not enough, affects SSSF pairing in ZZ
    //TVar::CandidateDecayMode defaultHDecayMode = PDGHelpers::HDecayMode;
    //if (fabs(charge)>0.01) PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_ZW);
    if (
      PDGHelpers::HDecayMode!=TVar::CandidateDecay_WW
      && PDGHelpers::HDecayMode!=TVar::CandidateDecay_ZZ
      && PDGHelpers::HDecayMode!=TVar::CandidateDecay_ZW
      ){
      MELAerr << "TUtil::ConvertVectorFormat: PDGHelpers::HDecayMode = " << PDGHelpers::HDecayMode << " is set incorrectly for Ndaughters=" << daughters.size() << ". There is no automatic mechnism for this scenario. ";
      MELAerr << "Suggestion: Call PDGHelpers::setCandidateDecayMode(TVar::CandidateDecay_XY) before this call." << endl;
      assert(0);
    }
    cand->sortDaughters();
    //PDGHelpers::setCandidateDecayMode(defaultHDecayMode);
  }

  /***** Adaptation of LHEAnalyzer::Event::addVVCandidateMother *****/
  if (mothers.size()>0){ // ==2
    for (auto& part:mothers) cand->addMother(part);
    if (isGen) cand->setGenStatus(-1); // Candidate is a gen. particle!
  }
  /***** Adaptation of LHEAnalyzer::Event::addVVCandidateAppendages *****/
  if (aparticles.size()>0){
    for (auto& part:aparticles){
      const int& partId = part->id;
      if (PDGHelpers::isALepton(partId)) cand->addAssociatedLeptons(part);
      else if (PDGHelpers::isANeutrino(partId)) cand->addAssociatedNeutrinos(part); // Be careful: Neutrinos are neutrinos, but also "leptons" in MELACandidate!
      else if (PDGHelpers::isAPhoton(partId)) cand->addAssociatedPhotons(part);
      else if (PDGHelpers::isAJet(partId)) cand->addAssociatedJets(part);
    }
    cand->addAssociatedVs(); // For the VH topology
  }

  if (candList!=0 && cand!=0) candList->push_back(cand);
  return cand;
}

// Convert the vector of top daughters (as simple particles) to MELAParticles and create a MELATopCandidate
// The output lists could be members of TEvtProb directly.
MELATopCandidate* TUtil::ConvertTopCandidate(
  // Input
  SimpleParticleCollection_t* TopDaughters,
  // Outputs
  std::vector<MELAParticle*>* particleList,
  std::vector<MELATopCandidate*>* topCandList
  ){
  MELATopCandidate* cand=0;

  if (TopDaughters==0){ MELAerr << "TUtil::ConvertTopCandidate: No daughters!" << endl; return cand; }
  else if (TopDaughters->size()==0){ MELAerr << "TUtil::ConvertTopCandidate: Daughter size==0!" << endl; return cand; }
  else if (!(TopDaughters->size()==1 || TopDaughters->size()==3)){ MELAerr << "TUtil::ConvertVectorFormat: Daughter size " << TopDaughters->size() << "!=1 or 3 is not supported!" << endl; return cand; }

  if (TopDaughters->size()==1){
    if (abs((TopDaughters->at(0)).first)==6 || (TopDaughters->at(0)).first==0){
      cand = new MELATopCandidate((TopDaughters->at(0)).first, (TopDaughters->at(0)).second);
      topCandList->push_back(cand);
    }
  }
  else if (TopDaughters->size()==3){
    MELAParticle* bottom = new MELAParticle((TopDaughters->at(0)).first, (TopDaughters->at(0)).second);
    MELAParticle* Wf = new MELAParticle((TopDaughters->at(1)).first, (TopDaughters->at(1)).second);
    MELAParticle* Wfb = new MELAParticle((TopDaughters->at(2)).first, (TopDaughters->at(2)).second);

    if (Wf->id<0 || Wfb->id>0){
      MELAParticle* parttmp = Wf;
      Wf=Wfb;
      Wfb=parttmp;
    }

    particleList->push_back(bottom);
    particleList->push_back(Wf);
    particleList->push_back(Wfb);

    cand = new MELATopCandidate(bottom, Wf, Wfb);
    topCandList->push_back(cand);
  }
  return cand;
}

void TUtil::PrintCandidateSummary(MELACandidate* cand){
  MELAout << "***** TUtil::PrintCandidateSummary *****" << endl;
  MELAout << "Candidate: " << cand << endl;
  if (cand) MELAout << *cand;
}

void TUtil::PrintCandidateSummary(simple_event_record* cand){
  MELAout << "***** TUtil::PrintCandidateSummary (Simple Event Record) *****" << endl;
  MELAout << "Candidate: " << cand << endl;
  if (cand!=0){
    MELAout << "\tAssociationCode: " << cand->AssociationCode << endl;
    MELAout << "\tAssociationVCompatibility: " << cand->AssociationVCompatibility << endl;
    MELAout << "\tnRequested_AssociatedJets: " << cand->nRequested_AssociatedJets << endl;
    MELAout << "\tnRequested_AssociatedLeptons: " << cand->nRequested_AssociatedLeptons << endl;
    MELAout << "\tnRequested_AssociatedPhotons: " << cand->nRequested_AssociatedPhotons << endl;
    MELAout << "\tnRequested_Tops: " << cand->nRequested_Tops << endl;
    MELAout << "\tnRequested_Antitops: " << cand->nRequested_Antitops << endl;
    MELAout << "\tHas " << cand->pMothers.size() << " mothers" << endl;
    for (unsigned int ip=0; ip<cand->pMothers.size(); ip++){
      SimpleParticle_t* part = &(cand->pMothers.at(ip));
      MELAout
        << "\t\tV" << ip << " (" << part->first << ") (X,Y,Z,T)=( "
        << part->second.X() << " , "
        << part->second.Y() << " , "
        << part->second.Z() << " , "
        << part->second.T() << " )" << endl;
    }
    MELAout << "\tHas " << cand->intermediateVid.size() << " sorted daughter Vs" << endl;
    for (unsigned int iv=0; iv<cand->intermediateVid.size(); iv++) MELAout << "\t\tV" << iv << " (" << cand->intermediateVid.at(iv) << ")" << endl;
    MELAout << "\tHas " << cand->pDaughters.size() << " daughters" << endl;
    for (unsigned int ip=0; ip<cand->pDaughters.size(); ip++){
      SimpleParticle_t* part = &(cand->pDaughters.at(ip));
      MELAout
        << "\t\tDau[" << ip << "] (" << part->first << ") (X,Y,Z,T)=( "
        << part->second.X() << " , "
        << part->second.Y() << " , "
        << part->second.Z() << " , "
        << part->second.T() << " )" << endl;
    }
    MELAout << "\tHas " << cand->pAssociated.size() << " associated particles" << endl;
    for (unsigned int ip=0; ip<cand->pAssociated.size(); ip++){
      SimpleParticle_t* part = &(cand->pAssociated.at(ip));
      MELAout
        << "\t\tAPart[" << ip << "] (" << part->first << ") (X,Y,Z,T)=( "
        << part->second.X() << " , "
        << part->second.Y() << " , "
        << part->second.Z() << " , "
        << part->second.T() << " )" << endl;
    }
    MELAout << "\tHas " << cand->pStableTops.size() << " stable tops" << endl;
    for (unsigned int ip=0; ip<cand->pStableTops.size(); ip++){
      SimpleParticle_t* part = &(cand->pStableTops.at(ip));
      MELAout
        << "\t\tAPart[" << ip << "] (" << part->first << ") (X,Y,Z,T)=( "
        << part->second.X() << " , "
        << part->second.Y() << " , "
        << part->second.Z() << " , "
        << part->second.T() << " )" << endl;
    }
    MELAout << "\tHas " << cand->pStableAntitops.size() << " stable antitops" << endl;
    for (unsigned int ip=0; ip<cand->pStableAntitops.size(); ip++){
      SimpleParticle_t* part = &(cand->pStableAntitops.at(ip));
      MELAout
        << "\t\tAPart[" << ip << "] (" << part->first << ") (X,Y,Z,T)=( "
        << part->second.X() << " , "
        << part->second.Y() << " , "
        << part->second.Z() << " , "
        << part->second.T() << " )" << endl;
    }

    MELAout << "\tHas " << cand->pTopDaughters.size() << " unstable tops" << endl;
    for (unsigned int ip=0; ip<cand->pTopDaughters.size(); ip++){
      MELAout << "\t\tTop[" << ip << "] daughters:" << endl;
      for (unsigned int jp=0; jp<cand->pTopDaughters.at(ip).size(); jp++){
        SimpleParticle_t* part = &(cand->pTopDaughters.at(ip).at(jp));
        MELAout
          << "\t\t- Top daughter[" << ip << jp << "] (" << part->first << ") (X,Y,Z,T)=( "
          << part->second.X() << " , "
          << part->second.Y() << " , "
          << part->second.Z() << " , "
          << part->second.T() << " )" << endl;
      }
    }
    MELAout << "\tHas " << cand->pAntitopDaughters.size() << " unstable antitops" << endl;
    for (unsigned int ip=0; ip<cand->pAntitopDaughters.size(); ip++){
      MELAout << "\t\tAntitop[" << ip << "] daughters:" << endl;
      for (unsigned int jp=0; jp<cand->pAntitopDaughters.at(ip).size(); jp++){
        SimpleParticle_t* part = &(cand->pAntitopDaughters.at(ip).at(jp));
        MELAout
          << "\t\t- Antitop daughter[" << ip << jp << "] (" << part->first << ") (X,Y,Z,T)=( "
          << part->second.X() << " , "
          << part->second.Y() << " , "
          << part->second.Z() << " , "
          << part->second.T() << " )" << endl;
      }
    }
  }
}
