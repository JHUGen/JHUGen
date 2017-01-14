#include "Mela.h"
#include "newZZMatrixElement.h"

#ifdef __CINT__


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ class MELAParticle;
#pragma link C++ class std::vector<MELAParticle*>;
#pragma link C++ class MELATopCandidate;
#pragma link C++ class MELACandidate;
#pragma link C++ class MelaIO;
#pragma link C++ namespace TVar;

#pragma link C++ namespace TUtil;
#pragma link C++ function TUtil::computeAngles;
#pragma link C++ function TUtil::computeAnglesCS;
#pragma link C++ function TUtil::computeVBFangles;
#pragma link C++ function TUtil::computeVHangles;
#pragma link C++ function TUtil::scaleMomentumToEnergy;
#pragma link C++ function TUtil::constrainedRemovePairMass;
#pragma link C++ function TUtil::removeMassFromPair;
#pragma link C++ function TUtil::adjustTopDaughters;
#pragma link C++ function TUtil::computeFakeJet;
//
#pragma link C++ function TUtil::SetAlphaS;
#pragma link C++ function TUtil::CheckPartonMomFraction;
#pragma link C++ function TUtil::ComputePDF;
#pragma link C++ function TUtil::SumMEPDF;
//
#pragma link C++ function TUtil::GetBoostedParticleVectors;
#pragma link C++ function TUtil::ConvertVectorFormat;
#pragma link C++ function TUtil::ConvertTopCandidate;
//
#pragma link C++ function TUtil::SumMatrixElementPDF;
#pragma link C++ function TUtil::JHUGenMatEl;
#pragma link C++ function TUtil::HJJMatEl;
#pragma link C++ function TUtil::VHiggsMatEl;
#pragma link C++ function TUtil::TTHiggsMatEl;
#pragma link C++ function TUtil::BBHiggsMatEl;

#pragma link C++ class SpinZeroCouplings;
#pragma link C++ class SpinOneCouplings;
#pragma link C++ class SpinTwoCouplings;
#pragma link C++ class TEvtProb;
#pragma link C++ class Mela;
#pragma link C++ class newZZMatrixElement;


#endif
