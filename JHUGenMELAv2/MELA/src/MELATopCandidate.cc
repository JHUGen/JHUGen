#include "MELATopCandidate.h"

using namespace std;

MELATopCandidate::MELATopCandidate(
  MELAParticle* lightQuark_,
  MELAParticle* Wferm_,
  MELAParticle* Wfermbar_
  ) : MELAParticle(),
  lightQuark(lightQuark_),
  Wferm(Wferm_),
  Wfermbar(Wfermbar_)
{
  if (lightQuark!=0){
    p4 = p4 + lightQuark->p4;
    addDaughter(lightQuark);
  }
  double Wcharge=0;
  if (Wferm!=0){
    p4 = p4 + Wferm->p4;
    Wcharge += Wferm->charge();
    addDaughter(Wferm);
  }
  if (Wfermbar!=0){
    p4 = p4 + Wfermbar->p4;
    Wcharge += Wfermbar->charge();
    addDaughter(Wfermbar);
  }
  if (
    (lightQuark!=0 && PDGHelpers::isDownTypeQuark(lightQuark->id) && lightQuark->id>0)
    ||
    (Wcharge>0)
    ) id = 6;
  else if (
    (lightQuark!=0 && PDGHelpers::isDownTypeQuark(lightQuark->id) && lightQuark->id<0)
    ||
    (Wcharge<0)
    ) id = -6;
  else id=0;
}

void MELATopCandidate::setLightQuark(MELAParticle* myParticle){ lightQuark=myParticle; if (lightQuark!=0) addDaughter(lightQuark); }
void MELATopCandidate::setWFermion(MELAParticle* myParticle){ Wferm=myParticle; if (Wferm!=0) addDaughter(Wferm); }
void MELATopCandidate::setWAntifermion(MELAParticle* myParticle){ Wfermbar=myParticle; if (Wfermbar!=0) addDaughter(Wfermbar); }

