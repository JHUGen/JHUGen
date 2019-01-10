#include "boost.h"
#include <stdlib.h>




// boostMom is the mom.vec. to be boosted into the direction of mom.vec. refMom (with mass refMass)
inline void boost_(double *boostMom, double *refMom, double *refMass){
double EuklidSP,spacialFact;

   EuklidSP    = refMom[1]*boostMom[1] + refMom[2]*boostMom[2] + refMom[3]*boostMom[3];
   spacialFact = (EuklidSP/((refMom[0]) + (*refMass)) + boostMom[0])/(*refMass);
   boostMom[0] = (refMom[0]*boostMom[0] + EuklidSP)/(*refMass);
   boostMom[1] = boostMom[1] + refMom[1]*spacialFact;
   boostMom[2] = boostMom[2] + refMom[2]*spacialFact;
   boostMom[3] = boostMom[3] + refMom[3]*spacialFact;
};


