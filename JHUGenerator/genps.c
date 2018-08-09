#include "genps.h"
#include <stdlib.h>

#define CHECK_OSREL 0
#define CHECK_EMCON 0
#define OSREL_PREC 1e-13
#define EMCON_PREC 1e-13


inline double lsp(double *Mom1, double *Mom2){
   return Mom1[0]*Mom2[0] - Mom1[1]*Mom2[1] - Mom1[2]*Mom2[2] - Mom1[3]*Mom2[3];
};


inline void swapmom_(double *Mom1, double *Mom2){
double tmp[4];

   tmp[0]=Mom1[0];
   tmp[1]=Mom1[1];
   tmp[2]=Mom1[2];
   tmp[3]=Mom1[3];

   Mom1[0]=Mom2[0];
   Mom1[1]=Mom2[1];
   Mom1[2]=Mom2[2];
   Mom1[3]=Mom2[3];

   Mom2[0]=tmp[0];
   Mom2[1]=tmp[1];
   Mom2[2]=tmp[2];
   Mom2[3]=tmp[3];
};

void swapdbl_(double *r1, double *r2){ //"inline" taken out due to error
double tmp;

   tmp=(*r1);
   (*r1)=(*r2);
   (*r2)=tmp;
};




void gensing_(int *NOut, double *Energy, double *Mass, double Mom[][4],
              int *Pcol1, int *Pcol2, double *SingDepth, int *NGrid){
static double *XRan;
double PSWgt,xJac,MInv1,mu=0.0;
double SingScale,eta;
const int Nminus1=(*NOut)-1;
const int Nminus2=(*NOut)-2;
int i;
static int ncall=0;
//  conventions for Pcol1, Pcol2:   ( i,j = 1 or 2 )
//  if Pcoli=0,1 then Pcolj becomes collinear to the beam axis
//  if Pcoli,Pcolj >=2 then Pcoli becomes collinear to Pcolj, Pcol2 denotes the first outgoing momentum
//  if Pcol1==Pcol2 then Pcol1 becomes soft


    if ( (*Pcol1)>(*Pcol2) ){ int tmp=(*Pcol1);(*Pcol1)=(*Pcol2);(*Pcol2)=tmp; }

// sort masses
    swapdbl_(&Mass[0],&Mass[(*Pcol2)-2]);
    if ( (*Pcol1)>1 && (*Pcol1)!=(*Pcol2) ) swapdbl_(&Mass[1],&Mass[(*Pcol1)-2]);

// initialize
   if ( ncall==0 ){
      XRan = malloc( (3*(*NOut)-4+1)*sizeof(double) );
      srand(time(NULL));
      for(i=0; i<=3*(*NOut)-4; i++) XRan[i]=rand()/((double)RAND_MAX + 1);
   }
   for(i=1; i<=Nminus1; i++ ) mu+=Mass[i];
   xJac = (*Energy)-Mass[0] - mu;
   MInv1 = xJac*XRan[0] + mu;


   SingScale = pow(10.0, log10(0.1) + ncall*( log10((*SingDepth))-log10(0.1) )/(*NGrid) );
   ncall++;
// generate collinear limit
   if ( (*Pcol1)!=(*Pcol2) ){
      if ( (*Pcol1)==0 ) XRan[Nminus1] = 1.0 - 0.5*SingScale;
      if ( (*Pcol1)==1 ) XRan[Nminus1] = 0.0 + 0.5*SingScale;
      if ( (*Pcol1)>1  ) {
         eta=sqrt(SingScale*MInv1/(*Energy)*2.0*(1.0-XRan[Nminus1])
            /(16.0*sq(M_PI*XRan[Nminus2]*(1.0-XRan[Nminus1]))+1.0)/XRan[Nminus1]);
         eta=sqrt(SingScale*2.0*(1.0-XRan[Nminus1])/(16.0*sq(M_PI*XRan[Nminus2]*(1.0-XRan[Nminus1]))+1.0)/XRan[Nminus1]);
         XRan[Nminus2+2]=(1.0-eta)*XRan[Nminus2];
         XRan[Nminus1+2]=(1.0-eta)*XRan[Nminus1];
      };
   };
// generate soft limit
   if ( (*Pcol1)==(*Pcol2) ){
      XRan[0]=(*Energy)/xJac*( sqrt(1.0-2.0*SingScale)-mu/(*Energy) );
   };


// generate phase space
   genps_(NOut,Energy,XRan,Mass,Mom,&PSWgt);

// sort masses
    swapdbl_(&Mass[0],&Mass[(*Pcol2)-2]);
    if ( (*Pcol1)>1 && (*Pcol1)!=(*Pcol2)) swapdbl_(&Mass[1],&Mass[(*Pcol1)-2]);
// sort momenta
    swapmom_(Mom[0],Mom[(*Pcol2)-2]);
    if ( (*Pcol1)>1 && (*Pcol1)!=(*Pcol2) ) swapmom_(Mom[1],Mom[(*Pcol1)-2]);
};






// TO DO:  special psgen for 2 particle final state,
//         generate initial state momenta
// inline void rocky_(int *NOut, double *Energy, double *XRan, double *Mass, double Mom[][4], double *PSWgt ){
//
// //    if ( (*NOut)==2 ){ genps2out_(Energy,XRan,Mass,Mom,PSWgt);  }
// //    else             { genps_(NOut,Energy,XRan,Mass,Mom,PSWgt); }
//
//    genps_(NOut,Energy,XRan,Mass,Mom,PSWgt);
//
// };






// generate phase space (using sequential splitting)
// missing normalization (Byckling norm * linear mapping)= (2*Pi)^-(3n-4) * (4*Pi)^(n-1)
void genps_(int *NOut, double *Energy, double *XRan, double *Mass, double Mom[][4], double *PSWgt ){
double MInv[10];
double EIntMom[10];
double xJac,mu=0.0;
const int Nminus1=(*NOut)-1;
const int Nminus2=(*NOut)-2;
int Levl,NMom;
double SingScale,eta;

// set invariant masses: MInv
// set Jacobian factors: xJac
   MInv[0] = (*Energy);
   MInv[Nminus1] = Mass[Nminus1];
   (*PSWgt)=0.5/MInv[0];

   for(NMom=1; NMom<=Nminus1; NMom++ ) mu+=Mass[NMom];
   for(Levl=1; Levl<=Nminus2; Levl++ ) {
      xJac = MInv[Levl-1]-Mass[Levl-1] - mu;
      MInv[Levl] = xJac*XRan[Levl-1] + mu;
      mu-=Mass[Levl];
      (*PSWgt)*= xJac;
   }


// split MInv[Levl] --> MInv[Levl+1],Mass[Levl]
   for(Levl=0; Levl<=Nminus2; Levl++ ) {
      (*PSWgt)*= split( &MInv[Levl], &XRan[Nminus2+2*Levl], &Mass[Levl], &MInv[Levl+1], Mom[Levl] );
      EIntMom[Levl] = MInv[Levl]-Mom[Levl][0];
   }


// set last momentum vector
   Mom[Nminus1][0] =  EIntMom[Nminus2];
   Mom[Nminus1][1] = -Mom[Nminus2][1];
   Mom[Nminus1][2] = -Mom[Nminus2][2];
   Mom[Nminus1][3] = -Mom[Nminus2][3];

// boost momenta to partonic CMS
   for(Levl=Nminus2; Levl>=1; Levl-- ) {
      if( MInv[Levl]<=1e-14 ){ (*PSWgt)*= 0.0; return; }    // NEW!
      for(NMom=Levl; NMom<=Nminus1; NMom++ ) {
          boostback( Mom[NMom],Mom[Levl-1],&MInv[Levl],&EIntMom[Levl-1] );
      }
   }


#if CHECK_OSREL
   checkOSrel(NOut,Energy,Mass,Mom,OSREL_PREC);
#endif
#if CHECK_EMCON
   checkEMcon(NOut,Energy,Mass,Mom,EMCON_PREC);
#endif

};





static inline double split(double *Energy, double *XRan, double *Mass, double *MInv, double *Mom){
double cosT,sinT,cosP,sinP;
double absVec;
const double TwoPi=2.0*M_PI;


   if( (*Energy) <= 1e-14 ){ return 0.0; }   // NEW!
   cosT  = 2.0*XRan[1]-1.0;
   sinT = 2.0*sqrt( XRan[1]*(1.0-XRan[1]) );

   //sincos(TwoPi*XRan[0], &sinP, &cosP);
   sinP=sin(TwoPi*XRan[0]);
   cosP=cos(TwoPi*XRan[0]);


   absVec = kappa( *Energy, (*Mass), *MInv );

   Mom[0] = sqrt( sq((*Mass)) + sq(absVec) );
   Mom[1] = absVec * ( cosP * sinT );
   Mom[2] = absVec * ( sinP * sinT );
   Mom[3] = absVec * (     cosT    );

   return 0.5*absVec;
};




// refMass,refEnergy: corresp. to reference vector
// refMomX: second vector from splitting, oppsosite to reference vector but with different energy
inline void boostback(double *boostMom, double *refMomX, double *refMass, double *refEnergy){
double EuklidSP,spacialFact;

   EuklidSP    = -refMomX[1]*boostMom[1] - refMomX[2]*boostMom[2] - refMomX[3]*boostMom[3];
   spacialFact = (EuklidSP/((*refEnergy) + (*refMass)) + boostMom[0])/(*refMass);
   boostMom[0] = ((*refEnergy)*boostMom[0] + EuklidSP)/(*refMass);
   boostMom[1] = boostMom[1] - refMomX[1]*spacialFact;
   boostMom[2] = boostMom[2] - refMomX[2]*spacialFact;
   boostMom[3] = boostMom[3] - refMomX[3]*spacialFact;

};





// kappa(x,y,z) = 1/(2*x) * SQRT( lambda(x^2,y^2,z^2) )
static inline double kappa(double x, double y, double z){
double la;

   if( x==0.0 ) {return 0.0;}
   else {
          la = sq(x) -2.0*(sq(y)+sq(z)) + sq((sq(y)-sq(z))/x);
          if ( fabs(la/x) <=1e-13 ) {return 0.0;}
          else return 0.5*sqrt( la );
   }
};


static inline double sq(double x){ return x*x; }


static void checkOSrel(int *NOut, double *Energy, double *Mass, double Mom[][4], double limit){
double OSrel;
int i;

   for(i=0; i<=(*NOut)-1; i++ ){
      OSrel = (sq(Mass[i])-(sq(Mom[i][0])-sq(Mom[i][1])-sq(Mom[i][2])-sq(Mom[i][3])))/sq(*Energy);
      if( (fabs(OSrel)>(limit)) ) {
              printf("OS relation for p%d violated by: %4.2e \n",i,OSrel);
      }
   }
}


static void checkEMcon(int *NOut, double *Energy, double *Mass, double Mom[][4], double limit){
double EMcon[4]={(*Energy), 0.0, 0.0, 0.0};
int i;

   for(i=0; i<=(*NOut)-1; i++ ){
      EMcon[0] -= Mom[i][0];
      EMcon[1] += Mom[i][1];
      EMcon[2] += Mom[i][2];
      EMcon[3] += Mom[i][3];
   }
   EMcon[0]=EMcon[0]/(*Energy);
   EMcon[1]=EMcon[1]/(*Energy);
   EMcon[2]=EMcon[2]/(*Energy);
   EMcon[3]=EMcon[3]/(*Energy);

   if( (fabs(EMcon[0])>(limit)) || (fabs(EMcon[1])>(limit)) || (fabs(EMcon[2])>(limit)) || (fabs(EMcon[3])>(limit)) ){
      printf("violated EM conservation: \n");
      printf("E: %4.2e \n",EMcon[0]);
      printf("x: %4.2e \n",EMcon[1]);
      printf("y: %4.2e \n",EMcon[2]);
      printf("z: %4.2e \n",EMcon[3]);
   };
};

