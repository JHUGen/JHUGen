* eval propagator functions
id PropDenom(pDumX?,mX?) = PropDenom(pDumX.pDumX-mX^2);
argument PropDenom;
#call simplify();
endargument;



*#call encapsulateColor();
*id Color(?args) = 1;

#call colorAlgebra2();
id SUND(Col2,Col1) = SUND(Col1,Col2);
id SUND(Col1,Col2) = ICol(1);




* expand Ga(pi-pj)=Ga(pi)-Ga(pj)
id Ga(DumStr1?,LorX?)=Ga(DumStr1,LorX);



* move Chir to the left
#call moveChirL();
#call simplify();




* Dirac algebra
repeat;
   #call applyDiracEq(1);
   #call commuteToLeft(p2);
   #call commuteToRight(p1);
endrepeat;
#call simplify();




* apply transversality condition for polarization vectors
#call applyTransversality();




* pull out loop momenta q1
id once Ga(DumStr1?,q1) = Ga(DumStr1,DumLor1)*q1(DumLor1);
id once Ga(DumStr1?,q1) = Ga(DumStr1,DumLor2)*q1(DumLor2);
id once Ga(DumStr1?,q1) = Ga(DumStr1,DumLor3)*q1(DumLor3);
id once Ga(DumStr1?,q1) = Ga(DumStr1,DumLor4)*q1(DumLor4);
id once pDumX?FVec.q1 = pDumX(DumLor5)*q1(DumLor5);
id once pDumX?FVec.q1 = pDumX(DumLor6)*q1(DumLor6);
id once pDumX?FVec.q1 = pDumX(DumLor7)*q1(DumLor7);
id once pDumX?FVec.q1 = pDumX(DumLor8)*q1(DumLor8);
id once LeviCiv(q1,LorX?,LorY?,LorZ?) = LeviCiv(DumLor9,LorX,LorY,LorZ)*q1(DumLor9);


* identify tensor integrals (Denner's conventions)
multiply SIntDummy;
id q1(LorW?DumId)*q1(LorX?DumId)*q1(LorY?DumId)*q1(LorZ?DumId) = LoopMom(LorW,LorX,LorY,LorZ)/SIntDummy;
id q1(LorW?DumId)*q1(LorX?DumId)*q1(LorY?DumId)                = LoopMom(LorW,LorX,LorY)     /SIntDummy;
id q1(LorW?DumId)*q1(LorX?DumId)                               = LoopMom(LorW,LorX)          /SIntDummy;
id q1(LorW?DumId)                                              = LoopMom(LorW)               /SIntDummy;
id SIntDummy                                                   = LoopMom(0);


id LoopDenom( q1,m0?, pDumW?,m1? ) * LoopMom(?args)                   
   = i_*Pi^2 * TI(2,?args,pDumW-q1,m0,m1);
id LoopDenom( q1,m0?, pDumW?,m1?, pDumX?,m2? ) * LoopMom(?args) 
   = i_*Pi^2 * TI(3,?args,pDumW-q1,pDumX-q1,m0,m1,m2);
id LoopDenom( q1,m0?, pDumW?,m1?, pDumX?,m2?, pDumY?,m3? ) * LoopMom(?args) 
   = i_*Pi^2 * TI(4,?args,pDumW-q1,pDumX-q1,pDumY-q1,m0,m1,m2,m3);

id TI(2,0,pDumW?,m0?,m1?) 
 = SI(2,pDumW,m0,m1);
id TI(3,0,pDumW?,pDumX?,m0?,m1?,m2?) 
 = SI(3,pDumW,pDumX,m0,m1,m2);
id TI(4,0,pDumW?,pDumX?,pDumY?,m0?,m1?,m2?,m3?) 
 = SI(4,pDumW,pDumX,pDumY,m0,m1,m2,m3);

* remove LoopMom(0) from the LO contribution
id LoopMom(0)=1;



* insert tensor decomposition ( scalar integrals must be identified first )
#include  /afs/cern.ch/user/m/maschulz/lib/FeynArtsToForm/LorDec.frm




sum DumLor1;
sum DumLor2;
sum DumLor3;
sum DumLor4;
sum DumLor5;
sum DumLor6;
sum DumLor7;
sum DumLor8;
sum DumLor9;

#call simplify();
.sort 





#call FORMTrace(2);


id Ga(1,p4) = Ga(1,p1) + Ga(1,p2) - Ga(1,p3);

* Dirac algebra
repeat;
   #call applyDiracEq(1);
   #call commuteToLeft(p2);
   #call commuteToRight(p1);
endrepeat;
#call simplify();






id i_=cI;
argument;
   id i_=cI;
endargument;



***#call insertMandelstam2to2();



id MW = CW*MZ;
id 1/MW = 1/CW/MZ;


id ASpi(1,p2,0)*Chir(1,+1)*Ga(1,cep3)*Spi(1,p1,0) = SME(1,+1);
id ASpi(1,p2,0)*Chir(1,-1)*Ga(1,cep3)*Spi(1,p1,0) = SME(1,-1);

id ASpi(1,p2,0)*Chir(1,+1)*Ga(1,p3)*Spi(1,p1,0) = SME(2,+1);
id ASpi(1,p2,0)*Chir(1,-1)*Ga(1,p3)*Spi(1,p1,0) = SME(2,-1);






