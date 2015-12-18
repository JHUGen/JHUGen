

* eval propagator functions
id PropDenom(pDumX?,mX?) = PropDenom(pDumX.pDumX-mX^2);
argument PropDenom;
#call simplify();
endargument;


* remove sum over internal colors
id SumOver(?args) = 1;




* expand Ga(pi-pj)=Ga(pi)-Ga(pj)
id Ga(DumStr1?,LorX?)=Ga(DumStr1,LorX);


* move Chir to the left
#call moveChirL();
#call simplify();



id Ga(1,p4) = Ga(1,p1) + Ga(1,p2) - Ga(1,p3);



* Dirac algebra
repeat;
   #call applyDiracEq(1);
   #call commuteToLeft(p2);
   #call commuteToRight(p1);
endrepeat;
#call simplify();


* move Ga(ep3) to the left
repeat;
   #call commuteToLeft(cep3);
   #call simplify();
   #call applyTransversality();
endrepeat;


* apply transversality condition for polarization vectors
#call applyTransversality();


* move Chir to the left
#call moveChirR();


#call simplify();


id i_=cI;
argument;
   id i_=cI;
endargument;


id MW = CW*MZ;
id 1/MW = 1/CW/MZ;

id ASpi(1,p2,0)*Ga(1,cep3)*Chir(1,+1)*Spi(1,p1,0) = SME(1,+1);
id ASpi(1,p2,0)*Ga(1,cep3)*Chir(1,-1)*Spi(1,p1,0) = SME(1,-1);

id ASpi(1,p2,0)*Ga(1,p3)*Chir(1,+1)*Spi(1,p1,0)*p4.cep3 = SME(2,+1);
id ASpi(1,p2,0)*Ga(1,p3)*Chir(1,-1)*Spi(1,p1,0)*p4.cep3 = SME(2,-1);

id ASpi(1,p2,0)*Ga(1,Lor2)*Chir(1,+1)*Spi(1,p1,0)*LeviCiv(p3,p4,cep3,Lor2) = SME(3,+1);
id ASpi(1,p2,0)*Ga(1,Lor2)*Chir(1,-1)*Spi(1,p1,0)*LeviCiv(p3,p4,cep3,Lor2) = SME(3,-1);





