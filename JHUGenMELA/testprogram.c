#include <stdio.h>
#include <math.h>
#include <stdlib.h>



int main(void){
// input unit = GeV/100 such that 125GeV is 1.25 in the code
  double MReso = 125.0/100.0;
  double GaReso= 0.1/100.0;
  double P[6][4];
  double MatElSq;
  int MYIDUP[4];
  double Hggcoupl[3][2];
  double Hvvcoupl[4][2];
  double Zqqcoupl[2][2];
  double Zvvcoupl[2][2];
  double Gqqcoupl[2][2];
  double Gggcoupl[5][2];
  double Gvvcoupl[10][2];

  Hggcoupl[0][0]=1.0;  Hggcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
  Hggcoupl[1][0]=0.0;  Hggcoupl[1][1]=0.0;  
  Hggcoupl[2][0]=0.0;  Hggcoupl[2][1]=0.0;  
  Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0;  
  Hvvcoupl[1][0]=0.0;  Hvvcoupl[1][1]=0.0;  
  Hvvcoupl[2][0]=0.0;  Hvvcoupl[2][1]=0.0;  
  Hvvcoupl[3][0]=0.0;  Hvvcoupl[3][1]=0.0;  

  Zqqcoupl[0][0]=1.0;  Zqqcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
  Zqqcoupl[1][0]=1.0;  Zqqcoupl[1][1]=0.0;
  
  Zvvcoupl[0][0]=0.0;  Zvvcoupl[0][1]=0.0;
  Zvvcoupl[1][0]=1.0;  Zvvcoupl[1][1]=0.0;
  
  Gqqcoupl[0][0]=1.0;  Gqqcoupl[0][1]=0.0;  
  Gqqcoupl[1][0]=1.0;  Gqqcoupl[1][1]=0.0; 
  Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0;
  Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
  Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
  Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
  Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;
  Gvvcoupl[0][0]=1.0;  Gvvcoupl[0][1]=0.0;
  Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
  Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
  Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
  Gvvcoupl[4][0]=1.0;  Gvvcoupl[4][1]=0.0;
  Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
  Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
  Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
  Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
  Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

  
// particle ID: +7=e+,  -7=e-,  +8=mu+,  -8=mu-
  MYIDUP[0]=+7;
  MYIDUP[1]=-7;
  MYIDUP[2]=+7;
  MYIDUP[3]=-7;
  
// p(i,0:3) = (E(i),px(i),py(i),pz(i)) in units of 100GeV
// i=0,1: glu1,glu2 (outgoing convention)
// i=2,3: correspond to MY_IDUP(1),MY_IDUP(0)
// i=4,5: correspond to MY_IDUP(3),MY_IDUP(2)
  P[0][0]=-0.0645264033200954;
  P[0][1]=0.0;
  P[0][2]=0.0;
  P[0][3]=-0.0645264033200954;
  
  P[1][0]=-6.0537234356954572;
  P[1][1]=0.0;
  P[1][2]=0.0;
  P[1][3]=6.0537234356954572;
  
  P[2][0]=4.7598878302889442;
  P[2][1]=0.4544925597087586;
  P[2][2]=0.0917597774970785;
  P[2][3]=-4.7372511874858247;
  
  P[3][0]=0.8177159853499207;
  P[3][1]=-0.2768889802512220;
  P[3][2]=-0.0202805643015353;
  P[3][3]=-0.7691427851991084;
  
  P[4][0]=0.3191706236713395;
  P[4][1]=-0.0470102494218651;
  P[4][2]=-0.0079466854602927;
  P[4][3]=-0.3155895651859249;
  
  P[5][0]=0.2214753997053482;
  P[5][1]=-0.1305933300356715;
  P[5][2]=-0.0635325277352505;
  P[5][3]=-0.1672134945045036;
  

 __modhiggs_MOD_evalamp_gg_h_vv(P, &MReso,  &GaReso, Hggcoupl, Hvvcoupl, MYIDUP, &MatElSq);
 printf("\n ");
 printf("Matr.el. squared (spin-0): %20.17e \n ",MatElSq);
 printf("result should be (spin-0): %20.17e \n ",0.0045682366425370);
 printf("ratio: %20.17e \n ",MatElSq/0.0045682366425370);

 __modzprime_MOD_evalamp_qqb_zprime_vv(P, &MReso,  &GaReso, Zqqcoupl, Zvvcoupl, MYIDUP, &MatElSq);
 printf("\n ");
 printf("Matr.el. squared (spin-1): %20.17e \n ",MatElSq);
 printf("result should be (spin-1): %20.17e \n ",0.0020357978978982);
 printf("ratio: %20.17e \n ",MatElSq/0.0020357978978982);

  __modgraviton_MOD_evalamp_gg_g_vv(P, &MReso,  &GaReso, Gggcoupl, Gvvcoupl, MYIDUP, &MatElSq);
 printf("\n ");
 printf("Matr.el. squared (gg spin-2): %20.17e \n ",MatElSq);
 printf("result should be (gg spin-2): %20.17e \n ",0.0307096869320374);
 printf("ratio: %20.17e \n ",MatElSq/0.0307096869320374);

  __modgraviton_MOD_evalamp_qqb_g_vv(P, &MReso,  &GaReso, Gqqcoupl, Gvvcoupl, MYIDUP, &MatElSq);
 printf("\n ");
 printf("Matr.el. squared (qq spin-2): %20.17e \n ",MatElSq);
 printf("result should be (qq spin-2): %20.17e \n ",0.0004838377647021);
 printf("ratio: %20.17e \n ",MatElSq/0.0004838377647021);






  P[0][0]=-1.25;
  P[0][1]=0.0;
  P[0][2]=0.0;
  P[0][3]=0.0;
  
  P[1][0]=0.0;
  P[1][1]=0.0;
  P[1][2]=0.0;
  P[1][3]=0.0;
  
  P[2][0]=0.5723124900092045;
  P[2][1]=0.2083223910863685;
  P[2][2]=0.1441561361468674;
  P[2][3]=-0.5131884410270751;
  
  P[3][0]=0.3270364795109065;
  P[3][1]=-0.1343151455361828;
  P[3][2]=0.0368904244639839;
  P[3][3]=0.2958908535141779;
  
  P[4][0]=0.0584259632074581;
  P[4][1]=-0.0206108978979021;
  P[4][2]=-0.0424960639382021;
  P[4][3]=-0.0343928570247043;
  
  P[5][0]=0.2922250672724307;
  P[5][1]=-0.0533963476522836;
  P[5][2]=-0.1385504966726492;
  P[5][3]=0.2516904445376015;


 __modzprime_MOD_evalamp_zprime_vv(P, &MReso,  &GaReso, Zvvcoupl, MYIDUP, &MatElSq);
 printf("\n ");
 printf("no production dynamics\n ");
 printf("Matr.el. squared (spin-1): %20.17e \n ",MatElSq);
 printf("result should be (spin-1): %20.17e \n ",5.25259158248293146e-10);
 printf("ratio: %20.17e \n ",MatElSq/5.25259158248293146e-10);

  __modgraviton_MOD_evalamp_g_vv(P, &MReso,  &GaReso, Gvvcoupl, MYIDUP, &MatElSq);
 printf("\n ");
 printf("no production dynamics\n ");
 printf("Matr.el. squared (spin-2): %20.17e \n ",MatElSq);
 printf("result should be (spin-2): %20.17e \n ",3.50330723427981412e-09);
 printf("ratio: %20.17e \n ",MatElSq/3.50330723427981412e-09);




 
return 0;
};

