#include <iostream>
using namespace std;

// Example, run in ROOT:
// [] .L analyticPartonLuminosities.C
// [] main()

double pdfs_qqbar(double m, double Y, double sqrts);
double pdfs_gg(double m, double Y, double sqrts);
double pdfs_qqbarPrime(double m, double Y, double sqrts);

int main(){

    std::cout << "hello world" << std::endl;
    sqrtS0 = 14000;
    
    double testPoint1_qq  = pdfs_qqbar(10, 0, sqrtS0);
    double testPoint1_gg  = pdfs_gg(10, 0, sqrtS0);
    double testPoint1_qqp = pdfs_qqbarPrime(10, 0, sqrtS0);

    double testPoint2_qq  = pdfs_qqbar(100, 0, sqrtS0);
    double testPoint2_gg  = pdfs_gg(100, 0, sqrtS0);
    double testPoint2_qqp = pdfs_qqbarPrime(100, 0, sqrtS0);

    double testPoint3_qq  = pdfs_qqbar(1000, 0, sqrtS0);
    double testPoint3_gg  = pdfs_gg(1000, 0, sqrtS0);
    double testPoint3_qqp = pdfs_qqbarPrime(1000, 0, sqrtS0);

    std::cout << "m = 10 GeV, Y = 0..." << std::endl;
    std::cout << "testPoint1_qq = " << testPoint1_qq << ", testPoint1_gg = " << testPoint1_gg << ", testPoint1_qqp = " << testPoint1_qqp << std::endl;
    std::cout << "m = 100 GeV, Y = 0..." << std::endl;
    std::cout << "testPoint2_qq = " << testPoint2_qq << ", testPoint2_gg = " << testPoint2_gg << ", testPoint2_qqp = " << testPoint2_qqp << std::endl;
    std::cout << "m = 1000 GeV, Y = 0..." << std::endl;
    std::cout << "testPoint3_qq = " << testPoint3_qq << ", testPoint3_gg = " << testPoint3_gg << ", testPoint3_qqp = " << testPoint3_qqp << std::endl;
    
    return 0;
    
}



////////////////////////////////////////////////////////////////////////////
// Q-QBAR
////////////////////////////////////////////////////////////////////////////
double pdfs_qqbar(double m, double Y, double sqrts)
{
    
    double YVal = Y;
    double sqrtsVal = sqrts;
    double mVal = m;
    
    double weightu = 0.5;
	double weightd = 0.5;
	double weightc = 1.0;
	double weights = 1.0;
	double weightb = 1.0;
    
    double s0 = sqrtsVal*sqrtsVal;
	double s = mVal*mVal;
	double Q = mVal;
	double xa = exp(YVal)*sqrt(s/s0);
	double xb = exp(-YVal)*sqrt(s/s0);
    
	// PDF parameters
	// up params
	double u0par0 = 0.03134; double u0par1 =-2.068e-05; double u0par2 = 1.283e-08; 
	double u1par0 = 0.9; double u1par1 =-0.0004307; double u1par2 = 2.458e-07;
	double u2par0 =-0.1369; double u2par1 = 0.003423; double u2par2 =-2.155e-06;
	double u3par0 =-0.4013; double u3par1 =-0.0002574; double u3par2 = 1.561e-07;
	double u4par0 = 0.5782; double u4par1 =-0.004728; double u4par2 = 2.906e-06;
	double ubar0par0 = 0.02856; double ubar0par1 =-2.112e-05; double ubar0par2 = 1.272e-08;
	double ubar1par0 =-0.06822; double ubar1par1 = 3.172e-05; double ubar1par2 =-2.008e-08;
	double ubar2par0 = 0.1967; double ubar2par1 =-0.000118; double ubar2par2 = 6.871e-08;
	double ubar3par0 =-0.2251; double ubar3par1 = 0.0001295; double ubar3par2 =-7.181e-08;
	double ubar4par0 =-0.4068; double ubar4par1 =-0.0002956; double ubar4par2 = 1.783e-07;
	double ubar5par0 =-2.251; double ubar5par1 =-0.0001699; double ubar5par2 = 1.492e-07;
	// down params
	double d0par0 = 0.03278; double d0par1 =-2.915e-05; double d0par2 = 1.809e-08;
	double d1par0 = 0.479; double d1par1 =-0.0002559; double d1par2 = 1.557e-07;
	double d2par0 =-0.5972; double d2par1 = 0.0003118; double d2par2 =-1.905e-07;
	double d3par0 =-0.3892; double d3par1 =-0.000317; double d3par2 = 1.944e-07;
	double d4par0 = 0.5007; double d4par1 =-0.001665; double d4par2 = 9.895e-07;
	double dbar0par0 = 0.02328; double dbar0par1 =-1.367e-05; double dbar0par2 = 8.246e-09;
	double dbar1par0 = 0.09422; double dbar1par1 =-0.0001019; double dbar1par2 = 6.375e-08;
	double dbar2par0 =-0.5296; double dbar2par1 = 0.000466; double dbar2par2 =-2.896e-07;
	double dbar3par0 = 0.5354; double dbar3par1 =-0.0004404; double dbar3par2 = 2.728e-07;
	double dbar4par0 =-0.4386; double dbar4par1 =-0.0002605; double dbar4par2 = 1.582e-07;
	double dbar5par0 =-1.289; double dbar5par1 =-0.001618; double dbar5par2 = 9.601e-07;
	// charm, strange, bottom params
	double c0par0 = 0.01829; double c0par1 =-6.93e-06; double c0par2 = 3.796e-09;
	double c1par0 = 0.03081; double c1par1 = 4.325e-05; double c1par2 =-3.95e-08;
	double c2par0 = 0.5398; double c2par1 =-4.284e-05; double c2par2 =-1.362e-08;
	double c3par0 =-0.5986; double c3par1 = 0.002565; double c3par2 =-1.937e-06;
	double c4par0 =-0.4534; double c4par1 =-0.0002329; double c4par2 = 1.343e-07;
	double c5par0 =-8.657; double c5par1 =-0.005157; double c5par2 = 3.68e-06;
	double s0par0 = 0.01312; double s0par1 =-3.743e-06; double s0par2 = 2.076e-09;
	double s1par0 =-0.001416; double s1par1 =-7.649e-06; double s1par2 = 4.757e-09;
	double s2par0 = 0.2864; double s2par1 =-6.693e-05; double s2par2 = 3.566e-08;
	double s3par0 =-0.4857; double s3par1 =-0.000253; double s3par2 = 1.541e-07;
	double s4par0 =-10.33; double s4par1 =-0.001601; double s4par2 = 9.718e-07;
	double b0par0 = 0.005934; double b0par1 = 2.516e-06; double b0par2 =-1.828e-09;
	double b1par0 =-0.003063; double b1par1 =-6.761e-06; double b1par2 = 4.298e-09;
	double b2par0 = 0.1174; double b2par1 = 3.752e-05; double b2par2 =-2.863e-08;
	double b3par0 =-0.5549; double b3par1 =-0.0002205; double b3par2 = 1.334e-07;
	double b4par0 =-10.18; double b4par1 =-0.001136; double b4par2 = 6.931e-07;
    
	
	// PDF definition
	double up0 = u0par0 + u0par1*Q + u0par2*Q*Q;
	double up1 = u1par0 + u1par1*Q + u1par2*Q*Q;
	double up2 = u2par0 + u2par1*Q + u2par2*Q*Q;
	double up3 = u3par0 + u3par1*Q + u3par2*Q*Q;
	double up4 = u4par0 + u4par1*Q + u4par2*Q*Q;
	double antiup0 = ubar0par0 + ubar0par1*Q + ubar0par2*Q*Q;
	double antiup1 = ubar1par0 + ubar1par1*Q + ubar1par2*Q*Q;
	double antiup2 = ubar2par0 + ubar2par1*Q + ubar2par2*Q*Q;
	double antiup3 = ubar3par0 + ubar3par1*Q + ubar3par2*Q*Q;
	double antiup4 = ubar4par0 + ubar4par1*Q + ubar4par2*Q*Q;
	double antiup5 = ubar5par0 + ubar5par1*Q + ubar5par2*Q*Q;
	double down0 = d0par0 + d0par1*Q + d0par2*Q*Q;
	double down1 = d1par0 + d1par1*Q + d1par2*Q*Q;
	double down2 = d2par0 + d2par1*Q + d2par2*Q*Q;
	double down3 = d3par0 + d3par1*Q + d3par2*Q*Q;
	double down4 = d4par0 + d4par1*Q + d4par2*Q*Q;
	double antidown0 = dbar0par0 + dbar0par1*Q + dbar0par2*Q*Q;
	double antidown1 = dbar1par0 + dbar1par1*Q + dbar1par2*Q*Q;
	double antidown2 = dbar2par0 + dbar2par1*Q + dbar2par2*Q*Q;
	double antidown3 = dbar3par0 + dbar3par1*Q + dbar3par2*Q*Q;
	double antidown4 = dbar4par0 + dbar4par1*Q + dbar4par2*Q*Q;
	double antidown5 = dbar5par0 + dbar5par1*Q + dbar5par2*Q*Q;
	double charm0 = c0par0 + c0par1*Q + c0par2*Q*Q;
	double charm1 = c1par0 + c1par1*Q + c1par2*Q*Q;
	double charm2 = c2par0 + c2par1*Q + c2par2*Q*Q;
	double charm3 = c3par0 + c3par1*Q + c3par2*Q*Q;
	double charm4 = c4par0 + c4par1*Q + c4par2*Q*Q;
	double charm5 = c5par0 + c5par1*Q + c5par2*Q*Q;
	double strange0 = s0par0 + s0par1*Q + s0par2*Q*Q;
	double strange1 = s1par0 + s1par1*Q + s1par2*Q*Q;
	double strange2 = s2par0 + s2par1*Q + s2par2*Q*Q;
	double strange3 = s3par0 + s3par1*Q + s3par2*Q*Q;
	double strange4 = s4par0 + s4par1*Q + s4par2*Q*Q;
	double bottom0 = b0par0 + b0par1*Q + b0par2*Q*Q;
	double bottom1 = b1par0 + b1par1*Q + b1par2*Q*Q;
	double bottom2 = b2par0 + b2par1*Q + b2par2*Q*Q;
	double bottom3 = b3par0 + b3par1*Q + b3par2*Q*Q;
	double bottom4 = b4par0 + b4par1*Q + b4par2*Q*Q;
    
	double FuncAu1 = (up0+up1*xa+up2*pow(xa,2))*pow((1-xa),4)*pow(xa,up3)*exp(1.0+up4*xa);
	double FuncBu1 = (antiup0+antiup1*xb+antiup2*pow(xb,2)+antiup3*pow(xb,3))*pow((1-xb),4)*pow(xb,antiup4)*exp(1.0+antiup5*xb);
	double FuncAu2 = (up0+up1*xb+up2*pow(xb,2))*pow((1-xb),4)*pow(xb,up3)*exp(1.0+up4*xb);
	double FuncBu2 = (antiup0+antiup1*xa+antiup2*pow(xa,2)+antiup3*pow(xa,3))*pow((1-xa),4)*pow(xa,antiup4)*exp(1.0+antiup5*xa);
	double FuncABu = FuncAu1/xa*FuncBu1/xb+FuncAu2/xa*FuncBu2/xb;
    
	
	double FuncAd1 = (down0+down1*xa+down2*pow(xa,2))*pow((1-xa),4)*pow(xa,down3)*exp(1.0+down4*xa);
	double FuncBd1 = (antidown0+antidown1*xb+antidown2*pow(xb,2)+antidown3*pow(xb,3))*pow((1-xb),4)*pow(xb,antidown4)*exp(1.0+antidown5*xb);
	double FuncAd2 = (down0+down1*xb+down2*pow(xb,2))*pow((1-xb),4)*pow(xb,down3)*exp(1.0+down4*xb);
	double FuncBd2 = (antidown0+antidown1*xa+antidown2*pow(xa,2)+antidown3*pow(xa,3))*pow((1-xa),4)*pow(xa,antidown4)*exp(1.0+antidown5*xa);
	double FuncABd = FuncAd1/xa*FuncBd1/xb+FuncAd2/xa*FuncBd2/xb;
    
	double Funcca = (charm0+charm1*xa+charm2*pow(xa,2)+charm3*pow(xa,3))*pow((1-xa),4)*pow(xa,charm4)*exp(1.0+charm5*xa);
	double Funccb = (charm0+charm1*xb+charm2*pow(xb,2)+charm3*pow(xb,3))*pow((1-xb),4)*pow(xb,charm4)*exp(1.0+charm5*xb);
	double Funcsa = (strange0+strange1*xa+strange2*pow(xa,2))*pow((1-xa),4)*pow(xa,strange3)*exp(1.0+strange4*xa);
	double Funcsb = (strange0+strange1*xb+strange2*pow(xb,2))*pow((1-xb),4)*pow(xb,strange3)*exp(1.0+strange4*xb);
	double Funcba = (bottom0+bottom1*xa+bottom2*pow(xa,2))*pow((1-xa),4)*pow(xa,bottom3)*exp(1.0+bottom4*xa);
	double Funcbb = (bottom0+bottom1*xb+bottom2*pow(xb,2))*pow((1-xb),4)*pow(xb,bottom3)*exp(1.0+bottom4*xb);
	double FuncABc = Funcsa*Funcsb/xa/xb;
	double FuncABs = Funcca*Funccb/xa/xb;
	double FuncABb = Funcba*Funcbb/xa/xb;
    
    
	double totSec = 2*mVal*(
                              (FuncABu)*weightu
                              +(FuncABd)*weightd
                              +(FuncABc)*weightc
                              +(FuncABs)*weights
                              +(FuncABb)*weightb
                              );
    
	if(( mVal <= 600. && fabs(YVal) > 20*pow(float(mVal),float(-0.32))) || ( mVal > 600. && fabs(YVal) > 21*pow(float(mVal),float(-0.34))))
        {
	    //Find totSec when mZZ, YVal=0
	    double xa0 = sqrt(s/s0); //at YVal=0 xa=xb
        
	    //up
	    //if xa=xb then FuncAu1=FuncAu2 and FuncBu1=FuncBu2
	    FuncAu1 = (up0+up1*xa0+up2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,up3)*exp(1.0+up4*xa0);
	    FuncBu1 = (antiup0+antiup1*xa0+antiup2*pow(xa0,2)+antiup3*pow(xa0,3))*pow((1-xa0),4)*pow(xa0,antiup4)*exp(1.0+antiup5*xa0);
	    FuncABu = 2*(FuncAu1/xa0*FuncBu1/xa0);
        
	    //down
	    //if xa=xb then FuncAd1=FuncAd2 and FuncBd1=FuncBd2
	    FuncAd1 = (down0+down1*xa0+down2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,down3)*exp(1.0+down4*xa0);
	    FuncBd1 = (antidown0+antidown1*xa0+antidown2*pow(xa0,2)+antidown3*pow(xa0,3))*pow((1-xa0),4)*pow(xa0,antidown4)*exp(1.0+antidown5*xa0);
	    FuncABd = 2*(FuncAd1/xa0*FuncBd1/xa0);
        
	    //sea
	    Funcca = (charm0+charm1*xa0+charm2*pow(xa0,2)+charm3*pow(xa0,3))*pow((1-xa0),4)*pow(xa0,charm4)*exp(1.0+charm5*xa0); //Funcca=Funccb
	    Funcsa = (strange0+strange1*xa0+strange2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,strange3)*exp(1.0+strange4*xa0); //Funcsa=Funcsb
	    Funcba = (bottom0+bottom1*xa0+bottom2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,bottom3)*exp(1.0+bottom4*xa0); //Funcba=Funcbb
	    FuncABc = Funcsa*Funcsa/xa0/xa0;
	    FuncABs = Funcca*Funcca/xa0/xa0;
	    FuncABb = Funcba*Funcba/xa0/xa0;
	    double totSec0 = 2*mVal*(
                                   (FuncABu)*weightu
                                   +(FuncABd)*weightd
                                   +(FuncABc)*weightc
                                   +(FuncABs)*weights
                                   +(FuncABb)*weightb
                                   );
	    totSec = 1.e-5*totSec0;
	    }
    
    if (totSec <= 0.) totSec = 0.00001;
    
    totSec*=2.; // correcting for the half weight factor
	return totSec;
}
////////////////////////////////////////////////////////////////////////////
// GLU-GLU
////////////////////////////////////////////////////////////////////////////
double pdfs_gg(double m, double Y, double sqrtsVal)
{
    
    double s0 = sqrtsVal*sqrtsVal;
	double s = m*m;
	double Q = m;
	double xa = exp(Y)*sqrt(s/s0);
	double xb = exp(-Y)*sqrt(s/s0);
    
	double weightg = 1.0;
    
	//gluon params
	double g0par0 = 0.2282; double g0par1 = -0.0002252; double g0par2 = 1.383e-07;
	double g1par0 = 0.01968; double g1par1 = -0.0002993; double g1par2 = 1.986e-07;
	double g2par0 = 3.624; double g2par1 = -0.003164; double g2par2 = 1.941e-06;
	double g3par0 = -0.578; double g3par1 = -0.0003; double g3par2 = 1.828e-07;
	double g4par0 = -7.515; double g4par1 = -0.001355; double g4par2 = 8.199e-07;
    
	double gluon0 = g0par0 + g0par1*Q + g0par2*Q*Q;
	double gluon1 = g1par0 + g1par1*Q + g1par2*Q*Q;
	double gluon2 = g2par0 + g2par1*Q + g2par2*Q*Q;
	double gluon3 = g3par0 + g3par1*Q + g3par2*Q*Q;
	double gluon4 = g4par0 + g4par1*Q + g4par2*Q*Q;
    
	double Funcga = (gluon0+gluon1*xa+gluon2*pow(xa,2))*pow((1-xa),4)*pow(xa,gluon3)*exp(1.0+gluon4*xa);
	double Funcgb = (gluon0+gluon1*xb+gluon2*pow(xb,2))*pow((1-xb),4)*pow(xb,gluon3)*exp(1.0+gluon4*xb);
	double FuncABg = Funcga*Funcgb/xa/xb;
	
	double totSec = 2*m*((FuncABg)*weightg);
    
	if(( m <= 600. && fabs(Y) > 20*pow(m,-0.32)) || ( m > 600. && fabs(Y) > 21*pow(m,-0.34)))
        {
	    //Find totSec when mZZ, Y=0
	    double xa0 = sqrt(s/s0); //at Y=0 xa=xb
                                   //if xa=xb then Funcga=Funcgb
	    Funcga = (gluon0+gluon1*xa0+gluon2*pow(xa0,2))*pow((1-xa0),4)*pow(xa0,gluon3)*exp(1.0+gluon4*xa0);
	    FuncABg = Funcga*Funcga/xa0/xa0;
	    double totSec0 = 2*m*((FuncABg)*weightg);
	    totSec = 1.e-5*totSec0;
        }
    
    if (totSec <= 0.) totSec = 0.00001; 
    totSec/=10.;
    
    totSec*=2.; // correcting for the half weight factor
    return totSec;	                                        
}
////////////////////////////////////////////////////////////////////////////
// Q-QBARPRIME
////////////////////////////////////////////////////////////////////////////
double pdfs_qqbarPrime(double m, double Y, double sqrts)
{
    
    double YVal = Y;
    double sqrtsVal = sqrts;
    double mVal = m;
    
    double s0 = sqrtsVal*sqrtsVal;
	double s = mVal*mVal;
	double Q = mVal;
	double xa = exp(YVal)*sqrt(s/s0);
	double xb = exp(-YVal)*sqrt(s/s0);
    
	
    // missing CKM mixing for the moment
	double weightu = 0.5;
	double weightd = 0.5;
	double weightc = 1.0;
	double weights = 1.0;
	//double weightb = 0.;
	
    
	// PDF parameters
	// up params
	double u0par0 = 0.03134; double u0par1 =-2.068e-05; double u0par2 = 1.283e-08; 
	double u1par0 = 0.9; double u1par1 =-0.0004307; double u1par2 = 2.458e-07;
	double u2par0 =-0.1369; double u2par1 = 0.003423; double u2par2 =-2.155e-06;
	double u3par0 =-0.4013; double u3par1 =-0.0002574; double u3par2 = 1.561e-07;
	double u4par0 = 0.5782; double u4par1 =-0.004728; double u4par2 = 2.906e-06;
	double ubar0par0 = 0.02856; double ubar0par1 =-2.112e-05; double ubar0par2 = 1.272e-08;
	double ubar1par0 =-0.06822; double ubar1par1 = 3.172e-05; double ubar1par2 =-2.008e-08;
	double ubar2par0 = 0.1967; double ubar2par1 =-0.000118; double ubar2par2 = 6.871e-08;
	double ubar3par0 =-0.2251; double ubar3par1 = 0.0001295; double ubar3par2 =-7.181e-08;
	double ubar4par0 =-0.4068; double ubar4par1 =-0.0002956; double ubar4par2 = 1.783e-07;
	double ubar5par0 =-2.251; double ubar5par1 =-0.0001699; double ubar5par2 = 1.492e-07;
	// down params
	double d0par0 = 0.03278; double d0par1 =-2.915e-05; double d0par2 = 1.809e-08;
	double d1par0 = 0.479; double d1par1 =-0.0002559; double d1par2 = 1.557e-07;
	double d2par0 =-0.5972; double d2par1 = 0.0003118; double d2par2 =-1.905e-07;
	double d3par0 =-0.3892; double d3par1 =-0.000317; double d3par2 = 1.944e-07;
	double d4par0 = 0.5007; double d4par1 =-0.001665; double d4par2 = 9.895e-07;
	double dbar0par0 = 0.02328; double dbar0par1 =-1.367e-05; double dbar0par2 = 8.246e-09;
	double dbar1par0 = 0.09422; double dbar1par1 =-0.0001019; double dbar1par2 = 6.375e-08;
	double dbar2par0 =-0.5296; double dbar2par1 = 0.000466; double dbar2par2 =-2.896e-07;
	double dbar3par0 = 0.5354; double dbar3par1 =-0.0004404; double dbar3par2 = 2.728e-07;
	double dbar4par0 =-0.4386; double dbar4par1 =-0.0002605; double dbar4par2 = 1.582e-07;
	double dbar5par0 =-1.289; double dbar5par1 =-0.001618; double dbar5par2 = 9.601e-07;
	// charm, strange, bottom params
	double c0par0 = 0.01829; double c0par1 =-6.93e-06; double c0par2 = 3.796e-09;
	double c1par0 = 0.03081; double c1par1 = 4.325e-05; double c1par2 =-3.95e-08;
	double c2par0 = 0.5398; double c2par1 =-4.284e-05; double c2par2 =-1.362e-08;
	double c3par0 =-0.5986; double c3par1 = 0.002565; double c3par2 =-1.937e-06;
	double c4par0 =-0.4534; double c4par1 =-0.0002329; double c4par2 = 1.343e-07;
	double c5par0 =-8.657; double c5par1 =-0.005157; double c5par2 = 3.68e-06;
	double s0par0 = 0.01312; double s0par1 =-3.743e-06; double s0par2 = 2.076e-09;
	double s1par0 =-0.001416; double s1par1 =-7.649e-06; double s1par2 = 4.757e-09;
	double s2par0 = 0.2864; double s2par1 =-6.693e-05; double s2par2 = 3.566e-08;
	double s3par0 =-0.4857; double s3par1 =-0.000253; double s3par2 = 1.541e-07;
	double s4par0 =-10.33; double s4par1 =-0.001601; double s4par2 = 9.718e-07;
	double b0par0 = 0.005934; double b0par1 = 2.516e-06; double b0par2 =-1.828e-09;
	double b1par0 =-0.003063; double b1par1 =-6.761e-06; double b1par2 = 4.298e-09;
	double b2par0 = 0.1174; double b2par1 = 3.752e-05; double b2par2 =-2.863e-08;
	double b3par0 =-0.5549; double b3par1 =-0.0002205; double b3par2 = 1.334e-07;
	double b4par0 =-10.18; double b4par1 =-0.001136; double b4par2 = 6.931e-07;
    
	
	// PDF definition
	double up0 = u0par0 + u0par1*Q + u0par2*Q*Q;
	double up1 = u1par0 + u1par1*Q + u1par2*Q*Q;
	double up2 = u2par0 + u2par1*Q + u2par2*Q*Q;
	double up3 = u3par0 + u3par1*Q + u3par2*Q*Q;
	double up4 = u4par0 + u4par1*Q + u4par2*Q*Q;
	double antiup0 = ubar0par0 + ubar0par1*Q + ubar0par2*Q*Q;
	double antiup1 = ubar1par0 + ubar1par1*Q + ubar1par2*Q*Q;
	double antiup2 = ubar2par0 + ubar2par1*Q + ubar2par2*Q*Q;
	double antiup3 = ubar3par0 + ubar3par1*Q + ubar3par2*Q*Q;
	double antiup4 = ubar4par0 + ubar4par1*Q + ubar4par2*Q*Q;
	double antiup5 = ubar5par0 + ubar5par1*Q + ubar5par2*Q*Q;
	double down0 = d0par0 + d0par1*Q + d0par2*Q*Q;
	double down1 = d1par0 + d1par1*Q + d1par2*Q*Q;
	double down2 = d2par0 + d2par1*Q + d2par2*Q*Q;
	double down3 = d3par0 + d3par1*Q + d3par2*Q*Q;
	double down4 = d4par0 + d4par1*Q + d4par2*Q*Q;
	double antidown0 = dbar0par0 + dbar0par1*Q + dbar0par2*Q*Q;
	double antidown1 = dbar1par0 + dbar1par1*Q + dbar1par2*Q*Q;
	double antidown2 = dbar2par0 + dbar2par1*Q + dbar2par2*Q*Q;
	double antidown3 = dbar3par0 + dbar3par1*Q + dbar3par2*Q*Q;
	double antidown4 = dbar4par0 + dbar4par1*Q + dbar4par2*Q*Q;
	double antidown5 = dbar5par0 + dbar5par1*Q + dbar5par2*Q*Q;
	double charm0 = c0par0 + c0par1*Q + c0par2*Q*Q;
	double charm1 = c1par0 + c1par1*Q + c1par2*Q*Q;
	double charm2 = c2par0 + c2par1*Q + c2par2*Q*Q;
	double charm3 = c3par0 + c3par1*Q + c3par2*Q*Q;
	double charm4 = c4par0 + c4par1*Q + c4par2*Q*Q;
	double charm5 = c5par0 + c5par1*Q + c5par2*Q*Q;
	double strange0 = s0par0 + s0par1*Q + s0par2*Q*Q;
	double strange1 = s1par0 + s1par1*Q + s1par2*Q*Q;
	double strange2 = s2par0 + s2par1*Q + s2par2*Q*Q;
	double strange3 = s3par0 + s3par1*Q + s3par2*Q*Q;
	double strange4 = s4par0 + s4par1*Q + s4par2*Q*Q;
	double bottom0 = b0par0 + b0par1*Q + b0par2*Q*Q;
	double bottom1 = b1par0 + b1par1*Q + b1par2*Q*Q;
	double bottom2 = b2par0 + b2par1*Q + b2par2*Q*Q;
	double bottom3 = b3par0 + b3par1*Q + b3par2*Q*Q;
	double bottom4 = b4par0 + b4par1*Q + b4par2*Q*Q;
    
	double FuncAu1 = (up0+up1*xa+up2*pow(xa,2))*pow((1-xa),4)*pow(xa,up3)*exp(1.0+up4*xa);
	double FuncBu1 = (antiup0+antiup1*xb+antiup2*pow(xb,2)+antiup3*pow(xb,3))*pow((1-xb),4)*pow(xb,antiup4)*exp(1.0+antiup5*xb);
	double FuncAu2 = (up0+up1*xb+up2*pow(xb,2))*pow((1-xb),4)*pow(xb,up3)*exp(1.0+up4*xb);
	double FuncBu2 = (antiup0+antiup1*xa+antiup2*pow(xa,2)+antiup3*pow(xa,3))*pow((1-xa),4)*pow(xa,antiup4)*exp(1.0+antiup5*xa);
	double FuncABu = FuncAu1/xa*FuncBu1/xb+FuncAu2/xa*FuncBu2/xb; // uubar + sym
    
	
	double FuncAd1 = (down0+down1*xa+down2*pow(xa,2))*pow((1-xa),4)*pow(xa,down3)*exp(1.0+down4*xa);
	double FuncBd1 = (antidown0+antidown1*xb+antidown2*pow(xb,2)+antidown3*pow(xb,3))*pow((1-xb),4)*pow(xb,antidown4)*exp(1.0+antidown5*xb);
	double FuncAd2 = (down0+down1*xb+down2*pow(xb,2))*pow((1-xb),4)*pow(xb,down3)*exp(1.0+down4*xb);
	double FuncBd2 = (antidown0+antidown1*xa+antidown2*pow(xa,2)+antidown3*pow(xa,3))*pow((1-xa),4)*pow(xa,antidown4)*exp(1.0+antidown5*xa);
	double FuncABd = FuncAd1/xa*FuncBd1/xb+FuncAd2/xa*FuncBd2/xb;  // ddbar + sym
    
	double Funcca = (charm0+charm1*xa+charm2*pow(xa,2)+charm3*pow(xa,3))*pow((1-xa),4)*pow(xa,charm4)*exp(1.0+charm5*xa);
	double Funccb = (charm0+charm1*xb+charm2*pow(xb,2)+charm3*pow(xb,3))*pow((1-xb),4)*pow(xb,charm4)*exp(1.0+charm5*xb);
	double Funcsa = (strange0+strange1*xa+strange2*pow(xa,2))*pow((1-xa),4)*pow(xa,strange3)*exp(1.0+strange4*xa);
	double Funcsb = (strange0+strange1*xb+strange2*pow(xb,2))*pow((1-xb),4)*pow(xb,strange3)*exp(1.0+strange4*xb);
	double Funcba = (bottom0+bottom1*xa+bottom2*pow(xa,2))*pow((1-xa),4)*pow(xa,bottom3)*exp(1.0+bottom4*xa);
	double Funcbb = (bottom0+bottom1*xb+bottom2*pow(xb,2))*pow((1-xb),4)*pow(xb,bottom3)*exp(1.0+bottom4*xb);
    //double FuncABc = Funcsa*Funcsb/xa/xb; // ccbar + sym
	//double FuncABs = Funcca*Funccb/xa/xb; // ssbar + sym
	//double FuncABb = Funcba*Funcbb/xa/xb; // bbbar + sym
    
    // For W+'s
    // u-dbar 
    double FuncAB_udbar = FuncAu1/xa*FuncBd1/xb + FuncAu2/xa*FuncBd2/xb;
    // c-sbar
    double FuncAB_csbar = Funcca*Funcsb/xa/xb;
    
    // For W+'s
    // d-ubar 
    double FuncAB_dubar = FuncAd1/xa*FuncBu1/xb + FuncAd2/xa*FuncBu2/xb;
    // s-cbar
    double FuncAB_scbar = Funccb*Funcsa/xa/xb;    
    
	double totSec = 2*mVal*(
                              (FuncAB_udbar)*weightu
                              +(FuncAB_dubar)*weightd
                              +(FuncAB_csbar)*weightc
                              +(FuncAB_scbar)*weights
                              );
    
    if (totSec <= 0.) totSec = 0.00001;
    
    totSec*=2.; // correcting for the half weight factor
	return totSec;
}

