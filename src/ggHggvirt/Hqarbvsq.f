      double precision function Hqarbvsq(i1,i2,i3,i4)
      implicit none
C-----Author: Keith Ellis
C-----July  2005
C----Matrix element squared for the virtual corrections to the process
C----q(p1)+r(p3) -> q(p2)+r(p4)+Higgs
C----summed over incoming/outgoing colors and spins.
C    with a factor of gsq*ason2pi*Asq removed.
      include 'constants.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'scheme.f'
      double complex L0,L1,Lsm1,Lsm1_2me,lnrat
      double precision s123,s124,s134,s234,s12,s13,s14,s23,s24,s34
      double precision fun1,fun2,fun3,fun4,fun5,fun6,fun7,
     . fun8,fun9,fun10,fun11,fun12
      double precision Hqarb,Hqarbsq,mhsq
      integer i1,i2,i3,i4
      logical CheckEGZ
      common/CheckEGZ/CheckEGZ
!$omp threadprivate(/CheckEGZ/)
 
c      double complex function L0(x,y)
c      double complex function L1(x,y)
c      double complex function Ls0(x1,y1,x2,y2)
c      double complex function Ls1(x1,y1,x2,y2)
c      double complex function Lsm1(x1,y1,x2,y2)
c      double complex function Lsm1_2mh(s,t,m1sq,m2sq)
c      double complex function Lsm1_2me(s,t,m1sq,m3sq)

      fun1(i1,i2,i3,i4)=
     . -s(i1,i2)*s(i3,i4)/(2d0*s(i1,i3)**2)
     . -s(i2,i4)**2/(2d0*s(i1,i2)*s(i3,i4))
     . +(3d0*s(i1,i3)*s(i2,i4)-s(i2,i3)**2+s(i1,i4)*s(i2,i3)-s(i1,i4)**2
     . -s(i1,i3)**2)/(s(i1,i2)*s(i3,i4))
     . -s(i1,i4)**2*s(i2,i3)**2/(2d0*s(i1,i2)*s(i1,i3)**2*s(i3,i4))
     .   -2d0*(s(i1,i3)*s(i2,i4)-s(i1,i4)*s(i2,i3))**2
     .  /(s(i1,i2)**2*s(i3,i4)**2)+s(i2,i4)/s(i1,i3)
     .  +s(i1,i4)*s(i2,i3)/s(i1,i3)**2-2d0
 
      fun2(i1,i2,i3,i4) = 
     . (s(i1,i2)*s(i3,i4)*(s(i1,i2)*s(i3,i4)+s(i2,i3)
     . *(s(i2,i4)+2d0*s(i2,i3)-s(i1,i4)))+s(i2,i3)**2
     . *(s(i2,i4)+s(i1,i4))**2)
     .  /(2*s(i1,i2)**3*s(i3,i4))


      fun3(i1,i2,i3,i4) = 
     . s(i3,i4)/(2d0*s(i2,i3))+s(i2,i3)*(s(i2,i4)+s(i1,i4))
     . *(s(i2,i4)+4d0*s(i2,i3)+3d0*s(i1,i4))/(2d0*s(i1,i2)**2*s(i3,i4))
     . +(3d0*s(i2,i4)+4d0*s(i2,i3))/(2d0*s(i1,i2))

      fun4(i1,i2,i3,i4) =  
     . -2d0*s(i3,i4)/s(i2,i3)-s(i2,i3)*(s(i2,i4)+s(i1,i4))
     . *(s(i2,i4)+2d0*s(i2,i3)+5d0*s(i1,i4))/(2d0*s(i1,i2)**2*s(i3,i4))
     . -(4d0*s(i2,i4)+6*s(i2,i3)-3d0*s(i1,i4))/(2d0*s(i1,i2))


      fun5(i1,i2,i3,i4) =   
     . -(s(i2,i4)**2*(s(i3,i4)**2-s(i2,i4)*s(i3,i4)
     . +s(i1,i3)*s(i3,i4)-s(i1,i3)*s(i2,i4)+s(i1,i4)*s(i2,i3))
     . +s(i1,i4)*s(i2,i3)*(2d0*s(i2,i4)*s(i3,i4)
     . -s(i1,i3)*s(i2,i4)+s(i1,i4)*s(i2,i3)))
     . /(s(i1,i4)*s(i2,i4)**2*s(i3,i4))


      fun6(i1,i2,i3,i4)=
     . s(i1,i2)*s(i3,i4)/(2d0*s(i1,i3)*s(i2,i3))
     . +(4d0*s(i2,i3)*s(i2,i4)+2d0*s(i1,i4)*s(i2,i4)
     . -3d0*s(i1,i3)*s(i2,i4)+3d0*s(i1,i4)*s(i2,i3))
     . /(2d0*s(i1,i2)*s(i3,i4))+s(i1,i4)**2*s(i2,i3)
     . /(2d0*s(i1,i2)*s(i1,i3)*s(i3,i4))-s(i1,i4)/s(i1,i3)+0.5d0
 

      fun7(i1,i2,i3,i4)=
     . (s(i1,i4)**2*s(i2,i3)**2+s(i1,i3)*s(i3,i4)**2*s(i1,i2)
     . -s(i1,i3)*s(i1,i4)*s(i2,i3)*s(i2,i4))
     . /2/s(i1,i3)**2/s(i3,i4)/s(i1,i2)

      fun8(i1,i2,i3,i4)=
     . (s(i1,i4)*s(i2,i3)-s(i1,i3)*s(i2,i4))/s(i1,i2)/s(i3,i4)/2d0
 

      fun9(i1,i2,i3,i4)= 
     . +s(i1,i3)*s(i2,i4)**2/(s(i1,i2)*s(i2,i3)*s(i3,i4))
     . -s(i1,i2)*s(i3,i4)/(s(i1,i3)*s(i2,i3))-(2d0*s(i2,i4)**2
     . +2d0*s(i2,i3)*s(i2,i4)+5d0*s(i1,i4)*s(i2,i4)
     . -5d0*s(i1,i3)*s(i2,i4)+5d0*s(i1,i4)*s(i2,i3))
     . /(2d0*s(i1,i2)*s(i3,i4))
     . -s(i1,i4)**2*s(i2,i3)/(s(i1,i2)*s(i1,i3)*s(i3,i4))
     . -2d0*s(i2,i4)/s(i2,i3)+2d0*s(i1,i4)/s(i1,i3)-1d0


      fun10(i1,i2,i3,i4)= 
     . (s(i1,i2)*s(i2,i3)*s(i3,i4)**2+s(i1,i3)**2*s(i2,i4)**2
     . -s(i1,i3)*s(i1,i4)*s(i2,i3)*s(i2,i4))
     . /(s(i1,i2)*s(i2,i3)**2*s(i3,i4))

      fun11(i1,i2,i3,i4)=2d0*s(i1,i3)*s(i2,i4)/s(i1,i2)/s(i3,i4) 

      fun12(i1,i2,i3,i4)= 
     . +s(i1,i3)*(s(i1,i3)*(s(i1,i4)-s(i2,i4))
     . +2d0*s(i1,i4)*s(i2,i3))/2/s(i1,i2)**2
     . /s(i3,i4)+s(i1,i4)/s(i1,i2)/2d0

C-----End of statement functions


      s123=s(i1,i2)+s(i1,i3)+s(i2,i3)
      s124=s(i1,i2)+s(i1,i4)+s(i2,i4)
      s234=s(i2,i3)+s(i2,i4)+s(i3,i4)
      s134=s(i1,i3)+s(i1,i4)+s(i3,i4)

      s12=s(i1,i2)      
      s13=s(i1,i3)      
      s14=s(i1,i4)      
      s23=s(i2,i3)      
      s24=s(i2,i4)      
      s34=s(i3,i4)
      mhsq=s12+s13+s14+s23+s24+s34

      Hqarb=Hqarbsq(i1,i2,i3,i4)


      Hqarbvsq=epinv**2*(2.D0/xn-2.D0*xn)

      Hqarbvsq=Hqarbvsq+epinv * (
     &     + 3.D0/xn-3.D0*xn
     &     - dble(lnrat(-s34,musq))/xn
     &     - dble(lnrat(-s12,musq))/xn
     &     + 2.D0*dble(lnrat(-s13,musq))/xn
     &     + 2.D0*dble(lnrat(-s24,musq))/xn
     &     - 2.D0*dble(lnrat(-s14,musq))/xn
     &     + dble(lnrat(-s14,musq))*xn
     &     - 2.D0*dble(lnrat(-s23,musq))/xn
     &     + dble(lnrat(-s23,musq))*xn)

      Hqarbvsq = Hqarbvsq + 11.D0
     &     + 8.D0/xn
     &     + 80.D0/9.D0*xn
     &     - 20.D0/9.D0*nf
     &     + dble(Lsm1_2me(s134,s234,s34,mhsq))/xn
     &     - 2.D0*dble(Lsm1_2me(s123,s134,s13,mhsq))/xn
     &     + 2.D0*dble(Lsm1_2me(s123,s234,s23,mhsq))/xn
     &     - dble(Lsm1_2me(s123,s234,s23,mhsq))*xn
     &     + dble(Lsm1_2me(s123,s124,s12,mhsq))/xn
     &     + 2.D0*dble(Lsm1_2me(s124,s134,s14,mhsq))/xn
     &     - dble(Lsm1_2me(s124,s134,s14,mhsq))*xn
     &     - 2.D0*dble(Lsm1_2me(s124,s234,s24,mhsq))/xn
     &     - 3.D0/2.D0*dble(lnrat(-s34,musq))/xn
     &     - 13.D0/6.D0*dble(lnrat(-s34,musq))*xn
     &     + 2.D0/3.D0*dble(lnrat( - s34,musq))*nf
     &

      Hqarbvsq = Hqarbvsq 
     &     - 3.D0/2.D0*dble(lnrat(-s12,musq))/xn
     &     - 13.D0/6.D0*dble(lnrat(-s12,musq))*xn
     &     +  2.D0/3.D0*dble(lnrat(-s12,musq))*nf
     &     - dble(lnrat(-s13,musq)**2)/xn
     &     - dble(lnrat(-s24,musq)**2)/xn
     &     + dble(lnrat(-s14,musq)**2)/xn
     &     + dble(lnrat(-s23,musq)**2)/xn
     &     + 1.D0/2.D0*dble(lnrat(-s34,musq)**2)/xn
     &     + 1.D0/2.D0*dble(lnrat(-s12,musq)**2)/xn
     &     - 1.D0/2.D0*dble(lnrat(-s14,musq)**2)*xn
     &     - 1.D0/2.D0*dble(lnrat(-s23,musq)**2)*xn
     &

      if (checkEGZ) then
c--- remove finite renormalization
      Hqarbvsq = Hqarbvsq - 11.D0     
      endif

      Hqarbvsq=Hqarb*Hqarbvsq

C---end of terms proportional to lowest order      
      
      Hqarbvsq=Hqarbvsq+V/xn*(
     . +Dble(Lsm1(-s12,-s123,-s13,-s123))*fun1(i2,i1,i3,i4)
     . +Dble(Lsm1(-s12,-s124,-s24,-s124))*fun1(i1,i2,i4,i3)
     . +Dble(Lsm1(-s34,-s234,-s24,-s234))*fun1(i2,i1,i3,i4)
     . +Dble(Lsm1(-s34,-s134,-s13,-s134))*fun1(i1,i2,i4,i3))
     .  +V*(xn-2d0/xn)/2d0
     . *(Dble(Lsm1(-s12,-s123,-s23,-s123))*fun1(i1,i2,i3,i4)
     .  +Dble(Lsm1(-s12,-s124,-s14,-s124))*fun1(i2,i1,i4,i3)
     .  +Dble(Lsm1(-s34,-s234,-s23,-s234))*fun1(i2,i1,i4,i3)
     .  +Dble(Lsm1(-s34,-s134,-s14,-s134))*fun1(i1,i2,i3,i4))

      Hqarbvsq=Hqarbvsq
     . +fun2(i1,i2,i3,i4)*dble(L1(-s123,-s12))*(1d0/xn+xn)*0.5d0*V
     . +fun2(i2,i1,i4,i3)*dble(L1(-s124,-s12))*(1d0/xn+xn)*0.5d0*V  
     . +fun2(i4,i3,i2,i1)*dble(L1(-s234,-s34))*(1d0/xn+xn)*0.5d0*V     
     . +fun2(i3,i4,i1,i2)*dble(L1(-s134,-s34))*(1d0/xn+xn)*0.5d0*V

     . +fun3(i1,i2,i3,i4)*dble(L0(-s123,-s12))*xn*0.5d0*V
     . +fun3(i2,i1,i4,i3)*dble(L0(-s124,-s12))*xn*0.5d0*V  
     . +fun3(i4,i3,i2,i1)*dble(L0(-s234,-s34))*xn*0.5d0*V     
     . +fun3(i3,i4,i1,i2)*dble(L0(-s134,-s34))*xn*0.5d0*V

     . +fun4(i1,i2,i3,i4)*dble(L0(-s123,-s12))/xn*0.5d0*V
     . +fun4(i2,i1,i4,i3)*dble(L0(-s124,-s12))/xn*0.5d0*V  
     . +fun4(i4,i3,i2,i1)*dble(L0(-s234,-s34))/xn*0.5d0*V     
     . +fun4(i3,i4,i1,i2)*dble(L0(-s134,-s34))/xn*0.5d0*V

     . +fun5(i1,i2,i3,i4)*dble(L0(-s124,-s14))*0.5d0*V*(-xn/2+1d0/xn)
     . +fun5(i2,i1,i4,i3)*dble(L0(-s123,-s23))*0.5d0*V*(-xn/2+1d0/xn)
     . +fun5(i4,i3,i2,i1)*dble(L0(-s134,-s14))*0.5d0*V*(-xn/2+1d0/xn)
     . +fun5(i3,i4,i1,i2)*dble(L0(-s234,-s23))*0.5d0*V*(-xn/2+1d0/xn)

     . -fun5(i1,i2,i4,i3)*dble(L0(-s123,-s13))/xn*0.5d0*V
     . -fun5(i2,i1,i3,i4)*dble(L0(-s124,-s24))/xn*0.5d0*V  
     . -fun5(i4,i3,i1,i2)*dble(L0(-s234,-s24))/xn*0.5d0*V     
     . -fun5(i3,i4,i2,i1)*dble(L0(-s134,-s13))/xn*0.5d0*V

      Hqarbvsq=Hqarbvsq
     . +0.5d0*V*xn*dble(lnrat(-s123,-s12))*fun6(i1,i2,i3,i4)
     . +0.5d0*V*xn*dble(lnrat(-s124,-s12))*fun6(i2,i1,i4,i3)
     . +0.5d0*V*xn*dble(lnrat(-s234,-s34))*fun6(i4,i3,i2,i1)
     . +0.5d0*V*xn*dble(lnrat(-s134,-s34))*fun6(i3,i4,i1,i2)

     . +0.5d0*V*xn*dble(lnrat(-s123,-s23))*fun7(i1,i2,i3,i4)
     . +0.5d0*V*xn*dble(lnrat(-s124,-s14))*fun7(i2,i1,i4,i3)
     . +0.5d0*V*xn*dble(lnrat(-s234,-s23))*fun7(i4,i3,i2,i1)
     . +0.5d0*V*xn*dble(lnrat(-s134,-s14))*fun7(i3,i4,i1,i2)

     . +0.5d0*V*xn*dble(lnrat(-s12,-s14))*fun8(i1,i2,i3,i4)
     . +0.5d0*V*xn*dble(lnrat(-s12,-s14))*fun8(i2,i1,i4,i3)
     . +0.5d0*V*xn*dble(lnrat(-s34,-s23))*fun8(i4,i3,i2,i1)
     . +0.5d0*V*xn*dble(lnrat(-s34,-s23))*fun8(i3,i4,i1,i2)


      Hqarbvsq=Hqarbvsq
     . +0.5d0*V/xn*dble(lnrat(-s123,-s12))*fun9(i1,i2,i3,i4)
     . +0.5d0*V/xn*dble(lnrat(-s124,-s12))*fun9(i2,i1,i4,i3)
     . +0.5d0*V/xn*dble(lnrat(-s234,-s34))*fun9(i4,i3,i2,i1)
     . +0.5d0*V/xn*dble(lnrat(-s134,-s34))*fun9(i3,i4,i1,i2)

     . +0.5d0*V/xn*dble(lnrat(-s123,-s13))*fun10(i1,i2,i3,i4)
     . +0.5d0*V/xn*dble(lnrat(-s124,-s24))*fun10(i2,i1,i4,i3)
     . +0.5d0*V/xn*dble(lnrat(-s234,-s24))*fun10(i4,i3,i2,i1)
     . +0.5d0*V/xn*dble(lnrat(-s134,-s13))*fun10(i3,i4,i1,i2)

     . +0.5d0*V/xn*dble(lnrat(-s12,-s13))*fun11(i1,i2,i3,i4)
     . +0.5d0*V/xn*dble(lnrat(-s12,-s23))*fun11(i2,i1,i4,i3)
     . +0.5d0*V/xn*dble(lnrat(-s34,-s24))*fun11(i4,i3,i2,i1)
     . +0.5d0*V/xn*dble(lnrat(-s34,-s14))*fun11(i3,i4,i1,i2)


      Hqarbvsq=Hqarbvsq

     . -0.5d0*V/xn*dble(lnrat(-s123,-s23))*fun10(i2,i1,i3,i4)
     . -0.5d0*V/xn*dble(lnrat(-s124,-s14))*fun10(i1,i2,i4,i3)
     . -0.5d0*V/xn*dble(lnrat(-s234,-s23))*fun10(i3,i4,i2,i1)
     . -0.5d0*V/xn*dble(lnrat(-s134,-s14))*fun10(i4,i3,i1,i2)

     . -0.5d0*V/xn*dble(lnrat(-s12,-s13))*fun11(i2,i1,i3,i4)
     . -0.5d0*V/xn*dble(lnrat(-s12,-s14))*fun11(i1,i2,i4,i3)
     . -0.5d0*V/xn*dble(lnrat(-s34,-s23))*fun11(i3,i4,i2,i1)
     . -0.5d0*V/xn*dble(lnrat(-s34,-s24))*fun11(i4,i3,i1,i2)


     . +fun12(i1,i2,i3,i4)*0.5d0*V*(xn+1d0/xn)
     . +fun12(i2,i1,i4,i3)*0.5d0*V*(xn+1d0/xn)
     . +fun12(i4,i3,i2,i1)*0.5d0*V*(xn+1d0/xn)
     . +fun12(i3,i4,i1,i2)*0.5d0*V*(xn+1d0/xn)

c--- translation between schemes, according to Eq. (54) of EGZ
      if (scheme .eq. 'dred') then
      Hqarbvsq=Hqarbvsq+(2d0*Cf)*Hqarb
      endif

      return
      end
