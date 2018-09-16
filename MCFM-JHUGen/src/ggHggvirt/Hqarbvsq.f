      function Hqarbvsq(i1,i2,i3,i4)
      implicit none
      include 'types.f'
      real(dp):: Hqarbvsq
      
C-----Author: Keith Ellis
C-----July  2005
C----Matrix element squared for the virtual corrections to the process
C----q(p1)+r(p3) -> q(p2)+r(p4)+Higgs
C----summed over incoming/outgoing colors and spins.
C    with a factor of gsq*ason2pi*Asq removed.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'scheme.f'
      complex(dp):: L0,L1,Lsm1,Lsm1_2me,lnrat
      real(dp):: s123,s124,s134,s234,s12,s13,s14,s23,s24,s34
      real(dp):: fun1,fun2,fun3,fun4,fun5,fun6,fun7,
     & fun8,fun9,fun10,fun11,fun12
      real(dp):: Hqarb,Hqarbsq,mhsq
      integer:: i1,i2,i3,i4
      logical:: CheckEGZ
      common/CheckEGZ/CheckEGZ
!$omp threadprivate(/CheckEGZ/)
 
c      function L0(x,y)
c      implicit none
c      include 'types.f'
c      complex(dp):: L0
c      function L1(x,y)
c      implicit none
c      include 'types.f'
c      complex(dp):: L1
c      function Ls0(x1,y1,x2,y2)
c      implicit none
c      include 'types.f'
c      complex(dp):: Ls0
c      function Ls1(x1,y1,x2,y2)
c      implicit none
c      include 'types.f'
c      complex(dp):: Ls1
c      function Lsm1(x1,y1,x2,y2)
c      implicit none
c      include 'types.f'
c      complex(dp):: Lsm1
c      function Lsm1_2mh(s,t,m1sq,m2sq)
c      implicit none
c      include 'types.f'
c      complex(dp):: Lsm1_2mh
c      function Lsm1_2me(s,t,m1sq,m3sq)
c      implicit none
c      include 'types.f'
c      complex(dp):: Lsm1_2me

      fun1(i1,i2,i3,i4)=
     & -s(i1,i2)*s(i3,i4)/(2._dp*s(i1,i3)**2)
     & -s(i2,i4)**2/(2._dp*s(i1,i2)*s(i3,i4))
     & +(3._dp*s(i1,i3)*s(i2,i4)-s(i2,i3)**2+s(i1,i4)*s(i2,i3)-s(i1,i4)**2
     & -s(i1,i3)**2)/(s(i1,i2)*s(i3,i4))
     & -s(i1,i4)**2*s(i2,i3)**2/(2._dp*s(i1,i2)*s(i1,i3)**2*s(i3,i4))
     &   -2._dp*(s(i1,i3)*s(i2,i4)-s(i1,i4)*s(i2,i3))**2
     &  /(s(i1,i2)**2*s(i3,i4)**2)+s(i2,i4)/s(i1,i3)
     &  +s(i1,i4)*s(i2,i3)/s(i1,i3)**2-2._dp
 
      fun2(i1,i2,i3,i4) = 
     & (s(i1,i2)*s(i3,i4)*(s(i1,i2)*s(i3,i4)+s(i2,i3)
     & *(s(i2,i4)+2._dp*s(i2,i3)-s(i1,i4)))+s(i2,i3)**2
     & *(s(i2,i4)+s(i1,i4))**2)
     &  /(2*s(i1,i2)**3*s(i3,i4))


      fun3(i1,i2,i3,i4) = 
     & s(i3,i4)/(2._dp*s(i2,i3))+s(i2,i3)*(s(i2,i4)+s(i1,i4))
     & *(s(i2,i4)+4._dp*s(i2,i3)+3._dp*s(i1,i4))/(2._dp*s(i1,i2)**2*s(i3,i4))
     & +(3._dp*s(i2,i4)+4._dp*s(i2,i3))/(2._dp*s(i1,i2))

      fun4(i1,i2,i3,i4) =  
     & -2._dp*s(i3,i4)/s(i2,i3)-s(i2,i3)*(s(i2,i4)+s(i1,i4))
     & *(s(i2,i4)+2._dp*s(i2,i3)+5._dp*s(i1,i4))/(2._dp*s(i1,i2)**2*s(i3,i4))
     & -(4._dp*s(i2,i4)+6*s(i2,i3)-3._dp*s(i1,i4))/(2._dp*s(i1,i2))


      fun5(i1,i2,i3,i4) =   
     & -(s(i2,i4)**2*(s(i3,i4)**2-s(i2,i4)*s(i3,i4)
     & +s(i1,i3)*s(i3,i4)-s(i1,i3)*s(i2,i4)+s(i1,i4)*s(i2,i3))
     & +s(i1,i4)*s(i2,i3)*(2._dp*s(i2,i4)*s(i3,i4)
     & -s(i1,i3)*s(i2,i4)+s(i1,i4)*s(i2,i3)))
     & /(s(i1,i4)*s(i2,i4)**2*s(i3,i4))


      fun6(i1,i2,i3,i4)=
     & s(i1,i2)*s(i3,i4)/(2._dp*s(i1,i3)*s(i2,i3))
     & +(4._dp*s(i2,i3)*s(i2,i4)+2._dp*s(i1,i4)*s(i2,i4)
     & -3._dp*s(i1,i3)*s(i2,i4)+3._dp*s(i1,i4)*s(i2,i3))
     & /(2._dp*s(i1,i2)*s(i3,i4))+s(i1,i4)**2*s(i2,i3)
     & /(2._dp*s(i1,i2)*s(i1,i3)*s(i3,i4))-s(i1,i4)/s(i1,i3)+0.5_dp
 

      fun7(i1,i2,i3,i4)=
     & (s(i1,i4)**2*s(i2,i3)**2+s(i1,i3)*s(i3,i4)**2*s(i1,i2)
     & -s(i1,i3)*s(i1,i4)*s(i2,i3)*s(i2,i4))
     & /2/s(i1,i3)**2/s(i3,i4)/s(i1,i2)

      fun8(i1,i2,i3,i4)=
     & (s(i1,i4)*s(i2,i3)-s(i1,i3)*s(i2,i4))/s(i1,i2)/s(i3,i4)/2._dp
 

      fun9(i1,i2,i3,i4)= 
     & +s(i1,i3)*s(i2,i4)**2/(s(i1,i2)*s(i2,i3)*s(i3,i4))
     & -s(i1,i2)*s(i3,i4)/(s(i1,i3)*s(i2,i3))-(2._dp*s(i2,i4)**2
     & +2._dp*s(i2,i3)*s(i2,i4)+5._dp*s(i1,i4)*s(i2,i4)
     & -5._dp*s(i1,i3)*s(i2,i4)+5._dp*s(i1,i4)*s(i2,i3))
     & /(2._dp*s(i1,i2)*s(i3,i4))
     & -s(i1,i4)**2*s(i2,i3)/(s(i1,i2)*s(i1,i3)*s(i3,i4))
     & -2._dp*s(i2,i4)/s(i2,i3)+2._dp*s(i1,i4)/s(i1,i3)-1._dp


      fun10(i1,i2,i3,i4)= 
     & (s(i1,i2)*s(i2,i3)*s(i3,i4)**2+s(i1,i3)**2*s(i2,i4)**2
     & -s(i1,i3)*s(i1,i4)*s(i2,i3)*s(i2,i4))
     & /(s(i1,i2)*s(i2,i3)**2*s(i3,i4))

      fun11(i1,i2,i3,i4)=2._dp*s(i1,i3)*s(i2,i4)/s(i1,i2)/s(i3,i4) 

      fun12(i1,i2,i3,i4)= 
     & +s(i1,i3)*(s(i1,i3)*(s(i1,i4)-s(i2,i4))
     & +2._dp*s(i1,i4)*s(i2,i3))/2/s(i1,i2)**2
     & /s(i3,i4)+s(i1,i4)/s(i1,i2)/2._dp

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


      Hqarbvsq=epinv**2*(2._dp/xn-2._dp*xn)

      Hqarbvsq=Hqarbvsq+epinv * (
     &     + 3._dp/xn-3._dp*xn
     &     - real(lnrat(-s34,musq))/xn
     &     - real(lnrat(-s12,musq))/xn
     &     + 2._dp*real(lnrat(-s13,musq))/xn
     &     + 2._dp*real(lnrat(-s24,musq))/xn
     &     - 2._dp*real(lnrat(-s14,musq))/xn
     &     + real(lnrat(-s14,musq))*xn
     &     - 2._dp*real(lnrat(-s23,musq))/xn
     &     + real(lnrat(-s23,musq))*xn)

      Hqarbvsq = Hqarbvsq + 11._dp
     &     + 8._dp/xn
     &     + 80._dp/9._dp*xn
     &     - 20._dp/9._dp*nf
     &     + real(Lsm1_2me(s134,s234,s34,mhsq))/xn
     &     - 2._dp*real(Lsm1_2me(s123,s134,s13,mhsq))/xn
     &     + 2._dp*real(Lsm1_2me(s123,s234,s23,mhsq))/xn
     &     - real(Lsm1_2me(s123,s234,s23,mhsq))*xn
     &     + real(Lsm1_2me(s123,s124,s12,mhsq))/xn
     &     + 2._dp*real(Lsm1_2me(s124,s134,s14,mhsq))/xn
     &     - real(Lsm1_2me(s124,s134,s14,mhsq))*xn
     &     - 2._dp*real(Lsm1_2me(s124,s234,s24,mhsq))/xn
     &     - 3._dp/2._dp*real(lnrat(-s34,musq))/xn
     &     - 13._dp/6._dp*real(lnrat(-s34,musq))*xn
     &     + 2._dp/3._dp*real(lnrat( - s34,musq))*nf
     &

      Hqarbvsq = Hqarbvsq 
     &     - 3._dp/2._dp*real(lnrat(-s12,musq))/xn
     &     - 13._dp/6._dp*real(lnrat(-s12,musq))*xn
     &     +  2._dp/3._dp*real(lnrat(-s12,musq))*nf
     &     - real(lnrat(-s13,musq)**2)/xn
     &     - real(lnrat(-s24,musq)**2)/xn
     &     + real(lnrat(-s14,musq)**2)/xn
     &     + real(lnrat(-s23,musq)**2)/xn
     &     + 1._dp/2._dp*real(lnrat(-s34,musq)**2)/xn
     &     + 1._dp/2._dp*real(lnrat(-s12,musq)**2)/xn
     &     - 1._dp/2._dp*real(lnrat(-s14,musq)**2)*xn
     &     - 1._dp/2._dp*real(lnrat(-s23,musq)**2)*xn
     &

      if (checkEGZ) then
c--- remove finite renormalization
      Hqarbvsq = Hqarbvsq - 11._dp     
      endif

      Hqarbvsq=Hqarb*Hqarbvsq

C---end of terms proportional to lowest order      
      
      Hqarbvsq=Hqarbvsq+V/xn*(
     & +Dble(Lsm1(-s12,-s123,-s13,-s123))*fun1(i2,i1,i3,i4)
     & +Dble(Lsm1(-s12,-s124,-s24,-s124))*fun1(i1,i2,i4,i3)
     & +Dble(Lsm1(-s34,-s234,-s24,-s234))*fun1(i2,i1,i3,i4)
     & +Dble(Lsm1(-s34,-s134,-s13,-s134))*fun1(i1,i2,i4,i3))
     &  +V*(xn-2._dp/xn)/2._dp
     & *(Dble(Lsm1(-s12,-s123,-s23,-s123))*fun1(i1,i2,i3,i4)
     &  +Dble(Lsm1(-s12,-s124,-s14,-s124))*fun1(i2,i1,i4,i3)
     &  +Dble(Lsm1(-s34,-s234,-s23,-s234))*fun1(i2,i1,i4,i3)
     &  +Dble(Lsm1(-s34,-s134,-s14,-s134))*fun1(i1,i2,i3,i4))

      Hqarbvsq=Hqarbvsq
     & +fun2(i1,i2,i3,i4)*real(L1(-s123,-s12))*(1._dp/xn+xn)*0.5_dp*V
     & +fun2(i2,i1,i4,i3)*real(L1(-s124,-s12))*(1._dp/xn+xn)*0.5_dp*V  
     & +fun2(i4,i3,i2,i1)*real(L1(-s234,-s34))*(1._dp/xn+xn)*0.5_dp*V     
     & +fun2(i3,i4,i1,i2)*real(L1(-s134,-s34))*(1._dp/xn+xn)*0.5_dp*V

     & +fun3(i1,i2,i3,i4)*real(L0(-s123,-s12))*xn*0.5_dp*V
     & +fun3(i2,i1,i4,i3)*real(L0(-s124,-s12))*xn*0.5_dp*V  
     & +fun3(i4,i3,i2,i1)*real(L0(-s234,-s34))*xn*0.5_dp*V     
     & +fun3(i3,i4,i1,i2)*real(L0(-s134,-s34))*xn*0.5_dp*V

     & +fun4(i1,i2,i3,i4)*real(L0(-s123,-s12))/xn*0.5_dp*V
     & +fun4(i2,i1,i4,i3)*real(L0(-s124,-s12))/xn*0.5_dp*V  
     & +fun4(i4,i3,i2,i1)*real(L0(-s234,-s34))/xn*0.5_dp*V     
     & +fun4(i3,i4,i1,i2)*real(L0(-s134,-s34))/xn*0.5_dp*V

     & +fun5(i1,i2,i3,i4)*real(L0(-s124,-s14))*0.5_dp*V*(-xn/2+1._dp/xn)
     & +fun5(i2,i1,i4,i3)*real(L0(-s123,-s23))*0.5_dp*V*(-xn/2+1._dp/xn)
     & +fun5(i4,i3,i2,i1)*real(L0(-s134,-s14))*0.5_dp*V*(-xn/2+1._dp/xn)
     & +fun5(i3,i4,i1,i2)*real(L0(-s234,-s23))*0.5_dp*V*(-xn/2+1._dp/xn)

     & -fun5(i1,i2,i4,i3)*real(L0(-s123,-s13))/xn*0.5_dp*V
     & -fun5(i2,i1,i3,i4)*real(L0(-s124,-s24))/xn*0.5_dp*V  
     & -fun5(i4,i3,i1,i2)*real(L0(-s234,-s24))/xn*0.5_dp*V     
     & -fun5(i3,i4,i2,i1)*real(L0(-s134,-s13))/xn*0.5_dp*V

      Hqarbvsq=Hqarbvsq
     & +0.5_dp*V*xn*real(lnrat(-s123,-s12))*fun6(i1,i2,i3,i4)
     & +0.5_dp*V*xn*real(lnrat(-s124,-s12))*fun6(i2,i1,i4,i3)
     & +0.5_dp*V*xn*real(lnrat(-s234,-s34))*fun6(i4,i3,i2,i1)
     & +0.5_dp*V*xn*real(lnrat(-s134,-s34))*fun6(i3,i4,i1,i2)

     & +0.5_dp*V*xn*real(lnrat(-s123,-s23))*fun7(i1,i2,i3,i4)
     & +0.5_dp*V*xn*real(lnrat(-s124,-s14))*fun7(i2,i1,i4,i3)
     & +0.5_dp*V*xn*real(lnrat(-s234,-s23))*fun7(i4,i3,i2,i1)
     & +0.5_dp*V*xn*real(lnrat(-s134,-s14))*fun7(i3,i4,i1,i2)

     & +0.5_dp*V*xn*real(lnrat(-s12,-s14))*fun8(i1,i2,i3,i4)
     & +0.5_dp*V*xn*real(lnrat(-s12,-s14))*fun8(i2,i1,i4,i3)
     & +0.5_dp*V*xn*real(lnrat(-s34,-s23))*fun8(i4,i3,i2,i1)
     & +0.5_dp*V*xn*real(lnrat(-s34,-s23))*fun8(i3,i4,i1,i2)


      Hqarbvsq=Hqarbvsq
     & +0.5_dp*V/xn*real(lnrat(-s123,-s12))*fun9(i1,i2,i3,i4)
     & +0.5_dp*V/xn*real(lnrat(-s124,-s12))*fun9(i2,i1,i4,i3)
     & +0.5_dp*V/xn*real(lnrat(-s234,-s34))*fun9(i4,i3,i2,i1)
     & +0.5_dp*V/xn*real(lnrat(-s134,-s34))*fun9(i3,i4,i1,i2)

     & +0.5_dp*V/xn*real(lnrat(-s123,-s13))*fun10(i1,i2,i3,i4)
     & +0.5_dp*V/xn*real(lnrat(-s124,-s24))*fun10(i2,i1,i4,i3)
     & +0.5_dp*V/xn*real(lnrat(-s234,-s24))*fun10(i4,i3,i2,i1)
     & +0.5_dp*V/xn*real(lnrat(-s134,-s13))*fun10(i3,i4,i1,i2)

     & +0.5_dp*V/xn*real(lnrat(-s12,-s13))*fun11(i1,i2,i3,i4)
     & +0.5_dp*V/xn*real(lnrat(-s12,-s23))*fun11(i2,i1,i4,i3)
     & +0.5_dp*V/xn*real(lnrat(-s34,-s24))*fun11(i4,i3,i2,i1)
     & +0.5_dp*V/xn*real(lnrat(-s34,-s14))*fun11(i3,i4,i1,i2)


      Hqarbvsq=Hqarbvsq

     & -0.5_dp*V/xn*real(lnrat(-s123,-s23))*fun10(i2,i1,i3,i4)
     & -0.5_dp*V/xn*real(lnrat(-s124,-s14))*fun10(i1,i2,i4,i3)
     & -0.5_dp*V/xn*real(lnrat(-s234,-s23))*fun10(i3,i4,i2,i1)
     & -0.5_dp*V/xn*real(lnrat(-s134,-s14))*fun10(i4,i3,i1,i2)

     & -0.5_dp*V/xn*real(lnrat(-s12,-s13))*fun11(i2,i1,i3,i4)
     & -0.5_dp*V/xn*real(lnrat(-s12,-s14))*fun11(i1,i2,i4,i3)
     & -0.5_dp*V/xn*real(lnrat(-s34,-s23))*fun11(i3,i4,i2,i1)
     & -0.5_dp*V/xn*real(lnrat(-s34,-s24))*fun11(i4,i3,i1,i2)


     & +fun12(i1,i2,i3,i4)*0.5_dp*V*(xn+1._dp/xn)
     & +fun12(i2,i1,i4,i3)*0.5_dp*V*(xn+1._dp/xn)
     & +fun12(i4,i3,i2,i1)*0.5_dp*V*(xn+1._dp/xn)
     & +fun12(i3,i4,i1,i2)*0.5_dp*V*(xn+1._dp/xn)

c--- translation between schemes, according to Eq. (54) of EGZ
      if (scheme == 'dred') then
      Hqarbvsq=Hqarbvsq+(2._dp*Cf)*Hqarb
      endif

      return
      end
