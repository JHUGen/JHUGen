      function Hqaqavsq(i1,i2,i3,i4)
      implicit none
      include 'types.f'
      real(dp):: Hqaqavsq
      
C     Amplitude squared with certain factors removed
C     for the process H --> q(p1)+a(p2)+q(p3)+a(p4)
 
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
      real(dp):: fun1,fun2,fun3,fun4,fun5,fun6,fun7,fun8
      real(dp):: Hqaqa,Hqaqasq,mhsq
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


      fun1(s12,s13,s14,s23,s24,s34)=(
     &  - 2._dp*s12*s13*s24*s34+ s12*s13**2*s34
     &  - 4._dp*s12*s14*s23*s34+ s12*s24**2*s34
     &  + 2._dp*s12**2*s34**2- 2._dp*s13*s14*s23*s24
     &  - s13*s24**3+ s13**2*s14*s23
     &  - s13**3*s24+ s14*s23*s24**2+ 2._dp*s14**2*s23**2)
     & /(4._dp*s12*s23*s14*s34)

      fun2(s12,s13,s14,s23,s24,s34)=(
     &  - 2._dp*s12*s13*s14*s23*s24*s34
     &  + 2._dp*s12*s13**2*s14*s23*s34
     &  - 4._dp*s12*s13**2*s24**2*s34
     &  + 3._dp*s12*s13**3*s24*s34
     &  - 2._dp*s12*s13**4*s34
     &  + s12*s14**2*s23**2*s34
     &  + 3._dp*s12**2*s13*s24*s34**2
     &  - 3._dp*s12**2*s13**2*s34**2
     &  + s12**2*s14*s23*s34**2
     &  - s12**3*s34**3
     &  + 3._dp*s13*s14**2*s23**2*s24
     &  - 4._dp*s13**2*s14*s23*s24**2
     &  - 3._dp*s13**2*s14**2*s23**2
     &  + 3._dp*s13**3*s14*s23*s24
     &  + 2._dp*s13**3*s24**3
     &  - 2._dp*s13**4*s14*s23
     &  + 2._dp*s13**5*s24
     &  - s14**3*s23**3)
     & /(8._dp*s12*s13**2*s23*s14*s34)

      fun3(s12,s13,s14,s23,s24,s34)=(
     &  - 2._dp*s12*s14*s23*s34
     &  + s12*s23*s24*s34
     &  + s12*s23**2*s34
     &  + s12**2*s34**2
     &  + s14*s23**2*s24
     &  + s14*s23**3
     &  + s14**2*s23**2
     &  + s23**3*s24)
     & /(8._dp*s12**2*s14*s23)

      fun4(s12,s13,s14,s23,s24,s34)=(
     &  - 2._dp*s12*s23*s34+ 5._dp*s12*s24*s34- 3._dp*s14*s23*s24
     &  + 6._dp*s14*s23**2+ 4._dp*s23*s24**2- s23**2*s24)
     & /(8._dp*s14*s12*s23)

      fun5(s12,s13,s14,s23,s24,s34)=(
     &  + s12*s14*s23**2*s24*s34- 6._dp*s12*s14*s23**3*s34
     &  + 2._dp*s12*s14**2*s23**2*s34- 2._dp*s12*s23**2*s24**2*s34
     &  + 3._dp*s12*s23**3*s24*s34+ 2._dp*s12**2*s14*s23*s34**2
     &  - 3._dp*s12**2*s23*s24*s34**2+ 2._dp*s12**2*s23**2*s34**2
     &  - 2._dp*s12**3*s34**3- 2._dp*s14*s23**3*s24**2
     &  - 4._dp*s14**2*s23**3*s24- 2._dp*s14**3*s23**3)
     & /(8._dp*s12**2*s23**2*s14*s34)


      fun6(s12,s13,s14,s23,s24,s34)=(- 3._dp*s12*s34+ 3._dp*s13*s24
     &  + 5._dp*s14*s23- 2._dp*s23*s24+ 4._dp*s24**2)
     & /(8._dp*s14*s23)

      fun7(s12,s13,s14,s23,s24,s34)=(
     &  - 7._dp*s12*s13*s14*s23**2*s34
     &  - 2._dp*s12*s13*s23*s24**2*s34
     &  + 2._dp*s12*s13*s23**2*s24*s34
     &  - s12*s13**2*s23*s24*s34
     &  + 2._dp*s12*s14**2*s23**2*s34
     &  + s12**2*s13*s23*s34**2
     &  + 2._dp*s12**2*s14*s23*s34**2
     &  - 2._dp*s12**3*s34**3
     &  - 2._dp*s13*s14*s23**2*s24**2
     &  - 2._dp*s14**3*s23**3)
     & /(8._dp*s12*s13*s23**2*s14*s34)

      fun8(s12,s13,s14,s23,s24,s34)=(
     &  - 2._dp*s12*s13*s14*s24- 2._dp*s12*s13*s23*s24
     &  + s12*s13*s23**2+ s12*s13*s24**2
     &  + s12*s13**2*s24+ 2._dp*s12*s14*s23**2
     &  + 2._dp*s12*s14*s34**2+ 2._dp*s12*s14**2*s23
     &  + s12*s14**2*s24+ 2._dp*s12*s23*s34**2+s12**2*s13*s23
     &  + s12**2*s14*s24+ 2._dp*s12**2*s14*s34
     &  + 2._dp*s12**2*s23*s34- 2._dp*s13*s14*s24*s34+ s13*s14*s24**2
     &  + s13*s14*s34**2 + s13*s14**2*s34
     &  - 2._dp*s13*s23*s24*s34+ s13*s23*s24**2
     &  + s13*s24**2*s34+ s13**2*s14*s24
     &  + s13**2*s23*s24+ s13**2*s24*s34
     &  + 2._dp*s14*s23**2*s34+ 2._dp*s14**2*s23*s34
     &  + s23*s24*s34**2  + s23**2*s24*s34)
     &  /(64._dp*s12*s14*s23*s34)


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

      Hqaqa=Hqaqasq(i1,i2,i3,i4)


      Hqaqavsq = + epinv**2*(+2._dp/xn-2._dp*xn)
      Hqaqavsq = Hqaqavsq + epinv * (
     &     + 3._dp/xn- 3._dp*xn
     &     - 8._dp/3._dp*0.5_dp*nf
     &     + 4._dp/3._dp*nf
     &     - real(lnrat( - s14,musq))/xn
     &     - real(lnrat( - s23,musq))/xn
     &     - real(lnrat( - s12,musq))/xn
     &     - real(lnrat( - s34,musq))/xn
     &     + real(lnrat( - s24,musq))*(xn+1._dp/xn)
     &     + real(lnrat( - s13,musq))*(xn+1._dp/xn)
     &     )
      Hqaqavsq = Hqaqavsq + 11._dp
     &     + 8._dp/xn
     &     + 80._dp/9._dp*xn
     &     - 40._dp/9._dp*0.5_dp*nf
     &     + real(Lsm1_2me(s123,s234,s23,mhsq))/xn
     &     + real(Lsm1_2me(s123,s124,s12,mhsq))/xn
     &     + real(Lsm1_2me(s134,s234,s34,mhsq))/xn
     &     + real(Lsm1_2me(s124,s134,s14,mhsq))/xn
     &     - real(Lsm1_2me(s123,s134,s13,mhsq))*(xn+1._dp/xn)
     &     - real(Lsm1_2me(s124,s234,s24,mhsq))*(xn+1._dp/xn)
     &     - 3._dp/4._dp*real(lnrat( - s14,musq))/xn
     &     - 13._dp/12._dp*real(lnrat( - s14,musq))*xn
     &     + 2._dp/3._dp*real(lnrat( - s14,musq))*0.5_dp*nf
     &
      Hqaqavsq = Hqaqavsq 
     &     + 1._dp/2._dp*real(lnrat( - s14,musq)**2)/xn
     &     - 3._dp/4._dp*real(lnrat( - s23,musq))/xn
     &     - 13._dp/12._dp*real(lnrat( - s23,musq))*xn
     &     + 2._dp/3._dp*real(lnrat( - s23,musq))*0.5_dp*nf
     &     + 1._dp/2._dp*real(lnrat( - s23,musq)**2)/xn
     &     - 3._dp/4._dp*real(lnrat( - s12,musq))/xn
     &     - 13._dp/12._dp*real(lnrat( - s12,musq))*xn
     &     + 2._dp/3._dp*real(lnrat( - s12,musq))*0.5_dp*nf
     &     + 1._dp/2._dp*real(lnrat( - s12,musq)**2)/xn
     &     - 3._dp/4._dp*real(lnrat( - s34,musq))/xn
     &     - 13._dp/12._dp*real(lnrat( - s34,musq))*xn
     &     + 2._dp/3._dp*real(lnrat( - s34,musq))*0.5_dp*nf
     &     + 1._dp/2._dp*real(lnrat( - s34,musq)**2)/xn
     &     - 1._dp/2._dp*real(lnrat( - s24,musq)**2)*(xn+1._dp/xn)
     &     - 1._dp/2._dp*real(lnrat( - s13,musq)**2)*(xn+1._dp/xn)

      if (checkEGZ) then
c--- remove finite renormalization
      Hqaqavsq = Hqaqavsq - 11._dp     
      endif
    
      Hqaqavsq =Hqaqavsq*Hqaqa
      
      Hqaqavsq = Hqaqavsq+V*(1._dp+1._dp/xn**2)*( 
     & -real(Lsm1(-s12,-s124,-s24,-s124))*fun1(s12,s13,s14,s23,s24,s34)
     & -real(Lsm1(-s23,-s234,-s24,-s234))*fun1(s23,s13,s34,s12,s24,s14)
     & -real(Lsm1(-s14,-s124,-s24,-s124))*fun1(s14,s13,s12,s34,s24,s23)
     & -real(Lsm1(-s34,-s234,-s24,-s234))*fun1(s34,s13,s23,s14,s24,s12)
     & -real(Lsm1(-s34,-s134,-s13,-s134))*fun1(s12,s13,s14,s23,s24,s34)
     & -real(Lsm1(-s14,-s134,-s13,-s134))*fun1(s23,s13,s34,s12,s24,s14)
     & -real(Lsm1(-s12,-s123,-s13,-s123))*fun1(s14,s13,s12,s34,s24,s23)
     & -real(Lsm1(-s23,-s123,-s13,-s123))*fun1(s34,s13,s23,s14,s24,s12))

c--- NB: in sub-formula below, I have replaced the line:
c     & -real(Lsm1(-s12,-s124,-s24,-s124))*fun2(s12,s24,s23,s14,s13,s34)
c---     with the one:
c     & -real(Lsm1(-s12,-s124,-s14,-s124))*fun2(s12,s24,s23,s14,s13,s34)
c--- This is as prescribed by Eq. (37) of EGZ and restores 1<>3 symmetry

      Hqaqavsq = Hqaqavsq+V/xn**2*( 
     & -real(Lsm1(-s12,-s123,-s23,-s123))*fun2(s12,s13,s14,s23,s24,s34)
     & -real(Lsm1(-s23,-s123,-s12,-s123))*fun2(s23,s13,s34,s12,s24,s14)
c     & -real(Lsm1(-s12,-s124,-s24,-s124))*fun2(s12,s24,s23,s14,s13,s34)
     & -real(Lsm1(-s12,-s124,-s14,-s124))*fun2(s12,s24,s23,s14,s13,s34)
     & -real(Lsm1(-s14,-s124,-s12,-s124))*fun2(s14,s24,s34,s12,s13,s23)
     & -real(Lsm1(-s34,-s234,-s23,-s234))*fun2(s23,s24,s12,s34,s13,s14)
     & -real(Lsm1(-s23,-s234,-s34,-s234))*fun2(s34,s24,s14,s23,s13,s12)
     & -real(Lsm1(-s34,-s134,-s14,-s134))*fun2(s34,s13,s23,s14,s24,s12)
     & -real(Lsm1(-s34,-s134,-s14,-s134))*fun2(s14,s13,s12,s34,s24,s23))


      Hqaqavsq = Hqaqavsq+V*(1._dp+1._dp/xn**2)*( 
     & + real(L1(-s134,-s34))*fun3(s34,s13,s23,s14,s24,s12)
     & + real(L1(-s134,-s14))*fun3(s14,s13,s12,s34,s24,s23)
     & + real(L1(-s234,-s34))*fun3(s34,s24,s14,s23,s13,s12)
     & + real(L1(-s124,-s14))*fun3(s14,s24,s34,s12,s13,s23)
     & + real(L1(-s123,-s12))*fun3(s12,s13,s14,s23,s24,s34)
     & + real(L1(-s123,-s23))*fun3(s23,s13,s34,s12,s24,s14)
     & + real(L1(-s124,-s12))*fun3(s12,s24,s23,s14,s13,s34)
     & + real(L1(-s234,-s23))*fun3(s23,s24,s12,s34,s13,s14))

      Hqaqavsq = Hqaqavsq+V*( 
     & + real(L0(-s134,-s34))*fun4(s34,s13,s23,s14,s24,s12)
     & + real(L0(-s134,-s14))*fun4(s14,s13,s12,s34,s24,s23)
     & + real(L0(-s234,-s34))*fun4(s34,s24,s14,s23,s13,s12)
     & + real(L0(-s124,-s14))*fun4(s14,s24,s34,s12,s13,s23)
     & + real(L0(-s123,-s12))*fun4(s12,s13,s14,s23,s24,s34)
     & + real(L0(-s123,-s23))*fun4(s23,s13,s34,s12,s24,s14)
     & + real(L0(-s124,-s12))*fun4(s12,s24,s23,s14,s13,s34)
     & + real(L0(-s234,-s23))*fun4(s23,s24,s12,s34,s13,s14))


      Hqaqavsq = Hqaqavsq+V/xn**2*( 
     & + real(L0(-s134,-s34))*fun5(s34,s13,s23,s14,s24,s12)
     & + real(L0(-s134,-s14))*fun5(s14,s13,s12,s34,s24,s23)
     & + real(L0(-s234,-s34))*fun5(s34,s24,s14,s23,s13,s12)
     & + real(L0(-s124,-s14))*fun5(s14,s24,s34,s12,s13,s23)
     & + real(L0(-s123,-s12))*fun5(s12,s13,s14,s23,s24,s34)
     & + real(L0(-s123,-s23))*fun5(s23,s13,s34,s12,s24,s14)
     & + real(L0(-s124,-s12))*fun5(s12,s24,s23,s14,s13,s34)
     & + real(L0(-s234,-s23))*fun5(s23,s24,s12,s34,s13,s14))

      Hqaqavsq = Hqaqavsq+V*( 
     & +real(lnrat(-s134,-s34))*fun6(s34,s13,s23,s14,s24,s12)
     & +real(lnrat(-s134,-s14))*fun6(s14,s13,s12,s34,s24,s23)
     & +real(lnrat(-s234,-s34))*fun6(s34,s24,s14,s23,s13,s12)
     & +real(lnrat(-s124,-s14))*fun6(s14,s24,s34,s12,s13,s23)
     & +real(lnrat(-s123,-s12))*fun6(s12,s13,s14,s23,s24,s34)
     & +real(lnrat(-s123,-s23))*fun6(s23,s13,s34,s12,s24,s14)
     & +real(lnrat(-s124,-s12))*fun6(s12,s24,s23,s14,s13,s34)
     & +real(lnrat(-s234,-s23))*fun6(s23,s24,s12,s34,s13,s14))


      Hqaqavsq = Hqaqavsq+V/xn**2*( 
     & +real(lnrat(-s134,-s34))*fun7(s34,s13,s23,s14,s24,s12)
     & +real(lnrat(-s134,-s14))*fun7(s14,s13,s12,s34,s24,s23)
     & +real(lnrat(-s234,-s34))*fun7(s34,s24,s14,s23,s13,s12)
     & +real(lnrat(-s124,-s14))*fun7(s14,s24,s34,s12,s13,s23)
     & +real(lnrat(-s123,-s12))*fun7(s12,s13,s14,s23,s24,s34)
     & +real(lnrat(-s123,-s23))*fun7(s23,s13,s34,s12,s24,s14)
     & +real(lnrat(-s124,-s12))*fun7(s12,s24,s23,s14,s13,s34)
     & +real(lnrat(-s234,-s23))*fun7(s23,s24,s12,s34,s13,s14))


      Hqaqavsq = Hqaqavsq+V*(1._dp+1._dp/xn**2)*( 
     & +fun8(s34,s13,s23,s14,s24,s12)
     & +fun8(s14,s13,s12,s34,s24,s23)
     & +fun8(s34,s24,s14,s23,s13,s12)
     & +fun8(s14,s24,s34,s12,s13,s23)
     & +fun8(s12,s13,s14,s23,s24,s34)
     & +fun8(s23,s13,s34,s12,s24,s14)
     & +fun8(s12,s24,s23,s14,s13,s34)
     & +fun8(s23,s24,s12,s34,s13,s14))

c--- translation between schemes, according to Eq. (54) of EGZ
      if (scheme == 'dred') then
      Hqaqavsq=Hqaqavsq+(2._dp*Cf)*Hqaqa
      endif


c      write(6,*) 'Hqaqa:born',Hqaqa
c      write(6,*) 'final:Hqaqavsq',Hqaqavsq

      return
      end
