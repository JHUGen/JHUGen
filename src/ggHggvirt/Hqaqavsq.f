      double precision function Hqaqavsq(i1,i2,i3,i4)
      implicit none
C     Amplitude squared with certain factors removed
C     for the process H --> q(p1)+a(p2)+q(p3)+a(p4)
 
      include 'constants.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'scheme.f'
      double complex L0,L1,Lsm1,Lsm1_2me,lnrat
      double precision s123,s124,s134,s234,s12,s13,s14,s23,s24,s34
      double precision fun1,fun2,fun3,fun4,fun5,fun6,fun7,fun8
      double precision Hqaqa,Hqaqasq,mhsq
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


      fun1(s12,s13,s14,s23,s24,s34)=(
     &  - 2.D0*s12*s13*s24*s34+ s12*s13**2*s34
     &  - 4.D0*s12*s14*s23*s34+ s12*s24**2*s34
     &  + 2.D0*s12**2*s34**2- 2.D0*s13*s14*s23*s24
     &  - s13*s24**3+ s13**2*s14*s23
     &  - s13**3*s24+ s14*s23*s24**2+ 2.D0*s14**2*s23**2)
     & /(4d0*s12*s23*s14*s34)

      fun2(s12,s13,s14,s23,s24,s34)=(
     &  - 2.D0*s12*s13*s14*s23*s24*s34
     &  + 2.D0*s12*s13**2*s14*s23*s34
     &  - 4.D0*s12*s13**2*s24**2*s34
     &  + 3.D0*s12*s13**3*s24*s34
     &  - 2.D0*s12*s13**4*s34
     &  + s12*s14**2*s23**2*s34
     &  + 3.D0*s12**2*s13*s24*s34**2
     &  - 3.D0*s12**2*s13**2*s34**2
     &  + s12**2*s14*s23*s34**2
     &  - s12**3*s34**3
     &  + 3.D0*s13*s14**2*s23**2*s24
     &  - 4.D0*s13**2*s14*s23*s24**2
     &  - 3.D0*s13**2*s14**2*s23**2
     &  + 3.D0*s13**3*s14*s23*s24
     &  + 2.D0*s13**3*s24**3
     &  - 2.D0*s13**4*s14*s23
     &  + 2.D0*s13**5*s24
     &  - s14**3*s23**3)
     & /(8d0*s12*s13**2*s23*s14*s34)

      fun3(s12,s13,s14,s23,s24,s34)=(
     &  - 2.D0*s12*s14*s23*s34
     &  + s12*s23*s24*s34
     &  + s12*s23**2*s34
     &  + s12**2*s34**2
     &  + s14*s23**2*s24
     &  + s14*s23**3
     &  + s14**2*s23**2
     &  + s23**3*s24)
     & /(8d0*s12**2*s14*s23)

      fun4(s12,s13,s14,s23,s24,s34)=(
     &  - 2.D0*s12*s23*s34+ 5.D0*s12*s24*s34- 3.D0*s14*s23*s24
     &  + 6.D0*s14*s23**2+ 4.D0*s23*s24**2- s23**2*s24)
     & /(8d0*s14*s12*s23)

      fun5(s12,s13,s14,s23,s24,s34)=(
     &  + s12*s14*s23**2*s24*s34- 6.D0*s12*s14*s23**3*s34
     &  + 2.D0*s12*s14**2*s23**2*s34- 2.D0*s12*s23**2*s24**2*s34
     &  + 3.D0*s12*s23**3*s24*s34+ 2.D0*s12**2*s14*s23*s34**2
     &  - 3.D0*s12**2*s23*s24*s34**2+ 2.D0*s12**2*s23**2*s34**2
     &  - 2.D0*s12**3*s34**3- 2.D0*s14*s23**3*s24**2
     &  - 4.D0*s14**2*s23**3*s24- 2.D0*s14**3*s23**3)
     & /(8d0*s12**2*s23**2*s14*s34)


      fun6(s12,s13,s14,s23,s24,s34)=(- 3.D0*s12*s34+ 3.D0*s13*s24
     &  + 5.D0*s14*s23- 2.D0*s23*s24+ 4.D0*s24**2)
     & /(8d0*s14*s23)

      fun7(s12,s13,s14,s23,s24,s34)=(
     &  - 7.D0*s12*s13*s14*s23**2*s34
     &  - 2.D0*s12*s13*s23*s24**2*s34
     &  + 2.D0*s12*s13*s23**2*s24*s34
     &  - s12*s13**2*s23*s24*s34
     &  + 2.D0*s12*s14**2*s23**2*s34
     &  + s12**2*s13*s23*s34**2
     &  + 2.D0*s12**2*s14*s23*s34**2
     &  - 2.D0*s12**3*s34**3
     &  - 2.D0*s13*s14*s23**2*s24**2
     &  - 2.D0*s14**3*s23**3)
     & /(8d0*s12*s13*s23**2*s14*s34)

      fun8(s12,s13,s14,s23,s24,s34)=(
     &  - 2.D0*s12*s13*s14*s24- 2.D0*s12*s13*s23*s24
     &  + s12*s13*s23**2+ s12*s13*s24**2
     &  + s12*s13**2*s24+ 2.D0*s12*s14*s23**2
     &  + 2.D0*s12*s14*s34**2+ 2.D0*s12*s14**2*s23
     &  + s12*s14**2*s24+ 2.D0*s12*s23*s34**2+s12**2*s13*s23
     &  + s12**2*s14*s24+ 2.D0*s12**2*s14*s34
     &  + 2.D0*s12**2*s23*s34- 2.D0*s13*s14*s24*s34+ s13*s14*s24**2
     &  + s13*s14*s34**2 + s13*s14**2*s34
     &  - 2.D0*s13*s23*s24*s34+ s13*s23*s24**2
     &  + s13*s24**2*s34+ s13**2*s14*s24
     &  + s13**2*s23*s24+ s13**2*s24*s34
     &  + 2.D0*s14*s23**2*s34+ 2.D0*s14**2*s23*s34
     &  + s23*s24*s34**2  + s23**2*s24*s34)
     &  /(64d0*s12*s14*s23*s34)


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


      Hqaqavsq = + epinv**2*(+2.D0/xn-2.D0*xn)
      Hqaqavsq = Hqaqavsq + epinv * (
     &     + 3.D0/xn- 3.D0*xn
     &     - 8.D0/3.D0*0.5d0*nf
     &     + 4.D0/3.D0*nf
     &     - dble(lnrat( - s14,musq))/xn
     &     - dble(lnrat( - s23,musq))/xn
     &     - dble(lnrat( - s12,musq))/xn
     &     - dble(lnrat( - s34,musq))/xn
     &     + dble(lnrat( - s24,musq))*(xn+1d0/xn)
     &     + dble(lnrat( - s13,musq))*(xn+1d0/xn)
     &     )
      Hqaqavsq = Hqaqavsq + 11.D0
     &     + 8.D0/xn
     &     + 80.D0/9.D0*xn
     &     - 40.D0/9.D0*0.5d0*nf
     &     + dble(Lsm1_2me(s123,s234,s23,mhsq))/xn
     &     + dble(Lsm1_2me(s123,s124,s12,mhsq))/xn
     &     + dble(Lsm1_2me(s134,s234,s34,mhsq))/xn
     &     + dble(Lsm1_2me(s124,s134,s14,mhsq))/xn
     &     - dble(Lsm1_2me(s123,s134,s13,mhsq))*(xn+1d0/xn)
     &     - dble(Lsm1_2me(s124,s234,s24,mhsq))*(xn+1d0/xn)
     &     - 3.D0/4.D0*dble(lnrat( - s14,musq))/xn
     &     - 13.D0/12.D0*dble(lnrat( - s14,musq))*xn
     &     + 2.D0/3.D0*dble(lnrat( - s14,musq))*0.5d0*nf
     &
      Hqaqavsq = Hqaqavsq 
     &     + 1.D0/2.D0*dble(lnrat( - s14,musq)**2)/xn
     &     - 3.D0/4.D0*dble(lnrat( - s23,musq))/xn
     &     - 13.D0/12.D0*dble(lnrat( - s23,musq))*xn
     &     + 2.D0/3.D0*dble(lnrat( - s23,musq))*0.5d0*nf
     &     + 1.D0/2.D0*dble(lnrat( - s23,musq)**2)/xn
     &     - 3.D0/4.D0*dble(lnrat( - s12,musq))/xn
     &     - 13.D0/12.D0*dble(lnrat( - s12,musq))*xn
     &     + 2.D0/3.D0*dble(lnrat( - s12,musq))*0.5d0*nf
     &     + 1.D0/2.D0*dble(lnrat( - s12,musq)**2)/xn
     &     - 3.D0/4.D0*dble(lnrat( - s34,musq))/xn
     &     - 13.D0/12.D0*dble(lnrat( - s34,musq))*xn
     &     + 2.D0/3.D0*dble(lnrat( - s34,musq))*0.5d0*nf
     &     + 1.D0/2.D0*dble(lnrat( - s34,musq)**2)/xn
     &     - 1.D0/2.D0*dble(lnrat( - s24,musq)**2)*(xn+1d0/xn)
     &     - 1.D0/2.D0*dble(lnrat( - s13,musq)**2)*(xn+1d0/xn)

      if (checkEGZ) then
c--- remove finite renormalization
      Hqaqavsq = Hqaqavsq - 11.D0     
      endif
    
      Hqaqavsq =Hqaqavsq*Hqaqa
      
      Hqaqavsq = Hqaqavsq+V*(1d0+1d0/xn**2)*( 
     & -dble(Lsm1(-s12,-s124,-s24,-s124))*fun1(s12,s13,s14,s23,s24,s34)
     & -dble(Lsm1(-s23,-s234,-s24,-s234))*fun1(s23,s13,s34,s12,s24,s14)
     & -dble(Lsm1(-s14,-s124,-s24,-s124))*fun1(s14,s13,s12,s34,s24,s23)
     & -dble(Lsm1(-s34,-s234,-s24,-s234))*fun1(s34,s13,s23,s14,s24,s12)
     & -dble(Lsm1(-s34,-s134,-s13,-s134))*fun1(s12,s13,s14,s23,s24,s34)
     & -dble(Lsm1(-s14,-s134,-s13,-s134))*fun1(s23,s13,s34,s12,s24,s14)
     & -dble(Lsm1(-s12,-s123,-s13,-s123))*fun1(s14,s13,s12,s34,s24,s23)
     & -dble(Lsm1(-s23,-s123,-s13,-s123))*fun1(s34,s13,s23,s14,s24,s12))

c--- NB: in sub-formula below, I have replaced the line:
c     & -dble(Lsm1(-s12,-s124,-s24,-s124))*fun2(s12,s24,s23,s14,s13,s34)
c---     with the one:
c     & -dble(Lsm1(-s12,-s124,-s14,-s124))*fun2(s12,s24,s23,s14,s13,s34)
c--- This is as prescribed by Eq. (37) of EGZ and restores 1<>3 symmetry

      Hqaqavsq = Hqaqavsq+V/xn**2*( 
     & -dble(Lsm1(-s12,-s123,-s23,-s123))*fun2(s12,s13,s14,s23,s24,s34)
     & -dble(Lsm1(-s23,-s123,-s12,-s123))*fun2(s23,s13,s34,s12,s24,s14)
c     & -dble(Lsm1(-s12,-s124,-s24,-s124))*fun2(s12,s24,s23,s14,s13,s34)
     & -dble(Lsm1(-s12,-s124,-s14,-s124))*fun2(s12,s24,s23,s14,s13,s34)
     & -dble(Lsm1(-s14,-s124,-s12,-s124))*fun2(s14,s24,s34,s12,s13,s23)
     & -dble(Lsm1(-s34,-s234,-s23,-s234))*fun2(s23,s24,s12,s34,s13,s14)
     & -dble(Lsm1(-s23,-s234,-s34,-s234))*fun2(s34,s24,s14,s23,s13,s12)
     & -dble(Lsm1(-s34,-s134,-s14,-s134))*fun2(s34,s13,s23,s14,s24,s12)
     & -dble(Lsm1(-s34,-s134,-s14,-s134))*fun2(s14,s13,s12,s34,s24,s23))


      Hqaqavsq = Hqaqavsq+V*(1d0+1d0/xn**2)*( 
     & + dble(L1(-s134,-s34))*fun3(s34,s13,s23,s14,s24,s12)
     & + dble(L1(-s134,-s14))*fun3(s14,s13,s12,s34,s24,s23)
     & + dble(L1(-s234,-s34))*fun3(s34,s24,s14,s23,s13,s12)
     & + dble(L1(-s124,-s14))*fun3(s14,s24,s34,s12,s13,s23)
     & + dble(L1(-s123,-s12))*fun3(s12,s13,s14,s23,s24,s34)
     & + dble(L1(-s123,-s23))*fun3(s23,s13,s34,s12,s24,s14)
     & + dble(L1(-s124,-s12))*fun3(s12,s24,s23,s14,s13,s34)
     & + dble(L1(-s234,-s23))*fun3(s23,s24,s12,s34,s13,s14))

      Hqaqavsq = Hqaqavsq+V*( 
     & + dble(L0(-s134,-s34))*fun4(s34,s13,s23,s14,s24,s12)
     & + dble(L0(-s134,-s14))*fun4(s14,s13,s12,s34,s24,s23)
     & + dble(L0(-s234,-s34))*fun4(s34,s24,s14,s23,s13,s12)
     & + dble(L0(-s124,-s14))*fun4(s14,s24,s34,s12,s13,s23)
     & + dble(L0(-s123,-s12))*fun4(s12,s13,s14,s23,s24,s34)
     & + dble(L0(-s123,-s23))*fun4(s23,s13,s34,s12,s24,s14)
     & + dble(L0(-s124,-s12))*fun4(s12,s24,s23,s14,s13,s34)
     & + dble(L0(-s234,-s23))*fun4(s23,s24,s12,s34,s13,s14))


      Hqaqavsq = Hqaqavsq+V/xn**2*( 
     & + dble(L0(-s134,-s34))*fun5(s34,s13,s23,s14,s24,s12)
     & + dble(L0(-s134,-s14))*fun5(s14,s13,s12,s34,s24,s23)
     & + dble(L0(-s234,-s34))*fun5(s34,s24,s14,s23,s13,s12)
     & + dble(L0(-s124,-s14))*fun5(s14,s24,s34,s12,s13,s23)
     & + dble(L0(-s123,-s12))*fun5(s12,s13,s14,s23,s24,s34)
     & + dble(L0(-s123,-s23))*fun5(s23,s13,s34,s12,s24,s14)
     & + dble(L0(-s124,-s12))*fun5(s12,s24,s23,s14,s13,s34)
     & + dble(L0(-s234,-s23))*fun5(s23,s24,s12,s34,s13,s14))

      Hqaqavsq = Hqaqavsq+V*( 
     & +dble(lnrat(-s134,-s34))*fun6(s34,s13,s23,s14,s24,s12)
     & +dble(lnrat(-s134,-s14))*fun6(s14,s13,s12,s34,s24,s23)
     & +dble(lnrat(-s234,-s34))*fun6(s34,s24,s14,s23,s13,s12)
     & +dble(lnrat(-s124,-s14))*fun6(s14,s24,s34,s12,s13,s23)
     & +dble(lnrat(-s123,-s12))*fun6(s12,s13,s14,s23,s24,s34)
     & +dble(lnrat(-s123,-s23))*fun6(s23,s13,s34,s12,s24,s14)
     & +dble(lnrat(-s124,-s12))*fun6(s12,s24,s23,s14,s13,s34)
     & +dble(lnrat(-s234,-s23))*fun6(s23,s24,s12,s34,s13,s14))


      Hqaqavsq = Hqaqavsq+V/xn**2*( 
     & +dble(lnrat(-s134,-s34))*fun7(s34,s13,s23,s14,s24,s12)
     & +dble(lnrat(-s134,-s14))*fun7(s14,s13,s12,s34,s24,s23)
     & +dble(lnrat(-s234,-s34))*fun7(s34,s24,s14,s23,s13,s12)
     & +dble(lnrat(-s124,-s14))*fun7(s14,s24,s34,s12,s13,s23)
     & +dble(lnrat(-s123,-s12))*fun7(s12,s13,s14,s23,s24,s34)
     & +dble(lnrat(-s123,-s23))*fun7(s23,s13,s34,s12,s24,s14)
     & +dble(lnrat(-s124,-s12))*fun7(s12,s24,s23,s14,s13,s34)
     & +dble(lnrat(-s234,-s23))*fun7(s23,s24,s12,s34,s13,s14))


      Hqaqavsq = Hqaqavsq+V*(1d0+1d0/xn**2)*( 
     & +fun8(s34,s13,s23,s14,s24,s12)
     & +fun8(s14,s13,s12,s34,s24,s23)
     & +fun8(s34,s24,s14,s23,s13,s12)
     & +fun8(s14,s24,s34,s12,s13,s23)
     & +fun8(s12,s13,s14,s23,s24,s34)
     & +fun8(s23,s13,s34,s12,s24,s14)
     & +fun8(s12,s24,s23,s14,s13,s34)
     & +fun8(s23,s24,s12,s34,s13,s14))

c--- translation between schemes, according to Eq. (54) of EGZ
      if (scheme .eq. 'dred') then
      Hqaqavsq=Hqaqavsq+(2d0*Cf)*Hqaqa
      endif


c      write(6,*) 'Hqaqa:born',Hqaqa
c      write(6,*) 'final:Hqaqavsq',Hqaqavsq

      return
      end
