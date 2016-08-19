      subroutine Dfill_recur3(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N0)
      implicit none
C     Implements the calculation of the formfactors
C     for small momenta, as in DD Eq.5.64-5.70 etc
C     N0 is the offset in the common block

C--- Currently: calculates up to rank 4 with at least one recursion
c---            calculates rank 5 with no recursion
c---            calculates D00iiii components of rank 6

      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      integer C234,C134,C124,C123,np,ep,N0,j,k,pvCcache,
     . i1,i2,i3,i4,i5,step,kmax
      parameter(np=3)
      double precision p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,f(np),
     . Gr(np,np),DetGr
      double complex S0000(-2:0),S0000i(np,-2:0),
     . Shat3zz(np,-2:0),Shat4zz(np,z1max,-2:0),
     . Shat5zz(np,z2max,-2:0),Shat6zz(np,z3max,-2:0),
     . Shat5zzzz(np,-2:0),Shat6zzzz(np,z1max,-2:0),
     . Shat7zz(np,z4max,-2:0),Shat7zzzzzz(np,-2:0),
     . Shat1(np,-2:0),Shat2(np,z1max,-2:0),
     . Shat3(np,z2max,-2:0),Shat4(np,z3max,-2:0),Shat5(np,z4max,-2:0),
     . Shat6(np,z5max,-2:0),Shat7(np,z6max,-2:0)
      double complex csum0(-2:0),csum1(-2:0),csum2(-2:0),
     . csum11(-2:0),csum00(-2:0),csum12(-2:0),csum22(-2:0),
     . csum111(-2:0),csum112(-2:0),csum122(-2:0),csum222(-2:0),
     . csum1111(-2:0),csum1112(-2:0),csum1122(-2:0),
     . csum1222(-2:0),csum2222(-2:0),
     . csum11111(-2:0),csum11112(-2:0),csum11122(-2:0),
     . csum11222(-2:0),csum12222(-2:0),csum22222(-2:0),
     . csum001(-2:0),csum002(-2:0),csum0011(-2:0),
     . csum0012(-2:0),csum0022(-2:0),csum0000(-2:0),
     . csum00111(-2:0),csum00112(-2:0),csum00122(-2:0),csum00222(-2:0),
     . Czero5(z5max,-2:0),Czero4(z4max,-2:0),
     . Czero3(z3max,-2:0),Czero2(z2max,-2:0),Czero1(z1max,-2:0),
     . Czero0(-2:0)
      logical,save:: first=.true.
!$omp threadprivate(first)

      if (first) then
        first=.false.
c--- These lines must be uncommented
c        call Array3dim
c        call DArraysetup
        write(6,*) 'Need to uncomment lines in Dfill_recur3'
        stop
      endif

c--- Not necessary, routine upgraded now
c      if ((m1 .ne. 0d0).or.(m2 .ne. 0d0).or.(m3 .ne. 0d0)) then
c      write(6,*) 'Dfill: nonzero internal masses not foreseen'
c      stop
c      endif

      C234=pvCcache(p2,p3,p2p3,m2,m3,m4)
      C134=pvCcache(p1p2,p3,p4,m1,m3,m4)
      C124=pvCcache(p1,p2p3,p4,m1,m2,m4)
      C123=pvCcache(p1,p2,p1p2,m1,m2,m3)
       
C----We have changed the sign of fi (different from Dfill) to agree
C----with notation of Denner-Dittmaier
      f(1) = -m2 + m1 + p1
      f(2) = -m3 + m1 + p1p2
      f(3) = -m4 + m1 + p4

      Gr(1,1) = 2*(p1)
      Gr(2,2) = 2*(p1p2)
      Gr(3,3) = 2*(p4)
      Gr(1,2) = (p1+p1p2 - p2)
      Gr(2,1) = Gr(1,2)
      Gr(1,3) = (p1+p4 - p2p3)
      Gr(3,1) = Gr(1,3)
      Gr(2,3) = (p1p2 - p3+p4)
      Gr(3,2) = Gr(2,3)
      call determinant(3,np,Gr,DetGr)
      write(6,*) 'small P: 3x3 DetGr = ',DetGr
      
      do ep=-2,0
      csum0(ep)=Cv(cc0+C234,ep)+Cv(cc1+C234,ep)+Cv(cc2+C234,ep)
      csum1(ep)=Cv(cc1+C234,ep)+Cv(cc11+C234,ep)+Cv(cc12+C234,ep)
      csum2(ep)=Cv(cc2+C234,ep)+Cv(cc12+C234,ep)+Cv(cc22+C234,ep)
      enddo
      do ep=-2,0
      csum00(ep)=Cv(cc00+C234,ep)+Cv(cc001+C234,ep)+Cv(cc002+C234,ep)
      csum11(ep)=Cv(cc11+C234,ep)+Cv(cc111+C234,ep)+Cv(cc112+C234,ep)
      csum12(ep)=Cv(cc12+C234,ep)+Cv(cc112+C234,ep)+Cv(cc122+C234,ep)
      csum22(ep)=Cv(cc22+C234,ep)+Cv(cc122+C234,ep)+Cv(cc222+C234,ep)
      csum111(ep)=
     & Cv(cc111+C234,ep)+Cv(cc1111+C234,ep)+Cv(cc1112+C234,ep)
      csum112(ep)=
     & Cv(cc112+C234,ep)+Cv(cc1112+C234,ep)+Cv(cc1122+C234,ep)
      csum122(ep)=
     & Cv(cc122+C234,ep)+Cv(cc1122+C234,ep)+Cv(cc1222+C234,ep)
      csum222(ep)=
     & Cv(cc222+C234,ep)+Cv(cc1222+C234,ep)+Cv(cc2222+C234,ep)

      csum001(ep)=
     & Cv(cc001+C234,ep)+Cv(cc0011+C234,ep)+Cv(cc0012+C234,ep)
      csum002(ep)=
     & Cv(cc002+C234,ep)+Cv(cc0012+C234,ep)+Cv(cc0022+C234,ep)

      csum0000(ep)=Cv(cc0000+C234,ep)+Cv(cc00001+C234,ep)
     & +Cv(cc00002+C234,ep)
      csum0011(ep)=
     & Cv(cc0011+C234,ep)+Cv(cc00111+C234,ep)+Cv(cc00112+C234,ep)
      csum0012(ep)=
     & Cv(cc0012+C234,ep)+Cv(cc00112+C234,ep)+Cv(cc00122+C234,ep)
      csum0022(ep)=
     & Cv(cc0022+C234,ep)+Cv(cc00122+C234,ep)+Cv(cc00222+C234,ep)

      csum1111(ep)=
     & Cv(cc1111+C234,ep)+Cv(cc11111+C234,ep)+Cv(cc11112+C234,ep)
      csum1112(ep)=
     & Cv(cc1112+C234,ep)+Cv(cc11112+C234,ep)+Cv(cc11122+C234,ep)
      csum1122(ep)=
     & Cv(cc1122+C234,ep)+Cv(cc11122+C234,ep)+Cv(cc11222+C234,ep)
      csum1222(ep)=
     & Cv(cc1222+C234,ep)+Cv(cc11222+C234,ep)+Cv(cc12222+C234,ep)
      csum2222(ep)=
     & Cv(cc2222+C234,ep)+Cv(cc12222+C234,ep)+Cv(cc22222+C234,ep)

      csum11111(ep)=
     & Cv(cc11111+C234,ep)+Cv(cc111111+C234,ep)+Cv(cc111112+C234,ep)
      csum11112(ep)=
     & Cv(cc11112+C234,ep)+Cv(cc111112+C234,ep)+Cv(cc111122+C234,ep)
      csum11122(ep)=
     & Cv(cc11122+C234,ep)+Cv(cc111122+C234,ep)+Cv(cc111222+C234,ep)
      csum11222(ep)=
     & Cv(cc11222+C234,ep)+Cv(cc111222+C234,ep)+Cv(cc112222+C234,ep)
      csum12222(ep)=
     & Cv(cc12222+C234,ep)+Cv(cc112222+C234,ep)+Cv(cc122222+C234,ep)
      csum22222(ep)=
     & Cv(cc22222+C234,ep)+Cv(cc122222+C234,ep)+Cv(cc222222+C234,ep)

      csum00111(ep)=
     & Cv(cc00111+C234,ep)+Cv(cc001111+C234,ep)+Cv(cc001112+C234,ep)
      csum00112(ep)=
     & Cv(cc00112+C234,ep)+Cv(cc001112+C234,ep)+Cv(cc001122+C234,ep)
      csum00122(ep)=
     & Cv(cc00122+C234,ep)+Cv(cc001122+C234,ep)+Cv(cc001222+C234,ep)
      csum00222(ep)=
     & Cv(cc00222+C234,ep)+Cv(cc001222+C234,ep)+Cv(cc002222+C234,ep)
      enddo

      do ep=-2,0
      include 'Shat.f'
      enddo

c--- Note: these functions are not used in this recursion
c--- definitions of the S00 functions from Dfill_recur
c      include 'S00_def.f'
c      include 'S00i_def.f'
c      include 'S0000_def.f'
c      include 'S00ii_def.f'
c      include 'S0000i_def.f'
c      include 'S00iii_def.f'
c      include 'S000000_def.f'
c      include 'S0000ii_def.f'
c      include 'S00iiii_def.f'

      do ep=-2,0

      Shat3zz(1,ep)=Cv(cc00+C134,ep)-Cv(cc00+C234,ep)
      Shat3zz(2,ep)=Cv(cc00+C124,ep)-Cv(cc00+C234,ep)
      Shat3zz(3,ep)=Cv(cc00+C123,ep)-Cv(cc00+C234,ep)

      Shat4zz(1,1,ep)=Csum00(ep)
      Shat4zz(2,1,ep)=Csum00(ep)+Cv(cc001+C124,ep)
      Shat4zz(3,1,ep)=Csum00(ep)+Cv(cc001+C123,ep)

      Shat4zz(1,2,ep)=+Cv(cc001+C134,ep)-Cv(cc001+C234,ep)
      Shat4zz(2,2,ep)=-Cv(cc001+C234,ep)
      Shat4zz(3,2,ep)=-Cv(cc001+C234,ep)+Cv(cc002+C123,ep)

      Shat4zz(1,3,ep)=+Cv(cc002+C134,ep)-Cv(cc002+C234,ep)
      Shat4zz(2,3,ep)=+Cv(cc002+C124,ep)-Cv(cc002+C234,ep)
      Shat4zz(3,3,ep)=-Cv(cc002+C234,ep)


      Shat5zz(1,z2(1,1),ep)=-Csum00(ep)-Csum001(ep)-Csum002(ep)
      Shat5zz(2,z2(1,1),ep)=-Csum00(ep)-Csum001(ep)-Csum002(ep)
     & +Cv(cc0011+C124,ep)
      Shat5zz(3,z2(1,1),ep)=-Csum00(ep)-Csum001(ep)-Csum002(ep)
     & +Cv(cc0011+C123,ep)

      Shat5zz(1,z2(1,2),ep)=Csum001(ep)
      Shat5zz(2,z2(1,2),ep)=Csum001(ep)
      Shat5zz(3,z2(1,2),ep)=Csum001(ep)+ Cv(cc0012+C123,ep)

      Shat5zz(1,z2(2,2),ep)=Cv(cc0011+C134,ep)- Cv(cc0011+C234,ep)
      Shat5zz(2,z2(2,2),ep)=- Cv(cc0011+C234,ep)
      Shat5zz(3,z2(2,2),ep)=Cv(cc0022+C123,ep)- Cv(cc0011+C234,ep)

      Shat5zz(1,z2(2,3),ep)=Cv(cc0012+C134,ep)- Cv(cc0012+C234,ep)
      Shat5zz(2,z2(2,3),ep)=-Cv(cc0012+C234,ep)
      Shat5zz(3,z2(2,3),ep)=-Cv(cc0012+C234,ep)

      Shat5zz(1,z2(1,3),ep)=Csum002(ep)
      Shat5zz(2,z2(1,3),ep)=Csum002(ep)+ Cv(cc0012+C124,ep)
      Shat5zz(3,z2(1,3),ep)=+Csum002(ep)

      Shat5zz(1,z2(3,3),ep)=Cv(cc0022+C134,ep)- Cv(cc0022+C234,ep)
      Shat5zz(2,z2(3,3),ep)=Cv(cc0022+C124,ep)- Cv(cc0022+C234,ep)
      Shat5zz(3,z2(3,3),ep)=-Cv(cc0022+C234,ep)

      Shat6zz(1,z3(1,1,1),ep)=Csum0011(ep)+Csum0022(ep)+2D0*Csum0012(ep)
     & +Csum00(ep) + 2D0*Csum001(ep)+ 2D0*Csum002(ep)
      Shat6zz(2,z3(1,1,1),ep)=Csum0011(ep)+Csum0022(ep)+2D0*Csum0012(ep)
     & +2.D0*Csum002(ep)+2.D0*Csum001(ep)+Csum00(ep)
     & +Cv(cc00111+C124,ep)
      Shat6zz(3,z3(1,1,1),ep)=Csum0011(ep)+Csum0022(ep)+2D0*Csum0012(ep)
     & +2.D0*Csum002(ep)+2.D0*Csum001(ep)+Csum00(ep)
     & +Cv(cc00111+C123,ep)

      Shat6zz(1,z3(1,1,2),ep)=-Csum0011(ep)-Csum0012(ep)-Csum001(ep)
      Shat6zz(2,z3(1,1,2),ep)=-Csum0011(ep)-Csum0012(ep)-Csum001(ep)
      Shat6zz(3,z3(1,1,2),ep)=-Csum0011(ep)-Csum0012(ep)-Csum001(ep)
     & +Cv(cc00112+C123,ep)

      Shat6zz(1,z3(1,1,3),ep)=-Csum0022(ep)-Csum0012(ep)-Csum002(ep)
      Shat6zz(2,z3(1,1,3),ep)=-Csum0022(ep)-Csum0012(ep)
     & -Csum002(ep)+Cv(cc00112+C124,ep)
      Shat6zz(3,z3(1,1,3),ep)=-Csum0022(ep)-Csum0012(ep)-Csum002(ep)

      Shat6zz(1,z3(1,2,2),ep)=+Csum0011(ep)
      Shat6zz(2,z3(1,2,2),ep)=+Csum0011(ep)
      Shat6zz(3,z3(1,2,2),ep)=+Csum0011(ep)+Cv(cc00122+C123,ep)

      Shat6zz(1,z3(1,2,3),ep)=+Csum0012(ep)
      Shat6zz(2,z3(1,2,3),ep)=+Csum0012(ep)
      Shat6zz(3,z3(1,2,3),ep)=+Csum0012(ep)

      Shat6zz(1,z3(1,3,3),ep)=+Csum0022(ep)
      Shat6zz(2,z3(1,3,3),ep)=+Csum0022(ep)+Cv(cc00122+C124,ep)
      Shat6zz(3,z3(1,3,3),ep)=+Csum0022(ep)

      Shat6zz(1,z3(2,2,2),ep)=+Cv(cc00111+C134,ep)-Cv(cc00111+C234,ep)
      Shat6zz(2,z3(2,2,2),ep)=-Cv(cc00111+C234,ep)
      Shat6zz(3,z3(2,2,2),ep)=-Cv(cc00111+C234,ep)+Cv(cc00222+C123,ep)

      Shat6zz(1,z3(2,2,3),ep)=+Cv(cc00112+C134,ep)-Cv(cc00112+C234,ep)
      Shat6zz(2,z3(2,2,3),ep)=-Cv(cc00112+C234,ep)
      Shat6zz(3,z3(2,2,3),ep)=-Cv(cc00112+C234,ep)

      Shat6zz(1,z3(3,3,2),ep)=+Cv(cc00122+C134,ep)-Cv(cc00122+C234,ep)
      Shat6zz(2,z3(3,3,2),ep)=-Cv(cc00122+C234,ep)
      Shat6zz(3,z3(3,3,2),ep)=-Cv(cc00122+C234,ep)

      Shat6zz(1,z3(3,3,3),ep)=+Cv(cc00222+C134,ep)-Cv(cc00222+C234,ep)
      Shat6zz(2,z3(3,3,3),ep)=+Cv(cc00222+C124,ep)-Cv(cc00222+C234,ep)
      Shat6zz(3,z3(3,3,3),ep)=-Cv(cc00222+C234,ep)

      Shat5zzzz(1,ep)=Cv(cc0000+C134,ep)-Cv(cc0000+C234,ep)
      Shat5zzzz(2,ep)=Cv(cc0000+C124,ep)-Cv(cc0000+C234,ep)
      Shat5zzzz(3,ep)=Cv(cc0000+C123,ep)-Cv(cc0000+C234,ep)

      Shat6zzzz(1,1,ep)=Csum0000(ep)
      Shat6zzzz(2,1,ep)=Csum0000(ep)+Cv(cc00001+C124,ep)
      Shat6zzzz(3,1,ep)=Csum0000(ep)+Cv(cc00001+C123,ep)

      Shat6zzzz(1,2,ep)=+Cv(cc00001+C134,ep)-Cv(cc00001+C234,ep)
      Shat6zzzz(2,2,ep)=-Cv(cc00001+C234,ep)
      Shat6zzzz(3,2,ep)=-Cv(cc00001+C234,ep)+Cv(cc00002+C123,ep)

      Shat6zzzz(1,3,ep)=+Cv(cc00002+C134,ep)-Cv(cc00002+C234,ep)
      Shat6zzzz(2,3,ep)=+Cv(cc00002+C124,ep)-Cv(cc00002+C234,ep)
      Shat6zzzz(3,3,ep)=-Cv(cc00002+C234,ep)

      Shat7zz(1,z4(1,1,1,1),ep) =
     &  - Csum00(ep)
     &  - 3.D0*Csum001(ep)
     &  - 3.D0*Csum002(ep)
     &  - 3.D0*Csum0011(ep)
     &  - 3.D0*Csum0022(ep)
     &  - 6.D0*Csum0012(ep)
     &  - Csum00111(ep)
     &  - Csum00222(ep)
     &  - 3.D0*Csum00122(ep)
     &  - 3.D0*Csum00112(ep)

      Shat7zz(2,z4(1,1,1,1),ep) =
     &  - Csum00(ep)
     &  - 3.D0*Csum001(ep)
     &  - 3.D0*Csum002(ep)
     &  - 3.D0*Csum0011(ep)
     &  - 3.D0*Csum0022(ep)
     &  - 6.D0*Csum0012(ep)
     &  - Csum00111(ep)
     &  - Csum00222(ep)
     &  - 3.D0*Csum00122(ep)
     &  - 3.D0*Csum00112(ep)
     &  + Cv(cc001111 + C124,ep)

      Shat7zz(3,z4(1,1,1,1),ep) =
     &  - Csum00(ep)
     &  - 3.D0*Csum001(ep)
     &  - 3.D0*Csum002(ep)
     &  - 3.D0*Csum0011(ep)
     &  - 3.D0*Csum0022(ep)
     &  - 6.D0*Csum0012(ep)
     &  - Csum00111(ep)
     &  - Csum00222(ep)
     &  - 3.D0*Csum00122(ep)
     &  - 3.D0*Csum00112(ep)
     &  + Cv(cc001111 + C123,ep)

      Shat7zz(1,z4(1,1,1,2),ep) =
     &  + Csum001(ep)
     &  + 2.D0*Csum0011(ep)
     &  + 2.D0*Csum0012(ep)
     &  + Csum00111(ep)
     &  + Csum00122(ep)
     &  + 2.D0*Csum00112(ep)

      Shat7zz(2,z4(1,1,1,2),ep) =
     &  + Csum001(ep)
     &  + 2.D0*Csum0011(ep)
     &  + 2.D0*Csum0012(ep)
     &  + Csum00111(ep)
     &  + Csum00122(ep)
     &  + 2.D0*Csum00112(ep)

      Shat7zz(3,z4(1,1,1,2),ep) =
     &  + Csum001(ep)
     &  + 2.D0*Csum0011(ep)
     &  + 2.D0*Csum0012(ep)
     &  + Csum00111(ep)
     &  + Csum00122(ep)
     &  + 2.D0*Csum00112(ep)
     &  + Cv(cc001112 + C123,ep)

      Shat7zz(1,z4(1,1,1,3),ep) =
     &  + Csum002(ep)
     &  + 2.D0*Csum0022(ep)
     &  + 2.D0*Csum0012(ep)
     &  + Csum00222(ep)
     &  + 2.D0*Csum00122(ep)
     &  + Csum00112(ep)

      Shat7zz(2,z4(1,1,1,3),ep) =
     &  + Csum002(ep)
     &  + 2.D0*Csum0022(ep)
     &  + 2.D0*Csum0012(ep)
     &  + Csum00222(ep)
     &  + 2.D0*Csum00122(ep)
     &  + Csum00112(ep)
     &  + Cv(cc001112 + C124,ep)

      Shat7zz(3,z4(1,1,1,3),ep) =
     &  + Csum002(ep)
     &  + 2.D0*Csum0022(ep)
     &  + 2.D0*Csum0012(ep)
     &  + Csum00222(ep)
     &  + 2.D0*Csum00122(ep)
     &  + Csum00112(ep)

      Shat7zz(1,z4(1,1,2,2),ep) =
     &  - Csum0011(ep)
     &  - Csum00111(ep)
     &  - Csum00112(ep)

      Shat7zz(2,z4(1,1,2,2),ep) =
     &  - Csum0011(ep)
     &  - Csum00111(ep)
     &  - Csum00112(ep)

      Shat7zz(3,z4(1,1,2,2),ep) =
     &  - Csum0011(ep)
     &  - Csum00111(ep)
     &  - Csum00112(ep)
     &  + Cv(cc001122 + C123,ep)

      Shat7zz(1,z4(1,1,2,3),ep) =
     &  - Csum0012(ep)
     &  - Csum00122(ep)
     &  - Csum00112(ep)

      Shat7zz(2,z4(1,1,2,3),ep) =
     &  - Csum0012(ep)
     &  - Csum00122(ep)
     &  - Csum00112(ep)

      Shat7zz(3,z4(1,1,2,3),ep) =
     &  - Csum0012(ep)
     &  - Csum00122(ep)
     &  - Csum00112(ep)

      Shat7zz(1,z4(1,1,3,3),ep) =
     &  - Csum0022(ep)
     &  - Csum00222(ep)
     &  - Csum00122(ep)

      Shat7zz(2,z4(1,1,3,3),ep) =
     &  - Csum0022(ep)
     &  - Csum00222(ep)
     &  - Csum00122(ep)
     &  + Cv(cc001122 + C124,ep)

      Shat7zz(3,z4(1,1,3,3),ep) =
     &  - Csum0022(ep)
     &  - Csum00222(ep)
     &  - Csum00122(ep)

      Shat7zz(1,z4(1,2,2,2),ep) =
     &  + Csum00111(ep)

      Shat7zz(2,z4(1,2,2,2),ep) =
     &  + Csum00111(ep)

      Shat7zz(3,z4(1,2,2,2),ep) =
     &  + Csum00111(ep)
     &  + Cv(cc001222 + C123,ep)

      Shat7zz(1,z4(1,2,2,3),ep) =
     &  + Csum00112(ep)

      Shat7zz(2,z4(1,2,2,3),ep) =
     &  + Csum00112(ep)

      Shat7zz(3,z4(1,2,2,3),ep) =
     &  + Csum00112(ep)

      Shat7zz(1,z4(1,2,3,3),ep) =
     &  + Csum00122(ep)

      Shat7zz(2,z4(1,2,3,3),ep) =
     &  + Csum00122(ep)

      Shat7zz(3,z4(1,2,3,3),ep) =
     &  + Csum00122(ep)

      Shat7zz(1,z4(1,3,3,3),ep) =
     &  + Csum00222(ep)

      Shat7zz(2,z4(1,3,3,3),ep) =
     &  + Csum00222(ep)
     &  + Cv(cc001222 + C124,ep)

      Shat7zz(3,z4(1,3,3,3),ep) =
     &  + Csum00222(ep)

      Shat7zz(1,z4(2,2,2,2),ep) =
     &  + Cv(cc001111 + C134,ep)
     &  - Cv(cc001111 + C234,ep)

      Shat7zz(2,z4(2,2,2,2),ep) =
     &  - Cv(cc001111 + C234,ep)

      Shat7zz(3,z4(2,2,2,2),ep) =
     &  - Cv(cc001111 + C234,ep)
     &  + Cv(cc002222 + C123,ep)

      Shat7zz(1,z4(2,2,2,3),ep) =
     &  + Cv(cc001112 + C134,ep)
     &  - Cv(cc001112 + C234,ep)

      Shat7zz(2,z4(2,2,2,3),ep) =
     &  - Cv(cc001112 + C234,ep)

      Shat7zz(3,z4(2,2,2,3),ep) =
     &  - Cv(cc001112 + C234,ep)

      Shat7zz(1,z4(2,2,3,3),ep) =
     &  + Cv(cc001122 + C134,ep)
     &  - Cv(cc001122 + C234,ep)

      Shat7zz(2,z4(2,2,3,3),ep) =
     &  - Cv(cc001122 + C234,ep)

      Shat7zz(3,z4(2,2,3,3),ep) =
     &  - Cv(cc001122 + C234,ep)

      Shat7zz(1,z4(2,3,3,3),ep) =
     &  + Cv(cc001222 + C134,ep)
     &  - Cv(cc001222 + C234,ep)

      Shat7zz(2,z4(2,3,3,3),ep) =
     &  - Cv(cc001222 + C234,ep)

      Shat7zz(3,z4(2,3,3,3),ep) =
     &  - Cv(cc001222 + C234,ep)

      Shat7zz(1,z4(3,3,3,3),ep) =
     &  + Cv(cc002222 + C134,ep)
     &  - Cv(cc002222 + C234,ep)

      Shat7zz(2,z4(3,3,3,3),ep) =
     &  + Cv(cc002222 + C124,ep)
     &  - Cv(cc002222 + C234,ep)

      Shat7zz(3,z4(3,3,3,3),ep) =
     &  - Cv(cc002222 + C234,ep)

c--- definitions of Shat7zzzz have not been completed yet, so current
c--- implementation of D000000i will not be correct

c      Shat7zzzz(1,z2(1,1),ep) =
c     &  - Csum00001(ep)
c     &  - Csum00002(ep)
c     &  - Csum0000(ep)
c     &  - 4.D0*D0000001(P,K,L,m0,m1,m2,m3)

c      Shat7zzzz(2,z2(1,1),ep) =
c     &  - Csum00001(ep)
c     &  - Csum00002(ep)
c     &  - Csum0000(ep)
c     &  + Cv(cc000011 + C124,ep)

c      Shat7zzzz(3,z2(1,1),ep) =
c     &  - Csum00001(ep)
c     &  - Csum00002(ep)
c     &  - Csum0000(ep)
c     &  + Cv(cc000011 + C123,ep)

c      Shat7zzzz(1,z2(1,2),ep) =
c     &  + Csum00001(ep)
c     &  - 2.D0*D0000002(P,K,L,m0,m1,m2,m3)

c      Shat7zzzz(2,z2(1,2),ep) =
c     &  + Csum00001(ep)
c     &  - 2.D0*D0000001(P,K,L,m0,m1,m2,m3)

c      Shat7zzzz(3,z2(1,2),ep) =
c     &  + Csum00001(ep)
c     &  + Cv(cc000012 + C123,ep)

c      Shat7zzzz(1,z2(1,3),ep) =
c     &  + Csum00002(ep)
c     &  - 2.D0*D0000003(P,K,L,m0,m1,m2,m3)

c      Shat7zzzz(2,z2(1,3),ep) =
c     &  + Csum00002(ep)
c     &  + Cv(cc000012 + C124,ep)

c      Shat7zzzz(3,z2(1,3),ep) =
c     &  + Csum00002(ep)
c     &  - 2.D0*D0000001(P,K,L,m0,m1,m2,m3)

c      Shat7zzzz(1,z2(2,2),ep) =
c     &  + Cv(cc000011 + C134,ep)
c     &  - Cv(cc000011 + C234,ep)

c      Shat7zzzz(2,z2(2,2),ep) =
c     &  - 4.D0*D0000002(P,K,L,m0,m1,m2,m3)
c     &  - Cv(cc000011 + C234,ep)

c      Shat7zzzz(3,z2(2,2),ep) =
c     &  - Cv(cc000011 + C234,ep)
c     &  + Cv(cc000022 + C123,ep)

c      Shat7zzzz(1,z2(2,3),ep) =
c     &  + Cv(cc000012 + C134,ep)
c     &  - Cv(cc000012 + C234,ep)

c      Shat7zzzz(2,z2(2,3),ep) =
c     &  - 2.D0*D0000003(P,K,L,m0,m1,m2,m3)
c     &  - Cv(cc000012 + C234,ep)

c      Shat7zzzz(3,z2(2,3),ep) =
c     &  - 2.D0*D0000002(P,K,L,m0,m1,m2,m3)
c     &  - Cv(cc000012 + C234,ep)

c      Shat7zzzz(1,z2(3,3),ep) =
c     &  + Cv(cc000022 + C134,ep)
c     &  - Cv(cc000022 + C234,ep)

c      Shat7zzzz(2,z2(3,3),ep) =
c     &  + Cv(cc000022 + C124,ep)
c     &  - Cv(cc000022 + C234,ep)

c      Shat7zzzz(3,z2(3,3),ep) =
c     &  - 4.D0*D0000003(P,K,L,m0,m1,m2,m3)
c     &  - Cv(cc000022 + C234,ep)

      Shat7zzzzzz(1,ep) =
     &  + Cv(cc000000 + C134,ep)
     &  - Cv(cc000000 + C234,ep)

      Shat7zzzzzz(2,ep) =
     &  + Cv(cc000000 + C124,ep)
     &  - Cv(cc000000 + C234,ep)

      Shat7zzzzzz(3,ep) =
     &  + Cv(cc000000 + C123,ep)
     &  - Cv(cc000000 + C234,ep)

      enddo

      do ep=-2,0

c--- note: these are the triangle parts of the S00 functions that
c---       are defined above (and commented out), except that these
c---       are a factor of two smaller

c--- definitions of Czero5 in include file for now
      include 'Czero5.f'

      Czero4(z4(1,1,1,1),ep)=+Csum0(ep)+3d0*Csum1(ep)
     . +3d0*Csum2(ep)+3d0*Csum11(ep)+6d0*Csum12(ep)+3d0*Csum22(ep)
     . +Csum111(ep)+Csum222(ep)+3d0*Csum112(ep)+3d0*Csum122(ep)
      Czero4(z4(1,1,1,2),ep)=-(Csum1(ep)+2d0*Csum11(ep)
     . +2d0*Csum12(ep)+Csum111(ep)+2d0*Csum112(ep)+Csum122(ep))
      Czero4(z4(1,1,1,3),ep)=-(Csum2(ep)+2d0*Csum12(ep)
     . +2d0*Csum22(ep)+Csum112(ep)+2d0*Csum122(ep)+Csum222(ep))
      Czero4(z4(1,1,2,2),ep)=+Csum11(ep)+Csum111(ep)+Csum112(ep)
      Czero4(z4(1,1,2,3),ep)=+Csum12(ep)+Csum112(ep)+Csum122(ep)
      Czero4(z4(1,1,3,3),ep)=+Csum22(ep)+Csum122(ep)+Csum222(ep)
      Czero4(z4(1,2,2,2),ep)=-Csum111(ep)
      Czero4(z4(1,2,2,3),ep)=-Csum112(ep)
      Czero4(z4(1,2,3,3),ep)=-Csum122(ep)
      Czero4(z4(1,3,3,3),ep)=-Csum222(ep)
      Czero4(z4(2,2,2,2),ep)=+Cv(cc1111+C234,ep)
      Czero4(z4(2,2,2,3),ep)=+Cv(cc1112+C234,ep)
      Czero4(z4(2,2,3,3),ep)=+Cv(cc1122+C234,ep)
      Czero4(z4(2,3,3,3),ep)=+Cv(cc1222+C234,ep)
      Czero4(z4(3,3,3,3),ep)=+Cv(cc2222+C234,ep)

      Czero3(z3(1,1,1),ep)=
     &  - Csum0(ep)
     &  - 2.D0*Csum1(ep)
     &  - 2.D0*Csum2(ep)
     &  - Csum11(ep)
     &  - Csum22(ep)
     &  - 2.D0*Csum12(ep)

      Czero3(z3(2,2,2),ep)=
     &  + Cv(cc111 + C234,ep)
      Czero3(z3(3,3,3),ep)=
     &  + Cv(cc222 + C234,ep)

      Czero3(z3(1,1,2),ep)=
     &  + Csum1(ep)
     &  + Csum11(ep)
     &  + Csum12(ep)

      Czero3(z3(1,1,3),ep)=
     &  + Csum2(ep)
     &  + Csum22(ep)
     &  + Csum12(ep)

      Czero3(z3(1,2,2),ep)=
     &  - Csum11(ep)

      Czero3(z3(1,3,3),ep)=
     &  - Csum22(ep)

      Czero3(z3(2,2,3),ep)=
     &  + Cv(cc112 + C234,ep)

      Czero3(z3(2,3,3),ep)=
     &  + Cv(cc122 + C234,ep)

      Czero3(z3(1,2,3),ep)=
     &  - Csum12(ep)

      Czero2(z2(1,1),ep)=
     &  + Csum0(ep)
     &  + Csum1(ep)
     &  + Csum2(ep)

      Czero2(z2(2,2),ep)=
     &  + Cv(cc11 + C234,ep)

      Czero2(z2(3,3),ep)=
     &  + Cv(cc22 + C234,ep)

      Czero2(z2(1,2),ep)=
     &  - Csum1(ep)

      Czero2(z2(1,3),ep)=
     &  - Csum2(ep)

      Czero2(z2(2,3),ep)=
     &  + Cv(cc12 + C234,ep)

      Czero1(1,ep)=
     &  - Csum0(ep)

      Czero1(2,ep)=
     &  + Cv(cc1 + C234,ep)

      Czero1(3,ep)=
     &  + Cv(cc2 + C234,ep)

      Czero0(ep)=
     &  + Cv(cc0 + C234,ep)

      enddo

    
c--- find the largest f(k) for recursion relations
      kmax=1
      do k=2,np
      if (abs(f(k)) .ge. abs(f(kmax))) kmax=k
      enddo

      write(6,*) 'f(kmax) =',f(kmax)

C----Begin the iteration scheme

C--- Set all the Dv to zero
      do ep=-2,0
      do j=1,Ndd
      Dv(j+N0,ep)=czip
      enddo
      enddo

      do step=0,5
      if (step .eq. 5) goto 105
      if (step .eq. 4) goto 104
      if (step .eq. 3) goto 103
      if (step .eq. 2) goto 102
      if (step .eq. 1) goto 101
      if (step .eq. 0) goto 100

C----step 5: calculate D00iiii, Diiiii 
 105  continue

C--- Fixes D00iiii using extension of 5.69
C--- knowing Diiii, with corrections of order Gr(i,j)*Diiiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call runP_00iiii(i1,i2,i3,i4,m1,Gr,Czero4,N0)
      enddo
      enddo
      enddo
      enddo

C--- Fixes Diiiii using extension of 5.70
c--- knowing D00iiii, with a correction of order Gr(i,j)*Diiiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      do i5=i4,np
      call runP_iiiii(kmax,i1,i2,i3,i4,i5,f,Gr,Shat6,N0)
      enddo
      enddo
      enddo
      enddo
      enddo

C----step 4: calculate D00iii, Diiii, D0000i
 104  continue

C--- Fixes D00iii using extension of 5.69
C--- knowing Diii, with corrections of order Gr(i,j)*Diiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call runP_00iii(i1,i2,i3,m1,Gr,Czero3,N0)
      enddo
      enddo
      enddo

C--- Fixes Diiii using extension of 5.70
c--- knowing D00iii, with a correction of order Gr(i,j)*Diiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call runP_iiii(kmax,i1,i2,i3,i4,f,Gr,Shat5,N0)
      enddo
      enddo
      enddo
      enddo

c--- calculate S0000i (needs D00i) - required for D0000
      include 'S0000i_def.f'

C--- Fixes D0000i, with corrections of order Gr(i,j)*D00iii
      do i1=1,np
      call runP_0000i(i1,Gr,S0000i,N0)
      enddo
      
C----step 3: calculate D00ii, Diii, D0000 
 103  continue

C--- Fixes D00ii using 5.69
C--- knowing Dii, with corrections of order Gr(i,j)*Diiii
      do i1=1,np
      do i2=i1,np
      call runP_00ii(i1,i2,m1,Gr,Czero2,N0)
      enddo
      enddo

C--- Fixes Diii using 5.70
c--- knowing D00ii, with a correction of order Gr(i,j)*Diiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call runP_iii(kmax,i1,i2,i3,f,Gr,Shat4,N0)
      enddo
      enddo
      enddo

c--- calculate S0000 (needs D00) - required for D0000
      include 'S0000_def.f'

C--- Fixes D0000, with corrections of order Gr(i,j)*D00ii
      call runP_0000(Gr,S0000,N0)
     
C---- step 2: calculate D00i and Dii
 102  continue
C--- Fixes D00i according to 5.67 Denner-Dittmaier
C--- knowing Di, with correction of order Gr(i,j)*Diii
      do i1=1,np
      call runP_00i(i1,m1,Gr,Czero1,N0)
      enddo

C--- Fixes Dii using Eq. 5.68 Denner-Dittmaier
C--- knowing D00i with correction of order Gr(i,j)*Diii
      do i1=1,np
      do i2=i1,np
      call runP_ii(kmax,i1,i2,f,Gr,Shat3,N0)
      enddo
      enddo

C----step 1: calculate D00 and Di
 101  continue
C--- Fixes D00 according to 5.65 Denner-Dittmaier
C--- knowing D0 with corrections of order Gr(i,j)*Dii  
      call runP_00(m1,Gr,Czero0,N0)

C--- Fixes Di according to 5.66 Denner-Dittmaier
C--- knowing D00 with corrections of order Gr(i,j)*Dii
      do i1=1,np
      call runP_i(kmax,i1,f,Gr,Shat2,N0)
      enddo

C----step 0: calculate D0
 100  continue
C--- Fixes D0 according to 5.64 Denner-Dittmaier
C--- with corrections of order Gr(i,j)*Di
      call runP_0(kmax,f,Gr,Shat1,N0)

      enddo


c--- check the contents of box array    
c      write(6,*) 'D array'
c      do ip=1,Ndd
c      if (abs(Dsing(ip,p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4)).ne.0d0) then
c      write(6,'(i3,2f20.15)') ip,Dv(ip+N0,-1)
c     .    /Dsing(ip,p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4)
c      endif
c      enddo
c      pause
      

c   77 format(a3,i2,a5,3('(',e13.6,',',e13.6,') '))
    
      return
      end


