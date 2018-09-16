      subroutine Dfill_recur2(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N0,
     .                         exceptional)
      implicit none
C     Implements the calculation of the formfactors
C     for small Gram Determinant and small Y, as in DD Eq.5.54-5.61 etc
C     N0 is the offset in the common block

C--- Currently: calculates up to rank 3 with at least one recursion
c---            calculates ranks 4 and 5 with no recursion
c---            calculates metric tensor components of ranks 6 and 7

c--- JC: 11/22/2012 added an extra level of recursion. No additional
c---     identities are used, but the extra loop improves the
c---     numerical precision

      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'pvverbose.f'
      integer C234,C134,C124,C123,np,ep,N0,i,j,k,l,pvCcache,
     . i1,i2,i3,i4,i5,n,m,jx,step,kgt,lgt,ixt,jxt
      parameter(np=3)
      double precision p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,f(np),
     . Gtwiddle(np,np),Xtwiddle0(np),Gr(np,np),DetGr,Gtt(np,np,np,np),
     . Xtwiddle(0:np,0:np),Y(4,4),DetY,Xtmax
      double complex 
     . Shat3zz(np,-2:0),Shat4zz(np,z1max,-2:0),
     . Shat5zz(np,z2max,-2:0),Shat6zz(np,z3max,-2:0),
     . Shat5zzzz(np,-2:0),Shat6zzzz(np,z1max,-2:0),
     . Shat7zz(np,z4max,-2:0),Shat7zzzz(np,z2max,-2:0),
     . Shat7zzzzzz(np,-2:0),
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
      logical exceptional
      logical,save::first=.true.
!$omp threadprivate(first)

      if (first) then
        first=.false.
        call Array3dim
        call DArraysetup
      endif

      exceptional=.false.

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
      if (pvverbose) write(6,*) 'small Y: 3x3 DetGr = ',DetGr
      
c--- commented out: this is now tested in pvDfill
cC     Y(i,j)=mi^2+mj^2-(q_i-q_j)^2
cC     where q_1=0,  q_2=p1,  q_3=p_1+p_2, q_4=p_1+p_2+p_3;
      Y(1,1) = (2*m1)
      Y(1,2) = (m1 + m2 - p1)
      Y(2,1) = Y(1,2)
      Y(1,3) = (m1 + m3 - p1p2)
      Y(3,1) = Y(1,3)
      Y(1,4) = (m1 + m4 - p4)
      Y(4,1) = Y(1,4)
      Y(2,2) = (2*m2)
      Y(2,3) = (m2 + m3 - p2)
      Y(3,2) = Y(2,3)
      Y(2,4) = (m2 + m4 - p2p3)
      Y(4,2) = Y(2,4)
      Y(3,3) = (2*m3)
      Y(3,4) = (m3 + m4 - p3)
      Y(4,3) = Y(3,4)
      Y(4,4) = (2*m4)
      call determinant(4,4,Y,DetY)
      if (pvverbose) write(6,*) 'DetY = ',DetY

c      Ysing=pvGramsing(Y,4)

      Gtwiddle(1,1)=Gr(2,2)*Gr(3,3)-Gr(2,3)**2
      Gtwiddle(2,2)=Gr(1,1)*Gr(3,3)-Gr(1,3)**2
      Gtwiddle(3,3)=Gr(1,1)*Gr(2,2)-Gr(1,2)**2
      Gtwiddle(1,2)=-(Gr(1,2)*Gr(3,3)-Gr(1,3)*Gr(2,3))
      Gtwiddle(2,1)=Gtwiddle(1,2)
      Gtwiddle(2,3)=-(Gr(1,1)*Gr(2,3)-Gr(1,2)*Gr(1,3))
      Gtwiddle(3,2)=Gtwiddle(2,3)
      Gtwiddle(1,3)=Gr(1,2)*Gr(2,3)-Gr(1,3)*Gr(2,2)
      Gtwiddle(3,1)=Gtwiddle(1,3)

      do j=1,3
      Xtwiddle0(j)=
     . -Gtwiddle(j,1)*f(1)-Gtwiddle(j,2)*f(2)-Gtwiddle(j,3)*f(3)
      Xtwiddle(0,j)=
     . -Gtwiddle(j,1)*f(1)-Gtwiddle(j,2)*f(2)-Gtwiddle(j,3)*f(3)
      Xtwiddle(j,0)=Xtwiddle(0,j)
      enddo
      Xtwiddle(0,0)=DetGr
c      if (abs(Xtwiddle0(2)) .gt. abs(Xtwiddle0(jx))) jx=2
c      if (abs(Xtwiddle0(3)) .gt. abs(Xtwiddle0(jx))) jx=3

C----setup Gtt
      Gtt(1,2,1,2)=-Gr(3,3)
      Gtt(1,2,1,3)=+Gr(3,2)
      Gtt(1,2,2,3)=-Gr(3,1)
      Gtt(1,3,1,2)=+Gr(2,3)
      Gtt(1,3,1,3)=-Gr(2,2)
      Gtt(1,3,2,3)=+Gr(2,1)
      Gtt(2,3,1,2)=-Gr(1,3)
      Gtt(2,3,1,3)=+Gr(1,2)
      Gtt(2,3,2,3)=-Gr(1,1)

      do i=1,3
      do k=i,3
      do j=1,3
      do l=j,3
      if ((i .eq. k) .or. (j .eq. l)) then
      Gtt(i,k,j,l)=0d0
      endif
      Gtt(k,i,j,l)=-Gtt(i,k,j,l)
      Gtt(i,k,l,j)=-Gtt(i,k,j,l)
      Gtt(k,i,l,j)=+Gtt(i,k,j,l)
      enddo
      enddo
      enddo
      enddo

c--- setup Xtwiddle
      do i=1,3
      do j=1,3
      Xtwiddle(i,j)=2d0*m1*Gtwiddle(i,j)
      do n=1,3
      do m=1,3
      Xtwiddle(i,j)=Xtwiddle(i,j)+Gtt(i,n,j,m)*f(n)*f(m)
      enddo
      enddo
      enddo
      enddo

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

      jx=1
      do j=2,np
      if (abs(Xtwiddle0(j)) .ge. abs(Xtwiddle0(jx))) jx=j
      enddo

      kgt=1
      lgt=1
      do k=1,np
      do l=k,np
      if (abs(Gtwiddle(k,l)) .ge. abs(Gtwiddle(kgt,lgt))) then
      kgt=k
      lgt=l
      endif
      enddo
      enddo
      
      ixt=1
      jxt=1
      do i=1,np
      do j=i,np
      if (abs(Xtwiddle(i,j)) .ge. abs(Xtwiddle(ixt,jxt))) then
      ixt=i
      jxt=j
      endif
      enddo
      enddo
      
c      do k=1,np
c      do l=k,np
c      write(6,*) k,l,Gtwiddle(k,l),Gtwiddle(kgt,lgt)
c      enddo
c      enddo
c      pause

c      write(6,*) '   Xtwiddle(0,0)', Xtwiddle(0,0)
c      do j=1,np
c      write(6,*) 'j, Xtwiddle(0,j)',j, Xtwiddle(0,j)
c      enddo
      
c      do j=1,np
c      do k=1,np
c      write(6,*) 'j, k, Xtwiddle(j,k)',j,k, Xtwiddle(j,k)
c      enddo
c      enddo

      Xtmax=abs(Gr(1,1))
      do j=1,np
      do k=1,np
      Xtmax=max(abs(Xtwiddle(j,k)),Xtmax)
      enddo
      enddo
c      write(6,*) 'Xtmax=',Xtmax

c--- check for exceptional case where none of the Xtwiddle(i,j) elements
c---  are large compared to DetGr (Xtwiddle(0,0)) or Xtwiddle(0,j)
c---  [see note at end of Sec.5.5 in DD].
      if ( (Xtmax/abs(Xtwiddle(0,0)) .lt. 1d1)
     . .or.(Xtmax/abs(Xtwiddle0(jx)) .lt. 1d1)) then
        if (pvverbose) then
          write(6,*) 'EXCEPTIONAL CASE'
        write(6,*) 'Maximum Xtwiddle(i,j) = ',Xtmax
        write(6,*) '        Xtwiddle(0,0) = ',Xtwiddle(0,0)
        write(6,*) 'maximum Xtwiddle(0,j) = ',Xtwiddle(0,jx)
      endif
      exceptional=.true.
      return
      endif
    
C----Begin the iteration scheme


C--- Set all the Dv to zero
      do ep=-2,0
      do j=1,Ndd
      Dv(j+N0,ep)=czip
      enddo
      enddo


      do step=0,3
      if (step .eq. 3) goto 103
      if (step .eq. 2) goto 102
      if (step .eq. 1) goto 101
      if (step .eq. 0) goto 100

C--- step 3
 103  continue

C--- step 2: calculate D00iiii, D00iiiii, Diiii, Diiiii,
c---                   D0000ii, D0000iii, D000000,D000000i
 102  continue

C--- a) Calculate D00iiii
C---    Small terms of order Xtwiddle(0,k)*Diiiii,Xtwiddle(0,0)*Diiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00llll(kgt,lgt,Xtwiddle,Gtwiddle,Shat6,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate D00llli, requires D00llll
C---    Small terms of order Xtwiddle(0,k)*Diiiii,Xtwiddle(0,0)*Diiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00llli(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat6,N0)
      endif
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then
C---    Calculate D00lli1i2, requires D00llli1
C---    Small terms of order Xtwiddle(0,k)*Diiiii,Xtwiddle(0,0)*Diiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00lli1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat6,N0)
      endif
      enddo
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt) .and. (i3 .ne. lgt)) then
C---    Calculate D00li1i2i3, requires D00lli1i2
C---    Small terms of order Xtwiddle(0,k)*Diiiii,Xtwiddle(0,0)*Diiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00li1i2i3(kgt,lgt,i1,i2,i3,Xtwiddle,Gtwiddle,Shat6,N0)
      endif
      enddo
      enddo
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      if (   (i1 .ne. lgt) .and. (i2 .ne. lgt)
     . .and. (i3 .ne. lgt) .and. (i4 .ne. lgt)) then
C---    Calculate D00i1i2i3i4, requires D00li1i2i3
C---    Small terms of order Xtwiddle(0,k)*Diiiii,Xtwiddle(0,0)*Diiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00i1i2i3i4(kgt,lgt,i1,i2,i3,i4,
     .                     Xtwiddle,Gtwiddle,Shat6,N0)
      endif
      enddo
      enddo
      enddo
      enddo

C--- b) Calculate D00iiiii
C---    Small terms of order Xtwiddle(0,k)*Diiiiii,Xtwiddle(0,0)*Diiiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00lllll(kgt,lgt,Xtwiddle,Gtwiddle,Shat7,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate D00lllli, requires D00lllll
C---    Small terms of order Xtwiddle(0,k)*Diiiiii,Xtwiddle(0,0)*Diiiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00lllli(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then
C---    Calculate D00llli1i2, requires D00lllli1
C---    Small terms of order Xtwiddle(0,k)*Diiiiii,Xtwiddle(0,0)*Diiiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00llli1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt) .and. (i3 .ne. lgt)) then
C---    Calculate D00lli1i2i3, requires D00llli1i2
C---    Small terms of order Xtwiddle(0,k)*Diiiiii,Xtwiddle(0,0)*Diiiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00lli1i2i3(kgt,lgt,i1,i2,i3,Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo
      enddo
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      if (   (i1 .ne. lgt) .and. (i2 .ne. lgt)
     . .and. (i3 .ne. lgt) .and. (i4 .ne. lgt)) then
C---    Calculate D00li1i2i3i4, requires D00lli1i2i3
C---    Small terms of order Xtwiddle(0,k)*Diiiiii,Xtwiddle(0,0)*Diiiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00li1i2i3i4(kgt,lgt,i1,i2,i3,i4,
     .                      Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo
      enddo
      enddo
      enddo
      
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      do i5=i4,np
      if (   (i1 .ne. lgt) .and. (i2 .ne. lgt)
     . .and. (i3 .ne. lgt) .and. (i4 .ne. lgt) .and. (i5 .ne. lgt)) then
C---    Calculate D00i1i2i3i4i5, requires D00li1i2i3i4
C---    Small terms of order Xtwiddle(0,k)*Diiiiii,Xtwiddle(0,0)*Diiiiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00i1i2i3i4i5(kgt,lgt,i1,i2,i3,i4,i5,
     .                       Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo
      enddo
      enddo
      enddo
      enddo

C--- c) Calculate Diiiii, requires D00iiii,D00iiiii
C---    Small terms of order Xtwiddle(0,j)*Diiiiii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      do i5=i4,np
      call runY_i1i2i3i4i5(ixt,jxt,i1,i2,i3,i4,i5,
     .                     f,Xtwiddle,Gtt,Gtwiddle,Shat6,Czero5,N0)
      enddo
      enddo
      enddo
      enddo
      enddo

C--- d) Calculate Diiii, requires D00iii,D00iiii
C---    Small terms of order Xtwiddle(0,j)*Diiiii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call runY_i1i2i3i4(ixt,jxt,i1,i2,i3,i4,
     . f,Xtwiddle,Gtt,Gtwiddle,Shat5,Czero4,N0)
      enddo
      enddo
      enddo
      enddo

C--- e) Calculate D0000ii
C---    Small terms of order Xtwiddle(0,k)*Dzziii,Xtwiddle(0,0)*Dzziiii
C---    Denominator Gtwiddle(k,l)
      call runY_0000ll(kgt,lgt,Xtwiddle,Gtwiddle,Shat6zz,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate D0000li, requires D0000ll
C---    Small terms of order Xtwiddle(0,k)*D00iii,Xtwiddle(0,0)*D00iiii
C---    Denominator Gtwiddle(k,l)
      call runY_0000li(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat6zz,N0)
      endif 
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1.ne.lgt) .and. (i2 .ne. lgt)) then  
C---    Calculate D0000i1i2, requires D0000li1,D0000li2
C---    Small terms of order Xtwiddle(0,k)*D00iii,Xtwiddle(0,0)*D00iiii
C---    Denominator Gtwiddle(k,l)
      call runY_0000i1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat6zz,N0)
      endif
      enddo
      enddo
      
C--- f) Calculate D0000iii
C---    Small terms of order Xtwiddle(0,k)*Dzziiii,Xtwiddle(0,0)*Dzziiiii
C---    Denominator Gtwiddle(k,l)
      call runY_0000lll(kgt,lgt,Xtwiddle,Gtwiddle,Shat7zz,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate D0000lli, requires D0000lll
C---    Small terms of order Xtwiddle(0,k)*D00iiii,Xtwiddle(0,0)*D00iiiii
C---    Denominator Gtwiddle(k,l)
      call runY_0000lli(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat7zz,N0)
      endif
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then
C---    Calculate D0000li1i2, requires D0000lli1
C---    Small terms of order Xtwiddle(0,k)*D00iiii,Xtwiddle(0,0)*D00iiiii
C---    Denominator Gtwiddle(k,l)
      call runY_0000li1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat7zz,N0)
      endif
      enddo
      enddo

C---    Calculate D0000i1i2i3, requires D0000li1i2,D0000li2i3,D0000li3i1
C---    Small terms of order Xtwiddle(0,k)*D00iiii,Xtwiddle(0,0)*D00iiiii
C---    Denominator Gtwiddle(k,l)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      if ((i1 .ne. lgt) .and.(i2 .ne. lgt) .and.(i3 .ne. lgt)) then
      call runY_0000i1i2i3(kgt,lgt,i1,i2,i3,
     .                     Xtwiddle,Gtwiddle,Shat7zz,N0)
      endif
      enddo
      enddo
      enddo

C--- g) Calculate D000000
C---    Small terms of order Xtwiddle(0,k)*D0000i,Xtwiddle(0,0)*D0000ii
C---    Denominator Gtwiddle(k,l)
      call runY_000000(kgt,lgt,Xtwiddle,Gtwiddle,Shat6zzzz,N0)

C--- h) Calculate D000000i
C---    Small terms of order Xtwiddle(0,k)*D0000ii,Xtwiddle(0,0)*D0000iii
C---    Denominator Gtwiddle(k,l)
      call runY_000000l(kgt,lgt,Xtwiddle,Gtwiddle,Shat7zzzz,N0)

      do i1=1,np
      if (i1 .ne. lgt)
     . call runY_000000i(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat7zzzz,N0)
      enddo

C--- step 1: calculate D00ii, D00iii, Dii, Diii, D0000, D0000i
 101  continue

C--- a) Calculate D00ii
C---    Small terms of order Xtwiddle(0,k)*Diii,Xtwiddle(0,0)*Diiii
C---    Denominator Gtwiddle(k,l)
      call runY_00ll(kgt,lgt,Xtwiddle,Gtwiddle,Shat4,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate D00li, requires D00ll
C---    Small terms of order Xtwiddle(0,k)*Diii,Xtwiddle(0,0)*Diiii
C---    Denominator Gtwiddle(k,l)
      call runY_00li(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat4,N0)
      endif 
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then  
C---    Calculate D00i1i2, requires D00li1,D00li2
C---    Small terms of order Xtwiddle(0,k)*Diii,Xtwiddle(0,0)*Diiii
C---    Denominator Gtwiddle(k,l)
      call runY_00i1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat4,N0)
      endif
      enddo
      enddo
      
c--- b) Calculate D00iii
C---    Small terms of order Xtwiddle(0,k)*Diiii,Xtwiddle(0,0)*Diiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00lll(kgt,lgt,Xtwiddle,Gtwiddle,Shat5,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate D00lli, requires D00lll
C---    Small terms of order Xtwiddle(0,k)*Diiii,Xtwiddle(0,0)*Diiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00lli(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat5,N0)
      endif
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then
C---    Calculate D00li1i2, requires D00lli1
C---    Small terms of order Xtwiddle(0,k)*Diiii,Xtwiddle(0,0)*Diiiii
C---    Denominator Gtwiddle(k,l)
      call runY_00li1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat5,N0)
      endif
      enddo
      enddo

C---    Calculate D00i1i2i3, requires D00li1i2,D00li2i3,D00li3i1
C---    Small terms of order Xtwiddle(0,k)*Diiii,Xtwiddle(0,0)*Diiiii
C---    Denominator Gtwiddle(k,l)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      if ((i1 .ne. lgt) .and.(i2 .ne. lgt) .and.(i3 .ne. lgt)) then
      call runY_00i1i2i3(kgt,lgt,i1,i2,i3,Xtwiddle,Gtwiddle,Shat5,N0)
      endif
      enddo
      enddo
      enddo

C--- c) Calculate Diii, requires D00ii,D00iii
C---    Small terms of order Xtwiddle(0,j)*Diiii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call runY_i1i2i3(ixt,jxt,i1,i2,i3,f,Xtwiddle,Gtt,Gtwiddle,Shat4,
     . Czero3,N0)
      enddo
      enddo
      enddo

C--- d) Calculate Dii, requires D00i,D00ii
C---    Small terms of order Xtwiddle(0,j)*Diii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      do i2=i1,np
      call runY_i1i2(ixt,jxt,i1,i2,f,Xtwiddle,Gtt,Gtwiddle,Shat3,
     . Czero2,N0)
      enddo
      enddo

C--- e) Calculate D0000
C---    Small terms of order Xtwiddle(0,k)*D00i,Xtwiddle(0,0)*D00ii
C---    Denominator Gtwiddle(k,l)
      call runY_0000(kgt,lgt,Xtwiddle,Gtwiddle,Shat4zz,N0)

C--- f) Calculate D0000l
C---    Small terms of order Xtwiddle(0,k)*D00ii,Xtwiddle(0,0)*D00iii
C---    Denominator Gtwiddle(k,l)
      call runY_0000l(kgt,lgt,Xtwiddle,Gtwiddle,Shat5zz,N0)

      do i1=1,np
      if (i1 .ne. lgt)
     . call runY_0000i(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat5zz,N0)
      enddo


C--- step 0: calculate D00,D00i,D0 and Di
 100  continue
C--- a) Calculate D00
C---    Small terms of order Xtwiddle(0,k)*Di,Xtwiddle(0,0)*Dii
C---    Denominator Gtwiddle(kgt,lgt)
      call runY_00(kgt,lgt,Xtwiddle,Gtwiddle,Shat2,N0)

C--- b) Calculate D00l, D00i
C---    Small terms of order Xtwiddle(0,k)*Dii,Xtwiddle(0,0)*Diii
C---    Denominator Gtwiddle(k,l)
      call runY_00l(kgt,lgt,Xtwiddle,Gtwiddle,Shat3,N0)
C---    Calculate D00i1, requires D00l
C---    Small terms of order Xtwiddle(0,k)*Dli1,Xtwiddle(0,0)*Dkli1
C---    Denominator Gtwiddle(k,l)
      do i1=1,np
      if (i1 .ne. lgt) then
        call runY_00i(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat3,N0)
      endif
      enddo

C--- c) Calculate Di, requires D00i
C---    Small terms of order Xtwiddle(0,j)*Dii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      call runY_i(ixt,jxt,i1,f,Xtwiddle,Gtt,Gtwiddle,Shat2,Czero1,N0)
      enddo

C--- d) Calculates D0
C---    Requires D00, small terms of order Xtwiddle(0,j)*Di
C---    Denominator Xtwiddle(i,j)
      call runY_0(ixt,jxt,f,Xtwiddle,Gtwiddle,Gtt,Shat1,Czero0,N0)

c--- check the contents of box array    
c      write(6,*) 'recur2: D array'
c      do ip=1,24
c        write(6,'(i3,2e20.12)') ip,Dv(ip+N0,0)
c      enddo
c      pause

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

c      write(6,*) 'Leaving small Y recursion'
c      pause

c   77 format(a3,i2,a5,3('(',e13.6,',',e13.6,') '))
    
      return
      end


