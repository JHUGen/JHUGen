      subroutine Dfill_recur(p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,N0)
      implicit none
C     Implements the calculation of the formfactors
C     for small Gram Determinant as in DD Eq.5.41-5.48 etc
C     N0 is the offset in the common block

C--- Currently: calculates up to rank 4 with at least one recursion
c---            calculates rank 5 with no recursion
c---            calculates metric tensor components of rank 6

      include 'TRconstants.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'Darraydef.f'
      include 'pvverbose.f'
      integer C234,C134,C124,C123,np,ep,N0,i,j,k,l,pvCcache,
     . i1,i2,i3,i4,i5,step,jmax,kmax,lmax
      parameter(np=3)
      double precision p1,p2,p3,p4,p1p2,p2p3,m1,m2,m3,m4,f(np),
     . Gtwiddle(np,np),Xtwiddle0(np),Gr(np,np),DetGr,Gtt(np,np,np,np)
      double complex S00(-2:0),S0000(-2:0),S000000(-2:0),
     . S0000i(np,-2:0),S0000ii(z2max,-2:0),
     . S00i(np,-2:0),S00ii(z2max,-2:0),
     . S00iii(z3max,-2:0),S00iiii(z4max,-2:0),
     . Shat3zz(np,-2:0),Shat4zz(np,z1max,-2:0),
     . Shat5zz(np,z2max,-2:0),Shat6zz(np,z3max,-2:0),
     . Shat5zzzz(np,-2:0),Shat6zzzz(np,z1max,-2:0),
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
     . csum0012(-2:0),csum0022(-2:0),csum0000(-2:0)
      logical,save:: first=.true.
!$omp threadprivate(first)


      if (first) then
        first=.false.
        call Array3dim
        call DArraysetup
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
      
      if (pvverbose) write(6,*) 'small G: 3x3 DetGr = ',DetGr
      
      Gtwiddle(1,1)=Gr(2,2)*Gr(3,3)-Gr(2,3)**2
      Gtwiddle(2,2)=Gr(1,1)*Gr(3,3)-Gr(1,3)**2
      Gtwiddle(3,3)=Gr(1,1)*Gr(2,2)-Gr(1,2)**2
      Gtwiddle(1,2)=-(Gr(1,2)*Gr(3,3)-Gr(1,3)*Gr(2,3))
      Gtwiddle(2,1)=Gtwiddle(1,2)
      Gtwiddle(2,3)=-(Gr(1,1)*Gr(2,3)-Gr(1,2)*Gr(1,3))
      Gtwiddle(3,2)=Gtwiddle(2,3)
      Gtwiddle(1,3)=Gr(1,2)*Gr(2,3)-Gr(1,3)*Gr(2,2)
      Gtwiddle(3,1)=Gtwiddle(1,3)

      do j=1,np
      Xtwiddle0(j)=
     . -Gtwiddle(j,1)*f(1)-Gtwiddle(j,2)*f(2)-Gtwiddle(j,3)*f(3)
      enddo


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

      do i=1,np
      do k=i,np
      do j=1,np
      do l=j,np
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

      enddo

      do ep=-2,0
      include 'Shat.f'
      enddo

      do ep=-2,0
c--- These are now calculated in the recursion, which is required when
c--- the first propagator has a non-zero mass
c      include 'S00.f'

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

      enddo


      jmax=1
      do j=2,np
      if (abs(Xtwiddle0(j)) .ge. abs(Xtwiddle0(jmax))) jmax=j
      enddo
      
      kmax=1
      lmax=1
      do k=1,np
      do l=k,np
      if (abs(Gtwiddle(k,l)) .ge. abs(Gtwiddle(kmax,lmax))) then
      kmax=k
      lmax=l
      endif
      enddo
      enddo
      
C----Begin the iteration scheme

C Set all the Dv to zero
      do ep=-2,0
      do j=1,Ndd
      Dv(j+N0,ep)=czip
      enddo
      enddo


      do step=0,5
c      if (step .eq. 6) goto 106
      if (step .eq. 5) goto 105
      if (step .eq. 4) goto 104
      if (step .eq. 3) goto 103
      if (step .eq. 2) goto 102
      if (step .eq. 1) goto 101
      if (step .eq. 0) goto 100


C----step 6 
c 106  continue

C---step 5
 105  continue
C--- Fixes D000000 according to extension of Denner-Dittmaier 
C--- knowing D0000 with correction of order Delta*D0000ii

c--- calculate S000000 (needs D0000)
      do ep=-2,0
      S000000(ep)=2d0*Cv(cc0000+C234,ep)+2d0*m1*Dv(dd0000+N0,ep)
      enddo

      call run_000000(kmax,lmax,DetGr,f,Gtwiddle,Gtt,
     . Shat5zzzz,Shat6zzzz,S000000,N0)
C--- Fixes D0000ii according to extension of Denner-Dittmaier
C--- knowing D00ii,D0000i,D000000 with correction of order Delta*D00iiii

c--- calculate S0000ii (needs D00ii)
      do ep=-2,0
      S0000ii(z2(1,1),ep)=+2d0*(Csum00(ep)+Csum001(ep)+Csum002(ep))
     .                    +2d0*m1*Dv(dd0011+N0,ep)
      S0000ii(z2(1,2),ep)=-2d0*Csum001(ep)
     .                    +2d0*m1*Dv(dd0012+N0,ep)
      S0000ii(z2(1,3),ep)=-2d0*Csum002(ep)
     .                    +2d0*m1*Dv(dd0013+N0,ep)
      S0000ii(z2(2,2),ep)=+2d0*Cv(cc0011+C234,ep)
     .                    +2d0*m1*Dv(dd0022+N0,ep)
      S0000ii(z2(2,3),ep)=+2d0*Cv(cc0012+C234,ep)
     .                    +2d0*m1*Dv(dd0023+N0,ep)
      S0000ii(z2(3,3),ep)=+2d0*Cv(cc0022+C234,ep)
     .                    +2d0*m1*Dv(dd0033+N0,ep)
      enddo

      do i1=1,np
      do i2=i1,np
      call run_0000ii(kmax,lmax,i1,i2,DetGr,f,Gtwiddle,Gtt, 
     . Shat5zz,Shat6zzzz,S0000ii,Shat6zz,N0) 
      enddo
      enddo
C--- Fixes D00iiii according to extension of Denner-Dittmaier
C--- knowing Diiii,D00iii,D0000ii, corrections of order Delta*Diiiiii

c--- calculate S00iiii (needs Diiii)
      do ep=-2,0
      S00iiii(z4(1,1,1,1),ep)=+2d0*(+Csum0(ep)+3d0*Csum1(ep)
     . +3d0*Csum2(ep)+3d0*Csum11(ep)+6d0*Csum12(ep)+3d0*Csum22(ep)
     . +Csum111(ep)+Csum222(ep)+3d0*Csum112(ep)+3d0*Csum122(ep))
     .                        +2d0*m1*Dv(dd1111+N0,ep)
      S00iiii(z4(1,1,1,2),ep)=-2d0*(Csum1(ep)+2d0*Csum11(ep)
     . +2d0*Csum12(ep)+Csum111(ep)+2d0*Csum112(ep)+Csum122(ep))
     .                        +2d0*m1*Dv(dd1112+N0,ep)
      S00iiii(z4(1,1,1,3),ep)=-2d0*(Csum2(ep)+2d0*Csum12(ep)
     . +2d0*Csum22(ep)+Csum112(ep)+2d0*Csum122(ep)+Csum222(ep))
     .                        +2d0*m1*Dv(dd1113+N0,ep)
      S00iiii(z4(1,1,2,2),ep)=+2d0*(Csum11(ep)+Csum111(ep)+Csum112(ep))
     .                        +2d0*m1*Dv(dd1122+N0,ep)
      S00iiii(z4(1,1,2,3),ep)=+2d0*(Csum12(ep)+Csum112(ep)+Csum122(ep))
     .                        +2d0*m1*Dv(dd1123+N0,ep)
      S00iiii(z4(1,1,3,3),ep)=+2d0*(Csum22(ep)+Csum122(ep)+Csum222(ep))
     .                        +2d0*m1*Dv(dd1133+N0,ep)
      S00iiii(z4(1,2,2,2),ep)=-2d0*Csum111(ep)
     .                        +2d0*m1*Dv(dd1222+N0,ep)
      S00iiii(z4(1,2,2,3),ep)=-2d0*Csum112(ep)
     .                        +2d0*m1*Dv(dd1223+N0,ep)
      S00iiii(z4(1,2,3,3),ep)=-2d0*Csum122(ep)
     .                        +2d0*m1*Dv(dd1233+N0,ep)
      S00iiii(z4(1,3,3,3),ep)=-2d0*Csum222(ep)
     .                        +2d0*m1*Dv(dd1333+N0,ep)
      S00iiii(z4(2,2,2,2),ep)=+2d0*Cv(cc1111+C234,ep)
     .                        +2d0*m1*Dv(dd2222+N0,ep)
      S00iiii(z4(2,2,2,3),ep)=+2d0*Cv(cc1112+C234,ep)
     .                        +2d0*m1*Dv(dd2223+N0,ep)
      S00iiii(z4(2,2,3,3),ep)=+2d0*Cv(cc1122+C234,ep)
     .                        +2d0*m1*Dv(dd2233+N0,ep)
      S00iiii(z4(2,3,3,3),ep)=+2d0*Cv(cc1222+C234,ep)
     .                        +2d0*m1*Dv(dd2333+N0,ep)
      S00iiii(z4(3,3,3,3),ep)=+2d0*Cv(cc2222+C234,ep)
     .                        +2d0*m1*Dv(dd3333+N0,ep)
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call run_00iiii(kmax,lmax,i1,i2,i3,i4,DetGr,f,Gtwiddle,Gtt,
     . Shat5,Shat6,S00iiii,Shat6zz,N0)
      enddo
      enddo
      enddo
      enddo
C--- Fixes Diiiii according to extension of Denner-Dittmaier
c--- knowing D00iiii with a correction of order Delta*Diiiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      do i5=i4,np
      call run_iiiii(jmax,i1,i2,i3,i4,i5,DetGr,Xtwiddle0,Gtwiddle,
     . Shat6,N0)
      enddo
      enddo
      enddo
      enddo
      enddo

C----step 4 Calculate Diiii
 104  continue
C--- Fixes D0000i according to extension of Denner-Dittmaier
C--- knowing D00i,D0000 with corrections of order Delta*D00iii
 
c--- calculate S0000i (needs D00i)
      do ep=-2,0
      S0000i(1,ep)=-2d0*Csum00(ep)+2d0*m1*Dv(dd001+N0,ep)
      S0000i(2,ep)=+2d0*Cv(cc001+C234,ep)+2d0*m1*Dv(dd002+N0,ep)
      S0000i(3,ep)=+2d0*Cv(cc002+C234,ep)+2d0*m1*Dv(dd003+N0,ep)
      enddo

      do i1=1,np
      call run_0000i(kmax,lmax,i1,DetGr,f,Gtwiddle,Gtt, 
     . Shat4zz,Shat5zzzz,S0000i,Shat5zz,N0) 
      enddo
C--- Fixes D00iii according to extension of Denner-Dittmaier
c--- knowing Diii,D00ii,D0000i with corrections of order Delta*Diiiii

c--- calculate S00iii (needs Diii)
      do ep=-2,0
      S00iii(z3(1,1,1),ep)=-2d0*(Csum0(ep)+2d0*Csum1(ep)+2d0*Csum2(ep)
     .                          +Csum11(ep)+2d0*Csum12(ep)+Csum22(ep))
     .                     +2d0*m1*Dv(dd111+N0,ep)
      S00iii(z3(1,1,2),ep)=+2d0*(Csum1(ep)+Csum11(ep)+Csum12(ep))
     .                     +2d0*m1*Dv(dd112+N0,ep)
      S00iii(z3(1,1,3),ep)=+2d0*(Csum2(ep)+Csum12(ep)+Csum22(ep))
     .                     +2d0*m1*Dv(dd113+N0,ep)
      S00iii(z3(1,2,2),ep)=-2d0*Csum11(ep)
     .                     +2d0*m1*Dv(dd122+N0,ep)
      S00iii(z3(1,2,3),ep)=-2d0*Csum12(ep)
     .                     +2d0*m1*Dv(dd123+N0,ep)
      S00iii(z3(1,3,3),ep)=-2d0*Csum22(ep)
     .                     +2d0*m1*Dv(dd133+N0,ep)
      S00iii(z3(2,2,2),ep)=+2d0*Cv(cc111+C234,ep)
     .                     +2d0*m1*Dv(dd222+N0,ep)
      S00iii(z3(2,2,3),ep)=+2d0*Cv(cc112+C234,ep)
     .                     +2d0*m1*Dv(dd223+N0,ep)
      S00iii(z3(2,3,3),ep)=+2d0*Cv(cc122+C234,ep)
     .                     +2d0*m1*Dv(dd233+N0,ep)
      S00iii(z3(3,3,3),ep)=+2d0*Cv(cc222+C234,ep)
     .                     +2d0*m1*Dv(dd333+N0,ep)
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call run_00iii(kmax,lmax,i1,i2,i3,DetGr,f,Gtwiddle,Gtt,
     . Shat4,Shat5,S00iii,Shat5zz,N0)
      enddo
      enddo
      enddo
C--- Fixes Diiii according to extension of Denner-Dittmaier
C--- knowing D00iii and dropping terms of order Delta*Diiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call run_iiii(jmax,i1,i2,i3,i4,DetGr,Xtwiddle0,Gtwiddle,Shat5,N0)
      enddo
      enddo
      enddo
      enddo

C----step 3 
 103  continue
C--- Fixes D0000 using 5.46
C--- knowing D00 with correction of order Delta*D00ii

c--- calculate S0000 (needs D00)
      do ep=-2,0
      S0000(ep)=2d0*Cv(cc00+C234,ep)+2d0*m1*Dv(dd00+N0,ep)
      enddo

      call run_0000(kmax,lmax,DetGr,f,Gtwiddle,Gtt,Shat3zz,Shat4zz,
     . S0000,N0)
C--- Fixes D00ii using 5.47
C--- knowing Dii,D00i,D0000 with corrections of order Delta*Diiii

c--- calculate S00ii (needs Dii)
      do ep=-2,0
      S00ii(z2(1,1),ep)=+2d0*(Csum0(ep)+Csum1(ep)+Csum2(ep))
     .                  +2d0*m1*Dv(dd11+N0,ep)
      S00ii(z2(1,2),ep)=-2d0*Csum1(ep)+2d0*m1*Dv(dd12+N0,ep)
      S00ii(z2(1,3),ep)=-2d0*Csum2(ep)+2d0*m1*Dv(dd13+N0,ep)
      S00ii(z2(2,2),ep)=+2d0*Cv(cc11+C234,ep)+2d0*m1*Dv(dd22+N0,ep)
      S00ii(z2(2,3),ep)=+2d0*Cv(cc12+C234,ep)+2d0*m1*Dv(dd23+N0,ep)
      S00ii(z2(3,3),ep)=+2d0*Cv(cc22+C234,ep)+2d0*m1*Dv(dd33+N0,ep)
      enddo

      do i1=1,np
      do i2=i1,np
      call run_00ii(kmax,lmax,i1,i2,DetGr,f,Gtwiddle,Gtt,Shat3,Shat4,
     . S00ii,Shat4zz,N0)
      enddo
      enddo
C--- Fixes Diii using 5.48
c--- knowing D00ii with a correction of order Delta*Diiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call run_iii(jmax,i1,i2,i3,DetGr,Xtwiddle0,Gtwiddle,Shat4,N0)
      enddo
      enddo
      enddo

C---- step 2 calculate D00i and Dii
 102  continue
C--- Fixes D00i according to 5.44 Denner-Dittmaier
C--- knowing Di, D00 with correction of order Delta*Diii

c--- calculate S00i (needs Di)
      do ep=-2,0
      S00i(1,ep)=-2d0*Csum0(ep)+2d0*m1*Dv(dd1+N0,ep) 
      S00i(2,ep)=+2d0*Cv(cc1+C234,ep)+2d0*m1*Dv(dd2+N0,ep)
      S00i(3,ep)=+2d0*Cv(cc2+C234,ep)+2d0*m1*Dv(dd3+N0,ep)
      enddo

      do i3=1,np
      call run_00i(kmax,lmax,i3,DetGr,f,Gtwiddle,Gtt,Shat2,Shat3,
     . Shat3zz,S00i,N0)
      enddo
C--- Fixes Dii using Eq. 5.45 Denner-Dittmaier
C--- knowing D00i with correction of order Delta*Diii
      do i2=1,np
      do i3=i2,np
      call run_ii(jmax,i2,i3,DetGr,Xtwiddle0,Gtwiddle,Shat3,N0)
      enddo
      enddo

C----step 1 calculate D00
 101  continue
C--- Fixes D00 according to 5.42 Denner-Dittmaier
C--- knowing D0 with corrections of order Delta*Dii

c--- calculate S00 (needs D0)
      do ep=-2,0
      S00(ep)=2d0*Cv(cc0+C234,ep)+2d0*m1*Dv(dd0+N0,ep)
      enddo

      call run_00(kmax,lmax,DetGr,f,Gtwiddle,Gtt,Shat1,Shat2,S00,N0)
C--- Fixes Di according to 5.43 Denner-Dittmaier
C--- knowing D00 with corrections of order Delta*Dii
      do i2=1,np
      call run_i(jmax,i2,DetGr,Xtwiddle0,Gtwiddle,Shat2,N0)
      enddo

C----step 0 calculate D0
 100  continue
      call run_0(jmax,DetGr,Xtwiddle0,Gtwiddle,Shat1,N0)

c--- check the contents of box array    
c      write(6,*) 'recur: D array'
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
      
c--- check the contents of triangle arrays   
c      write(6,*) p2,p3,p2p3,m2,m3,m4 
c      do ip=1,Ncc
c      if (abs(Csing(ip,p2p3,p2,p3,m2,m3,m4)).ne.0d0) then
c      write(6,'(i3,4f20.15)') ip,Cv(ip+C234,-1),Cv(ip+C234,-1)
c     .    /Csing(ip,p2p3,p2,p3,m2,m3,m4)
c      endif
c      enddo
c      pause

c      write(6,*) p1p2,p3,p4,m1,m3,m4
c      do ip=1,Ncc
c      if (abs(Csing(ip,p4,p1p2,p3,m1,m3,m4)).ne.0d0) then
c      write(6,'(i3,4f20.15)') ip,Cv(ip+C134,-1),Cv(ip+C134,-1)
c     .    /Csing(ip,p4,p1p2,p3,m1,m3,m4)
c      endif
c      enddo
c      pause

c      write(6,*) p1,p2p3,p4,m1,m2,m4
c      do ip=1,Ncc
c      if (abs(Csing(ip,p4,p1,p2p3,m1,m2,m4)).ne.0d0) then
c      write(6,'(i3,4f20.15)') ip,Cv(ip+C124,-1),Cv(ip+C124,-1)
c     .    /Csing(ip,p4,p1,p2p3,m1,m2,m4)
c      endif
c      enddo
c      pause

c      write(6,*) p1,p2,p1p2,m1,m2,m3
c      do ip=1,Ncc
c      if (abs(Csing(ip,p1p2,p1,p2,m1,m2,m3)).ne.0d0) then
c      write(6,'(i3,4f20.15)') ip,Cv(ip+C123,-1),Cv(ip+C123,-1)
c     .    /Csing(ip,p1p2,p1,p2,m1,m2,m3)
c      endif
c      enddo
c      pause

c   77 format(a3,i2,a5,3('(',e13.6,',',e13.6,') '))
    
      return
      end


