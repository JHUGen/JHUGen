      subroutine pvDtensor(q1,q2,q3,m1s,m2s,m3s,m4s,
     . FD0,FD1,FD2,FD3,FD4,FD5,FD6)
      implicit none
C     q1,q2,q3 are the loop offset momenta
C     m1s,m2s,m3s,m4s are the squares of the masses in the propagators
      include 'pvDnames.f'
      include 'pvDv.f'
      include 'TRydef.f'
      include 'TRmaxindex.f'
      include 'pvrecurflags.f'
      include 'TRbadpoint.f'
      include 'pvforcerecalc.f'
      include 'pvverbose.f'
      include 'pvDitry.f'
      include 'TRmetric.f'
      double complex FD0(-2:0),FD1(y1max,-2:0),FD2(y2max,-2:0),
     . FD3(y3max,-2:0),FD4(y4max,-2:0),FD5(y5max,-2:0),FD6(y6max,-2:0)
      double precision p1(4),p2(4),p3(4),p4(4),p1Dp1,p2Dp2,p3Dp3,p4Dp4,
     . s12,s23,q1(4),q2(4),q3(4),p12(4),p23(4),m1s,m2s,m3s,m4s
      double precision pvSPK,pvSPKL,pvSPKK,pvSDDP,
     . pvSPKKK,pvSPPKK,pvSPPKL,pvSDDPP,pvSDDPK,pvSDDDD,
     . pvSDDDDP,pvSDDPPP,pvSDDPPK,pvSDDPKL,
     . pvSPPPKK,pvSPPPPK,pvSPPPKL,pvSPPKKL,
     . pvSDDDDDD,pvSDDDDPP,pvSDDDDPK,pvSDDPPPP,pvSDDPPPK,pvSDDPPKK,
     . pvSDDPPKL,pvSDDPKKK,pvSPPPPPK,pvSPPPPKK,pvSPPPKKK,pvSPPPPKL,
     . pvSPPPKKL,pvSPPKKLL
      integer nu,n1,n2,n3,n4,n5,n6,ep,pvDcache,D01,itry
      logical failed
      double precision q1save(4),q2save(4),q3save(4)
      common/q123save/q1save,q2save,q3save  
      logical,save:: first=.true.
!$omp threadprivate(first,/q123save/)

c--- value of itry specifies which calculation to use for tensor coefficients:
c---    itry=0            regular PV recursion
c---    itry=1            small G recursion
c---    itry=2            small G and small Y recursion

      itry=0      ! Regular PV recursion is the default
      doGsing =.false.
      doGYsing=.false.
      doPsing =.false.
      doPFsing=.false.
      
      if (first) then
      first=.false.
      call pvarraysetup
      endif
    
      q1save(:)=q1(:)
      q2save(:)=q2(:)
      q3save(:)=q3(:)    
    
      call pvYcalc(q1,q2,q3,m1s,m2s,m3s,m4s)


C     p1,p2,p3,p4 are the external momenta
      do nu=1,4
      p1(nu)=+q1(nu)
      p2(nu)=+q2(nu)-q1(nu)
      p3(nu)=+q3(nu)-q2(nu)
      p4(nu)=-q3(nu)
      p23(nu)=p2(nu)+p3(nu)
      p12(nu)=p1(nu)+p2(nu)
      enddo

      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      p2Dp2=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      p3Dp3=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p4Dp4=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s23=p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2

c--- point from which to continue
   11 continue

      pvforcerecalc=.false.
      if     (itry .eq. 1) then
        doGsing=.true.
        pvforcerecalc=.true.
      elseif (itry .eq. 2) then
        doGsing=.false.
        doGYsing=.true.
        pvforcerecalc=.true.
      elseif (itry .eq. 3) then
c--- cannot compute point: set flag
        pvbadpoint=.true.   
      if (pvverbose) write(6,*) 'flag: badpoint  set in pvDtensor'
      return  
      endif
      
      D01=pvDcache(p1Dp1,p2Dp2,p3Dp3,p4Dp4,s12,s23,m1s,m2s,m3s,m4s)
      
      do ep=-2,0
      FD0(ep)=Dv(D01+dd0,ep)
      enddo
c      write(6,*) 'Dtensor,FD0',FD0,D01

c   Id,FD1(n1?,q1?,q2?,q3?,?x)=
c    +d_(n1,q1)*D1(q1,q2,q3,?x)
c    +d_(n1,q2)*D2(q1,q2,q3,?x)
c    +d_(n1,q3)*D3(q1,q2,q3,?x);

      do ep=-2,0
      do n1=1,4
      FD1(n1,ep)=
     . +Dv(D01+dd1,ep)*q1(n1)
     . +Dv(D01+dd2,ep)*q2(n1)
     . +Dv(D01+dd3,ep)*q3(n1)

      enddo
      enddo

c   Id,FD2(n1?,n2?,q1?,q2?,q3?,?x)=
c    +d_(n1,q1)*d_(n2,q1)*D11(q1,q2,q3,?x)
c    +d_(n1,q2)*d_(n2,q2)*D22(q1,q2,q3,?x)
c    +d_(n1,q3)*d_(n2,q3)*D33(q1,q2,q3,?x)
c    +pvSPK(n1,n2,q1,q2)*D12(q1,q2,q3,?x)
c    +pvSPK(n1,n2,q1,q3)*D13(q1,q2,q3,?x)
c    +pvSPK(n1,n2,q2,q3)*D23(q1,q2,q3,?x)
c    +d_(n1,n2)*D00(q1,q2,q3,?x);

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      FD2(y2(n1,n2),ep)=
     . +q1(n1)*q1(n2)*Dv(D01+dd11,ep)
     . +q2(n1)*q2(n2)*Dv(D01+dd22,ep)
     . +q3(n1)*q3(n2)*Dv(D01+dd33,ep)
     . +pvSPK(n1,n2,q1,q2)*Dv(D01+dd12,ep)
     . +pvSPK(n1,n2,q1,q3)*Dv(D01+dd13,ep)
     . +pvSPK(n1,n2,q2,q3)*Dv(D01+dd23,ep)
     . +g(n1,n2)*Dv(D01+dd00,ep)

      enddo
      enddo
      enddo

      if (maxdindex .eq. 2) then
        call pvDcheck(2,q1,q2,q3,m1s,m2s,m3s,m4s,
     &   FD0,FD1,FD2,FD3,FD4,FD5,FD6,failed)
      if (failed) then
        itry=itry+1
        goto 11
      endif
        return
      endif

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      FD3(y3(n1,n2,n3),ep)=
     . +q1(n1)*q1(n2)*q1(n3)*Dv(D01+dd111,ep)
     . +q2(n1)*q2(n2)*q2(n3)*Dv(D01+dd222,ep)
     . +q3(n1)*q3(n2)*q3(n3)*Dv(D01+dd333,ep)
     . +pvSPKK(n1,n2,n3,q2,q1)*Dv(D01+dd112,ep)
     . +pvSPKK(n1,n2,n3,q3,q1)*Dv(D01+dd113,ep)
     . +pvSPKK(n1,n2,n3,q1,q2)*Dv(D01+dd122,ep)
     . +pvSPKK(n1,n2,n3,q1,q3)*Dv(D01+dd133,ep)
     . +pvSPKK(n1,n2,n3,q3,q2)*Dv(D01+dd223,ep)
     . +pvSPKK(n1,n2,n3,q2,q3)*Dv(D01+dd233,ep)
     . +pvSPKL(n1,n2,n3,q1,q2,q3)*Dv(D01+dd123,ep)
     . +pvSDDP(n1,n2,n3,q1)*Dv(D01+dd001,ep)
     . +pvSDDP(n1,n2,n3,q2)*Dv(D01+dd002,ep)
     . +pvSDDP(n1,n2,n3,q3)*Dv(D01+dd003,ep)

      enddo
      enddo
      enddo
      enddo

      if (maxdindex .eq. 3) then
        if (pvDitry(D01) .eq. -1) then
        call pvDcheck(3,q1,q2,q3,m1s,m2s,m3s,m4s,
     &   FD0,FD1,FD2,FD3,FD4,FD5,FD6,failed)
      if (failed) then
c        write(6,*) 'recursion in pvDtensor'
        itry=itry+1
        goto 11
      endif
      pvDitry(D01)=itry
c      write(6,*) 'pvDtensor perfect'
        endif
        return
      endif

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      FD4(y4(n1,n2,n3,n4),ep)=
     . +q1(n1)*q1(n2)*q1(n3)*q1(n4)*Dv(D01+dd1111,ep)
     . +q2(n1)*q2(n2)*q2(n3)*q2(n4)*Dv(D01+dd2222,ep)
     . +q3(n1)*q3(n2)*q3(n3)*q3(n4)*Dv(D01+dd3333,ep)
     . +pvSPKKK(n1,n2,n3,n4,q2,q1)*Dv(D01+dd1112,ep)
     . +pvSPKKK(n1,n2,n3,n4,q3,q1)*Dv(D01+dd1113,ep)
     . +pvSPKKK(n1,n2,n3,n4,q1,q2)*Dv(D01+dd1222,ep)
     . +pvSPKKK(n1,n2,n3,n4,q3,q2)*Dv(D01+dd2223,ep)
     . +pvSPKKK(n1,n2,n3,n4,q1,q3)*Dv(D01+dd1333,ep)
     . +pvSPKKK(n1,n2,n3,n4,q2,q3)*Dv(D01+dd2333,ep)

     . +pvSPPKK(n1,n2,n3,n4,q1,q2)*Dv(D01+dd1122,ep)
     . +pvSPPKK(n1,n2,n3,n4,q1,q3)*Dv(D01+dd1133,ep)
     . +pvSPPKK(n1,n2,n3,n4,q2,q3)*Dv(D01+dd2233,ep)

     . +pvSPPKL(n1,n2,n3,n4,q1,q2,q3)*Dv(D01+dd1123,ep)
     . +pvSPPKL(n1,n2,n3,n4,q2,q1,q3)*Dv(D01+dd1223,ep)
     . +pvSPPKL(n1,n2,n3,n4,q3,q1,q2)*Dv(D01+dd1233,ep)

     . +pvSDDPP(n1,n2,n3,n4,q1)*Dv(D01+dd0011,ep)
     . +pvSDDPP(n1,n2,n3,n4,q2)*Dv(D01+dd0022,ep)
     . +pvSDDPP(n1,n2,n3,n4,q3)*Dv(D01+dd0033,ep)

     . +pvSDDPK(n1,n2,n3,n4,q1,q2)*Dv(D01+dd0012,ep)
     . +pvSDDPK(n1,n2,n3,n4,q2,q3)*Dv(D01+dd0023,ep)
     . +pvSDDPK(n1,n2,n3,n4,q1,q3)*Dv(D01+dd0013,ep)

     . +pvSDDDD(n1,n2,n3,n4)*Dv(D01+dd0000,ep)

      enddo
      enddo
      enddo
      enddo
      enddo

      if (maxdindex .eq. 4) then
        if (pvDitry(D01) .eq. -1) then
        call pvDcheck(4,q1,q2,q3,m1s,m2s,m3s,m4s,
     &   FD0,FD1,FD2,FD3,FD4,FD5,FD6,failed)
      if (failed) then
c        write(6,*) 'recursion in pvDtensor'
        itry=itry+1
        goto 11
      endif
      pvDitry(D01)=itry
c      write(6,*) 'pvDtensor perfect'
        endif
        return
      endif
      
      if (maxdindex .eq. 4) return

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      FD5(y5(n1,n2,n3,n4,n5),ep)=
     . +Dv(D01+dd00001,ep)*pvSDDDDP(n1,n2,n3,n4,n5,q1)
     . +Dv(D01+dd00002,ep)*pvSDDDDP(n1,n2,n3,n4,n5,q2)
     . +Dv(D01+dd00003,ep)*pvSDDDDP(n1,n2,n3,n4,n5,q3)
     . +Dv(D01+dd00111,ep)*pvSDDPPP(n1,n2,n3,n4,n5,q1)
     . +Dv(D01+dd00112,ep)*pvSDDPPK(n1,n2,n3,n4,n5,q1,q2)
     . +Dv(D01+dd00113,ep)*pvSDDPPK(n1,n2,n3,n4,n5,q1,q3)
     . +Dv(D01+dd00122,ep)*pvSDDPPK(n1,n2,n3,n4,n5,q2,q1)
     . +Dv(D01+dd00123,ep)*pvSDDPKL(n1,n2,n3,n4,n5,q1,q2,q3)
     . +Dv(D01+dd00133,ep)*pvSDDPPK(n1,n2,n3,n4,n5,q3,q1)
     . +Dv(D01+dd00222,ep)*pvSDDPPP(n1,n2,n3,n4,n5,q2)
     . +Dv(D01+dd00223,ep)*pvSDDPPK(n1,n2,n3,n4,n5,q2,q3)
     . +Dv(D01+dd00233,ep)*pvSDDPPK(n1,n2,n3,n4,n5,q3,q2)
     . +Dv(D01+dd00333,ep)*pvSDDPPP(n1,n2,n3,n4,n5,q3)
     . +Dv(D01+dd11111,ep)*q1(n1)*q1(n2)*q1(n3)*q1(n4)*q1(n5)
     . +Dv(D01+dd11112,ep)*pvSPPPPK(n1,n2,n3,n4,n5,q1,q2)
     . +Dv(D01+dd11113,ep)*pvSPPPPK(n1,n2,n3,n4,n5,q1,q3)
     . +Dv(D01+dd11122,ep)*pvSPPPKK(n1,n2,n3,n4,n5,q1,q2)
     . +Dv(D01+dd11123,ep)*pvSPPPKL(n1,n2,n3,n4,n5,q1,q2,q3)
     . +Dv(D01+dd11133,ep)*pvSPPPKK(n1,n2,n3,n4,n5,q1,q3)
     . +Dv(D01+dd11222,ep)*pvSPPPKK(n1,n2,n3,n4,n5,q2,q1)
     . +Dv(D01+dd11223,ep)*pvSPPKKL(n1,n2,n3,n4,n5,q1,q2,q3)
     . +Dv(D01+dd11233,ep)*pvSPPKKL(n1,n2,n3,n4,n5,q1,q3,q2)
     . +Dv(D01+dd11333,ep)*pvSPPPKK(n1,n2,n3,n4,n5,q3,q1)
     . +Dv(D01+dd12222,ep)*pvSPPPPK(n1,n2,n3,n4,n5,q2,q1)
     . +Dv(D01+dd12223,ep)*pvSPPPKL(n1,n2,n3,n4,n5,q2,q1,q3)
     . +Dv(D01+dd12233,ep)*pvSPPKKL(n1,n2,n3,n4,n5,q2,q3,q1)
     . +Dv(D01+dd12333,ep)*pvSPPPKL(n1,n2,n3,n4,n5,q3,q1,q2)
     . +Dv(D01+dd13333,ep)*pvSPPPPK(n1,n2,n3,n4,n5,q3,q1)
     . +Dv(D01+dd22222,ep)*q2(n1)*q2(n2)*q2(n3)*q2(n4)*q2(n5)
     . +Dv(D01+dd22223,ep)*pvSPPPPK(n1,n2,n3,n4,n5,q2,q3)
     . +Dv(D01+dd22233,ep)*pvSPPPKK(n1,n2,n3,n4,n5,q2,q3)
     . +Dv(D01+dd22333,ep)*pvSPPPKK(n1,n2,n3,n4,n5,q3,q2)
     . +Dv(D01+dd23333,ep)*pvSPPPPK(n1,n2,n3,n4,n5,q3,q2)
     . +Dv(D01+dd33333,ep)*q3(n1)*q3(n2)*q3(n3)*q3(n4)*q3(n5) 
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      if (maxdindex .eq. 5) return

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      FD6(y6(n1,n2,n3,n4,n5,n6),ep)=
     .  +Dv(D01+dd000000,ep)*pvSDDDDDD(n1,n2,n3,n4,n5,n6)
     .  +Dv(D01+dd000011,ep)*pvSDDDDPP(n1,n2,n3,n4,n5,n6,q1)
     .  +Dv(D01+dd000012,ep)*pvSDDDDPK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Dv(D01+dd000013,ep)*pvSDDDDPK(n1,n2,n3,n4,n5,n6,q1,q3)
     .  +Dv(D01+dd000022,ep)*pvSDDDDPP(n1,n2,n3,n4,n5,n6,q2)
     .  +Dv(D01+dd000023,ep)*pvSDDDDPK(n1,n2,n3,n4,n5,n6,q2,q3)
     .  +Dv(D01+dd000033,ep)*pvSDDDDPP(n1,n2,n3,n4,n5,n6,q3)
     .  +Dv(D01+dd001111,ep)*pvSDDPPPP(n1,n2,n3,n4,n5,n6,q1)
     .  +Dv(D01+dd001112,ep)*pvSDDPPPK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Dv(D01+dd001113,ep)*pvSDDPPPK(n1,n2,n3,n4,n5,n6,q1,q3)
     .  +Dv(D01+dd001122,ep)*pvSDDPPKK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Dv(D01+dd001123,ep)*pvSDDPPKL(n1,n2,n3,n4,n5,n6,q1,q2,q3)
     .  +Dv(D01+dd001133,ep)*pvSDDPPKK(n1,n2,n3,n4,n5,n6,q1,q3)
     .  +Dv(D01+dd001222,ep)*pvSDDPKKK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Dv(D01+dd001223,ep)*pvSDDPPKL(n1,n2,n3,n4,n5,n6,q2,q1,q3)
     .  +Dv(D01+dd001233,ep)*pvSDDPPKL(n1,n2,n3,n4,n5,n6,q3,q1,q2)
     .  +Dv(D01+dd001333,ep)*pvSDDPPPK(n1,n2,n3,n4,n5,n6,q3,q1)
     .  +Dv(D01+dd002222,ep)*pvSDDPPPP(n1,n2,n3,n4,n5,n6,q2)
     .  +Dv(D01+dd002223,ep)*pvSDDPPPK(n1,n2,n3,n4,n5,n6,q2,q3)
     .  +Dv(D01+dd002233,ep)*pvSDDPPKK(n1,n2,n3,n4,n5,n6,q2,q3)
     .  +Dv(D01+dd002333,ep)*pvSDDPPPK(n1,n2,n3,n4,n5,n6,q3,q2)
     .  +Dv(D01+dd003333,ep)*pvSDDPPPP(n1,n2,n3,n4,n5,n6,q3)
     .  +Dv(D01+dd111111,ep)*q1(n1)*q1(n2)*q1(n3)*q1(n4)*q1(n5)*q1(n6)
     .  +Dv(D01+dd111112,ep)*pvSPPPPPK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Dv(D01+dd111113,ep)*pvSPPPPPK(n1,n2,n3,n4,n5,n6,q1,q3)
     .  +Dv(D01+dd111122,ep)*pvSPPPPKK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Dv(D01+dd111123,ep)*pvSPPPPKL(n1,n2,n3,n4,n5,n6,q1,q2,q3)
     .  +Dv(D01+dd111133,ep)*pvSPPPPKK(n1,n2,n3,n4,n5,n6,q1,q3)
     .  +Dv(D01+dd111222,ep)*pvSPPPKKK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Dv(D01+dd111223,ep)*pvSPPPKKL(n1,n2,n3,n4,n5,n6,q1,q2,q3)
     .  +Dv(D01+dd111233,ep)*pvSPPPKKL(n1,n2,n3,n4,n5,n6,q1,q3,q2)
     .  +Dv(D01+dd111333,ep)*pvSPPPKKK(n1,n2,n3,n4,n5,n6,q1,q3)
     .  +Dv(D01+dd112222,ep)*pvSPPPPKK(n1,n2,n3,n4,n5,n6,q2,q1)
     .  +Dv(D01+dd112223,ep)*pvSPPPKKL(n1,n2,n3,n4,n5,n6,q2,q1,q3)
     .  +Dv(D01+dd112233,ep)*pvSPPKKLL(n1,n2,n3,n4,n5,n6,q1,q2,q3)
     .  +Dv(D01+dd112333,ep)*pvSPPPKKL(n1,n2,n3,n4,n5,n6,q3,q1,q2)
     .  +Dv(D01+dd113333,ep)*pvSPPPPKK(n1,n2,n3,n4,n5,n6,q3,q1)
     .  +Dv(D01+dd122222,ep)*pvSPPPPPK(n1,n2,n3,n4,n5,n6,q2,q1)
     .  +Dv(D01+dd122223,ep)*pvSPPPPKL(n1,n2,n3,n4,n5,n6,q2,q1,q3)
     .  +Dv(D01+dd122233,ep)*pvSPPPKKL(n1,n2,n3,n4,n5,n6,q2,q3,q1)
     .  +Dv(D01+dd122333,ep)*pvSPPPKKL(n1,n2,n3,n4,n5,n6,q3,q2,q1)
     .  +Dv(D01+dd123333,ep)*pvSPPPPKL(n1,n2,n3,n4,n5,n6,q3,q1,q2)
     .  +Dv(D01+dd133333,ep)*pvSPPPPPK(n1,n2,n3,n4,n5,n6,q3,q1)
     .  +Dv(D01+dd222222,ep)*q2(n1)*q2(n2)*q2(n3)*q2(n4)*q2(n5)*q2(n6)
     .  +Dv(D01+dd222223,ep)*pvSPPPPPK(n1,n2,n3,n4,n5,n6,q2,q3)
     .  +Dv(D01+dd222233,ep)*pvSPPPPKK(n1,n2,n3,n4,n5,n6,q2,q3)
     .  +Dv(D01+dd222333,ep)*pvSPPPKKK(n1,n2,n3,n4,n5,n6,q2,q3)
     .  +Dv(D01+dd223333,ep)*pvSPPPPKK(n1,n2,n3,n4,n5,n6,q3,q2)
     .  +Dv(D01+dd233333,ep)*pvSPPPPPK(n1,n2,n3,n4,n5,n6,q3,q2)
     .  +Dv(D01+dd333333,ep)*q3(n1)*q3(n2)*q3(n3)*q3(n4)*q3(n5)*q3(n6)

      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end

