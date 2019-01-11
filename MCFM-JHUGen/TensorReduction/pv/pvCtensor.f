      subroutine pvCtensor(q1,q2,m1s,m2s,m3s,
     . FC0,FC1,FC2,FC3,FC4,FC5,FC6)
      implicit none
C     q1,q2 are the momenta in the propagators
C     m1s,m2s,m3s are the squares of the masses in the propagators
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'TRydef.f'
      include 'TRmaxindex.f'
      include 'pvRespectmaxcindex.f'
      include 'pvrecurflags.f'
      include 'TRbadpoint.f'
      include 'pvforcerecalc.f'
      include 'pvverbose.f'
      include 'pvCitry.f'
      include 'TRmetric.f'
      double complex FC0(-2:0),FC1(y1max,-2:0),FC2(y2max,-2:0),
     . FC3(y3max,-2:0),FC4(y4max,-2:0),FC5(y5max,-2:0),FC6(y6max,-2:0),
     . FC7(y7max,-2:0)
      double precision p1Dp1,p2Dp2,p3Dp3,q1(4),q2(4)
      double precision m1s,m2s,m3s,
     . pvSPK,
     . pvSPKK,pvSDDP,
     . pvSPKKK,pvSPPKK,
     . pvSDDPP,pvSDDPK,pvSDDDD,
     . pvSDDDDP,pvSDDPPP,pvSDDPPK,pvSPPPKK,pvSPPPPK,
     . pvSDDDDDD,pvSDDDDPP,pvSDDDDPK,pvSDDPPPP,pvSDDPPPK,pvSDDPPKK,
     . pvSDDPKKK,pvSPPPPPK,pvSPPPPKK,pvSPPPKKK
      double precision q1save(4),q2save(4)
      integer n1,n2,n3,n4,n5,n6,n7,ep,C0i,pvCcache,itry
      logical failed
      common/q12save/q1save,q2save
      logical,save:: first=.true.
!$omp threadprivate(first,/q12save/)

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

c      p1(nu)=q1(nu)
c      p2(nu)=q2(nu)-q1(nu)

      p1Dp1=q1(4)**2-q1(1)**2-q1(2)**2-q1(3)**2
      p3Dp3=q2(4)**2-q2(1)**2-q2(2)**2-q2(3)**2
      p2Dp2=p1Dp1+p3Dp3
     . -2d0*(q1(4)*q2(4)-q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3))

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
      if (pvverbose) write(6,*) 'flag: badpoint  set in pvCtensor'
      return  
      endif
      
      C0i=pvCcache(p1Dp1,p2Dp2,p3Dp3,m1s,m2s,m3s)
c      write(6,*) 'C cache for ',p1Dp1,p2Dp2,p3Dp3,m1s,m2s,m3s
c      do n1=1,Ncc
c      write(6,'(i4,2f28.16)') n1,Cv(C0i+n1,0)
c      enddo
c      write(6,*)
      
      do ep=-2,0
      FC0(ep)=Cv(C0i+cc0,ep)
      enddo

      do ep=-2,0
      do n1=1,4
      FC1(n1,ep)=+Cv(C0i+cc1,ep)*q1(n1)+Cv(C0i+cc2,ep)*q2(n1)
      enddo
      enddo


      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      FC2(y2(n1,n2),ep)=
     . +q1(n1)*q1(n2)*Cv(C0i+cc11,ep)
     . +q2(n1)*q2(n2)*Cv(C0i+cc22,ep)
     . +pvSPK(n1,n2,q1,q2)*Cv(C0i+cc12,ep)
     . +g(n1,n2)*Cv(C0i+cc00,ep)
      enddo
      enddo
      enddo

      if ((maxcindex .eq. 2) .and. (pvRespectmaxcindex)) then
        if (pvCitry(C0i) .eq. -1) then
        call pvCcheck(2,q1,q2,m1s,m2s,m3s,
     &   FC0,FC1,FC2,FC3,FC4,FC5,FC6,failed)
      if (failed) then
c          write(6,*) 'recursion in pvCtensor'
        itry=itry+1
        goto 11
      endif
      pvCitry(C0i)=itry
c        if (itry .gt. 0) then
c        write(6,*) 'pvCtensor success: itry=',itry
c        endif
        endif
        return
      endif

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      FC3(y3(n1,n2,n3),ep)=
     . +q1(n1)*q1(n2)*q1(n3)*Cv(C0i+cc111,ep)
     . +q2(n1)*q2(n2)*q2(n3)*Cv(C0i+cc222,ep)
     . +pvSPKK(n1,n2,n3,q2,q1)*Cv(C0i+cc112,ep)
     . +pvSPKK(n1,n2,n3,q1,q2)*Cv(C0i+cc122,ep)
     . +pvSDDP(n1,n2,n3,q1)*Cv(C0i+cc001,ep)
     . +pvSDDP(n1,n2,n3,q2)*Cv(C0i+cc002,ep)
      enddo
      enddo
      enddo
      enddo

c--- This code can be used to check that the calculated tensor (e.g.
c---  by one of the recursion methods) satisfies the usual PV relations
c      call PVchecktri(q1,q2,m1s,m2s,m3s,FC0,FC1,FC2,FC3,FC4,FC5,FC6,failed)
c---

      if ((maxcindex .eq. 3) .and. (pvRespectmaxcindex)) then
        call pvCcheck(3,q1,q2,m1s,m2s,m3s,
     &   FC0,FC1,FC2,FC3,FC4,FC5,FC6,failed)
      if (failed) then
c        write(6,*) 'recursion in pvCtensor'
        itry=itry+1
        goto 11
      endif
        return
      endif

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      FC4(y4(n1,n2,n3,n4),ep)=
     . +q1(n1)*q1(n2)*q1(n3)*q1(n4)*Cv(C0i+cc1111,ep)
     . +q2(n1)*q2(n2)*q2(n3)*q2(n4)*Cv(C0i+cc2222,ep)
     . +pvSPKKK(n1,n2,n3,n4,q2,q1)*Cv(C0i+cc1112,ep)
     . +pvSPKKK(n1,n2,n3,n4,q1,q2)*Cv(C0i+cc1222,ep)
     . +pvSPPKK(n1,n2,n3,n4,q1,q2)*Cv(C0i+cc1122,ep)
     . +pvSDDPP(n1,n2,n3,n4,q1)*Cv(C0i+cc0011,ep)
     . +pvSDDPP(n1,n2,n3,n4,q2)*Cv(C0i+cc0022,ep)
     . +pvSDDPK(n1,n2,n3,n4,q1,q2)*Cv(C0i+cc0012,ep)
     . +pvSDDDD(n1,n2,n3,n4)*Cv(C0i+cc0000,ep)


      enddo
      enddo
      enddo
      enddo
      enddo

      if ((maxcindex .eq. 4) .and. (pvRespectmaxcindex)) then
        call pvCcheck(4,q1,q2,m1s,m2s,m3s,
     &   FC0,FC1,FC2,FC3,FC4,FC5,FC6,failed)
      if (failed) then
c        write(6,*) 'recursion in pvCtensor'
        itry=itry+1
        goto 11
      endif
        return
      endif

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4

      FC5(y5(n1,n2,n3,n4,n5),ep)=
     . +Cv(C0i+cc00001,ep)*pvSDDDDP(n1,n2,n3,n4,n5,q1)
     . +Cv(C0i+cc00002,ep)*pvSDDDDP(n1,n2,n3,n4,n5,q2)
     . +Cv(C0i+cc00111,ep)*pvSDDPPP(n1,n2,n3,n4,n5,q1)
     . +Cv(C0i+cc00112,ep)*pvSDDPPK(n1,n2,n3,n4,n5,q1,q2)
     . +Cv(C0i+cc00122,ep)*pvSDDPPK(n1,n2,n3,n4,n5,q2,q1)
     . +Cv(C0i+cc00222,ep)*pvSDDPPP(n1,n2,n3,n4,n5,q2)
     . +Cv(C0i+cc11111,ep)*q1(n1)*q1(n2)*q1(n3)*q1(n4)*q1(n5)
     . +Cv(C0i+cc11112,ep)*pvSPPPPK(n1,n2,n3,n4,n5,q1,q2)
     . +Cv(C0i+cc11122,ep)*pvSPPPKK(n1,n2,n3,n4,n5,q1,q2)
     . +Cv(C0i+cc11222,ep)*pvSPPPKK(n1,n2,n3,n4,n5,q2,q1)
     . +Cv(C0i+cc12222,ep)*pvSPPPPK(n1,n2,n3,n4,n5,q2,q1)
     . +Cv(C0i+cc22222,ep)*q2(n1)*q2(n2)*q2(n3)*q2(n4)*q2(n5)


      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      if ((maxcindex .eq. 5) .and. (pvRespectmaxcindex)) then
        call pvCcheck(5,q1,q2,m1s,m2s,m3s,
     &   FC0,FC1,FC2,FC3,FC4,FC5,FC6,failed)
      if (failed) then
c        write(6,*) 'recursion in pvCtensor'
        itry=itry+1
        goto 11
      endif
        return
      endif

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      FC6(y6(n1,n2,n3,n4,n5,n6),ep)=
     .  +Cv(C0i+cc000000,ep)*pvSDDDDDD(n1,n2,n3,n4,n5,n6)
     .  +Cv(C0i+cc000011,ep)*pvSDDDDPP(n1,n2,n3,n4,n5,n6,q1)
     .  +Cv(C0i+cc000012,ep)*pvSDDDDPK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Cv(C0i+cc000022,ep)*pvSDDDDPP(n1,n2,n3,n4,n5,n6,q2)
     .  +Cv(C0i+cc001111,ep)*pvSDDPPPP(n1,n2,n3,n4,n5,n6,q1)
     .  +Cv(C0i+cc001112,ep)*pvSDDPPPK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Cv(C0i+cc001122,ep)*pvSDDPPKK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Cv(C0i+cc001222,ep)*pvSDDPKKK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Cv(C0i+cc002222,ep)*pvSDDPPPP(n1,n2,n3,n4,n5,n6,q2)
     .  +Cv(C0i+cc111111,ep)*q1(n1)*q1(n2)*q1(n3)*q1(n4)*q1(n5)*q1(n6)
     .  +Cv(C0i+cc111112,ep)*pvSPPPPPK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Cv(C0i+cc111122,ep)*pvSPPPPKK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Cv(C0i+cc111222,ep)*pvSPPPKKK(n1,n2,n3,n4,n5,n6,q1,q2)
     .  +Cv(C0i+cc112222,ep)*pvSPPPPKK(n1,n2,n3,n4,n5,n6,q2,q1)
     .  +Cv(C0i+cc122222,ep)*pvSPPPPPK(n1,n2,n3,n4,n5,n6,q2,q1)
     .  +Cv(C0i+cc222222,ep)*q2(n1)*q2(n2)*q2(n3)*q2(n4)*q2(n5)*q2(n6)


      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      if ((maxcindex .eq. 6) .and. (pvRespectmaxcindex)) then
        if (pvCitry(C0i) .eq. -1) then
        call pvCcheck(6,q1,q2,m1s,m2s,m3s,
     &   FC0,FC1,FC2,FC3,FC4,FC5,FC6,failed)
      if (failed) then
c          write(6,*) 'recursion in pvCtensor'
        itry=itry+1
        goto 11
      endif
      pvCitry(C0i)=itry
c        if (itry .gt. 0) then
c        write(6,*) 'pvCtensor success: itry=',itry
c        endif
        endif
        return
      endif

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      do n7=n6,4
c--- crude implementation of FC7 (used for checking only)
      include 'pvFC7.f'
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      
      return
      end

