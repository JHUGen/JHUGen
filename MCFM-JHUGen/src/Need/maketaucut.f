      subroutine maketaucut(pparton,pjet,jets,isub,passed)
c---- J. Campbell, April 2015
c----
c---- Adapted from code written by R. Boughezal et al. for
c---- 1-jettiness calculation in W+jet events
c----
c---- Implements cut on 1-jettiness for processes with nqcdjets=1
c---- Implements cut on 0-jettiness for processes with nqcdjets=0
c----
c---  Only keeps events with tau > taucut (passed via common block)
c---  tau < taucut will be treated separately using SCET
c----
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'nqcdjets.f'
      include 'taucut.f'
      include 'plabel.f'
      include 'first.f'
      logical passed
      integer i,j,isub,ihard,jets
      real(dp) pparton(mxpart,4),pjet(mxpart,4),pttwo,
     & ptjet(3),nn1(4),nn2(4),nn3(4),tau,taua,taub,tauj
      integer, save :: ipp
!$omp threadprivate(ipp)
      logical failedndp 
      common/failedndp/failedndp 
!$omp threadprivate(/failedndp/) 

c--- determine beginning of parton entries in plabel
      if (first) then
        first=.false.
        ipp=3
        do while ((ipp < mxpart ) .and. (plabel(ipp) .ne. 'pp'))
          ipp=ipp+1
        enddo
        if (ipp == mxpart) then
          write(6,*) 'Could not identify partons in maketaucut.f'
          stop
        endif
c        write(6,*) 'found ipp=',ipp  
      endif
      
      passed=.false.

c--- safety cut on pt(H)
c      if (ntau == 1) then
c        if (pttwo(3,4,pparton) < 1._dp) return
c      endif

!---  Demand at least (ntau) hard jets after implementing the jet
!---  algorithm, otherwise the event is removed
      if (jets .lt. ntau) return
      
c--- if there are not enough partons (i.e. we're doing virtual) then return
      if ((npart .eq. ipp-3+ntau) .and. (isub .eq. 0)) then
        passed=.true.
        return
      endif

!--- routine is not designed for more than 3 jets
      if (jets .gt. 3) then
        write(6,*) 'Error: >3 jets found in maketaucut.f'
        write(6,*) ' jets=',jets
        stop
      endif

!--- routine is not designed for nqcdjets>2
      if (nqcdjets .gt. 2) then
        write(6,*) 'Error: unimplemented nqcdjets in maketaucut.f'
        write(6,*) ' nqcdjets=',nqcdjets
        stop
      endif

c--- reference directions for beam axes
      nn1(:)=(/0._dp,0._dp, 1._dp,1._dp/)
      nn2(:)=(/0._dp,0._dp,-1._dp,1._dp/)

c---  compute hard jet reference direction for 1-jettiness
      if (ntau > 0) then
        ptjet(:)=0._dp
        do i=0,jets-1
          ptjet(i+1)=sqrt(pjet(ipp+i,1)**2+pjet(ipp+i,2)**2)
        enddo
        ihard=ipp
        if (jets > 1) then
          if (ptjet(2) .gt. max(ptjet(1),ptjet(3))) ihard=ipp+1
        endif
        if (jets > 2) then
          if (ptjet(3) .gt. max(ptjet(1),ptjet(2))) ihard=ipp+2
        endif
        nn3(1:3)=pjet(ihard,1:3)
     &    /sqrt(pjet(ihard,1)**2+pjet(ihard,2)**2+pjet(ihard,3)**2)
        nn3(4)=1.0_dp
      endif

      tau=0._dp
c--- compute N-jettiness tau, by dotting the final state QCD parton
c--- 4-momentum into ni, for i=1,2 and 3 (ntau>0 only)
      do j = ipp,npart+2-isub     
      taua =pparton(j,4)*nn1(4)-pparton(j,1)*nn1(1)
     &     -pparton(j,2)*nn1(2)-pparton(j,3)*nn1(3)
      taub=pparton(j,4)*nn2(4)-pparton(j,1)*nn2(1)
     &    -pparton(j,2)*nn2(2)-pparton(j,3)*nn2(3)
      if (ntau .eq. 0) then
        tau=tau+min(abs(taua),abs(taub))
      else
        tauj=pparton(j,4)*nn3(4)-pparton(j,1)*nn3(1)
     &      -pparton(j,2)*nn3(2)-pparton(j,3)*nn3(3)
        tau=tau+min(abs(taua),abs(taub),abs(tauj))
      endif
      enddo

c--- check to make sure no NaN
      if (tau .ne. tau) then
        call writeout(pparton)
        write(6,*) 'maketaucut.f:  tau=',tau
        stop
      endif
c      write(6,*) 'tau',tau
      
c--- if tau is too small, implement the jettiness cut to remove the event
      if (tau .lt. taucut) return

c--- at this point cut has been passed
      passed=.true.
c      if (failedndp) then 
c         write(6,*) 'WARNING, phase space point has passed, '
c         write(6,*) 'both tau cut AND has failed ndp check'
c      endif

      return
      end
      
