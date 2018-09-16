      subroutine maketaucut_bb(pparton,jets,isub,passed)
!====== C Williams Oct 2015
!======Special version of the routine maketaucuts designed for
!=====Diboson style processes with hadronic decays at LO (i.e. not included
!=====in the tau stuff, some of the jet checks from maketaucut.f are not
!=====ideal.
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'taucut.f'
      include 'plabel.f'
      include 'first.f'
      logical passed
      integer i,j,isub,ihard,jets
      real(dp) :: pparton(mxpart,4),pttwo,
     &     ptjet(3),nn1(4),nn2(4),nn3(4),tau,taua,taub,tauj
      integer, save :: ipp
!$omp threadprivate(ipp)
      logical failedndp 
      common/failedndp/failedndp 
!$omp threadprivate(/failedndp/) 


      if(ntau.ne.0) then
         write(6,*) 'maketaucut_bb.f is designed only for'
         write(6,*) 'use when ntau=0 found ntau=',ntau
         stop
      endif
     
      
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

      
!---  Demand at least (ntau) hard jets after implementing the jet
!---  algorithm, otherwise the event is removed
      if (jets .lt. ntau) return

      
c--- if there are not enough partons (i.e. we're doing virtual) then return
      if ((npart .eq. ipp-3+ntau) .and. (isub .eq. 0)) then
        passed=.true.
        return
      endif

!=====note b quarks not included as reference directions since
!=====LO decay
      
c--- reference directions for beam axes
      nn1(:)=(/0._dp,0._dp, 1._dp,1._dp/)
      nn2(:)=(/0._dp,0._dp,-1._dp,1._dp/)

      
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
     &        -pparton(j,2)*nn3(2)-pparton(j,3)*nn3(3)
        tau=tau+min(abs(taua),abs(taub),abs(tauj))
      endif
      enddo
c---  check to make sure no NaN
      if (tau .ne. tau) then
         call writeout(pparton)
        write(6,*) 'maketaucut.f:  tau=',tau
        stop
      endif
c     write(6,*) 'tau',tau
      
c---  if tau is too small, implement the jettiness cut to remove the event
      if (tau .lt. taucut) return

c---  at this point cut has been passed
      passed=.true.
c     if (failedndp) then 
c     write(6,*) 'WARNING, phase space point has passed, '
c         write(6,*) 'both tau cut AND has failed ndp check'
c     endif
      
      return
      end
      
