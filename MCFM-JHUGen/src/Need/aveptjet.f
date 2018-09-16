      function aveptjet(p)
      implicit none
      include 'types.f'
      real(dp):: aveptjet

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      include 'jetlabel.f'
      integer:: j,countjet,isub,oldjets
      logical:: is_hadronic
      real(dp):: p(mxpart,4),pjet(mxpart,4),pt,rcut
      common/rcut/rcut

      aveptjet=0._dp

      if (abs(p(npart+2,4)) > 1.e-8_dp) then
        isub=0  ! real term
      else
        isub=1  ! subtraction term
      endif

c-- cluster jets but make sure recorded number of jets is not changed
      oldjets=jets
      call genclust2(p,rcut,pjet,isub)

      countjet=0
      do j=3,npart+2
        if (countjet == jets) goto 99
        if (is_hadronic(j)) then
          countjet=countjet+1
          aveptjet=aveptjet+pt(j,pjet)
        endif
      enddo

   99 continue

c--- restore old value of jets
      jets=oldjets

c--- dummy value returned if countjet=0, since this process
c--- must have nqcdjets > 0 - so this point will be dumped anyway
      if (countjet == 0) then
        aveptjet=10._dp
        return
      endif

      aveptjet=aveptjet/real(countjet,dp)

      return
      end


