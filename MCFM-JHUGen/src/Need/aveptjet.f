      double precision function aveptjet(p)
      implicit none
      include 'constants.f'
      include 'npart.f'
      include 'jetlabel.f'
      integer j,countjet,isub,oldjets
      logical is_hadronic
      double precision p(mxpart,4),pjet(mxpart,4),pt,rcut
      common/rcut/rcut
      
      aveptjet=0d0

      if (abs(p(npart+2,4)) .gt. 1d-8) then
        isub=0  ! real term
      else
        isub=1  ! subtraction term
      endif
      
c-- cluster jets but make sure recorded number of jets is not changed
      oldjets=jets     
      call genclust2(p,rcut,pjet,isub)
      
      countjet=0
      do j=3,npart+2
        if (countjet .eq. jets) goto 99
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
      if (countjet .eq. 0) then
        aveptjet=10d0
        return
      endif
     
      aveptjet=aveptjet/dfloat(countjet)
      
      return
      end
      
      
