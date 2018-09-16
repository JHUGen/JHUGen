      ! functions for counting particles 

      function count_photo()
       implicit none
      include 'types.f'
      integer:: count_photo
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j
      logical:: is_photon
      external is_photon

      count_photo=0
      do j=1,mxpart
         if (is_photon(j)) then 
            count_photo=count_photo+1
         endif
      enddo 

      return
      end

      
      function count_jets()
       implicit none
      include 'types.f'
      integer:: count_jets
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      integer:: j
      logical:: is_hadronic
      external is_hadronic

c---- Count final state jets
      count_jets=0

      do j=3,npart+2
         if (is_hadronic(j)) then 
            count_jets=count_jets+1
         endif
      enddo 

      return 
      end
