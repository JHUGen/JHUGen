      subroutine scaleset_HT(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c--- partonic HT (i.e. scalar sum of particle pt's, whether or not they pass cuts)
      implicit none
      include 'constants.f'
      include 'npart.f'
      integer j
      double precision p(mxpart,4),mu0,pt

      mu0=0d0
      do j=3,npart+2
      mu0=mu0+pt(j,p)
      enddo

      return
      end
      
