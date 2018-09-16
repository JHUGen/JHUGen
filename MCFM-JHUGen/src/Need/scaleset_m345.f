      subroutine scaleset_m345(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3, 4 and 5
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      real(dp):: p(mxpart,4),mu0

      if((kcase==kWgamma) .or.
     &   (kcase==kZgamma) .or.
     &   (kcase==kWgajet) .or.
     &   (kcase==ktrigam) .or.
     &   (kcase==kZgajet)) then
        mu0=(p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
     &     -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2       
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale m(345) not supported for this process.'
        stop
      endif
      
      return
      end
      
