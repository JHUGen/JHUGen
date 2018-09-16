      subroutine scaleset_Msqpt345sq(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt345^2), where M is the mass of the particle (345)
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'breit.f'
      double precision p(mxpart,4),mu0

      if((case .eq. 'tt_bbl') .or.
     &   (case .eq. 'tt_bbu') .or.
     &   (case .eq. 'tt_bbh')) then
        mu0=(p(3,4)+p(4,4)+p(5,4))**2-(p(3,3)+p(4,3)+p(5,3))**2
        mu0=dsqrt(dabs(mu0))
      else
        write(6,*) 'dynamicscale sqrt(M^2+pt345^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
