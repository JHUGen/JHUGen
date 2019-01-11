      subroutine scaleset_shat(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3 and 4
      implicit none
      include 'constants.f'
      include 'part.f'
      double precision p(mxpart,4),mu0

      if (part .eq. 'lord') then
        mu0=(p(1,4)+p(2,4))**2-(p(1,1)+p(2,1))**2
     &     -(p(1,2)+p(2,2))**2-(p(1,3)+p(2,3))**2       
        mu0=dsqrt(dabs(mu0))
      else
        write(6,*) 'dynamicscale s-hat not supported beyond LO.'
        stop
      endif
      
      return
      end
      
