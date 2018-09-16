      subroutine scaleset_Msqpt5sq(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt5^2), where M is the mass of the particle (34)
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'breit.f'
      double precision p(mxpart,4),mu0,pt

      if((case .eq. 'Wgamma') .or.
     &   (case .eq. 'Zgamma') .or.
     &   (case .eq. 'Zgajet') .or.
     &   (case .eq. 'Zga2jt')) then
        mu0=mass3**2+pt(5,p)**2
        mu0=dsqrt(dabs(mu0))
      else
        write(6,*) 'dynamicscale sqrt(M^2+pt5^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      

