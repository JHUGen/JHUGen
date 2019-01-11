      subroutine scaleset_Msqpt34sq(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt34^2), where M is the mass of the particle (34)
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'breit.f'
      double precision p(mxpart,4),mu0,pttwo

      if((case .eq. 'W_only') .or.
     &   (case .eq. 'Z_only') .or.
     &   (case .eq. 'W_1jet') .or.
     &   (case .eq. 'W_2jet') .or.
     &   (case .eq. 'W_3jet') .or.
     &   (case .eq. 'Z_1jet') .or.
     &   (case .eq. 'Z_2jet') .or.
     &   (case .eq. 'Z_3jet') .or.
     &   (case .eq. 'Wbbbar') .or.
     &   (case .eq. 'Wbbmas') .or.
     &   (case .eq. 'Zbbbar') .or.
     &   (case .eq. 'ggfus0') .or.
     &   (case .eq. 'ggfus1') .or.
     &   (case .eq. 'ggfus2') .or.
     &   (case .eq. 'ggfus3') .or.
     &   (case .eq. 'gagajj') .or.
     &   (case .eq. 'httjet') .or.
     &   (case .eq. 'Higaga') .or.
     &   (case .eq. 'Hgagaj')) then
        mu0=mass3**2+pttwo(3,4,p)**2
        mu0=dsqrt(dabs(mu0))
      else
        write(6,*) 'dynamicscale sqrt(M^2+pt34^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
