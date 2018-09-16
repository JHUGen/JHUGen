      subroutine scaleset_Msqpt34sq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt34^2), where M is the mass of the particle (34)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'breit.f'
      real(dp):: p(mxpart,4),mu0,pttwo

      if((kcase==kW_only) .or.
     &   (kcase==kZ_only) .or.
     &   (kcase==kW_1jet) .or.
     &   (kcase==kW_2jet) .or.
     &   (kcase==kW_3jet) .or.
     &   (kcase==kZ_1jet) .or.
     &   (kcase==kZ_2jet) .or.
     &   (kcase==kZ_3jet) .or.
     &   (kcase==kWbbbar) .or.
     &   (kcase==kWbbmas) .or.
     &   (kcase==kZbbbar) .or.
     &   (kcase==kggfus0) .or.
     &   (kcase==kggfus1) .or.
     &   (kcase==kggfus2) .or.
     &   (kcase==kggfus3) .or.
     &   (kcase==kgagajj) .or.
     &   (kcase==khttjet) .or.
     &   (kcase==kHigaga) .or.
     &   (kcase==kHgagaj)) then
        mu0=mass3**2+pttwo(3,4,p)**2
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale sqrt(M^2+pt34^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
