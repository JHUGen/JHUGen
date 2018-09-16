      subroutine scaleset_Msqpt5sq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt5^2), where M is the mass of the particle (34)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'breit.f'
      real(dp):: p(mxpart,4),mu0,pt

      if((kcase==kWgamma) .or.
     &   (kcase==kZgamma) .or.
     &   (kcase==kZgajet) .or.
     &   (kcase==kZga2jt)) then
        mu0=mass3**2+pt(5,p)**2
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale sqrt(M^2+pt5^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      

