      subroutine scaleset_ptphoton(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  pt(photon)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      real(dp):: p(mxpart,4),mu0,pt

      if    ((kcase==kWgamma) .or.
     &       (kcase==kZgamma)) then
        mu0=pt(5,p)
      elseif((kcase==kdirgam)) then
        mu0=pt(3,p)
      else
        write(6,*) 'dynamicscale pt(photon)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
