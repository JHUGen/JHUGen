      subroutine scaleset_ptphoton(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c---  pt(photon)
      implicit none
      include 'constants.f'
      include 'process.f'
      double precision p(mxpart,4),mu0,pt

      if    ((case .eq. 'Wgamma') .or.
     &       (case .eq. 'Zgamma')) then
        mu0=pt(5,p)
      elseif((case .eq. 'dirgam')) then
        mu0=pt(3,p)
      else
        write(6,*) 'dynamicscale pt(photon)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
