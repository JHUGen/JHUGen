      double precision function H4prenorm()
c--- This function returns the appropriate renormalization factor for
c--- the Higgs + 4 parton amplitudes
c--- it includes:    a) strong coupling renormalization
c---                 b) finite renormalization of Hgg effective coupling
c---                 c) finite renormalization of alpha-s in dred scheme
c---
c--- Note that this function returns zero when checking the
c--- (unrenormalized) results in the EGZ paper
      implicit none
      include 'constants.f'
      include 'b0.f'
      include 'epinv.f'
      include 'scheme.f'
      logical CheckEGZ
      common/CheckEGZ/CheckEGZ
!$omp threadprivate(/CheckEGZ/)

      if (CheckEGZ) then
        H4prenorm=0d0
      else
        H4prenorm=(-4d0*b0/xn*epinv+11d0/xn)
        if (scheme .eq. 'dred') then
        H4prenorm=H4prenorm+2d0/3d0
        endif
      endif  
      
      return
      end
