      subroutine fill_APPLgrid(p)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4)
c--- This is a hook to call the APPLgrid routine that should
c--- write out the grids
c---  e.g. call fill_grid(p) in v1.2.6 of APPLgrid

c      call fill_grid(p)
      write(6,*) 'No APPLgrid write-out routine available!'
      write(6,*) 'You must edit src/User/fill_APPLgrid appropriately.'
      stop

      return
      end

