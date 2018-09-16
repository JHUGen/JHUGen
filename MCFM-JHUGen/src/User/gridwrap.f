      subroutine write_grid(xx)
      implicit none
      include 'types.f'
      real(dp):: xx
      write(6,*) 'write_grid:You have to link to gridwrap.cxx'
      write(6,*) 'write_grid:for applgrid'
      write(6,*) 'write_grid:Program stopped'
      stop
      end

      subroutine fill_grid(p)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4)
      write(6,*) 'fill_grid:You have to link to gridwrap.cxx'
      write(6,*) 'fill_grid:for applgrid'
      write(6,*) 'fill_grid:Program stopped'
      stop
      end

c      subroutine book_grid()
c      implicit none
c      include 'types.f'
c      write(6,*) 'book_grid:You have to link to gridwrap.cxx'
c      write(6,*) 'book_grid:for applgrid'
c      write(6,*) 'book_grid:Program stopped'
c      stop
c      end

