      subroutine write_grid(xx)
      implicit none
      double precision xx
      write(6,*) 'write_grid:You have to link to gridwrap.cxx'
      write(6,*) 'write_grid:for applgrid'
      write(6,*) 'write_grid:Program stopped'
      stop
      end

      subroutine fill_grid(p)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4)
      write(6,*) 'fill_grid:You have to link to gridwrap.cxx'
      write(6,*) 'fill_grid:for applgrid'
      write(6,*) 'fill_grid:Program stopped'
      stop
      end

c      subroutine book_grid()
c      write(6,*) 'book_grid:You have to link to gridwrap.cxx'
c      write(6,*) 'book_grid:for applgrid'
c      write(6,*) 'book_grid:Program stopped'
c      stop
c      end

