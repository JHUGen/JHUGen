      double precision function Li3(x)
      implicit none
C     returns Li_3(x) for real x, minf < x < 1
      double precision x
      double complex wgplg
      
      if (x .gt. 1d0) then
        write(6,*) 'x>1 in Li3 function, src/Lib/Li3.f'
        stop
      endif
      
      Li3 = dble(wgplg(2,1,x))
      end
