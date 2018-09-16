      double precision function Li4(x)
C     returns Li_4(x) for real x, minf < x < 1
      implicit none
      double precision x
      double complex wgplg
      
      if (x .gt. 1d0) then
        write(6,*) 'x>1 in Li4 function, src/Lib/Li4.f'
        stop
      endif
      
      Li4 = dble(wgplg(3,1,x))
      end
