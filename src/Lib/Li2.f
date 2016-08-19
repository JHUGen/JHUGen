      double precision function Li2(x)
      implicit none
C     returns Li_2(x) for real x, minf < x < 1
      double precision x
      double complex wgplg
      
      if (x .gt. 1d0) then
        write(6,*) 'x>1 in Li2 function, src/Lib/Li2.f'
        stop
      endif
      
      Li2 = dble(wgplg(1,1,x))
      end
