      function Li4(x)
      implicit none
      include 'types.f'
      real(dp):: Li4
C     returns Li_4(x) for real x, minf < x < 1

      real(dp):: x
      complex(dp):: wgplg

      if (x > 1._dp) then
        write(6,*) 'x>1 in Li4 function, src/Lib/Li4.f'
        stop
      endif

      Li4 = real(wgplg(3,1,x))
      end
