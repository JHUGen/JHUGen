      function Li2(x)
      implicit none
      include 'types.f'
      real(dp):: Li2

C     returns Li_2(x) for real x, minf < x < 1
      real(dp):: x
      complex(dp):: wgplg

      if (x > 1._dp) then
        write(6,*) 'x>1 in Li2 function, src/Lib/Li2.f'
        stop
      endif

      Li2 = real(wgplg(1,1,x))
      end
