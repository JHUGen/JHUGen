      function Bdiff(s34,s12,msq)
      implicit none
      include 'types.f'
      complex(dp):: Bdiff

      include 'scale.f'
      include 'first.f'
      complex(dp):: qlI2
      real(dp):: s34,s12,msq
      if (first) then
         call qlinit
         first=.false.
      endif
      Bdiff=qlI2(s34,msq,msq,musq,0)-qlI2(s12,msq,msq,musq,0)
      return
      end
