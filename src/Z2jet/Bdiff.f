      double complex function Bdiff(s34,s12,msq)
      implicit none
      include 'scale.f'
      include 'first.f'
      double complex qlI2
      double precision s34,s12,msq
      if (first) then
         call qlinit
         first=.false.
      endif
      Bdiff=qlI2(s34,msq,msq,musq,0)-qlI2(s12,msq,msq,musq,0)
      return
      end
