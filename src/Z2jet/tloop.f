      double complex function tloop(s23,mtsq)
      implicit none
      include 'scale.f'
      include 'first.f'
      double precision s23,mtsq
      double complex qlI3,Bdiff

      if (first) then
      call qlinit
      first=.false.
      endif

      tloop=-1d0-2d0*mtsq/s23*(
     & 3d0*s23*qlI3(0d0,0d0,s23,mtsq,mtsq,mtsq,musq,0)
     & +6d0*Bdiff(s23,0d0,mtsq)+12d0)
      return
      end
