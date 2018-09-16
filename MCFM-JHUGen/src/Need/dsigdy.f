      function dsigdy(x)
      implicit none
      include 'types.f'
      include 'nyy.f'
      real(dp)::dsigdy
      real(dp)::sig(nyy),xyy(nyy),err
      logical:: first
      data first/.true./
      save sig,xyy
      if (first) then
      first=.false.
      open(unit=47,file='outw+.dat',status='old')
      do ny=1,nyy
      read(47,*) xyy(ny),sig(ny),err
c      write(6,*) xyy(ny),sig(ny)
      enddo
      close(unit=47)
      endif
      mpot=3
      dsigdy=1d3*ddvdif(sig,xyy,nyy,x,mpot)
      return
      end
