      real*8 function dsigdy(x)      
      implicit real*8 (a-h,o-z)
      include 'nyy.f'
      real*8 sig(nyy),xyy(nyy),err
      logical first
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
