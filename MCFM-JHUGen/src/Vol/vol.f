      function vol(s,n)
      implicit none
      include 'types.f'
      real(dp):: vol
c     volume of massless phase space including all (2*pi)'s
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: s
      integer:: n,fac
      vol=twopi**(4-3*n)*(pi/two)**(n-1)*s**(n-2)
      vol=vol/real(fac(n-1)*fac(n-2),dp)
c      if (n == 6) vol=vol/(5._dp*4._dp*3._dp*2._dp)/(4._dp*3._dp*2._dp)
c      if (n == 5) vol=vol/(4._dp*3._dp*2._dp)    /(3._dp*2._dp)
c      if (n == 4) vol=vol/(3._dp*2._dp)        /(2._dp)
c      if (n == 3) vol=vol/(2._dp)            /(1._dp)
      return
      end
      
      function fac(n)
       implicit none
      include 'types.f'
      integer:: fac
      integer:: j,n
      j=1
      fac=1
 10   continue
      if (j <= n) then
      fac=fac*j
      j=j+1
      else
      return
      endif
      go to 10
      end
