      double complex function F31m(s)
      implicit none
      include 'epinv.f'
      include 'scale.f'
      double precision s
c      double complex qlI3
      double complex lnrat
c      integer ep
c      F31m=czip
c      do ep=-2,0
c      F31m=F31m+s*epinv**(-ep)*qlI3(0d0,0d0,s,0d0,0d0,0d0,musq,ep)
c      enddo
      
c--- NOTE: checked on 8/31/09 that this agrees with the expression above
      F31m=epinv**2-epinv*lnrat(-s,musq)+0.5d0*lnrat(-s,musq)**2
      
      return
      end

