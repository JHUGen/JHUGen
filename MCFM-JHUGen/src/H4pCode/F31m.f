      function F31m(s)
      implicit none
      include 'types.f'
      complex(dp):: F31m

      include 'epinv.f'
      include 'scale.f'
      real(dp):: s
c      complex(dp):: qlI3
      complex(dp):: lnrat
c      integer:: ep
c      F31m=czip
c      do ep=-2,0
c      F31m=F31m+s*epinv**(-ep)*qlI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,ep)
c      enddo

c--- NOTE: checked on 8/31/09 that this agrees with the expression above
      F31m=epinv**2-epinv*lnrat(-s,musq)+0.5_dp*lnrat(-s,musq)**2

      return
      end

