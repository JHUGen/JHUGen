      function F41m(psq,s,t)
      implicit none
      include 'types.f'
      complex(dp):: F41m
c--- note: ordering of arguments to function is taken from, e.g.
c---       arXiV:0804.4149v3 (App. B) and not hep-ph/0607139 (Eq. 21)

c      include 'epinv.f'
c      include 'scale.f'
      real(dp):: s,t,psq
c      real(dp):: den
c      complex(dp):: qlI4
      complex(dp):: F31m,F41mF
c      integer:: ep
c      den=s*t
c      F41m=czip
c      do ep=-2,0
c      F41m=F41m+den*epinv**(-ep)
c     & *qlI4(0._dp,0._dp,0._dp,psq,s,t,0._dp,0._dp,0._dp,0._dp,musq,ep)
c      enddo

c--- NOTE: checked on 8/30/09 that this agrees with the expression above
      F41m=2._dp*(F31m(s)+F31m(t)-F31m(psq)+F41mF(psq,s,t))

      return
      end

