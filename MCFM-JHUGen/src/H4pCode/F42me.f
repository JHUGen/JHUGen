      function F42me(psq,qsq,s,t)
      implicit none
      include 'types.f'
      complex(dp):: F42me
c--- note: ordering of arguments to function is taken from, e.g.
c---       arXiV:0804.4149v3 (App. B) and not hep-ph/0607139 (Eq. 22)

c      include 'epinv.f'
c      include 'scale.f'
      real(dp):: s,t,psq,qsq
c      real(dp):: den
c      complex(dp):: qlI4
      complex(dp):: F31m,F42meF
c      integer:: ep
c      den=s*t-psq*qsq
c      F42me=czip
c      do ep=-2,0
c      F42me=F42me
c     & +den*epinv**(-ep)
c     & *qlI4(0._dp,psq,0._dp,qsq,s,t,0._dp,0._dp,0._dp,0._dp,musq,ep)
c      enddo

c--- NOTE: checked on 8/30/09 that this agrees with the expression above
      F42me=2._dp*(F31m(s)+F31m(t)-F31m(psq)-F31m(qsq)
     &          +F42meF(psq,qsq,s,t))

      return
      end

