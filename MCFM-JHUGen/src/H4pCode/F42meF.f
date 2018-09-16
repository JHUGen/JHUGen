      function F42meF(psq,qsq,s,t)
      implicit none
      include 'types.f'
      complex(dp):: F42meF
c--- note: ordering of arguments to function is taken from, e.g.
c---       arXiV:0804.4149v3 (App. B) and not hep-ph/0607139 (Eq. 22)

c      include 'scale.f'
      real(dp):: s,t,psq,qsq
c      real(dp):: den
c      complex(dp):: qlI4,qlI3
      complex(dp):: Lsm1_2me
c      den=s*t-psq*qsq
c--- note: added a factor of 1/2 here
c      F42meF=
c     & +den*qlI4(0._dp,psq,0._dp,qsq,s,t,0._dp,0._dp,0._dp,0._dp,musq,0)/2._dp
c     & -  s*qlI3(0._dp,0._dp,s,0._dp,0._dp,0._dp,musq,0)
c     & -  t*qlI3(0._dp,0._dp,t,0._dp,0._dp,0._dp,musq,0)
c     & +psq*qlI3(0._dp,0._dp,psq,0._dp,0._dp,0._dp,musq,0)
c     & +qsq*qlI3(0._dp,0._dp,qsq,0._dp,0._dp,0._dp,musq,0)

c--- NOTE: checked on 8/30/09 that Lsm1_2me == (expression above)
      F42meF=Lsm1_2me(s,t,psq,qsq)

      return
      end

