      double complex function F42me(psq,qsq,s,t)
c--- note: ordering of arguments to function is taken from, e.g.
c---       arXiV:0804.4149v3 (App. B) and not hep-ph/0607139 (Eq. 22)
      implicit none
c      include 'epinv.f'
c      include 'scale.f'
      double precision s,t,psq,qsq
c      double precision den
c      double complex qlI4
      double complex F31m,F42meF
c      integer ep
c      den=s*t-psq*qsq
c      F42me=czip
c      do ep=-2,0
c      F42me=F42me
c     . +den*epinv**(-ep)
c     . *qlI4(0d0,psq,0d0,qsq,s,t,0d0,0d0,0d0,0d0,musq,ep)
c      enddo

c--- NOTE: checked on 8/30/09 that this agrees with the expression above
      F42me=2d0*(F31m(s)+F31m(t)-F31m(psq)-F31m(qsq)
     .          +F42meF(psq,qsq,s,t))
     
      return
      end

