      double complex function F41m(psq,s,t)
c--- note: ordering of arguments to function is taken from, e.g.
c---       arXiV:0804.4149v3 (App. B) and not hep-ph/0607139 (Eq. 21)
      implicit none
c      include 'epinv.f'
c      include 'scale.f'
      double precision s,t,psq
c      double precision den
c      double complex qlI4
      double complex F31m,F41mF
c      integer ep
c      den=s*t
c      F41m=czip
c      do ep=-2,0
c      F41m=F41m+den*epinv**(-ep)
c     . *qlI4(0d0,0d0,0d0,psq,s,t,0d0,0d0,0d0,0d0,musq,ep)
c      enddo
      
c--- NOTE: checked on 8/30/09 that this agrees with the expression above
      F41m=2d0*(F31m(s)+F31m(t)-F31m(psq)+F41mF(psq,s,t))
      
      return
      end

