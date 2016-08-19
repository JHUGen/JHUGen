      double complex function F42meF(psq,qsq,s,t)
c--- note: ordering of arguments to function is taken from, e.g.
c---       arXiV:0804.4149v3 (App. B) and not hep-ph/0607139 (Eq. 22)
      implicit none
c      include 'scale.f'
      double precision s,t,psq,qsq
c      double precision den
c      double complex qlI4,qlI3
      double complex Lsm1_2me
c      den=s*t-psq*qsq
c--- note: added a factor of 1/2 here
c      F42meF=
c     . +den*qlI4(0d0,psq,0d0,qsq,s,t,0d0,0d0,0d0,0d0,musq,0)/2d0
c     . -  s*qlI3(0d0,0d0,s,0d0,0d0,0d0,musq,0)
c     . -  t*qlI3(0d0,0d0,t,0d0,0d0,0d0,musq,0)
c     . +psq*qlI3(0d0,0d0,psq,0d0,0d0,0d0,musq,0)
c     . +qsq*qlI3(0d0,0d0,qsq,0d0,0d0,0d0,musq,0)
      
c--- NOTE: checked on 8/30/09 that Lsm1_2me == (expression above)
      F42meF=Lsm1_2me(s,t,psq,qsq)
     
      return
      end

