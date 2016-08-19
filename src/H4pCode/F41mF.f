      double complex function F41mF(psq,s,t)
c--- note: ordering of arguments to function is taken from, e.g.
c---       arXiV:0804.4149v3 (App. B) and not hep-ph/0607139 (Eq. 21)
      implicit none
c      include 'scale.f'
      double precision s,t,psq
c      double precision den
c      double complex qlI4,qlI3
      double complex Lsm1

c      den=s*t
c--- MODIFIED: added factor of 1/2 and removed finite pieces from poles
c      F41mF=
c     . +den*qlI4(0d0,0d0,0d0,psq,s,t,0d0,0d0,0d0,0d0,musq,0)/2d0
c     . -   s*qlI3(0d0,0d0,s,0d0,0d0,0d0,musq,0)
c     . -   t*qlI3(0d0,0d0,t,0d0,0d0,0d0,musq,0)
c     . + psq*qlI3(0d0,0d0,psq,0d0,0d0,0d0,musq,0)
     
c--- NOTE: checked on 8/30/09 that Lsm1 == (expression above)
      F41mF=Lsm1(-s,-psq,-t,-psq)
      
      return
      end

