      double complex function F42mhF(psq,qsq,s,t)
      implicit none
c      include 'scale.f'
c      include 'epinv.f'
      double precision s,t,psq,qsq
c      double precision den
c      double complex qlI4,qlI3
      double complex Lsm1_2mht
c      integer ep

c      den=s*t
c--- note: added a factor of 1/2 here
c      F42mhF=
c     . +den*qlI4(0d0,0d0,psq,qsq,s,t,0d0,0d0,0d0,0d0,musq,0)/2d0
c     . -  s*qlI3(0d0,0d0,s,0d0,0d0,0d0,musq,0)/2d0
c     . -  t*qlI3(0d0,0d0,t,0d0,0d0,0d0,musq,0)
c     . +psq*qlI3(0d0,0d0,psq,0d0,0d0,0d0,musq,0)/2d0
c     . +qsq*qlI3(0d0,0d0,qsq,0d0,0d0,0d0,musq,0)/2d0

c--- NOTE: checked on 8/30/09 that Lsm1_2mht == (expression above)
      F42mhF=Lsm1_2mht(s,t,psq,qsq)
      
      return
      end

