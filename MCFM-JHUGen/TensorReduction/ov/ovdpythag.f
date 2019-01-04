C  (C) Copr. 1986-92 Numerical Recipes Software ]2w.1,r1..

      FUNCTION ovdpythag(a,b)
      implicit none
      DOUBLE PRECISION a,b,ovdpythag
      DOUBLE PRECISION absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        ovdpythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          ovdpythag=0.0d0
        else
          ovdpythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END
