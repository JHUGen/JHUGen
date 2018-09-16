C  (C) Copr. 1986-92 Numerical Recipes Software ]2w.1,r1..

      FUNCTION pvdpythag(a,b)
      implicit none
      DOUBLE PRECISION a,b,pvdpythag
      DOUBLE PRECISION absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pvdpythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          pvdpythag=0.0d0
        else
          pvdpythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END
