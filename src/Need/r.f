      double precision function r(p,i,j)
c----calculate the jets separation between p(i) and p(j)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),r1,r2,dely,delphi,ei,ej,pti2,ptj2
      integer i,j

      pti2=p(i,1)**2+p(i,2)**2
      ptj2=p(j,1)**2+p(j,2)**2

      ei=dsqrt(pti2+p(i,3)**2)
      ej=dsqrt(ptj2+p(j,3)**2)

      r1= (ei+p(i,3))*(ej-p(j,3))/
     &   ((ej+p(j,3))*(ei-p(i,3)))
      dely=0.5d0*dlog(r1)

      r2= (p(i,1)*p(j,1)+p(i,2)*p(j,2))/dsqrt(pti2*ptj2)
      if (r2 .gt. +0.9999999D0) r2=+1D0
      if (r2 .lt. -0.9999999D0) r2=-1D0
      delphi=dacos(r2)

      r=dsqrt(dely**2+delphi**2)
      
      return
      end
      
