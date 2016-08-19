      double precision function Rgen(pjet,i,p,j)
c--- Calculate the angular separation between pjet(i) and p(j)
c--- This routine is a generalization of R.f: Rgen(p,i,p,j) == R(p,i,j) 
      implicit none
      include 'constants.f'
      double precision pjet(mxpart,4),p(mxpart,4),
     & r1,r2,dely,delphi,ei,ej
      integer i,j

      ei=dsqrt(pjet(i,1)**2+pjet(i,2)**2+pjet(i,3)**2)
      ej=dsqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)

      r1= (ei+pjet(i,3))*(ej-p(j,3))/
     .     ((ej+p(j,3))*(ei-pjet(i,3)))
      dely=0.5d0*dlog(r1)

      r2= (pjet(i,1)*p(j,1)+pjet(i,2)*p(j,2))
     .     /dsqrt((pjet(i,1)**2+pjet(i,2)**2)*(p(j,1)**2+p(j,2)**2))
      if (r2 .gt. +0.9999999D0) r2=+1D0
      if (r2 .lt. -0.9999999D0) r2=-1D0
      delphi=dacos(r2)

      Rgen=dsqrt(dely**2+delphi**2)
      
      return
      end
      
