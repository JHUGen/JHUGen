      function Rgen(pjet,i,p,j)
      implicit none
      include 'types.f'
      real(dp):: Rgen
c--- Calculate the angular separation between pjet(i) and p(j)
c--- This routine is a generalization of R.f: Rgen(p,i,p,j) == R(p,i,j)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: pjet(mxpart,4),p(mxpart,4),
     & r1,r2,dely,delphi,ei,ej
      integer:: i,j

      ei=sqrt(pjet(i,1)**2+pjet(i,2)**2+pjet(i,3)**2)
      ej=sqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)

      r1= (ei+pjet(i,3))*(ej-p(j,3))/
     &     ((ej+p(j,3))*(ei-pjet(i,3)))
      dely=0.5_dp*log(r1)

      r2= (pjet(i,1)*p(j,1)+pjet(i,2)*p(j,2))
     &     /sqrt((pjet(i,1)**2+pjet(i,2)**2)*(p(j,1)**2+p(j,2)**2))
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)

      Rgen=sqrt(dely**2+delphi**2)

      return
      end

