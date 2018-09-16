      function r(p,i,j)
      implicit none
      include 'types.f'
      real(dp):: r
c----calculate the jets separation between p(i) and p(j)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp):: p(mxpart,4),r1,r2,dely,delphi,ei,ej,pti2,ptj2
      integer:: i,j

      pti2=p(i,1)**2+p(i,2)**2
      ptj2=p(j,1)**2+p(j,2)**2

      ei=sqrt(pti2+p(i,3)**2)
      ej=sqrt(ptj2+p(j,3)**2)

      r1= (ei+p(i,3))*(ej-p(j,3))/
     &   ((ej+p(j,3))*(ei-p(i,3)))
      dely=0.5_dp*log(r1)

      r2= (p(i,1)*p(j,1)+p(i,2)*p(j,2))/sqrt(pti2*ptj2)
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)

      r=sqrt(dely**2+delphi**2)

      return
      end

