      subroutine dotem(N,p,s)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),s(mxpart,mxpart)
      integer j,k,N
c---returns 2*piDpj for massless particles
      do j=1,N
      s(j,j)=0d0
      do k=j+1,N
      s(j,k)=2d0*
     & (p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
      s(k,j)=s(j,k)
      enddo
      enddo
      return
      end
