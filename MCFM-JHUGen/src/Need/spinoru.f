      subroutine spinoru(N,p,za,zb)
      implicit none
      include 'types.f'
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),rt(mxpart)
      complex(dp):: c23(mxpart),f(mxpart)
      integer:: i,j,N
      
c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=czip
         zb(j,j)=za(j,j)

C-----positive energy case
         if (p(j,4) > zip) then
            rt(j)=sqrt(p(j,4)+p(j,1))
            c23(j)=cplx2(p(j,3),-p(j,2))
            f(j)=cone
         else
C-----negative energy case
            rt(j)=sqrt(-p(j,4)-p(j,1))
            c23(j)=cplx2(-p(j,3),p(j,2))
            f(j)=im
         endif
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &              -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)
     &   *(c23(i)*rt(j)/rt(i)-c23(j)*rt(i)/rt(j))

         if (abs(s(i,j))<1e-5_dp) then
         zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
         else
         zb(i,j)=-s(i,j)/za(i,j)
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

      return
      end
