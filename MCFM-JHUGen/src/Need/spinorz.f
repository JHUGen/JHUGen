      subroutine spinorz(N,p,za,zb)
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double precision p(mxpart,4)
      double complex c12(mxpart),f(mxpart),rt(mxpart)
      integer i,j,N
      
c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=czip
         zb(j,j)=za(j,j)

C-----positive energy case
            rt(j)=dcmplx(p(j,4)+p(j,3))
            if (dble(rt(j)) .eq. 0d0) then
            write(6,*) 'spinorz:j',j
            write(6,*) 'spinorz fails for momenta directed along z axis'
            stop
            endif
            rt(j)=cdsqrt(rt(j))
            c12(j)=dcmplx(p(j,1),p(j,2))
            f(j)=cone
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &              -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)
     &   *(c12(i)*(rt(j)/rt(i))-c12(j)*(rt(i)/rt(j)))

         if (abs(s(i,j)).lt.1d-9) then
         zb(i,j)=-(f(i)*f(j))**2*dconjg(za(i,j))
         else
         zb(i,j)=-dcmplx(s(i,j))/za(i,j)
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

      return
      end
