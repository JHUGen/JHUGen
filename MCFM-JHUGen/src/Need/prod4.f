      subroutine prod4(p,za,zb,s)
c  extended to deal with negative energies ie with all momenta outgoing
c---all particle assumed massless
      implicit none
      double precision p(6,4),s(6,6)
      double complex za(6,6),zb(6,6),c23(6),f(6),one,im
      double precision rt(6)
      integer i,j
      parameter(one=(1d0,0d0),im=(0d0,1d0))
      
      do j=1,6
         za(j,j)=(0D0,0D0)
         zb(j,j)=za(j,j)
         s(j,j)=0d0
         if (p(j,4) .gt. 0.d0) then
         rt(j)=dsqrt(p(j,4)-p(j,1))
         c23(j)=dcmplx(p(j,2),-p(j,3))
         f(j)=one
         else
         rt(j)=dsqrt(-p(j,4)+p(j,1))
         c23(j)=dcmplx(-p(j,2),p(j,3))
         f(j)=im
         endif
      enddo
       do i=2,6
         do j=1,i-1
         s(i,j)=2d0*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &   -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)*(c23(i)*rt(j)/rt(i)-c23(j)*rt(i)/rt(j))
         za(j,i)=-za(i,j)
         zb(i,j)=s(i,j)/za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
       enddo
      enddo

      return
      end

      subroutine prod2(p,za,zb)
c  extended to deal with negative energies ie with all momenta outgoing
      implicit none
      double precision p(4,4),dot(4,4)
      double complex za(4,4),zb(4,4),c23(4),f(4),one,im
      double precision rt(4)
      integer i,j
      parameter(one=(1.d0,0.d0),im=(0.d0,1.d0))
      
      do j=1,4
         za(j,j)=(0.D0,0.D0)
         zb(j,j)=za(j,j)
         if (p(j,4) .gt. 0.d0) then
         rt(j)=dsqrt(p(j,4)-p(j,1))
         c23(j)=dcmplx(p(j,2),-p(j,3))
         f(j)=one
         else
         rt(j)=dsqrt(-p(j,4)+p(j,1))
         c23(j)=dcmplx(-p(j,2),p(j,3))
         f(j)=im
         endif
      enddo
       do i=2,4
         do j=1,i-1
         dot(i,j)=2d0*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &   -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)*(c23(i)*rt(j)/rt(i)-c23(j)*rt(i)/rt(j))
         za(j,i)=-za(i,j)
         zb(i,j)=dot(i,j)/za(i,j)
         zb(j,i)=-zb(i,j)
       enddo
      enddo
      return
      end


      subroutine prod3(p,za,zb,dot)
c  extended to deal with negative energies ie with all momenta outgoing
      implicit none
      double precision p(5,4),dot(5,5)
      double complex za(5,5),zb(5,5),c23(5),f(5),one,im
      double precision rt(5)
      integer i,j
      parameter(one=(1.d0,0.d0),im=(0.d0,1.d0))
      do j=1,5
         za(j,j)=(0.D0,0.D0)
         zb(j,j)=za(j,j)
         if (p(j,4) .gt. 0.d0) then
         rt(j)=dsqrt(p(j,4)-p(j,1))
         c23(j)=dcmplx(p(j,2),-p(j,3))
         f(j)=one
         else
         rt(j)=dsqrt(-p(j,4)+p(j,1))
         c23(j)=dcmplx(-p(j,2),p(j,3))
         f(j)=im
         endif
      enddo
       do i=2,5
         do j=1,i-1
         dot(i,j)=2d0*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &   -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)*(c23(i)*rt(j)/rt(i)-c23(j)*rt(i)/rt(j))
         za(j,i)=-za(i,j)
         zb(i,j)=dot(i,j)/za(i,j)
         zb(j,i)=-zb(i,j)
       enddo
      enddo
      return
      end

