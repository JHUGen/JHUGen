      subroutine spinorcurr(N,p,za,zb,zab,zba)
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double precision p(mxpart,4),rt(mxpart)
      double complex c23(mxpart),f(mxpart),fac,
     & zab(mxpart,4,mxpart),zba(mxpart,4,mxpart)
      integer i,j,N
      
c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=czip
         zb(j,j)=za(j,j)
         zab(j,:,j)=dcmplx(2d0*p(j,:))
         zba(j,:,j)=zab(j,:,j)
C-----positive energy case
         if (p(j,4) .gt. 0d0) then
            rt(j)=dsqrt(p(j,4)+p(j,1))
            c23(j)=dcmplx(p(j,3),-p(j,2))
            f(j)=cone
         else
C-----negative energy case
            rt(j)=dsqrt(-p(j,4)-p(j,1))
            c23(j)=dcmplx(-p(j,3),p(j,2))
            f(j)=im
         endif
      enddo

      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)
     &              -p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)
     &   *(c23(i)*dcmplx(rt(j)/rt(i))-c23(j)*dcmplx(rt(i)/rt(j)))

      fac=f(i)*f(j)
C      elseif(h .eq. +1) then
      zba(i,3,j)=fac
     . *(dcmplx(rt(i)/rt(j))*c23(j)+dcmplx(rt(j)/rt(i))*dconjg(c23(i)))
      zba(i,2,j)=fac*im
     . *(dcmplx(rt(i)/rt(j))*c23(j)-dcmplx(rt(j)/rt(i))*dconjg(c23(i)))
      zba(i,1,j)=fac
     . *(dcmplx(rt(i)*rt(j))-dconjg(c23(i))*c23(j)/dcmplx(rt(i)*rt(j)))
      zba(i,4,j)=fac
     . *(dcmplx(rt(i)*rt(j))+dconjg(c23(i))*c23(j)/dcmplx(rt(i)*rt(j)))

C      elseif(h .eq. -1) then
      zab(i,3,j)=+fac
     . *(dcmplx(rt(i)/rt(j))*dconjg(c23(j))+dcmplx(rt(j)/rt(i))*c23(i))
      zab(i,2,j)=fac*im
     . *(dcmplx(rt(j)/rt(i))*c23(i)-dcmplx(rt(i)/rt(j))*dconjg(c23(j)))
      zab(i,1,j)=fac
     . *(dcmplx(rt(i)*rt(j))-c23(i)*dconjg(c23(j))/dcmplx(rt(i)*rt(j)))
      zab(i,4,j)=fac
     . *(dcmplx(rt(i)*rt(j))+c23(i)*dconjg(c23(j))/dcmplx(rt(i)*rt(j)))

      if (abs(s(i,j)).lt.1d-5) then
      zb(i,j)=-(f(i)*f(j))**2*dconjg(za(i,j))
      else
      zb(i,j)=-dcmplx(s(i,j))/za(i,j)
      endif

      za(j,i)=-za(i,j)
      zb(j,i)=-zb(i,j)
      s(j,i)=s(i,j)
      zab(j,:,i)=zba(i,:,j)
      zba(j,:,i)=zab(i,:,j)
      enddo
      enddo
      return
      end
