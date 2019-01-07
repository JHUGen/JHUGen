      subroutine spinork(N,p,zabikj,zbaikj,k)
c---Calculate spinor products dotted in with a vector k
c---extended to deal with negative energies ie with all momenta outgoing
C   zabikj=<i-|k|j-> zbaikj=<i+|k|j+> 
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl, 
c---za(i,j)*zb(j,i)=s(i,j)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),rt(mxpart),k(4),kp,km,flip(mxpart)
      double complex pr(mxpart),pl(mxpart),f(mxpart),kr,kl,
     & zabikj(mxpart,mxpart),zbaikj(mxpart,mxpart)
      integer i,j,N
      
C--setup components for vector which is contracted in
      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=dcmplx(+k(3),-k(2))
      kl=dcmplx(+k(3),+k(2))

c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
            zabikj(j,j)=czip
            zbaikj(j,j)=czip
            
C-----positive energy case
            if (p(j,4) .gt. 0d0) then
            flip(j)=1d0
            f(j)=cone
            else
            flip(j)=-1d0
            f(j)=im
            endif
            rt(j)=dsqrt(flip(j)*(p(j,4)+p(j,1)))
            pr(j)=dcmplx(flip(j)*p(j,3),-flip(j)*p(j,2))
            pl(j)=Dconjg(pr(j))
      enddo
      do i=1,N
         do j=1,i
         zabikj(i,j)=f(i)*f(j)
     & *(pr(i)*pl(j)*dcmplx(kp/(rt(i)*rt(j)))
     &    -pr(i)*kl*dcmplx(rt(j)/rt(i))
     &    -dcmplx(rt(i)/rt(j))*kr*pl(j)+dcmplx(rt(i)*rt(j)*km))
         zbaikj(j,i)=zabikj(i,j) 
         zabikj(j,i)=flip(i)*flip(j)*Dconjg(zabikj(i,j))
         zbaikj(i,j)=flip(i)*flip(j)*Dconjg(zbaikj(j,i))

      enddo
      enddo

      return
      end
