      subroutine qqb_w_cjet_massless(p,msq)
      implicit none
c----Matrix element for W production
C----averaged over initial colours and spins
C for nwz=+1
c     u(-p1)+dbar(-p2)--> W^+(n(p3)+e^+(p4))   + c(p5)(massless)
C For nwz=-1
c     d(-p1)+ubar(-p2)--> W^-(e^-(p3)+nbar(p4))+ c(p5)(massless) 
c---
      include 'constants.f'
      include 'ckm.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double precision qgWq,qbgWqb,gqbWqb,gqWq,w1jet
     

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(5,p,s)
c---protect from soft and collinear singularities
c      if  ((-s(1,5) .lt. cutoff) .or. (-s(2,5) .lt. cutoff)) return
c---calculate the propagator
      fac=gwsq**2*gsq*V

c      qqbWg= +aveqq*fac*w1jet(1,2,3,4,5)
      gqbWqb=-aveqg*fac*w1jet(5,2,3,4,1)
      qgWq=  -aveqg*fac*w1jet(1,5,3,4,2)
      
c      qbqWg= +aveqq*fac*w1jet(2,1,3,4,5)
      qbgWqb=-aveqg*fac*w1jet(5,1,3,4,2)
      gqWq=  -aveqg*fac*w1jet(2,5,3,4,1)

      do j=-(nf-2),(nf-2)
      do k=-(nf-2),(nf-2)
      if ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=Vsq(j,-4)*qgWq
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=Vsq(j,+4)*qbgWqb
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=Vsq(-4,k)*gqWq
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=Vsq(+4,k)*gqbWqb
      endif

      enddo
      enddo
      return
      end
 
