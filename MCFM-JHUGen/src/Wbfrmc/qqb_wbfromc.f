      subroutine qqb_wbfromc(p,msq)
      implicit none
c----Matrix element for W production
C----averaged over initial colours and spins
C for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + bbar(p5)
C For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ b(p5) 
c---
      include 'constants.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'nflav.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double precision qgWq,qbgWqb,gqbWqb,gqWq,w1cjet

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
C--setup s products through common block
      call dotem(5,p,s)

      fac=gwsq**2*gsq*V

      qgWq=0d0
      gqWq=0d0
      qbgWqb=0d0
      gqbWqb=0d0
      if (nwz .eq. 1) then
      qgWq=  -aveqg*fac*w1cjet(1,5,3,4,2)
      gqWq=  -aveqg*fac*w1cjet(2,5,3,4,1)
      elseif (nwz .eq. -1) then
      gqbWqb=-aveqg*fac*w1cjet(2,5,4,3,1)
      qbgWqb=-aveqg*fac*w1cjet(1,5,4,3,2)
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav
      if (((j .eq. 2) .or. (j .eq. 4)) .and. (k .eq. 0)) then
          msq(j,k)=Vsq(j,-5)*qgWq
      elseif ((j .eq. 0) .and. ((k .eq. +2).or.(k .eq. +4))) then
          msq(j,k)=Vsq(-5,k)*gqWq
      elseif (((j .eq. -2).or.(j .eq. -4)) .and. (k .eq. 0))then
          msq(j,k)=Vsq(j,+5)*qbgWqb
      elseif ((j .eq. 0) .and. ((k .eq. -2).or.(k .eq. -4))) then
          msq(j,k)=Vsq(+5,k)*gqbWqb
      endif
      enddo
      enddo
      return
      end
 


