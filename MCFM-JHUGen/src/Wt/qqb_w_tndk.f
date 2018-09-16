      subroutine qqb_w_tndk(p,msq)
      implicit none
      include 'types.f'

c----Matrix element for W production
C----averaged over initial colours and spins
C for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4))   + tbar(p5)
C For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ t(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'nwz.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      real(dp):: qgWq,qbgWqb,gqbWqb,gqWq,w1cjet

      msq(:,:)=zero

C--setup s products through common block
      call dotem(5,p,s)

      fac=gwsq**2*gsq*V

      if     (nwz == -1) then
c---- basic process is g+b -> W- + t
        qgWq=  -aveqg*fac*w1cjet(1,5,3,4,2)
        gqWq=  -aveqg*fac*w1cjet(2,5,3,4,1)
      elseif (nwz == +1) then
c---- basic process is g+b~ -> W+ + t~
        gqbWqb=-aveqg*fac*w1cjet(2,5,4,3,1)
        qbgWqb=-aveqg*fac*w1cjet(1,5,4,3,2)
      else
        write(6,*) 'Problem with nwz in qqb_w_tndk.f: nwz=',nwz
        stop
      endif

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j == +5) .and. (k == 0) .and. (nwz == -1)) then
          msq(j,k)=qgWq
      elseif ((j == -5) .and. (k == 0) .and. (nwz == +1)) then
          msq(j,k)=qbgWqb
      elseif ((j == 0) .and. (k == +5) .and. (nwz == -1)) then
          msq(j,k)=gqWq
      elseif ((j == 0) .and. (k == -5) .and. (nwz == +1)) then
          msq(j,k)=gqbWqb
      endif

      enddo
      enddo
      return
      end

