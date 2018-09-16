      subroutine qqb_w_cjet_gvec(p,n,in,msq)
      implicit none
      include 'types.f'

c----Matrix element for W production in association with a top quark
C----averaged over initial colours and spins
c    contracted with the vector v(mu)
C For nwz=+1
c     c(-p1)+sbar(-p2)--> g(p5)+ W^+(n(p3)+e^+(p4))
C For nwz=-1
c     s(-p1)+cbar(-p2)--> g(p5)+ W^-(e^-(p3)+nbar(p4))
c---
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'nwz.f'
      include 'ckm.f'
      include 'nflav.f'
      integer:: j,k,in
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: wcjetn,p1p2(-1:1,-1:1),n(4)

      real(dp):: FAC

      do j=-nf,nf
      do k=-nf,nf
        msq(j,k)=0._dp
      enddo
      enddo

      do j=-1,+1
      do k=-1,+1
        p1p2(j,k)=0._dp
      enddo
      enddo

      fac=2._dp*gsq*V*gwsq**2
      call dotem(5,p,s)

      if     (in == 1) then
        if (nwz == +1) p1p2(0,-1)=-aveqg*fac*wcjetn(5,2,1,p,n)
        if (nwz == -1) p1p2(0,+1)=-aveqg*fac*wcjetn(2,5,1,p,n)
      elseif (in == 2) then
        if (nwz == -1) p1p2(+1,0)=-aveqg*fac*wcjetn(1,5,2,p,n)
        if (nwz == +1) p1p2(-1,0)=-aveqg*fac*wcjetn(5,1,2,p,n)
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav
      if     (((j == +1).or.(j == +3)) .and. (k == 0)) then
          msq(j,k)=Vsq(j,-4)*p1p2(+1,0)
      elseif (((j == -1).or.(j == -3)) .and. (k == 0)) then
          msq(j,k)=Vsq(j,+4)*p1p2(-1,0)
      elseif ((j == 0) .and. ((k == +1).or.(k == +3))) then
          msq(j,k)=Vsq(-4,k)*p1p2(0,+1)
      elseif ((j == 0) .and. ((k == -1).or.(k == -3))) then
          msq(j,k)=Vsq(+4,k)*p1p2(0,-1)
      endif

      enddo
      enddo

      return
      end
