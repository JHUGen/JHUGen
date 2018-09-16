      subroutine qqb_w_tndk_gvec(p,n,in,msq)
      implicit none
      include 'types.f'

c----Matrix element for W production in association with a top quark
C----averaged over initial colours and spins
c    contracted with the vector v(mu)
C For nwz=+1
c     t(-p1)+bbar(-p2)--> g(p5)+ W^+(n(p3)+e^+(p4))
C For nwz=-1
c     b(-p1)+tbar(-p2)--> g(p5)+ W^-(e^-(p3)+nbar(p4))
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
      integer:: j,k,in
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: wcjetn,p1p2(-1:1,-1:1),n(4)

      real(dp):: FAC

      do j=-nf,nf
      do k=-nf,nf
        msq(j,k)=zero
      enddo
      enddo

      do j=-1,+1
      do k=-1,+1
        p1p2(j,k)=zero
      enddo
      enddo

      fac=two*gsq*V*gwsq**2
      call dotem(5,p,s)

      if     (in == 1) then
        if (nwz == +1) p1p2(0,-1)=-aveqg*fac*wcjetn(5,2,1,p,n)
        if (nwz == -1) p1p2(0,+1)=-aveqg*fac*wcjetn(2,5,1,p,n)
      elseif (in == 2) then
        if (nwz == -1) p1p2(+1,0)=-aveqg*fac*wcjetn(1,5,2,p,n)
        if (nwz == +1) p1p2(-1,0)=-aveqg*fac*wcjetn(5,1,2,p,n)
      endif

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j == +5) .and. (k == 0) .and. (nwz == -1)) then
          msq(j,k)=p1p2(+1,0)
      elseif ((j == -5) .and. (k == 0) .and. (nwz == +1)) then
          msq(j,k)=p1p2(-1,0)
      elseif ((j == 0) .and. (k == +5) .and. (nwz == -1)) then
          msq(j,k)=p1p2(0,+1)
      elseif ((j == 0) .and. (k == -5) .and. (nwz == +1)) then
          msq(j,k)=p1p2(0,-1)
      endif

      enddo
      enddo

      return
      end
