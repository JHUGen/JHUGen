      subroutine qqb_w_tndk_gvec(p,n,in,msq)
      implicit none
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
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'nwz.f'
      integer j,k,in
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision wcjetn,p1p2(-1:1,-1:1),n(4)
     
      double precision FAC

      do j=-nf,nf
      do k=-nf,nf
        msq(j,k)=0d0
      enddo
      enddo

      do j=-1,+1
      do k=-1,+1
        p1p2(j,k)=0d0
      enddo
      enddo

      fac=2d0*gsq*V*gwsq**2
      call dotem(5,p,s)

      if     (in .eq. 1) then
        if (nwz .eq. +1) p1p2(0,-1)=-aveqg*fac*wcjetn(5,2,1,p,n)
        if (nwz .eq. -1) p1p2(0,+1)=-aveqg*fac*wcjetn(2,5,1,p,n)
      elseif (in .eq. 2) then
        if (nwz .eq. -1) p1p2(+1,0)=-aveqg*fac*wcjetn(1,5,2,p,n)
        if (nwz .eq. +1) p1p2(-1,0)=-aveqg*fac*wcjetn(5,1,2,p,n)
      endif

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      if     ((j .eq. +5) .and. (k .eq. 0) .and. (nwz .eq. -1)) then
          msq(j,k)=p1p2(+1,0)
      elseif ((j .eq. -5) .and. (k .eq. 0) .and. (nwz .eq. +1)) then
          msq(j,k)=p1p2(-1,0)
      elseif ((j .eq. 0) .and. (k .eq. +5) .and. (nwz .eq. -1)) then
          msq(j,k)=p1p2(0,+1)
      elseif ((j .eq. 0) .and. (k .eq. -5) .and. (nwz .eq. +1)) then
          msq(j,k)=p1p2(0,-1)
      endif

      enddo
      enddo
      
      return
      end
