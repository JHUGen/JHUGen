      subroutine qqb_w_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
c----Matrix element for W production
C----averaged over initial colours and spins
c    contracted with the vector v(mu)
C For nwz=+1
c     u(-p1)+dbar(-p2)--> g(p5)+ W^+(n(p3)+e^+(p4))
C For nwz=-1
c     d(-p1)+ubar(-p2)--> g(p5)+ W^-(e^-(p3)+nbar(p4))
c---
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'ckm.f'
      include 'masses.f'
      integer:: j,k,in
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: w1jetn,p1p2(-1:1,-1:1),n(4)

      real(dp):: fac,prop

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=zip
      enddo
      enddo

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=zip
      enddo
      enddo

      call dotem(5,p,s)
c---calculate the propagator
      prop=((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      fac=two*gsq*V*gwsq**2/prop

      if (in == 1) then
      p1p2(0,-1)=-aveqg*fac*w1jetn(5,2,3,4,1,p,n)
      p1p2(0,+1)=-aveqg*fac*w1jetn(2,5,3,4,1,p,n)
      elseif (in == 2) then
      p1p2(+1,0)=-aveqg*fac*w1jetn(1,5,3,4,2,p,n)
      p1p2(-1,0)=-aveqg*fac*w1jetn(5,1,3,4,2,p,n)
      elseif (in == 5) then
      p1p2(1,-1)=+aveqq*fac*w1jetn(1,2,3,4,5,p,n)
      p1p2(-1,1)=+aveqq*fac*w1jetn(2,1,3,4,5,p,n)
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*p1p2(1,-1)
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*p1p2(-1,1)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*p1p2(+1,0)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*p1p2(-1,0)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*p1p2(0,+1)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*p1p2(0,-1)
      endif

      enddo
      enddo

      return
      end


      function w1jetn(j1,j2,j3,j4,j5,p,n)
      implicit none
      include 'types.f'
      real(dp):: w1jetn

C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> W(l(p3)+a(p4)+g(p5)
c   contracted with the vector n(mu)
c---rewritten so that momentum conservation is not used
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5
      real(dp):: p(mxpart,4),n(4),nDn,
     &                 nDp1,nDp2,nDp3,nDp4

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDp3=n(4)*p(j3,4)-n(3)*p(j3,3)-n(2)*p(j3,2)-n(1)*p(j3,1)
      nDp4=n(4)*p(j4,4)-n(3)*p(j4,3)-n(2)*p(j4,2)-n(1)*p(j4,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)

c      w1jetn=(nDp1*s(j2,j3)/s(j1,j5)-nDp2*s(j1,j4)/s(j2,j5))**2
c     & +(s(j2,j3)*nDp1/s(j1,j5)-s(j1,j4)*nDp2/s(j2,j5))
c     & *(nDp2+nDp3-nDp4-nDp1)
c     & -(s(j1,j4)-s(j2,j3))**2*s(j3,j4)*nDn/four/s(j1,j5)/s(j2,j5)
c     & -(nDp1+nDp4)*(nDp2+nDp3)

      w1jetn=
     &  s(j2,j3)*nDp1**2/s(j1,j5)**2*(s(j1,j4)+s(j4,j5))
     & +s(j1,j4)*nDp2**2/s(j2,j5)**2*(s(j2,j3)+s(j3,j5))

     & +(s(j1,j2)*s(j3,j5)*s(j4,j5)*nDn/four
     & -nDp1*nDp2*(two*s(j1,j4)*s(j2,j3)+s(j1,j4)*s(j3,j5)
     & +s(j2,j3)*s(j4,j5)+s(j3,j5)*s(j4,j5)))/s(j1,j5)/s(j2,j5)

     & +(s(j1,j4)*nDp1*nDp3+s(j4,j5)*nDp1*nDp3-s(j2,j3)*nDp1*nDp4
     & -(s(j1,j3)+s(j2,j3))*s(j4,j5)*nDn/four)/s(j1,j5)

     & +(s(j2,j3)*nDp2*nDp4+s(j3,j5)*nDp2*nDp4-s(j1,j4)*nDp2*nDp3
     & -(s(j1,j4)+s(j2,j4))*s(j3,j5)*nDn/four)/s(j2,j5)
     & +s(j3,j4)*nDn/four-nDp3*nDp4


      return
      end

