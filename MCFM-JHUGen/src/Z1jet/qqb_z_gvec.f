      subroutine qqb_z_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
C***********************************************************************
c     Author: R.K. Ellis                                               *
c     September, 1999.                                                 *
c     Matrix element for Z production                                  *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p5)              *
c     u(-p1)+dbar(-p2)--> g(p5)+ Z^+(l(p3)+a(p4))                      *
C***********************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'nflav.f'
      integer:: j,k,in
C--in is the label of the parton dotted with n
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: z1jetn,fac,p1p2(-1:1,-1:1),n(4)
      complex(dp):: prop

      msq(:,:)=zip

      call dotem(5,p,s)

C-----Protect from photon pole by cutting off at some value about 10 GeV
c      if (s(3,4) < 4._dp*mbsq) return

      fac=16._dp*cf*xn*esq**2*gsq
      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

      p1p2(:,:)=zip

      if (in == 1) then
      p1p2(0,-1)=-aveqg*fac*z1jetn(5,2,1,p,n)
      p1p2(0,+1)=-aveqg*fac*z1jetn(2,5,1,p,n)
      elseif (in == 2) then
      p1p2(+1,0)=-aveqg*fac*z1jetn(1,5,2,p,n)
      p1p2(-1,0)=-aveqg*fac*z1jetn(5,1,2,p,n)
      elseif (in == 5) then
      p1p2(-1,1)=+aveqq*fac*z1jetn(2,1,5,p,n)
      p1p2(1,-1)=+aveqq*fac*z1jetn(1,2,5,p,n)
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=zip
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=+(abs(Q(j)*q1+L(j)*l1*prop)**2
     &              +abs(Q(j)*q1+R(j)*r1*prop)**2)*p1p2(1,-1)
     &             +(abs(Q(j)*q1+L(j)*r1*prop)**2
     &              +abs(Q(j)*q1+R(j)*l1*prop)**2)*p1p2(-1,1)
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=+(abs(Q(k)*q1+L(k)*l1*prop)**2
     &              +abs(Q(k)*q1+R(k)*r1*prop)**2)*p1p2(-1,1)
     &             +(abs(Q(k)*q1+L(k)*r1*prop)**2
     &              +abs(Q(k)*q1+R(k)*l1*prop)**2)*p1p2(1,-1)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=+(abs(Q(j)*q1+L(j)*l1*prop)**2
     &              +abs(Q(j)*q1+R(j)*r1*prop)**2)*p1p2(+1,0)
     &             +(abs(Q(j)*q1+L(j)*r1*prop)**2
     &              +abs(Q(j)*q1+R(j)*l1*prop)**2)*p1p2(-1,0)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=+(abs(Q(-j)*q1+L(-j)*l1*prop)**2
     &              +abs(Q(-j)*q1+R(-j)*r1*prop)**2)*p1p2(-1,0)
     &             +(abs(Q(-j)*q1+L(-j)*r1*prop)**2
     &              +abs(Q(-j)*q1+R(-j)*l1*prop)**2)*p1p2(+1,0)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=+(abs(Q(k)*q1+L(k)*l1*prop)**2
     &              +abs(Q(k)*q1+R(k)*r1*prop)**2)*p1p2(0,+1)
     &             +(abs(Q(k)*q1+L(k)*r1*prop)**2
     &              +abs(Q(k)*q1+R(k)*l1*prop)**2)*p1p2(0,-1)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=+(abs(Q(-k)*q1+L(-k)*l1*prop)**2
     &              +abs(Q(-k)*q1+R(-k)*r1*prop)**2)*p1p2(0,-1)
     &             +(abs(Q(-k)*q1+L(-k)*r1*prop)**2
     &              +abs(Q(-k)*q1+R(-k)*l1*prop)**2)*p1p2(0,+1)
      endif

   19 continue
      enddo
      enddo

      return
      end

      function z1jetn(j1,j2,j5,p,n)
      implicit none
      include 'types.f'
      real(dp):: z1jetn

C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> Z(l(p3)+a(p4))+g(p5)
c   contracted with the vector n(mu)
c   before spin/color average
c---overall factor of 16 gs**2*gw**4*xw**2*CF*xn removed
c--note QED propagator included.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'

      integer:: j1,j2,j3,j4,j5
      real(dp):: n(4),p(mxpart,4),nDn,nDp1,nDp2,nDp3,nDp4
      j3=3
      j4=4

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDp3=n(4)*p(j3,4)-n(3)*p(j3,3)-n(2)*p(j3,2)-n(1)*p(j3,1)
      nDp4=n(4)*p(j4,4)-n(3)*p(j4,3)-n(2)*p(j4,2)-n(1)*p(j4,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)


c      z1jetn=((nDp1*s(j2,j3)/s(j1,j5)-nDp2*s(j1,j4)/s(j2,j5))**2
c     & +two*(s(j2,j3)*nDp1/s(j1,j5)-s(j1,j4)*nDp2/s(j2,j5))*(nDp2+nDp3)
c     & -(s(j1,4)-s(j2,3))**2*s(j3,j4)*nDn/(four*s(j1,j5)*s(j2,j5))
c     & +(nDp2+nDp3)**2)/s(j3,j4)**2

      z1jetn=
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

      z1jetn=z1jetn/s(j3,j4)**2

      return
      end





