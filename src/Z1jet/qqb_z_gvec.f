      subroutine qqb_z_gvec(p,n,in,msq)
C*********************************************************************** 
c     Author: R.K. Ellis                                               *
c     September, 1999.                                                 *
c     Matrix element for Z production                                  *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p5)              *
c     u(-p1)+dbar(-p2)--> g(p5)+ Z^+(l(p3)+a(p4))                      *
C*********************************************************************** 
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'nflav.f'
      integer j,k,in
C--in is the label of the parton dotted with n
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision z1jetn,fac,p1p2(-1:1,-1:1),n(4)
      double complex prop


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(5,p,s)

C-----Protect from photon pole by cutting off at some value about 10 GeV
c      if (s(3,4) .lt. 4d0*mbsq) return

      fac=16d0*cf*xn*esq**2*gsq
      prop=s(3,4)/dcmplx((s(3,4)-zmass**2),zmass*zwidth)

      do j=-1,1
      do k=-1,1
      p1p2(j,k)=0d0
      enddo
      enddo
       
      if (in .eq. 1) then
      p1p2(0,-1)=-aveqg*fac*z1jetn(5,2,1,p,n)
      p1p2(0,+1)=-aveqg*fac*z1jetn(2,5,1,p,n)
      elseif (in .eq. 2) then
      p1p2(+1,0)=-aveqg*fac*z1jetn(1,5,2,p,n)
      p1p2(-1,0)=-aveqg*fac*z1jetn(5,1,2,p,n)
      elseif (in .eq. 5) then      
      p1p2(-1,1)=+aveqq*fac*z1jetn(2,1,5,p,n)
      p1p2(1,-1)=+aveqq*fac*z1jetn(1,2,5,p,n)
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=0d0
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=+(cdabs(Q(j)*q1+L(j)*l1*prop)**2
     .              +cdabs(Q(j)*q1+R(j)*r1*prop)**2)*p1p2(1,-1)
     .             +(cdabs(Q(j)*q1+L(j)*r1*prop)**2
     .              +cdabs(Q(j)*q1+R(j)*l1*prop)**2)*p1p2(-1,1)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=+(cdabs(Q(k)*q1+L(k)*l1*prop)**2
     .              +cdabs(Q(k)*q1+R(k)*r1*prop)**2)*p1p2(-1,1)
     .             +(cdabs(Q(k)*q1+L(k)*r1*prop)**2
     .              +cdabs(Q(k)*q1+R(k)*l1*prop)**2)*p1p2(1,-1)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+(cdabs(Q(j)*q1+L(j)*l1*prop)**2
     .              +cdabs(Q(j)*q1+R(j)*r1*prop)**2)*p1p2(+1,0)
     .             +(cdabs(Q(j)*q1+L(j)*r1*prop)**2
     .              +cdabs(Q(j)*q1+R(j)*l1*prop)**2)*p1p2(-1,0)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=+(cdabs(Q(-j)*q1+L(-j)*l1*prop)**2
     .              +cdabs(Q(-j)*q1+R(-j)*r1*prop)**2)*p1p2(-1,0)
     .             +(cdabs(Q(-j)*q1+L(-j)*r1*prop)**2
     .              +cdabs(Q(-j)*q1+R(-j)*l1*prop)**2)*p1p2(+1,0)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=+(cdabs(Q(k)*q1+L(k)*l1*prop)**2
     .              +cdabs(Q(k)*q1+R(k)*r1*prop)**2)*p1p2(0,+1)
     .             +(cdabs(Q(k)*q1+L(k)*r1*prop)**2
     .              +cdabs(Q(k)*q1+R(k)*l1*prop)**2)*p1p2(0,-1)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=+(cdabs(Q(-k)*q1+L(-k)*l1*prop)**2
     .              +cdabs(Q(-k)*q1+R(-k)*r1*prop)**2)*p1p2(0,-1)
     .             +(cdabs(Q(-k)*q1+L(-k)*r1*prop)**2
     .              +cdabs(Q(-k)*q1+R(-k)*l1*prop)**2)*p1p2(0,+1)
      endif

   19 continue
      enddo
      enddo

      return
      end
 
      double precision function z1jetn(j1,j2,j5,p,n)
      implicit none 
C---calculates the amplitude squared for the process 
c   q(p1)+qbar(p2) --> Z(l(p3)+a(p4))+g(p5)
c   contracted with the vector n(mu) 
c   before spin/color average
c---overall factor of 16 gs**2*gw**4*xw**2*CF*xn removed
c--note QED propagator included.
      include 'constants.f'
      include 'sprods_com.f'

      integer j1,j2,j3,j4,j5
      double precision n(4),p(mxpart,4),nDn,nDp1,nDp2,nDp3
      j3=3
      j4=4

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDp3=n(4)*p(j3,4)-n(3)*p(j3,3)-n(2)*p(j3,2)-n(1)*p(j3,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)

      z1jetn=((nDp1*s(j2,j3)/s(j1,j5)-nDp2*s(j1,j4)/s(j2,j5))**2
     . +two*(s(j2,j3)*nDp1/s(j1,j5)-s(j1,j4)*nDp2/s(j2,j5))*(nDp2+nDp3)
     . -(s(j1,4)-s(j2,3))**2*s(j3,j4)*nDn/(four*s(j1,j5)*s(j2,j5))
     . +(nDp2+nDp3)**2)/s(j3,j4)**2

      return
      end




