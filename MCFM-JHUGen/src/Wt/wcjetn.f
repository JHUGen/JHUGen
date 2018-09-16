      function wcjetn(p1,p2,p5,p,n)
      implicit none
      include 'types.f'
      real(dp):: wcjetn

C---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> W(l(p3)+a(p4)+g(p5)
c   contracted with the vector n(mu)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p5
      real(dp):: p(mxpart,4),n(4),nDn,prop,
     &                 nDp1,nDp2,nDp3,nDp4
      p3=3
      p4=4
      nDp1=n(4)*p(p1,4)-n(3)*p(p1,3)-n(2)*p(p1,2)-n(1)*p(p1,1)
      nDp2=n(4)*p(p2,4)-n(3)*p(p2,3)-n(2)*p(p2,2)-n(1)*p(p2,1)
      nDp3=n(4)*p(p3,4)-n(3)*p(p3,3)-n(2)*p(p3,2)-n(1)*p(p3,1)
      nDp4=n(4)*p(p4,4)-n(3)*p(p4,3)-n(2)*p(p4,2)-n(1)*p(p4,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,p5)

c---calculate the propagator
      prop=((s(p3,p4)-wmass**2)**2+(wmass*wwidth)**2)

      wcjetn=-nDp3*nDp4+0.25_dp*s(p3,p4)*nDn
     & +s(p1,p4)*nDp2**2*(s(p3,p2)+s(p3,p5))/s(p2,p5)**2
     & +s(p3,p2)*nDp1**2*(s(p1,p4)+s(p4,p5))/s(p1,p5)**2

     & -(nDp2*(nDp3*s(p1,p4)-(s(p3,p2)+s(p3,p5))*nDp4)
     & +s(p3,p5)*0.25_dp*nDn*(s(p4,p2)+s(p1,p4)))/s(p2,p5)

     & -(nDp1*(s(p3,p2)*nDp4-(s(p1,p4)+s(p4,p5))*nDp3)
     & +s(p4,p5)*0.25_dp*nDn*(s(p3,p2)+s(p3,p1)))/s(p1,p5)

     & -(nDp1*nDp2*(s(p3,p2)*s(p4,p5)+s(p3,p5)*s(p1,p4)
     & +s(p3,p5)*s(p4,p5)+2_dp*s(p3,p2)*s(p1,p4))
     & -s(p3,p5)*s(p4,p5)*s(p1,p2)*0.25_dp*nDn)/s(p1,p5)/s(p2,p5)

c--- massless expression
c      wcjetn=(nDp1*s(p2,p3)/s(p1,p5)-nDp2*s(p1,p4)/s(p2,p5))**2
c     & +(s(p2,p3)*nDp1/s(p1,p5)-s(p1,p4)*nDp2/s(p2,p5))
c     & *(nDp2+nDp3-nDp4-nDp1)
c     & -(s(p1,p4)-s(p2,p3))**2*s(p3,p4)*nDn/4_dp/s(p1,p5)/s(p2,p5)
c     & -(nDp1+nDp4)*(nDp2+nDp3)

      wcjetn=wcjetn/prop

      return
      end
