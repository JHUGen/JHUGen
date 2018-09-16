      subroutine phase4m(r,p1,p2,p3,p4,p5,p6,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'debug.f'
c---- generate phase space for 2-->4 process
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^8)
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p12(4),p56(4),p456(4),s3min
      real(dp):: wt,wt3456,wt456,wt56,wt0
      integer:: j,iflag
      parameter(wt0=1._dp/twopi**2)
      data iflag/0/
      save iflag 
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      s3min=mt
      if (iflag == 1) then
      iflag=0
      call phi1_2m(mt,r(1),r(2),r(3),s3min,p12,p3,p456,wt3456,*99)
      call phi1_2m(mt,r(4),r(5),r(6),s3min,p456,p4,p56,wt456,*99)
      else
      iflag=1
      call phi1_2m(mt,r(1),r(2),r(3),s3min,p12,p4,p456,wt3456,*99)
      call phi1_2m(mt,r(4),r(5),r(6),s3min,p456,p3,p56,wt456,*99)
      endif 

c p56 is the b-bbar system
C--decay 56 into massless b quarks
      call phi3m0(r(7),r(8),p56,p5,p6,wt56,*99)
      wt=wt0*wt3456*wt456*wt56
      if (debug) write(6,*) 'wt in phase4',wt
      return
 99   wt=0._dp
      return 1
      end

