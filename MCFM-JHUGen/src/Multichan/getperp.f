      subroutine getperp(p,i1,i2,v1,v2,*)
      implicit none
      include 'types.f'
************************************************************************
*                                                                      *
*     Given two momenta, p(i1,..) and p(i2, ..), construct two         *
*     vectors v1 and v2 that satisfy:                                  *
*                                                                      *
*      v1.p1 = 0, v1.p2 = 0                                            *
*      v2.p1 = 0, v2.p2 = 0                                            *
*      v1.v2 = 0                                                       *
*      v1.v1 = -1, v2.v2 = -1                                          *
*                                                                      *
*     Author: J.M.Campbell                                             *
*       Date: 19th March 2009                                          *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),p1(4),p2(4),v1(4),v2(4),dotpr,p1Dp2,
     & offset,vsq,vDp1,vDp2,v1Dv2,tiny
      integer:: i1,i2,nu
      parameter(tiny=1.e-15_dp)

c      p(i1,1)=   0.0000000000000000._dp
c      p(i1,2)=  0.0000000000000000._dp
c      p(i1,3)=624.49264065711168._dp
c      p(i1,4)= -624.49264065711168._dp
c      p(i2,1)= -8.42527279039759410e-006_dp
c      p(i2,2)= -1.62606796446384395e-006_dp
c      p(i2,3)= -421.23910530482470._dp
c      p(i2,4)=  421.23910530482482._dp

c--- set up vectors p1,p2 and normalize
      do nu=1,4
        p1(nu)=p(i1,nu)/p(i1,4)
        p2(nu)=p(i2,nu)/p(i2,4)
      enddo

      p1Dp2=dotpr(p1,p2)

c--- construct momentum v1 orthogonal to p1 and p2

c--- offset is to avoid exceptional configurations, usually it is zero
      offset=0._dp

   22 continue

      v1(4)=0._dp
      v1(1)=5._dp+offset
      v1(2)=3._dp-offset
      v1(3)=4._dp-offset

      vsq=dotpr(v1,v1)
      vDp1=dotpr(v1,p1)
      vDp2=dotpr(v1,p2)

c--- avoid the case where we construct a momentum v1 with v1.v1=0
      if (abs(vsq-two*vDp1*vDp2/p1Dp2) < tiny) then
        offset=offset+0.7_dp
        goto 22
      endif

      do nu=1,4
        v1(nu)=v1(nu)-(vDp1*p2(nu)+vDp2*p1(nu))/p1Dp2
      enddo

c--- normalize so that v1.v1=-1
      vsq=dotpr(v1,v1)
      do nu=1,4
        v1(nu)=v1(nu)/sqrt(-vsq)
      enddo

c--- construct momentum v2 orthogonal to p1 and p2

c--- offset is to avoid exceptional configurations, usually it is zero
      offset=0._dp

   33 continue

      v2(4)=0._dp
      v2(1)=1._dp-offset
      v2(2)=2._dp+offset
      v2(3)=3._dp+offset

      vsq=dotpr(v2,v2)
      vDp1=dotpr(v2,p1)
      vDp2=dotpr(v2,p2)

c--- avoid the case where we construct a momentum v2 with v2.v2=0
      if (abs(vsq-two*vDp1*vDp2/p1Dp2) < tiny) then
        offset=offset+0.6_dp
        goto 33
      endif

      do nu=1,4
        v2(nu)=v2(nu)-(vDp1*p2(nu)+vDp2*p1(nu))/p1Dp2
      enddo

c--- now ensure that v1.v2=0
      v1Dv2=dotpr(v1,v2)
      vsq=dotpr(v1,v1)
      do nu=1,4
        v2(nu)=v2(nu)-v1(nu)*v1Dv2/vsq
      enddo

c--- normalize so that v2.v2=-1
      vsq=dotpr(v2,v2)
c--- avoid the case where we construct a momentum v2 with v2.v2=0
c--- ( -vsq should be positive)
      if (-vsq < tiny) then
        offset=offset+0.6_dp
        goto 33
      endif

      do nu=1,4
        v2(nu)=v2(nu)/sqrt(-vsq)
      enddo

c--- for checking
c      write(6,*) 'v1.v1 = ',dotpr(v1,v1)
c      write(6,*) 'v2.v2 = ',dotpr(v2,v2)
c      write(6,*) 'v1.v2 = ',dotpr(v1,v2)
c      write(6,*) 'v1.p1 = ',dotpr(v1,p1)
c      write(6,*) 'v1.p2 = ',dotpr(v1,p2)
c      write(6,*) 'v2.p1 = ',dotpr(v2,p1)
c      write(6,*) 'v2.p2 = ',dotpr(v2,p2)

      if ((v1(4) .ne. v1(4)) .or. (v2(4) .ne. v2(4))) then
!        write(6,*) 'warning: could not generate momenta in getperp.f'
        return 1
c        write(6,*) 'i1,i2',i1,i2
c        write(6,*) 'p(i1,:)',p(i1,:)
c        write(6,*) 'p(i2,:)',p(i2,:)
c        write(6,*) 'p1',p1
c        write(6,*) 'p2',p2
c        write(6,*) 'p1Dp2',p1Dp2
c        write(6,*) 'vDp1',vDp1
c        write(6,*) 'vDp2',vDp2
c        write(6,*) 'v1',v1
c        write(6,*) 'v2',v2
c        write(6,*) 'offset',offset
c        stop
      endif

      return
      end


      function dotpr(p,q)
      implicit none
      include 'types.f'
      real(dp):: dotpr
c--- returns the dot product of the two four vectors p(4) and q(4)

      real(dp):: p(4),q(4)

      dotpr=p(4)*q(4)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)

      return
      end

