      subroutine getperp(p,i1,i2,v1,v2)
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
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),p1(4),p2(4),v1(4),v2(4),dotpr,p1Dp2,
     . offset,vsq,vDp1,vDp2,v1Dv2,tiny
      integer i1,i2,nu      
      parameter(tiny=1d-15)
      
c--- set up vectors p1,p2 and normalize
      do nu=1,4
        p1(nu)=p(i1,nu)/p(i1,4)
        p2(nu)=p(i2,nu)/p(i2,4)
      enddo
      
      p1Dp2=dotpr(p1,p2)

c--- construct momentum v1 orthognal to p1 and p2
      
c--- offset is to avoid exceptional configurations, usually it is zero
      offset=0d0
      
   22 continue
   
      v1(4)=0d0
      v1(1)=5d0+offset
      v1(2)=3d0-offset
      v1(3)=4d0-offset
      
      vsq=dotpr(v1,v1)
      vDp1=dotpr(v1,p1)
      vDp2=dotpr(v1,p2)

c--- avoid the case where we construct a momentum v1 with v1.v1=0      
      if (abs(vsq-2d0*vDp1*vDp2/p1Dp2) .lt. tiny) then
        offset=offset+0.7d0
        goto 22
      endif
      
      do nu=1,4
        v1(nu)=v1(nu)-(vDp1*p2(nu)+vDp2*p1(nu))/p1Dp2
      enddo

c--- normalize so that v1.v1=-1
      vsq=dotpr(v1,v1)
      do nu=1,4
        v1(nu)=v1(nu)/dsqrt(-vsq)
      enddo

c--- construct momentum v2 orthogonal to p1 and p2
      
c--- offset is to avoid exceptional configurations, usually it is zero
      offset=0d0
      
   33 continue
   
      v2(4)=0d0
      v2(1)=1d0-offset
      v2(2)=2d0+offset
      v2(3)=3d0+offset
      
      vsq=dotpr(v2,v2)
      vDp1=dotpr(v2,p1)
      vDp2=dotpr(v2,p2)

c--- avoid the case where we construct a momentum v2 with v2.v2=0 
      if (abs(vsq-2d0*vDp1*vDp2/p1Dp2) .lt. tiny) then
        offset=offset+0.6d0
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
      if (-vsq .lt. tiny) then
        offset=offset+0.6d0
        goto 33
      endif
      
      do nu=1,4
        v2(nu)=v2(nu)/dsqrt(-vsq)
      enddo

c--- for checking
c      write(6,*) 'v1.v1 = ',dotpr(v1,v1)
c      write(6,*) 'v2.v2 = ',dotpr(v2,v2)
c      write(6,*) 'v1.v2 = ',dotpr(v1,v2)
c      write(6,*) 'v1.p1 = ',dotpr(v1,p1)
c      write(6,*) 'v1.p2 = ',dotpr(v1,p2)
c      write(6,*) 'v2.p1 = ',dotpr(v2,p1)
c      write(6,*) 'v2.p2 = ',dotpr(v2,p2)

      return
      end
      
      
      double precision function dotpr(p,q)
c--- returns the dot product of the two four vectors p(4) and q(4)
      implicit none
      double precision p(4),q(4)
      
      dotpr=p(4)*q(4)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      
      return
      end
      
