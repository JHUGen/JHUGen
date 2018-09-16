      subroutine phi3m0(xth,xphi,p0,p1,p2,wt,*)
      implicit none
      include 'types.f'
c     massive particle p0 in rest frame
c     decaying into p1 fixed mass 0 and p2 fixed mass 0.
c     vectors returned p1 and p2 are in the frame in which
C     p0 is supplied
c result is 1/8/pi * 2|p|/sqrts  * domega/(4*pi)
c     factor of (2*pi)^4 included in definition of phase space
c     Expression evaluated is
c     d^4 p1 d^4 p2 (2 pi)^4 delta(p0-p1-p2)/(2 pi)^6
c     delta(p2^2) delta(p3^2)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'kprocess.f'
      real(dp):: p0(4),p1(4),p2(4),p1cm(4)
      real(dp):: xth,xphi,phi,s,roots,costh,sinth
      real(dp):: wt0,wt
      integer:: j
      parameter(wt0=one/eight/pi)
      wt=0._dp

      s=p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2
      if (s < zip) then
c       if (case(1:5) .ne. 'vlchk') then
!        write(6,*) 's<0 in phi3m0',s
c       endif
       return 1
      endif

      roots=sqrt(s)
      costh=two*xth-one
      sinth=sqrt(one-costh**2)
      phi=twopi*xphi

      wt=wt0

      p1cm(4)=roots/two
      p1cm(1)=roots/two*sinth*sin(phi)
      p1cm(2)=roots/two*sinth*cos(phi)
      p1cm(3)=roots/two*costh

c      write(6,*) 'e',roots/two*(s+m1sq-m2sq)/s
c      write(6,*) 'p',roots/two*lambda/s

c      write(6,*) 'sinth',sinth
c      write(6,*) 'costh',costh
c      write(6,*) 'p1cm**2',p1cm(4)**2-p1cm(1)**2-p1cm(2)**2-p1cm(3)**2
c      pause

      call boost(roots,p0,p1cm,p1)
      do j=1,4
      p2(j)=p0(j)-p1(j)
      enddo

      if (  (p0(4) < 0._dp)
     & .or. (p1(4) < 0._dp)
     & .or. (p2(4) < 0._dp)) then
      write(6,*) 'p0',p0(4),p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2,s
      write(6,*) 'p1',p1(4),p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      write(6,*) 'p2',p2(4),p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      write(6,*) 'in phi3m0'
      endif
      return
      end

