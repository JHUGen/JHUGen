      subroutine phi1_2m_bw(m2,x3,xth,xphi,s3min,p1,p2,p3,
     . bwmass,bwwidth,wt,*)
c     This routine differs from phi1_2m in that s3 is always generated
c     according to a Breit-Wigner described by the additional
c     arguments bwmass and bwwidth
c-----------------------------------------------------------------------
c     massive particle p1 decaying into p2 mass m2 and p3 mass-squared s3.
c     with invariant mass of particle three s3 integrated over.
c     s3min is the minimum value of s3.
c     Vectors returned p2 and p3 are in the same frame as p1 is supplied.
c     Expression evaluated is
c     ds3 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6
c     delta(p2^2-m2) delta(p3^2-s3)
      implicit none
      include 'constants.f'
      include 'debug.f'
      double precision p1(4),p2(4),p3(4),p3cm(4)
      double precision x3,xth,xphi,costh,sinth,phi,cphi,sphi
      double precision wt,wt0,w3
      double precision s3max,s3min,xx,xexp
      double precision m1,m2,m3,s1,s2,s3,lambda,bwmass,bwwidth
      integer j
      parameter(wt0=one/8d0/pi)

      wt=0d0
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      if (s1 .lt. 0d0) return 1
      m1=dsqrt(s1)
      s2=m2**2
      s3max=(m2-m1)**2
      if (s3min .gt. s3max) return 1

      xx=1d0
      call breitw(x3,s3min,s3max,bwmass,bwwidth,s3,w3)

      m3=dsqrt(s3)

      costh=two*xth-one
      phi=twopi*xphi
      sinth=dsqrt(one-costh**2)
      cphi=dcos(phi)
      sphi=dsin(phi)
      lambda=((s1-s2-s3)**2-4d0*s2*s3)

      if ((lambda .lt. 0d0) .or. debug) then
      write(6,*) 'lambda in phi1_2m_bw',lambda
      write(6,*) 's1 in phi1_2m_bw',s1
      write(6,*) 's2 in phi1_2m_bw',s2
      write(6,*) 's3 in phi1_2m_bw',s3
      write(6,*) 'm1 in phi1_2m_bw',m1
      write(6,*) 'm2 in phi1_2m_bw',m2
      write(6,*) 'm3 in phi1_2m_bw',m3
      write(6,*) 'm1-m2-m3 in phi1_2m_bw',m1-m2-m3
      write(6,*) 'xx in phi1_2m_bw',xx
      write(6,*) 'x3 in phi1_2m_bw',x3
      write(6,*) 'bwmass in phi1_2m',bwmass
      write(6,*) 'bwwidth in phi1_2m',bwwidth
      return 1
      endif
      lambda=dsqrt(lambda)

      wt=wt0*w3*lambda/s1

      if(debug) write(6,*) 'wt in phi1_2m_bw',wt

      p3cm(4)=m1/two*(s1+s3-s2)/s1
      p3cm(1)=m1/two*lambda/s1*sinth*sphi
      p3cm(2)=m1/two*lambda/s1*sinth*cphi
      p3cm(3)=m1/two*lambda/s1*costh


      call boost(m1,p1,p3cm,p3)
      do j=1,4
      p2(j)=p1(j)-p3(j)
      enddo


      if (  (p1(4) .lt. 0d0)
     & .or. (p2(4) .lt. 0d0)
     & .or. (p3(4) .lt. 0d0)) then
      return 1
      endif
      return
      end



