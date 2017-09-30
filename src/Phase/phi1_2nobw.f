      subroutine phi1_2nobw(x1,x2,x3,x4,p1,p2,p3,wt,*)
c     massive particle p1 decaying into p2 mass m2 and p3 mass m3.
c     with invariant mass of particle two s2 and particle three s3 integrated over.
c
c     Modified version of phi1_2 such that neither p2 nor p3 decays using BW
c
c     vectors returned p2 and p3 are in the same frame as p1 is supplied
c     Expression evaluated is 
c     ds2 ds3 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6 
c     delta(p2^2-s2) delta(p3^2-s3)
c
      implicit none
      include 'constants.f'
      include 'heavyflav.f'
      include 'masses.f'
      include 'process.f'
      include 'zerowidth.f'
      include 'verbose.f'
      include 'breit.f'
      include 'limits.f'
      include 'dm_params.f' 
      include 'first.f' 
      double precision p1(4),p2(4),p3(4),p3cm(4)
      double precision x1,x2,x3,x4,costh,sinth,phi,cphi,sphi
      double precision wt,wt0,w2,w3
      double precision s2max,s2min,s3max,s3min
      double precision m1,m2,s1,s2,s3,lambda
      integer j
      logical oldzerowidth
      common/lambda/lambda,s1,s2,s3
      parameter(wt0=one/8d0/pi)
!$omp threadprivate(/lambda/)

      if (verbose) then
        if(first) then
        write(6,*) 'phase space using phi1_2nobw'
        first=.false.
        endif
      endif

      wt=0d0
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      if (s1 .lt. 0d0) return 1
      m1=dsqrt(s1)

      if (m1 .lt. 1d-8) return 1 ! for safety

      s2min=1d-8
      s2max=s1

c--- adjust limits to match input file if it is more restrictive
      s2min=max(s2min,m3456min**2)
      s2max=min(s2max,m3456max**2)
      if (s2max .lt. s2min)  return 1 ! for safety
      
      w2=s2max-s2min
      s2=s2max*x1+s2min*(1d0-x1)
      
      m2=dsqrt(s2)
      s3min=1d-8
      s3max=(m2-m1)**2

      if (s3max .lt. s3min) return 1 ! for safety

      w3=s3max-s3min
      s3=s3max*x2+s3min*(1d0-x2)

      costh=two*x3-one      
      phi=twopi*x4
      sinth=dsqrt(one-costh**2)
      cphi=dcos(phi)
      sphi=dsin(phi)
      lambda=((s1-s2-s3)**2-4d0*s2*s3)

      if (lambda .lt. 0d0) then
        return 1
      endif

      lambda=dsqrt(lambda)
      wt=wt0*w2*w3*lambda/s1

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
       if (case(1:5) .ne. 'vlchk') then 
        write(6,*) '   m1=',m1
        write(6,*) 's2min=',s2min
        write(6,*) 's2max=',s2max
        write(6,*) 's3min=',s3min
        write(6,*) 's3max=',s3max
        write(6,*) 'p1',p1(4),p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2,s1
        write(6,*) 'p2',p2(4),p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2,s2
        write(6,*) 'p3',p3(4),p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2,s3
        write(6,*) 'n2,n3',n2,n3
        write(6,*) 'in phi1_2bw.f'
       endif
       return 1
      endif

      return
      end



