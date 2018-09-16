      subroutine phi1_2bis(x1,x2,x3,x4,p1,p2,p3,wt,*)
      implicit none
      include 'types.f'
c     massive particle p1 decaying into p2 mass m2 and p3 mass m3.
c     with invariant mass
c     of particle two s2 and particle three s3 integrated over.
c     vectors returned p2 and p3 are in the same frame as p1 is supplied
c     Expression evaluate is
c     ds2 ds3 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6
c     delta(p2^2-s2) delta(p3^2-s3)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'heavyflav.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'zerowidth.f'
      include 'verbose.f'
      include 'breit.f'
      include 'first.f'
      real(dp):: p1(4),p2(4),p3(4),p3cm(4)
      real(dp):: x1,x2,x3,x4,costh,sinth,phi,cphi,sphi
      real(dp):: wt,w2,w3
      real(dp):: s2max,s2min,s3max,s3min
      real(dp):: m1,m2,s1,s2,s3,lambda
      real(dp), parameter::wt0=one/8._dp/pi
      integer:: j
      logical:: oldzerowidth
      common/lambda/lambda,s1,s2,s3
!$omp threadprivate(/lambda/)

      if (verbose) then
      if(first) then
c      if (n2 == 1) write(6,*) 'generating phase space with bw,n2=',n2
c      if (n3 == 1) write(6,*) 'generating phase space with bw,n3=',n3
      first=.false.
      endif
      endif

      wt=0._dp
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      if (s1 < 0._dp) return 1
      m1=sqrt(s1)

c--- if both particles are produced on-shell, reject if m1 too small
      if (
     & zerowidth
     & .and. (m1 < mass2*real(n2,dp)+mass3*real(n3,dp))
     &    ) return 1

c--- top is on-shell for W+t processes, so reject if m1 too small
      if ( ((kcase==kW_twdk) .or. (kcase==kWtdkay)
     &  .or.(kcase==kW_cwdk) .or. (kcase==kWtbwdk)
     &  .or.(kcase==kqq_tth) .or. (kcase==ktth_ww)
     &  .or.(kcase==kqq_ttz) .or. (kcase==kqqtthz)
     & .or. (kcase==kqq_ttw) .or. (kcase==kttwldk))
     & .and. (m1 < mass2) ) return 1
c      s2min=bbsqmin
c      s2max=min(s1,bbsqmax)
      s2min=1.e-15_dp
      s2max=s1
      if (((kcase==kWbbmas) .and. (flav == 5))
     ..or.((kcase==kZbbmas) .and. (flav == 5))
     ..or.(kcase==kZccmas) .or. (kcase==kvlchkm)
     ..or.(kcase==kWbbjet) .or. (kcase==kWbbjem)
     ..or.(kcase==kW_bjet) ) then
        s2min=4._dp*mb**2
      elseif ((kcase==kWbbmas) .and. (flav == 4)) then
        s2min=4._dp*mc**2
      elseif (((kcase==kZbbmas) .and. (flav == 6))
     &    .or. (kcase==kqq_ttz) .or. (kcase==kqqtthz)) then
        s2min=4._dp*mt**2
      elseif (kcase==kW_cjet) then
        s2min=mc**2
      elseif (kcase==kW_tndk)  then
        s2min=mt**2
      elseif (kcase==kWtbndk) then
        s2min=(mt+mb)**2
      elseif ((kcase==kqq_tbg) .or. (kcase==kqqtbgg)) then
        s2min=mt**2
      elseif ((kcase==ktt_ldk) .or. (kcase==ktt_hdk)
     &   .or. (kcase==ktt_udk) .or. (kcase==ktthWdk)
     &   .or. (kcase==kttdkay) .or. (kcase==ktdecay)) then
        s2min=mb**2
      elseif ((kcase==k4ftwdk) .or. (kcase==k4ftjet)
     &   .or. (kcase==kdk_4ft)) then
        s2min=(mt+mb)**2
      elseif ((kcase==kW_twdk) .or. (kcase==kWtdkay)
     &   .or. (kcase==kW_cwdk) .or. (kcase==kWtbwdk)
     &   .or. (kcase==kqq_tth)
     &   .or. (kcase==kqq_ttz) .or. (kcase==kqqtthz)
     &   .or. (kcase==kqq_ttw) .or. (kcase==ktth_ww)) then
        oldzerowidth=zerowidth
        zerowidth=.true.
      endif
      if (s2min > s2max) return 1
      if (n2 == 0) then
         w2=s2max-s2min
         s2=s2max*x1+s2min*(1._dp-x1)
      elseif (n2 == 1) then
         call breitw(x1,s2min,s2max,mass2,width2,s2,w2)
      endif

      if ((kcase==kW_twdk) .or. (kcase==kWtdkay)
     & .or.(kcase==kW_cwdk) .or. (kcase==kWtbwdk)
     & .or.(kcase==kqq_tth)
     & .or. (kcase==kqq_ttz) .or. (kcase==kqqtthz)
     & .or. (kcase==kqq_ttw).or. (kcase==ktth_ww))  then
        zerowidth=oldzerowidth
      endif

      m2=sqrt(s2)
      s3min=1.e-15_dp
      if ((kcase==kqq_tbg) .or. (kcase==kqqtbgg)) s3min=mb**2
      if ((kcase==kqq_tth)
     &  .or. (kcase==kqq_ttw)
     &  .or. (kcase==ktth_ww)
     &  .or. (kcase==kqq_ttz) .or. (kcase==kqqtthz))
     &  s3min=4._dp*mb**2
c      s3min=mb**2 ! DEBUG: hack for s36 small
c      s3min=mt**2 ! DEBUG: hack for s46 small
      s3max=(m2-m1)**2
      if (s3max < s3min) return 1 ! for safety
c      if (s3max-s3min < 1.e-9_dp) return 1
      if (n3 == 0) then
         w3=s3max-s3min
         s3=s3max*x2+s3min*(1._dp-x2)
      elseif (n3 == 1) then
         call breitw(x2,s3min,s3max,mass3,width3,s3,w3)
      endif

      costh=two*x3-one
      phi=twopi*x4
      sinth=sqrt(one-costh**2)
      cphi=cos(phi)
      sphi=sin(phi)
      lambda=((s1-s2-s3)**2-4._dp*s2*s3)

      if (lambda < 0._dp) then
c      write(6,*) '(lambda < 0) in phi1_2.f',lambda
c      write(6,*) 'sqrt(s1)',sqrt(s1)
c      write(6,*) 'sqrt(s2)',sqrt(s2)
c      write(6,*) 'sqrt(s3)',sqrt(s3)
c      write(6,*) s3min,s3,s3max,m1,m2,sqrt(s1),sqrt(s2)
      return 1
      endif
      lambda=sqrt(lambda)
      wt=wt0*w2*w3*lambda/s1


      p3cm(4)=m1/two*(s1+s3-s2)/s1
      p3cm(1)=m1/two*lambda/s1*sinth*sphi
      p3cm(2)=m1/two*lambda/s1*sinth*cphi
      p3cm(3)=m1/two*lambda/s1*costh
      call boost(m1,p1,p3cm,p3)
      do j=1,4
      p2(j)=p1(j)-p3(j)
      enddo
      if (  (p1(4) < 0._dp)
     & .or. (p2(4) < 0._dp)
     & .or. (p3(4) < 0._dp)) then
c       if (case(1:5) .ne. 'vlchk') then
        write(6,*) '   m1=',m1
        write(6,*) 's2min=',s2min
        write(6,*) 's2max=',s2max
        write(6,*) 's3min=',s3min
        write(6,*) 's3max=',s3max
        write(6,*) 'p1',p1(4),p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2,s1
        write(6,*) 'p2',p2(4),p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2,s2
        write(6,*) 'p3',p3(4),p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2,s3
        write(6,*) 'n2,n3',n2,n3
        write(6,*) 'in phi1_2.f'
c       endif
       return 1
      endif

      return
      end



