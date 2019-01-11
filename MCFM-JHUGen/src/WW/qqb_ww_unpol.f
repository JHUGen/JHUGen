      subroutine qqb_ww_unpol(p,msq)
      implicit none
C----Author R.K.Ellis December 2003
c----Matrix element for WW production
c----in the notation of DKS
C----averaged over initial colours and spins
C----massless final state particles
c     q(-p1)+qbar(-p2)-->q'(p5)+bar{q'}(p6)+n(p3)+ebar(p4)
c--- note that non-leptonic W decays do not include scattering diagrams
c--- unpolarized W-pairs
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zerowidth.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'plabel.f'

      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),qdks(mxpart,4)
      double complex prop12,prop34,prop56,cprop,cs_z(2,2),cs_g(2,2)
      double complex AWW(2,-1:1,-1:1)
      double complex a123456(-1:1,-1:1),b123456(-1:1,-1:1)
      double complex a213456(-1:1,-1:1),b213456(-1:1,-1:1)
      double complex a126543(-1:1,-1:1),b126543(-1:1,-1:1)
      double complex a216543(-1:1,-1:1),b216543(-1:1,-1:1)

      double precision fac
      integer j,k,jk,hel34,hel56,n34max,n34min,n56max,n56min
      integer, parameter :: minus=1,mplus=2
      double precision, parameter :: mp(nf)=(/-1d0,+1d0,-1d0,+1d0,-1d0/)
      fac=gw**8*xn*aveqq
C---multiply by factor for c-sbar+u-dbar hadronic decay
      if (plabel(5) .eq. 'qj') fac=2d0*xn*fac
      
      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
      enddo
      enddo

C----Change the momenta to DKS notation 
C   swapped possibility if we want to swap momenta for hadronic case
c   We have --- f(p1) + f'(p2)-->mu^-(p3)+nubar(p4)+e^+(p6)+nu(p5)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)
c----
C   or normal configuration
c   We have --- f(p1) + f'(p2)-->mu^-(p5)+nubar(p6)+e^+(p4)+nu(p3)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)

      if ((plabel(5) .eq. 'qj') .and. (plabel(3) .eq. 'el')) then 
C ---- swapped case
      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(3,j)
      qdks(4,j)=p(4,j)
      qdks(5,j)=p(6,j)
      qdks(6,j)=p(5,j)
      enddo
      else
      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(5,j)
      qdks(4,j)=p(6,j)
      qdks(5,j)=p(4,j)
      qdks(6,j)=p(3,j)
      enddo
      endif


      call spinoru(6,qdks,za,zb)
c--   s returned from sprod (common block) is 2*dot product
      
      if (zerowidth .neqv. .true.) then
      write(6,*) 'Please set zerowidth .eqv. .true.'
      stop
      endif

      prop12=s(1,2)/dcmplx(s(1,2)-zmass**2,zmass*zwidth)
      prop34=s(3,4)/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
      prop56=s(5,6)/dcmplx(s(5,6)-wmass**2,wmass*wwidth)
      cprop=dcmplx(1d0)

      
c-- couplings with or without photon pole
      do j=1,2
      cs_z(minus,j)=+mp(j)*l(j)*sin2w*prop12
      cs_z(mplus,j)=-mp(j)*2d0*Q(j)*xw*prop12
      cs_g(minus,j)=+mp(j)*2d0*Q(j)*xw
      cs_g(mplus,j)=+mp(j)*2d0*Q(j)*xw
      enddo


c---case dbar-d and d-dbar
   
      call susana(1,2,6,5,4,3,za,zb,a126543,b126543)
      call susana(2,1,6,5,4,3,za,zb,a216543,b216543)
      call susana(1,2,3,4,5,6,za,zb,a123456,b123456)
      call susana(2,1,3,4,5,6,za,zb,a213456,b213456)

     
      do j=-nf,nf
      k=-j
c--Exclude gluon-gluon initial state
      if (j.eq.0) go to 20
      jk=max(j,k)

c--assign values
c---Remember that base process is ub-u so this has the natural 123456 order

      n34min=-1
      n34max=+1
      n56min=-1
      n56max=+1
      if (j .gt. 0) then
          if         (tau(jk) .eq. +1d0) then
            do hel34=n34min,n34max
            do hel56=n56min,n56max
            AWW(minus,hel34,hel56)=prop56*prop34*(a213456(hel34,hel56)
     .      +(cs_z(minus,2)+cs_g(minus,2))*b213456(hel34,hel56))
            AWW(mplus,hel34,hel56)=(cs_z(mplus,2)+cs_g(mplus,2))
     .      *b123456(hel34,hel56)*prop56*prop34
            enddo
            enddo
          elseif     (tau(jk) .eq. -1d0) then
            do hel34=n34min,n34max
            do hel56=n56min,n56max
            AWW(minus,hel34,hel56)=prop56*prop34*(a216543(hel34,hel56)
     .      +(cs_z(minus,1)+cs_g(minus,1))*b216543(hel34,hel56))
            AWW(mplus,hel34,hel56)=(cs_z(mplus,1)+cs_g(mplus,1))
     .      *b126543(hel34,hel56)*prop56*prop34
            enddo
            enddo
          endif
      elseif (j .lt. 0) then
          if     (tau(jk) .eq. +1d0) then
            do hel34=n34min,n34max
            do hel56=n56min,n56max
            AWW(minus,hel34,hel56)=prop56*prop34*(a123456(hel34,hel56)
     .      +(cs_z(minus,2)+cs_g(minus,2))*b123456(hel34,hel56))
            AWW(mplus,hel34,hel56)=(cs_z(mplus,2)+cs_g(mplus,2))
     .      *b213456(hel34,hel56)*prop56*prop34
            enddo
            enddo
          elseif (tau(jk) .eq. -1d0) then
            do hel34=n34min,n34max
            do hel56=n56min,n56max
            AWW(minus,hel34,hel56)=prop56*prop34*(a126543(hel34,hel56)
     .      +(cs_z(minus,1)+cs_g(minus,1))*b126543(hel34,hel56)) 
            AWW(mplus,hel34,hel56)=(cs_z(mplus,1)+cs_g(mplus,1))
     .      *b216543(hel34,hel56)*prop56*prop34
            enddo
            enddo
          endif
      endif
     
C-- Inclusion of width for W's a la Baur and Zeppenfeld with cprop.
      do hel34=n34min,n34max
      do hel56=n56min,n56max
      AWW(minus,hel34,hel56)=cprop*AWW(minus,hel34,hel56)
      AWW(mplus,hel34,hel56)=cprop*AWW(mplus,hel34,hel56)

      msq(j,k)=msq(j,k)
     & +fac*(abs(AWW(minus,hel34,hel56))**2
     &      +abs(AWW(mplus,hel34,hel56))**2)
     & /dfloat(n34max-n34min+1)/dfloat(n56max-n56min+1)
      enddo
      enddo
 20   continue

      enddo
      
      return
      end


