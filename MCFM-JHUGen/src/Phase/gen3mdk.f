      subroutine gen3mdk(r,p,m3,m4,m5,pswt,*)
      implicit none
      include 'types.f'
c--- this routine is an extension of gen3m to include the decay
c--- of one of the heavy particles


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'masses.f'
      include 'zerowidth.f'
      include 'limits.f'
      include 'breit.f'
      real(dp):: r(mxdim)
      real(dp):: p(mxpart,4),m3,m4,m5,pswt,wtbw,wten,smin
      real(dp):: pt(4),pb(4),pw(4),pe(4),pn(4),pbbar(4),pq(4)
      integer:: nu

c--- these must be set this way for gen3m
      mass3=m4 ! =mb
      n3=0

c--- first call gen3m, which uses r(1) ... r(7)
      call gen3m(r,p,m3,m4,m5,pswt,*999)
c--- this returns p3 -> p345, p4 -> p6, p5 -> p7
      do nu=1,4
      pt(nu)=p(3,nu)
      pbbar(nu)=p(4,nu)
      pq(nu)=p(5,nu)
      enddo

c--- these must be set this way for phi1_2m
      mass3=wmass
      n3=1

c--- set up minimum invariant mass for the W in the top decay
      if (zerowidth) then
        smin=wmass**2
      else
        smin=wsqmin
      endif

c--- decay top -> b W
      call phi1_2m(mb,r(8),r(9),r(10),smin,pt,pb,pw,wtbw,*999)
c--- decay W -> e n
      call phi3m0(r(11),r(12),pw,pe,pn,wten,*999)

c--- compute new weight
      pswt=pswt/twopi**2*wtbw*wten*pi*mt*twidth

c--- put vectors in their proper places
      do nu=1,4
      p(3,nu)=pn(nu)
      p(4,nu)=pe(nu)
      p(5,nu)=pb(nu)
      p(6,nu)=pbbar(nu)
      p(7,nu)=pq(nu)
      enddo

      return

  999 continue
      pswt=0._dp
      return 1

      end

