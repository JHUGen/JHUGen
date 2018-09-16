      subroutine gen5mdk(r,p,pswt,*)
      implicit none
      include 'types.f'
c--- this routine is an extension of gen4 to include the decay
c--- of one of the heavy particles


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'masses.f'
      include 'zerowidth.f'
      include 'limits.f'
      include 'kprocess.f'
      real(dp):: r(mxdim)
      real(dp):: p(mxpart,4),pswt,smin
      real(dp):: p1(4),p2(4),p8(4),p9(4),p7(4),p6(4)
      real(dp):: p567(4),p56(4),p5(4),p3(4),p4(4)
      real(dp):: wt567,wt56
      integer:: nu
      real(dp):: wt0
      parameter(wt0=1._dp/twopi**2)

      if ((kcase==kZ_tdkj) .or. (kcase==kH_tdkj)
     &.or.(kcase==kZtdk2j) ) then
      call gen5(r,p,pswt,*99)

      p1(:)=p(1,:)
      p2(:)=p(2,:)
      p3(:)=p(3,:)
      p4(:)=p(4,:)
      p567(:)=p(5,:)
      p8(:)=p(6,:)
      p9(:)=p(7,:)

c--- set up minimum invariant mass for the W in the top decay
      if (zerowidth) then
        smin=wmass**2
      else
        smin=bbsqmin
      endif

c--- decay top -> b W
      call phi1_2m_bw(mb,r(14),r(15),r(16),smin,p567,p7,p56,
     & wmass,wwidth,wt567,*99)
c--- decay W -> e n
      call phi3m0(r(17),r(18),p56,p5,p6,wt56,*99)

c--- compute new weight
      pswt=wt0*pswt*wt567*wt56*pi*mt*twidth

      else

      write(6,*) 'Case not foreseen in gen5mdk.f'
      stop

      endif

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=p7(nu)
      p(8,nu)=p8(nu)
      p(9,nu)=p9(nu)
      enddo

      return

   99 continue
      return 1

      end

