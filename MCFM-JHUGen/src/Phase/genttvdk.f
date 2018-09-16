      subroutine genttvdk(r,p,pswt,*)
      implicit none
      include 'types.f'
c--- this routine is an extension of gen4 to include the decay
c--- of the top quarks for ttV processes
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'mxdim.f'
      include 'breit.f'
      include 'masses.f'
      include 'limits.f'
      include 'zerowidth.f'
      integer:: nu
      real(dp):: r(mxdim)
      real(dp):: p(mxpart,4),pswt,
     & p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),p9(4),p10(4),
     & p34(4),p345(4),p78(4),p678(4),
     & smin,wt34,wt345,wt78,wt678

      call gen4(r,p,pswt,*99)

      p1(:)=p(1,:)
      p2(:)=p(2,:)
      p9(:)=p(3,:)
      p10(:)=p(4,:)
      p345(:)=p(5,:)
      p678(:)=p(6,:)

c--- set up minimum invariant mass for the W in the top decay
      if (zerowidth) then
        smin=wmass**2
      else
        smin=bbsqmin
      endif

c--- decay top -> b W
      call phi1_2m_bw(mb,r(11),r(12),r(13),smin,p345,p5,p34,
     & wmass,wwidth,wt345,*99)
c--- decay W -> e n
      call phi3m0(r(14),r(15),p34,p3,p4,wt34,*99)
      pswt=pswt/twopi**2*wt345*wt34*pi*mt*twidth

c--- decay antitop -> b W
      call phi1_2m_bw(mb,r(16),r(17),r(18),smin,p678,p6,p78,
     & wmass,wwidth,wt678,*99)
c--- decay W -> e n
      call phi3m0(r(19),r(20),p78,p7,p8,wt78,*99)
      pswt=pswt/twopi**2*wt678*wt78*pi*mt*twidth

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
      p(10,nu)=p10(nu)
      enddo

      return

   99 continue
      return 1

      end


