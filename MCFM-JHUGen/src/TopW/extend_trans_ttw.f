      subroutine extend_trans_ttw(pold,p,ptrans,pext)
      implicit none
      include 'types.f'
      
************************************************************************
*     Authors: J.M. Campbell and R.K. Ellis                            *
*     August, 2008.                                                    *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: pold(mxpart,4),p(mxpart,4),ptrans(mxpart,4),
     & pext(mxpart,4),pt(4),ptt(4),ptb(4),pttb(4),
     & p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),
     & p3out(4),p4out(4),p5out(4),p6out(4),p7out(4),p8out(4)
      integer:: nu
      
      do nu=1,4
        pt(nu)=p(5,nu)
        ptt(nu)=ptrans(5,nu)
        p3(nu)=pold(3,nu)
        p4(nu)=pold(4,nu)
        p5(nu)=pold(5,nu)
      enddo
c--- Boost input vector pn to output vector pnout using the same
c--- transformation as required to boost massive vector pt to ptt
      call boostx(p3,pt,ptt,p3out)
      call boostx(p4,pt,ptt,p4out)
      call boostx(p5,pt,ptt,p5out)

      do nu=1,4
        ptb(nu)=p(6,nu)
        pttb(nu)=ptrans(6,nu)
        p6(nu)=pold(6,nu)
        p7(nu)=pold(7,nu)
        p8(nu)=pold(8,nu)
      enddo
      
      call boostx(p6,ptb,pttb,p6out)
      call boostx(p7,ptb,pttb,p7out)
      call boostx(p8,ptb,pttb,p8out)


      do nu=1,4
        pext(1,nu)=ptrans(1,nu)
        pext(2,nu)=ptrans(2,nu)
        pext(3,nu)=p3out(nu)
        pext(4,nu)=p4out(nu)
        pext(5,nu)=p5out(nu)
        pext(6,nu)=p6out(nu)
        pext(7,nu)=p7out(nu)
        pext(8,nu)=p8out(nu)
        pext(9,nu)=ptrans(3,nu)
        pext(10,nu)=ptrans(4,nu)
        pext(11,nu)=ptrans(7,nu)
      enddo
            
      return
      end
      
