      subroutine extend_trans_wt(pold,p,ptrans,pext)
      implicit none
c--- take vector 
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: pold(mxpart,4),p(mxpart,4),ptrans(mxpart,4),
     & pext(mxpart,4),p5(4),p6(4),p7(4),pt(4),ptt(4),
     & p5out(4),p6out(4),p7out(4)
      integer:: j,nu
      
      do nu=1,4
        pt(nu)=p(5,nu)
        ptt(nu)=ptrans(5,nu)
        p5(nu)=pold(5,nu)
        p6(nu)=pold(6,nu)
        p7(nu)=pold(7,nu)
      enddo
      
      call boostx(p5,pt,ptt,p5out)
      call boostx(p6,pt,ptt,p6out)
      call boostx(p7,pt,ptt,p7out)

      do nu=1,4
        pext(1,nu)=ptrans(1,nu)
        pext(2,nu)=ptrans(2,nu)
        pext(3,nu)=ptrans(3,nu)
        pext(4,nu)=ptrans(4,nu)
        pext(5,nu)=p5out(nu)
        pext(6,nu)=p6out(nu)
        pext(7,nu)=p7out(nu)
        do j=8,mxpart
        pext(j,nu)=zero
        enddo
      enddo
            
      return
      end
      
