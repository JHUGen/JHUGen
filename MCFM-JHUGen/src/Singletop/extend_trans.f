      subroutine extend_trans(pold,p,ptrans,pext)
      implicit none
      include 'types.f'
      
c--- take vector 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: pold(mxpart,4),p(mxpart,4),ptrans(mxpart,4),
     & pext(mxpart,4),p3(4),p4(4),p5(4),pt(4),ptt(4),
     & p3out(4),p4out(4),p5out(4)
      integer:: nu
      
      do nu=1,4
        pt(nu)=p(3,nu)
        ptt(nu)=ptrans(3,nu)
        p3(nu)=pold(3,nu)
        p4(nu)=pold(4,nu)
        p5(nu)=pold(5,nu)
      enddo
      
      call boostx(p3,pt,ptt,p3out)
      call boostx(p4,pt,ptt,p4out)
      call boostx(p5,pt,ptt,p5out)

      do nu=1,4
        pext(1,nu)=ptrans(1,nu)
        pext(2,nu)=ptrans(2,nu)
        pext(3,nu)=p3out(nu)
        pext(4,nu)=p4out(nu)
        pext(5,nu)=p5out(nu)
        pext(6,nu)=ptrans(4,nu)
      enddo
            
      return
      end
      
