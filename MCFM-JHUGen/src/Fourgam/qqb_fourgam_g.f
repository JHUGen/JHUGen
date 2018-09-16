      subroutine qqb_fourgam_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*    Author: C Williams                                                *
*    Sept, 2014.                                                      *
*    LO matrix element squared, averaged over initial colors         *
*    and spins (plus all crossings)                                    *
c     q(-p1)+qbar(-p2) --> gam(p3) + gam(p4) + gam(p5) + gam(p6) +g(p7)        *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      integer:: j
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & qqb,qbq,qg,gq,qbg,gqb,fac,cfac
      real(dp):: symfac 
      parameter (symfac=one/24d0)
      real(dp):: real_aaaaj
      external real_aaaaj

      call spinoru(7,p,za,zb)

c--- overall coupling, color and identical particle factors
      fac=16d0*esq**4*xn*symfac*gsq*Cf*2d0
      
      qqb=aveqq*fac*real_aaaaj(1,2,3,4,5,6,7,za,zb)
      qg=aveqg*fac*real_aaaaj(1,7,3,4,5,6,2,za,zb)
      gq=aveqg*fac*real_aaaaj(2,7,3,4,5,6,1,za,zb)

      qbq=aveqq*fac*real_aaaaj(2,1,3,4,5,6,7,za,zb)
      qbg=aveqg*fac*real_aaaaj(7,1,3,4,5,6,2,za,zb)
      gqb=aveqg*fac*real_aaaaj(7,2,3,4,5,6,1,za,zb)

      msq(:,:)=0d0
      do j=1,nf
        cfac=Q(j)**8
        msq(j,-j)=cfac*qqb
        msq(-j,j)=cfac*qbq
        msq(j,0)=cfac*qg
        msq(0,j)=cfac*gq
        msq(-j,0)=cfac*qbg
        msq(0,-j)=cfac*gqb
      enddo
      

      
      return
      end
      
      
