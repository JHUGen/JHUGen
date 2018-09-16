      subroutine qqb_fourgam(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*    Author: C Williams                                                *
*    Sept, 2014.                                                       *
*    LO matrix element squared, averaged over initial colors           *
*    and spins (plus all crossings)                                    *
c     q(-p1)+qbar(-p2) --> gam(p3) + gam(p4) + gam(p5) + gam(p6)       *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'zprods_decl.f'
      integer:: j
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & ampsq_3gam1g,qqb,qbq,fac,cfac
      real(dp):: symfac 
      parameter (symfac=one/24d0)

      call spinoru(6,p,za,zb)

c--- overall coupling, color and identical particle factors
      fac=16d0*esq**4*xn*symfac
      qqb=aveqq*fac*ampsq_3gam1g(3,4,5,6,1,2,za,zb)
      qbq=aveqq*fac*ampsq_3gam1g(3,4,5,6,2,1,za,zb)

      msq(:,:)=0d0
      do j=1,5
        cfac=Q(j)**8
        msq(j,-j)=cfac*qqb
        msq(-j,j)=cfac*qbq
      enddo
      
      return
      end
      
      
