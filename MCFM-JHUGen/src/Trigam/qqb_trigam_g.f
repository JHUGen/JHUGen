      subroutine qqb_trigam_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*    Author: J. M. Campbell                                            *
*    March, 2013.                                                      *
*    Real matrix element squared, averaged over initial colors         *
*    and spins (plus all crossings)                                    *
c     q(-p1)+qbar(-p2) --> gam(p3) + gam(p4) + gam(p5) + g(p6)         *
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
     & ampsq_3gam1g,qqb,qbq,qg,gq,qbg,gqb,fac,cfac
      
      call spinoru(6,p,za,zb)

c--- overall coupling, color and identical particle factors
      fac=16._dp*esq**3*gsq*xn*Cf/6._dp
      
      qqb=aveqq*fac*ampsq_3gam1g(6,3,4,5,1,2,za,zb)
      qg=aveqg*fac*ampsq_3gam1g(2,3,4,5,1,6,za,zb)
      gq=aveqg*fac*ampsq_3gam1g(1,3,4,5,2,6,za,zb)
      qbq=qqb
      qbg=qg
      gqb=gq
      
c      qbq=aveqq*fac*ampsq_3gam1g(6,3,4,5,2,1,za,zb)
c      qbg=aveqg*fac*ampsq_3gam1g(2,3,4,5,6,1,za,zb)
c      gqb=aveqg*fac*ampsq_3gam1g(1,3,4,5,6,2,za,zb)

      msq(:,:)=0._dp
      do j=1,nf
        cfac=Q(j)**6
        msq(j,-j)=cfac*qqb
        msq(-j,j)=cfac*qbq
        msq(j,0)=cfac*qg
        msq(0,j)=cfac*gq
        msq(-j,0)=cfac*qbg
        msq(0,-j)=cfac*gqb
      enddo
      
      return
      end
      
      
