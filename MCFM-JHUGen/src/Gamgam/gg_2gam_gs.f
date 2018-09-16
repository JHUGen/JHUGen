      subroutine gg_2gam_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     October, 2002.                                                   *
*     Modified by CW to Include photon fragmentation dipoles Feb 11    *
*    Matrix element SUBTRACTION squared averag'd over init'l colors    *
*    and spins                                                         *
*     g(-p1) + g(-p2) -->  gamma(p3) + gamma(p4) + g(p5)               *
************************************************************************
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'frag.f'
      include 'ewcharge.f'
      include 'phot_dip.f'
      integer:: j,k,nd
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),sub15_2v,sub25_1v,
     & msq15_2v(-nf:nf,-nf:nf),msq25_1v(-nf:nf,-nf:nf)
      external gg_2gam,gg_2gam_gvec

      ndmax=2

c---- calculate both initial-initial dipoles
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & gg_2gam,gg_2gam_gvec)
      call dips(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & gg_2gam,gg_2gam_gvec)

      msq(:,:,:)=0._dp

      msq(1,0,0)=2._dp*xn*(sub15_2(gg)*msq15_2(0,0)
     &                    +sub15_2v*msq15_2v(0,0))
      msq(2,0,0)=2._dp*xn*(sub25_1(gg)*msq25_1(0,0)
     &                    +sub25_1v*msq25_1v(0,0))

      do j=1,nf
      msq(1,+j,0)=(aveqg/avegg)*(sub15_2(gq)*msq15_2(0,0)
     &                          +sub15_2v*msq15_2v(0,0))
      msq(1,-j,0)=(aveqg/avegg)*(sub15_2(gq)*msq15_2(0,0)
     &                          +sub15_2v*msq15_2v(0,0))
      msq(2,0,+j)=(aveqg/avegg)*(sub25_1(gq)*msq25_1(0,0)
     &                          +sub25_1v*msq25_1v(0,0))
      msq(2,0,-j)=(aveqg/avegg)*(sub25_1(gq)*msq25_1(0,0)
     &                          +sub25_1v*msq25_1v(0,0))
      enddo
      
      return
      end
      
