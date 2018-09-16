      subroutine qq_tbg_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Real subtraction for s-channel single top + jet                  *
*                                                                      *
*     q(p1) + q(p2) -> t(p3) + b(p4) + g(p5) + g(p6)                   *
*                                                                      *
*      (and related crossings and MEs)                                 *
*                                                                      *
*     Author: J. Campbell, June 23, 2008                               *
*                                                                      *
************************************************************************
*                                                                      *
*     IMPORTANT NOTE!                                                  *
*                                                                      *
*     For now, we only include radiation from heavy quark line         *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ptilde.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'nflav.f'
      include 'stopscales.f'
      include 'breit.f'
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),corrH
      real(dp)::
     & msq35_4(-nf:nf,-nf:nf),sub35_4(4),
     & msq36_4(-nf:nf,-nf:nf),sub36_4(4),
     & msq35_6(-nf:nf,-nf:nf),sub35_6(4),
     & msq36_5(-nf:nf,-nf:nf),sub36_5(4),
     & msq45_3(-nf:nf,-nf:nf),sub45_3(4),
     & msq46_3(-nf:nf,-nf:nf),sub46_3(4),
     & msq45_6(-nf:nf,-nf:nf),sub45_6(4),
     & msq46_5(-nf:nf,-nf:nf),sub46_5(4),
     & msq56_3(-nf:nf,-nf:nf),sub56_3(4),
     & msq56_4(-nf:nf,-nf:nf),sub56_4(4),
     & msq56_3v(-nf:nf,-nf:nf),sub56_3v_gg,sub56_3v_gq,
     & msq56_4v(-nf:nf,-nf:nf),sub56_4v_gg,sub56_4v_gq,
     & dummyv(-nf:nf,-nf:nf),dsubv,
     & subv_gg,subv_gq
      integer:: j,k,nd
      real(dp):: oldmass2
      common/subv_ff/subv_gg,subv_gq
!$omp threadprivate(/subv_ff/)
      external qq_tbg,qq_tbg_gvec,donothing_gvec

      ndmax=10

c--- initialize
      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

      qqproc=.true.
      qgproc=.false.
      gqproc=.false.
      ggproc=.false.

c--- save original value of mass2
      oldmass2=mass2

c--- subtractions for massive line
c--- note: mass of final state particle for final-final dipoles
c---       is passed in (clumsily) via mass2
      mass2=mt
      call dips_mass( 1,p,3,5,4,sub35_4,dsubv,msq35_4,dummyv,
     & qq_tbg,donothing_gvec)
      call dips_mass( 2,p,3,6,4,sub36_4,dsubv,msq36_4,dummyv,
     & qq_tbg,donothing_gvec)
      call dips_mass( 3,p,3,5,6,sub35_6,dsubv,msq35_6,dummyv,
     & qq_tbg,donothing_gvec)
      call dips_mass( 4,p,3,6,5,sub36_5,dsubv,msq36_5,dummyv,
     & qq_tbg,donothing_gvec)
      mass2=mb
      call dips_mass( 5,p,4,5,3,sub45_3,dsubv,msq45_3,dummyv,
     & qq_tbg,donothing_gvec)
      call dips_mass( 6,p,4,6,3,sub46_3,dsubv,msq46_3,dummyv,
     & qq_tbg,donothing_gvec)
      call dips_mass( 7,p,4,5,6,sub45_6,dsubv,msq45_6,dummyv,
     & qq_tbg,donothing_gvec)
      call dips_mass( 8,p,4,6,5,sub46_5,dsubv,msq46_5,dummyv,
     & qq_tbg,donothing_gvec)

c--- these dipoles are for g->gg or g->qq~ (mq=0) splittings, so that
c--- the dipole transformations are identical and ME's are the same
      qqproc=.false.
      qgproc=.false.
      gqproc=.true.
      ggproc=.true.
      mass2=0._dp
      call dips_mass( 9,p,5,6,3,sub56_3,dsubv,msq56_3,msq56_3v,
     & qq_tbg,qq_tbg_gvec)
      sub56_3v_gg=subv_gg
      sub56_3v_gq=subv_gq
      call dips_mass(10,p,5,6,4,sub56_4,dsubv,msq56_4,msq56_4v,
     & qq_tbg,qq_tbg_gvec)
      sub56_4v_gg=subv_gg
      sub56_4v_gq=subv_gq

c--- reset mass2 to original value
      mass2=oldmass2

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

c--- correction factors to replace the "gsq" that appears in the sub...
c--- expressions with the correct gsq for that line
c      corrL=as_L/as
      corrH=as_H/as

c--- leading colour (Ca) goes in (0,0)
c--- subleading colour (2*Cf-Ca) goes in (0,1)
c--- Tr terms go in (1,0)

c--- leading contributions
      msq( 9,0,0)=(half*ca*(
     &               msq56_3(0,0)*sub56_3(gg)+msq56_3v(0,0)*sub56_3v_gg)
     &              )*corrH
      msq( 9,1,0)=(half*2._dp*Tr*real(nflav,dp)*(
     &               msq56_3(0,0)*sub56_3(gq)-msq56_3v(0,0)*sub56_3v_gq)
     &              )*corrH
      msq(10,0,0)=(half*ca*(
     &               msq56_4(0,0)*sub56_4(gg)+msq56_4v(0,0)*sub56_4v_gg)
     &              )*corrH
      msq(10,1,0)=(half*2._dp*Tr*real(nflav,dp)*(
     &               msq56_4(0,0)*sub56_4(gq)-msq56_4v(0,0)*sub56_4v_gq)
     &              )*corrH
      msq( 3,0,0)=half*ca*(msq35_6(0,0)*sub35_6(qq))*corrH
      msq( 4,0,0)=half*ca*(msq36_5(0,0)*sub36_5(qq))*corrH
      msq( 7,0,0)=half*ca*(msq45_6(0,0)*sub45_6(qq))*corrH
      msq( 8,0,0)=half*ca*(msq46_5(0,0)*sub46_5(qq))*corrH
c--- subleading contributions
      msq( 1,0,1)=half*(2._dp*cf-ca)*(msq35_4(0,0)*sub35_4(qq))*corrH
      msq( 2,0,1)=half*(2._dp*cf-ca)*(msq36_4(0,0)*sub36_4(qq))*corrH
      msq( 5,0,1)=half*(2._dp*cf-ca)*(msq45_3(0,0)*sub45_3(qq))*corrH
      msq( 6,0,1)=half*(2._dp*cf-ca)*(msq46_3(0,0)*sub46_3(qq))*corrH

c      do j=-4,4
c      do k=-4,4
c
c      if     ((j > 0) .and. (k < 0)) then
cc--- leading contributions
c        msq( 9,j,k)=(half*ca*(
c     &               msq56_3(j,k)*sub56_3(gg)+msq56_3v(j,k)*sub56_3v_gg)
c     &              +half*real(nflav,dp)*(
c     &               msq56_3(j,k)*sub56_3(gq)-msq56_3v(j,k)*sub56_3v_gq)
c     &              )*corrH
c        msq(10,j,k)=(half*ca*(
c     &               msq56_4(j,k)*sub56_4(gg)+msq56_4v(j,k)*sub56_4v_gg)
c     &              +half*real(nflav,dp)*(
c     &               msq56_4(j,k)*sub56_4(gq)-msq56_4v(j,k)*sub56_4v_gq)
c     &              )*corrH
c        msq( 3,j,k)=half*ca*(msq35_6(j,k)*sub35_6(qq))*corrH
c        msq( 4,j,k)=half*ca*(msq36_5(j,k)*sub36_5(qq))*corrH
c        msq( 7,j,k)=half*ca*(msq45_6(j,k)*sub45_6(qq))*corrH
c        msq( 8,j,k)=half*ca*(msq46_5(j,k)*sub46_5(qq))*corrH
cc--- subleading contributions
c      msq( 1,j,k)=half*(2._dp*cf-ca)*(msq35_4(j,k)*sub35_4(qq))*corrH
c      msq( 2,j,k)=half*(2._dp*cf-ca)*(msq36_4(j,k)*sub36_4(qq))*corrH
c      msq( 5,j,k)=half*(2._dp*cf-ca)*(msq45_3(j,k)*sub45_3(qq))*corrH
c      msq( 6,j,k)=half*(2._dp*cf-ca)*(msq46_3(j,k)*sub46_3(qq))*corrH
c      elseif ((j < 0) .and. (k > 0)) then
cc--- leading contributions
c        msq( 9,j,k)=(half*ca*(
c     &               msq56_3(j,k)*sub56_3(gg)+msq56_3v(j,k)*sub56_3v_gg)
c     &              +half*real(nflav,dp)*(
c     &               msq56_3(j,k)*sub56_3(gq)-msq56_3v(j,k)*sub56_3v_gq)
c     &              )*corrH
c        msq(10,j,k)=(half*ca*(
c     &               msq56_4(j,k)*sub56_4(gg)+msq56_4v(j,k)*sub56_4v_gg)
c     &              +half*real(nflav,dp)*(
c     &               msq56_4(j,k)*sub56_4(gq)-msq56_4v(j,k)*sub56_4v_gq)
c     &              )*corrH
c        msq( 3,j,k)=half*ca*(msq35_6(j,k)*sub35_6(qq))*corrH
c        msq( 4,j,k)=half*ca*(msq36_5(j,k)*sub36_5(qq))*corrH
c        msq( 7,j,k)=half*ca*(msq45_6(j,k)*sub45_6(qq))*corrH
c        msq( 8,j,k)=half*ca*(msq46_5(j,k)*sub46_5(qq))*corrH
cc--- subleading contributions
c      msq( 1,j,k)=half*(2._dp*cf-ca)*(msq35_4(j,k)*sub35_4(qq))*corrH
c      msq( 2,j,k)=half*(2._dp*cf-ca)*(msq36_4(j,k)*sub36_4(qq))*corrH
c      msq( 5,j,k)=half*(2._dp*cf-ca)*(msq45_3(j,k)*sub45_3(qq))*corrH
c      msq( 6,j,k)=half*(2._dp*cf-ca)*(msq46_3(j,k)*sub46_3(qq))*corrH
c      endif
c
c   99 continue
c
c      enddo
c      enddo

      return
      end

