      subroutine qg_tbqdk_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Real subtraction for t-channel single top, with explicit b-quark *
*                                                                      *
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5) + g(p6)                  *
*                                                                      *
*      (and related crossings and MEs)                                 *
*                                                                      *
*     Author: J. Campbell, March 18, 2008                              *
*                         (added decay May 2011)                       *
*                                                                      *
************************************************************************
*                                                                      *
*    Modified so that off-diagonal subtractions are initial-final      *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ckm.f'
      include 'ptilde.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'noglue.f'
      include 'stopscales.f'
      include 'breit.f'
      real(dp):: corrL,corrH
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp)::
     & msq16_5(-nf:nf,-nf:nf),sub16_5(4),
     & msq56_1(-nf:nf,-nf:nf),sub56_1(4),
     & msq26_4(-nf:nf,-nf:nf),sub26_4(4),
     & msq46_2(-nf:nf,-nf:nf),sub46_2(4),
     & msq26_3(-nf:nf,-nf:nf),sub26_3(4),
     & msq36_2(-nf:nf,-nf:nf),sub36_2(4),
     & msq36_4(-nf:nf,-nf:nf),sub36_4(4),
     & msq46_3(-nf:nf,-nf:nf),sub46_3(4),
     & msq26_4v(-nf:nf,-nf:nf),sub26_4v,
     & msq26_3v(-nf:nf,-nf:nf),sub26_3v,
     & msq26_5(-nf:nf,-nf:nf),sub26_5(4),
     & msq56_2(-nf:nf,-nf:nf),sub56_2(4),
     & msq16_4(-nf:nf,-nf:nf),sub16_4(4),
     & msq46_1(-nf:nf,-nf:nf),sub46_1(4),
     & msq16_3(-nf:nf,-nf:nf),sub16_3(4),
     & msq36_1(-nf:nf,-nf:nf),sub36_1(4),
     & msq16_4v(-nf:nf,-nf:nf),sub16_4v,
     & msq16_3v(-nf:nf,-nf:nf),sub16_3v,
c     & msq16_2(-nf:nf,-nf:nf),sub16_2(4),
c     & msq26_1(-nf:nf,-nf:nf),sub26_1(4),
c     & msq16_2v(-nf:nf,-nf:nf),sub16_2v,
c     & msq26_1v(-nf:nf,-nf:nf),sub26_1v,
c     & msq15_2v(-nf:nf,-nf:nf),sub15_2v,
c     & msq25_1v(-nf:nf,-nf:nf),sub25_1v,
     & msq16_5v(-nf:nf,-nf:nf),sub16_5v,
     & msq26_5v(-nf:nf,-nf:nf),sub26_5v,
     & msq15_6(-nf:nf,-nf:nf),sub15_6(4),
     & msq25_6(-nf:nf,-nf:nf),sub25_6(4),
     & msq25_6v(-nf:nf,-nf:nf),sub25_6v,
     & msq15_6v(-nf:nf,-nf:nf),sub15_6v,
     & dummyv(-nf:nf,-nf:nf),dsubv
      integer:: j,k,nd
      real(dp):: oldmass2
      external qg_tbqdk,qg_tbqdk_gvec,donothing_gvec

c--- Note that the subtractions here must separate the initial-final
c--- and final-initial dipoles because of the top/bottom masses
      ndmax=14

      qqproc=.true.
      qgproc=.false.
      gqproc=.false.
      ggproc=.false.

c--- save original value of mass2
      oldmass2=mass2

c--- NOTE: to pass to the "decay" (dk) version of this file, the labels
c---       of the dipoles have not been changed.
c---       For the calls to dips_mass, the momenta
c---       to which they refer will be transformed according to:
c---       3 -> (3+4+5), 4 -> 6, 5 -> 7, 6 -> 8

c--- DIPOLES REQUIRED FOR OFF-DIAGONAL SUBTRACTIONS
c--- subtractions for massless line
      mass2=zip
      call dips_mass(1,p,1,6,5,sub16_5,sub16_5v,msq16_5,msq16_5v,
     & qg_tbqdk,qg_tbqdk_gvec)
      call dips_mass(1,p,5,6,1,sub56_1,dsubv,msq56_1,dummyv,
     & qg_tbqdk,donothing_gvec)

c--- subtractions for massless line
      call dips_mass(2,p,2,6,5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     & qg_tbqdk,qg_tbqdk_gvec)
      call dips_mass(2,p,5,6,2,sub56_2,dsubv,msq56_2,dummyv,
     & qg_tbqdk,donothing_gvec)

      call dips_mass(3,p,1,5,6,sub15_6,sub15_6v,msq15_6,msq15_6v,
     & qg_tbqdk,qg_tbqdk_gvec)
      call dips_mass(4,p,2,5,6,sub25_6,sub25_6v,msq25_6,msq25_6v,
     & qg_tbqdk,qg_tbqdk_gvec)


c--- implement shortcuts for checking
      if ((noglue) .or. (ggonly)) then
      ndmax=4 ! only need the first four dipoles
      else

c--- DIPOLES FOR QG AND QBARG
c--- subtractions for massive line
c--- note: mass of final state particle for final-initial and
c---  final-final dipoles is passed in (clumsily) via mass2
      mass2=mb
      call dips_mass(5,p,2,6,4,sub26_4,sub26_4v,msq26_4,msq26_4v,
     & qg_tbqdk,qg_tbqdk_gvec)
      call dips_mass(6,p,4,6,2,sub46_2,dsubv,msq46_2,dummyv,
     & qg_tbqdk,donothing_gvec)
      mass2=mt
      call dips_mass(7,p,2,6,3,sub26_3,sub26_3v,msq26_3,msq26_3v,
     & qg_tbqdk,qg_tbqdk_gvec)
      call dips_mass(8,p,3,6,2,sub36_2,dsubv,msq36_2,dummyv,
     & qg_tbqdk,donothing_gvec)
      mass2=mt
      call dips_mass(9,p,3,6,4,sub36_4,dsubv,msq36_4,dummyv,
     & qg_tbqdk,donothing_gvec)
      mass2=mb
      call dips_mass(10,p,4,6,3,sub46_3,dsubv,msq46_3,dummyv,
     & qg_tbqdk,donothing_gvec)

c--- EXTRA DIPOLES FOR GQ AND GQBAR
c--- subtractions for massive line
c--- note: mass of final state particle for final-initial and
c---  final-final dipoles is passed in (clumsily) via mass2
      mass2=mb
      call dips_mass(11,p,1,6,4,sub16_4,sub16_4v,msq16_4,msq16_4v,
     & qg_tbqdk,qg_tbqdk_gvec)
      call dips_mass(12,p,4,6,1,sub46_1,dsubv,msq46_1,dummyv,
     & qg_tbqdk,donothing_gvec)
      mass2=mt
      call dips_mass(13,p,1,6,3,sub16_3,sub16_3v,msq16_3,msq16_3v,
     & qg_tbqdk,qg_tbqdk_gvec)
      call dips_mass(14,p,3,6,1,sub36_1,dsubv,msq36_1,dummyv,
     & qg_tbqdk,donothing_gvec)

      endif

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
      corrL=as_L/as
      corrH=as_H/as

c      write(6,*) 'msq(5,2,0)',
c     & msq26_4(2,0),sub26_4(gg),msq26_4v(2,0),sub26_4v

      do j=-4,4
      do k=-4,4

      if ((noglue) .and. (j*k == 0)) goto 99

      if ((ggonly) .and. ((j .ne. 0) .or. (k .ne. 0))) goto 99

      if     ((j .ne. 0) .and. (k == 0)) then
c--- subtractions for qg and qbarg
        msq(1,j,k)=2._dp*cf*msq16_5(j,k)*(sub16_5(qq)+sub56_1(qq))*corrL
        msq(5,j,k)=xn*(msq26_4(j,k)*sub26_4(gg)+msq26_4v(j,k)*sub26_4v)
     &               *corrH
        msq(6,j,k)=xn*(msq46_2(j,k)*sub46_2(qq))*corrH
        msq(7,j,k)=xn*(msq26_3(j,k)*sub26_3(gg)+msq26_3v(j,k)*sub26_3v)
     &               *corrH
        msq(8,j,k)=xn*(msq36_2(j,k)*sub36_2(qq))*corrH
        msq(9,j,k)=-(msq36_4(j,k)*sub36_4(qq))/xn*corrH
        msq(10,j,k)=-(msq46_3(j,k)*sub46_3(qq))/xn*corrH
      elseif ((j == 0) .and. (k .ne. 0)) then
c--- subtractions for gq and gqbar
        msq( 2,j,k)=2._dp*cf*msq26_5(j,k)*(sub26_5(qq)+sub56_2(qq))*corrL
        msq(11,j,k)=xn*(msq16_4(j,k)*sub16_4(gg)+msq16_4v(j,k)*sub16_4v)
     &                *corrH
        msq(12,j,k)=xn*(msq46_1(j,k)*sub46_1(qq))*corrH
        msq(13,j,k)=xn*(msq16_3(j,k)*sub16_3(gg)+msq16_3v(j,k)*sub16_3v)
     &                *corrH
        msq(14,j,k)=xn*(msq36_1(j,k)*sub36_1(qq))*corrH
        msq( 9,j,k)=-(msq36_4(j,k)*sub36_4(qq))/xn*corrH
        msq(10,j,k)=-(msq46_3(j,k)*sub46_3(qq))/xn*corrH
      elseif ((j == 0) .and. (k == 0)) then
c--- subtractions for gg
       msq(1,j,k)=2._dp*tr*(msq16_5(4,k)+msq16_5(3,k)
     &                   +msq16_5(2,k)+msq16_5(1,k))*sub16_5(qg)*corrL
       msq(2,j,k)=2._dp*tr*(msq26_5(j,4)+msq26_5(j,3)
     &                   +msq26_5(j,2)+msq26_5(j,1))*sub26_5(qg)*corrL
       msq(3,j,k)=2._dp*tr*(msq15_6(-4,k)+msq15_6(-3,k)
     &                   +msq15_6(-2,k)+msq15_6(-1,k))*sub15_6(qg)*corrL
       msq(4,j,k)=2._dp*tr*(msq25_6(j,-4)+msq25_6(j,-3)
     &                   +msq25_6(j,-2)+msq25_6(j,-1))*sub25_6(qg)*corrL
      elseif (j*k < 0) then
c--- subtractions for qqbar/qbarq
       msq(1,j,k)=2._dp*cf*(msq16_5(0,k)*sub16_5(gq)
     &                   +msq16_5v(0,k)*sub16_5v)*corrH
       msq(2,j,k)=2._dp*cf*(msq26_5(j,0)*sub26_5(gq)
     &                   +msq26_5v(j,0)*sub26_5v)*corrH
      elseif (j*k > 0) then
c--- subtractions for qq/qbarqbar
       msq(1,j,k)=(2._dp*Vsum(k)-Vsq(-j,k))*cf*(
     &       msq16_5(0,k)*sub16_5(gq)+msq16_5v(0,k)*sub16_5v)*corrH
       msq(2,j,k)=(2._dp*Vsum(j)-Vsq(j,-k))*cf*(
     &       msq26_5(j,0)*sub26_5(gq)+msq26_5v(j,0)*sub26_5v)*corrH
       msq(3,j,k)=Vsq(-j,k)*cf*(
     &       msq15_6(0,k)*sub15_6(gq)+msq15_6v(0,k)*sub15_6v)*corrH
       msq(4,j,k)=Vsq(j,-k)*cf*(
     &       msq25_6(j,0)*sub25_6(gq)+msq25_6v(j,0)*sub25_6v)*corrH
      endif

   99 continue

      enddo
      enddo

c      call qg_tbqdk(p,dummyv)
c
c      e25=cf*gsq*(2._dp*dot(p,2,5))
c     & /dot(p,2,6)/dot(p,5,6)*dummyv(0,4)
c      e13=xn/2._dp*gsq*dummyv(0,4)*(
c     & 2._dp*dot(p,1,3)/dot(p,1,6)-mt**2/dot(p,3,6))/dot(p,3,6)
c      e14=xn/2._dp*gsq*dummyv(0,4)*(
c     & 2._dp*dot(p,1,4)/dot(p,1,6)-mb**2/dot(p,4,6))/dot(p,4,6)
c      e34=-1._dp/2._dp/xn*gsq*dummyv(0,4)*(
c     & 2._dp*dot(p,3,4)/dot(p,3,6)/dot(p,4,6)
c     & -mb**2/dot(p,4,6)**2-mt**2/dot(p,3,6)**2)
c      write(6,*) '(0,4)'
c      write(6,*) 'eikonal 25',e25
c      write(6,*) 'eikonal 13',e13
c      write(6,*) 'eikonal 14',e14
c      write(6,*) 'eikonal 34',e34
c      write(6,*) '   sum    ',e25+e13+e14+e34

c      e15=cf*gsq*(2._dp*dot(p,1,5))
c     & /dot(p,1,6)/dot(p,5,6)*dummyv(4,0)
c      e23=xn/2._dp*gsq*dummyv(4,0)*(
c     & 2._dp*dot(p,2,3)/dot(p,2,6)-mt**2/dot(p,3,6))/dot(p,3,6)
c      e24=xn/2._dp*gsq*dummyv(4,0)*(
c     & 2._dp*dot(p,2,4)/dot(p,2,6)-mb**2/dot(p,4,6))/dot(p,4,6)
c      e34=-1._dp/2._dp/xn*gsq*dummyv(4,0)*(
c     & 2._dp*dot(p,3,4)/dot(p,3,6)/dot(p,4,6)
c     & -mb**2/dot(p,4,6)**2-mt**2/dot(p,3,6)**2)
c      write(6,*) '(4,0)'
c      write(6,*) 'eikonal 15',e15
c      write(6,*) 'eikonal 23',e23
c      write(6,*) 'eikonal 24',e24
c      write(6,*) 'eikonal 34',e34
c      write(6,*) '   sum    ',e15+e23+e24+e34

      return
      end

