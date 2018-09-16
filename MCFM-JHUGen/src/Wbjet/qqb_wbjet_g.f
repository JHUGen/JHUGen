      subroutine qqb_wbjet_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J. M. Campbell                                           *
*     January, 2004.                                                   *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+b(-p2) -->  W^+ + b(p5)+f(p6)+g(p7)
c                        |
c                         --> nu(p3)+e^+(p4)
c -- or
c     q(-p1)+b(-p2) -->  W^- + b(p5)+f(p6)+g(p7)
c                        |
c                         --> e^-(p3)+nu~(p4)
c
c--- Extended to include charm quark production via the variable "flav"


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'heavyflav.f'
      include 'zprods_com.f'
      include 'nwz.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf),
     & Vsm(-nf:nf),nlightf,
     & QQb_bb,QbQ_bb,
     & QQ_us_sd,QQ_su_sd,
     & QbQb_ds_su,QbQb_sd_su,
     & QQb_us_sd,QQb_su_sd,
     & QbQ_us_sd,QbQ_su_sd,
     & QG_d_dsc,QG_u_dcc,
     & GQ_d_dsc,GQ_u_dcc,
     & QbG_u_ucs,QbG_d_ucc,
     & GQb_u_ucs,GQb_d_ucc

      integer:: isub
      common/isub/isub

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(7,p,za,zb)

************************************************************************
*     Calculate contributions from the QQBQQB matrix elements          *
************************************************************************

      if (isub == 1) then
c--- basic Q-Q amplitudes
        call addhel_wbj(1,6,2,5,7,3,4,QQ_us_sd)
        call addhel_wbj(2,6,1,5,7,3,4,QQ_su_sd)
c--- basic Qb-Qb amplitudes
        call addhel_wbj(6,1,5,2,7,3,4,QbQb_ds_su)
        call addhel_wbj(6,2,5,1,7,3,4,QbQb_sd_su)
c--- basic Q-Qb amplitudes
c--- scattering
        call addhel_wbj(1,6,5,2,7,3,4,QQb_us_sd)
        call addhel_wbj(6,2,1,5,7,3,4,QQb_su_sd)
c--- basic Qb-Q amplitudes
c--- scattering
        call addhel_wbj(2,6,5,1,7,3,4,QbQ_su_sd)
        call addhel_wbj(6,1,2,5,7,3,4,QbQ_us_sd)
c--- basic Q-G amplitudes
        call addhel_wbj(7,6,1,5,2,3,4,QG_d_dsc)
c--- basic QB-G amplitudes
        call addhel_wbj(6,7,5,1,2,3,4,QbG_u_ucs)
c--- basic G-Q amplitudes
        call addhel_wbj(7,6,2,5,1,3,4,GQ_d_dsc)
c--- basic G-QB amplitudes
        call addhel_wbj(6,7,5,2,1,3,4,GQb_u_ucs)
      endif

      if (isub == 2) then
c--- gq amplitudes
        call addhel_wbj(1,7,5,6,2,3,4,QG_u_dcc)
        call addhel_wbj(7,1,6,5,2,3,4,QbG_d_ucc)
        call addhel_wbj(2,7,6,5,1,3,4,GQ_u_dcc)
        call addhel_wbj(7,2,5,6,1,3,4,GQb_d_ucc)
c--- qqbar amplitudes
        call addhel_wbj(1,2,6,5,7,3,4,QQb_bb)
        call addhel_wbj(2,1,6,5,7,3,4,QbQ_bb)
      endif

c--- This is a consistency check
c      call addhel_wbj(1,7,6,5,2,3,4,QG_u_dcc)
c      write(6,*) 'addhel_wbj results'
c      write(6,*) 'QG_u_dcc',QG_u_dcc
c      write(6,*) 'wbbgamp result'
c      call wbbgamp(1,7,2,5,6,3,4,qgWbbq)
c      write(6,*) 'qgWbbq',gsq**3*gw**4/4._dp*aveqq*qgWbbq
c      write(6,*) 'ratio',gsq**3*gw**4/4._dp*aveqq*qgWbbq/QG_u_dcc
c      write(6,*) 'redsqm result',
c     & gsq**3*gw**4/4._dp*aveqq*redmsqm(1,7,2,5,6,3,4,zip)
c      pause

c--- set up auxiliary array

      do j=-nf,nf
        Vsm(j)=Vsum(j)
        if (abs(j) >= flav) Vsm(j)=0._dp
c--- make sure that elements are either one or zero
        if (Vsm(j) > 0._dp) Vsm(j)=1._dp
      enddo
c      pause

c--- This sum represents the number of light flavours that may be
c--- attached to the W
c---  for the bottom quark, this is 2
c---  for the charm quark, it is 1 (no extra charm in the final state)
      if (nwz == +1) then
        nlightf=Vsm(+1)+Vsm(+2)+Vsm(+3)+Vsm(+4)
      else
        nlightf=Vsm(-1)+Vsm(-2)+Vsm(-3)+Vsm(-4)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp

      if ((abs(j) .ne. flav) .and. (abs(k) .ne. flav)) goto 97
      if ((abs(j) == flav) .and. (abs(k) == flav)) goto 97
c--- so that either abs(j) or abs(k) = flav (but not both).

c--- contributions with 1 b-quark in the final state
      if (isub == 1) then
        if ((j > 0) .and. (k < 0)) then
c--- e.g.  b d~ -> b u~ or u b~ -> b~ d
             msq(j,k)=+(Vsm(j)*QQb_us_sd+Vsm(k)*QQb_su_sd)
        elseif ((j < 0) .and. (k > 0)) then
c--- e.g.  d~ b -> b u~ or b~ u -> b~ d
             msq(j,k)=+(Vsm(k)*QbQ_su_sd+Vsm(j)*QbQ_us_sd)
        elseif ((j > 0) .and. (k > 0)) then
c--- e.g.  b u -> b d  or u b b d
             msq(j,k)=+(Vsm(k)*QQ_su_sd+Vsm(j)*QQ_us_sd)
        elseif ((j < 0) .and. (k < 0)) then
c--- e.g.  b~ d~ -> b~ u~ or d~ b~ -> b~ u~
             msq(j,k)=+(Vsm(j)*QbQb_ds_su+Vsm(k)*QbQb_sd_su)
        elseif ((j == 0) .and. (k > 0)) then
             msq(j,k)=nlightf*(aveqg/aveqq)*GQ_d_dsc
        elseif ((j == 0) .and. (k < 0)) then
             msq(j,k)=nlightf*(aveqg/aveqq)*GQb_u_ucs
        elseif ((j > 0) .and. (k == 0)) then
             msq(j,k)=nlightf*(aveqg/aveqq)*QG_d_dsc
        elseif ((j < 0) .and. (k == 0)) then
             msq(j,k)=nlightf*(aveqg/aveqq)*QbG_u_ucs
        endif
      endif

 97   continue

      if (isub == 2) then
c--- contributions with 2 b-quarks in the final state
c--- Note interchange of b and bbar below.
         if ((j > 0) .and. (k == 0)) then
c---     q G --> l lb Qb Q q'
              msq(j,k)=msq(j,k)+Vsm(j)*(aveqg/aveqq)*QG_u_dcc
         elseif ((j < 0) .and. (k == 0)) then
c---     qb G --> l lb Q Qb qb'
              msq(j,k)=msq(j,k)+Vsm(j)*(aveqg/aveqq)*QbG_d_ucc
         elseif ((j == 0) .and. (k > 0)) then
c---     G q --> l lb Q Qb q'
              msq(j,k)=msq(j,k)+Vsm(k)*(aveqg/aveqq)*GQ_u_dcc
         elseif ((j == 0) .and. (k < 0)) then
c---     G qb --> l lb Qb Q qb'
              msq(j,k)=msq(j,k)+Vsm(k)*(aveqg/aveqq)*GQb_d_ucc
c---     q qb' --> l lb Q Qb
         elseif ((j > 0) .and. (k < 0)) then
              msq(j,k)=msq(j,k)+Vsq(j,k)*QQb_bb
c---     qb' q --> l lb Q Qb
         elseif ((j < 0) .and. (k > 0)) then
              msq(j,k)=msq(j,k)+Vsq(j,k)*QbQ_bb
         endif
      endif

      enddo
      enddo


      return
      end

