************************************************************************
*     This is the W- routine                                           *
************************************************************************
      subroutine qqb_wm2jet_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J. M. Campbell                                           *
*     July, 2001.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  g*  + W^+ + g(p7)
c                           |     |
c                           |     --> nu(p3)+e^+(p4)
c                           |
c                           ---> f(p5)+f(p6)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'flags.f'
      include 'lc.f'
      include 'first.f'
      include 'mpicommon.f'
      integer:: j,k,n1,n2
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: mmsq_gg,mmsq_qqb,mmsq_qbq,mmsq_qg
      real(dp):: mmsq_gq,mmsq_gqb,mmsq_qbg
      real(dp):: fac,Vfac
      real(dp):: propsq
      real(dp)::
     & QQ_ud_dd,QQ_us_ds,QQ_uu_du,QQ_du_dd,QQ_su_sd,QQ_uu_ud,
     .QbQb_ds_us,QbQb_du_uu,QbQb_dd_ud,QbQb_sd_su,QbQb_ud_uu,QbQb_dd_du,
     & RRb_dd_du,RRb_ss_du,RRb_uu_du,RRb_ud_dd,RRb_ud_ss,RRb_ud_uu,
     & QQb_ud_dd,QQb_us_ds,QQb_ud_uu,QQb_du_dd,QQb_su_sd,QQb_uu_ud,
     & RbR_dd_du,RbR_ss_du,RbR_uu_du,RbR_ud_dd,RbR_ud_ss,RbR_ud_uu,
     & QbQ_ud_dd,QbQ_us_ds,QbQ_ud_uu,QbQ_du_dd,QbQ_su_sd,QbQ_uu_ud,
     & QG_d_ddu,QG_d_dsc,QG_u_ddd,QG_u_dcc,QG_u_duu,QG_u_udu,
     & GQ_d_ddu,GQ_d_dsc,GQ_u_ddd,GQ_u_dcc,GQ_u_duu,GQ_u_udu,
     & QbG_d_ddu,QbG_u_ucs,QbG_u_uud,QbG_d_uuu,QbG_d_ucc,QbG_d_udd,
     & GQb_d_ddu,GQb_u_ucs,GQb_u_uud,GQb_d_udd,GQb_d_ucc,GQb_d_uuu
      integer,parameter :: jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter :: kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      if (first) then
      first=.false.
        if (rank == 0) then
        if ((Gflag) .or. (QandGflag)) then
          write(*,*) 'Using QQGG+G (REAL) matrix elements'
!          write(*,*) '[LC is     N   ]'
!          write(*,*) '[SLC is   1/N  ]'
!          write(*,*) '[SSLC is 1/N**3]'
        endif
        if ((Qflag) .or. (QandGflag)) then
          write(*,*) 'Using QQBQQB+G (REAL) matrix elements'
!          write(*,*) '[LC is   1 ]'
!          write(*,*) '[SLC is 1/N]'
        endif
        if     (colourchoice == 1) then
          write(*,*) 'Leading colour only in REAL'
        elseif (colourchoice == 2) then
          write(*,*) 'Sub-leading colour only in REAL'
        elseif (colourchoice == 3) then
          write(*,*) 'Sub-sub-leading colour only in REAL'
        elseif (colourchoice == 0) then
          write(*,*) 'Total of all colour structures in REAL'
        else
          write(*,*) 'Bad colourchoice'
          stop
        endif
        endif
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(7,p,za,zb)
      propsq=s(3,4)**2/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)

      if (Gflag) then
************************************************************************
*     Calculate contributions from the QQGG matrix elements            *
************************************************************************

      call spinoru(7,p,za,zb)
      call xwqqggg(5,1,2,7,6,3,4,mmsq_gg)
      call xwqqggg(1,5,6,7,2,3,4,mmsq_qqb)
      call xwqqggg(2,5,6,7,1,3,4,mmsq_qbq)
      call xwqqggg(1,2,6,7,5,3,4,mmsq_qg)
      call xwqqggg(2,1,6,7,5,3,4,mmsq_gq)
      call xwqqggg(5,1,6,7,2,3,4,mmsq_gqb)
      call xwqqggg(5,2,6,7,1,3,4,mmsq_qbg)

      endif

      if (Qflag) then
************************************************************************
*     Calculate contributions from the QQBQQB matrix elements          *
************************************************************************

c--- basic Q-Q amplitudes
      call addhel(1,5,2,6,7,3,4,QQ_ud_dd,QQ_us_ds,QQ_uu_du)
      call addhel(2,6,1,5,7,3,4,QQ_du_dd,QQ_su_sd,QQ_uu_ud)
c--- basic Qb-Qb amplitudes
      call addhel(5,1,6,2,7,3,4,QbQb_dd_ud,QbQb_ds_us,QbQb_du_uu)
      call addhel(6,2,5,1,7,3,4,QbQb_dd_du,QbQb_sd_su,QbQb_ud_uu)
c--- basic Q-Qb amplitudes
c--- annihilation
      call addhel(6,5,1,2,7,3,4,RRb_dd_du,RRb_ss_du,RRb_uu_du)
      call addhel(1,2,6,5,7,3,4,RRb_ud_dd,RRb_ud_ss,RRb_ud_uu)
c--- scattering
      call addhel(1,5,6,2,7,3,4,QQb_ud_dd,QQb_us_ds,QQb_ud_uu)
      call addhel(6,2,1,5,7,3,4,QQb_du_dd,QQb_su_sd,QQb_uu_ud)
c--- basic Qb-Q amplitudes
c--- annihilation
      call addhel(6,5,2,1,7,3,4,RbR_dd_du,RbR_ss_du,RbR_uu_du)
      call addhel(2,1,6,5,7,3,4,RbR_ud_dd,RbR_ud_ss,RbR_ud_uu)
c--- scattering
      call addhel(2,5,6,1,7,3,4,QbQ_ud_dd,QbQ_us_ds,QbQ_ud_uu)
      call addhel(6,1,2,5,7,3,4,QbQ_du_dd,QbQ_su_sd,QbQ_uu_ud)
c--- basic Q-G amplitudes
      call addhel(7,6,1,5,2,3,4,QG_d_ddu,QG_d_dsc,QG_u_udu)
      call addhel(1,5,7,6,2,3,4,QG_u_ddd,QG_u_dcc,QG_u_duu)
c--- Notes added 28/5/09 (for checking u+g processes):
c---     QG_d_dsc corresponds to:     d + g -> d + c + sb
c---     QG_u_udu corresponds to:     d + g -> d + u + db
c---     QG_u_ddd corresponds to:     d + g -> u + u + ub
c---     QG_u_dcc corresponds to:     d + g -> u + c + cb

c--- basic QB-G amplitudes
      call addhel(6,7,5,1,2,3,4,QbG_d_ddu,QbG_u_ucs,QbG_u_uud)
      call addhel(5,1,6,7,2,3,4,QbG_d_udd,QbG_d_ucc,QbG_d_uuu)
c--- Notes added 28/5/09 (for checking db+g processes):
c---     QbG_u_ucs corresponds to:     ub + g -> ub + sb + c
c---     QbG_u_uud is NOT NEEDED:           db + g -> db + db + u
c---     QbG_d_ddu is used instead:    ub + g -> ub + db + u
c---     QbG_d_udd is NOT NEEDED:           ub + g -> db + ub + u
c---     QbG_d_uuu is used instead:    ub + g -> db + db + d
c---     QbG_d_ucc corresponds to:     ub + g -> db + sb + s

c--- basic G-Q amplitudes
      call addhel(7,5,2,6,1,3,4,GQ_d_ddu,GQ_d_dsc,GQ_u_udu)
      call addhel(2,6,7,5,1,3,4,GQ_u_ddd,GQ_u_dcc,GQ_u_duu)
c--- Notes added 28/5/09 (for checking g+u processes):
c---     GQ_d_dsc corresponds to:     g + d -> c + d + sb
c---     GQ_u_udu corresponds to:     g + d -> u + d + db
c---     GQ_u_ddd corresponds to:     g + d -> u + u + ub
c---     GQ_u_dcc corresponds to:     g + d -> c + u + cb

c--- basic G-QB amplitudes
      call addhel(5,7,6,2,1,3,4,GQb_d_udd,GQb_u_ucs,GQb_u_uud)
      call addhel(6,2,5,7,1,3,4,GQb_d_ddu,GQb_d_ucc,GQb_d_uuu)
c--- Notes added 28/5/09 (for checking g+db processes):
c---     GQb_u_ucs corresponds to:     g + ub -> sb + ub + c
c---     GQb_u_uud is NOT NEEDED     :     g + db -> db + db + u
c---     GQb_d_udd is used instead:    g + ub -> db + ub + u
c---     GQb_d_ddu is NOT NEEDED:          g + ub -> ub + db + u
c---     GQB_d_uuu is used instead:    g + ub -> db + db + d
c---     GQb_d_ucc corresponds to:     g + ub -> sb + db + s

      endif

c--- note the factor of 4._dp*xw**2 relative to wbb
      fac=4._dp*gsq**3*(gwsq/2._dp)**2
c--- extra factor of 2**3=8 to compensate for Ta normalization
      fac=fac*8._dp

      do j=-nf,nf
      do k=-nf,nf

      msq(j,k)=0._dp

      if (Gflag) then
************************************************************************
*     Sum the contributions from the QQGG matrix elements              *
************************************************************************

c--- note the identical particle factor of 1/6 for the
c--- q-qb initial states, due to 3 gluons in the final state
      if     ((j == 0) .and. (k == 0)) then
        Vfac=0._dp
        do n1=1,nf
          do n2=-nf,-1
            Vfac=Vfac+Vsq(n1,n2)
          enddo

        enddo
        msq(j,k)=propsq*mmsq_gg*Vfac*(gwsq**2/4._dp/esq**2)
      elseif ((j > 0) .and. (k < 0)) then
        msq(j,k)=propsq*mmsq_qqb*Vsq(j,k)
     &            *(aveqq/avegg)*(gwsq**2/4._dp/esq**2)/6._dp
      elseif ((j < 0) .and. (k > 0)) then
        msq(j,k)=propsq*mmsq_qbq*Vsq(j,k)
     &            *(aveqq/avegg)*(gwsq**2/4._dp/esq**2)/6._dp
      elseif ((j > 0) .and. (k == 0)) then
        msq(j,k)=half*propsq*mmsq_qg
     &            *(aveqg/avegg)*(gwsq**2/4._dp/esq**2)
     &            *(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))
      elseif ((j < 0) .and. (k == 0)) then
        msq(j,k)=half*propsq*mmsq_qbg
     &            *(aveqg/avegg)*(gwsq**2/4._dp/esq**2)
     &            *(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))
      elseif ((j == 0) .and. (k > 0)) then
        msq(j,k)=half*propsq*mmsq_gq
     &            *(aveqg/avegg)*(gwsq**2/4._dp/esq**2)
     &            *(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))
      elseif ((j == 0) .and. (k < 0)) then
        msq(j,k)=half*propsq*mmsq_gqb
     &            *(aveqg/avegg)*(gwsq**2/4._dp/esq**2)
     &            *(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))
      endif
      endif

      if (Qflag) then
************************************************************************
*     Sum the contributions from the QQBQQB matrix elements            *
************************************************************************

c--- note the factor of 4._dp*xw**2 relative to wbb
      fac=gsq**3*gwsq**2
c--- extra factor of 2**3=8 to compensate for Ta normalization
      fac=fac*8._dp

      if     ((j == 0) .and. (k == 0)) then
c--- no glue-glue contribution here
      elseif ((j > 0) .and. (k > 0)) then
c--- Q Q --> Q Q
        if ((jj(j) == 1) .and. (kk(k) == 2)) then
          msq(j,k)=msq(j,k)+Vsq(j,-k)*QQ_ud_dd*0.5_dp
     &           +(Vsum(j)-Vsq(j,-k))*QQ_us_ds
        elseif ((jj(j) == 2) .and. (kk(k) == 1)) then
          msq(j,k)=msq(j,k)+Vsq(k,-j)*QQ_du_dd*0.5_dp
     &           +(Vsum(k)-Vsq(k,-j))*QQ_su_sd
        elseif ((jj(j) == 1) .and. (kk(k) == 1)) then
          if (j == k) msq(j,k)=msq(j,k)+Vsum(j)*QQ_uu_du
          if (j .ne. k) msq(j,k)=msq(j,k)+Vsum(j)*QQ_us_ds
     &                                   +Vsum(k)*QQ_su_sd
        endif
      elseif ((j < 0) .and. (k < 0)) then
c--- Qb Qb --> Qb Qb
        if ((jj(j) == -1) .and. (kk(k) == -2)) then
          msQ(j,k)=msq(j,k)+Vsq(k,-j)*QbQb_ud_uu*0.5_dp
     &           +(Vsum(k)-Vsq(k,-j))*QbQb_sd_su
        elseif ((jj(j) == -2) .and. (kk(k) == -1)) then
          msq(j,k)=msq(j,k)+Vsq(j,-k)*QbQb_du_uu*0.5_dp
     &           +(Vsum(j)-Vsq(j,-k))*QbQb_ds_us
        elseif ((jj(j) == -2) .and. (kk(k) == -2)) then
          if (j == k) msq(j,k)=msq(j,k)+Vsum(j)*QbQb_dd_ud
          if (j .ne. k) msq(j,k)=msq(j,k)+Vsum(j)*QbQb_ds_us
     &                                   +Vsum(k)*QbQb_sd_su
        endif
      elseif ((j > 0) .and. (k < 0)) then
c--- Q Qb --> Q Qb
        if ((jj(j)==1) .and. (kk(k)==-2)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*(RRb_ud_dd+RRb_ud_uu
     &                               +real(nf-2,dp)*RRb_ud_ss)
     &           +(Vsum(j)-Vsq(j,k))*QQb_us_ds
     &           +(Vsum(k)-Vsq(j,k))*QQb_su_sd
        elseif ((jj(j)==1) .and. (kk(k)==-1)) then
          if (j==-k) then
            Vfac=0._dp
            do n1=1,nf
            do n2=-nf,-1
            if ((n1 .ne. j) .and. (n2 .ne. k)) then
              Vfac=Vfac+Vsq(n1,n2)
            endif
            enddo
            enddo
            msq(j,k)=Vsum(j)*RRb_uu_du
     &                 +Vfac*RRb_ss_du
          else
            msq(j,k)=msq(j,k)+Vsum(j)*QQb_us_ds
          endif
        elseif ((jj(j)==2) .and. (kk(k)==-2)) then
          if (j==-k) then
            Vfac=0._dp
            do n1=1,nf
            do n2=-nf,-1
            if ((n1 .ne. j) .and. (n2 .ne. k)) then
              Vfac=Vfac+Vsq(n1,n2)
            endif
            enddo
            enddo
            msq(j,k)=Vsum(k)*RRb_dd_du
     &                 +Vfac*RRb_ss_du
          else
            msq(j,k)=msq(j,k)+Vsum(k)*QQb_su_sd
          endif
        endif
      elseif ((j < 0) .and. (k > 0)) then
c--- Qb Q --> Q Qb
        if ((jj(j)==-2) .and. (kk(k)==1)) then
          msq(j,k)=msq(j,k)+Vsq(j,k)*(RbR_ud_dd+RbR_ud_uu
     &                               +real(nf-2,dp)*RbR_ud_ss)
     &           +(Vsum(k)-Vsq(j,k))*QbQ_us_ds
     &           +(Vsum(j)-Vsq(j,k))*QbQ_su_sd
        elseif ((jj(j)==-1) .and. (kk(k)==1)) then
          if (j==-k) then
            Vfac=0._dp
            do n1=-nf,-1
            do n2=1,nf
            if ((n1 .ne. j) .and. (n2 .ne. k)) then
              Vfac=Vfac+Vsq(n1,n2)
            endif
            enddo
            enddo
            msq(j,k)=Vsum(k)*RbR_uu_du
     &                 +Vfac*RbR_ss_du
          else
            msq(j,k)=msq(j,k)+Vsum(k)*QbQ_us_ds
          endif
        elseif ((jj(j)==-2) .and. (kk(k)==2)) then
          if (j==-k) then
            Vfac=0._dp
            do n1=-nf,-1
            do n2=1,nf
            if ((n1 .ne. j) .and. (n2 .ne. k)) then
              Vfac=Vfac+Vsq(n1,n2)
            endif
            enddo
            enddo
            msq(j,k)=Vsum(j)*RbR_dd_du
     &                 +Vfac*RbR_ss_du
          else
            msq(j,k)=msq(j,k)+Vsum(j)*QbQ_su_sd
          endif
        endif
      elseif ((j > 0) .and. (k == 0)) then
c--- Q G --> Q q qb
        if( jj(j) == 1) then
          if (j== 5) then
            Vfac=0._dp
          else
            Vfac=Vsq(j,-j-1)
          endif
              msq(j,k)=msq(j,k)+Vfac*(
     &                   (aveqg/aveqq)*QG_u_udu
     &                  +real(nf-2,dp)*(aveqg/aveqq)*QG_u_dcc
     &                  +(aveqg/aveqq)*QG_u_ddd*0.5_dp)
     &                         +(2._dp-Vfac)*
     &                   (aveqg/aveqq)*QG_d_dsc
        elseif( jj(j) == 2) then
         msq(j,k)=msq(j,k)+Vsq(-j,j-1)*(aveqg/aveqq)*QG_d_ddu*0.5_dp
     &                    +(2._dp-Vsq(-j,j-1))*(aveqg/aveqq)*QG_d_dsc
        endif
      elseif ((j < 0) .and. (k == 0)) then
c--- QB G --> QB qb q
        if( jj(j) == -2) then
         msq(j,k)=msq(j,k)+Vsq(j,-j-1)*(
     &              (aveqg/aveqq)*QbG_d_ddu
     &             +real(nf-2,dp)*(aveqg/aveqq)*QbG_d_ucc
     &             +(aveqg/aveqq)*QbG_d_uuu*0.5_dp)
     &                    +(2._dp-Vsq(j,-j-1))*
     &              (aveqg/aveqq)*QbG_u_ucs
        elseif( jj(j) == -1) then
          if (j== -5) then
            Vfac=0._dp
          else
            Vfac=Vsq(-j,j-1)
          endif
         msq(j,k)=msq(j,k)+Vfac*(aveqg/aveqq)*QbG_u_uud*0.5_dp
     &                    +(2._dp-Vfac)*(aveqg/aveqq)*QbG_u_ucs
        endif
c--- G Q --> Q q qb
      elseif ((j == 0) .and. (k > 0)) then
        if( kk(k) == 1) then
         if (k == 5) then
           Vfac=0._dp
         else
           Vfac=Vsq(k,-k-1)
         endif
         msq(j,k)=msq(j,k)+Vfac*(
     &             +real(nf-2,dp)*(aveqg/aveqq)*GQ_u_dcc
     &             +(aveqg/aveqq)*GQ_u_ddd*0.5_dp
     &             +(aveqg/aveqq)*GQ_u_udu)
     &                         +(2._dp-Vfac)*
     &            (aveqg/aveqq)*GQ_d_dsc
        elseif( kk(k) == 2) then
         msq(j,k)=msq(j,k)+Vsq(-k,k-1)*(aveqg/aveqq)*GQ_d_ddu*0.5_dp
     &                    +(2._dp-Vsq(-k,k-1))*(aveqg/aveqq)*GQ_d_dsc
        endif
      elseif ((j == 0) .and. (k < 0)) then
c--- G QB --> QB qb q
        if( kk(k) == -2) then
              msq(j,k)=msq(j,k)+Vsq(k,-k-1)*(
     &                   (aveqg/aveqq)*GQb_d_udd
     &                  +real(nf-2,dp)*(aveqg/aveqq)*GQb_d_ucc
     &                  +(aveqg/aveqq)*GQb_d_uuu*0.5_dp)
     &                         +(2._dp-Vsq(k,-k-1))*
     &                   (aveqg/aveqq)*GQb_u_ucs
        elseif( kk(k) == -1) then
         if (k == -5) then
           Vfac=0._dp
         else
           Vfac=Vsq(-k,k-1)
         endif
         msq(j,k)=msq(j,k)+Vfac*(aveqg/aveqq)*GQb_u_uud*0.5_dp
     &                    +(2._dp-Vfac)*(aveqg/aveqq)*GQb_u_ucs
        endif
      endif

      endif

      enddo
      enddo

      return
      end
