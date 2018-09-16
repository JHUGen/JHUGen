      subroutine qqb_w2jet_v(p,msqv)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     January 2001.                                                    *
*     Additions Aug. 2001, for the 4Q piece                            *
*      (singularities cancelled)                                       *
*                                                                      *
*     Calculate the virtual matrix element squared and                 *
*     subtraction terms for the process                                *
*                                                                      *
*     q(-p1) + qbar(-p2) --> W + j(p5) + j(p6)                         *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
*                                                                      *
*     where the partons are either q(p5) and qbar(p6) [Qflag = .true.] *
*                               or g(p5) and g(p6)    [Gflag = .true.] *
*                                                                      *
*     The value of COLOURCHOICE determines which colour structures     *
*     are included in the terms for the QQGG piece                     *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'flags.f'
      include 'msq_cs.f'
      include 'lc.f'
      include 'noglue.f'
      include 'first.f'
      include 'ppmax.f'
      include 'mpicommon.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),fac,
     & mmsq_qqb,mmsq_qbq,mmsq_gq,mmsq_gqb,mmsq_qg,mmsq_qbg,mmsq_gg,
     & p(mxpart,4),q(mxpart,4),pswap(mxpart,4)
      complex(dp):: prop
      integer:: nu,j,k,n1,n2,cs,rvcolourchoice
      real(dp):: subuv(0:2)
      real(dp):: Vfac
      real(dp):: mqq(0:2,fn:nf,fn:nf)
      real(dp):: ppmsqx(0:2,ppmax)
      real(dp):: msqx_cs(0:2,-nf:nf,-nf:nf)
      real(dp):: qqb_ijkk,qqb_iikl,qqb_ijkj,qqb_ijik,
     &                 qqb_ijii,qqb_ijjj,qqb_iiij,qqb_iiji,
     &                 qbq_ijkk,qbq_iikl,qbq_ijkj,qbq_ijik,
     &                 qbq_ijii,qbq_ijjj,qbq_iiij,qbq_iiji,
     &                 qq_ijkk,qq_iikl,qq_ijkj,qq_ijik,
     &                 qq_ijii,qq_ijjj,qq_iiij,qq_iiji,
     &                 qbqb_ijkk,qbqb_iikl,qbqb_ijkj,qbqb_ijik,
     &                 qbqb_ijii,qbqb_ijjj,qbqb_iiij,qbqb_iiji
c      character*30 runstring
c      common/runstring/runstring
      common/mqq/mqq
      common/rvcolourchoice/rvcolourchoice
!$omp threadprivate(/mqq/,/rvcolourchoice/)

      if (first) then
        first=.false.
        if (rank == 0) then
        if ((Gflag) .or. (QandGflag)) then
          write(*,*) 'Using QQGG (VIRTUAL) matrix elements'
!          write(*,*) '[LC is     N   ]'
!          write(*,*) '[SLC is   1/N  ]'
!          write(*,*) '[SSLC is 1/N**3]'
        endif
        if ((Qflag) .or. (QandGflag)) then
          write(*,*) 'Using QQBQQB (VIRTUAL) matrix elements'
!          write(*,*) '[LC is   1 ]'
!          write(*,*) '[SLC is 1/N]'
        endif
        if     (rvcolourchoice == 1) then
          write(*,*) 'Leading colour only in VIRTUAL'
        elseif (rvcolourchoice == 2) then
          write(*,*) 'Sub-leading colour only in VIRTUAL'
        elseif (rvcolourchoice == 3) then
          write(*,*) 'Sub-sub-leading colour only in VIRTUAL'
        elseif (rvcolourchoice == 0) then
          write(*,*) 'Total of all colour structures in VIRTUAL'
        else
          write(*,*) 'Bad colourchoice'
          stop
        endif
        endif
      endif

      scheme='dred'

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

c--- calculate the lowest order matrix element and fill the
c--- common block twopij with s_{ij}
      call qqb_wp2jetx_new(p,msq,mqq,ppmsqx,msqx_cs)

      prop=s(3,4)/cplx2(s(3,4)-wmass**2,wmass*wwidth)

************************************************************************
*     Contributions from QQGG matrix elements                          *
************************************************************************
      if (Gflag) then
c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      if     (colourchoice == 1) then
        subuv(1)=2._dp*xn
     &  *(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
        subuv(2)=subuv(1)
      elseif (colourchoice == 2) then
        subuv(0)=2._dp*xn
     &  *(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
      elseif (colourchoice == 3) then
c--- all zero already
      elseif (colourchoice == 0) then
        subuv(1)=2._dp*xn
     &  *(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
        subuv(2)=subuv(1)
        subuv(0)=subuv(1)
      endif

c--- transfer lowest order matrix elements
c--- NB: this breaks the routine if Qflag = Gflag = .true.

      do cs=0,2
        do j=-nf,nf
        do k=-nf,nf
        msq_cs(cs,j,k)=msqx_cs(cs,j,k)
        enddo
        enddo
      enddo

c--- when testing alpha-dependence, we do not need loop contribution
c      if (runstring(1:5) == 'alpha') then
c        return
c      endif

c--- Now calculate the relevant lowest-order matrix elements
c--- for each possible initial state from the QQGG contribution

c---  calculate the qqb terms
CALL    0--> q(p2)+g(p5)+g(p6)+qb(p1)+l(p3)+lbar(p4)
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
      do nu=1,4
      pswap(1,nu)=p(2,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(5,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qqb)

c---  calculate the qbq terms
CALL    0--> q(p1)+g(p5)+g(p6)+qb(p2)+l(p3)+lbar(p4)
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(1,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(5,nu)
      pswap(4,nu)=p(2,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qbq)

c---  calculate the gq terms
CALL    0--> q(p5)+g(p1)+g(p6)+qb(p2)+l(p3)+lbar(p4)
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(1,nu)
      pswap(4,nu)=p(2,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_gq)

c---  calculate the qg terms
CALL    0--> q(p5)+g(p2)+g(p6)+qb(p1)+l(p3)+lbar(p4)
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qg)

c---  calculate the gqb terms
CALL    0--> q(p2)+g(p1)+g(p6)+qb(p5)+l(p3)+lbar(p4)
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(2,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(1,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_gqb)

c---  calculate the qbg terms
CALL    0--> q(p1)+g(p2)+g(p6)+qb(p5)+l(p3)+lbar(p4)
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l+(p6)
      do nu=1,4
      pswap(1,nu)=p(1,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_qbg)

c--- calculate the gg terms
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(6,nu)
      pswap(5,nu)=p(3,nu)
      pswap(6,nu)=p(4,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xwqqgg_v(mmsq_gg)
      endif

************************************************************************
*     Contributions from QQQQ matrix elements                          *
************************************************************************
      if (Qflag) then
c--- UV counter-term is already included in a6routine.f
      subuv(1)=0._dp
      subuv(2)=subuv(1)
      subuv(0)=subuv(1)

c--- transfer lowest order matrix elements
c--- NB: this breaks the routine if Qflag = Gflag = .true.

      do cs=0,2
        do j=-nf,nf
        do k=-nf,nf
        msq_cs(cs,j,k)=mqq(cs,j,k)
        enddo
        enddo
      enddo

c--- when testing alpha-dependence, we do not need loop contribution
c      if (runstring(1:5) == 'alpha') then
c        return
c      endif

c--- early return if there's no contribution here
      if ((gqonly) .or. (ggonly)) return

c---  Now transform momenta into a notation
c---  suitable for calling the BDKW function with notation which is
c---  q-(-p4)+Q+(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
      do nu=1,4
      q(1,nu)=p(2,nu)
      q(2,nu)=p(6,nu)
      q(3,nu)=p(5,nu)
      q(4,nu)=p(1,nu)
      q(5,nu)=p(4,nu)
      q(6,nu)=p(3,nu)
      enddo
      call spinoru(6,q,za,zb)
      fac=V*xn*gw**4*gsq**2*ason2pi

c--- set-up qqb matrix elements
      call qqbw2j_loop(1,2,3,4,5,6,qqb_ijkk,qqb_iikl,qqb_ijkj,qqb_ijik,
     &                             qqb_ijii,qqb_ijjj,qqb_iiij,qqb_iiji)
c--- qbq
      call qqbw2j_loop(4,2,3,1,5,6,qbq_ijkk,qbq_iikl,qbq_ijkj,qbq_ijik,
     &                             qbq_ijii,qbq_ijjj,qbq_iiij,qbq_iiji)
c--- qq (note that roles of iiij and ijii are reversed)
      call qqbw2j_loop(2,1,3,4,5,6,qq_ijkk,qq_iikl,qq_ijkj,qq_ijik,
     &                             qq_iiij,qq_ijjj,qq_ijii,qq_iiji)
c--- qbqb (note that roles of ijjj and iiji are reversed)
      call qqbw2j_loop(1,2,4,3,5,6,qbqb_ijkk,qbqb_iikl,qbqb_ijkj,
     &         qbqb_ijik,qbqb_ijii,qbqb_iiji,qbqb_iiij,qbqb_ijjj)

      endif

c--- Add VIRTUAL terms

      do j=-nf,nf
      do k=-nf,nf

************************************************************************
*     Contributions from QQGG matrix elements                          *
************************************************************************
      if (Gflag) then
      if     ((j > 0) .and. (k < 0)) then
        msqv(j,k)=msqv(j,k)+Vsq(j,k)*mmsq_qqb*abs(prop)**2
     &           *half*(aveqq/avegg)*(gwsq**2/4._dp/esq**2)
      elseif ((j < 0) .and. (k > 0)) then
        msqv(j,k)=msqv(j,k)+Vsq(j,k)*mmsq_qbq*abs(prop)**2
     &           *half*(aveqq/avegg)*(gwsq**2/4._dp/esq**2)
      elseif ((j > 0) .and. (k == 0)) then
        msqv(j,k)=msqv(j,k)+mmsq_qg*abs(prop)**2*
     &           (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))
     &           *(aveqg/avegg)*(gwsq**2/4._dp/esq**2)
      elseif ((j < 0) .and. (k == 0)) then
        msqv(j,k)=msqv(j,k)+mmsq_qbg*abs(prop)**2*
     &           (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))
     &           *(aveqg/avegg)*(gwsq**2/4._dp/esq**2)
      elseif ((j == 0) .and. (k > 0)) then
        msqv(j,k)=msqv(j,k)+mmsq_gq*abs(prop)**2*
     &           (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))
     &           *(aveqg/avegg)*(gwsq**2/4._dp/esq**2)
      elseif ((j == 0) .and. (k < 0)) then
        msqv(j,k)=msqv(j,k)+mmsq_gqb*abs(prop)**2*
     &           (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))
     &           *(aveqg/avegg)*(gwsq**2/4._dp/esq**2)
      elseif ((j == 0) .and. (k == 0)) then
        Vfac=0._dp
        do n1=1,nf
          do n2=-nf,-1
            Vfac=Vfac+Vsq(n1,n2)
          enddo
        enddo
        msqv(j,k)=msqv(j,k)
     &           +mmsq_gg*Vfac*abs(prop)**2*(gwsq**2/4._dp/esq**2)
      endif
      endif

      if (Qflag) then
      if     ((j > 0) .and. (k < 0)) then
        if (j .ne. -k) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &      Vsq(j,k)*(qqb_ijii+qqb_ijjj)
     &     +Vsq(j,k)*real(nf-2,dp)*qqb_ijkk
     &     +(Vsum(j)-Vsq(j,k))*qqb_ijkj
     &     +(Vsum(k)-Vsq(j,k))*qqb_ijik)
        else
          Vfac=0._dp
          do n1=1,nf
          do n2=-nf,-1
          if ((n1 .ne. j) .and. (n2 .ne. k)) then
            Vfac=Vfac+Vsq(n1,n2)
          endif
          enddo
          enddo
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &      Vsum(k)*qqb_iiij
     &     +Vsum(j)*qqb_iiji
     &     +Vfac*qqb_iikl)
        endif
      elseif ((j < 0) .and. (k > 0)) then
        if (j .ne. -k) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &      Vsq(j,k)*(qbq_ijii+qbq_ijjj)
     &     +Vsq(j,k)*real(nf-2,dp)*qbq_ijkk
     &     +(Vsum(k)-Vsq(j,k))*qbq_ijkj
     &     +(Vsum(j)-Vsq(j,k))*qbq_ijik)
        else
          Vfac=0._dp
          do n1=-nf,-1
          do n2=1,nf
          if ((n1 .ne. j) .and. (n2 .ne. k)) then
            Vfac=Vfac+Vsq(n1,n2)
          endif
          enddo
          enddo
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &      Vsum(j)*qbq_iiij
     &     +Vsum(k)*qbq_iiji
     &     +Vfac*qbq_iikl)
        endif
      elseif ((j > 0) .and. (k > 0)) then
        if (j .ne. k) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &      Vsq(j,-k)*half*qq_ijjj
     &     +Vsq(k,-j)*half*qq_ijii
     &     +(Vsum(j)-Vsq(j,-k))*qq_ijkj
     &     +(Vsum(k)-Vsq(k,-j))*qq_ijik)
        else
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &      Vsum(j)*qq_iiji)
        endif
      elseif ((j < 0) .and. (k < 0)) then
        if (j .ne. k) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &      Vsq(j,-k)*half*qbqb_ijjj
     &     +Vsq(k,-j)*half*qbqb_ijii
     &     +(Vsum(j)-Vsq(j,-k))*qbqb_ijkj
     &     +(Vsum(k)-Vsq(k,-j))*qbqb_ijik)
        else
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &      Vsum(j)*qbqb_iiji)
        endif
      endif
      endif

      enddo
      enddo

************************************************************************
*     UV contributions are included here                               *
************************************************************************
      do j=-nf,nf
      do k=-nf,nf

      do cs=0,2
      msqv(j,k)=msqv(j,k)+
     &  ason2pi*(-subuv(cs))*msq_cs(cs,j,k)
      enddo

      enddo
      enddo

      return
      end



