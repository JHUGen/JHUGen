      subroutine qqb_z2jet_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2001.                                                     *
*     Tested by JMC, 6/7/01 and found to agree with Madgraph           *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  Z^0 + f(p5)+f(p6)+g(p7)
c                           |
c                            --> e^-(p3)+e^+(p4)
c

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'flags.f'
      include 'lc.f'
      include 'first.f'
      include 'mpicommon.f'
      integer:: j,k,nquark
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: mmsq_gg(2,2),mmsq_qqb(2,2),mmsq_qbq(2,2),
     & mmsq_qg(2,2),mmsq_gq(2,2),mmsq_gqb(2,2),mmsq_qbg(2,2)
      real(dp):: fac
      real(dp):: msqi_qq(2),msqi_qbqb(2),
     &                 msqi_qqb(2),msqi_qbq(2),
     &                 msqi_qqbs(2),msqi_qbqs(2),
     &                 msqi_qg(2),msqi_qbg(2),
     &                 msqi_gqb(2),msqi_gq(2),ggtemp
      real(dp):: msqn_qq(2,2),msqn_qbqb(2,2),
     &                 msqn_qqb(2,2),msqn_qbq(2,2),
     &                 msqn_qqbs(2,2),msqn_qbqs(2,2),
     &                 msqn_qg(2,2),msqn_qbg(2,2),
     &                 msqn_gqb(2,2),msqn_gq(2,2)
      complex(dp):: prop
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      include 'cplx.h'
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

C---call spinor routine and load common block twopij
      call spinoru(7,p,za,zb)

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

      if (Gflag) then
************************************************************************
*     Calculate contributions from the QQGGG matrix elements            *
************************************************************************

c-- matrix elements for gg -> qbq
      call xzqqggg(5,1,2,7,6,3,4,mmsq_gg)
c-- matrix elements for qqb -> gg
      call xzqqggg(1,5,6,7,2,3,4,mmsq_qqb)
c-- matrix elements for qbq -> gg
      call xzqqggg(2,5,6,7,1,3,4,mmsq_qbq)
c-- matrix elements for qg -> qg
      call xzqqggg(1,2,6,7,5,3,4,mmsq_qg)
c-- matrix elements for gq -> gq
      call xzqqggg(2,1,6,7,5,3,4,mmsq_gq)
c-- matrix elements for gqb -> gqb
      call xzqqggg(5,1,6,7,2,3,4,mmsq_gqb)
c-- matrix elements for qbg -> qbg
      call xzqqggg(5,2,6,7,1,3,4,mmsq_qbg)

      do j=-nf,nf
      do k=-nf,nf

      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
           ggtemp=0._dp
           do nquark=1,nf

           ggtemp=ggtemp
     &            +abs(Q(nquark)*q1+L(nquark)*l1*prop)**2*mmsq_gg(1,1)
     &            +abs(Q(nquark)*q1+R(nquark)*r1*prop)**2*mmsq_gg(2,2)
     &            +abs(Q(nquark)*q1+L(nquark)*r1*prop)**2*mmsq_gg(1,2)
     &            +abs(Q(nquark)*q1+R(nquark)*l1*prop)**2*mmsq_gg(2,1)
          enddo
          msq(j,k)=ggtemp
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=+abs(Q(j)*q1+L(j)*l1*prop)**2*mmsq_qqb(1,1)
     &             +abs(Q(j)*q1+R(j)*r1*prop)**2*mmsq_qqb(2,2)
     &             +abs(Q(j)*q1+L(j)*r1*prop)**2*mmsq_qqb(1,2)
     &             +abs(Q(j)*q1+R(j)*l1*prop)**2*mmsq_qqb(2,1)

c---Statistical factor
          msq(j,k)=aveqq/avegg*msq(j,k)/6._dp
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=+abs(Q(k)*q1+L(k)*l1*prop)**2*mmsq_qbq(1,1)
     &             +abs(Q(k)*q1+R(k)*r1*prop)**2*mmsq_qbq(2,2)
     &             +abs(Q(k)*q1+L(k)*r1*prop)**2*mmsq_qbq(1,2)
     &             +abs(Q(k)*q1+R(k)*l1*prop)**2*mmsq_qbq(2,1)
c---Statistical factor
          msq(j,k)=aveqq/avegg*msq(j,k)/6._dp
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=+abs(Q(j)*q1+L(j)*l1*prop)**2*mmsq_qg(1,1)
     &             +abs(Q(j)*q1+R(j)*r1*prop)**2*mmsq_qg(2,2)
     &             +abs(Q(j)*q1+L(j)*r1*prop)**2*mmsq_qg(1,2)
     &             +abs(Q(j)*q1+R(j)*l1*prop)**2*mmsq_qg(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=+abs(Q(-j)*q1+L(-j)*l1*prop)**2*mmsq_qbg(1,1)
     &             +abs(Q(-j)*q1+R(-j)*r1*prop)**2*mmsq_qbg(2,2)
     &             +abs(Q(-j)*q1+L(-j)*r1*prop)**2*mmsq_qbg(1,2)
     &             +abs(Q(-j)*q1+R(-j)*l1*prop)**2*mmsq_qbg(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=+abs(Q(k)*q1+L(k)*l1*prop)**2*mmsq_gq(1,1)
     &             +abs(Q(k)*q1+R(k)*r1*prop)**2*mmsq_gq(2,2)
     &             +abs(Q(k)*q1+L(k)*r1*prop)**2*mmsq_gq(1,2)
     &             +abs(Q(k)*q1+R(k)*l1*prop)**2*mmsq_gq(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=+abs(Q(-k)*q1+L(-k)*l1*prop)**2*mmsq_gqb(1,1)
     &             +abs(Q(-k)*q1+R(-k)*r1*prop)**2*mmsq_gqb(2,2)
     &             +abs(Q(-k)*q1+L(-k)*r1*prop)**2*mmsq_gqb(1,2)
     &             +abs(Q(-k)*q1+R(-k)*l1*prop)**2*mmsq_gqb(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
      endif
   19 continue
      enddo
      enddo
      endif

      if (Qflag) then
c--- note the factor of 4._dp*xw**2 relative to wbb
      fac=4._dp*gsq**3*esq**2*8._dp
c--- extra factor of 2**3=8 to compensate for Ta normalization

c--- note: the following two arrays end up being overall 1<->2 symmetric
      call msq_ZqqQQg(5,1,2,6,7,4,3,msqn_qqbs,msqi_qqbs)
      call msq_ZqqQQg(1,6,5,2,7,4,3,msqn_qbqs,msqi_qbqs)

      call msq_ZqqQQg(2,1,5,6,7,4,3,msqn_qqb,msqi_qqb)
      call msq_ZqqQQg(1,2,5,6,7,4,3,msqn_qbq,msqi_qbq)

      call msq_ZqqQQg(1,5,2,6,7,4,3,msqn_qbqb,msqi_qbqb)
      call msq_ZqqQQg(5,1,6,2,7,4,3,msqn_qq,msqi_qq)


c      call msq_ZqqQQg(5,1,6,7,2,4,3,msqn_qg,msqi_qg)
      call msq_ZqqQQg(7,1,5,6,2,4,3,msqn_qg,msqi_qg)
      call msq_ZqqQQg(2,7,5,6,1,4,3,msqn_gqb,msqi_gqb)
      call msq_ZqqQQg(7,2,5,6,1,4,3,msqn_gq,msqi_gq)
      call msq_ZqqQQg(1,7,5,6,2,4,3,msqn_qbg,msqi_qbg)

      do j=-nf,nf
      do k=-nf,nf

      if ((j > 0) .and. (k < 0)) then
c-qqb
           if (k==-j) then
            msq(j,k)=msq(j,k)+fac*aveqq*(msqi_qqbs(jj(j))
     &      +real(1+jj(j),dp)*msqn_qqb(jj(j),1)
     &      +real(3-jj(j),dp)*msqn_qqb(jj(j),2))
           else
            msq(j,k)=msq(j,k)+fac*aveqq*msqn_qqbs(jj(j),-kk(k))
           endif
      elseif ((j < 0) .and. (k > 0)) then
c-qbq
          if (j ==-k) then
            msq(j,k)=msq(j,k)+fac*aveqq*(msqi_qbqs(kk(k))
     &      +real(1+kk(k),dp)*msqn_qbq(kk(k),1)
     &      +real(3-kk(k),dp)*msqn_qbq(kk(k),2))
          else
           msq(j,k)=msq(j,k)+fac*aveqq*msqn_qbqs(-jj(j),kk(k))
          endif
      elseif ((j > 0) .and. (k == 0)) then
c-qg
          msq(j,k)=msq(j,k)
     &     +fac*aveqg*(half*msqi_qg(jj(j))
     &      +real(1+jj(j),dp)*msqn_qg(jj(j),1)
     &      +real(3-jj(j),dp)*msqn_qg(jj(j),2)
     & )
c-qbg
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=msq(j,k)
     &     +fac*aveqg*(half*msqi_qbg(-jj(j))
     &      +real(1-jj(j),dp)*msqn_qbg(-jj(j),1)
     &      +real(3+jj(j),dp)*msqn_qbg(-jj(j),2)
     & )
      elseif ((j == 0) .and. (k > 0)) then
c-gq
          msq(j,k)=msq(j,k)
     &      +fac*aveqg*(half*msqi_gq(kk(k))
     &      +real(1+kk(k),dp)*msqn_gq(kk(k),1)
     &      +real(3-kk(k),dp)*msqn_gq(kk(k),2)
     & )
      elseif ((j == 0) .and. (k < 0)) then
c-gqb
          msq(j,k)=msq(j,k)
     &      +fac*aveqg*(half*msqi_gqb(-kk(k))
     &      +real(1-kk(k),dp)*msqn_gqb(-kk(k),1)
     &      +real(3+kk(k),dp)*msqn_gqb(-kk(k),2)
     & )
      elseif ((j > 0) .and. (k > 0)) then
c-qq
          if (j==k) then
          msq(j,k)=msq(j,k)+half*fac*aveqq*msqi_qq(jj(j))
          else
          msq(j,k)=msq(j,k)+fac*aveqq*msqn_qq(jj(j),kk(k))
         endif
      elseif ((j < 0) .and. (k < 0)) then
c-qbqb
          if (j==k) then
          msq(j,k)=msq(j,k)+half*fac*aveqq*msqi_qbqb(-jj(j))
          else
          msq(j,k)=msq(j,k)+fac*aveqq*msqn_qbqb(-jj(j),-kk(k))
          endif
      endif


      enddo
      enddo

      endif


      return
      end
