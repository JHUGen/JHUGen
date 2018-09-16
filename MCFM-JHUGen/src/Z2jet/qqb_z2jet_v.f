      subroutine qqb_z2jet_v(p,msqv)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     September 2001.                                                  *
*                                                                      *
*     Calculate the virtual matrix element squared and                 *
*     subtraction terms for the process                                *
*                                                                      *
*     q(-p1) + qbar(-p2) --> Z + j(p5) + j(p6)                         *
*                            |                                         *
*                            --> e^-(p3) + e^+(p4)                     *
*                                                                      *
*     where the partons are either q(p5) and qbar(p6) [Qflag = .true.] *
*                               or g(p5) and g(p6)    [Gflag = .true.] *
*                                                                      *
*     The value of COLOURCHOICE determines which colour structures     *
*     are included in the terms for the QQGG piece                     *
*                                                                      *
************************************************************************
*                                                                      *
*     By setting compare=.true. it is possible to test the squared     *
*     tree amplitudes in this routine with those from qqb_z2jet.f      *
*      (checked on 10/22/01 by JMC)                                    *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'flags.f'
      include 'msq_cs.f'
      include 'lc.f'
      include 'first.f'
      include 'ppmax.f'
      include 'mpicommon.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     & mmsq_qqb(2,2),mmsq_qbq(2,2),mmsq_gq(2,2),mmsq_qg(2,2),
     & mmsq_qbg(2,2),mmsq_gqb(2,2),mmsq_gg(2,2),
     & p(mxpart,4),pswap(mxpart,4),fac
      complex(dp):: mmsq_qqb_ax(2,2),mmsq_qbq_ax(2,2),
     & mmsq_gq_ax(2,2),mmsq_qg_ax(2,2),mmsq_qbg_ax(2,2),
     & mmsq_gqb_ax(2,2),mmsq_gg_ax(2,2),
     & mmsq_qqb_vec(2,2),mmsq_qbq_vec(2,2),
     & mmsq_gq_vec(2,2),mmsq_qg_vec(2,2),mmsq_qbg_vec(2,2),
     & mmsq_gqb_vec(2,2),mmsq_gg_vec(2,2)
      complex(dp):: atreez,a61z,a62z,a63z,prop,vcouple(2),
     & tamp,lamp,tampup,lampup,tampdo,lampdo,tamps,lamps,lampx,lampsx
      integer:: nu,j,k,cs,polq,polb,polz
      integer:: nup,ndo,rvcolourchoice
      integer, parameter :: jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      real(dp):: subuv(0:2)
      real(dp):: faclo,v2(2),vQ(nf,2),idfac
      real(dp):: mqq(0:2,fn:nf,fn:nf)
      real(dp):: ppmsqx(0:2,ppmax)
      real(dp):: msqx_cs(0:2,-nf:nf,-nf:nf)
      complex(dp):: atreez_526143(2,2,2),atreez_251643(2,2,2),
     &               atreez_625143(2,2,2),atreez_261543(2,2,2),
     &               atreez_162543(2,2,2),atreez_615243(2,2,2),
     &               atreez_152643(2,2,2),atreez_516243(2,2,2),
     &               atreez_562143(2,2,2),atreez_651243(2,2,2),
     &               atreez_265143(2,2,2),atreez_621543(2,2,2),
     &               atreez_561243(2,2,2),atreez_652143(2,2,2),
     &               atreez_165243(2,2,2),atreez_612543(2,2,2),
     &               a61z_526143(2,2,2),a61z_251643(2,2,2),
     &               a61z_625143(2,2,2),a61z_261543(2,2,2),
     &               a61z_162543(2,2,2),a61z_615243(2,2,2),
     &               a61z_152643(2,2,2),a61z_516243(2,2,2),
     &               a61z_562143(2,2,2),a61z_651243(2,2,2),
     &               a61z_265143(2,2,2),a61z_621543(2,2,2),
     &               a61z_561243(2,2,2),a61z_652143(2,2,2),
     &               a61z_165243(2,2,2),a61z_612543(2,2,2),
     &               a62z_625143(2,2,2),a62z_261543(2,2,2),
     &               a62z_526143(2,2,2),a62z_251643(2,2,2),
     &               a62z_152643(2,2,2),a62z_516243(2,2,2),
     &               a62z_162543(2,2,2),a62z_615243(2,2,2),
     &               a62z_265143(2,2,2),a62z_621543(2,2,2),
     &               a62z_562143(2,2,2),a62z_651243(2,2,2),
     &               a62z_165243(2,2,2),a62z_612543(2,2,2),
     &               a62z_561243(2,2,2),a62z_652143(2,2,2),
     &               a63z_526143(2,2,2),a63z_625143(2,2,2),
     &               a63z_162543(2,2,2),a63z_152643(2,2,2),
     &               a63z_562143(2,2,2),a63z_265143(2,2,2),
     &               a63z_561243(2,2,2),a63z_165243(2,2,2)

      logical:: compare,checkvector,checkaxial
      common/mqq/mqq
      common/rvcolourchoice/rvcolourchoice
      parameter (nup=2,ndo=nf-nup)
      parameter (compare=.false.)
c--- These parameters allow for a point-by-point comparison of the
c--- vector and axial pieces with MadLoop. They should normally
c--- both be set to .false.
      parameter(checkvector=.false.,checkaxial=.false.)
! discardvecax: if true, discard diagrams where Z couples directly
! to a closed quark loop, through vector and axial-vector couplings
      logical, parameter:: discardvecax=.true.
!$omp threadprivate(/mqq/,/rvcolourchoice/)
      include 'cplx.h'

      if (Qflag .and. Gflag) then
        write(6,*) 'Both Qflag and Gflag cannot be true'
        write(6,*) 'They are set in file options.DAT'
        write(6,*) 'Failed in qqb_z2jet_v.f'
        stop
      endif

      if (first) then
        first=.false.
        if (rank == 0) then
        if (discardvecax) then
          write(6,*) 'No vector or axial-vector closed quark-loop Z+2 jet diagrams'
          write(6,*)
        endif
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

c--- get special phase space point for MadLoop comparison
      if (checkvector .or. checkaxial) then
        call ps_check(p,1)
      endif


      scheme='dred'
c--- calculate the lowest order matrix element and fill the
c--- common block twopij with s_{ij}
      call qqb_z2jetx_new(p,msq,mqq,ppmsqx,msqx_cs)

c--- write out ug Born amplitude when checking
      if (checkvector .or. checkaxial) then
        write(6,*) 'Madloop check: ug Born',msq(2,0)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

      prop=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth)

      v2(1)=l1
      v2(2)=r1

      do j=1,nf
        vQ(j,1)=L(j)
        vQ(j,2)=R(j)
      enddo

c---- DEBUG: check vector piece
      if (checkvector) then
        prop=czip
      endif

************************************************************************
*     Endpoint contributions from QQGG matrix elements                 *
*     Loop matrix elements are also initialized here                   *
************************************************************************
      if (Gflag) then
c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      if     (colourchoice == 1) then
        subuv(1)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
        subuv(2)=subuv(1)
      elseif (colourchoice == 2) then
        subuv(0)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
      elseif (colourchoice == 3) then
c--- all zero already
      elseif (colourchoice == 0) then
        subuv(1)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
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

c--- Now calculate the relevant lowest-order matrix elements
c--- for each possible initial state from the QQGG contribution

c---  calculate the qqb terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p2)+g(p6)+g(p5)+qb(p1)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(2,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(5,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq_qqb,mmsq_qqb_vec,mmsq_qqb_ax)

c--- obtain qbq from qqb by symmetry
      do polq=1,2
      do polz=1,2
        mmsq_qbq(polq,polz)=mmsq_qqb(3-polq,polz)
        mmsq_qbq_vec(polq,polz)=mmsq_qqb_vec(3-polq,polz)
        mmsq_qbq_ax(polq,polz)=-mmsq_qqb_ax(3-polq,polz)
      enddo
      enddo
c---  calculate the qbq terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p1)+g(p6)+g(p5)+qb(p2)+lbar(p4)+l(p3)
c      do nu=1,4
c      pswap(1,nu)=p(1,nu)
c      pswap(2,nu)=p(6,nu)
c      pswap(3,nu)=p(5,nu)
c      pswap(4,nu)=p(2,nu)
c      pswap(5,nu)=p(4,nu)
c      pswap(6,nu)=p(3,nu)
c      enddo
c      call spinoru(6,pswap,za,zb)
c      call xzqqgg_v(mmsq_qbq,mmsq_qbq_vec,mmsq_qbq_ax)
c      do polq=1,2
c      do polz=1,2
c        write(6,*) mmsq_qbq(polq,polz)-mmsq_qqb(3-polq,polz)
c        write(6,*) mmsq_qbq_vec(polq,polz)-mmsq_qqb_vec(3-polq,polz)
c        write(6,*) mmsq_qbq_ax(polq,polz)+mmsq_qqb_ax(3-polq,polz)
c      enddo
c      enddo
c      pause

c---  calculate the gq terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p5)+g(p6)+g(p1)+qb(p2)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(1,nu)
      pswap(4,nu)=p(2,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq_gq,mmsq_gq_vec,mmsq_gq_ax)

c--- obtain gqb from gq by symmetry
      do polq=1,2
      do polz=1,2
        mmsq_gqb(polq,polz)=mmsq_gq(3-polq,polz)
        mmsq_gqb_vec(polq,polz)=mmsq_gq_vec(3-polq,polz)
        mmsq_gqb_ax(polq,polz)=-mmsq_gq_ax(3-polq,polz)
      enddo
      enddo
c---  calculate the gqb terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p2)+g(p6)+g(p1)+qb(p5)+lbar(p4)+l(p3)
c      do nu=1,4
c      pswap(1,nu)=p(2,nu)
c      pswap(2,nu)=p(6,nu)
c      pswap(3,nu)=p(1,nu)
c      pswap(4,nu)=p(5,nu)
c      pswap(5,nu)=p(4,nu)
c      pswap(6,nu)=p(3,nu)
c      enddo
c      call spinoru(6,pswap,za,zb)
c      call xzqqgg_v(mmsq_gqb,mmsq_gqb_vec,mmsq_gqb_ax)

c---  calculate the qg terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p5)+g(p6)+g(p2)+qb(p1)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(5,nu)
      pswap(2,nu)=p(6,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(1,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq_qg,mmsq_qg_vec,mmsq_qg_ax)

c--- obtain qbg from qg by symmetry
      do polq=1,2
      do polz=1,2
        mmsq_qbg(polq,polz)=mmsq_qg(3-polq,polz)
        mmsq_qbg_vec(polq,polz)=mmsq_qg_vec(3-polq,polz)
        mmsq_qbg_ax(polq,polz)=-mmsq_qg_ax(3-polq,polz)
      enddo
      enddo
c---  calculate the qbg terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p1)+g(p6)+g(p2)+qb(p5)+lbar(p4)+l(p3)
c      do nu=1,4
c      pswap(1,nu)=p(1,nu)
c      pswap(2,nu)=p(6,nu)
c      pswap(3,nu)=p(2,nu)
c      pswap(4,nu)=p(5,nu)
c      pswap(5,nu)=p(4,nu)
c      pswap(6,nu)=p(3,nu)
c      enddo
c      call spinoru(6,pswap,za,zb)
c      call xzqqgg_v(mmsq_qbg,mmsq_qbg_vec,mmsq_qbg_ax)

c--- calculate the gg terms
c-BDKW  0--> q(p1)+g(p2)+g(p3)+qb(p4)+lbar(p5)+l(p6)
CALL    0--> q(p6)+g(p1)+g(p2)+qb(p5)+lbar(p4)+l(p3)
      do nu=1,4
      pswap(1,nu)=p(6,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq_gg,mmsq_gg_vec,mmsq_gg_ax)
      endif

************************************************************************
*     Endpoint contributions from QQQQ matrix elements                 *
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

      endif

      if (discardvecax) then
        mmsq_qqb_ax(:,:)=zip
        mmsq_qbq_ax(:,:)=zip
        mmsq_gq_ax(:,:)=zip
        mmsq_qg_ax(:,:)=zip
        mmsq_qbg_ax(:,:)=zip
        mmsq_gqb_ax(:,:)=zip
        mmsq_gg_ax(:,:)=zip
        mmsq_qqb_vec(:,:)=zip
        mmsq_qbq_vec(:,:)=zip
        mmsq_gq_vec(:,:)=zip
        mmsq_qg_vec(:,:)=zip
        mmsq_qbg_vec(:,:)=zip
        mmsq_gqb_vec(:,:)=zip
        mmsq_gg_vec(:,:)=zip
      endif

c--- Add VIRTUAL terms

************************************************************************
*     Include loop contributions from QQGG matrix elements             *
************************************************************************
      if (Gflag) then

c--- compute correct vector-like coupling for diagrams with Z coupled to a loop
      vcouple(1)=czip
      vcouple(2)=czip
      do j=1,nf
      do polz=1,2
      vcouple(polz)=vcouple(polz)
     & +Q(j)*q1+0.5_dp*(vQ(j,1)+vQ(j,2))*v2(polz)*prop
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
      do polq=1,2
      do polz=1,2

c--- quark-antiquark
      if ((j > 0) .and. (k < 0)) then
        if (j == -k)
     &  msqv(j,k)=msqv(j,k)+half*(aveqq/avegg)*(mmsq_qqb(polq,polz)*(
     &          abs(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_qqb_vec(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_qqb_ax(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

c--- antiquark-quark
      elseif ((j < 0) .and. (k > 0)) then
        if (j == -k)
     &  msqv(j,k)=msqv(j,k)+half*(aveqq/avegg)*(mmsq_qbq(polq,polz)*(
     &          abs(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_qbq_vec(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_qbq_ax(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))


c--- quark-gluon
      elseif ((j > 0) .and. (k == 0)) then

************************** BEGIN MADLOOP CHECKING CODE ************************
        if (checkvector) then
c--- MadLoop check: vector couplings only - no Z coupling to leptons and quarks
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(
     &                   +real(mmsq_qg_vec(polq,polz)
     &          *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*vcouple(polz))
     &                   )
        endif

        if (checkaxial) then
c--- MadLoop check: axial couplings only - photon contribution removed
c------- full Z in Born (no photon)
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(
     &              +real(mmsq_qg_ax(polq,polz)
     &        *conjg(Q(j)*q1*0+vQ(j,polq)*v2(polz)*prop)
     &         *(v2(polz)*prop)/sin2w))

c------- axial Z everywehere (Born and virt, quarks and leptons)
c        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(
c     &              +real(mmsq_qg_ax(polq,polz)
c     &         *conjg(Q(j)*q1*0+
c     &           (vQ(j,polq)-vQ(j,3-polq))*(v2(polz)-v2(3-polz))*prop)
c     &         *((v2(polz)-v2(3-polz))*prop)/sin2w))/8._dp
      endif
*************************** END MADLOOP CHECKING CODE *************************

      if ((checkvector.eqv..false.).and.(checkaxial.eqv..false.)) then
c--- normal case
      msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_qg(polq,polz)*(
     &        abs(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_qg_vec(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_qg_ax(polq,polz)
     &        *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))


        endif

c--- antiquark-gluon
      elseif ((j < 0) .and. (k == 0)) then

        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_qbg(polq,polz)*(
     &          abs(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_qbg_vec(polq,polz)
     &        *conjg(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_qbg_ax(polq,polz)
     &        *conjg(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

c--- gluon-quark
      elseif ((j == 0) .and. (k > 0)) then
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_gq(polq,polz)*(
     &          abs(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_gq_vec(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_gq_ax(polq,polz)
     &        *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

c--- gluon-antiquark
      elseif ((j == 0) .and. (k < 0)) then
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_gqb(polq,polz)*(
     &          abs(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_gqb_vec(polq,polz)
     &        *conjg(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_gqb_ax(polq,polz)
     &        *conjg(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

c--- gluon-gluon
      elseif ((j == 0) .and. (k == 0)) then
        msqv(j,k)=msqv(j,k)+real(ndo,dp)*(mmsq_gg(polq,polz)*(
     &             abs(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_gg_vec(polq,polz)
     &        *conjg(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_gg_ax(polq,polz)
     &        *conjg(Q(1)*q1+vQ(1,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))
        msqv(j,k)=msqv(j,k)+real(nup,dp)*(mmsq_gg(polq,polz)*(
     &          abs(Q(2)*q1+vQ(2,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_gg_vec(polq,polz)
     &        *conjg(Q(2)*q1+vQ(2,polq)*v2(polz)*prop)*vcouple(polz))
     &                 +real(mmsq_gg_ax(polq,polz)
     &        *conjg(Q(2)*q1+vQ(2,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))
      endif
      enddo
      enddo
      enddo
      enddo

      endif

************************************************************************
*     Include loop contributions from QQQQ matrix elements             *
************************************************************************
      if (Qflag) then

      call spinoru(6,p,za,zb)

      faclo=4._dp*V*aveqq*esq**2*gsq**2
      fac=faclo*xn*0.5_dp*ason2pi

c--- Set-up the desired amplitudes, in order to minimize number of calls
      do polq=1,2
      do polz=1,2
      do polb=1,2
c--- atreez
      atreez_526143(polq,polb,polz)=
     &       atreez(polq,polb,polz,5,2,6,1,4,3,za,zb)
      atreez_251643(polq,polb,polz)=
     &       atreez(polq,polb,polz,2,5,1,6,4,3,za,zb)
      atreez_625143(polq,polb,polz)=
     &       atreez(polq,polb,polz,6,2,5,1,4,3,za,zb)
      atreez_261543(polq,polb,polz)=
     &       atreez(polq,polb,polz,2,6,1,5,4,3,za,zb)
      atreez_162543(polq,polb,polz)=
     &       atreez(polq,polb,polz,1,6,2,5,4,3,za,zb)
      atreez_615243(polq,polb,polz)=
     &       atreez(polq,polb,polz,6,1,5,2,4,3,za,zb)
      atreez_152643(polq,polb,polz)=
     &       atreez(polq,polb,polz,1,5,2,6,4,3,za,zb)
      atreez_516243(polq,polb,polz)=
     &       atreez(polq,polb,polz,5,1,6,2,4,3,za,zb)
      atreez_562143(polq,polb,polz)=
     &       atreez(polq,polb,polz,5,6,2,1,4,3,za,zb)
      atreez_651243(polq,polb,polz)=
     &       atreez(polq,polb,polz,6,5,1,2,4,3,za,zb)
      atreez_265143(polq,polb,polz)=
     &       atreez(polq,polb,polz,2,6,5,1,4,3,za,zb)
      atreez_621543(polq,polb,polz)=
     &       atreez(polq,polb,polz,6,2,1,5,4,3,za,zb)
      atreez_561243(polq,polb,polz)=
     &       atreez(polq,polb,polz,5,6,1,2,4,3,za,zb)
      atreez_652143(polq,polb,polz)=
     &       atreez(polq,polb,polz,6,5,2,1,4,3,za,zb)
      atreez_165243(polq,polb,polz)=
     &       atreez(polq,polb,polz,1,6,5,2,4,3,za,zb)
      atreez_612543(polq,polb,polz)=
     &       atreez(polq,polb,polz,6,1,2,5,4,3,za,zb)
c--- a61z
      a61z_526143(polq,polb,polz)=a61z(polq,polb,polz,5,2,6,1,4,3,za,zb)
      a61z_251643(polq,polb,polz)=a61z(polq,polb,polz,2,5,1,6,4,3,za,zb)
      a61z_625143(polq,polb,polz)=a61z(polq,polb,polz,6,2,5,1,4,3,za,zb)
      a61z_261543(polq,polb,polz)=a61z(polq,polb,polz,2,6,1,5,4,3,za,zb)
      a61z_162543(polq,polb,polz)=a61z(polq,polb,polz,1,6,2,5,4,3,za,zb)
      a61z_615243(polq,polb,polz)=a61z(polq,polb,polz,6,1,5,2,4,3,za,zb)
      a61z_152643(polq,polb,polz)=a61z(polq,polb,polz,1,5,2,6,4,3,za,zb)
      a61z_516243(polq,polb,polz)=a61z(polq,polb,polz,5,1,6,2,4,3,za,zb)
      a61z_562143(polq,polb,polz)=a61z(polq,polb,polz,5,6,2,1,4,3,za,zb)
      a61z_651243(polq,polb,polz)=a61z(polq,polb,polz,6,5,1,2,4,3,za,zb)
      a61z_265143(polq,polb,polz)=a61z(polq,polb,polz,2,6,5,1,4,3,za,zb)
      a61z_621543(polq,polb,polz)=a61z(polq,polb,polz,6,2,1,5,4,3,za,zb)
      a61z_561243(polq,polb,polz)=a61z(polq,polb,polz,5,6,1,2,4,3,za,zb)
      a61z_652143(polq,polb,polz)=a61z(polq,polb,polz,6,5,2,1,4,3,za,zb)
      a61z_165243(polq,polb,polz)=a61z(polq,polb,polz,1,6,5,2,4,3,za,zb)
      a61z_612543(polq,polb,polz)=a61z(polq,polb,polz,6,1,2,5,4,3,za,zb)
c--- a62z
      a62z_625143(polq,polb,polz)=a62z(polq,polb,polz,6,2,5,1,4,3,za,zb)
      a62z_261543(polq,polb,polz)=a62z(polq,polb,polz,2,6,1,5,4,3,za,zb)
      a62z_526143(polq,polb,polz)=a62z(polq,polb,polz,5,2,6,1,4,3,za,zb)
      a62z_251643(polq,polb,polz)=a62z(polq,polb,polz,2,5,1,6,4,3,za,zb)
      a62z_152643(polq,polb,polz)=a62z(polq,polb,polz,1,5,2,6,4,3,za,zb)
      a62z_516243(polq,polb,polz)=a62z(polq,polb,polz,5,1,6,2,4,3,za,zb)
      a62z_162543(polq,polb,polz)=a62z(polq,polb,polz,1,6,2,5,4,3,za,zb)
      a62z_615243(polq,polb,polz)=a62z(polq,polb,polz,6,1,5,2,4,3,za,zb)
      a62z_265143(polq,polb,polz)=a62z(polq,polb,polz,2,6,5,1,4,3,za,zb)
      a62z_621543(polq,polb,polz)=a62z(polq,polb,polz,6,2,1,5,4,3,za,zb)
      a62z_562143(polq,polb,polz)=a62z(polq,polb,polz,5,6,2,1,4,3,za,zb)
      a62z_651243(polq,polb,polz)=a62z(polq,polb,polz,6,5,1,2,4,3,za,zb)
      a62z_165243(polq,polb,polz)=a62z(polq,polb,polz,1,6,5,2,4,3,za,zb)
      a62z_612543(polq,polb,polz)=a62z(polq,polb,polz,6,1,2,5,4,3,za,zb)
      a62z_561243(polq,polb,polz)=a62z(polq,polb,polz,5,6,1,2,4,3,za,zb)
      a62z_652143(polq,polb,polz)=a62z(polq,polb,polz,6,5,2,1,4,3,za,zb)
c--- a63z
      if (discardvecax) then
      a63z_526143(polq,polb,polz)=czip
      a63z_625143(polq,polb,polz)=czip
      a63z_162543(polq,polb,polz)=czip
      a63z_152643(polq,polb,polz)=czip
      a63z_562143(polq,polb,polz)=czip
      a63z_265143(polq,polb,polz)=czip
      a63z_561243(polq,polb,polz)=czip
      a63z_165243(polq,polb,polz)=czip
      else
      a63z_526143(polq,polb,polz)=a63z(polq,polb,polz,5,2,6,1,4,3,za,zb)
      a63z_625143(polq,polb,polz)=a63z(polq,polb,polz,6,2,5,1,4,3,za,zb)
      a63z_162543(polq,polb,polz)=a63z(polq,polb,polz,1,6,2,5,4,3,za,zb)
      a63z_152643(polq,polb,polz)=a63z(polq,polb,polz,1,5,2,6,4,3,za,zb)
      a63z_562143(polq,polb,polz)=a63z(polq,polb,polz,5,6,2,1,4,3,za,zb)
      a63z_265143(polq,polb,polz)=a63z(polq,polb,polz,2,6,5,1,4,3,za,zb)
      a63z_561243(polq,polb,polz)=a63z(polq,polb,polz,5,6,1,2,4,3,za,zb)
      a63z_165243(polq,polb,polz)=a63z(polq,polb,polz,1,6,5,2,4,3,za,zb)
      endif
      enddo
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf

c---Desired formula =(Att*(A61+A61o+(A63-A63s/n+A62s+A62os)/n)
c                     +Atts*(A61s+A61os+(A63s-A63/n +A62+A62o)/n))
c                   =tamp*(lamp+lampx)
c                   +tamps*(lamps+lampsx)
c   where lamp,lampsx are the pieces without the "s" (5<->6 swap)
c     and lamps,lampx are the pieces that have the "s"
      do polq=1,2
      do polz=1,2
      do polb=1,2
        tamp=czip
        lamp=czip
        lampx=czip
        tamps=czip
        lamps=czip
        lampsx=czip
c--- Q-Q
        if ((j > 0) .and. (k > 0)) then
          idfac=1._dp
          tamp=atreez_526143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -atreez_251643(3-polb,3-polq,polz)
     &         *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
          lamp=a61z_526143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a61z_251643(3-polb,3-polq,polz)
     &         *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
     &        +a63z_526143(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop
          if (j == k) then
c--- statistical factor
          idfac=0.5_dp
          tamps=-(atreez_625143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &         -atreez_261543(3-polb,3-polq,polz)
     &         *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop))
          lampx=-(
     &        -a63z_625143(polq,polb,polz)/xn**2
     &         *v2(polz)/sin2w*prop
     &        +a62z_625143(polq,polb,polz)/xn
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a62z_261543(3-polb,3-polq,polz)/xn
     &         *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop))
          lamps=-(a61z_625143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a61z_261543(3-polb,3-polq,polz)
     &         *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
     &        +a63z_625143(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop)
          lampsx=
     &        -a63z_526143(polq,polb,polz)/xn**2
     &         *v2(polz)/sin2w*prop
     &        +a62z_526143(polq,polb,polz)/xn
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a62z_251643(3-polb,3-polq,polz)/xn
     &         *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
          endif
          if (compare) then
            msqv(j,k)=msqv(j,k)
     &      +idfac*faclo*real(tamp*conjg(tamp))
     &      +idfac*faclo*real(tamps*conjg(tamps))
            if (polq == polb) msqv(j,k)=msqv(j,k)
     &      +idfac*faclo*real(tamp*conjg(tamps))*(-2._dp/xn)
          else
            msqv(j,k)=msqv(j,k)
     &      +idfac*fac*2._dp*real(tamp*conjg(lamp))
     &      +idfac*fac*2._dp*real(tamps*conjg(lamps))
            if (polq == polb) msqv(j,k)=msqv(j,k)
     &      +idfac*fac*2._dp*real(tamp*conjg(lampx))
     &      +idfac*fac*2._dp*real(tamps*conjg(lampsx))
          endif
c--- Qbar-Qbar
        elseif ((j < 0) .and. (k < 0)) then
          idfac=1._dp
          tamp=atreez_162543(polq,polb,polz)
     &         *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &        -atreez_615243(3-polb,3-polq,polz)
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
          lamp=a61z_162543(polq,polb,polz)
     &         *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &        -a61z_615243(3-polb,3-polq,polz)
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     &        +a63z_162543(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop
          if (j == k) then
c--- statistical factor
          idfac=0.5_dp
          tamps=-(atreez_152643(polq,polb,polz)
     &         *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &        -atreez_516243(3-polb,3-polq,polz)
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop))
          lampx=-(
     &        -a63z_152643(polq,polb,polz)/xn**2
     &         *v2(polz)/sin2w*prop
     &         +a62z_152643(polq,polb,polz)/xn
     &         *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &        -a62z_516243(3-polb,3-polq,polz)/xn
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop))
          lamps=-(a61z_152643(polq,polb,polz)
     &         *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &        -a61z_516243(3-polb,3-polq,polz)
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     &        +a63z_152643(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop)
          lampsx=
     &         -a63z_162543(polq,polb,polz)/xn**2
     &         *v2(polz)/sin2w*prop
     &         +a62z_162543(polq,polb,polz)/xn
     &         *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &         -a62z_615243(3-polb,3-polq,polz)/xn
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
          endif
          if (compare) then
            msqv(j,k)=msqv(j,k)
     &      +idfac*faclo*real(tamp*conjg(tamp))
     &      +idfac*faclo*real(tamps*conjg(tamps))
            if (polq == polb) msqv(j,k)=msqv(j,k)
     &      +idfac*faclo*real(tamp*conjg(tamps))*(-2._dp/xn)
          else
            msqv(j,k)=msqv(j,k)
     &      +idfac*fac*2._dp*real(tamp*conjg(lamp))
     &      +idfac*fac*2._dp*real(tamps*conjg(lamps))
            if (polq == polb) msqv(j,k)=msqv(j,k)
     &      +idfac*fac*2._dp*real(tamp*conjg(lampx))
     &      +idfac*fac*2._dp*real(tamps*conjg(lampsx))
          endif
c--- Q-Qbar
        elseif ((j > 0) .and. (k < 0)) then
          tamp=atreez_562143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -atreez_651243(3-polb,3-polq,polz)
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
          lamp=a61z_562143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a61z_651243(3-polb,3-polq,polz)
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     &        +a63z_562143(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop
          if (j == -k) then
          tamps=-(atreez_265143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -atreez_621543(3-polb,3-polq,polz)
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop))
          lampx=-(
     &        -a63z_265143(polq,polb,polz)/xn**2
     &         *v2(polz)/sin2w*prop
     &        +a62z_265143(polq,polb,polz)/xn
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a62z_621543(3-polb,3-polq,polz)/xn
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop))
          lamps=-(a61z_265143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a61z_621543(3-polb,3-polq,polz)
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     &        +a63z_265143(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop)
          lampsx=
     &        -a63z_562143(polq,polb,polz)/xn**2
     &         *v2(polz)/sin2w*prop
     &        +a62z_562143(polq,polb,polz)/xn
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a62z_651243(3-polb,3-polq,polz)/xn
     &         *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
          tampdo=atreez_265143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -atreez_621543(3-polb,3-polq,polz)
     &         *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
          lampdo=a61z_265143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a61z_621543(3-polb,3-polq,polz)
     &         *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
     &        +a63z_265143(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop
          tampup=atreez_265143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -atreez_621543(3-polb,3-polq,polz)
     &         *(Q(2)*q1+vQ(2,polb)*v2(polz)*prop)
          lampup=a61z_265143(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a61z_621543(3-polb,3-polq,polz)
     &         *(Q(2)*q1+vQ(2,polb)*v2(polz)*prop)
     &        +a63z_265143(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop
          endif
          if (compare) then
            msqv(j,k)=msqv(j,k)
     &      +faclo*real(tamp*conjg(tamp))
     &      +faclo*real(tamps*conjg(tamps))
            if (polq == polb) msqv(j,k)=msqv(j,k)
     &      +faclo*real(tamp*conjg(tamps))*(-2._dp/xn)
          else
            msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &        +real(tamp*conjg(lamp))
     &        +real(tamps*conjg(lamps)))
            if (polq == polb) msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &        +real(tamp*conjg(lampx))
     &        +real(tamps*conjg(lampsx)))
          endif
c--- add additional annihilation diagrams if necessary
          if (j == -k) then
            if (jj(j) == 1) then
              if (compare) then
                msqv(j,k)=msqv(j,k)+faclo*(
     &          +real(tampup*conjg(tampup))*real(nup,dp)
     &          +real(tampdo*conjg(tampdo))*real(ndo-1,dp))
              else
                msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &          +real(tampup*conjg(lampup))*real(nup,dp)
     &          +real(tampdo*conjg(lampdo))*real(ndo-1,dp))
              endif
            elseif (jj(j) == 2) then
              if (compare) then
                msqv(j,k)=msqv(j,k)+faclo*(
     &          +real(tampup*conjg(tampup))*real(nup-1,dp)
     &          +real(tampdo*conjg(tampdo))*real(ndo,dp))
              else
                msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &          +real(tampup*conjg(lampup))*real(nup-1,dp)
     &          +real(tampdo*conjg(lampdo))*real(ndo,dp))
              endif
            endif
          endif
c--- Qbar-Q
        elseif ((j < 0) .and. (k > 0)) then
          tamp=atreez_561243(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -atreez_652143(3-polb,3-polq,polz)
     &         *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
          lamp=a61z_561243(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -a61z_652143(3-polb,3-polq,polz)
     &         *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
     &        +a63z_561243(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop
          if (j == -k) then
          tamps=-(atreez_165243(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -atreez_612543(3-polb,3-polq,polz)
     &         *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop))
          lampx=-(
     &        -a63z_165243(polq,polb,polz)/xn**2
     &         *v2(polz)/sin2w*prop
     &        +a62z_165243(polq,polb,polz)/xn
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -a62z_612543(3-polb,3-polq,polz)/xn
     &         *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop))
          lamps=-(a61z_165243(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -a61z_612543(3-polb,3-polq,polz)
     &         *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
     &        +a63z_165243(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop)
          lampsx=
     &        -a63z_561243(polq,polb,polz)/xn**2
     &         *v2(polz)/sin2w*prop
     &        +a62z_561243(polq,polb,polz)/xn
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -a62z_652143(3-polb,3-polq,polz)/xn
     &         *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
          tampdo=atreez_165243(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -atreez_612543(3-polb,3-polq,polz)
     &         *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
          lampdo=a61z_165243(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -a61z_612543(3-polb,3-polq,polz)
     &         *(Q(1)*q1+vQ(1,polb)*v2(polz)*prop)
     &        +a63z_165243(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop
          tampup=atreez_165243(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -atreez_612543(3-polb,3-polq,polz)
     &         *(Q(2)*q1+vQ(2,polb)*v2(polz)*prop)
          lampup=a61z_165243(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -a61z_612543(3-polb,3-polq,polz)
     &         *(Q(2)*q1+vQ(2,polb)*v2(polz)*prop)
     &        +a63z_165243(polq,polb,polz)/xn
     &         *v2(polz)/sin2w*prop
          endif
          if (compare) then
            msqv(j,k)=msqv(j,k)
     &      +faclo*real(tamp*conjg(tamp))
     &      +faclo*real(tamps*conjg(tamps))
            if (polq == polb) msqv(j,k)=msqv(j,k)
     &      +faclo*real(tamp*conjg(tamps))*(-2._dp/xn)
          else
            msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &        +real(tamp*conjg(lamp))
     &        +real(tamps*conjg(lamps)))
            if (polq == polb) msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &        +real(tamp*conjg(lampx))
     &        +real(tamps*conjg(lampsx)))
          endif
c--- add additional annihilation diagrams if necessary
          if (j == -k) then
            if (jj(j) == -1) then
              if (compare) then
                msqv(j,k)=msqv(j,k)+faclo*(
     &          +real(tampup*conjg(tampup))*real(nup,dp)
     &          +real(tampdo*conjg(tampdo))*real(ndo-1,dp))
              else
                msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &          +real(tampup*conjg(lampup))*real(nup,dp)
     &          +real(tampdo*conjg(lampdo))*real(ndo-1,dp))
              endif
            elseif (jj(j) == -2) then
              if (compare) then
                msqv(j,k)=msqv(j,k)+faclo*(
     &          +real(tampup*conjg(tampup))*real(nup-1,dp)
     &          +real(tampdo*conjg(tampdo))*real(ndo,dp))
              else
                msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &          +real(tampup*conjg(lampup))*real(nup-1,dp)
     &          +real(tampdo*conjg(lampdo))*real(ndo,dp))
              endif
            endif
          endif
        endif
      enddo
      enddo
      enddo

      enddo
      enddo

      endif

c--- write out ug virtual amplitude when checking
      if (checkvector .or. checkaxial) then
        write(6,*) 'Madloop check: ug Virt',msqv(2,0)
c      pause
      endif

************************************************************************
*     UV contributions are included here                               *
C     This is the correction to put the answer in UV renormalized      *
C     dred scheme with msbar coupling                                  *
************************************************************************
      do j=-nf,nf
      do k=-nf,nf

      do cs=0,2
      msqv(j,k)=msqv(j,k)-ason2pi*subuv(cs)*msq_cs(cs,j,k)
      enddo

      enddo
      enddo

      return
      end



