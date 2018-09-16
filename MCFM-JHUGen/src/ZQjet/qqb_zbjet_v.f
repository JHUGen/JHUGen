      subroutine qqb_zbjet_v(p,msqv)
      implicit none
      include 'types.f'

************************************************************************
*     Author: J. Campbell                                              *
*     August 2005.                                                     *
*                                                                      *
*     Calculate the virtual matrix element squared and                 *
*     subtraction terms for the process                                *
*                                                                      *
*     q(-p1) + qbar(-p2) --> Z + b(p5) + j(p6)                         *
*                            |                                         *
*                            --> e^-(p3) + e^+(p4)                     *
*                                                                      *
************************************************************************
*                                                                      *
*     Adapted from the generic Z+2 jet routine                         *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'msq_cs.f'
      include 'b0.f'
      include 'nflav.f'
      include 'heavyflav.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     & mmsq_gq(2,2),mmsq_qg(2,2),mmsq_qbg(2,2),mmsq_gqb(2,2),
     & p(mxpart,4),pswap(mxpart,4),fac
      complex(dp):: atreez,a61z,a63z,prop,vcouple(2),tamp,lamp,
     & mmsq_gq_vec(2,2),mmsq_gq_ax(2,2),
     & mmsq_qg_vec(2,2),mmsq_qg_ax(2,2),
     & mmsq_qbg_vec(2,2),mmsq_qbg_ax(2,2),
     & mmsq_gqb_vec(2,2),mmsq_gqb_ax(2,2)
      integer:: nu,j,k,cs,polq,polb,polz
      real(dp):: subuv(0:2)
      real(dp):: faclo,v2(2),vQ(nf,2)
      complex(dp):: atreez_526143(2,2,2),atreez_251643(2,2,2),
     &               atreez_625143(2,2,2),atreez_261543(2,2,2),
     &               atreez_162543(2,2,2),atreez_615243(2,2,2),
     &               atreez_152643(2,2,2),atreez_516243(2,2,2),
     &               atreez_562143(2,2,2),atreez_651243(2,2,2),
     &               atreez_561243(2,2,2),atreez_652143(2,2,2),
     &               a61z_526143(2,2,2),a61z_251643(2,2,2),
     &               a61z_625143(2,2,2),a61z_261543(2,2,2),
     &               a61z_162543(2,2,2),a61z_615243(2,2,2),
     &               a61z_152643(2,2,2),a61z_516243(2,2,2),
     &               a61z_562143(2,2,2),a61z_651243(2,2,2),
     &               a61z_561243(2,2,2),a61z_652143(2,2,2),
     &               a63z_526143(2,2,2),a63z_625143(2,2,2),
     &               a63z_162543(2,2,2),a63z_152643(2,2,2),
     &               a63z_562143(2,2,2),a63z_561243(2,2,2),
     &               a63z_652143(2,2,2),a63z_651243(2,2,2)

      logical:: compare
      parameter (compare=.false.)

      scheme='dred'
c--- calculate the lowest order matrix element and fill the
c--- common block twopij with s_{ij}
      call qqb_zbjet(p,msq)

c--- initialize the matrix element squared
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

c--- compute correct vector-like coupling for diagrams with Z coupled to a loop
      vcouple(1)=czip
      vcouple(2)=czip
      do j=1,nf
      do polz=1,2
      vcouple(polz)=vcouple(polz)
     & +Q(j)*q1+0.5_dp*(vQ(j,1)+vQ(j,2))*v2(polz)*prop
      enddo
      enddo

************************************************************************
*     Endpoint contributions from QQGG matrix elements                 *
*     Loop matrix elements are also initialized here                   *
************************************************************************
c----UV counterterm contains the finite renormalization to arrive
c----at the MS bar scheme.
c      subuv(1)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
c--- This version should be more robust wrt nflav. (8/11/05, JMC)
      subuv(1)=2._dp*xn*(epinv*b0/xn-1._dp/6._dp)
      subuv(2)=subuv(1)
      subuv(0)=subuv(1)

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

************************************************************************
*     Endpoint contributions from QQQQ matrix elements                 *
************************************************************************
c--- UV counter-term is already included in a6routine.f

************************************************************************
*     Include loop contributions from QQGG matrix elements             *
************************************************************************

      do j=-nflav,nflav,nflav
      do k=-nflav,nflav,nflav

      do polq=1,2
      do polz=1,2
      if     ((j ==  +flav) .and. (k == 0)) then
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_qg(polq,polz)*(
     &             abs(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_qg_vec(polq,polz)
     &          *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)*vcouple(polz))
     &                     +real(mmsq_qg_ax(polq,polz)
     &          *conjg(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))

      elseif ((j == -flav) .and. (k == 0)) then
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_qbg(polq,polz)*(
     &             abs(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_qbg_vec(polq,polz)
     &        *conjg(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)*vcouple(polz))
     &                     +real(mmsq_qbg_ax(polq,polz)
     &        *conjg(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &        *(v2(polz)*prop)/sin2w))

      elseif ((j == 0) .and. (k == +flav)) then
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_gq(polq,polz)*(
     &             abs(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_gq_vec(polq,polz)
     &          *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)*vcouple(polz))
     &                     +real(mmsq_gq_ax(polq,polz)
     &          *conjg(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &          *(v2(polz)*prop)/sin2w))
      elseif ((j == 0) .and. (k == -flav)) then
        msqv(j,k)=msqv(j,k)+(aveqg/avegg)*(mmsq_gqb(polq,polz)*(
     &             abs(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)**2)
     &                     +real(mmsq_gqb_vec(polq,polz)
     &        *conjg(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)*vcouple(polz))
     &                     +real(mmsq_gqb_ax(polq,polz)
     &        *conjg(Q(-k)*q1+vQ(-k,polq)*v2(polz)*prop)
     &        *(v2(polz)*prop)/sin2w))
      endif

      enddo
      enddo

************************************************************************
*     UV contributions are included here                               *
************************************************************************
      do cs=0,2
      msqv(j,k)=msqv(j,k)+
     &  ason2pi*(-subuv(cs))*msq_cs(cs,j,k)
      enddo

      enddo
      enddo

************************************************************************
*     Include loop contributions from QQQQ matrix elements             *
************************************************************************

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
      atreez_561243(polq,polb,polz)=
     &       atreez(polq,polb,polz,5,6,1,2,4,3,za,zb)
      atreez_652143(polq,polb,polz)=
     &       atreez(polq,polb,polz,6,5,2,1,4,3,za,zb)
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
      a61z_561243(polq,polb,polz)=a61z(polq,polb,polz,5,6,1,2,4,3,za,zb)
      a61z_652143(polq,polb,polz)=a61z(polq,polb,polz,6,5,2,1,4,3,za,zb)
c--- a63z
      a63z_526143(polq,polb,polz)=a63z(polq,polb,polz,5,2,6,1,4,3,za,zb)
      a63z_625143(polq,polb,polz)=a63z(polq,polb,polz,6,2,5,1,4,3,za,zb)
      a63z_162543(polq,polb,polz)=a63z(polq,polb,polz,1,6,2,5,4,3,za,zb)
      a63z_152643(polq,polb,polz)=a63z(polq,polb,polz,1,5,2,6,4,3,za,zb)
      a63z_562143(polq,polb,polz)=a63z(polq,polb,polz,5,6,2,1,4,3,za,zb)
      a63z_561243(polq,polb,polz)=a63z(polq,polb,polz,5,6,1,2,4,3,za,zb)
      a63z_652143(polq,polb,polz)=a63z(polq,polb,polz,6,5,2,1,4,3,za,zb)
      a63z_651243(polq,polb,polz)=a63z(polq,polb,polz,6,5,1,2,4,3,za,zb)
      enddo
      enddo
      enddo

      do j=-nflav,nflav
      do k=-nflav,nflav

      if ((abs(j) .ne. flav) .and. (abs(k) .ne. flav)) goto 99
      if ((abs(j) == flav) .and. (abs(k) == flav)) goto 99
c--- so that either abs(j) or abs(k) = flav (but not both).

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
c--- Q-Q
        if      ((j > 0) .and. (k > 0)) then
          if     (j == +flav) then
            tamp=atreez_526143(polq,polb,polz)
     &           *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          -atreez_251643(3-polb,3-polq,polz)
     &           *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
            lamp=a61z_526143(polq,polb,polz)
     &           *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          -a61z_251643(3-polb,3-polq,polz)
     &           *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
     &          +a63z_526143(polq,polb,polz)/xn
     &           *v2(polz)/sin2w*prop
          elseif (k == +flav) then
            tamp=atreez_625143(polq,polb,polz)
     &           *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          -atreez_261543(3-polb,3-polq,polz)
     &           *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
            lamp=a61z_625143(polq,polb,polz)
     &           *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          -a61z_261543(3-polb,3-polq,polz)
     &           *(Q(k)*q1+vQ(k,polb)*v2(polz)*prop)
     &          +a63z_625143(polq,polb,polz)/xn
     &           *v2(polz)/sin2w*prop
          endif
          if (compare) then
            msqv(j,k)=msqv(j,k)
     &      +faclo*real(tamp*conjg(tamp))
          else
            msqv(j,k)=msqv(j,k)
     &      +fac*2._dp*real(tamp*conjg(lamp))
          endif
c--- Qbar-Qbar
        elseif ((j < 0) .and. (k < 0)) then
          if     (j == -flav) then
            tamp=atreez_162543(polq,polb,polz)
     &           *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &          -atreez_615243(3-polb,3-polq,polz)
     &           *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
            lamp=a61z_162543(polq,polb,polz)
     &           *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &          -a61z_615243(3-polb,3-polq,polz)
     &           *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     &          +a63z_162543(polq,polb,polz)/xn
     &           *v2(polz)/sin2w*prop
          elseif (k == -flav) then
            tamp=atreez_152643(polq,polb,polz)
     &           *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &          -atreez_516243(3-polb,3-polq,polz)
     &           *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
            lamp=a61z_152643(polq,polb,polz)
     &           *(Q(-j)*q1+vQ(-j,polq)*v2(polz)*prop)
     &          -a61z_516243(3-polb,3-polq,polz)
     &           *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     &          +a63z_152643(polq,polb,polz)/xn
     &           *v2(polz)/sin2w*prop
          endif
          if (compare) then
            msqv(j,k)=msqv(j,k)
     &      +faclo*real(tamp*conjg(tamp))
          else
            msqv(j,k)=msqv(j,k)
     &      +fac*2._dp*real(tamp*conjg(lamp))
          endif
c--- Q-Qbar
        elseif ((j > 0) .and. (k < 0)) then
          if     (j == +flav) then
            tamp=atreez_562143(polq,polb,polz)
     &           *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          -atreez_651243(3-polb,3-polq,polz)
     &           *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
            lamp=a61z_562143(polq,polb,polz)
     &           *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          -a61z_651243(3-polb,3-polq,polz)
     &           *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     &          +a63z_562143(polq,polb,polz)/xn
     &           *v2(polz)/sin2w*prop
          elseif (k == -flav) then
            tamp=atreez_652143(polq,polb,polz)
     &           *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          -atreez_561243(3-polb,3-polq,polz)
     &           *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
            lamp=a61z_652143(polq,polb,polz)
     &           *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &          -a61z_561243(3-polb,3-polq,polz)
     &           *(Q(-k)*q1+vQ(-k,polb)*v2(polz)*prop)
     &          +a63z_652143(polq,polb,polz)/xn
     &           *v2(polz)/sin2w*prop
          endif
          if (compare) then
            msqv(j,k)=msqv(j,k)
     &      +faclo*real(tamp*conjg(tamp))
          else
            msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &        +real(tamp*conjg(lamp)))
          endif
c--- Qbar-Q
        elseif ((j < 0) .and. (k > 0)) then
          if     (j == -flav) then
            tamp=atreez_651243(polq,polb,polz)
     &           *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &          -atreez_562143(3-polb,3-polq,polz)
     &           *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
            lamp=a61z_651243(polq,polb,polz)
     &           *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &          -a61z_562143(3-polb,3-polq,polz)
     &           *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
     &          +a63z_651243(polq,polb,polz)/xn
     &           *v2(polz)/sin2w*prop
          elseif (k == +flav) then
            tamp=atreez_561243(polq,polb,polz)
     &           *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &          -atreez_652143(3-polb,3-polq,polz)
     &           *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
            lamp=a61z_561243(polq,polb,polz)
     &           *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &          -a61z_652143(3-polb,3-polq,polz)
     &           *(Q(-j)*q1+vQ(-j,polb)*v2(polz)*prop)
     &          +a63z_561243(polq,polb,polz)/xn
     &           *v2(polz)/sin2w*prop
          endif
          if (compare) then
            msqv(j,k)=msqv(j,k)
     &      +faclo*real(tamp*conjg(tamp))
          else
            msqv(j,k)=msqv(j,k)+fac*2._dp*(
     &        +real(tamp*conjg(lamp)))
          endif
        endif
      enddo
      enddo
      enddo

   99 continue

      enddo
      enddo

      return
      end



