      subroutine qqb_zbb_v(p,msqv)
      implicit none
      include 'types.f'

************************************************************************
*     Author: J.M. Campbell                                            *
*     March, 1999.                                                     *
*     Updated February, 2000                                           *
*     calculate the virtual matrix element squared and subtraction     *
*     terms for the process                                            *
*     q(-p1)+qb(-p2) --> e^-(p3)+e^+(p4)+b(p5)+bb(p6)                  *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'noglue.f'
      include 'b0.f'
      include 'heavyflav.f'
      include 'nflav.f'
      include 'first.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     & p(mxpart,4),q_bdkw(mxpart,4),faclo,subuv,scalesq,
     & fac,v2(2),vQ(nf,2),mmsq(2,2),pswap(mxpart,4)
      complex(dp):: tamp,lamp,atreez,a61z,prop,
     & atreez_123456(2,2,2),atreez_214356(2,2,2),
     & atreez_423156(2,2,2),atreez_241356(2,2,2),
     & a61z_123456(2,2,2),a61z_214356(2,2,2),
     & a61z_423156(2,2,2),a61z_241356(2,2,2),
     & mmsq_vec(2,2),mmsq_ax(2,2),vcouple(2)
      integer:: nu,j,k,polq,polb,polz
      save scalesq

      scheme='dred'

      if (first) then
       if     (flav == 5) then
         scalesq=mbsq
       elseif (flav == 4) then
         scalesq=mcsq
       else
         write(6,*) 'Invalid flav in qqb_zbb_v.f, flav=',flav
       endif
       first=.false.
      endif

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

c--- shortcut if we're doing gqonly
      if (gqonly) return

c---calculate the lowest order matrix element and fill the common block
c---twopij with s_{ij} (in rke notation)
      call qqb_zbb(p,msq)

      if (
     &      (s(5,6) < four*scalesq)
     & .or. (s(1,5)*s(2,5)/s(1,2) < scalesq)
     & .or. (s(1,6)*s(2,6)/s(1,2) < scalesq) ) return

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

c--- calculate the gg terms
c ---Call the two gluon process which is defined in xzqqgg_v
C ---in the notation
C     0 ---> q(p1)+g(p2)+g(p3)+qbar(p4)+a(p5)+  l(p6)
C ---compared with ours which is:-
c     0 ---> b(p6)+g(p1)+g(p2)+bb(p5)+e^+(p4)+e^-(p3)
      do nu=1,4
      pswap(1,nu)=p(6,nu)
      pswap(2,nu)=p(1,nu)
      pswap(3,nu)=p(2,nu)
      pswap(4,nu)=p(5,nu)
      pswap(5,nu)=p(4,nu)
      pswap(6,nu)=p(3,nu)
      enddo
      call spinoru(6,pswap,za,zb)
      call xzqqgg_v(mmsq,mmsq_vec,mmsq_ax)

c---  Now transform momenta into a notation
c---  suitable for calling the BDKW function with notation which is
c---    q-(-p4)+Q+(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
compared to ours which is (see above)
c--     q (-p1)+b (-p5)+l-(-p4) ---> q+(p2)+b (p6)+e-(p3)

      do nu=1,4
      q_bdkw(4,nu)=p(1,nu)
      q_bdkw(1,nu)=p(2,nu)
      q_bdkw(6,nu)=p(3,nu)
      q_bdkw(5,nu)=p(4,nu)
      q_bdkw(2,nu)=p(5,nu)
      q_bdkw(3,nu)=p(6,nu)
      enddo
      call spinoru(6,q_bdkw,za,zb)

      faclo=4._dp*V*aveqq*esq**2*gsq**2
      fac=faclo*xn*0.5_dp*ason2pi

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

c--- set-up amplitudes first, to improve efficiency
      do polq=1,2
      do polz=1,2
      do polb=1,2
      atreez_123456(polq,polb,polz)
     &      =atreez(polq,polb,polz,1,2,3,4,5,6,za,zb)
      atreez_214356(polq,polb,polz)
     &      =atreez(polq,polb,polz,2,1,4,3,5,6,za,zb)
      atreez_423156(polq,polb,polz)
     &      =atreez(polq,polb,polz,4,2,3,1,5,6,za,zb)
      atreez_241356(polq,polb,polz)
     &      =atreez(polq,polb,polz,2,4,1,3,5,6,za,zb)
      a61z_123456(polq,polb,polz)=a61z(polq,polb,polz,1,2,3,4,5,6,za,zb)
      a61z_214356(polq,polb,polz)=a61z(polq,polb,polz,2,1,4,3,5,6,za,zb)
      a61z_423156(polq,polb,polz)=a61z(polq,polb,polz,4,2,3,1,5,6,za,zb)
      a61z_241356(polq,polb,polz)=a61z(polq,polb,polz,2,4,1,3,5,6,za,zb)
      enddo
      enddo
      enddo

      do j=-nflav,nflav
      k=-j

      do polq=1,2
      do polz=1,2
      do polb=1,2
        if     ((j == 0) .and. (k == 0)) then
          tamp=0._dp
          lamp=0._dp
        elseif ((j > 0) .and. (k < 0)) then
          tamp=atreez_123456(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -atreez_214356(3-polb,3-polq,polz)
     &         *(Q(flav)*q1+vQ(flav,polb)*v2(polz)*prop)
          lamp=a61z_123456(polq,polb,polz)
     &         *(Q(j)*q1+vQ(j,polq)*v2(polz)*prop)
     &        -a61z_214356(3-polb,3-polq,polz)
     &         *(Q(flav)*q1+vQ(flav,polb)*v2(polz)*prop)
        elseif ((j < 0) .and. (k > 0)) then
          tamp=atreez_423156(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -atreez_241356(3-polb,3-polq,polz)
     &         *(Q(flav)*q1+vQ(flav,polb)*v2(polz)*prop)
          lamp=a61z_423156(polq,polb,polz)
     &         *(Q(k)*q1+vQ(k,polq)*v2(polz)*prop)
     &        -a61z_241356(3-polb,3-polq,polz)
     &         *(Q(flav)*q1+vQ(flav,polb)*v2(polz)*prop)
        endif
        msqv(j,k)=msqv(j,k)+fac*2._dp*real(tamp*conjg(lamp))
      enddo
      if ((j == 0) .and. (k == 0)) then
        msqv(j,k)=msqv(j,k)+mmsq(polq,polz)*(
     &    abs(Q(flav)*q1+vQ(flav,polq)*v2(polz)*prop)**2)
     &            +real(mmsq_vec(polq,polz)
     &   *conjg(Q(flav)*q1+vQ(flav,polq)*v2(polz)*prop)*vcouple(polz))
     &            +real(mmsq_ax(polq,polz)
     &   *conjg(Q(flav)*q1+vQ(flav,polq)*v2(polz)*prop)
     &   *(v2(polz)*prop)/sin2w)
      endif
      enddo
      enddo

      enddo

c--- add in UV counter-term for gg sub-process here
c--- (UV subtraction occurs for the other pieces in a6routine.f)
c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      subuv=2._dp*(epinv*b0-xn/6._dp)
      msqv(0,0)=msqv(0,0)-ason2pi*subuv*msq(0,0)

      return
      end

