      subroutine qqb_wbjet_v(p,msqv)
      implicit none
      include 'types.f'

************************************************************************
*     Author: J. M. Campbell                                           *
*     January 22nd, 2004.                                              *
*                                                                      *
*     Calculate the virtual matrix element squared and                 *
*     subtraction terms for the process                                *
*                                                                      *
*     q(-p1) + b(-p2) --> W + b(p5) + f(p6)                            *
*                         |                                            *
*                         --> nu(p3) + e^+(p4)                         *
*                                                                      *
*    This process requires a subset of the matrix elements that are    *
*    needed for W+2 jet production, cf. qqb_w2jet_v.f                  *
*                                                                      *
*    Extended to include charm quark production via the variable "flav"*
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'heavyflav.f'
      include 'ckm.f'
      include 'zprods_com.f'
      include 'scheme.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),fac,
     & p(mxpart,4),q(mxpart,4),Vsm(-nf:nf)
      integer:: nu,j,k
      real(dp):: qqb_jkji,qbq_jkji,qq_jkji,qbqb_jkji,
     &                 qqbxjkik,qbqxjkki,qq_jkki,qbqb_jkki

      scheme='dred'
c--- calculate the lowest order matrix element and fill the
c--- common block twopij with s_{ij}
      call qqb_wbjet(p,msq)

c--- initialize matrix elements
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

c      prop=s(3,4)/(s(3,4)-wmass**2+im*wmass*wwidth)

************************************************************************
*     Contributions from QQQQ matrix elements                          *
************************************************************************

c--- UV counter-term is already included in a6routine.f

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

c--- set-up matrix elements
      call qqbwbj_loop(1,2,3,4,5,6,
     &                       qqb_jkji,qbq_jkji,qq_jkji,qbqb_jkji,
     &                       qqbxjkik,qbqxjkki,qq_jkki,qbqb_jkki)

c--- set up auxiliary array
      do j=-nf,nf
        Vsm(j)=Vsum(j)
        if (abs(j) >= flav) Vsm(j)=0._dp
c--- make sure that elements are either one or zero
        if (Vsm(j) > 0._dp) Vsm(j)=1._dp
      enddo

c--- Add VIRTUAL terms
      do j=-nf,nf
      do k=-nf,nf

      if ((abs(j) .ne. flav) .and. (abs(k) .ne. flav)) goto 99
      if ((abs(j) == flav) .and. (abs(k) == flav)) goto 99
c--- so that either abs(j) or abs(k) = flav (but not both).

      if     ((j > 0) .and. (k < 0)) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &     +Vsm(j)*qqbxjkik
     &     +Vsm(k)*qqb_jkji)
      elseif ((j < 0) .and. (k > 0)) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &     +Vsm(j)*qbqxjkki
     &     +Vsm(k)*qbq_jkji)
      elseif ((j > 0) .and. (k > 0)) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &     +Vsm(j)*qq_jkki
     &     +Vsm(k)*qq_jkji)
      elseif ((j < 0) .and. (k < 0)) then
          msqv(j,k)=msqv(j,k)+fac*aveqq*(
     &     +Vsm(j)*qbqb_jkki
     &     +Vsm(k)*qbqb_jkji)
      endif

   99 continue

c      write(6,*) '_v', j,k,msq(j,k)
c      write(6,*) 'epinv', epinv
c      write(6,*)

      enddo
      enddo

      return
      end


      subroutine qqbwbj_loop(i1,i2,i3,i4,i5,i6,
     &                       qqb_jkji,qbq_jkji,qq_jkji,qbqb_jkji,
     &                       qqbxjkik,qbqxjkki,qq_jkki,qbqb_jkki)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp)::a61_1432(2),a61_2314(2),a61_4132(2),a61_2341(2),
     &               a61_2431(2),a61_2134(2),a61_1342(2),a61_4312(2),
     &               atr_1432(2),atr_2314(2),atr_4132(2),atr_2341(2),
     &               atr_2431(2),atr_2134(2),atr_1342(2),atr_4312(2)
      complex(dp):: atrLLL,atrLRL,a61LLL,a61LRL
      real(dp):: qqb_jkji,qbq_jkji,qq_jkji,qbqb_jkji,
     &                 qqbxjkik,qbqxjkki,qq_jkki,qbqb_jkki
c--- this routine returns various pieces of the tree-loop interference
c--- for the process 0 -> q(1) qb(4) Q(3) Qb(2) + W (-> l(5) + lbar(6))
c--- the returned pieces are labelled according to the qqb -> qqb
c--- crossing and are calculated in the FORM program sq4q, with
c--- reference to the lowest order qqb_w2jetx.f
c--- there are 2 helicity amplitudes, we label LLL as 1, LRL as 2
C--- Note that a61LLL(i1,i2,i3,i4,i5,i6,za,zb) corresponds to all outgoing
*     q(j1,-)+Q(j3,-)+e(j6,-)+q~(j4)+Q~(j2)+e~(j5)
C--- Note that a61LRL(i1,i2,i3,i4,i5,i6,za,zb) corresponds to all outgoing
*     q(j1,-)+Q(j3,+)+e(j6,-)+q~(j4)+Q~(j2)+e~(j5)

c--- q-qb amplitudes
      a61_1432(1)=a61LLL(i1,i4,i3,i2,i5,i6,za,zb)
      a61_1432(2)=a61LRL(i1,i4,i3,i2,i5,i6,za,zb)
      a61_2314(1)=a61LLL(i2,i3,i1,i4,i5,i6,za,zb)
      a61_2314(2)=a61LRL(i2,i3,i1,i4,i5,i6,za,zb)

      atr_1432(1)=conjg(atrLLL(i1,i4,i3,i2,i5,i6,za,zb))
      atr_1432(2)=conjg(atrLRL(i1,i4,i3,i2,i5,i6,za,zb))
      atr_2314(1)=conjg(atrLLL(i2,i3,i1,i4,i5,i6,za,zb))
      atr_2314(2)=conjg(atrLRL(i2,i3,i1,i4,i5,i6,za,zb))

c--- qb-q amplitudes
      a61_4132(1)=a61LLL(i4,i1,i3,i2,i5,i6,za,zb)
      a61_4132(2)=a61LRL(i4,i1,i3,i2,i5,i6,za,zb)
      a61_2341(1)=a61LLL(i2,i3,i4,i1,i5,i6,za,zb)
      a61_2341(2)=a61LRL(i2,i3,i4,i1,i5,i6,za,zb)

      atr_4132(1)=conjg(atrLLL(i4,i1,i3,i2,i5,i6,za,zb))
      atr_4132(2)=conjg(atrLRL(i4,i1,i3,i2,i5,i6,za,zb))
      atr_2341(1)=conjg(atrLLL(i2,i3,i4,i1,i5,i6,za,zb))
      atr_2341(2)=conjg(atrLRL(i2,i3,i4,i1,i5,i6,za,zb))

c--- q-q amplitudes
      a61_2431(1)=a61LLL(i2,i4,i3,i1,i5,i6,za,zb)
      a61_2431(2)=a61LRL(i2,i4,i3,i1,i5,i6,za,zb)
      a61_2134(1)=a61LLL(i2,i1,i3,i4,i5,i6,za,zb)
      a61_2134(2)=a61LRL(i2,i1,i3,i4,i5,i6,za,zb)

      atr_2431(1)=conjg(atrLLL(i2,i4,i3,i1,i5,i6,za,zb))
      atr_2431(2)=conjg(atrLRL(i2,i4,i3,i1,i5,i6,za,zb))
      atr_2134(1)=conjg(atrLLL(i2,i1,i3,i4,i5,i6,za,zb))
      atr_2134(2)=conjg(atrLRL(i2,i1,i3,i4,i5,i6,za,zb))

c--- qb-qb amplitudes
      a61_1342(1)=a61LLL(i1,i3,i4,i2,i5,i6,za,zb)
      a61_1342(2)=a61LRL(i1,i3,i4,i2,i5,i6,za,zb)
      a61_4312(1)=a61LLL(i4,i3,i1,i2,i5,i6,za,zb)
      a61_4312(2)=a61LRL(i4,i3,i1,i2,i5,i6,za,zb)

      atr_1342(1)=conjg(atrLLL(i1,i3,i4,i2,i5,i6,za,zb))
      atr_1342(2)=conjg(atrLRL(i1,i3,i4,i2,i5,i6,za,zb))
      atr_4312(1)=conjg(atrLLL(i4,i3,i1,i2,i5,i6,za,zb))
      atr_4312(2)=conjg(atrLRL(i4,i3,i1,i2,i5,i6,za,zb))

* Note: in the expressions below, there are no LRL (2) interference
*       terms, since these would correspond to polarizations
*       that don't match, ie. (LL -> RR) x (LR -> LR)
C--transition (2-->1) (6->2 in original notation)
      qqb_jkji=real(a61_1432(1)*atr_1432(1)+a61_1432(2)*atr_1432(2))
C--transition (4-->2) (1->6 in original notation)
      qqbxjkik=real(a61_2314(1)*atr_2314(1)+a61_2314(2)*atr_2314(2))

C--transition (2-->4) (6->1 in original notation)
      qbqxjkki=real(a61_4132(1)*atr_4132(1)+a61_4132(2)*atr_4132(2))
C--transition (1-->2) (2->6 in original notation)
      qbq_jkji=real(a61_2341(1)*atr_2341(1)+a61_2341(2)*atr_2341(2))

C--transition (1-->2) (2->6 in original notation)
      qq_jkji=real(a61_2431(1)*atr_2431(1)+a61_2431(2)*atr_2431(2))
C--transition (4-->2) (1->6 in original notation)
      qq_jkki=real(a61_2134(1)*atr_2134(1)+a61_2134(2)*atr_2134(2))

C--transition (2-->1) (6->2 in original notation)
      qbqb_jkji=real(a61_1342(1)*atr_1342(1)+a61_1342(2)*atr_1342(2))
C--transition (2-->4) (6->1 in original notation)
      qbqb_jkki=real(a61_4312(1)*atr_4312(1)+a61_4312(2)*atr_4312(2))

      return
      end




