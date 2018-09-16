      subroutine qqb_wbjet(p,msq)
      implicit none
      include 'types.f'

c--- R.K. Ellis,  8/3/04
c--- matrix element squared and averaged over initial colours and spins
c     q(-p1) + b(-p2) --> W^+ + b(p5) + f(p6)
c                         |
c                         --> nu(p3) + e^+(p4)
c---or
c     q(-p1) + b(-p2) --> W^- + b(p5) + f(p6)
c                         |
c                         --> e^-(p3) + nu~(p4)
c--- all momenta are incoming
c--- all momenta are incoming
c--- Extended to include charm quark production via the variable
c--- "flav". Note that NLO corrections are not yet extended this way.

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'heavyflav.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'msq_cs.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     & facqq,prop,Vsm(-nf:nf),qQ_jkki,Qq_jkji,qbQb_jkki,Qbqb_jkji
      complex(dp):: aLLL
      integer:: j,k,i1,i2,i3,i4,i5,i6

C Statement function defines the sum of the following two diagrams

c     q5(L)----<----------q2            q5(L)------<--------q2
c                 0                             0
c                 0                             0
c                 0                             0
c     q6(L)------<--------q1            q6(L)------<--------q1
c             )                                         )
c            (                                         (
c             )                                         )
c     l3(L)-------<-------l4            l3(L)-------<-------l4

      aLLL(i1,i2,i3,i4,i5,i6)=-cone/(s(i3,i4)*s(i2,i5))*
     & (+za(i6,i3)*zb(i2,i1)/(s(i3,i4)+s(i4,i6)+s(i3,i6))
     & *(+zb(i4,i3)*za(i3,i5)+zb(i4,i6)*za(i6,i5))
     & +za(i6,i5)*zb(i4,i1)/(s(i2,i6)+s(i5,i6)+s(i2,i5))
     & *(+zb(i2,i6)*za(i6,i3)+zb(i2,i5)*za(i5,i3)))

c--- initialize matrix elements
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

c--- set up spinors (and dotproducts)
      call spinoru(6,p,za,zb)
      prop=s(3,4)**2/((s(3,4)-wmass**2)**2+wmass**2*wwidth**2)
      facqq=4._dp*V*gsq**2*(gwsq/2._dp)**2*aveqq*prop

c--- now square these amplitudes separating into color structures
c   1) Amplitude
c   2) Amplitude with (5<-->6)
c   0) Interference between above
c

c--- q(i) q(j) --> q(i) ( --> W q(k) ) q(j)
c--- eg u(1)b(2)->nu(3)e^+(4)b(5)d(6)
        qQ_jkki=abs(aLLL(1,2,3,4,5,6))**2+abs(aLLL(1,5,3,4,2,6))**2

c--- q(i) q(j) --> q(i) q(j) ( --> W q(k) ), eg
c--- eg b(1)u(2)->nu(3)e^+(4)b(5)d(6)
c---  exchange (1<-->2) in above
        Qq_jkji=abs(aLLL(2,1,3,4,5,6))**2+abs(aLLL(2,5,3,4,1,6))**2

c--- qb(i) qb(j) --> qb(i) qb(j) ( --> W qb(k) )
c--- eg b~(1)d~(2)->nu(3)e^+(4)b~(5)u~(6)
        Qbqb_jkji=abs(aLLL(6,1,3,4,5,2))**2+abs(aLLL(6,5,3,4,1,2))**2

c--- qb(i) qb(j) --> qb(i) ( --> W qb(k) ) qb(j)
C--- eg d~(1)b~(2)->nu(3)e^+(4)b~(5)u~(6)
C--- exchange (1<-->2) in above
        qbQb_jkki=abs(aLLL(6,2,3,4,5,1))**2+abs(aLLL(6,5,3,4,2,1))**2

c--- set up auxiliary array

      do j=-nf,nf
        Vsm(j)=Vsum(j)
        if (abs(j) >= flav) Vsm(j)=0._dp
c--- make sure that elements are either one or zero
        if (Vsm(j) > 0._dp) Vsm(j)=1._dp
      enddo

      do j=-nf,nf
      do k=-nf,nf

      msq_cs(0,j,k)=0._dp
      msq_cs(1,j,k)=0._dp
      msq_cs(2,j,k)=0._dp

      if ((abs(j) .ne. flav) .and. (abs(k) .ne. flav)) goto 99
      if ((abs(j) == flav) .and. (abs(k) == flav)) goto 99
c--- so that either abs(j) or abs(k) = flav (but not both).

c--- note that, for (q,qb) and (qb,q) contribution (2) refers to the
c--- straightforward case and (1) to the case with 5<->6
c--- This is reversed for (q,q) and (qb,qb)
C--- Note added RKE 8/17/04
C--- Actually this comment doesn't seem to be true for (qb,q)
C--- in order to get the poles to cancel
C--- There is presumably some interplay with the definition
C--- of the various color terms in qqb_wbjet_z.f,

      if ((j > 0) .and. (k < 0)) then
c--- e.g.  b d~ -> b u~
           msq_cs(2,j,k)=facqq*Vsm(k)*Qbqb_jkji
c--- e.g.  u b~ -> b~ d
           msq_cs(1,j,k)=facqq*Vsm(j)*qQ_jkki

      elseif ((j < 0) .and. (k > 0)) then
c--- e.g.  b~ u -> b~ d
           msq_cs(2,j,k)=facqq*Vsm(k)*Qq_jkji
c--- e.g.  d~ b -> b u~
           msq_cs(1,j,k)=facqq*Vsm(j)*qbQb_jkki

      elseif ((j > 0) .and. (k > 0)) then
c--- e.g.  b u -> b d
           msq_cs(1,j,k)=facqq*Vsm(k)*Qq_jkji
c--- e.g.  u b -> b d
           msq_cs(2,j,k)=facqq*Vsm(j)*qQ_jkki

      elseif ((j < 0) .and. (k < 0)) then
c--- e.g.  b~ d~ -> b~ u~
           msq_cs(1,j,k)=facqq*Vsm(k)*Qbqb_jkji
c--- e.g.  d~ b~ -> b~ u~
           msq_cs(2,j,k)=facqq*Vsm(j)*qbQb_jkki
      endif

      msq(j,k)=msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k)

   99 continue
      enddo
      enddo

      return
      end

