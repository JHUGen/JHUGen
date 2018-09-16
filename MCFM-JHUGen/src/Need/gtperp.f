      SUBROUTINE GTPERP(PTDQ,P,I,J,K,C,D)
      IMPLICIT NONE
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
C---Find the vectors perpendicular to P(I) and P(J)
C   C AND D are purely space-like vectors in the P(I)+P(J) CMF,
C   with C in the same plane as P(K) and D perpendicular to it,
C   both having length PTDQ*SQRT(2*DOT(P,I,J))

      real(dp):: PTDQ,P(mxpart,4),C(4),D(4),PTF,DIJ,DIK,DJK,DOT
      real(dp):: QI(4),QJ(4),EPS4
      integer:: I,J,K,L
      DIJ=DOT(P,I,J)
      DIK=DOT(P,I,K)
      DJK=DOT(P,J,K)
      PTF=PTDQ/SQRT(DIK*DJK)
      DO L=1,4
        C(L)=PTF*(DIJ*P(K,L)-DJK*P(I,L)-DIK*P(J,L))
        QI(L)=P(I,L)
        QJ(L)=P(J,L)
      ENDDO
      DO L=1,4
        D(L)=EPS4(L,QI,QJ,C)/DIJ
      ENDDO
      END

C-----------------------------------------------------------------------
      function EPS4(I,A,B,C)
      implicit none
      include 'types.f'
      real(dp):: EPS4
      real(dp):: EPS3,A(4),B(4),C(4),AA(3),BB(3),CC(3)
      integer::I,J,K
      integer,PARAMETER::S(4)=(/+1,-1,+1,+1/)
      J=1
      DO K=1,3
        IF (I.EQ.J) J=J+1
        AA(K)=A(J)
        BB(K)=B(J)
        CC(K)=C(J)
        J=J+1
      ENDDO
      EPS4=0._dp
      DO J=1,3
        EPS4=EPS4+CC(J)*EPS3(J,AA,BB)
      ENDDO
      EPS4=S(I)*EPS4
      END
C-----------------------------------------------------------------------
      function EPS3(I,A,B)
      implicit none
      include 'types.f'
      real(dp):: EPS3
      real(dp):: A(3),B(3),AA(2),BB(2)
      integer:: I,J,K
      integer,PARAMETER::S(3)=(/+1,-1,+1/)
      J=1
      DO K=1,2
        IF (I.EQ.J) J=J+1
        AA(K)=A(J)
        BB(K)=B(J)
        J=J+1
      ENDDO
      EPS3=S(I)*(AA(1)*BB(2)-AA(2)*BB(1))
      END
C-----------------------------------------------------------------------

