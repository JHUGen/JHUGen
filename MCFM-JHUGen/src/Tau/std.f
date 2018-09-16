      subroutine STD
      implicit none
      include 'types.f'
      include 'cplx.h'
      include 'constants.f'
      integer I,J,K
      complex(dp):: SPL(10,10),SMN(10,10),C23(10)
      real(dp):: ROOT(10),PLAB(4,10)
      COMMON/CSTD/SPL,SMN
      COMMON/MOM/PLAB
!$omp threadprivate(/MOM/,/CSTD/)


      DO 5 K=1,10
      SPL(K,K)=czip
      SMN(K,K)=SPL(K,K)
      ROOT(K)=sqrt(PLAB(4,K)-PLAB(1,K))
   5  C23(K)=cplx2(PLAB(2,K),PLAB(3,K))
      DO 10 I=2,10
      DO 10 J=1,I-1
      SPL(I,J)=C23(I)*ROOT(J)/ROOT(I) - C23(J)*ROOT(I)/ROOT(J)
      SPL(J,I)=-SPL(I,J)
      SMN(I,J)=-conjg(SPL(I,J))
      SMN(J,I)=-SMN(I,J)
  10  CONTINUE
      RETURN
      END
C
