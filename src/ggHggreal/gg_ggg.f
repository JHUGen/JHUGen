      DOUBLE PRECISION FUNCTION GG_GGG(IHEL)
C  
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR PROCESS : g g -> g g g (h)
C  
      IMPLICIT NONE
      include 'constants.f'
      integer ngluons
      parameter (ngluons=5)

      INTEGER IHEL(NGLUONS)
      double complex amp_h5g

      integer PERM(NGLUONS),IHELX(NGLUONS)
      DOUBLE COMPLEX Z(6),ZTEMP
      INTEGER I,J

c      LOGICAL CHECKS
c      DATA CHECKS/.false./ 
       
c--   permutations
      INTEGER IP(NGLUONS,6)

c--   color matrix
c     c_ij= (N^2-1)*N^3/4*CIJ/2 in the GM normalization
      INTEGER CIJ(6,6)

      DATA  (IP(I,1),i=1,5)/1,2,3,4,5/
      DATA  (IP(I,2),i=1,5)/1,2,4,3,5/
      DATA  (IP(I,3),i=1,5)/1,3,2,4,5/
      DATA  (IP(I,4),i=1,5)/1,3,4,2,5/
      DATA  (IP(I,5),i=1,5)/1,4,2,3,5/
      DATA  (IP(I,6),i=1,5)/1,4,3,2,5/

      DATA  (CIJ(I,1),i=1,6)/8,4,4,2,2,0/
      DATA  (CIJ(I,2),i=1,6)/4,8,2,0,4,2/
      DATA  (CIJ(I,3),i=1,6)/4,2,8,4,0,2/
      DATA  (CIJ(I,4),i=1,6)/2,0,4,8,2,4/
      DATA  (CIJ(I,5),i=1,6)/2,4,0,2,8,4/
      DATA  (CIJ(I,6),i=1,6)/0,2,2,4,4,8/

      save ip,cij
!$omp threadprivate(ip,cij)

C-----
C  BEGIN CODE
C-----

      GG_GGG = 0D0

         DO I=1,6

         PERM(1)=IP(1,I)
         PERM(2)=IP(2,I)
         PERM(3)=IP(3,I)
         PERM(4)=IP(4,I)
         PERM(5)=IP(5,I)


         call iperm(IHEL,PERM,IHELX,NGLUONS)


         Z(I)=amp_h5g(PERM,IHELX) 

c
c    check on the properties of the amplitudes
c
c         if (checks) then
c         call checker(z,i,perm,ihelx)
c         endif
c=============================================================

         ENDDO !sum over permutations

c
c   sum over color
c         
         DO J = 1, 6
            ZTEMP = (0.D0,0.D0)
            DO I = 1, 6
               ZTEMP = ZTEMP + Z(I)*CIJ(I,J)
            ENDDO
            GG_GGG =GG_GGG+DBLE(ZTEMP*DCONJG(Z(J)))  
         ENDDO   



c     Overall normalization of CIJ
         GG_GGG=GG_GGG*XN**3*(XN**2-1D0)/4d0 
         
      RETURN
      END


