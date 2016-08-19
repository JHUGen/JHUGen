      SUBROUTINE h5g(ANS)
C  
C RETURNS AMPLITUDE SQUARED SUMMED OVER COLORS
C AND HELICITIES FOR THE POINT IN PHASE SPACE P(mxpart,4)
C FOR PROCESS : g g -> g g g h 

C overall coupling factors have been removed
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL,   NCOMB              
      PARAMETER (NEXTERNAL=5, NCOMB= 16)

C  
C ARGUMENTS 
C  
      DOUBLE PRECISION ANS
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB)
      DOUBLE PRECISION T,GG_GGG
      INTEGER IHEL
      include 'hels.f'

C ----------
C BEGIN CODE
C ----------


      ANS=0D0

c--  sum over helicities
      DO IHEL=1,16
              T=GG_GGG(NHEL(1,IHEL))            
              ANS=ANS+T
      ENDDO

c--  Multiply by 2d0 to account for the other 16 helicity configs
      ANS=ANS*2d0

      RETURN
      END
       
       

