      SUBROUTINE h5g(ANS)
C
C RETURNS AMPLITUDE SQUARED SUMMED OVER COLORS
C AND HELICITIES FOR THE POINT IN PHASE SPACE P(mxpart,4)
C FOR PROCESS : g g -> g g g h

C overall coupling factors have been removed
C

      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
C
C CONSTANTS
C
      integer,parameter:: NEXTERNAL=5, NCOMB= 16

C
C ARGUMENTS
C
      real(dp)::ANS
C
C LOCAL VARIABLES
C
      integer::NHEL(NEXTERNAL,NCOMB)
      real(dp)::T,GG_GGG
      integer:: IHEL
      include 'hels.f'

C ----------
C BEGIN CODE
C ----------


      ANS=zip

c--  sum over helicities
      DO IHEL=1,16
              T=GG_GGG(NHEL(1,IHEL))
              ANS=ANS+T
      ENDDO

c--  Multiply by two to account for the other 16 helicity configs
      ANS=ANS*two

      RETURN
      END



