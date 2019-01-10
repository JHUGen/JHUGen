      DOUBLE PRECISION FUNCTION ALPHAS(Q,AMZ,NLOOP)
c--- this is simply a wrapper to the LHAPDF implementation of the
c--- running coupling, in the style of the native MCFM routine

c--- Note that the inputs AMZ and NLOOP are not used
      IMPLICIT NONE
      DOUBLE PRECISION Q,AMZ,alphasPDF
      INTEGER NLOOP
            
      ALPHAS=alphasPDF(Q)

      RETURN
      END


