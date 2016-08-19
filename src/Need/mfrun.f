C-----------------------------------------------------------------------------
C
      double precision function massfrun(mf,scale,asmz,nloop)
C
C-----------------------------------------------------------------------------
C
C       This function returns the 'nloop' value of a MSbar fermion mass
C       at a given scale.
C
C       INPUT: mf    = MSbar mass of fermion at MSbar fermion mass scale 
C              scale = scale at which the running mass is evaluated
C              asmz  = AS(MZ) : this is passed to alphas(scale,asmz,2)
C              nloop = # of loops in the evolutionC       
C
C       COMMON BLOCKS: COMMON/QMASS/CMASS,BMASS,TMASS
C                      contains the MS-bar masses of the heavy quarks.
C
C       EXTERNAL:      double precision alphas(scale,asmz,2)
C                      
C-----------------------------------------------------------------------------
C
      implicit none
      include 'masses.f'
C
C     ARGUMENTS
C
      double precision  mf,scale,asmz
      integer           nloop
C
C     LOCAL
C
      double precision  beta0, beta1,gamma0,gamma1
      double precision  A1,as,asmf,l2
      integer  nf
C
C     EXTERNAL
C
      double precision  alphas
      external          alphas
C
C     COMMON
C
c      real *8      cmass,bmass,tmass
c      COMMON/QMASS/CMASS,BMASS,TMASS
c
c     CONSTANTS
c
      double precision  One, Two, Three, Pi
      parameter( One = 1d0, Two = 2d0, Three = 3d0 )
      parameter( Pi = 3.14159265358979323846d0) 
cc
C
C

      if ( mf.gt.mt ) then
         nf = 6
      else
         nf = 5
      end if

      beta0 = ( 11d0 - Two/Three *nf )/4d0
      beta1 = ( 102d0  - 38d0/Three*nf )/16d0
      gamma0= one
      gamma1= ( 202d0/3d0  - 20d0/9d0*nf )/16d0
      A1    = -beta1*gamma0/beta0**2+gamma1/beta0
      as    = alphas(scale,asmz,nloop)
      asmf  = alphas(mf   ,asmz,nloop)
      l2    = (one+A1*as/Pi)/(one+A1*asmf/Pi)
      
      massfrun = mf * (as/asmf)**(gamma0/beta0)

      if(nloop.eq.2) massfrun=massfrun*l2
ccc
      return
      end


