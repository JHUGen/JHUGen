*  cpolylog_main.f 
*  By Tord Riemann
*  part of cpolylog.f
*  Latter is available at 
*  http://www-zeuthen.desy.de/theory/research/bhabha/bhabha.html

*  In case of any use, please, refer to the above webpage and to
*  S. Actis, M. Czakon, J. Gluza, T. Riemann
*  "Virtual Hadronic and Heavy-Fermion O(alpha**2) Corrections to Bhabha Scattering"
*  arXiv:0807.4691 [hep-ph], subm. to PRD
*  In the article, App. F, we give the mathematical background and references
*  for the formulae used

*================================================================================*)
*  cpolylog.f v.1.00 (26 Aug 2008)
*================================================================================*)
c 2008-06-09: tr added complex trilog function xtrilog(z) like the complex xspenz(z)
c
c function xspenz(z0)
c   transforms the argument of dilogarithm=li_2 with arbitrary complex argument z0 
c   to the region |z|<1, real(z)<0.5 and calls then the series xcdil(z)
c
c function xcdil(z)
c   a series expansion for complex dilogarithm=li_2 
c
c function xli3(z0)
c   transforms the argument of trilogarithm=li_3 with arbitrary complex argument z0 
c   to the region |z|<1, real(z)<0.5 and calls then the series xcli3(num,z)
c
c function xcli3(num,z)
c   a series expansion for complex trilogarithm=li_3 

c xspenz=li_2 is a copy of some auxiliary function from zfitter 6.30
c 1984 or so

c code tested with opensuse 10.3., gfortran: gcc version 4.2.1 (suse linux)
c test function: testcli3.f

C
      FUNCTION XSPENZ(Z)
C
C CALCULATES THE EULER DILOG FOR ARBITRARY COMPLEX ARGUMENT Z
C CALLS THE SERIES CALCULATION XCDIL 
C T. RIEMANN, ABOUT 1984, FOLLOWING:  
C t Hooft, G. and Veltman, M., SCALAR ONE LOOP INTEGRALS, 
C Nucl. Phys.B153, 365, App. A
C 

      IMPLICIT COMPLEX*16(X,Y)
      IMPLICIT REAL*8(A-H,O-W,Z)
      COMPLEX*16 Z,XCDIL,CSPENZ,LOG
      COMMON/CDZPIF/PI,F1
      EXTERNAL XCDIL
      DATA N/0/
C
      JPRINT=0
C
      N=N+1
      IF(N-1) 71,71,72
 71   PI=2.*(ACOS(0.D0)+ASIN(0.D0))
      F1=PI**2/6.
 72   CONTINUE
      N=1 ! added by JC to prevent overflow of N
      REZ=REAL(Z)
      AIZ=AIMAG(Z)
      AAZ=ABS(Z)
      IF(AAZ) 11,9,11
 9    CSPENZ=DCMPLX(0.D0,0.D0)
      GOTO 99
 11   IF(AAZ-1.) 6,4,1
 1    RE1=REAL(1./Z)
      IF(RE1-.5) 3,3,2
2     CONTINUE
      CSPENZ=XCDIL(1.-1./Z)-2.*F1-LOG(Z)*LOG(1.-1./Z)
     U      -.5*(LOG(-Z))**2
      GOTO 99
3     CONTINUE
      CSPENZ=-XCDIL(1./Z)-F1-.5*LOG(-Z)**2
      GOTO 99
 4    IF(REZ-1.) 6,5,1
 5    CSPENZ=DCMPLX(F1,0.D0)
      GOTO 99
 6    IF(REZ-.5) 7,7,8
7     CONTINUE
      CSPENZ=XCDIL(Z)
      GOTO 99
8     CONTINUE
      CSPENZ=-XCDIL(1.-Z)+F1-LOG(Z)*LOG(1.-Z)
 99   CONTINUE
      AAS= ABS(CSPENZ)
      RES=REAL(CSPENZ)
      AIS=AIMAG(CSPENZ)
      IF(JPRINT) 97,97,98
 98   CONTINUE
 97   CONTINUE
      XSPENZ=CSPENZ
      END

C================================================================ 
      FUNCTION XCDIL(Z)
C
C CALLED BY XSPENZ FOR THE CALCULATION OF THE EULER DILOG
C completely rewritten 13-06-2008
C evaluates Li_2(z) for complex z with |z|<1, Real(z)<1/2
C

      IMPLICIT COMPLEX*16(X,Y)
      IMPLICIT REAL*8(A-H,O-W,Z)
C
      COMPLEX*16 Z,Z1,CLZ,LOG
      COMMON/CDZPIF/PI,F1
c      DIMENSION ZETA(15)

      DATA N/0/,TCH/1.D-20/

      DOUBLE PRECISION C2fac(1:20)

      DATA c2fac(1)  / -0.027777777777777777778   d0 /
      DATA c2fac(2)  /  0.00027777777777777777778 d0 /
      DATA c2fac(3)  / -4.7241118669690098262     d-6/
      DATA c2fac(4)  /  9.1857730746619635509     d-8/
      DATA c2fac(5)  / -1.8978869988970999072 d-9    /
      DATA c2fac(6)  /  4.0647616451442255268 d-11   /
      DATA c2fac(7)  / -8.9216910204564525552 d-13   /
      DATA c2fac(8)  /  1.9939295860721075687 d-14   /
      DATA c2fac(9)  / -4.5189800296199181917 d-16   /
      DATA c2fac(10) /  1.0356517612181247014 d-17   /
      DATA c2fac(11) / -2.3952186210261867457 d-19   /
      DATA c2fac(12) /  5.5817858743250093363 d-21   /
      DATA c2fac(13) / -1.3091507554183212858 d-22   /
      DATA c2fac(14) /  3.0874198024267402932 d-24   /
      DATA c2fac(15) / -7.3159756527022034204 d-26   /
      DATA c2fac(16) /  1.740845657234000741  d-27   /
      DATA c2fac(17) / -4.1576356446138997196 d-29   /
      DATA c2fac(18) /  9.9621484882846221032 d-31   /
      DATA c2fac(19) / -2.3940344248961653005 d-32   /
      DATA c2fac(20) /  5.7683473553673900843 d-34   /
C
C ... Num = number of terms in sum xcdil(z)  for Li_2(z)
C     Num=12 gives more then 12 safe digits
C
      Num=12
C

      AAZ=ABS(Z)
      REZ=REAL(Z)
      IF(AAZ-1.) 4,2,3
 3    PRINT 1000
 1000 FORMAT(3X,6 (15HERROR MODULUS Z) )
 2    IF(REZ-.5) 4,4,3
 4    CONTINUE
      Z1=DCMPLX(1.D0,0.D0)-Z
      CLZ=LOG(Z1)
      xli2=-CLZ-.25*(CLZ)**2
       do 13 j=1,Num
        xli2 = xli2 + c2fac(j) * clz**(2*j+1)
 13      continue
      xcdil=xli2
      RETURN
      END
 


C================================================================ 
      FUNCTION xli3(z)
C
C T. RIEMANN, 06/2008, see also Mathematica notebook Li3.nb
C TRANSFORMS THE ARGUMENT OF trilogarithm WITH ARBITRARY COMPLEX ARGUMENT 
C TO THE REGION |Z|<1, Real(Z)<0.5 AND CALLS THEN A SERIES
C BASED ON SOME FORMULAE OF LAST CENTURIES TAKEN FROM MATHWORLD:
C http://mathworld.wolfram.com/Trilogarithm.html
C  
C THE Real ARGUMENT Z WITH Z>1 IS FORBIDDEN - THERE IS THE CUT
C (AND WE WOULD HAVE LOG OF NEGATIVE REAL ARGUMENT HERE)
C
C CALLS THE SERIES xcli3 
C 

      IMPLICIT COMPLEX*16(X,Y)
      IMPLICIT REAL*8(A-H,O-W,Z)
      COMPLEX*16 Z,xcli3,cli3,LOG
      COMMON/CDZPIF/PI,F1
      EXTERNAL xcli3
C
C ... Num = number of terms in sum xcli3(Num,z) for Li_3(z)
C     Num=12 gives about 12 safe digits
C     Num=20 gives about 18 safe digits
C
      Num=20      
C
      PI=3.1415926535897932385d0
      fzeta2=1.6449340668482264365d0
      F1=fzeta2
      fzeta3=1.2020569031595942854d0
C
      REZ=REAL(Z)
      AIZ=AIMAG(Z)
      AAZ=ABS(Z)
C
      IF(AIZ) 70,71,70       
 71   IF(REZ-1.D0) 70,70,72
 72   WRITE(*,*) 'FROM xli3(z): z is real with z>1, forbidden'
      stop 
 70   CONTINUE       
      IF(AAZ) 11,9,11
 9    cli3=DCMPLX(0.D0,0.D0)
      GOTO 99
 11   IF(AAZ-1.d0) 6,4,1
 1    RE1=REAL(1.d0/Z)
      IF(RE1-.5d0) 3,3,2
2     CONTINUE
C real(1/z)>1/2, |z|>1, this is my case 22
      cli3=-xcli3(Num,1.d0-1.d0/z)- xcli3(Num,1.d0-z) + LOG(1.d0/z)*F1 
     U   - 0.5d0*(LOG(1.d0/z))**2 * LOG(1.d0-1.d0/z)
     U   + (LOG(1.d0/z))**3/6.d0 - LOG(-z)*F1 - (LOG(-z))**3/6.d0 
     U      + fzeta3
      GOTO 99
3     CONTINUE
      cli3=xcli3(Num,1.d0/z) - F1 *LOG(-z) -(LOG(-z))**3/6.d0
      GOTO 99
 4    IF(REZ-1.d0) 6,5,1
 5    cli3=DCMPLX(fzeta3,0.D0)
      GOTO 99
 6    IF(REZ-.5d0) 7,7,8
7     CONTINUE
      cli3=xcli3(Num,z)
      GOTO 99
8     CONTINUE
C |z|<1, Real(z)>0.5, this is my case 12 
      cli3=-xcli3(Num,1.d0-1.d0/z) - xcli3(Num,1.d0-z) + LOG(z)*F1 
     U - 0.5d0*(LOG(z))**2 * LOG(1.d0-z) + (LOG(z))**3/6.d0 
     U      + fzeta3

 99   CONTINUE
      AAS= ABS(cli3)
      RES=REAL(cli3)
      AIS=AIMAG(cli3)
c      IF(JPRINT) 97,97,98
c 98   CONTINUE
c 97   CONTINUE
      xli3=cli3
      END

C================================================================ 

C================================================================ 
      FUNCTION xcli3(Num,z)
C
C T. RIEMANN, 2008 
C BASED ON: J. VOLLINGA, S. WEINZIERL, 
C Comput.Phys.Commun. 167 (2005) 177, hep-ph/0410259CPC, EQNS. 48,49
C CALLED BY xli3 FOR THE CALCULATION OF THE TRILOG WITH ARBITRARY
C COMPLEX ARGUMENT
C ----------------------------------------------------------------
C xcli3 IS A SERIES FOR COMPLEX Li_3(z) BASED ON BERNOULLI-RELATED 
C NUMBERS;      ARGUMENT z: |z|<1 AND Real(z)<0.5
C ----------------------------------------------------------------
C SOME REMARKS:
C Num=13 GAVE IN THE TESTS ALL THE SHOWN DIGITS CORRECT (DOUBLE COMPLEX) 
C THE Real(z)<0.5 IS IMPORTANT FOR THE GOOD CONVERGENCE
C THE |z|<1 IS LESS IMPORTANT, LOOK E.G. AT Z=0.5+10*I
C
C
      IMPLICIT COMPLEX*16(X,Y)
      IMPLICIT REAL*8(A-H,O-W,Z)
      include 'types.f'
      include 'cplx.h'

      COMPLEX*16 Z,LOG
      DOUBLE PRECISION C3fac(0:20)

      DATA c3fac(0)  / 1.                        d0    /
      DATA c3fac(1)  /-0.375                     d0    /
      DATA c3fac(2)  / 0.078703703703703703704   d0    /
      DATA c3fac(3)  /-0.0086805555555555555556  d0    /
      DATA c3fac(4)  / 0.00012962962962962962963 d0    /       
      DATA c3fac(5)  / 0.000081018518518518518519d0    /                          
      DATA c3fac(6)  /-3.4193571608537594932     d-6   /   
      DATA c3fac(7)  /-1.3286564625850340136     d-6   /                              
      DATA c3fac(8)  / 8.6608717561098513479     d-8   /                               
      DATA c3fac(9)  / 2.5260875955320399765     d-8   /                              
      DATA c3fac(10) /-2.1446944683640647609     d-9   /                            
      DATA c3fac(11) /-5.1401106220129789153     d-10  /                           
      DATA c3fac(12) / 5.2495821146008294364     d-11  /                            
      DATA c3fac(13) / 1.0887754406636318375     d-11  /                          
      DATA c3fac(14) /-1.2779396094493695306     d-12  /                               
      DATA c3fac(15) /-2.369824177308745209979778 d-13 /                            
      DATA c3fac(16) / 3.1043578879654622943      d-14 /                           
      DATA c3fac(17) / 5.261758629912506084131839 d-15 /                           
      DATA c3fac(18) /-7.538479549949265366       d-16 /                            
      DATA c3fac(19) /-1.18623225777522852530825  d-16 /                       
      DATA c3fac(20) / 1.8316979965491383382      d-17 /    

      IF(NUM-20) 2,2,3
 2     IF(ABS(Z)-1.00000001D0) 5,6,6
 5     IF(REAL(Z)-0.50000001D0) 4,7,7
 7     WRITE(*,*) 'FROM FUNCTION xcli3(Num,z): Real(Z)>0.5 NOT FORESEEN'
      GOTO 9
 6      WRITE(*,*) 'FROM FUNCTION xcli3(Num,z): |Z|>1 NOT FORESEEN'
      GOTO 9
 3      WRITE(*,*) 'FROM FUNCTION xcli3(Num,z): Num>20 NOT FORESEEN'
 9      STOP
 4      CONTINUE
C
      xrli3=cplx2(0._dp, 0._dp)

       do 10 j=0,Num
        xrli3 = xrli3 + c3fac(j)*(-Log(1.d0-z))**(j+1)
 10      continue
      xcli3=xrli3
      return
      end



c===================================================================
c===================================================================
c===================================================================
c===================================================================




