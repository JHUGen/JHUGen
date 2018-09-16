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
      implicit none
      include 'types.f'
      complex(dp)::xcli3
      complex(dp)::z,xrli3
      integer::j,Num
      real(dp),parameter::C3fac(0:20)=
     &  (/1.                      e0_dp,
     & -0.375                     e0_dp,
     &  0.078703703703703703704   e0_dp,
     & -0.0086805555555555555556  e0_dp,
     &  0.00012962962962962962963 e0_dp,
     &  0.000081018518518518518519e0_dp,
     & -3.4193571608537594932     e-6_dp,
     & -1.3286564625850340136     e-6_dp,
     &  8.6608717561098513479     e-8_dp,
     &  2.5260875955320399765     e-8_dp,
     & -2.1446944683640647609     e-9_dp,
     & -5.1401106220129789153     e-10_dp,
     &  5.2495821146008294364     e-11_dp,
     &  1.0887754406636318375     e-11_dp,
     & -1.2779396094493695306     e-12_dp,
     & -2.369824177308745209979778 e-13_dp,
     &  3.1043578879654622943      e-14_dp,
     &  5.261758629912506084131839 e-15_dp,
     & -7.538479549949265366       e-16_dp,
     & -1.18623225777522852530825  e-16_dp,
     &  1.8316979965491383382      e-17_dp/)

      IF(NUM-20) 2,2,3
 2    IF(ABS(Z)-1.00000001_dp) 5,6,6
 5    IF(REAL(Z)-0.50000001_dp) 4,7,7
 7    WRITE(*,*) 'FROM FUNCTION xcli3(Num,z): Real(Z)>0.5 NOT FORESEEN'
      GOTO 9
 6    WRITE(*,*) 'FROM FUNCTION xcli3(Num,z): |Z|>1 NOT FORESEEN'
      GOTO 9
 3    WRITE(*,*) 'FROM FUNCTION xcli3(Num,z): Num>20 NOT FORESEEN'
 9    STOP
 4    CONTINUE
C
      xrli3=0

      do 10 j=0,Num
        xrli3 = xrli3 + c3fac(j)*(-Log(1._dp-z))**(j+1)
 10     continue
      xcli3=xrli3
      return
      end


