      FUNCTION XCDIL(Z)
C
C CALLED BY XSPENZ FOR THE CALCULATION OF THE EULER DILOG
C completely rewritten 13-06-2008
C evaluates Li_2(z) for complex z with |z|<1, Real(z)<1/2
C
      implicit none
      include 'types.f'
      include 'constants.f'
      COMPLEX(dp)::XCDIL
      integer::j,Num
      COMPLEX(dp)::Z,Z1,CLZ,XLI2
      real(dp)::AAZ,REZ

      real(dp),parameter::C2fac(20)=
     & (/ -0.027777777777777777778e0_dp,
     &     0.00027777777777777777778e0_dp,
     &    -4.7241118669690098262 e-6_dp,
     &     9.1857730746619635509 e-8_dp,
     &    -1.8978869988970999072 e-9_dp,
     &     4.0647616451442255268 e-11_dp,
     &    -8.9216910204564525552 e-13_dp,
     &     1.9939295860721075687 e-14_dp,
     &    -4.5189800296199181917 e-16_dp,
     &     1.0356517612181247014 e-17_dp,
     &    -2.3952186210261867457 e-19_dp,
     &     5.5817858743250093363 e-21_dp,
     &    -1.3091507554183212858 e-22_dp,
     &     3.0874198024267402932 e-24_dp,
     &    -7.3159756527022034204 e-26_dp,
     &     1.740845657234000741  e-27_dp,
     &    -4.1576356446138997196 e-29_dp,
     &     9.9621484882846221032 e-31_dp,
     &    -2.3940344248961653005 e-32_dp,
     &     5.7683473553673900843 e-34_dp/)
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
 13	continue
	xcdil=xli2
      RETURN
      END



