      subroutine Wbb(q1,q2,q3,q4,q5,q6,q7,m,za,zb,h,qedi,qedf,qcda,qcdb)
      implicit none
      include 'types.f'
c-----Author R.K.Ellis
C     All contributions
C     1^-_q 2^+_qb 3^s3_b 4^s4_bb 5^s5_g 6^-_l 7^+_a
C     calculated in the Rodrigo scheme

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: m,s167,s267,s134,s234,s34,s345,beta,bm,bp,
     & p3Dq5,p4Dq5
      complex(dp):: zba2,helviol,tmp(2,2)
      complex(dp):: qedi(2,2,2),qedf(2,2,2),qcda(2,2,2),qcdb(2,2,2)
      complex(dp):: zba2153,zba2154,zba2157,zba6175,zba6174,zba6173
      complex(dp):: zba3267,zba4267,zba6243,zba5267,zba4351,zba3451

      integer:: q1,q2,q3,q4,q5,q6,q7,h
      zba2(q1,q2,q3,q4)=zb(q1,q2)*za(q2,q4)+zb(q1,q3)*za(q3,q4)
      s34=s(q3,q4)
      s167=s(q1,q6)+s(q1,q7)+s(q6,q7)
      s267=s(q2,q6)+s(q2,q7)+s(q6,q7)
      s134=s(q1,q3)+s(q1,q4)+s(q3,q4)
      s234=s(q2,q3)+s(q2,q4)+s(q3,q4)
      s345=s(q3,q4)+s(q3,q5)+s(q4,q5)
      beta=one-four*m**2/s34
      if (beta > zero) then
         beta=sqrt(beta)
      else
         beta=zero
      endif
      bp=half+half*beta
      bm=half-half*beta
      p3Dq5=+half*bp*s(q3,q5)+half*bm*s(q4,q5)
      p4Dq5=+half*bp*s(q4,q5)+half*bm*s(q3,q5)
      zba2153=zba2(q2,q1,q5,q3)
      zba2154=zba2(q2,q1,q5,q4)
      zba2157=zba2(q2,q1,q5,q7)
      zba6175=zba2(q6,q1,q7,q5)
      zba6174=zba2(q6,q1,q7,q4)
      zba6173=zba2(q6,q1,q7,q3)
      zba3267=zba2(q3,q2,q6,q7)
      zba4267=zba2(q4,q2,q6,q7)
      zba6243=zba2(q6,q2,q4,q3)
      zba5267=zba2(q5,q2,q6,q7)
      zba4351=zba2(q4,q3,q5,q1)
      zba3451=zba2(q3,q4,q5,q1)

C--   order of indices is gluon helicity s4,s5
C--This is qedi, for negative helicity gluon
      tmp(1,2)=(
     & -zba2153*zba4267*zb(q6,q2)/s267
     & -zba2157*zba6243*zb(q4,q2)/s234
     & +za(q1,q7)*zba6175*zb(q2,q4)**2*za(q4,q3)*zb(q5,q1)
     & /s167/s234
     & +za(q1,q3)*zba2(q4,q1,q3,q5)*zb(q2,q6)**2*za(q6,q7)*zb(q5,q1)
     & /s267/s134)
     & /(zb(q5,q1)*zb(q5,q2)*s34)

      tmp(2,1)=(
     & -zba2154*zba3267*zb(q6,q2)/s267
     & -zba2157*zba2(q6,q2,q3,q4)*zb(q3,q2)/s234
     & +za(q1,q7)*zba6175*zb(q2,q3)**2*za(q3,q4)*zb(q5,q1)
     & /s167/s234
     & +za(q1,q4)*zba2(q3,q1,q4,q5)*zb(q2,q6)**2*za(q6,q7)*zb(q5,q1)
     & /s267/s134)
     & /(zb(q5,q1)*zb(q5,q2)*s34)


      helviol=(
     & -zba2153*zba3267*zb(q6,q2)/s267-zba2157*zba6243*zb(q3,q2)/s234
     & -za(q1,q7)*zba6175*zb(q2,q4)*za(q4,q3)*zb(q3,q2)
     & *zb(q5,q1)/s167/s234
     & +za(q1,q3)*zba2(q3,q1,q4,q5)*zb(q2,q6)**2*za(q6,q7)*zb(q5,q1)
     & /s267/s134)
     & /(zb(q5,q2)*zb(q5,q1)*s34)

      tmp(1,1)=two*m/zb(q3,q4)*helviol
      tmp(2,2)=two*m/za(q3,q4)*helviol

      if (h == 1) then
      qedi(1,2,h)=+tmp(1,2)
      qedi(2,1,h)=+tmp(2,1)
      qedi(1,1,h)=+tmp(1,1)
      qedi(2,2,h)=+tmp(2,2)
      elseif (h == 2) then
      qedi(1,2,h)=-tmp(1,2)
      qedi(2,1,h)=-tmp(2,1)
      qedi(2,2,h)=-tmp(1,1)
      qedi(1,1,h)=-tmp(2,2)
      endif

C---This is qedf, for negative helicity gluon
      tmp(1,2)=+(bp/p3Dq5-bm/p4Dq5)/(two*zb(q4,q5)*s345)
     & *((zba4351*zba4267*zb(q2,q6)*za(q3,q5))/s267
     & +(zba6175*zb(q5,q4)+zba6173*zb(q3,q4))
     & *zb(q4,q2)*za(q3,q5)*za(q1,q7)/s167)

      tmp(2,1)=-(bp/p4Dq5-bm/p3Dq5)/(two*zb(q3,q5)*s345)
     & *((zba3451*zba3267*za(q4,q5)*zb(q2,q6))/s267
     & -(zb(q3,q5)*zba6175+zb(q3,q4)*zba6174)
     & *za(q4,q5)*za(q1,q7)*zb(q3,q2)/s167)

      tmp(2,2)=
     & -m/(two*zb(q4,q5)*za(q4,q3)*p3Dq5*s267*s345)*(
     & ((-za(q4,q3)*za(q1,q5)-za(q4,q5)*za(q3,q1))*bp*zb(q4,q3)
     & +zb(q4,q5)*za(q4,q5)*za(q1,q5))*zba4267*zb(q2,q6)
     & +za(q3,q5)*zba3267*zb(q4,q3)*zb(q2,q6)*za(q3,q1)*bp)

     & -m/(two*zb(q4,q5)*za(q4,q3)*p3Dq5*s167*s345)*(
     & +(zb(q4,q5)*za(q4,q5)-bp*zb(q4,q3)*za(q4,q3))
     & *zb(q4,q2)*zba6175*za(q1,q7)
     & +(zb(q4,q2)*za(q4,q5)-zb(q3,q2)*za(q3,q5))
     & *zba6173*zb(q4,q3)*za(q1,q7)*bp)

     & -m/(two*zb(q4,q5)*za(q4,q3)*p4Dq5*s267*s345)*(
     & +zb(q4,q3)*zba4267*zb(q2,q6)*za(q4,q1)*za(q3,q5)*bm
     & +(zb(q4,q5)*za(q1,q5)-bm*zb(q4,q3)*za(q3,q1))
     & *zba3267*zb(q2,q6)*za(q3,q5))

     & -m/(two*zb(q4,q5)*za(q4,q3)*p4Dq5*s167*s345)*(
     & +zb(q3,q2)*zba6175*zb(q4,q5)*za(q3,q5)*za(q1,q7)
     & +(zba6173*zb(q3,q2)-zba2(q6,q7,q1,q4)*zb(q4,q2))
     & *zb(q4,q3)*za(q3,q5)*za(q1,q7)*bm)

      tmp(1,1)=
     & -m/(two*zb(q4,q5)*zb(q4,q3)*p4Dq5*s267*s345)*(
     & -zb(q4,q5)*zba4267*zb(q2,q6)*za(q4,q5)*za(q1,q5)
     & +zb(q4,q3)*zba4267*zb(q2,q6)*bm*za(q4,q1)*za(q3,q5)
     & -zb(q4,q3)*zba3267*zb(q2,q6)*za(q3,q1)*za(q3,q5)*bm)

     & -m/(two*zb(q4,q5)*zb(q4,q3)*p4Dq5*s167*s345)*(
     & -zb(q4,q5)*zb(q4,q2)*zba6175*za(q4,q5)*za(q1,q7)
     & -zb(q4,q3)*zb(q4,q2)*zba6174*za(q1,q7)*za(q3,q5)*bm
     & +zb(q4,q3)*zb(q3,q2)*zba6173*za(q3,q5)*za(q1,q7)*bm)

     & -m/(two*zb(q4,q5)*zb(q4,q3)*p3Dq5*s267*s345)*(
     & -zb(q4,q5)*zba3267*zb(q2,q6)*za(q3,q5)*za(q1,q5)
     & -zb(q4,q3)*zba4267*zb(q2,q6)*za(q4,q1)*za(q3,q5)*bp
     & +zb(q4,q3)*zba3267*zb(q2,q6)*za(q3,q1)*za(q3,q5)*bp)

     & -m/(two*zb(q4,q5)*zb(q4,q3)*p3Dq5*s167*s345)*(
     & -zb(q4,q5)*zb(q3,q2)*zba6175*za(q3,q5)*za(q1,q7)
     & +zb(q4,q3)*zb(q4,q2)*zba6174*za(q3,q5)*za(q1,q7)*bp
     & -zb(q4,q3)*zb(q3,q2)*zba6173*za(q3,q5)*za(q1,q7)*bp)

      if (h == 1) then
      qedf(1,2,h)=+tmp(1,2)
      qedf(2,1,h)=+tmp(2,1)
      qedf(1,1,h)=+tmp(1,1)
      qedf(2,2,h)=+tmp(2,2)
      elseif (h == 2) then
      qedf(1,2,h)=-tmp(1,2)
      qedf(2,1,h)=-tmp(2,1)
      qedf(2,2,h)=-tmp(1,1)
      qedf(1,1,h)=-tmp(2,2)
      endif

C--   order of indices is gluon helicity s3,s4
c-- This is qcdb, for negative helicity gluon
      tmp(1,2)=
     & -(-zba4351*za(q3,q5)*bp*zba4267*zb(q2,q6))
     & /(two*zb(q4,q5)*p3Dq5*s267*s345)

     & -((zb(q4,q5)*zba6175+zb(q4,q3)*zba6173)
     & *zb(q4,q2)*za(q3,q5)*za(q1,q7)*bp)
     & /(zb(q4,q5)*two*p3Dq5*s167*s345)

     & -(zb(q4,q2)*zb(q4,q1)*zb(q2,q6)*za(q3,q1)*za(q1,q7))
     & /(zb(q4,q5)*zb(q2,q5)*s134*s34)

     & -(-zba4267*zb(q4,q1)*zb(q2,q6)*za(q3,q1)*za(q1,q5))
     & /(zb(q4,q5)*s134*s267*s34)

     & -(zba4267*zb(q2,q6)*za(q3,q5)*za(q1,q5))
     & /(s267*s345*s34)

     & -(zb(q4,q2)*zba6175*za(q3,q5)*za(q1,q7))
     & /(s167*s345*s34)

     & +(zb(q4,q2)*zb(q2,q6)*za(q3,q1)*za(q3,q7))
     & /(zb(q4,q5)*zb(q2,q5)*za(q4,q3)*s134)

     & -(zba4267*zb(q2,q6)*za(q3,q1)*za(q3,q5))
     & /(zb(q4,q5)*za(q4,q3)*s134*s267)

     & -(zba4267*zb(q2,q6)*za(q3,q1)*za(q3,q5))
     & /(zb(q4,q5)*za(q4,q3)*s267*s345)

     & -(zb(q4,q2)**2*zba6173*za(q1,q7))
     & /(zb(q4,q5)*zb(q2,q5)*s167*s34)

     & -(zb(q4,q2)*zb(q4,q3)*zba6173*za(q3,q5)*za(q1,q7))
     & /(zb(q4,q5)*s167*s345*s34)

      tmp(1,1)=
     & -m/(zb(q4,q5)*two*p3Dq5*s267*s345)*(
     & -zba4267*za(q4,q1)
     & +zba3267*za(q3,q1))*zb(q2,q6)*za(q3,q5)*bp

     & -m/(zb(q4,q3)*two*p3Dq5*s267*s345)*(
     &  -zba3267*zb(q2,q6)*za(q3,q5)*za(q1,q5))

     & -m/(zb(q4,q5)*two*p3Dq5*s167*s345)*(
     & +zb(q4,q2)*zba6174
     & -zb(q3,q2)*zba6173)*za(q3,q5)*za(q1,q7)*bp

     & -m/(zb(q4,q3)*two*p3Dq5*s167*s345)*(
     & -zb(q3,q2)*zba6175*za(q3,q5)*za(q1,q7))

     & -m/(zb(q4,q3)*zb(q4,q5)*zb(q2,q5)*s134*s34)*(
     &  -zb(q4,q2)*zb(q2,q6)*za(q1,q7)
     &  *(zb(q3,q1)*za(q3,q1)-zb(q4,q1)*za(q4,q1)))

     & -m/(zb(q4,q3)*zb(q4,q5)*zb(q2,q5)*s167*s34)*(
     & -zb(q4,q2)*za(q1,q7)
     &  *(zba6173*zb(q3,q2)-zba6174*zb(q4,q2)))

     & -m/(zb(q4,q3)*zb(q4,q5)*s134*s267*s34)*(
     & -zba4267*zb(q2,q6)*za(q1,q5)
     &  *(zb(q4,q1)*za(q4,q1)-zb(q3,q1)*za(q3,q1)))

     & -m/(zb(q4,q3)*zb(q4,q5)*s267*s345*s34)*(
     & -zba4267*zb(q3,q5)*zb(q2,q6)*za(q3,q5)*za(q1,q5))

     & -m/(zb(q4,q3)*zb(q4,q5)*s167*s345*s34)*(
     & -zb(q4,q2)*zb(q3,q5)*zba6175*za(q3,q5)*za(q1,q7))

     & -m/(zb(q4,q3)*s267*s345*s34)*(
     & +zba4267*zb(q2,q6)*za(q4,q5)*za(q1,q5))

     & -m/(zb(q4,q3)*s167*s345*s34)*(
     & +zb(q4,q2)*zba6175*za(q4,q5)*za(q1,q7))

     & -m/(zb(q4,q5)*zb(q2,q5)*s134*s34)*(
     & +zb(q4,q2)*zb(q2,q6)*(za(q4,q1)*za(q3,q7)+za(q4,q7)*za(q3,q1)))

     & -m/(zb(q5,q4)*s134*s267*s34)
     & *zba4267
     & *(za(q4,q1)*za(q3,q5)+za(q4,q5)*za(q3,q1))*zb(q2,q6)

     & -m/(two*zb(q4,q5)*s267*s345*s34)*(
     & -3._dp*za(q4,q1)*zba4267+za(q3,q1)*zba3267
     & +za(q5,q1)*zba5267)*zb(q2,q6)*za(q3,q5)

     & -m/(two*zb(q4,q5)*s167*s345*s34)*(
     & +3._dp*zb(q4,q2)*zba6174-zb(q3,q2)*zba6173
     & -zb(q5,q2)*zba6175)*za(q3,q5)*za(q1,q7)

      tmp(2,2)=
     & -m/(zb(q3,q5)*zb(q2,q5)*za(q4,q3)*s134*s34)*(
     & -(zb(q1,q4)*za(q4,q1)-zb(q1,q3)*za(q3,q1))
     & *zb(q3,q2)*zb(q2,q6)*za(q1,q7))

     & -m/(zb(q3,q5)*zb(q2,q5)*za(q4,q3)*s167*s34)*(
     & -(zb(q3,q2)*zba6173-zb(q4,q2)*zba6174)
     & *zb(q3,q2)*za(q1,q7))

     & -m/(zb(q3,q5)*za(q4,q3)*s134*s267*s34)*
     & (zb(q3,q1)*za(q3,q1)-zb(q4,q1)*za(q4,q1))
     & *zba3267*zb(q2,q6)*za(q1,q5)

     & -m/(zb(q3,q5)*za(q4,q3)*s267*s345*s34)*(
     & +zb(q4,q5)*zba3267*zb(q2,q6)*za(q4,q5)*za(q1,q5))

     & -m/(zb(q3,q5)*za(q4,q3)*s167*s345*s34)*(
     & +zb(q4,q5)*zb(q3,q2)*zba6175*za(q4,q5)*za(q1,q7))

     & -m/(za(q4,q3)*s267*s345*s34)*(
     & -zba3267*zb(q2,q6)*za(q3,q5)*za(q1,q5))

     & -m/(za(q4,q3)*s167*s345*s34)*(
     & -zb(q3,q2)*zba6175*za(q3,q5)*za(q1,q7))

     & -m/(zb(q3,q5)*zb(q2,q5)*za(q4,q3)**2*s134)*(
     & -zb(q3,q2)*zb(q2,q6)*(za(q4,q1)*za(q3,q7)+za(q4,q7)*za(q3,q1)))

     & -m/(zb(q3,q5)*za(q4,q3)*p3Dq5*s267*s345*two)*(
     & +(zba4267*za(q4,q1)-zba3267*za(q3,q1))
     & *zb(q2,q6)*zb(q4,q3)*za(q4,q5)*bm)

     & -m/(zb(q3,q5)*za(q4,q3)*p3Dq5*s167*s345*two)*(
     & (zb(q3,q2)*zba6173-zb(q4,q2)*zba6174)
     & *zb(q4,q3)*za(q4,q5)*za(q1,q7)*bm)

     & -m/(zb(q3,q5)*za(q4,q3)**2*s134*s267)*(
     & +zba3267*zb(q2,q6)
     & *(za(q4,q1)*za(q3,q5)+za(q4,q5)*za(q3,q1)))

     & -m/(two*zb(q3,q5)*za(q4,q3)**2*s267*s345)*(
     & (zba4267*za(q1,q4)-3._dp*zba3267*za(q1,q3)
     & +zba5267*za(q1,q5))*zb(q2,q6)*za(q4,q5))

     & -m/(two*zb(q3,q5)*za(q4,q3)**2*s167*s345)*(
     & (zb(q4,q2)*zba6174-3._dp*zb(q3,q2)*zba6173
     & +zb(q5,q2)*zba6175)*za(q4,q5)*za(q1,q7))

     & -m/(za(q4,q3)*p3Dq5*s267*s345*two)*(
     & +zba4267*zb(q2,q6)*za(q4,q5)*za(q1,q5))

     & -m/(za(q4,q3)*p3Dq5*s167*s345*two)*(
     & +zb(q4,q2)*zba6175*za(q4,q5)*za(q1,q7))

      tmp(2,1)=
     & -zba3267*zb(q2,q6)*za(q4,q5)*za(q1,q5)*bm
     & /(two*p3Dq5*s267*s345)

     & -zb(q3,q2)*zba6175*za(q4,q5)*za(q1,q7)*bm
     & /(two*p3Dq5*s167*s345)

     & -zb(q3,q2)*zb(q3,q1)*zb(q2,q6)*za(q4,q1)*za(q1,q7)
     & /(zb(q3,q5)*zb(q2,q5)*s134*s34)

     & -(zb(q3,q2)**2*zba6174*za(q1,q7))
     & /(zb(q3,q5)*zb(q2,q5)*s167*s34)

     & +zba3267*zb(q3,q1)*zb(q2,q6)*za(q4,q1)*za(q1,q5)
     & /(zb(q3,q5)*s134*s267*s34)

     & +zba3267*zb(q2,q6)*zba3451*za(q4,q5)
     & /(zb(q3,q5)*s267*s345*s34)

     & -zb(q4,q3)*zba3267*zb(q2,q6)*za(q4,q1)*za(q4,q5)*bm
     & /(zb(q3,q5)*two*p3Dq5*s267*s345)

     & +(zba6174*zb(q4,q3)+zba6175*zb(q5,q3))
     & *zb(q3,q2)*za(q4,q5)*za(q1,q7)
     & /(zb(q3,q5)*s167*s345*s34)

     & +zb(q4,q3)*zb(q3,q2)*zba6174*za(q4,q5)*za(q1,q7)*bm
     & /(zb(q3,q5)*two*p3Dq5*s167*s345)

     & -zb(q3,q2)*zb(q2,q6)*za(q4,q1)*za(q4,q7)
     & /(zb(q3,q5)*zb(q2,q5)*za(q4,q3)*s134)

     & +zba3267*zb(q2,q6)*za(q4,q1)*za(q4,q5)
     & /(zb(q3,q5)*za(q4,q3)*s134*s267)

      if (h == 1) then
      qcdb(1,1,h)=+tmp(1,1)
      qcdb(1,2,h)=+tmp(1,2)
      qcdb(2,1,h)=+tmp(2,1)
      qcdb(2,2,h)=+tmp(2,2)
      qcda(1,1,h)=qedi(1,1,h)+qedf(1,1,h)-qcdb(1,1,h)
      qcda(1,2,h)=qedi(1,2,h)+qedf(1,2,h)-qcdb(1,2,h)
      qcda(2,1,h)=qedi(2,1,h)+qedf(2,1,h)-qcdb(2,1,h)
      qcda(2,2,h)=qedi(2,2,h)+qedf(2,2,h)-qcdb(2,2,h)
      elseif (h == 2) then
      qcda(1,1,h)=-tmp(2,2)
      qcda(1,2,h)=-tmp(1,2)
      qcda(2,1,h)=-tmp(2,1)
      qcda(2,2,h)=-tmp(1,1)
      qcdb(1,1,h)=qedi(1,1,h)+qedf(1,1,h)-qcda(1,1,h)
      qcdb(1,2,h)=qedi(1,2,h)+qedf(1,2,h)-qcda(1,2,h)
      qcdb(2,1,h)=qedi(2,1,h)+qedf(2,1,h)-qcda(2,1,h)
      qcdb(2,2,h)=qedi(2,2,h)+qedf(2,2,h)-qcda(2,2,h)
      endif

      return

      end
