      subroutine box1(k1,k2,k3,k4,k5,k6,za,zb,app,apm,amp,amm)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'Higgsint.f'
      real(dp):: mtsq
      integer:: k1,k2,k3,k4,k5,k6
      complex(dp):: app,apm,amp,amm,boxamp,zab2
      real(dp):: s134
      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
c--- end statement functions

      mtsq=mt**2

c--- Box 1
      app=czip
      amm=czip
c      apm=czip
c      amp=czip

      if (Higgsint) return

c--- Improved expressions
      s134=s(k1,k3)+s(k1,k4)+s(k3,k4)
      boxamp =s(k1,k2)*zb(k2,k6)*za(k1,k3)
     & *(mtsq-s134)**3/zab2(k1,k3,k4,k2)**4*(mtsq*zb(k4,k2)*za(k1,k5)
     & +za(k1,k3)*zb(k3,k4)*za(k5,k6)*zb(k2,k6))
      apm=boxamp*(-half*im)/(s(k3,k4)*s(k5,k6))

      boxamp=(mtsq-s134)*s(k1,k2)/zab2(k2,k3,k4,k1)**4
     & *(mtsq*zb(k4,k1)*za(k2,k5)+zab2(k2,k1,k3,k4)*zab2(k5,k3,k4,k1))
     & *(mtsq*zb(k1,k6)-zb(k5,k6)*zab2(k5,k2,k6,k1))
     & *(mtsq*za(k2,k3)+za(k3,k4)*zab2(k2,k1,k3,k4))

      amp=boxamp*(-half*im)/(s(k3,k4)*s(k5,k6))

      return
      end
      
c--- FORM output
c      s134=za(k1,k3)*zb(k3,k1)+za(k1,k4)*zb(k4,k1)+za(k3,k4)*zb(k4,k3)
c      ga=-(s134-mtsq)/(za(k1,k3)*zb(k3,k2)+za(k1,k4)*zb(k4,k2))
c      de=-(s134-mtsq)/(za(k2,k3)*zb(k3,k1)+za(k2,k4)*zb(k4,k1))
c      boxamp =  + za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k4)*zb(
c     &    k2,k6)*ga**3 + za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k2,
c     &    k4)*zb(k2,k6)*ga**4 + za(k1,k2)*za(k1,k3)*za(k3,k5)*zb(k1,k2)
c     &    *zb(k2,k6)*zb(k3,k4)*ga**3

c      apm=boxamp*(-half*im)/(s(k3,k4)*s(k5,k6))
c
c      boxamp =  + za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k4)*zb(
c     &    k1,k6)*de**2 - za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,
c     &    k4)*zb(k2,k6)*de + za(k1,k2)*za(k1,k3)*za(k2,k5)*zb(k1,k2)*
c     &    zb(k1,k4)*zb(k1,k6)*de**3 - za(k1,k2)*za(k1,k3)*za(k2,k5)*zb(
c     &    k1,k2)*zb(k1,k4)*zb(k2,k6)*de**2 + za(k1,k2)*za(k1,k3)*za(k3,
c     &    k5)*zb(k1,k2)*zb(k1,k6)*zb(k3,k4)*de**2 - za(k1,k2)*za(k1,k3)
c     &    *za(k3,k5)*zb(k1,k2)*zb(k2,k6)*zb(k3,k4)*de + za(k1,k2)*za(k1
c     &    ,k5)*za(k2,k3)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*de**3 - za(k1,k2
c     &    )*za(k1,k5)*za(k2,k3)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6)*de**2 + 
c     &    za(k1,k2)*za(k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*
c     &    de**4 - za(k1,k2)*za(k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(
c     &    k2,k6)*de**3 + za(k1,k2)*za(k2,k3)*za(k3,k5)*zb(k1,k2)*zb(k1,
c     &    k6)*zb(k3,k4)*de**3 - za(k1,k2)*za(k2,k3)*za(k3,k5)*zb(k1,k2)
c     &    *zb(k2,k6)*zb(k3,k4)*de**2

c      amp=boxamp*(-half*im)/(s(k3,k4)*s(k5,k6))

