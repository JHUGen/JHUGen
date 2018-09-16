      subroutine triangle12(k1,k2,k3,k4,k5,k6,za,zb,app,apm)
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
      integer:: k1,k2,k3,k4,k5,k6,i1,i2
      complex(dp):: app,apm,triamp,iza,A134,A134c,s134h,s234h
c--- statement functions
      iza(i1,i2)=cone/za(i1,i2)
c--- end statement functions

      mtsq=mt**2
      apm=czip

c--- Triangle 12
      A134=-(za(k1,k3)*zb(k3,k2)+za(k1,k4)*zb(k4,k2))
      A134c=conjg(A134)
      s134h=s(k1,k3)+s(k1,k4)+s(k3,k4)+mtsq
      s234h=s(k2,k3)+s(k2,k4)+s(k3,k4)+mtsq

      triamp =  + iza(k1,k2)*mtsq*A134**(-2) * ( za(k1,k3)*za(k1,k5)*
     &    zb(k1,k2)*zb(k2,k4)*zb(k2,k6)*s234h + za(k1,k3)*za(k1,k5)*zb(
     &    k1,k2)*zb(k2,k4)*zb(k2,k6)*s134h )
      triamp = triamp + iza(k1,k2)*mtsq*A134**(-1) * ( za(k1,k3)*za(k1,
     &    k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6) + za(k1,k3)*za(k1,k5)*zb(k1
     &    ,k2)*zb(k1,k6)*zb(k2,k4) + za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(
     &    k2,k4)*zb(k2,k6) + za(k1,k5)*za(k2,k3)*zb(k1,k2)*zb(k2,k4)*
     &    zb(k2,k6) - 2._dp*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4)*zb(
     &    k4,k6) )
      triamp = triamp + iza(k1,k2)*mtsq*A134c**(-2) * ( za(k2,k3)*za(k2
     &    ,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*s234h + za(k2,k3)*za(k2,k5
     &    )*zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*s134h )
      triamp = triamp + iza(k1,k2)*mtsq*A134c**(-1) * ( za(k1,k3)*za(k2
     &    ,k5)*zb(k1,k2)*zb(k1,k4)*zb(k1,k6) + za(k1,k5)*za(k2,k3)*zb(
     &    k1,k2)*zb(k1,k4)*zb(k1,k6) + za(k2,k3)*za(k2,k5)*zb(k1,k2)*
     &    zb(k1,k4)*zb(k2,k6) + za(k2,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k6)
     &    *zb(k2,k4) - 2._dp*za(k2,k5)*za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(
     &    k4,k6) )

      app=triamp*(-half*im)/(s(k3,k4)*s(k5,k6))

      if (Higgsint) return

      triamp =  + A134**(-4) * (  - za(k1,k2)*za(k1,k3)*za(k1,k5)*zb(k1
     &    ,k2)*zb(k2,k4)*zb(k2,k6)*s234h**3 - za(k1,k2)*za(k1,k3)*za(k1
     &    ,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6)*s134h**3 )
      triamp = triamp + A134**(-3) * (  - za(k1,k2)*za(k1,k3)*za(k1,k5)
     &    *zb(k1,k2)*zb(k1,k4)*zb(k2,k6)*s134h**2 - za(k1,k2)*za(k1,k3)
     &    *za(k1,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4)*s134h**2 + za(k1,k2)
     &    *za(k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6)*s134h**2
     &     - za(k1,k2)*za(k1,k5)*za(k2,k3)*zb(k1,k2)*zb(k2,k4)*zb(k2,k6
     &    )*s234h**2 + za(k1,k2)*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4
     &    )*zb(k4,k6)*s234h**2 + za(k1,k2)*za(k1,k5)*za(k3,k4)*zb(k1,k2
     &    )*zb(k2,k4)*zb(k4,k6)*s134h**2 )
      triamp = triamp + A134**(-2) * (  - za(k1,k2)*za(k1,k3)*za(k1,k5)
     &    *zb(k1,k2)*zb(k1,k4)*zb(k1,k6)*s134h + za(k1,k2)*za(k1,k3)*
     &    za(k2,k5)*zb(k1,k2)*zb(k1,k4)*zb(k2,k6)*s134h + za(k1,k2)*za(
     &    k1,k3)*za(k2,k5)*zb(k1,k2)*zb(k1,k6)*zb(k2,k4)*s134h + za(k1,
     &    k2)*za(k1,k5)*za(k3,k4)*zb(k1,k2)*zb(k1,k4)*zb(k4,k6)*s134h
     &     - za(k1,k2)*za(k2,k5)*za(k3,k4)*zb(k1,k2)*zb(k2,k4)*zb(k4,k6
     &    )*s134h )
      triamp = triamp + A134**(-1) * ( za(k1,k2)*za(k1,k3)*za(k2,k5)*
     &    zb(k1,k2)*zb(k1,k4)*zb(k1,k6) - za(k1,k2)*za(k2,k5)*za(k3,k4)
     &    *zb(k1,k2)*zb(k1,k4)*zb(k4,k6) )
      triamp = triamp + mtsq*A134**(-3)*A134c * ( 2._dp*za(k1,k3)*za(k1,
     &    k5)*zb(k2,k4)*zb(k2,k6)*s234h + 2._dp*za(k1,k3)*za(k1,k5)*zb(
     &    k2,k4)*zb(k2,k6)*s134h )
      triamp = triamp + mtsq*A134**(-2) * (  - za(k1,k3)*za(k2,k5)*zb(
     &    k1,k4)*zb(k2,k6)*s234h - za(k1,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,
     &    k6)*s134h - za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4)*s234h - 
     &    za(k1,k5)*za(k2,k3)*zb(k1,k6)*zb(k2,k4)*s134h )
      triamp = triamp + mtsq*A134**(-2)*A134c * ( za(k1,k3)*za(k1,k5)*
     &    zb(k1,k4)*zb(k2,k6) + za(k1,k3)*za(k1,k5)*zb(k1,k6)*zb(k2,k4)
     &     - za(k1,k3)*za(k2,k5)*zb(k2,k4)*zb(k2,k6) + za(k1,k5)*za(k2,
     &    k3)*zb(k2,k4)*zb(k2,k6) - 2._dp*za(k1,k5)*za(k3,k4)*zb(k2,k4)*
     &    zb(k4,k6) )
      triamp = triamp + mtsq*A134**(-1) * (  - za(k1,k3)*za(k2,k5)*zb(
     &    k1,k4)*zb(k1,k6) - za(k1,k5)*za(k2,k3)*zb(k1,k4)*zb(k1,k6) - 
     &    za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k2,k6) + za(k2,k3)*za(k2,k5)
     &    *zb(k1,k6)*zb(k2,k4) + 2._dp*za(k2,k5)*za(k3,k4)*zb(k1,k4)*zb(
     &    k4,k6) )

      apm=triamp*(-half*im)/(s(k3,k4)*s(k5,k6))

      return
      end
