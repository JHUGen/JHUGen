      subroutine triangle9new(k1,k2,k3,k4,k5,k6,za,zb,app,apm)
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
      complex(dp):: app,apm,triamp,za56b,iza56b,
     & iza,izb,t2sum,t2prod,t1256sum
c--- statement functions

      iza(i1,i2)=cone/za(i1,i2)
      izb(i1,i2)=cone/zb(i1,i2)
      za56b(i1,i2)=za(i1,k5)*zb(k5,i2)+za(i1,k6)*zb(k6,i2)
      iza56b(i1,i2)=cone/za56b(i1,i2)
c--- end statement functions

      mtsq=mt**2
      apm=czip

c--- Triangle 9
      t2sum=-(-za(k1,k2)*zb(k2,k1)/za56b(k1,k1)*za(k5,k6)*zb(k6,k5)
     &        -za(k1,k2)*zb(k2,k1)/za56b(k1,k1)*mtsq
     &        +za56b(k2,k2))
      t2prod=-za(k1,k2)*zb(k2,k1)/za56b(k1,k1)*mtsq
     & *(za56b(k2,k2)
     &  -za(k5,k6)*zb(k6,k5)/za56b(k1,k1)*za(k1,k2)*zb(k2,k1))
      t1256sum=-(za(k1,k2)*zb(k2,k1)/za56b(k1,k1)*za(k5,k6)*zb(k6,k5)
     &          +za(k1,k2)*zb(k2,k1)/za56b(k1,k1)*mtsq
     &          +za(k1,k2)*zb(k2,k1))


      triamp =  + mtsq * (  - za(k1,k3)*za(k1,k5)**2*zb(k1,k2)*zb(k2,k4
     &    )*zb(k5,k6)*iza(k1,k2)**2*za56b(k2,k1)*iza56b(k1,k1)*iza56b(
     &    k1,k2) + za(k1,k3)*za(k1,k5)**2*zb(k1,k2)*zb(k2,k4)*zb(k5,k6)
     &    *iza(k1,k2)**2*iza56b(k1,k2)**2*t1256sum - za(k1,k3)*za(k1,k5
     &    )**2*zb(k1,k2)*zb(k5,k6)*iza(k1,k2)**3*za56b(k1,k4)*za56b(k2,
     &    k1)*iza56b(k1,k1)*iza56b(k1,k2) - za(k1,k3)*za(k1,k5)**2*zb(
     &    k1,k2)*zb(k5,k6)*iza(k1,k2)**3*za56b(k1,k4)*iza56b(k1,k2)**2*
     &    t2sum + za(k1,k3)*za(k1,k5)**2*zb(k1,k2)*zb(k5,k6)*iza(k1,k2)
     &    **3*za56b(k2,k4)*iza56b(k1,k2) - za(k1,k3)*za(k1,k5)**2*zb(k1
     &    ,k4)*zb(k5,k6)*iza(k1,k2)**3*za56b(k2,k2)*iza56b(k1,k2) - za(
     &    k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)**3*iza56b(
     &    k1,k2)*t1256sum - za(k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*
     &    iza(k1,k2)**3*iza56b(k1,k2)*t2sum + za(k1,k3)*za(k1,k5)**2*
     &    zb(k2,k4)*zb(k5,k6)*iza(k1,k2)**3*za56b(k1,k1)*iza56b(k1,k2)
     &    **2*t1256sum + za(k1,k3)*za(k1,k5)*za(k2,k5)*zb(k1,k2)*zb(k2,
     &    k4)*zb(k5,k6)*iza(k1,k2)**2*iza56b(k1,k2) )
      triamp = triamp + mtsq * ( za(k1,k3)*za(k1,k5)*za(k2,k5)*zb(k1,k4
     &    )*zb(k5,k6)*iza(k1,k2)**3 - za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(
     &    k1,k6)*zb(k2,k4)*iza(k1,k2)*iza56b(k1,k2) - za(k1,k3)*za(k1,
     &    k5)*zb(k1,k4)*zb(k1,k6)*iza(k1,k2)**2 - za(k1,k3)*za(k1,k5)*
     &    zb(k1,k6)*zb(k2,k4)*iza(k1,k2)**2*za56b(k2,k1)*iza56b(k1,k1)
     &     + za(k1,k3)*za(k1,k5)*zb(k1,k6)*zb(k2,k4)*iza(k1,k2)**2*
     &    za56b(k2,k2)*iza56b(k1,k2) + za(k1,k3)*za(k1,k5)*zb(k1,k6)*
     &    zb(k2,k4)*iza(k1,k2)**2*iza56b(k1,k2)*t1256sum - za(k1,k3)*
     &    za(k1,k5)*zb(k1,k6)*iza(k1,k2)**3*za56b(k1,k4)*za56b(k2,k1)*
     &    iza56b(k1,k1) - za(k1,k3)*za(k1,k5)*zb(k1,k6)*iza(k1,k2)**3*
     &    za56b(k1,k4)*iza56b(k1,k2)*t2sum + za(k1,k3)*za(k1,k5)*zb(k1,
     &    k6)*iza(k1,k2)**3*za56b(k2,k4) - za(k1,k5)**2*zb(k1,k2)*zb(k2
     &    ,k4)*zb(k5,k6)*iza(k1,k2)*za56b(k3,k1)*iza56b(k1,k1)*iza56b(
     &    k1,k2) - za(k1,k5)**2*zb(k1,k2)*zb(k5,k6)*iza(k1,k2)**2*
     &    za56b(k1,k4)*za56b(k3,k1)*iza56b(k1,k1)*iza56b(k1,k2) - za(k1
     &    ,k5)**2*zb(k2,k4)*zb(k5,k6)*iza(k1,k2)**2*za56b(k3,k1)*
     &    iza56b(k1,k2) )
      triamp = triamp + mtsq * (  - za(k1,k5)*zb(k1,k4)*zb(k1,k6)*iza(
     &    k1,k2)**2*izb(k1,k2)*za56b(k2,k2)*za56b(k3,k1)*iza56b(k2,k1)
     &     - za(k1,k5)*zb(k1,k4)*zb(k1,k6)*iza(k1,k2)**2*izb(k1,k2)*
     &    za56b(k3,k1)*iza56b(k2,k1)*t1256sum - za(k1,k5)*zb(k1,k4)*zb(
     &    k1,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k3,k1)*iza56b(k2,k1)*
     &    t2sum - 2._dp*za(k1,k5)*zb(k1,k6)*zb(k2,k4)*iza(k1,k2)*za56b(
     &    k3,k1)*iza56b(k1,k1) - 2._dp*za(k1,k5)*zb(k1,k6)*iza(k1,k2)**2
     &    *za56b(k1,k4)*za56b(k3,k1)*iza56b(k1,k1) + za(k1,k5)*zb(k1,k6
     &    )*iza(k1,k2)**2*za56b(k2,k4)*za56b(k3,k1)*iza56b(k2,k1) + za(
     &    k2,k5)*zb(k1,k4)*zb(k1,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,
     &    k1)*za56b(k3,k1)*iza56b(k2,k1)**2*t1256sum + za(k2,k5)*zb(k1,
     &    k4)*zb(k1,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,k1)*za56b(k3,
     &    k1)*iza56b(k2,k1)**2*t2sum + za(k2,k5)*zb(k1,k4)*zb(k1,k6)*
     &    iza(k1,k2)**2*izb(k1,k2)*za56b(k1,k2)*za56b(k3,k1)*iza56b(k2,
     &    k1) + za(k2,k5)*zb(k1,k6)*zb(k2,k4)*iza(k1,k2)*za56b(k3,k1)*
     &    iza56b(k2,k1) )
      triamp = triamp + mtsq**2 * (  - za(k1,k3)*za(k1,k5)*zb(k1,k2)*
     &    zb(k1,k6)*zb(k2,k4)*iza(k1,k2)*iza56b(k1,k1)*iza56b(k1,k2) - 
     &    za(k1,k3)*za(k1,k5)*zb(k1,k2)*zb(k1,k6)*iza(k1,k2)**2*za56b(
     &    k1,k4)*iza56b(k1,k1)*iza56b(k1,k2) - za(k1,k3)*za(k1,k5)*zb(
     &    k1,k6)*zb(k2,k4)*iza(k1,k2)**2*iza56b(k1,k2) )
      triamp = triamp + za(k1,k3)*za(k1,k5)**2*zb(k1,k2)*zb(k2,k4)*zb(
     & k5,k6)*iza(k1,k2)**2*za56b(k1,k1)*iza56b(k1,k2)**2*t1256sum - 
     &    za(k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)**3*
     &    za56b(k1,k1)*iza56b(k1,k2)*t2sum - za(k1,k3)*za(k1,k5)**2*zb(
     &    k2,k4)*zb(k5,k6)*iza(k1,k2)**3*za56b(k1,k1)*za56b(k2,k2)*
     &    iza56b(k1,k2)**2*t1256sum + za(k1,k3)*za(k1,k5)**2*zb(k2,k4)*
     &    zb(k5,k6)*iza(k1,k2)**3*za56b(k1,k1)*iza56b(k1,k2)**2*t2prod
     &     - za(k1,k3)*za(k1,k5)**2*zb(k2,k4)*zb(k5,k6)*iza(k1,k2)**3*
     &    za56b(k1,k1)*iza56b(k1,k2)**2*t1256sum**2 + za(k1,k3)*za(k1,
     &    k5)**2*zb(k2,k4)*zb(k5,k6)*iza(k1,k2)**3*za56b(k2,k1)*iza56b(
     &    k1,k2)*t1256sum + za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)
     &    **4*za56b(k1,k1)*za56b(k1,k4)*iza56b(k1,k2)**2*t2prod - za(k1
     &    ,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**4*za56b(k1,k1)*za56b(
     &    k1,k4)*iza56b(k1,k2)**2*t2sum**2 + za(k1,k3)*za(k1,k5)**2*zb(
     &    k5,k6)*iza(k1,k2)**4*za56b(k1,k1)*za56b(k2,k4)*iza56b(k1,k2)*
     &    t2sum
      triamp = triamp - za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**4*
     & za56b(k1,k4)*za56b(k2,k1)*iza56b(k1,k2)*t2sum - za(k1,k5)**2*zb(
     &    k1,k2)*zb(k2,k4)*zb(k5,k6)*iza(k1,k2)*za56b(k3,k1)*iza56b(k1,
     &    k2) - za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)**2*za56b(k3
     &    ,k1) - za(k1,k5)**2*zb(k2,k4)*zb(k5,k6)*iza(k1,k2)**2*za56b(
     &    k2,k1)*za56b(k3,k1)*iza56b(k1,k1) + za(k1,k5)**2*zb(k2,k4)*
     &    zb(k5,k6)*iza(k1,k2)**2*za56b(k2,k2)*za56b(k3,k1)*iza56b(k1,
     &    k2) + za(k1,k5)**2*zb(k2,k4)*zb(k5,k6)*iza(k1,k2)**2*za56b(k3
     &    ,k1)*iza56b(k1,k2)*t1256sum - za(k1,k5)**2*zb(k5,k6)*iza(k1,
     &    k2)**3*za56b(k1,k4)*za56b(k2,k1)*za56b(k3,k1)*iza56b(k1,k1)
     &     - za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**3*za56b(k1,k4)*za56b(k3
     &    ,k1)*iza56b(k1,k2)*t2sum + za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)
     &    **3*za56b(k2,k4)*za56b(k3,k1)

      app=triamp*(-half*im)/(s(k3,k4)*s(k5,k6))
      if (Higgsint) return
      triamp =  + mtsq * (  - za(k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6
     &    )*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,k1)*iza56b(k1,k2)**3*
     &    t1256sum**2 + za(k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(
     &    k1,k2)**2*izb(k1,k2)*za56b(k1,k1)*iza56b(k1,k2)**3*t2sum**2
     &     + za(k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)**2*
     &    izb(k1,k2)*za56b(k2,k1)**2*iza56b(k1,k1)*iza56b(k1,k2) + za(
     &    k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)**2*izb(k1,
     &    k2)*za56b(k2,k1)*iza56b(k1,k2)**2*t1256sum + 2._dp*za(k1,k3)*
     &    za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*
     &    za56b(k2,k1)*iza56b(k1,k2)**2*t2sum - za(k1,k3)*za(k1,k5)*za(
     &    k2,k5)*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,
     &    k1)*iza56b(k1,k2)**2*t1256sum - za(k1,k3)*za(k1,k5)*za(k2,k5)
     &    *zb(k1,k4)*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,k1)*
     &    iza56b(k1,k2)**2*t2sum - za(k1,k3)*za(k1,k5)*za(k2,k5)*zb(k1,
     &    k4)*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k2,k1)*iza56b(k1
     &    ,k2) )
      triamp = triamp + mtsq * ( za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k1,k6
     &    )*iza(k1,k2)*izb(k1,k2)*za56b(k1,k1)*iza56b(k1,k2)**2*
     &    t1256sum + za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k1,k6)*iza(k1,k2)
     &    *izb(k1,k2)*za56b(k1,k1)*iza56b(k1,k2)**2*t2sum + za(k1,k3)*
     &    za(k1,k5)*zb(k1,k4)*zb(k1,k6)*iza(k1,k2)*izb(k1,k2)*za56b(k2,
     &    k1)*iza56b(k1,k2) - za(k1,k3)*za(k1,k5)*zb(k1,k6)*iza(k1,k2)
     &    **2*izb(k1,k2)*za56b(k1,k1)*za56b(k1,k4)*iza56b(k1,k2)**3*
     &    t1256sum**2 + za(k1,k3)*za(k1,k5)*zb(k1,k6)*iza(k1,k2)**2*
     &    izb(k1,k2)*za56b(k1,k1)*za56b(k1,k4)*iza56b(k1,k2)**3*
     &    t2sum**2 - za(k1,k3)*za(k1,k5)*zb(k1,k6)*iza(k1,k2)**2*izb(k1
     &    ,k2)*za56b(k1,k1)*za56b(k2,k4)*iza56b(k1,k2)**2*t1256sum - 
     &    za(k1,k3)*za(k1,k5)*zb(k1,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(
     &    k1,k1)*za56b(k2,k4)*iza56b(k1,k2)**2*t2sum + za(k1,k3)*za(k1,
     &    k5)*zb(k1,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,k4)*za56b(k2,
     &    k1)**2*iza56b(k1,k1)*iza56b(k1,k2) + za(k1,k3)*za(k1,k5)*zb(
     &    k1,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,k4)*za56b(k2,k1)*
     &    iza56b(k1,k2)**2*t1256sum )
      triamp = triamp + mtsq * ( 2._dp*za(k1,k3)*za(k1,k5)*zb(k1,k6)*
     &    iza(k1,k2)**2*izb(k1,k2)*za56b(k1,k4)*za56b(k2,k1)*iza56b(k1,
     &    k2)**2*t2sum - za(k1,k3)*za(k1,k5)*zb(k1,k6)*iza(k1,k2)**2*
     &    izb(k1,k2)*za56b(k2,k1)*za56b(k2,k4)*iza56b(k1,k2) + 2._dp*za(
     &    k1,k5)**2*za(k2,k3)*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)*za56b(k2,
     &    k1)*iza56b(k1,k1)*iza56b(k1,k2) + za(k1,k5)**2*za(k2,k3)*zb(
     &    k1,k4)*zb(k5,k6)*iza(k1,k2)*iza56b(k1,k2)**2*t2sum + za(k1,k5
     &    )**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)*izb(k1,k2)*za56b(k2,k1)*
     &    za56b(k3,k1)*iza56b(k1,k1)*iza56b(k1,k2) + za(k1,k5)**2*zb(k1
     &    ,k4)*zb(k5,k6)*iza(k1,k2)*izb(k1,k2)*za56b(k3,k1)*iza56b(k1,
     &    k2)**2*t1256sum + za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)
     &    *izb(k1,k2)*za56b(k3,k1)*iza56b(k1,k2)**2*t2sum - za(k1,k5)*
     &    za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)*iza56b(k1,
     &    k2) + za(k1,k5)*za(k2,k3)*zb(k1,k4)*zb(k1,k6)*iza(k1,k2)*izb(
     &    k1,k2)*za56b(k2,k1)*iza56b(k1,k1) + za(k1,k5)*za(k2,k3)*zb(k1
     &    ,k4)*zb(k1,k6)*iza56b(k1,k2) )
      triamp = triamp + mtsq * ( 2._dp*za(k1,k5)*za(k2,k3)*zb(k1,k6)*
     &    iza(k1,k2)*za56b(k1,k4)*za56b(k2,k1)*iza56b(k1,k1)*iza56b(k1,
     &    k2) + za(k1,k5)*za(k2,k3)*zb(k1,k6)*iza(k1,k2)*za56b(k1,k4)*
     &    iza56b(k1,k2)**2*t2sum - za(k1,k5)*za(k2,k3)*zb(k1,k6)*iza(k1
     &    ,k2)*za56b(k2,k4)*iza56b(k1,k2) + za(k1,k5)*zb(k1,k4)*zb(k1,
     &    k6)*iza(k1,k2)*izb(k1,k2)**2*za56b(k2,k1)*za56b(k3,k1)*
     &    iza56b(k1,k1) + za(k1,k5)*zb(k1,k6)*iza(k1,k2)*izb(k1,k2)*
     &    za56b(k1,k4)*za56b(k2,k1)*za56b(k3,k1)*iza56b(k1,k1)*iza56b(
     &    k1,k2) + za(k1,k5)*zb(k1,k6)*iza(k1,k2)*izb(k1,k2)*za56b(k1,
     &    k4)*za56b(k3,k1)*iza56b(k1,k2)**2*t1256sum + za(k1,k5)*zb(k1,
     &    k6)*iza(k1,k2)*izb(k1,k2)*za56b(k1,k4)*za56b(k3,k1)*iza56b(k1
     &    ,k2)**2*t2sum - za(k2,k3)*za(k2,k5)*zb(k1,k4)*zb(k1,k6)*iza(
     &    k1,k2)*izb(k1,k2) - za(k2,k5)*zb(k1,k4)*zb(k1,k6)*iza(k1,k2)*
     &    izb(k1,k2)**2*za56b(k3,k1) )
      triamp = triamp + mtsq**2 * ( za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k1
     &    ,k6)*iza(k1,k2)*izb(k1,k2)*za56b(k2,k1)*iza56b(k1,k1)*iza56b(
     &    k1,k2) + za(k1,k3)*za(k1,k5)*zb(k1,k4)*zb(k1,k6)*iza(k1,k2)*
     &    izb(k1,k2)*iza56b(k1,k2)**2*t1256sum + za(k1,k3)*za(k1,k5)*
     &    zb(k1,k4)*zb(k1,k6)*iza(k1,k2)*izb(k1,k2)*iza56b(k1,k2)**2*
     &    t2sum + za(k1,k5)*za(k2,k3)*zb(k1,k4)*zb(k1,k6)*iza56b(k1,k1)
     &    *iza56b(k1,k2) )
      triamp = triamp - za(k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(
     & k1,k2)**2*izb(k1,k2)*za56b(k1,k1)**2*iza56b(k1,k2)**3*
     & t1256sum**2 + za(k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,
     &    k2)**2*izb(k1,k2)*za56b(k1,k1)**2*iza56b(k1,k2)**3*t2sum**2
     &     + za(k1,k3)*za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)**2*
     &    izb(k1,k2)*za56b(k1,k1)*za56b(k2,k1)*iza56b(k1,k2)**2*t2sum
     &     - 2._dp*za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**3*izb(k1
     &    ,k2)*za56b(k1,k1)**2*za56b(k1,k4)*iza56b(k1,k2)**4*t1256sum*
     &    t2prod + za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**3*izb(
     &    k1,k2)*za56b(k1,k1)**2*za56b(k1,k4)*iza56b(k1,k2)**4*
     &    t1256sum**3 - 2._dp*za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2
     &    )**3*izb(k1,k2)*za56b(k1,k1)**2*za56b(k1,k4)*iza56b(k1,k2)**4
     &    *t2sum*t2prod + za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**
     &    3*izb(k1,k2)*za56b(k1,k1)**2*za56b(k1,k4)*iza56b(k1,k2)**4*
     &    t2sum**3 + za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**3*
     &    izb(k1,k2)*za56b(k1,k1)**2*za56b(k2,k4)*iza56b(k1,k2)**3*
     &    t1256sum**2
      triamp = triamp - za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**3*
     & izb(k1,k2)*za56b(k1,k1)**2*za56b(k2,k4)*iza56b(k1,k2)**3*
     & t2sum**2 - za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**3*izb(k1
     &    ,k2)*za56b(k1,k1)*za56b(k1,k4)*za56b(k2,k1)*iza56b(k1,k2)**3*
     &    t2prod - za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**3*izb(
     &    k1,k2)*za56b(k1,k1)*za56b(k1,k4)*za56b(k2,k1)*iza56b(k1,k2)**
     &    3*t1256sum**2 + 2._dp*za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*iza(k1,
     &    k2)**3*izb(k1,k2)*za56b(k1,k1)*za56b(k1,k4)*za56b(k2,k1)*
     &    iza56b(k1,k2)**3*t2sum**2 - za(k1,k3)*za(k1,k5)**2*zb(k5,k6)*
     &    iza(k1,k2)**3*izb(k1,k2)*za56b(k1,k1)*za56b(k2,k1)*za56b(k2,
     &    k4)*iza56b(k1,k2)**2*t2sum + za(k1,k3)*za(k1,k5)**2*zb(k5,k6)
     &    *iza(k1,k2)**3*izb(k1,k2)*za56b(k1,k4)*za56b(k2,k1)**2*
     &    iza56b(k1,k2)**2*t2sum + za(k1,k5)**2*za(k2,k3)*zb(k1,k4)*zb(
     &    k5,k6)*iza(k1,k2)*za56b(k1,k1)*iza56b(k1,k2)**2*t2sum + za(k1
     &    ,k5)**2*za(k2,k3)*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)*za56b(k2,k1)
     &    *iza56b(k1,k2)
      triamp = triamp - za(k1,k5)**2*za(k2,k3)*zb(k5,k6)*iza(k1,k2)**2*
     & za56b(k1,k1)*za56b(k1,k4)*iza56b(k1,k2)**3*t2prod + za(k1,k5)**2
     &    *za(k2,k3)*zb(k5,k6)*iza(k1,k2)**2*za56b(k1,k1)*za56b(k1,k4)*
     &    iza56b(k1,k2)**3*t2sum**2 - za(k1,k5)**2*za(k2,k3)*zb(k5,k6)*
     &    iza(k1,k2)**2*za56b(k1,k1)*za56b(k2,k4)*iza56b(k1,k2)**2*
     &    t2sum + za(k1,k5)**2*za(k2,k3)*zb(k5,k6)*iza(k1,k2)**2*za56b(
     &    k1,k4)*za56b(k2,k1)**2*iza56b(k1,k1)*iza56b(k1,k2) + 2._dp*za(
     &    k1,k5)**2*za(k2,k3)*zb(k5,k6)*iza(k1,k2)**2*za56b(k1,k4)*
     &    za56b(k2,k1)*iza56b(k1,k2)**2*t2sum - za(k1,k5)**2*za(k2,k3)*
     &    zb(k5,k6)*iza(k1,k2)**2*za56b(k2,k1)*za56b(k2,k4)*iza56b(k1,
     &    k2) + za(k1,k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)*izb(k1,k2)*
     &    za56b(k1,k1)*za56b(k3,k1)*iza56b(k1,k2)**2*t1256sum + za(k1,
     &    k5)**2*zb(k1,k4)*zb(k5,k6)*iza(k1,k2)*izb(k1,k2)*za56b(k1,k1)
     &    *za56b(k3,k1)*iza56b(k1,k2)**2*t2sum + za(k1,k5)**2*zb(k1,k4)
     &    *zb(k5,k6)*iza(k1,k2)*izb(k1,k2)*za56b(k2,k1)*za56b(k3,k1)*
     &    iza56b(k1,k2)
      triamp = triamp - za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)
     & *za56b(k1,k1)*za56b(k1,k4)*za56b(k3,k1)*iza56b(k1,k2)**3*
     & t1256sum**2 + za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*
     &    za56b(k1,k1)*za56b(k1,k4)*za56b(k3,k1)*iza56b(k1,k2)**3*
     &    t2sum**2 - za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*
     &    za56b(k1,k1)*za56b(k2,k4)*za56b(k3,k1)*iza56b(k1,k2)**2*
     &    t1256sum - za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*
     &    za56b(k1,k1)*za56b(k2,k4)*za56b(k3,k1)*iza56b(k1,k2)**2*t2sum
     &     + za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,
     &    k4)*za56b(k2,k1)**2*za56b(k3,k1)*iza56b(k1,k1)*iza56b(k1,k2)
     &     + za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,
     &    k4)*za56b(k2,k1)*za56b(k3,k1)*iza56b(k1,k2)**2*t1256sum + 2._dp
     &    *za(k1,k5)**2*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k1,k4)
     &    *za56b(k2,k1)*za56b(k3,k1)*iza56b(k1,k2)**2*t2sum - za(k1,k5)
     &    **2*zb(k5,k6)*iza(k1,k2)**2*izb(k1,k2)*za56b(k2,k1)*za56b(k2,
     &    k4)*za56b(k3,k1)*iza56b(k1,k2)

      apm=triamp*(-half*im)/(s(k3,k4)*s(k5,k6))
      return
      end
