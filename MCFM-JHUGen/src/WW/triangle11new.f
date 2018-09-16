      subroutine triangle11new(k1,k2,k3,k4,k5,k6,za,zb,app,apm)
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
      integer:: k1,k2,k3,k4,k5,k6,i1,i2,k34h,k56h
      complex(dp):: app,apm,triamp,iza,izb,
     & s56,s34,s12,dot3456,delta,ga,x,y,t134sum,t234sum,
     & yp1,xm1,xpr,ymr
      parameter(k56h=11,k34h=12)

c--- statement functions
      iza(i1,i2)=cone/za(i1,i2)
      izb(i1,i2)=cone/zb(i1,i2)
c--- end statement functions

      mtsq=mt**2
      apm=czip

c--- Triangle 11
      s12=za(k1,k2)*zb(k2,k1)
      s56=za(k5,k6)*zb(k6,k5)
      s34=za(k3,k4)*zb(k4,k3)
      dot3456=(za(k5,k3)*zb(k3,k5)+za(k5,k4)*zb(k4,k5)
     &        +za(k6,k3)*zb(k3,k6)+za(k6,k4)*zb(k4,k6))/2._dp
      delta=sqrt(dot3456**2-s56*s34)
      ga=dot3456+delta

      x=(s34*(s56+ga)-mtsq*(ga+s34))/(s56*s34-ga**2)
      y=(-s56*(s34+ga)+mtsq*(s56+ga))/(s56*s34-ga**2)
      yp1=y+cone
      xm1=x-cone
      ymr=y-s56/ga
      xpr=x+s34/ga
      t134sum=(x*(ga+za(k1,k56h)*zb(k56h,k1))
     &        +y*(za(k3,k4)*zb(k4,k3)+za(k1,k34h)*zb(k34h,k1))
     &        +za(k3,k4)*zb(k4,k3)-mtsq+za(k1,k34h)*zb(k34h,k1)
     &        +za(k1,k56h)*zb(k56h,k1)*za(k3,k4)*zb(k4,k3)/ga)
      t234sum=(x*(ga+za(k2,k56h)*zb(k56h,k2))
     &        +y*(za(k3,k4)*zb(k4,k3)+za(k2,k34h)*zb(k34h,k2))
     &        +za(k3,k4)*zb(k4,k3)-mtsq+za(k2,k34h)*zb(k34h,k2)
     &        +za(k2,k56h)*zb(k56h,k2)*za(k3,k4)*zb(k4,k3)/ga)

      triamp =  + t234sum * (  - za(k1,k34h)*za(k2,k34h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h
     &    )*iza(k2,k56h)*yp1**2 - za(k1,k34h)*za(k3,k56h)*za(k5,k34h)*
     &    zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*yp1*ymr - 
     &    za(k1,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*
     &    zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,k34h)*yp1*xm1
     &     - za(k1,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k2,k56h)*zb(k4,k34h
     &    )*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k2,k34h)*yp1*xpr
     &     - za(k1,k56h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k4,k56h
     &    )*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*iza(k2,k34h)*xpr**2
     &     - za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h
     &    )*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)*xpr*ymr
     &     - za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,k56h
     &    )*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k2,k56h)*yp1*xpr
     &     - za(k1,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k56h
     &    )*iza(k1,k2)**2*iza(k1,k34h)*xm1*xpr )
      triamp = triamp + t234sum * (  - za(k2,k34h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h)*yp1*
     &    ymr - za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h)*izb(k1,k34h)*yp1
     &    *xm1 - za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*iza(k2,k56h)*yp1*xpr - za(k2,k56h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k2)**2*iza(k2,k34h)*izb(k1,k56h)*xpr*ymr - za(k2,k56h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*yp1*xpr - za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k2,k34h)*xm1*xpr - za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*izb(k1,k56h)*izb(k2,k56h)*yp1*ymr - za(k3
     &    ,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1
     &    ,k2)**2*izb(k1,k56h)*yp1*ymr - za(k3,k34h)*za(k5,k34h)*zb(k2,
     &    k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k2,k56h)*
     &    yp1**2 )
      triamp = triamp + t234sum * (  - za(k3,k34h)*za(k5,k56h)*zb(k2,
     &    k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k2,k56h)*yp1*
     &    xm1 - za(k3,k34h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,
     &    k2)**2*yp1*xm1 - za(k3,k56h)*za(k5,k34h)*zb(k2,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k34h)*xpr*ymr - za(k3,
     &    k56h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*xpr*
     &    ymr - za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k2,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k34h)*izb(k2,k34h)*xm1
     &    *xpr - za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*izb(k1,k34h)*xm1*xpr - za(k3,k56h)*za(k5,
     &    k56h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(
     &    k2,k34h)*xpr**2 )
      triamp = triamp + t234sum**2 * ( za(k1,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)**2*
     &    izb(k1,k34h)*yp1 + za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)**2*izb(k1,k56h)*
     &    xpr + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*iza(k1,k56h)*iza(k2,k56h)*izb(k1,k34h)*
     &    yp1 + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*iza(k1,k34h)*iza(k2,k34h)*izb(k1,k56h)*
     &    xpr + za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)**2*ymr + za(k3,
     &    k34h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k2)**2*iza(k1,k34h)*izb(k1,k56h)*izb(k2,k56h)*yp1 + za(k3,
     &    k34h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k1,k34h)*izb(k1,k56h)*yp1 + za(k3,k34h)*za(k5,k56h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)*xm1
     &     + za(k3,k56h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)
     &    **2*iza(k1,k56h)*izb(k1,k34h)*ymr )
      triamp = triamp + t234sum**2 * ( za(k3,k56h)*za(k5,k56h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(
     &    k1,k34h)**2*xm1 + za(k3,k56h)*za(k5,k56h)*zb(k2,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,k34h)*
     &    izb(k2,k34h)*xpr + za(k3,k56h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,k34h)*xpr )
      triamp = triamp + t234sum**3 * (  - za(k3,k34h)*za(k5,k34h)*zb(k4
     &    ,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)**2*izb(k1,k56h)
     &    **2 - za(k3,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k2)**2*iza(k1,k56h)**2*izb(k1,k34h)**2 )
      triamp = triamp + t134sum * ( za(k1,k34h)*za(k2,k34h)*za(k3,k56h)
     &    *za(k5,k56h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**
     &    2*iza(k1,k56h)**2*izb(k1,k34h)*yp1**2 + za(k1,k34h)*za(k2,
     &    k56h)*za(k3,k56h)*za(k5,k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*iza(k1,k56h)**2*izb(k1,k34h)*yp1*xm1 - 
     &    za(k1,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k34h)*
     &    iza(k1,k2)**2*iza(k1,k56h)*yp1*ymr - za(k1,k34h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2
     &    *iza(k1,k56h)*izb(k1,k34h)*yp1*xm1 - za(k1,k34h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h
     &    )*yp1*xpr + za(k1,k56h)*za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h
     &    )**2*izb(k1,k56h)*xpr*ymr + za(k1,k56h)*za(k2,k56h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k2)**2*iza(k1,k34h)**2*izb(k1,k56h)*xpr**2 - za(k1,k56h)*za(
     &    k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(
     &    k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)*xpr*ymr )
      triamp = triamp + t134sum * (  - za(k1,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*yp1*
     &    xpr - za(k1,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*iza(k1,k34h)*xm1*xpr + za(k2,k34h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)**2*yp1*ymr + 
     &    za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,k34h)*
     &    zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)*yp1*ymr
     &     + za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,k34h
     &    )*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,k34h)*yp1*ymr
     &     + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k2,k34h
     &    )*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,
     &    k34h)**2*yp1*xpr + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k1,k34h)*izb(k1,k56h)**2*yp1*xpr + za(k2,k56h)*za(k3,k34h)*
     &    za(k5,k56h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2
     &    *iza(k1,k34h)*izb(k1,k56h)*xm1*xpr )
      triamp = triamp + t134sum * ( za(k2,k56h)*za(k3,k56h)*za(k5,k56h)
     &    *zb(k1,k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**
     &    2*iza(k1,k56h)*izb(k1,k34h)**2*xm1*xpr + za(k2,k56h)*za(k3,
     &    k56h)*za(k5,k56h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k2)**2*iza(k1,k56h)*izb(k1,k34h)*xm1*xpr - za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(
     &    k1,k56h)*yp1**2 - za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k56h)*yp1*xm1 - za(k3,
     &    k34h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*yp1*
     &    xm1 - za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*izb(k1,k34h)*xpr*ymr - za(k3,k56h)*za(k5,
     &    k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*xpr*ymr - za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k2)**2*izb(k1,k34h)*xpr**2 )
      triamp = triamp + t134sum**2 * (  - za(k1,k34h)*za(k3,k56h)*za(k5
     &    ,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)**2*
     &    izb(k1,k34h)*yp1 - za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)**2*izb(k1,k56h)*
     &    xpr + za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)**2*izb(k1,k56h)
     &    **2*ymr + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k2,k34h)*zb(
     &    k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)**2*izb(k1,
     &    k34h)**2*yp1 + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k56h
     &    )*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)**2*izb(
     &    k1,k56h)**2*xpr + za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k2,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)**2*
     &    izb(k1,k34h)**2*xm1 - za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(
     &    k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)
     &    **2*yp1 - za(k3,k34h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k56h)*
     &    iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)*yp1 )
      triamp = triamp + t134sum**2 * (  - za(k3,k34h)*za(k5,k56h)*zb(k4
     &    ,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)*
     &    xm1 - za(k3,k56h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k2)**2*iza(k1,k56h)*izb(k1,k34h)*ymr - za(k3,k56h)*za(k5,k56h
     &    )*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,
     &    k56h)*izb(k1,k34h)**2*xpr - za(k3,k56h)*za(k5,k56h)*zb(k4,
     &    k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,k34h)*xpr
     &     )
      triamp = triamp + t134sum**3 * (  - za(k3,k34h)*za(k5,k34h)*zb(k4
     &    ,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)**2*izb(k1,k56h)
     &    **2 - za(k3,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k2)**2*iza(k1,k56h)**2*izb(k1,k34h)**2 )
      triamp = triamp + mtsq * ( za(k1,k5)*za(k2,k34h)*za(k3,k34h)*zb(
     &    k1,k2)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*
     &    izb(k1,k56h)*yp1 - za(k1,k5)*za(k2,k34h)*za(k3,k34h)*zb(k1,k4
     &    )*zb(k2,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,
     &    k56h)*yp1 + za(k1,k5)*za(k2,k34h)*za(k3,k56h)*zb(k1,k2)*zb(k4
     &    ,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,k34h)*
     &    yp1 - za(k1,k5)*za(k2,k34h)*za(k3,k56h)*zb(k1,k4)*zb(k2,k34h)
     &    *zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,k34h)*yp1 + 
     &    za(k1,k5)*za(k2,k56h)*za(k3,k34h)*zb(k1,k2)*zb(k4,k56h)*zb(k6
     &    ,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)*xpr - za(k1,k5
     &    )*za(k2,k56h)*za(k3,k34h)*zb(k1,k4)*zb(k2,k56h)*zb(k6,k56h)*
     &    iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)*xpr + za(k1,k5)*za(k2
     &    ,k56h)*za(k3,k56h)*zb(k1,k2)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k2)**2*iza(k1,k56h)*izb(k1,k34h)*xpr - za(k1,k5)*za(k2,k56h)*
     &    za(k3,k56h)*zb(k1,k4)*zb(k2,k56h)*zb(k6,k34h)*iza(k1,k2)**2*
     &    iza(k1,k56h)*izb(k1,k34h)*xpr )
      triamp = triamp + mtsq * (  - za(k1,k5)*za(k3,k34h)*zb(k1,k2)*zb(
     &    k1,k4)*zb(k6,k56h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k56h) - za(
     &    k1,k5)*za(k3,k56h)*zb(k1,k2)*zb(k1,k4)*zb(k6,k34h)*iza(k1,k2)
     &    *iza(k1,k56h)*izb(k1,k34h) + za(k1,k34h)*za(k2,k5)*za(k3,k56h
     &    )*zb(k1,k2)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h
     &    )*izb(k1,k34h)*yp1 + za(k1,k34h)*za(k2,k5)*za(k3,k56h)*zb(k2,
     &    k4)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*yp1 + za(k1,k56h)*
     &    za(k2,k5)*za(k3,k34h)*zb(k1,k2)*zb(k4,k56h)*zb(k6,k56h)*iza(
     &    k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)*xpr + za(k1,k56h)*za(k2,
     &    k5)*za(k3,k34h)*zb(k2,k4)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,
     &    k34h)*xpr + za(k2,k5)*za(k3,k34h)*zb(k1,k2)*zb(k2,k4)*zb(k6,
     &    k56h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k56h) + za(k2,k5)*za(k3,
     &    k34h)*zb(k1,k2)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,
     &    k56h)*yp1 + za(k2,k5)*za(k3,k34h)*zb(k1,k34h)*zb(k2,k4)*zb(k6
     &    ,k56h)*iza(k1,k2)**2*izb(k1,k56h)*yp1 + za(k2,k5)*za(k3,k56h)
     &    *zb(k1,k2)*zb(k2,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k56h)*izb(
     &    k1,k34h) )
      triamp = triamp + mtsq * ( za(k2,k5)*za(k3,k56h)*zb(k1,k2)*zb(k4,
     &    k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k34h)*xpr + za(k2,k5)*
     &    za(k3,k56h)*zb(k1,k56h)*zb(k2,k4)*zb(k6,k34h)*iza(k1,k2)**2*
     &    izb(k1,k34h)*xpr + za(k3,k34h)*za(k5,k34h)*zb(k1,k4)*zb(k2,
     &    k34h)*zb(k6,k56h)*iza(k1,k2)*iza(k1,k34h)*izb(k1,k56h)*ymr + 
     &    za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k2,k4)*zb(k6,k56h)*
     &    iza(k1,k2)*iza(k1,k34h)*izb(k1,k56h)*ymr + za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k4)*zb(k2,k56h)*zb(k6,k56h)*iza(k1,k2)*iza(k1,
     &    k34h)*izb(k1,k56h)*xm1 + za(k3,k34h)*za(k5,k56h)*zb(k2,k4)*
     &    zb(k6,k56h)*iza(k1,k2)*iza(k1,k34h)*xm1 + za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k4)*zb(k2,k34h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k56h)*izb(k1,k34h)*ymr + za(k3,k56h)*za(k5,k34h)*zb(k2,k4)*
     &    zb(k6,k34h)*iza(k1,k2)*iza(k1,k56h)*ymr + za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k4)*zb(k2,k56h)*zb(k6,k34h)*iza(k1,k2)*iza(k1,
     &    k56h)*izb(k1,k34h)*xm1 + za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*
     &    zb(k2,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k56h)*izb(k1,k34h)*
     &    xm1 )
      triamp = triamp + mtsq*t234sum * (  - za(k2,k5)*za(k3,k34h)*zb(k1
     &    ,k2)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(
     &    k1,k56h)**2 - za(k2,k5)*za(k3,k34h)*zb(k2,k4)*zb(k6,k56h)*
     &    iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h) - za(k2,k5)*za(k3,
     &    k56h)*zb(k1,k2)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,
     &    k56h)*izb(k1,k34h)**2 - za(k2,k5)*za(k3,k56h)*zb(k2,k4)*zb(k6
     &    ,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,k34h) - za(k3,k34h)*
     &    za(k5,k34h)*zb(k2,k4)*zb(k6,k56h)*iza(k1,k2)*iza(k1,k34h)**2*
     &    izb(k1,k56h) - za(k3,k56h)*za(k5,k56h)*zb(k2,k4)*zb(k6,k34h)*
     &    iza(k1,k2)*iza(k1,k56h)**2*izb(k1,k34h) )
      triamp = triamp + mtsq*t134sum * ( za(k1,k5)*za(k2,k34h)*za(k3,
     &    k34h)*zb(k1,k2)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,
     &    k34h)**2*izb(k1,k56h)**2 + za(k1,k5)*za(k2,k56h)*za(k3,k56h)*
     &    zb(k1,k2)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)
     &    **2*izb(k1,k34h)**2 + za(k1,k5)*za(k3,k34h)*zb(k1,k4)*zb(k6,
     &    k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h) + za(k1,k5)*za(
     &    k3,k56h)*zb(k1,k4)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*
     &    izb(k1,k34h) + za(k3,k34h)*za(k5,k34h)*zb(k1,k4)*zb(k2,k56h)*
     &    zb(k6,k56h)*iza(k1,k2)*iza(k1,k34h)**2*izb(k1,k56h)**2 + za(
     &    k3,k56h)*za(k5,k56h)*zb(k1,k4)*zb(k2,k34h)*zb(k6,k34h)*iza(k1
     &    ,k2)*iza(k1,k56h)**2*izb(k1,k34h)**2 )
      triamp = triamp + y * ( za(k1,k34h)*za(k2,k34h)*za(k3,k34h)*za(k5
     &    ,k56h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k2,k56h)*yp1**2 + za(k1,k34h)*za(k2,k34h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k1,k56h)*yp1**2 + za(k1,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*xpr*ymr + za(k1,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*yp1*xpr + za(k1,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*xm1*xpr + za(k1,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k1,k56h)*yp1*xm1 + za(k1,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(
     &    k2,k56h)*yp1*ymr + za(k1,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*yp1*ymr )
      triamp = triamp + y * (  - za(k1,k34h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*yp1**2 + 
     &    za(k1,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k56h)*zb(k2,k34h)*
     &    zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k56h)*yp1**2 + 
     &    za(k1,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k2,k56h)*
     &    zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k34h)*yp1*xpr - 
     &    za(k1,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*
     &    zb(k6,k34h)*iza(k1,k2)**2*yp1*xpr - za(k1,k34h)*za(k3,k34h)*
     &    za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2
     &    *yp1*xm1 + za(k1,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*
     &    zb(k2,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k56h
     &    )*yp1*xm1 + za(k1,k56h)*za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h
     &    )*yp1*ymr + za(k1,k56h)*za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k2,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k34h
     &    )*xpr*ymr )
      triamp = triamp + y * ( za(k1,k56h)*za(k2,k34h)*za(k3,k34h)*za(k5
     &    ,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k2,k56h)*yp1*xpr + za(k1,k56h)*za(k2,k34h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k2,k56h)*yp1*xm1 + za(k1,k56h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*xpr**2 + za(k1,k56h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k1,k34h)*xpr**2 + za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(
     &    k2,k34h)*xpr*ymr + za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k56h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(
     &    k2,k56h)*yp1*xpr - za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*xpr*ymr - za(k1,
     &    k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*yp1*xpr )
      triamp = triamp + y * ( za(k1,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1
     &    ,k34h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(
     &    k2,k34h)*xpr**2 + za(k1,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,
     &    k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(
     &    k2,k34h)*xm1*xpr - za(k1,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,
     &    k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*xpr**2 + za(k1,
     &    k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*xm1*xpr + za(k2,k34h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k2)**2*izb(k1,k56h)*yp1*ymr + 2._dp*za(k2,k34h)*za(k3,k34h)*
     &    za(k5,k34h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2
     &    *yp1*ymr + za(k2,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*
     &    zb(k2,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k34h
     &    )*yp1*xpr + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*
     &    zb(k2,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k56h
     &    )*yp1*xpr )
      triamp = triamp + y * ( za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1
     &    ,k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(
     &    k1,k34h)*xm1*xpr + 2._dp*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*
     &    zb(k2,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*xm1*xpr )
      triamp = triamp + y*t234sum * (  - za(k1,k34h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*yp1
     &     - za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k34h
     &    )*iza(k1,k2)**2*iza(k1,k34h)*xpr - za(k2,k34h)*za(k3,k34h)*
     &    za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h
     &    )*yp1 - za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k4,k56h)*zb(k6
     &    ,k34h)*iza(k1,k2)**2*iza(k2,k34h)*xpr - za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(
     &    k1,k56h)*ymr - za(k3,k34h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,k56h
     &    )*zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k56h)*yp1 - za(k3,k34h)*
     &    za(k5,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*ymr - za(k3
     &    ,k34h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*yp1
     &     - za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h
     &    )*iza(k1,k2)**2*izb(k1,k34h)*xm1 - za(k3,k34h)*za(k5,k56h)*
     &    zb(k2,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k34h
     &    )*xpr )
      triamp = triamp + y*t234sum * (  - za(k3,k34h)*za(k5,k56h)*zb(k4,
     &    k56h)*zb(k6,k34h)*iza(k1,k2)**2*xpr - za(k3,k34h)*za(k5,k56h)
     &    *zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*xm1 )
      triamp = triamp + y*t234sum**2 * ( za(k3,k34h)*za(k5,k34h)*zb(k4,
     &    k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h) + 
     &    za(k3,k34h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2
     &    *iza(k1,k56h)*izb(k1,k34h) )
      triamp = triamp + y*t134sum * (  - za(k1,k34h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*yp1
     &     - za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k34h
     &    )*iza(k1,k2)**2*iza(k1,k34h)*xpr + za(k2,k34h)*za(k3,k34h)*
     &    za(k5,k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2
     &    *iza(k1,k34h)*izb(k1,k56h)*ymr + za(k2,k34h)*za(k3,k34h)*za(
     &    k5,k56h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*
     &    iza(k1,k56h)*izb(k1,k34h)*yp1 + za(k2,k56h)*za(k3,k34h)*za(k5
     &    ,k34h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k1,k34h)*izb(k1,k56h)*xpr + za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k1,k56h)*izb(k1,k34h)*xm1 - za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k56h)*yp1
     &     - za(k3,k34h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)
     &    **2*ymr - za(k3,k34h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k34h)*
     &    iza(k1,k2)**2*yp1 )
      triamp = triamp + y*t134sum * (  - za(k3,k34h)*za(k5,k56h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k34h)*xpr
     &     - za(k3,k34h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)
     &    **2*xpr - za(k3,k34h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k34h)*
     &    iza(k1,k2)**2*xm1 )
      triamp = triamp + y*t134sum**2 * (  - za(k3,k34h)*za(k5,k34h)*zb(
     &    k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)
     &     - za(k3,k34h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)
     &    **2*iza(k1,k56h)*izb(k1,k34h) )
      triamp = triamp + y*mtsq * ( za(k1,k5)*za(k2,k34h)*za(k3,k34h)*
     &    zb(k1,k2)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k34h)*
     &    izb(k1,k56h) + za(k1,k5)*za(k2,k56h)*za(k3,k34h)*zb(k1,k2)*
     &    zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,
     &    k34h) + 2._dp*za(k1,k5)*za(k3,k34h)*zb(k1,k4)*zb(k6,k34h)*iza(
     &    k1,k2)**2 + za(k2,k5)*za(k3,k34h)*zb(k1,k2)*zb(k4,k34h)*zb(k6
     &    ,k34h)*iza(k1,k2)**2*izb(k1,k34h) + za(k2,k5)*za(k3,k34h)*zb(
     &    k1,k2)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k56h) + 2.
     &    _dp*za(k2,k5)*za(k3,k34h)*zb(k2,k4)*zb(k6,k34h)*iza(k1,k2)**2
     &     + za(k3,k34h)*za(k5,k34h)*zb(k1,k4)*zb(k2,k56h)*zb(k6,k34h)*
     &    iza(k1,k2)*iza(k1,k34h)*izb(k1,k56h) + za(k3,k34h)*za(k5,k34h
     &    )*zb(k2,k4)*zb(k6,k34h)*iza(k1,k2)*iza(k1,k34h) + za(k3,k34h)
     &    *za(k5,k56h)*zb(k1,k4)*zb(k2,k34h)*zb(k6,k34h)*iza(k1,k2)*
     &    iza(k1,k56h)*izb(k1,k34h) + za(k3,k34h)*za(k5,k56h)*zb(k2,k4)
     &    *zb(k6,k34h)*iza(k1,k2)*iza(k1,k56h) )
      triamp = triamp + x * ( za(k1,k34h)*za(k2,k34h)*za(k3,k56h)*za(k5
     &    ,k56h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k56h)*yp1**2 + za(k1,k34h)*za(k2,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k1,k56h)*yp1**2 + za(k1,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*xpr*ymr + za(k1,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*yp1*xpr + za(k1,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*xm1*xpr + za(k1,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k1,k56h)*yp1*xm1 + za(k1,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(
     &    k2,k56h)*yp1*ymr + za(k1,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*yp1*ymr )
      triamp = triamp + x * (  - za(k1,k34h)*za(k3,k56h)*za(k5,k34h)*
     &    zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*yp1**2 + 
     &    za(k1,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k2,k34h)*
     &    zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k2,k56h)*yp1**2 + 
     &    za(k1,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k2,k56h)*
     &    zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k2,k34h)*yp1*xpr - 
     &    za(k1,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*
     &    zb(k6,k56h)*iza(k1,k2)**2*yp1*xpr - za(k1,k34h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2
     &    *yp1*xm1 + za(k1,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*
     &    zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k2,k56h
     &    )*yp1*xm1 + za(k1,k56h)*za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*
     &    zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k2,k56h
     &    )*yp1*ymr + za(k1,k56h)*za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*
     &    zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h
     &    )*xpr*ymr )
      triamp = triamp + x * ( za(k1,k56h)*za(k2,k34h)*za(k3,k56h)*za(k5
     &    ,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k56h)*yp1*xpr + za(k1,k56h)*za(k2,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k56h)*yp1*xm1 + za(k1,k56h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*xpr**2 + za(k1,k56h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k1,k34h)*xpr**2 + za(k1,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(
     &    k2,k34h)*xpr*ymr + za(k1,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,
     &    k56h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(
     &    k2,k56h)*yp1*xpr - za(k1,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*xpr*ymr - za(k1,
     &    k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*yp1*xpr )
      triamp = triamp + x * ( za(k1,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1
     &    ,k34h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(
     &    k2,k34h)*xpr**2 + za(k1,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,
     &    k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(
     &    k2,k34h)*xm1*xpr - za(k1,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,
     &    k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*xpr**2 + za(k1,
     &    k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*xm1*xpr + za(k2,k34h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k2)**2*izb(k1,k56h)*yp1*ymr + 2._dp*za(k2,k34h)*za(k3,k56h)*
     &    za(k5,k34h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2
     &    *yp1*ymr + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*
     &    zb(k2,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k34h
     &    )*yp1*xpr + za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*
     &    zb(k2,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k56h
     &    )*yp1*xpr )
      triamp = triamp + x * ( za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1
     &    ,k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(
     &    k1,k34h)*xm1*xpr + 2._dp*za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*
     &    zb(k2,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*xm1*xpr )
      triamp = triamp + x*t234sum * (  - za(k1,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k56h)*yp1
     &     - za(k1,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k56h
     &    )*iza(k1,k2)**2*iza(k1,k34h)*xpr - za(k2,k34h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k2,k56h
     &    )*yp1 - za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k4,k56h)*zb(k6
     &    ,k56h)*iza(k1,k2)**2*iza(k2,k34h)*xpr - za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(
     &    k1,k56h)*ymr - za(k3,k56h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,k56h
     &    )*zb(k6,k56h)*iza(k1,k2)**2*izb(k2,k56h)*yp1 - za(k3,k56h)*
     &    za(k5,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*ymr - za(k3
     &    ,k56h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*yp1
     &     - za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h
     &    )*iza(k1,k2)**2*izb(k1,k34h)*xm1 - za(k3,k56h)*za(k5,k56h)*
     &    zb(k2,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k2,k34h
     &    )*xpr )
      triamp = triamp + x*t234sum * (  - za(k3,k56h)*za(k5,k56h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*xpr - za(k3,k56h)*za(k5,k56h)
     &    *zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*xm1 )
      triamp = triamp + x*t234sum**2 * ( za(k3,k56h)*za(k5,k34h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h) + 
     &    za(k3,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2
     &    *iza(k1,k56h)*izb(k1,k34h) )
      triamp = triamp + x*t134sum * (  - za(k1,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k56h)*yp1
     &     - za(k1,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k56h
     &    )*iza(k1,k2)**2*iza(k1,k34h)*xpr + za(k2,k34h)*za(k3,k56h)*
     &    za(k5,k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2
     &    *iza(k1,k34h)*izb(k1,k56h)*ymr + za(k2,k34h)*za(k3,k56h)*za(
     &    k5,k56h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*
     &    iza(k1,k56h)*izb(k1,k34h)*yp1 + za(k2,k56h)*za(k3,k56h)*za(k5
     &    ,k34h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k1,k34h)*izb(k1,k56h)*xpr + za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k1,k56h)*izb(k1,k34h)*xm1 - za(k3,k56h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k56h)*yp1
     &     - za(k3,k56h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)
     &    **2*ymr - za(k3,k56h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k56h)*
     &    iza(k1,k2)**2*yp1 )
      triamp = triamp + x*t134sum * (  - za(k3,k56h)*za(k5,k56h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k34h)*xpr
     &     - za(k3,k56h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)
     &    **2*xpr - za(k3,k56h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k56h)*
     &    iza(k1,k2)**2*xm1 )
      triamp = triamp + x*t134sum**2 * (  - za(k3,k56h)*za(k5,k34h)*zb(
     &    k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h)
     &     - za(k3,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)
     &    **2*iza(k1,k56h)*izb(k1,k34h) )
      triamp = triamp + x*mtsq * ( za(k1,k5)*za(k2,k34h)*za(k3,k56h)*
     &    zb(k1,k2)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*
     &    izb(k1,k56h) + za(k1,k5)*za(k2,k56h)*za(k3,k56h)*zb(k1,k2)*
     &    zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,
     &    k34h) + 2._dp*za(k1,k5)*za(k3,k56h)*zb(k1,k4)*zb(k6,k56h)*iza(
     &    k1,k2)**2 + za(k2,k5)*za(k3,k56h)*zb(k1,k2)*zb(k4,k34h)*zb(k6
     &    ,k56h)*iza(k1,k2)**2*izb(k1,k34h) + za(k2,k5)*za(k3,k56h)*zb(
     &    k1,k2)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k56h) + 2.
     &    _dp*za(k2,k5)*za(k3,k56h)*zb(k2,k4)*zb(k6,k56h)*iza(k1,k2)**2
     &     + za(k3,k56h)*za(k5,k34h)*zb(k1,k4)*zb(k2,k56h)*zb(k6,k56h)*
     &    iza(k1,k2)*iza(k1,k34h)*izb(k1,k56h) + za(k3,k56h)*za(k5,k34h
     &    )*zb(k2,k4)*zb(k6,k56h)*iza(k1,k2)*iza(k1,k34h) + za(k3,k56h)
     &    *za(k5,k56h)*zb(k1,k4)*zb(k2,k34h)*zb(k6,k56h)*iza(k1,k2)*
     &    iza(k1,k56h)*izb(k1,k34h) + za(k3,k56h)*za(k5,k56h)*zb(k2,k4)
     &    *zb(k6,k56h)*iza(k1,k2)*iza(k1,k56h) )
      triamp = triamp + x*y * ( za(k1,k34h)**2*za(k3,k56h)*za(k5,k56h)*
     &    zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h
     &    )*yp1 + za(k1,k34h)*za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1
     &    ,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h)*yp1
     &     - za(k1,k34h)*za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h
     &    )*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h)*yp1 + 
     &    za(k1,k34h)*za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k2,k34h)*
     &    zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*yp1 + za(
     &    k1,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(
     &    k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k2,k34h)*xpr + za(k1,
     &    k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,
     &    k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k34h)*xpr + za(k1,k34h
     &    )*za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k2,k56h)*zb(k4,k34h)
     &    *zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*yp1 + za(k1,k34h)*za(
     &    k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(
     &    k6,k56h)*iza(k1,k2)**2*izb(k2,k56h)*yp1 )
      triamp = triamp + x*y * ( za(k1,k34h)*za(k3,k34h)*za(k5,k56h)*zb(
     &    k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*yp1 + 2._dp*za(
     &    k1,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(
     &    k6,k34h)*iza(k1,k2)**2*ymr - za(k1,k34h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*yp1
     &     + za(k1,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k2,k34h
     &    )*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k56h)*yp1 + 
     &    za(k1,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k2,k56h)*
     &    zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k34h)*xpr + za(
     &    k1,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*
     &    zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k34h)*xpr - za(k1,k34h)*za(
     &    k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(
     &    k6,k34h)*iza(k1,k2)**2*izb(k2,k34h)*xpr + za(k1,k56h)**2*za(
     &    k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(
     &    k1,k2)**2*iza(k1,k34h)*xpr + za(k1,k56h)*za(k2,k34h)*za(k3,
     &    k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,
     &    k2)**2*iza(k2,k56h)*yp1 )
      triamp = triamp + x*y * ( za(k1,k56h)*za(k2,k34h)*za(k3,k34h)*za(
     &    k5,k56h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*
     &    iza(k1,k34h)*xpr + za(k1,k56h)*za(k2,k34h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k2,k56h)*yp1 - za(k1,k56h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h
     &    )*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k2,
     &    k34h)*xpr + za(k1,k56h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k2,k34h
     &    )*xpr + za(k1,k56h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k2
     &    ,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*xpr
     &     + za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)**2*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k56h)*yp1 - za(k1,k56h
     &    )*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k2,k34h)*zb(k4,k56h)
     &    *zb(k6,k56h)*iza(k1,k2)**2*izb(k2,k56h)*yp1 + za(k1,k56h)*za(
     &    k3,k34h)*za(k5,k34h)*zb(k1,k56h)*zb(k2,k34h)*zb(k4,k34h)*zb(
     &    k6,k56h)*iza(k1,k2)**2*izb(k2,k56h)*yp1 )
      triamp = triamp + x*y * ( za(k1,k56h)*za(k3,k34h)*za(k5,k56h)*zb(
     &    k1,k34h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*
     &    izb(k2,k34h)*xpr - za(k1,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*xpr + 2._dp*za(k1,
     &    k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*xm1 + za(k1,k56h)*za(k3,k56h)*za(k5,k34h)
     &    *zb(k1,k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**
     &    2*izb(k2,k34h)*xpr + za(k1,k56h)*za(k3,k56h)*za(k5,k34h)*zb(
     &    k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*xpr + za(k2,
     &    k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k2,k34h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k56h)*yp1 + 2._dp*za(k2
     &    ,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*yp1 + 2._dp*za(k2,k34h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*ymr
     &     + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k2,k34h
     &    )*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k34h)*xpr )
      triamp = triamp + x*y * ( za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(
     &    k1,k34h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*
     &    izb(k1,k56h)*yp1 + 2._dp*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*
     &    zb(k2,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*xm1 + za(k2
     &    ,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k2,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k34h)*xpr + 2._dp*za(k2
     &    ,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*xpr )
      triamp = triamp + x*y*t234sum * ( za(k1,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k1,k56h)*izb(k1,k34h) - za(k1,k34h)*za(k3,k56h)*za(k5,k56h)*
     &    zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h) + za(k1,
     &    k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,
     &    k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h) - za(k1,k56h)*
     &    za(k3,k34h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2
     &    *iza(k1,k34h) - za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k56h) - za(k3,k34h)*
     &    za(k5,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2 - za(k3,
     &    k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k2)**2*izb(k1,k34h) - za(k3,k56h)*za(k5,k34h)*zb(k4,k56h)*zb(
     &    k6,k34h)*iza(k1,k2)**2 )
      triamp = triamp + x*y*t134sum * ( za(k1,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(
     &    k1,k56h)*izb(k1,k34h) + za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h
     &    )*izb(k1,k56h) + za(k2,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k2,
     &    k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(
     &    k1,k56h) + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k2,k34h)*
     &    zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,
     &    k34h) + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k56h)*zb(k4
     &    ,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k1,k34h)*izb(k1,k56h) + 
     &    za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k2,k56h)*zb(k4,k34h)*
     &    zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*izb(k1,k34h) - za(k3,
     &    k34h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2 - za(
     &    k3,k56h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2 )
      triamp = triamp + x*y**2 * ( 2._dp*za(k1,k34h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2 + 2._dp
     &    *za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*
     &    zb(k6,k34h)*iza(k1,k2)**2 + 2._dp*za(k2,k34h)*za(k3,k34h)*za(
     &    k5,k56h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2 + 
     &    2._dp*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k2,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)**2 )
      triamp = triamp + x**2*y * ( 2._dp*za(k1,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2 + 2._dp
     &    *za(k1,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*
     &    zb(k6,k56h)*iza(k1,k2)**2 + 2._dp*za(k2,k34h)*za(k3,k56h)*za(
     &    k5,k56h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2 + 
     &    2._dp*za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k2,k56h)*zb(k4,
     &    k34h)*zb(k6,k56h)*iza(k1,k2)**2 )
      triamp = triamp - za(k1,k34h)**2*za(k3,k56h)*za(k5,k56h)*zb(k1,
     & k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*yp1**2*
     & xm1 + za(k1,k34h)*za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h
     &    )*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h)*yp1**2*
     &    ymr + za(k1,k34h)*za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k2,
     &    k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k1,k56h)*
     &    yp1**2*ymr + za(k1,k34h)*za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*
     &    zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h
     &    )*yp1**2*xm1 + za(k1,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h
     &    )*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(k2,
     &    k34h)*yp1*xpr*ymr + za(k1,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5
     &    ,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*yp1*xm1*xpr + za(k1,k34h)*za(k2,k56h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2
     &    *iza(k1,k56h)*yp1*xm1*xpr + za(k1,k34h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,
     &    k2)**2*izb(k2,k56h)*yp1**2*ymr
      triamp = triamp - za(k1,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)
     & *zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*yp1**2*xm1 + za(k1,k34h)*
     &    za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k2,k34h)*zb(k4,k34h)*
     &    zb(k6,k56h)*iza(k1,k2)**2*izb(k2,k56h)*yp1**2*xm1 + za(k1,
     &    k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k2,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k34h)*yp1*xpr*ymr - 
     &    za(k1,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*
     &    zb(k6,k34h)*iza(k1,k2)**2*yp1*xpr*ymr - za(k1,k34h)*za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(
     &    k1,k2)**2*izb(k1,k34h)*yp1*xm1*xpr + za(k1,k34h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k1,k56h)*zb(k2,k56h)*zb(k4,k34h)*zb(k6,k34h)*
     &    iza(k1,k2)**2*izb(k2,k34h)*yp1*xm1*xpr - za(k1,k56h)**2*za(k3
     &    ,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1
     &    ,k2)**2*iza(k1,k34h)*xpr**2*ymr + za(k1,k56h)*za(k2,k34h)*za(
     &    k3,k34h)*za(k5,k34h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(
     &    k1,k2)**2*iza(k1,k34h)*yp1*xpr*ymr
      triamp = triamp + za(k1,k56h)*za(k2,k34h)*za(k3,k56h)*za(k5,k34h)
     & *zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h)*
     & yp1*xpr*ymr + za(k1,k56h)*za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*
     &    zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k2)**2*iza(k2,k56h
     &    )*yp1*xm1*xpr + za(k1,k56h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*iza(
     &    k2,k34h)*xpr**2*ymr + za(k1,k56h)*za(k2,k56h)*za(k3,k34h)*za(
     &    k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*
     &    iza(k2,k34h)*xm1*xpr**2 + za(k1,k56h)*za(k2,k56h)*za(k3,k34h)
     &    *za(k5,k56h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**
     &    2*iza(k1,k34h)*xm1*xpr**2 - za(k1,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)**2*
     &    izb(k1,k56h)*yp1*xpr*ymr + za(k1,k56h)*za(k3,k34h)*za(k5,k34h
     &    )*zb(k1,k34h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k2)
     &    **2*izb(k2,k56h)*yp1*xpr*ymr + za(k1,k56h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k2,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k2)**2*izb(k2,k56h)*yp1*xm1*xpr
      triamp = triamp - za(k1,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)
     & *zb(k4,k34h)*zb(k6,k56h)*iza(k1,k2)**2*yp1*xm1*xpr + za(k1,k56h)
     &    *za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k2,k56h)*zb(k4,k56h)*
     &    zb(k6,k34h)*iza(k1,k2)**2*izb(k2,k34h)*xpr**2*ymr - za(k1,
     &    k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*xpr**2*ymr + za(k1,k56h)*za(k3,k56h)*za(
     &    k5,k56h)*zb(k1,k56h)*zb(k2,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(
     &    k1,k2)**2*izb(k2,k34h)*xm1*xpr**2 + za(k2,k34h)*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,k56h)*
     &    iza(k1,k2)**2*izb(k1,k56h)*yp1**2*ymr + za(k2,k34h)*za(k3,
     &    k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k2,k34h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k2)**2*izb(k1,k34h)*yp1*xpr*ymr + za(k2,k56h)*
     &    za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k2,k56h)*zb(k4,k56h)*
     &    zb(k6,k56h)*iza(k1,k2)**2*izb(k1,k56h)*yp1*xm1*xpr + za(k2,
     &    k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k2,k56h)*zb(k4,
     &    k56h)*zb(k6,k34h)*iza(k1,k2)**2*izb(k1,k34h)*xm1*xpr**2

      app=triamp*(-half*im)/(s(k3,k4)*s(k5,k6))

      if (Higgsint) return

      triamp =  + s12**(-1) * ( za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k34h)**2*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)*izb(k1,
     &    k56h)*yp1**2*ymr + za(k2,k34h)**2*za(k3,k34h)*za(k5,k56h)*zb(
     &    k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)*yp1**2*xm1 + 
     &    za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)**
     &    2*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*izb(k1,k56h)*yp1*xpr*
     &    ymr + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)*yp1*xpr*ymr + za(
     &    k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(
     &    k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*yp1*xm1*xpr + za(k2,k34h)*
     &    za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*
     &    zb(k6,k56h)*iza(k1,k34h)*yp1*xm1*xpr + za(k2,k34h)*za(k2,k56h
     &    )*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)
     &    *iza(k1,k56h)*yp1*xpr*ymr + za(k2,k34h)*za(k2,k56h)*za(k3,
     &    k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k56h)*yp1*xpr*ymr )
      triamp = triamp + s12**(-1) * ( za(k2,k34h)*za(k2,k56h)*za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(
     &    k1,k56h)*izb(k1,k34h)*yp1*xm1*xpr + za(k2,k34h)*za(k2,k56h)*
     &    za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*
     &    iza(k1,k56h)*yp1*xm1*xpr + za(k2,k34h)*za(k3,k34h)*za(k5,k34h
     &    )*zb(k1,k34h)**2*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k56h)*yp1**2*
     &    ymr + za(k2,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k56h)*yp1**2*xm1 + za(k2
     &    ,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(
     &    k6,k34h)*izb(k2,k34h)*yp1*xpr*ymr + za(k2,k34h)*za(k3,k56h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*
     &    izb(k2,k34h)*yp1*xpr*ymr + za(k2,k34h)*za(k3,k56h)*za(k5,k56h
     &    )*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k34h
     &    )*yp1*xm1*xpr + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,
     &    k56h)**2*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k34h)*yp1*xm1*xpr + 
     &    za(k2,k56h)**2*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k56h
     &    )*zb(k6,k34h)*iza(k1,k56h)*xpr**2*ymr )
      triamp = triamp + s12**(-1) * ( za(k2,k56h)**2*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k56h)**2*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)*
     &    izb(k1,k34h)*xm1*xpr**2 + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)
     &    *zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)*yp1*xpr*
     &    ymr + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k56h)*yp1*xpr*ymr + za(
     &    k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(
     &    k4,k56h)*zb(k6,k56h)*izb(k2,k56h)*yp1*xm1*xpr + za(k2,k56h)*
     &    za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k56h
     &    )*izb(k2,k56h)*yp1*xm1*xpr + za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(k2,
     &    k34h)*xpr**2*ymr + za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,
     &    k56h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k34h)*xm1*xpr**2 )
      triamp = triamp + s12**(-1)*t234sum * ( za(k2,k34h)**2*za(k3,k56h
     &    )*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)**2*yp1**2
     &     - 2._dp*za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4
     &    ,k34h)*zb(k6,k34h)*iza(k1,k56h)*izb(k2,k34h)*yp1*ymr - za(k2,
     &    k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,
     &    k34h)*iza(k1,k56h)*izb(k2,k34h)*yp1*xpr - za(k2,k34h)*za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k56h)*izb(k2,k34h)*yp1*xpr - 2._dp*za(k2,k34h)*za(k3,k56h)*za(
     &    k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*
     &    izb(k2,k34h)*yp1*xm1 + za(k2,k56h)**2*za(k3,k34h)*za(k5,k34h)
     &    *zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*xpr**2 - 2._dp*za(k2,
     &    k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,
     &    k56h)*iza(k1,k34h)*izb(k2,k56h)*xpr*ymr - za(k2,k56h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k34h)*izb(k2,k56h)*yp1*xpr - za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)*izb(k2
     &    ,k56h)*yp1*xpr )
      triamp = triamp + s12**(-1)*t234sum * (  - 2._dp*za(k2,k56h)*za(k3
     &    ,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1
     &    ,k34h)*izb(k2,k56h)*xm1*xpr + za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)**2*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)**2*yp1*ymr + 
     &    za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*
     &    zb(k6,k56h)*izb(k2,k56h)**2*yp1*ymr + za(k3,k34h)*za(k5,k34h)
     &    *zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k56h)
     &    **2*yp1**2 + za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*
     &    zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)**2*yp1*xm1 + za(k3,k34h)
     &    *za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k56h)*izb(k2,
     &    k56h)**2*yp1*xm1 + za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)**2*zb(
     &    k4,k56h)*zb(k6,k34h)*izb(k2,k34h)**2*xpr*ymr + za(k3,k56h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*
     &    izb(k2,k34h)**2*xpr*ymr + za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)
     &    *zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k34h)**2*xpr**2
     &     + za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h
     &    )*zb(k6,k34h)*izb(k2,k34h)**2*xm1*xpr )
      triamp = triamp + s12**(-1)*t234sum * ( za(k3,k56h)*za(k5,k56h)*
     &    zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k34h)**2*xm1*
     &    xpr )
      triamp = triamp + s12**(-1)*t234sum**2 * ( 2._dp*za(k2,k34h)*za(k3
     &    ,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)**2*
     &    izb(k2,k34h)*yp1 + 2._dp*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*izb(k2,k56h)*xpr - 
     &    za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*
     &    iza(k1,k34h)*izb(k2,k56h)**2*ymr - za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*izb(k2,k56h)
     &    **2*yp1 - za(k3,k34h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(
     &    k6,k56h)*iza(k1,k34h)*izb(k2,k56h)**2*yp1 - za(k3,k34h)*za(k5
     &    ,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*izb(
     &    k2,k56h)**2*xm1 - za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k56h)*izb(k2,k34h)**2*ymr - za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k56h)*izb(k2,k34h)**2*xpr - za(k3,k56h)*za(k5,k56h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*izb(k2,k34h)**2*
     &    xpr )
      triamp = triamp + s12**(-1)*t234sum**2 * (  - za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*izb(k2
     &    ,k34h)**2*xm1 )
      triamp = triamp + s12**(-1)*t234sum**3 * ( za(k3,k34h)*za(k5,k34h
     &    )*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*izb(k2,k56h)**2 + 
     &    za(k3,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)
     &    **2*izb(k2,k34h)**2 )
      triamp = triamp + s12**(-1)*t134sum * ( za(k2,k34h)**2*za(k3,k34h
     &    )*za(k5,k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k34h)**2*izb(k1,k56h)**2*yp1*ymr + za(k2,k34h)**2*za(k3,k34h)
     &    *za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)
     &    **2*izb(k1,k56h)*yp1*ymr + za(k2,k34h)**2*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)**2*
     &    izb(k1,k56h)*yp1**2 + za(k2,k34h)**2*za(k3,k34h)*za(k5,k56h)*
     &    zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*izb(k1,
     &    k56h)*yp1*xm1 + za(k2,k34h)**2*za(k3,k34h)*za(k5,k56h)*zb(k4,
     &    k34h)*zb(k6,k56h)*iza(k1,k34h)**2*yp1*xm1 + za(k2,k34h)**2*
     &    za(k3,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)
     &    **2*yp1**2 + 2._dp*za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*
     &    izb(k1,k56h)*xpr*ymr + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)
     &    **2*izb(k1,k56h)*yp1*xpr )
      triamp = triamp + s12**(-1)*t134sum * ( za(k2,k34h)*za(k2,k56h)*
     &    za(k3,k34h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)
     &    **2*yp1*xpr + 2._dp*za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*xm1*xpr + 2._dp*
     &    za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k4,k34h)*
     &    zb(k6,k34h)*iza(k1,k56h)**2*yp1*ymr + za(k2,k34h)*za(k2,k56h)
     &    *za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*
     &    iza(k1,k56h)**2*izb(k1,k34h)*yp1*xpr + 2._dp*za(k2,k34h)*za(k2
     &    ,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k56h)**2*izb(k1,k34h)*yp1*xm1 + za(k2,k34h)*za(
     &    k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(
     &    k1,k56h)**2*yp1*xpr + za(k2,k56h)**2*za(k3,k34h)*za(k5,k34h)*
     &    zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*xpr**2 + za(k2,k56h)
     &    **2*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k56h)**2*izb(k1,k34h)*xpr*ymr + za(k2,k56h)**2*
     &    za(k3,k56h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)
     &    **2*xpr*ymr )
      triamp = triamp + s12**(-1)*t134sum * ( za(k2,k56h)**2*za(k3,k56h
     &    )*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k56h)**2*izb(k1,k34h)**2*xm1*xpr + za(k2,k56h)**2*za(k3,k56h)
     &    *za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)
     &    **2*izb(k1,k34h)*xpr**2 + za(k2,k56h)**2*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)**2*
     &    izb(k1,k34h)*xm1*xpr )
      triamp = triamp + s12**(-1)*t134sum**2 * ( za(k2,k34h)**2*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k34h)**3*izb(k1,k56h)**2*ymr + za(k2,k34h)**2*za(k3,k34h)*za(
     &    k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**3*
     &    izb(k1,k56h)**2*yp1 + za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*
     &    zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)**3*izb(k1,k56h)*yp1 + 
     &    za(k2,k34h)**2*za(k3,k34h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k56h
     &    )*iza(k1,k34h)**3*izb(k1,k56h)*xm1 + 2._dp*za(k2,k34h)*za(k2,
     &    k56h)*za(k3,k34h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k34h)**3*izb(k1,k56h)*xpr + 2._dp*za(k2,k34h)*za(k2,k56h)*za(
     &    k3,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)**3*
     &    izb(k1,k34h)*yp1 + za(k2,k56h)**2*za(k3,k56h)*za(k5,k34h)*zb(
     &    k4,k34h)*zb(k6,k34h)*iza(k1,k56h)**3*izb(k1,k34h)*ymr + za(k2
     &    ,k56h)**2*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(
     &    k6,k34h)*iza(k1,k56h)**3*izb(k1,k34h)**2*xpr + za(k2,k56h)**2
     &    *za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*
     &    iza(k1,k56h)**3*izb(k1,k34h)**2*xm1 )
      triamp = triamp + s12**(-1)*t134sum**2 * ( za(k2,k56h)**2*za(k3,
     &    k56h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)**3*
     &    izb(k1,k34h)*xpr )
      triamp = triamp + s12**(-1)*t134sum**3 * ( za(k2,k34h)**2*za(k3,
     &    k34h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**4*
     &    izb(k1,k56h)**2 + za(k2,k56h)**2*za(k3,k56h)*za(k5,k56h)*zb(
     &    k4,k34h)*zb(k6,k34h)*iza(k1,k56h)**4*izb(k1,k34h)**2 )
      triamp = triamp + mtsq*s12**(-1) * (  - za(k2,k5)*za(k2,k34h)*za(
     &    k3,k34h)*zb(k1,k4)*zb(k1,k34h)*zb(k6,k56h)*iza(k1,k34h)*izb(
     &    k1,k56h)*yp1 - za(k2,k5)*za(k2,k56h)*za(k3,k56h)*zb(k1,k4)*
     &    zb(k1,k56h)*zb(k6,k34h)*iza(k1,k56h)*izb(k1,k34h)*xpr - za(k2
     &    ,k5)*za(k3,k34h)*zb(k1,k4)*zb(k1,k34h)*zb(k6,k56h)*izb(k2,
     &    k56h)*yp1 - za(k2,k5)*za(k3,k56h)*zb(k1,k4)*zb(k1,k56h)*zb(k6
     &    ,k34h)*izb(k2,k34h)*xpr )
      triamp = triamp + mtsq*s12**(-1)*t234sum * ( za(k2,k5)*za(k3,k34h
     &    )*zb(k1,k4)*zb(k6,k56h)*iza(k1,k34h)*izb(k2,k56h) + za(k2,k5)
     &    *za(k3,k56h)*zb(k1,k4)*zb(k6,k34h)*iza(k1,k56h)*izb(k2,k34h)
     &     )
      triamp = triamp + mtsq*s12**(-1)*t134sum * (  - za(k2,k5)*za(k2,
     &    k34h)*za(k3,k34h)*zb(k1,k4)*zb(k6,k56h)*iza(k1,k34h)**2*izb(
     &    k1,k56h) - za(k2,k5)*za(k2,k56h)*za(k3,k56h)*zb(k1,k4)*zb(k6,
     &    k34h)*iza(k1,k56h)**2*izb(k1,k34h) )
      triamp = triamp + y*s12**(-1) * ( za(k2,k34h)**2*za(k3,k34h)*za(
     &    k5,k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)*
     &    izb(k1,k56h)*yp1*ymr + za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)
     &    *zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k34h)*yp1*ymr + 
     &    za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h
     &    )*zb(k6,k34h)*iza(k1,k34h)*yp1**2 + za(k2,k34h)**2*za(k3,k34h
     &    )*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h
     &    )*yp1*xm1 + za(k2,k34h)**2*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h
     &    )*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k34h)*yp1*xm1 + 2._dp*za(k2,
     &    k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k56h)*yp1*ymr + 2._dp*za(k2,k34h)*za(
     &    k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(
     &    k6,k34h)*iza(k1,k34h)*xpr*ymr + za(k2,k34h)*za(k2,k56h)*za(k3
     &    ,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1
     &    ,k34h)*yp1*xpr + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k34h)*yp1*
     &    xpr )
      triamp = triamp + y*s12**(-1) * ( za(k2,k34h)*za(k2,k56h)*za(k3,
     &    k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k56h)*yp1*xpr + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*yp1*
     &    xpr + 2._dp*za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*
     &    zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*yp1*xm1 + 2.D
     &    0*za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)
     &    *zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)*xm1*xpr + 2._dp*za(k2,
     &    k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)**2*zb(k4,k34h)*zb(
     &    k6,k34h)*izb(k2,k34h)*yp1*ymr + za(k2,k34h)*za(k3,k34h)*za(k5
     &    ,k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k56h)*
     &    yp1*ymr + za(k2,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(
     &    k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k56h)*yp1*ymr + za(k2
     &    ,k34h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*izb(k2,k56h)*yp1**2 + za(k2,k34h)*za(k3,
     &    k34h)*za(k5,k56h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(
     &    k2,k34h)*yp1*xpr )
      triamp = triamp + y*s12**(-1) * ( za(k2,k34h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*izb(k2,
     &    k34h)*yp1*xpr + 2._dp*za(k2,k34h)*za(k3,k34h)*za(k5,k56h)*zb(
     &    k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k34h)*yp1
     &    *xm1 + za(k2,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k56h)*yp1*xm1 + za(k2,
     &    k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(
     &    k6,k34h)*izb(k2,k56h)*yp1*xm1 + za(k2,k56h)**2*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)*
     &    xpr*ymr + za(k2,k56h)**2*za(k3,k34h)*za(k5,k34h)*zb(k1,k56h)*
     &    zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*xpr*ymr + za(k2,k56h)**2
     &    *za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k56h)*izb(k1,k34h)*xm1*xpr + za(k2,k56h)**2*za(
     &    k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(
     &    k1,k56h)*xpr**2 + za(k2,k56h)**2*za(k3,k34h)*za(k5,k56h)*zb(
     &    k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)*xm1*xpr )
      triamp = triamp + y*s12**(-1) * ( za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k34h)*xpr
     &    *ymr + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k34h)*xpr*ymr + 2._dp*za(
     &    k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(
     &    k4,k56h)*zb(k6,k34h)*izb(k2,k56h)*xpr*ymr + za(k2,k56h)*za(k3
     &    ,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,
     &    k34h)*izb(k2,k56h)*yp1*xpr + za(k2,k56h)*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k56h)*yp1
     &    *xpr + za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k34h)*xpr**2 + za(k2,
     &    k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,
     &    k56h)*zb(k6,k34h)*izb(k2,k34h)*xm1*xpr + za(k2,k56h)*za(k3,
     &    k34h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*izb(
     &    k2,k34h)*xm1*xpr + 2._dp*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*
     &    zb(k1,k56h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k56h)*xm1*xpr )
      triamp = triamp + y*s12**(-1)*t234sum * (  - 2._dp*za(k2,k34h)*za(
     &    k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k34h)*iza(
     &    k1,k56h)*izb(k2,k34h)*yp1 - 2._dp*za(k2,k56h)*za(k3,k34h)*za(
     &    k5,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)*
     &    izb(k2,k56h)*xpr + za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)**2*zb(
     &    k4,k34h)*zb(k6,k34h)*izb(k2,k34h)**2*ymr + za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(k2,
     &    k56h)**2*ymr + za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h
     &    )*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k56h)**2*yp1 + za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*izb(k2,
     &    k56h)**2*yp1 + za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)**2*zb(k4,
     &    k56h)*zb(k6,k34h)*izb(k2,k34h)**2*xpr + za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*izb(k2,
     &    k34h)**2*xpr + za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h
     &    )*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k34h)**2*xm1 + za(k3,k34h)*
     &    za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,
     &    k56h)**2*xm1 )
      triamp = triamp + y*s12**(-1)*t234sum**2 * (  - za(k3,k34h)*za(k5
     &    ,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)*izb(
     &    k2,k56h)**2 - za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k34h)
     &    *zb(k6,k34h)*iza(k1,k56h)*izb(k2,k34h)**2 )
      triamp = triamp + y*s12**(-1)*t134sum * ( za(k2,k34h)**2*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k34h)**2*izb(k1,k56h)*ymr + za(k2,k34h)**2*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)**2*
     &    izb(k1,k56h)*yp1 + za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*zb(
     &    k4,k34h)*zb(k6,k34h)*iza(k1,k34h)**2*yp1 + za(k2,k34h)**2*za(
     &    k3,k34h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)**2*
     &    xm1 + 2._dp*za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)**2*xpr + 2._dp*za(k2,k34h
     &    )*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k34h)
     &    *iza(k1,k56h)**2*yp1 + za(k2,k56h)**2*za(k3,k34h)*za(k5,k34h)
     &    *zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)**2*ymr + za(k2,k56h)**2
     &    *za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*
     &    iza(k1,k56h)**2*izb(k1,k34h)*xpr + za(k2,k56h)**2*za(k3,k34h)
     &    *za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)
     &    **2*izb(k1,k34h)*xm1 )
      triamp = triamp + y*s12**(-1)*t134sum * ( za(k2,k56h)**2*za(k3,
     &    k34h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)**2*xpr
     &     )
      triamp = triamp + y*s12**(-1)*t134sum**2 * ( za(k2,k34h)**2*za(k3
     &    ,k34h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)**3*
     &    izb(k1,k56h) + za(k2,k56h)**2*za(k3,k34h)*za(k5,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k56h)**3*izb(k1,k34h) )
      triamp = triamp + y*mtsq*s12**(-1) * (  - za(k2,k5)*za(k2,k34h)*
     &    za(k3,k34h)*zb(k1,k4)*zb(k6,k34h)*iza(k1,k34h) - za(k2,k5)*
     &    za(k2,k56h)*za(k3,k34h)*zb(k1,k4)*zb(k6,k34h)*iza(k1,k56h) - 
     &    za(k2,k5)*za(k3,k34h)*zb(k1,k4)*zb(k1,k34h)*zb(k6,k34h)*izb(
     &    k2,k34h) - za(k2,k5)*za(k3,k34h)*zb(k1,k4)*zb(k1,k56h)*zb(k6,
     &    k34h)*izb(k2,k56h) )
      triamp = triamp + x*s12**(-1) * ( za(k2,k34h)**2*za(k3,k56h)*za(
     &    k5,k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*
     &    izb(k1,k56h)*yp1*ymr + za(k2,k34h)**2*za(k3,k56h)*za(k5,k34h)
     &    *zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)*yp1*ymr + 
     &    za(k2,k34h)**2*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k34h
     &    )*zb(k6,k56h)*iza(k1,k34h)*yp1**2 + za(k2,k34h)**2*za(k3,k56h
     &    )*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h
     &    )*yp1*xm1 + za(k2,k34h)**2*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h
     &    )*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)*yp1*xm1 + 2._dp*za(k2,
     &    k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,
     &    k34h)*zb(k6,k56h)*iza(k1,k56h)*yp1*ymr + 2._dp*za(k2,k34h)*za(
     &    k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(
     &    k6,k56h)*iza(k1,k34h)*xpr*ymr + za(k2,k34h)*za(k2,k56h)*za(k3
     &    ,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1
     &    ,k34h)*yp1*xpr + za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)*yp1*
     &    xpr )
      triamp = triamp + x*s12**(-1) * ( za(k2,k34h)*za(k2,k56h)*za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k56h)*yp1*xpr + za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k56h)*yp1*
     &    xpr + 2._dp*za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*
     &    zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k56h)*yp1*xm1 + 2.D
     &    0*za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)
     &    *zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*xm1*xpr + 2._dp*za(k2,
     &    k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)**2*zb(k4,k34h)*zb(
     &    k6,k56h)*izb(k2,k34h)*yp1*ymr + za(k2,k34h)*za(k3,k56h)*za(k5
     &    ,k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)*
     &    yp1*ymr + za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(
     &    k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k56h)*yp1*ymr + za(k2
     &    ,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,
     &    k34h)*zb(k6,k56h)*izb(k2,k56h)*yp1**2 + za(k2,k34h)*za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*izb(
     &    k2,k34h)*yp1*xpr )
      triamp = triamp + x*s12**(-1) * ( za(k2,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,
     &    k34h)*yp1*xpr + 2._dp*za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(
     &    k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k34h)*yp1
     &    *xm1 + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)*yp1*xm1 + za(k2,
     &    k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(
     &    k6,k56h)*izb(k2,k56h)*yp1*xm1 + za(k2,k56h)**2*za(k3,k56h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k56h)*
     &    xpr*ymr + za(k2,k56h)**2*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*
     &    zb(k4,k34h)*zb(k6,k56h)*iza(k1,k56h)*xpr*ymr + za(k2,k56h)**2
     &    *za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,
     &    k56h)*iza(k1,k56h)*izb(k1,k34h)*xm1*xpr + za(k2,k56h)**2*za(
     &    k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(
     &    k1,k56h)*xpr**2 + za(k2,k56h)**2*za(k3,k56h)*za(k5,k56h)*zb(
     &    k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k56h)*xm1*xpr )
      triamp = triamp + x*s12**(-1) * ( za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k34h)*xpr
     &    *ymr + za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k34h)*xpr*ymr + 2._dp*za(
     &    k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(
     &    k4,k56h)*zb(k6,k56h)*izb(k2,k56h)*xpr*ymr + za(k2,k56h)*za(k3
     &    ,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,
     &    k56h)*izb(k2,k56h)*yp1*xpr + za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k56h)*yp1
     &    *xpr + za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k34h)*xpr**2 + za(k2,
     &    k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,
     &    k56h)*zb(k6,k56h)*izb(k2,k34h)*xm1*xpr + za(k2,k56h)*za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k56h)*izb(
     &    k2,k34h)*xm1*xpr + 2._dp*za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*
     &    zb(k1,k56h)**2*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)*xm1*xpr )
      triamp = triamp + x*s12**(-1)*t234sum * (  - 2._dp*za(k2,k34h)*za(
     &    k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(
     &    k1,k56h)*izb(k2,k34h)*yp1 - 2._dp*za(k2,k56h)*za(k3,k56h)*za(
     &    k5,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*
     &    izb(k2,k56h)*xpr + za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)**2*zb(
     &    k4,k34h)*zb(k6,k56h)*izb(k2,k34h)**2*ymr + za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*izb(k2,
     &    k56h)**2*ymr + za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h
     &    )*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)**2*yp1 + za(k3,k56h)*
     &    za(k5,k34h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k56h)*izb(k2,
     &    k56h)**2*yp1 + za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)**2*zb(k4,
     &    k56h)*zb(k6,k56h)*izb(k2,k34h)**2*xpr + za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,
     &    k34h)**2*xpr + za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h
     &    )*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k34h)**2*xm1 + za(k3,k56h)*
     &    za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k56h)*zb(k6,k56h)*izb(k2,
     &    k56h)**2*xm1 )
      triamp = triamp + x*s12**(-1)*t234sum**2 * (  - za(k3,k56h)*za(k5
     &    ,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*izb(
     &    k2,k56h)**2 - za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k34h)
     &    *zb(k6,k56h)*iza(k1,k56h)*izb(k2,k34h)**2 )
      triamp = triamp + x*s12**(-1)*t134sum * ( za(k2,k34h)**2*za(k3,
     &    k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k34h)**2*izb(k1,k56h)*ymr + za(k2,k34h)**2*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*
     &    izb(k1,k56h)*yp1 + za(k2,k34h)**2*za(k3,k56h)*za(k5,k34h)*zb(
     &    k4,k34h)*zb(k6,k56h)*iza(k1,k34h)**2*yp1 + za(k2,k34h)**2*za(
     &    k3,k56h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*
     &    xm1 + 2._dp*za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*
     &    zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*xpr + 2._dp*za(k2,k34h
     &    )*za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k56h)
     &    *iza(k1,k56h)**2*yp1 + za(k2,k56h)**2*za(k3,k56h)*za(k5,k34h)
     &    *zb(k4,k34h)*zb(k6,k56h)*iza(k1,k56h)**2*ymr + za(k2,k56h)**2
     &    *za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*
     &    iza(k1,k56h)**2*izb(k1,k34h)*xpr + za(k2,k56h)**2*za(k3,k56h)
     &    *za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k56h)
     &    **2*izb(k1,k34h)*xm1 )
      triamp = triamp + x*s12**(-1)*t134sum * ( za(k2,k56h)**2*za(k3,
     &    k56h)*za(k5,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k56h)**2*xpr
     &     )
      triamp = triamp + x*s12**(-1)*t134sum**2 * ( za(k2,k34h)**2*za(k3
     &    ,k56h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**3*
     &    izb(k1,k56h) + za(k2,k56h)**2*za(k3,k56h)*za(k5,k56h)*zb(k4,
     &    k34h)*zb(k6,k56h)*iza(k1,k56h)**3*izb(k1,k34h) )
      triamp = triamp + x*mtsq*s12**(-1) * (  - za(k2,k5)*za(k2,k34h)*
     &    za(k3,k56h)*zb(k1,k4)*zb(k6,k56h)*iza(k1,k34h) - za(k2,k5)*
     &    za(k2,k56h)*za(k3,k56h)*zb(k1,k4)*zb(k6,k56h)*iza(k1,k56h) - 
     &    za(k2,k5)*za(k3,k56h)*zb(k1,k4)*zb(k1,k34h)*zb(k6,k56h)*izb(
     &    k2,k34h) - za(k2,k5)*za(k3,k56h)*zb(k1,k4)*zb(k1,k56h)*zb(k6,
     &    k56h)*izb(k2,k56h) )
      triamp = triamp + x*y*s12**(-1) * (  - 2._dp*za(k1,k34h)*za(k2,
     &    k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k56h)**2*yp1 - 2._dp*za(k1,k34h)*za(
     &    k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(
     &    k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*izb(k2,k34h)*yp1 - za(k1,
     &    k34h)*za(k2,k56h)**2*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(
     &    k4,k34h)*zb(k6,k34h)*iza(k1,k56h)**2*ymr - za(k1,k34h)*za(k2,
     &    k56h)**2*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*
     &    zb(k6,k34h)*iza(k1,k56h)**2*izb(k1,k34h)*xpr - za(k1,k34h)*
     &    za(k2,k56h)**2*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,
     &    k34h)*zb(k6,k34h)*iza(k1,k56h)**2*izb(k1,k34h)*xm1 - za(k1,
     &    k34h)*za(k2,k56h)**2*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(
     &    k4,k56h)*zb(k6,k34h)*iza(k1,k56h)**2*xpr + za(k1,k34h)*za(k3,
     &    k56h)*za(k5,k34h)*zb(k1,k34h)**2*zb(k1,k56h)*zb(k4,k34h)*zb(
     &    k6,k34h)*izb(k2,k34h)**2*ymr + za(k1,k34h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k1,k34h)**2*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(
     &    k2,k34h)**2*xpr )
      triamp = triamp + x*y*s12**(-1) * ( za(k1,k34h)*za(k3,k56h)*za(k5
     &    ,k56h)*zb(k1,k34h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*
     &    izb(k2,k34h)**2*xpr + za(k1,k34h)*za(k3,k56h)*za(k5,k56h)*zb(
     &    k1,k34h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k34h)
     &    **2*xm1 - za(k1,k56h)*za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*izb(k1
     &    ,k56h)*ymr - za(k1,k56h)*za(k2,k34h)**2*za(k3,k34h)*za(k5,
     &    k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*
     &    izb(k1,k56h)*yp1 - za(k1,k56h)*za(k2,k34h)**2*za(k3,k34h)*za(
     &    k5,k34h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)**2*
     &    yp1 - za(k1,k56h)*za(k2,k34h)**2*za(k3,k34h)*za(k5,k56h)*zb(
     &    k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*xm1 - 2._dp*
     &    za(k1,k56h)*za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*
     &    zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*xpr - 2._dp
     &    *za(k1,k56h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*
     &    zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*izb(k2,k56h)
     &    *xpr )
      triamp = triamp + x*y*s12**(-1) * ( za(k1,k56h)*za(k3,k34h)*za(k5
     &    ,k34h)*zb(k1,k34h)**2*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*
     &    izb(k2,k56h)**2*ymr + za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(
     &    k1,k34h)**2*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)
     &    **2*yp1 + za(k1,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(
     &    k1,k56h)**2*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k56h)**2*yp1 + za(
     &    k1,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)**2*
     &    zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)**2*xm1 + za(k2,k34h)**2*
     &    za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h
     &    )*iza(k1,k34h)*izb(k1,k56h)*yp1 + za(k2,k34h)**2*za(k3,k34h)*
     &    za(k5,k56h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)*
     &    yp1 + za(k2,k34h)**2*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(
     &    k4,k56h)*zb(k6,k34h)*iza(k1,k34h)*ymr + za(k2,k34h)**2*za(k3,
     &    k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k34h)*yp1 + za(k2,k34h)**2*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h
     &    )*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k34h)*yp1 )
      triamp = triamp + x*y*s12**(-1) * ( za(k2,k34h)**2*za(k3,k56h)*
     &    za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)*
     &    xm1 + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)**2*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*izb(k1,k56h)*
     &    ymr + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k34h)*ymr + 2._dp*za(k2,
     &    k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,
     &    k34h)*zb(k6,k56h)*iza(k1,k34h)*yp1 + 2._dp*za(k2,k34h)*za(k2,
     &    k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k34h)*zb(k6,
     &    k56h)*iza(k1,k56h)*yp1 + 2._dp*za(k2,k34h)*za(k2,k56h)*za(k3,
     &    k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k34h)*xpr + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*
     &    zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)*xm1 + za(k2,
     &    k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,
     &    k34h)*zb(k6,k56h)*iza(k1,k34h)*xm1 + za(k2,k34h)*za(k2,k56h)*
     &    za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*
     &    iza(k1,k56h)*ymr )
      triamp = triamp + x*y*s12**(-1) * ( za(k2,k34h)*za(k2,k56h)*za(k3
     &    ,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1
     &    ,k56h)*ymr + 2._dp*za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*yp1 + 
     &    2._dp*za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,
     &    k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)*xpr + za(k2,k34h)*
     &    za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h
     &    )*zb(k6,k34h)*iza(k1,k56h)*izb(k1,k34h)*xm1 + 2._dp*za(k2,k34h
     &    )*za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)
     &    *zb(k6,k34h)*iza(k1,k56h)*xpr + za(k2,k34h)*za(k2,k56h)*za(k3
     &    ,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1
     &    ,k56h)*xm1 + 2._dp*za(k2,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,
     &    k34h)**2*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k34h)*yp1 + za(k2,
     &    k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(
     &    k6,k56h)*izb(k2,k56h)*yp1 + za(k2,k34h)*za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,
     &    k56h)*yp1 )
      triamp = triamp + x*y*s12**(-1) * ( za(k2,k34h)*za(k3,k56h)*za(k5
     &    ,k34h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k34h)*
     &    ymr + za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k34h)*ymr + 2._dp*za(k2,
     &    k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*izb(k2,k34h)*yp1 + za(k2,k34h)*za(k3,k56h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*
     &    izb(k2,k56h)*ymr + za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,
     &    k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k56h)*yp1 + 
     &    za(k2,k34h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)**2*zb(k4,k34h
     &    )*zb(k6,k34h)*izb(k2,k56h)*yp1 + 2._dp*za(k2,k34h)*za(k3,k56h)
     &    *za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*
     &    izb(k2,k34h)*xpr + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,
     &    k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k34h)*xm1 + 
     &    za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h
     &    )*zb(k6,k34h)*izb(k2,k34h)*xm1 )
      triamp = triamp + x*y*s12**(-1) * ( za(k2,k34h)*za(k3,k56h)*za(k5
     &    ,k56h)*zb(k1,k56h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k56h)*
     &    xm1 + za(k2,k56h)**2*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(
     &    k4,k34h)*zb(k6,k56h)*iza(k1,k56h)*ymr + za(k2,k56h)**2*za(k3,
     &    k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k56h)*xpr + za(k2,k56h)**2*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h
     &    )*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k56h)*xpr + za(k2,k56h)**2*
     &    za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*
     &    iza(k1,k56h)*xm1 + za(k2,k56h)**2*za(k3,k56h)*za(k5,k34h)*zb(
     &    k1,k56h)**2*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*izb(k1,k34h)
     &    *xpr + za(k2,k56h)**2*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(
     &    k4,k56h)*zb(k6,k34h)*iza(k1,k56h)*xpr + za(k2,k56h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k34h)**2*zb(k4,k34h)*zb(k6,k56h)*izb(
     &    k2,k34h)*ymr + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h
     &    )**2*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)*ymr + za(k2,k56h)*
     &    za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*
     &    zb(k6,k56h)*izb(k2,k56h)*ymr )
      triamp = triamp + x*y*s12**(-1) * ( 2._dp*za(k2,k56h)*za(k3,k34h)*
     &    za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*
     &    izb(k2,k56h)*yp1 + za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,
     &    k34h)**2*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k34h)*xpr + za(k2,
     &    k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,
     &    k34h)*zb(k6,k56h)*izb(k2,k34h)*xpr + za(k2,k56h)*za(k3,k34h)*
     &    za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*
     &    izb(k2,k34h)*xm1 + 2._dp*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*
     &    zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)*
     &    xpr + za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,
     &    k56h)*zb(k4,k56h)*zb(k6,k56h)*izb(k2,k56h)*xm1 + za(k2,k56h)*
     &    za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)**2*zb(k4,k34h)*zb(k6,k56h
     &    )*izb(k2,k56h)*xm1 + za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(
     &    k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k34h)*xpr
     &     + za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)**2*zb(k4,
     &    k34h)*zb(k6,k34h)*izb(k2,k34h)*xpr )
      triamp = triamp + x*y*s12**(-1) * ( 2._dp*za(k2,k56h)*za(k3,k56h)*
     &    za(k5,k34h)*zb(k1,k56h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,
     &    k56h)*xpr )
      triamp = triamp + x*y*s12**(-1)*t234sum * (  - 2._dp*za(k1,k34h)*
     &    za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k34h)*
     &    zb(k6,k34h)*iza(k1,k56h)*izb(k2,k34h)**2 - 2._dp*za(k1,k56h)*
     &    za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*
     &    zb(k6,k56h)*iza(k1,k34h)*izb(k2,k56h)**2 - za(k2,k34h)*za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k56h)*izb(k2,k34h) - za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(
     &    k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)*izb(k2,k34h) - 
     &    za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,k56h)*
     &    zb(k6,k56h)*iza(k1,k34h)*izb(k2,k56h) - za(k2,k56h)*za(k3,
     &    k34h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,
     &    k34h)*izb(k2,k56h) + za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)**2*
     &    zb(k4,k34h)*zb(k6,k56h)*izb(k2,k34h)**2 + za(k3,k34h)*za(k5,
     &    k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*izb(k2,
     &    k56h)**2 + za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*
     &    zb(k4,k34h)*zb(k6,k34h)*izb(k2,k34h)**2 )
      triamp = triamp + x*y*s12**(-1)*t234sum * ( za(k3,k56h)*za(k5,
     &    k34h)*zb(k1,k56h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,k56h)**2
     &     )
      triamp = triamp + x*y*s12**(-1)*t134sum * (  - 2._dp*za(k1,k34h)*
     &    za(k2,k56h)**2*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h
     &    )*zb(k6,k34h)*iza(k1,k56h)**3*izb(k1,k34h) - 2._dp*za(k1,k56h)
     &    *za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k4,
     &    k56h)*zb(k6,k56h)*iza(k1,k34h)**3*izb(k1,k56h) + za(k2,k34h)
     &    **2*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,
     &    k56h)*iza(k1,k34h)**2*izb(k1,k56h) + za(k2,k34h)**2*za(k3,
     &    k56h)*za(k5,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)**2 + 
     &    za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*
     &    zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2*izb(k1,k56h) + za(k2,
     &    k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k4,k34h)*zb(k6,
     &    k56h)*iza(k1,k34h)**2 + za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h)
     &    **2*izb(k1,k34h) + za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,
     &    k56h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)**2 + za(k2,k56h)**
     &    2*za(k3,k34h)*za(k5,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k56h
     &    )**2 )
      triamp = triamp + x*y*s12**(-1)*t134sum * ( za(k2,k56h)**2*za(k3,
     &    k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,
     &    k56h)**2*izb(k1,k34h) )
      triamp = triamp + x*y**2*s12**(-1) * (  - za(k1,k34h)*za(k2,k56h)
     &    **2*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,
     &    k34h)*iza(k1,k56h)**2 + za(k1,k34h)*za(k3,k34h)*za(k5,k56h)*
     &    zb(k1,k34h)**2*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*izb(k2,
     &    k34h)**2 - za(k1,k56h)*za(k2,k34h)**2*za(k3,k34h)*za(k5,k34h)
     &    *zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h)**2 + za(k1,
     &    k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)**2*zb(
     &    k4,k56h)*zb(k6,k34h)*izb(k2,k56h)**2 + za(k2,k34h)**2*za(k3,
     &    k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,
     &    k34h) + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1
     &    ,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k34h) + za(k2,k34h)*za(
     &    k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(
     &    k6,k34h)*iza(k1,k34h) + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*
     &    za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k34h)*iza(k1,k56h)
     &     + za(k2,k34h)*za(k2,k56h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k56h
     &    )*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h) )
      triamp = triamp + x*y**2*s12**(-1) * ( za(k2,k34h)*za(k3,k34h)*
     &    za(k5,k56h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k34h)*izb(k2,
     &    k34h) + za(k2,k34h)*za(k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1
     &    ,k56h)*zb(k4,k34h)*zb(k6,k34h)*izb(k2,k34h) + za(k2,k34h)*za(
     &    k3,k34h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(
     &    k6,k34h)*izb(k2,k56h) + za(k2,k56h)**2*za(k3,k34h)*za(k5,k34h
     &    )*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k34h)*iza(k1,k56h) + za(k2,
     &    k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,
     &    k34h)*zb(k6,k34h)*izb(k2,k34h) + za(k2,k56h)*za(k3,k34h)*za(
     &    k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k34h)*izb(
     &    k2,k56h) + za(k2,k56h)*za(k3,k34h)*za(k5,k34h)*zb(k1,k56h)**2
     &    *zb(k4,k34h)*zb(k6,k34h)*izb(k2,k56h) )
      triamp = triamp + x**2*y*s12**(-1) * (  - za(k1,k34h)*za(k2,k56h)
     &    **2*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h)*zb(k4,k34h)*zb(k6,
     &    k56h)*iza(k1,k56h)**2 + za(k1,k34h)*za(k3,k56h)*za(k5,k56h)*
     &    zb(k1,k34h)**2*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,
     &    k34h)**2 - za(k1,k56h)*za(k2,k34h)**2*za(k3,k56h)*za(k5,k34h)
     &    *zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h)**2 + za(k1,
     &    k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)**2*zb(
     &    k4,k56h)*zb(k6,k56h)*izb(k2,k56h)**2 + za(k2,k34h)**2*za(k3,
     &    k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,
     &    k34h) + za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1
     &    ,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k34h) + za(k2,k34h)*za(
     &    k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)*zb(k4,k34h)*zb(
     &    k6,k56h)*iza(k1,k34h) + za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k1,k34h)*zb(k4,k56h)*zb(k6,k56h)*iza(k1,k56h)
     &     + za(k2,k34h)*za(k2,k56h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k56h
     &    )*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k56h) )
      triamp = triamp + x**2*y*s12**(-1) * ( za(k2,k34h)*za(k3,k56h)*
     &    za(k5,k56h)*zb(k1,k34h)**2*zb(k4,k56h)*zb(k6,k56h)*izb(k2,
     &    k34h) + za(k2,k34h)*za(k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k1
     &    ,k56h)*zb(k4,k34h)*zb(k6,k56h)*izb(k2,k34h) + za(k2,k34h)*za(
     &    k3,k56h)*za(k5,k56h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(
     &    k6,k56h)*izb(k2,k56h) + za(k2,k56h)**2*za(k3,k56h)*za(k5,k34h
     &    )*zb(k1,k56h)*zb(k4,k34h)*zb(k6,k56h)*iza(k1,k56h) + za(k2,
     &    k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,
     &    k34h)*zb(k6,k56h)*izb(k2,k34h) + za(k2,k56h)*za(k3,k56h)*za(
     &    k5,k34h)*zb(k1,k34h)*zb(k1,k56h)*zb(k4,k56h)*zb(k6,k56h)*izb(
     &    k2,k56h) + za(k2,k56h)*za(k3,k56h)*za(k5,k34h)*zb(k1,k56h)**2
     &    *zb(k4,k34h)*zb(k6,k56h)*izb(k2,k56h) )

      apm=triamp*(-half*im)/(s(k3,k4)*s(k5,k6))

      return
      end
