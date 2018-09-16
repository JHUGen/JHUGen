      function amp_dkqg_tbqdk_g(p1,p2,p3,p4,
     &  p7,p8,k5,k6,e6,mbdk)
      implicit none
      include 'types.f'
      real(dp):: amp_dkqg_tbqdk_g

c---- 4-flavour routine for gluon radiation in top decay
c---- averaged over initial colours and spins
c---- line is contracted with the vector n(mu)
c---- and routine allows for massive vectors p5,p6
c     q(-p1)+g(-p2)--> t(p3,p4,p5,p8)+bb(p6)+q'(p7)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'masses.f'
      integer:: p1,p2,p3,p4,p7,p8,b8,k5,k6,e6,i1,i2,b2,i,j,k,m
      complex(dp):: amp(2,2,2,2),za26b,za345b,za167b
      real(dp):: c6,prop26,prop167,s167,invtwoptDp8,prop58,
     & s17,s34,mbdk

c--- statement functions
      za167b(i1,i2)=za(i1,p1)*zb(p1,i2)+za(i1,p7)*zb(p7,i2)
     &             +za(i1,k6)*zb(k6,i2)+za(i1,e6)*zb(e6,i2)*c6
      za345b(i1,i2)=-za167b(i1,i2)-za(i1,p2)*zb(p2,i2)
     &                            -za(i1,p8)*zb(p8,i2)
      za26b(i1,i2)=za(i1,p2)*zb(p2,i2)
     &            +za(i1,k6)*zb(k6,i2)+za(i1,e6)*zb(e6,i2)*c6

      c6=mb**2/real(za(k6,e6)*zb(e6,k6))

      s167=real(za(p1,p7)*zb(p7,p1)
     &         +za(p1,k6)*zb(k6,p1)+c6*za(p1,e6)*zb(e6,p1)
     &         +za(p7,k6)*zb(k6,p7)+c6*za(p7,e6)*zb(e6,p7))
     &         +mb**2
      prop26=real(za(p2,k6)*zb(k6,p2)+c6*za(p2,e6)*zb(e6,p2))
      prop167=s167-mt**2
      invtwoptDp8=1._dp/real(za345b(p8,p8))
      prop58=real(za345b(p8,p8)-za(p8,p3)*zb(p3,p8)-za(p8,p4)*zb(p4,p8))

c--- choice of gauge vector
      b2=p1
c--- choice of gauge vector
      b8=p2

      amp(2,2,2,2)=czip

      amp(2,2,2,1)= + prop26**(-1)*prop58**(-1) * (  -1._dp/(za(p2,b2))/(
     &    za(p8,b8))*za(p3,k5)*za(p7,p8)*za(b8,k5)*zb(p2,k6)*zb(p4,p8)*
     &    zb(p8,k5)*za26b(b2,p1) +1._dp/(za(p2,b2))/(za(p8,b8))*za(p3,k5)*
     &    za(b8,k5)*zb(p2,k6)*zb(p8,k5)*za345b(p7,p4)*za26b(b2,p1) +
     &   1._dp/(za(p2,b2))/(za(p8,b8))/(zb(p3,k5))*za(p3,b8)*za(p7,p8)*zb(
     &    p2,k6)*zb(p3,p8)*zb(p4,p8)*za26b(b2,p1)*mbdk**2 -1._dp/(za(p2,b2
     &    ))/(za(p8,b8))/(zb(p3,k5))*za(p3,b8)*zb(p2,k6)*zb(p3,p8)*
     &    za345b(p7,p4)*za26b(b2,p1)*mbdk**2 )
      amp(2,2,2,1) = amp(2,2,2,1) + prop26**(-1)*invtwoptDp8 * (1._dp/(za(
     &    p2,b2))/(za(p8,b8))*za(p3,k5)*zb(p2,k6)*za345b(p7,p8)*za345b(
     &    b8,p4)*za26b(b2,p1) )
      amp(2,2,2,1) = amp(2,2,2,1) + prop167**(-1)*prop58**(-1) * (  -
     &   1._dp/(za(p2,b2))/(za(p8,b8))*za(p3,k5)*za(b2,p8)*za(b8,k5)*zb(p1
     &    ,k6)*zb(p4,p8)*zb(p8,k5)*za167b(p7,p2) +1._dp/(za(p2,b2))/(za(p8
     &    ,b8))*za(p3,k5)*za(b8,k5)*zb(p1,k6)*zb(p8,k5)*za167b(p7,p2)*
     &    za345b(b2,p4) +1._dp/(za(p2,b2))/(za(p8,b8))/(zb(p3,k5))*za(p3,
     &    b8)*za(b2,p8)*zb(p1,k6)*zb(p3,p8)*zb(p4,p8)*za167b(p7,p2)*
     &    mbdk**2 -1._dp/(za(p2,b2))/(za(p8,b8))/(zb(p3,k5))*za(p3,b8)*zb(
     &    p1,k6)*zb(p3,p8)*za167b(p7,p2)*za345b(b2,p4)*mbdk**2 )
      amp(2,2,2,1) = amp(2,2,2,1) + prop167**(-1)*invtwoptDp8 * (1._dp/(
     &    za(p2,b2))/(za(p8,b8))*za(p3,k5)*zb(p1,k6)*za167b(p7,p2)*
     &    za345b(b2,p8)*za345b(b8,p4) )
      amp(2,2,2,1) = amp(2,2,2,1) + mb**2*prop26**(-1)*prop58**(-1)
     &  * (  -1._dp/(za(p2,b2))/(za(p8,b8))/(za(k6,e6))*za(p3,k5)*za(b2,e6
     &    )*za(p7,p8)*za(b8,k5)*zb(p1,p2)*zb(p4,p8)*zb(p8,k5) +1._dp/(za(
     &    p2,b2))/(za(p8,b8))/(za(k6,e6))*za(p3,k5)*za(b2,e6)*za(b8,k5)
     &    *zb(p1,p2)*zb(p8,k5)*za345b(p7,p4) +1._dp/(za(p2,b2))/(za(p8,b8)
     &    )/(za(k6,e6))/(zb(p3,k5))*za(p3,b8)*za(b2,e6)*za(p7,p8)*zb(p1
     &    ,p2)*zb(p3,p8)*zb(p4,p8)*mbdk**2 -1._dp/(za(p2,b2))/(za(p8,b8))
     &    /(za(k6,e6))/(zb(p3,k5))*za(p3,b8)*za(b2,e6)*zb(p1,p2)*zb(p3,
     &    p8)*za345b(p7,p4)*mbdk**2 )
      amp(2,2,2,1) = amp(2,2,2,1) + mb**2*prop26**(-1)*invtwoptDp8 * (
     &   1._dp/(za(p2,b2))/(za(p8,b8))/(za(k6,e6))*za(p3,k5)*za(b2,e6)*zb(
     &    p1,p2)*za345b(p7,p8)*za345b(b8,p4) )
      amp(2,2,2,1) = amp(2,2,2,1) + mt**2*prop26**(-1)*invtwoptDp8 * (
     &     -1._dp/(za(p2,b2))/(za(p8,b8))*za(p3,k5)*za(p7,b8)*zb(p2,k6)*
     &    zb(p4,p8)*za26b(b2,p1) )
      amp(2,2,2,1) = amp(2,2,2,1) + mt**2*prop167**(-1)*prop58**(-1)
     &  * (1._dp/(za(p2,b2))/(za(p8,b8))*za(p3,k5)*za(b2,p7)*za(b8,k5)*zb(
     &    p1,k6)*zb(p2,p4)*zb(p8,k5) -1._dp/(za(p2,b2))/(za(p8,b8))/(zb(p3
     &    ,k5))*za(p3,b8)*za(b2,p7)*zb(p1,k6)*zb(p2,p4)*zb(p3,p8)*
     &    mbdk**2 )
      amp(2,2,2,1) = amp(2,2,2,1) + mt**2*prop167**(-1)*invtwoptDp8
     &  * (  -1._dp/(za(p2,b2))*za(p3,k5)*za(b2,p7)*zb(p1,k6)*zb(p2,p8)*
     &    zb(p4,p8) +1._dp/(za(p2,b2))/(za(p8,b8))*za(p3,k5)*za(b2,p7)*zb(
     &    p1,k6)*zb(p2,p8)*za345b(b8,p4) -1._dp/(za(p2,b2))/(za(p8,b8))*
     &    za(p3,k5)*za(b2,p7)*zb(p1,k6)*zb(p4,p8)*za345b(b8,p2) -1._dp/(
     &    za(p2,b2))/(za(p8,b8))*za(p3,k5)*za(b2,b8)*zb(p1,k6)*zb(p4,p8
     &    )*za167b(p7,p2) )
      amp(2,2,2,1) = amp(2,2,2,1) + mt**2*mb**2*prop26**(-1)*
     & invtwoptDp8 * (  -1._dp/(za(p2,b2))/(za(p8,b8))/(za(k6,e6))*za(p3,
     &    k5)*za(b2,e6)*za(p7,b8)*zb(p1,p2)*zb(p4,p8) )

      amp(2,2,1,2)= + prop26**(-1)*prop58**(-1) * (1._dp/(za(p2,b2))/(za(
     &    p3,k5))*za(p3,p8)**2*za(p7,p8)*zb(p2,k6)*zb(p4,p8)*za26b(b2,
     &    p1)*mbdk -1._dp/(za(p2,b2))/(za(p3,k5))*za(p3,p8)**2*zb(p2,k6)*
     &    za345b(p7,p4)*za26b(b2,p1)*mbdk )
      amp(2,2,1,2) = amp(2,2,1,2) + prop167**(-1)*prop58**(-1) * (1._dp/(
     &    za(p2,b2))/(za(p3,k5))*za(p3,p8)**2*za(b2,p8)*zb(p1,k6)*zb(p4
     &    ,p8)*za167b(p7,p2)*mbdk -1._dp/(za(p2,b2))/(za(p3,k5))*za(p3,p8)
     &    **2*zb(p1,k6)*za167b(p7,p2)*za345b(b2,p4)*mbdk )
      amp(2,2,1,2) = amp(2,2,1,2) + mb**2*prop26**(-1)*prop58**(-1)
     &  * (1._dp/(za(p2,b2))/(za(p3,k5))/(za(k6,e6))*za(p3,p8)**2*za(b2,e6
     &    )*za(p7,p8)*zb(p1,p2)*zb(p4,p8)*mbdk -1._dp/(za(p2,b2))/(za(p3,
     &    k5))/(za(k6,e6))*za(p3,p8)**2*za(b2,e6)*zb(p1,p2)*za345b(p7,
     &    p4)*mbdk )
      amp(2,2,1,2) = amp(2,2,1,2) + mt**2*prop167**(-1)*prop58**(-1)
     &  * (  -1._dp/(za(p2,b2))/(za(p3,k5))*za(p3,p8)**2*za(b2,p7)*zb(p1,
     &    k6)*zb(p2,p4)*mbdk )

      amp(2,2,1,1)= + prop26**(-1)*prop58**(-1) * (  -1._dp/(za(p2,b2))*
     &    za(p3,p8)*za(p7,p8)*za(p8,k5)*zb(p2,k6)*zb(p4,p8)*za26b(b2,p1
     &    ) +1._dp/(za(p2,b2))*za(p3,p8)*za(p8,k5)*zb(p2,k6)*za345b(p7,p4)
     &    *za26b(b2,p1) -1._dp/(za(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p3,
     &    p8)*za(p7,p8)*zb(p2,k6)*zb(p3,b8)*zb(p4,p8)*za26b(b2,p1)*
     &    mbdk**2 +1._dp/(za(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p3,p8)*zb(
     &    p2,k6)*zb(p3,b8)*za345b(p7,p4)*za26b(b2,p1)*mbdk**2 +1._dp/(za(
     &    p2,b2))/(zb(p8,b8))*za(p3,k5)*za(p7,p8)*za(p8,k5)*zb(p2,k6)*
     &    zb(p4,p8)*zb(b8,k5)*za26b(b2,p1) -1._dp/(za(p2,b2))/(zb(p8,b8))*
     &    za(p3,k5)*za(p8,k5)*zb(p2,k6)*zb(b8,k5)*za345b(p7,p4)*za26b(
     &    b2,p1) )
      amp(2,2,1,1) = amp(2,2,1,1) + prop26**(-1)*invtwoptDp8 * (  -1._dp/(
     &    za(p2,b2))*za(p3,k5)*za(p7,p8)*zb(p2,k6)*za345b(p8,p4)*za26b(
     &    b2,p1) -1._dp/(za(p2,b2))/(zb(p8,b8))*za(p3,k5)*zb(p2,k6)*
     &    za345b(p7,b8)*za345b(p8,p4)*za26b(b2,p1) )
      amp(2,2,1,1) = amp(2,2,1,1) + prop167**(-1)*prop58**(-1) * (  -
     &   1._dp/(za(p2,b2))*za(p3,p8)*za(b2,p8)*za(p8,k5)*zb(p1,k6)*zb(p4,
     &    p8)*za167b(p7,p2) +1._dp/(za(p2,b2))*za(p3,p8)*za(p8,k5)*zb(p1,
     &    k6)*za167b(p7,p2)*za345b(b2,p4) -1._dp/(za(p2,b2))/(zb(p3,k5))/(
     &    zb(p8,b8))*za(p3,p8)*za(b2,p8)*zb(p1,k6)*zb(p3,b8)*zb(p4,p8)*
     &    za167b(p7,p2)*mbdk**2 +1._dp/(za(p2,b2))/(zb(p3,k5))/(zb(p8,b8))
     &    *za(p3,p8)*zb(p1,k6)*zb(p3,b8)*za167b(p7,p2)*za345b(b2,p4)*
     &    mbdk**2 +1._dp/(za(p2,b2))/(zb(p8,b8))*za(p3,k5)*za(b2,p8)*za(p8
     &    ,k5)*zb(p1,k6)*zb(p4,p8)*zb(b8,k5)*za167b(p7,p2) -1._dp/(za(p2,
     &    b2))/(zb(p8,b8))*za(p3,k5)*za(p8,k5)*zb(p1,k6)*zb(b8,k5)*
     &    za167b(p7,p2)*za345b(b2,p4) )
      amp(2,2,1,1) = amp(2,2,1,1) + prop167**(-1)*invtwoptDp8 * (  -
     &   1._dp/(za(p2,b2))*za(p3,k5)*za(b2,p8)*zb(p1,k6)*za167b(p7,p2)*
     &    za345b(p8,p4) -1._dp/(za(p2,b2))/(zb(p8,b8))*za(p3,k5)*zb(p1,k6)
     &    *za167b(p7,p2)*za345b(b2,b8)*za345b(p8,p4) )
      amp(2,2,1,1) = amp(2,2,1,1) + mb**2*prop26**(-1)*prop58**(-1)
     &  * (  -1._dp/(za(p2,b2))/(za(k6,e6))*za(p3,p8)*za(b2,e6)*za(p7,p8)*
     &    za(p8,k5)*zb(p1,p2)*zb(p4,p8) +1._dp/(za(p2,b2))/(za(k6,e6))*za(
     &    p3,p8)*za(b2,e6)*za(p8,k5)*zb(p1,p2)*za345b(p7,p4) -1._dp/(za(p2
     &    ,b2))/(za(k6,e6))/(zb(p3,k5))/(zb(p8,b8))*za(p3,p8)*za(b2,e6)
     &    *za(p7,p8)*zb(p1,p2)*zb(p3,b8)*zb(p4,p8)*mbdk**2 +1._dp/(za(p2,
     &    b2))/(za(k6,e6))/(zb(p3,k5))/(zb(p8,b8))*za(p3,p8)*za(b2,e6)*
     &    zb(p1,p2)*zb(p3,b8)*za345b(p7,p4)*mbdk**2 +1._dp/(za(p2,b2))/(
     &    za(k6,e6))/(zb(p8,b8))*za(p3,k5)*za(b2,e6)*za(p7,p8)*za(p8,k5
     &    )*zb(p1,p2)*zb(p4,p8)*zb(b8,k5) -1._dp/(za(p2,b2))/(za(k6,e6))/(
     &    zb(p8,b8))*za(p3,k5)*za(b2,e6)*za(p8,k5)*zb(p1,p2)*zb(b8,k5)*
     &    za345b(p7,p4) )
      amp(2,2,1,1) = amp(2,2,1,1) + mb**2*prop26**(-1)*invtwoptDp8 * (
     &     -1._dp/(za(p2,b2))/(za(k6,e6))*za(p3,k5)*za(b2,e6)*za(p7,p8)*
     &    zb(p1,p2)*za345b(p8,p4) -1._dp/(za(p2,b2))/(za(k6,e6))/(zb(p8,b8
     &    ))*za(p3,k5)*za(b2,e6)*zb(p1,p2)*za345b(p7,b8)*za345b(p8,p4)
     &     )
      amp(2,2,1,1) = amp(2,2,1,1) + mt**2*prop26**(-1)*invtwoptDp8 * (
     &   1._dp/(za(p2,b2))/(zb(p8,b8))*za(p3,k5)*za(p7,p8)*zb(p2,k6)*zb(p4
     &    ,b8)*za26b(b2,p1) )
      amp(2,2,1,1) = amp(2,2,1,1) + mt**2*prop167**(-1)*prop58**(-1)
     &  * (1._dp/(za(p2,b2))*za(p3,p8)*za(b2,p7)*za(p8,k5)*zb(p1,k6)*zb(p2
     &    ,p4) +1._dp/(za(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p3,p8)*za(b2,
     &    p7)*zb(p1,k6)*zb(p2,p4)*zb(p3,b8)*mbdk**2 -1._dp/(za(p2,b2))/(
     &    zb(p8,b8))*za(p3,k5)*za(b2,p7)*za(p8,k5)*zb(p1,k6)*zb(p2,p4)*
     &    zb(b8,k5) )
      amp(2,2,1,1) = amp(2,2,1,1) + mt**2*prop167**(-1)*invtwoptDp8
     &  * (  -1._dp/(za(p2,b2))/(zb(p8,b8))*za(p3,k5)*za(b2,p7)*zb(p1,k6)*
     &    zb(p2,b8)*za345b(p8,p4) +1._dp/(za(p2,b2))/(zb(p8,b8))*za(p3,k5)
     &    *za(b2,p7)*zb(p1,k6)*zb(p4,b8)*za345b(p8,p2) +1._dp/(za(p2,b2))
     &    /(zb(p8,b8))*za(p3,k5)*za(b2,p8)*zb(p1,k6)*zb(p4,b8)*za167b(
     &    p7,p2) )
      amp(2,2,1,1) = amp(2,2,1,1) + mt**2*mb**2*prop26**(-1)*
     & invtwoptDp8 * (1._dp/(za(p2,b2))/(za(k6,e6))/(zb(p8,b8))*za(p3,k5)*
     &    za(b2,e6)*za(p7,p8)*zb(p1,p2)*zb(p4,b8) )

      amp(2,1,2,2)=czip

      amp(2,1,2,1)= + prop26**(-1)*prop58**(-1) * (1._dp/(za(p8,b8))/(zb(
     &    p2,b2))*za(p3,k5)*za(p7,p8)*za(b8,k5)*zb(p4,p8)*zb(b2,k6)*zb(
     &    p8,k5)*za26b(p2,p1) -1._dp/(za(p8,b8))/(zb(p2,b2))*za(p3,k5)*za(
     &    b8,k5)*zb(b2,k6)*zb(p8,k5)*za345b(p7,p4)*za26b(p2,p1) -1._dp/(
     &    za(p8,b8))/(zb(p2,b2))/(zb(p3,k5))*za(p3,b8)*za(p7,p8)*zb(p3,
     &    p8)*zb(p4,p8)*zb(b2,k6)*za26b(p2,p1)*mbdk**2 +1._dp/(za(p8,b8))
     &    /(zb(p2,b2))/(zb(p3,k5))*za(p3,b8)*zb(p3,p8)*zb(b2,k6)*
     &    za345b(p7,p4)*za26b(p2,p1)*mbdk**2 )
      amp(2,1,2,1) = amp(2,1,2,1) + prop26**(-1)*invtwoptDp8 * (  -1._dp/(
     &    za(p8,b8))/(zb(p2,b2))*za(p3,k5)*zb(b2,k6)*za345b(p7,p8)*
     &    za345b(b8,p4)*za26b(p2,p1) )
      amp(2,1,2,1) = amp(2,1,2,1) + prop167**(-1)*prop58**(-1) * (1._dp/(
     &    za(p8,b8))/(zb(p2,b2))*za(p2,p8)*za(p3,k5)*za(b8,k5)*zb(p1,k6
     &    )*zb(p4,p8)*zb(p8,k5)*za167b(p7,b2) -1._dp/(za(p8,b8))/(zb(p2,b2
     &    ))*za(p3,k5)*za(b8,k5)*zb(p1,k6)*zb(p8,k5)*za167b(p7,b2)*
     &    za345b(p2,p4) -1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(p3,k5))*za(p2,
     &    p8)*za(p3,b8)*zb(p1,k6)*zb(p3,p8)*zb(p4,p8)*za167b(p7,b2)*
     &    mbdk**2 +1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(p3,k5))*za(p3,b8)*zb(
     &    p1,k6)*zb(p3,p8)*za167b(p7,b2)*za345b(p2,p4)*mbdk**2 )
      amp(2,1,2,1) = amp(2,1,2,1) + prop167**(-1)*invtwoptDp8 * (  -
     &   1._dp/(za(p8,b8))/(zb(p2,b2))*za(p3,k5)*zb(p1,k6)*za167b(p7,b2)*
     &    za345b(p2,p8)*za345b(b8,p4) )
      amp(2,1,2,1) = amp(2,1,2,1) + mb**2*prop26**(-1)*prop58**(-1)
     &  * (1._dp/(za(p8,b8))/(za(k6,e6))/(zb(p2,b2))*za(p2,e6)*za(p3,k5)*
     &    za(p7,p8)*za(b8,k5)*zb(p1,b2)*zb(p4,p8)*zb(p8,k5) -1._dp/(za(p8,
     &    b8))/(za(k6,e6))/(zb(p2,b2))*za(p2,e6)*za(p3,k5)*za(b8,k5)*
     &    zb(p1,b2)*zb(p8,k5)*za345b(p7,p4) -1._dp/(za(p8,b8))/(za(k6,e6))
     &    /(zb(p2,b2))/(zb(p3,k5))*za(p2,e6)*za(p3,b8)*za(p7,p8)*zb(p1,
     &    b2)*zb(p3,p8)*zb(p4,p8)*mbdk**2 +1._dp/(za(p8,b8))/(za(k6,e6))/(
     &    zb(p2,b2))/(zb(p3,k5))*za(p2,e6)*za(p3,b8)*zb(p1,b2)*zb(p3,p8
     &    )*za345b(p7,p4)*mbdk**2 )
      amp(2,1,2,1) = amp(2,1,2,1) + mb**2*prop26**(-1)*invtwoptDp8 * (
     &     -1._dp/(za(p8,b8))/(za(k6,e6))/(zb(p2,b2))*za(p2,e6)*za(p3,k5)*
     &    zb(p1,b2)*za345b(p7,p8)*za345b(b8,p4) )
      amp(2,1,2,1) = amp(2,1,2,1) + mt**2*prop26**(-1)*invtwoptDp8 * (
     &   1._dp/(za(p8,b8))/(zb(p2,b2))*za(p3,k5)*za(p7,b8)*zb(p4,p8)*zb(b2
     &    ,k6)*za26b(p2,p1) )
      amp(2,1,2,1) = amp(2,1,2,1) + mt**2*prop167**(-1)*prop58**(-1)
     &  * (1._dp/(za(p8,b8))/(zb(p2,b2))*za(p2,p7)*za(p3,k5)*za(b8,k5)*zb(
     &    p1,k6)*zb(p4,b2)*zb(p8,k5) -1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(p3
     &    ,k5))*za(p2,p7)*za(p3,b8)*zb(p1,k6)*zb(p3,p8)*zb(p4,b2)*
     &    mbdk**2 )
      amp(2,1,2,1) = amp(2,1,2,1) + mt**2*prop167**(-1)*invtwoptDp8
     &  * (1._dp/(za(p8,b8))/(zb(p2,b2))*za(p2,p7)*za(p3,k5)*zb(p1,k6)*zb(
     &    p4,p8)*za345b(b8,b2) -1._dp/(za(p8,b8))/(zb(p2,b2))*za(p2,p7)*
     &    za(p3,k5)*zb(p1,k6)*zb(b2,p8)*za345b(b8,p4) +1._dp/(za(p8,b8))/(
     &    zb(p2,b2))*za(p2,b8)*za(p3,k5)*zb(p1,k6)*zb(p4,p8)*za167b(p7,
     &    b2) +1._dp/(zb(p2,b2))*za(p2,p7)*za(p3,k5)*zb(p1,k6)*zb(p4,p8)*
     &    zb(b2,p8) )
      amp(2,1,2,1) = amp(2,1,2,1) + mt**2*mb**2*prop26**(-1)*
     & invtwoptDp8 * (1._dp/(za(p8,b8))/(za(k6,e6))/(zb(p2,b2))*za(p2,e6)*
     &    za(p3,k5)*za(p7,b8)*zb(p1,b2)*zb(p4,p8) )

      amp(2,1,1,2)= + prop26**(-1)*prop58**(-1) * (  -1._dp/(za(p3,k5))/(
     &    zb(p2,b2))*za(p3,p8)**2*za(p7,p8)*zb(p4,p8)*zb(b2,k6)*za26b(
     &    p2,p1)*mbdk +1._dp/(za(p3,k5))/(zb(p2,b2))*za(p3,p8)**2*zb(b2,k6
     &    )*za345b(p7,p4)*za26b(p2,p1)*mbdk )
      amp(2,1,1,2) = amp(2,1,1,2) + prop167**(-1)*prop58**(-1) * (  -
     &   1._dp/(za(p3,k5))/(zb(p2,b2))*za(p2,p8)*za(p3,p8)**2*zb(p1,k6)*
     &    zb(p4,p8)*za167b(p7,b2)*mbdk +1._dp/(za(p3,k5))/(zb(p2,b2))*za(
     &    p3,p8)**2*zb(p1,k6)*za167b(p7,b2)*za345b(p2,p4)*mbdk )
      amp(2,1,1,2) = amp(2,1,1,2) + mb**2*prop26**(-1)*prop58**(-1)
     &  * (  -1._dp/(za(p3,k5))/(za(k6,e6))/(zb(p2,b2))*za(p2,e6)*za(p3,p8
     &    )**2*za(p7,p8)*zb(p1,b2)*zb(p4,p8)*mbdk +1._dp/(za(p3,k5))/(za(
     &    k6,e6))/(zb(p2,b2))*za(p2,e6)*za(p3,p8)**2*zb(p1,b2)*za345b(
     &    p7,p4)*mbdk )
      amp(2,1,1,2) = amp(2,1,1,2) + mt**2*prop167**(-1)*prop58**(-1)
     &  * (  -1._dp/(za(p3,k5))/(zb(p2,b2))*za(p2,p7)*za(p3,p8)**2*zb(p1,
     &    k6)*zb(p4,b2)*mbdk )

      amp(2,1,1,1)= + prop26**(-1)*prop58**(-1) * (1._dp/(zb(p2,b2))*za(p3
     &    ,p8)*za(p7,p8)*za(p8,k5)*zb(p4,p8)*zb(b2,k6)*za26b(p2,p1) -
     &   1._dp/(zb(p2,b2))*za(p3,p8)*za(p8,k5)*zb(b2,k6)*za345b(p7,p4)*
     &    za26b(p2,p1) +1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p3,p8
     &    )*za(p7,p8)*zb(p3,b8)*zb(p4,p8)*zb(b2,k6)*za26b(p2,p1)*
     &    mbdk**2 -1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p3,p8)*zb(
     &    p3,b8)*zb(b2,k6)*za345b(p7,p4)*za26b(p2,p1)*mbdk**2 -1._dp/(zb(
     &    p2,b2))/(zb(p8,b8))*za(p3,k5)*za(p7,p8)*za(p8,k5)*zb(p4,p8)*
     &    zb(b2,k6)*zb(b8,k5)*za26b(p2,p1) +1._dp/(zb(p2,b2))/(zb(p8,b8))*
     &    za(p3,k5)*za(p8,k5)*zb(b2,k6)*zb(b8,k5)*za345b(p7,p4)*za26b(
     &    p2,p1) )
      amp(2,1,1,1) = amp(2,1,1,1) + prop26**(-1)*invtwoptDp8 * (1._dp/(zb(
     &    p2,b2))*za(p3,k5)*za(p7,p8)*zb(b2,k6)*za345b(p8,p4)*za26b(p2,
     &    p1) +1._dp/(zb(p2,b2))/(zb(p8,b8))*za(p3,k5)*zb(b2,k6)*za345b(p7
     &    ,b8)*za345b(p8,p4)*za26b(p2,p1) )
      amp(2,1,1,1) = amp(2,1,1,1) + prop167**(-1)*prop58**(-1) * (1._dp/(
     &    zb(p2,b2))*za(p2,p8)*za(p3,p8)*za(p8,k5)*zb(p1,k6)*zb(p4,p8)*
     &    za167b(p7,b2) -1._dp/(zb(p2,b2))*za(p3,p8)*za(p8,k5)*zb(p1,k6)*
     &    za167b(p7,b2)*za345b(p2,p4) +1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(
     &    p8,b8))*za(p2,p8)*za(p3,p8)*zb(p1,k6)*zb(p3,b8)*zb(p4,p8)*
     &    za167b(p7,b2)*mbdk**2 -1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))
     &    *za(p3,p8)*zb(p1,k6)*zb(p3,b8)*za167b(p7,b2)*za345b(p2,p4)*
     &    mbdk**2 -1._dp/(zb(p2,b2))/(zb(p8,b8))*za(p2,p8)*za(p3,k5)*za(p8
     &    ,k5)*zb(p1,k6)*zb(p4,p8)*zb(b8,k5)*za167b(p7,b2) +1._dp/(zb(p2,
     &    b2))/(zb(p8,b8))*za(p3,k5)*za(p8,k5)*zb(p1,k6)*zb(b8,k5)*
     &    za167b(p7,b2)*za345b(p2,p4) )
      amp(2,1,1,1) = amp(2,1,1,1) + prop167**(-1)*invtwoptDp8 * (1._dp/(
     &    zb(p2,b2))*za(p2,p8)*za(p3,k5)*zb(p1,k6)*za167b(p7,b2)*
     &    za345b(p8,p4) +1._dp/(zb(p2,b2))/(zb(p8,b8))*za(p3,k5)*zb(p1,k6)
     &    *za167b(p7,b2)*za345b(p2,b8)*za345b(p8,p4) )
      amp(2,1,1,1) = amp(2,1,1,1) + mb**2*prop26**(-1)*prop58**(-1)
     &  * (1._dp/(za(k6,e6))/(zb(p2,b2))*za(p2,e6)*za(p3,p8)*za(p7,p8)*za(
     &    p8,k5)*zb(p1,b2)*zb(p4,p8) -1._dp/(za(k6,e6))/(zb(p2,b2))*za(p2,
     &    e6)*za(p3,p8)*za(p8,k5)*zb(p1,b2)*za345b(p7,p4) +1._dp/(za(k6,e6
     &    ))/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p2,e6)*za(p3,p8)*
     &    za(p7,p8)*zb(p1,b2)*zb(p3,b8)*zb(p4,p8)*mbdk**2 -1._dp/(za(k6,e6
     &    ))/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p2,e6)*za(p3,p8)*
     &    zb(p1,b2)*zb(p3,b8)*za345b(p7,p4)*mbdk**2 -1._dp/(za(k6,e6))/(
     &    zb(p2,b2))/(zb(p8,b8))*za(p2,e6)*za(p3,k5)*za(p7,p8)*za(p8,k5
     &    )*zb(p1,b2)*zb(p4,p8)*zb(b8,k5) +1._dp/(za(k6,e6))/(zb(p2,b2))/(
     &    zb(p8,b8))*za(p2,e6)*za(p3,k5)*za(p8,k5)*zb(p1,b2)*zb(b8,k5)*
     &    za345b(p7,p4) )
      amp(2,1,1,1) = amp(2,1,1,1) + mb**2*prop26**(-1)*invtwoptDp8 * (
     &   1._dp/(za(k6,e6))/(zb(p2,b2))*za(p2,e6)*za(p3,k5)*za(p7,p8)*zb(p1
     &    ,b2)*za345b(p8,p4) +1._dp/(za(k6,e6))/(zb(p2,b2))/(zb(p8,b8))*
     &    za(p2,e6)*za(p3,k5)*zb(p1,b2)*za345b(p7,b8)*za345b(p8,p4) )
      amp(2,1,1,1) = amp(2,1,1,1) + mt**2*prop26**(-1)*invtwoptDp8 * (
     &     -1._dp/(zb(p2,b2))/(zb(p8,b8))*za(p3,k5)*za(p7,p8)*zb(p4,b8)*
     &    zb(b2,k6)*za26b(p2,p1) )
      amp(2,1,1,1) = amp(2,1,1,1) + mt**2*prop167**(-1)*prop58**(-1)
     &  * (1._dp/(zb(p2,b2))*za(p2,p7)*za(p3,p8)*za(p8,k5)*zb(p1,k6)*zb(p4
     &    ,b2) +1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p2,p7)*za(p3,
     &    p8)*zb(p1,k6)*zb(p3,b8)*zb(p4,b2)*mbdk**2 -1._dp/(zb(p2,b2))/(
     &    zb(p8,b8))*za(p2,p7)*za(p3,k5)*za(p8,k5)*zb(p1,k6)*zb(p4,b2)*
     &    zb(b8,k5) )
      amp(2,1,1,1) = amp(2,1,1,1) + mt**2*prop167**(-1)*invtwoptDp8
     &  * (  -1._dp/(zb(p2,b2))/(zb(p8,b8))*za(p2,p7)*za(p3,k5)*zb(p1,k6)*
     &    zb(p4,b8)*za345b(p8,b2) +1._dp/(zb(p2,b2))/(zb(p8,b8))*za(p2,p7)
     &    *za(p3,k5)*zb(p1,k6)*zb(b2,b8)*za345b(p8,p4) -1._dp/(zb(p2,b2))
     &    /(zb(p8,b8))*za(p2,p8)*za(p3,k5)*zb(p1,k6)*zb(p4,b8)*za167b(
     &    p7,b2) )
      amp(2,1,1,1) = amp(2,1,1,1) + mt**2*mb**2*prop26**(-1)*
     & invtwoptDp8 * (  -1._dp/(za(k6,e6))/(zb(p2,b2))/(zb(p8,b8))*za(p2,
     &    e6)*za(p3,k5)*za(p7,p8)*zb(p1,b2)*zb(p4,b8) )

      amp(1,2,2,2)=czip

      amp(1,2,2,1)= + mb*prop26**(-1)*prop58**(-1) * (1._dp/(za(p2,b2))/(
     &    za(p8,b8))*za(p3,k5)*za(b2,k6)*za(p7,p8)*za(b8,k5)*zb(p1,p2)*
     &    zb(p4,p8)*zb(p8,k5) -1._dp/(za(p2,b2))/(za(p8,b8))*za(p3,k5)*za(
     &    b2,k6)*za(b8,k5)*zb(p1,p2)*zb(p8,k5)*za345b(p7,p4) -1._dp/(za(p2
     &    ,b2))/(za(p8,b8))/(zb(p3,k5))*za(p3,b8)*za(b2,k6)*za(p7,p8)*
     &    zb(p1,p2)*zb(p3,p8)*zb(p4,p8)*mbdk**2 +1._dp/(za(p2,b2))/(za(p8,
     &    b8))/(zb(p3,k5))*za(p3,b8)*za(b2,k6)*zb(p1,p2)*zb(p3,p8)*
     &    za345b(p7,p4)*mbdk**2 -1._dp/(za(p2,b2))/(za(p8,b8))/(zb(p3,k5))
     &    /(zb(k6,e6))*za(p3,b8)*za(p7,p8)*zb(p2,e6)*zb(p3,p8)*zb(p4,p8
     &    )*za26b(b2,p1)*mbdk**2 +1._dp/(za(p2,b2))/(za(p8,b8))/(zb(p3,k5)
     &    )/(zb(k6,e6))*za(p3,b8)*zb(p2,e6)*zb(p3,p8)*za345b(p7,p4)*
     &    za26b(b2,p1)*mbdk**2 +1._dp/(za(p2,b2))/(za(p8,b8))/(zb(k6,e6))*
     &    za(p3,k5)*za(p7,p8)*za(b8,k5)*zb(p2,e6)*zb(p4,p8)*zb(p8,k5)*
     &    za26b(b2,p1) -1._dp/(za(p2,b2))/(za(p8,b8))/(zb(k6,e6))*za(p3,k5
     &    )*za(b8,k5)*zb(p2,e6)*zb(p8,k5)*za345b(p7,p4)*za26b(b2,p1) )
      amp(1,2,2,1) = amp(1,2,2,1) + mb*prop26**(-1)*invtwoptDp8 * (  -
     &   1._dp/(za(p2,b2))/(za(p8,b8))*za(p3,k5)*za(b2,k6)*zb(p1,p2)*
     &    za345b(p7,p8)*za345b(b8,p4) -1._dp/(za(p2,b2))/(za(p8,b8))/(zb(
     &    k6,e6))*za(p3,k5)*zb(p2,e6)*za345b(p7,p8)*za345b(b8,p4)*
     &    za26b(b2,p1) )
      amp(1,2,2,1) = amp(1,2,2,1) + mb*prop167**(-1)*prop58**(-1) * (
     &     -1._dp/(za(p2,b2))/(za(p8,b8))/(zb(p3,k5))/(zb(k6,e6))*za(p3,b8
     &    )*za(b2,p8)*zb(p1,e6)*zb(p3,p8)*zb(p4,p8)*za167b(p7,p2)*
     &    mbdk**2 +1._dp/(za(p2,b2))/(za(p8,b8))/(zb(p3,k5))/(zb(k6,e6))*
     &    za(p3,b8)*zb(p1,e6)*zb(p3,p8)*za167b(p7,p2)*za345b(b2,p4)*
     &    mbdk**2 +1._dp/(za(p2,b2))/(za(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(
     &    b2,p8)*za(b8,k5)*zb(p1,e6)*zb(p4,p8)*zb(p8,k5)*za167b(p7,p2)
     &     -1._dp/(za(p2,b2))/(za(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(b8,k5)*
     &    zb(p1,e6)*zb(p8,k5)*za167b(p7,p2)*za345b(b2,p4) )
      amp(1,2,2,1) = amp(1,2,2,1) + mb*prop167**(-1)*invtwoptDp8 * (
     &     -1._dp/(za(p2,b2))/(za(p8,b8))/(zb(k6,e6))*za(p3,k5)*zb(p1,e6)*
     &    za167b(p7,p2)*za345b(b2,p8)*za345b(b8,p4) )
      amp(1,2,2,1) = amp(1,2,2,1) + mt**2*mb*prop26**(-1)*invtwoptDp8
     &  * (1._dp/(za(p2,b2))/(za(p8,b8))*za(p3,k5)*za(b2,k6)*za(p7,b8)*zb(
     &    p1,p2)*zb(p4,p8) +1._dp/(za(p2,b2))/(za(p8,b8))/(zb(k6,e6))*za(
     &    p3,k5)*za(p7,b8)*zb(p2,e6)*zb(p4,p8)*za26b(b2,p1) )
      amp(1,2,2,1) = amp(1,2,2,1) + mt**2*mb*prop167**(-1)*prop58**(-1)
     &  * (1._dp/(za(p2,b2))/(za(p8,b8))/(zb(p3,k5))/(zb(k6,e6))*za(p3,b8)
     &    *za(b2,p7)*zb(p1,e6)*zb(p2,p4)*zb(p3,p8)*mbdk**2 -1._dp/(za(p2,
     &    b2))/(za(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(b2,p7)*za(b8,k5)*
     &    zb(p1,e6)*zb(p2,p4)*zb(p8,k5) )
      amp(1,2,2,1) = amp(1,2,2,1) + mt**2*mb*prop167**(-1)*invtwoptDp8
     &  * (  -1._dp/(za(p2,b2))/(za(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(b2,p7
     &    )*zb(p1,e6)*zb(p2,p8)*za345b(b8,p4) +1._dp/(za(p2,b2))/(za(p8,b8
     &    ))/(zb(k6,e6))*za(p3,k5)*za(b2,p7)*zb(p1,e6)*zb(p4,p8)*
     &    za345b(b8,p2) +1._dp/(za(p2,b2))/(za(p8,b8))/(zb(k6,e6))*za(p3,
     &    k5)*za(b2,b8)*zb(p1,e6)*zb(p4,p8)*za167b(p7,p2) +1._dp/(za(p2,b2
     &    ))/(zb(k6,e6))*za(p3,k5)*za(b2,p7)*zb(p1,e6)*zb(p2,p8)*zb(p4,
     &    p8) )

      amp(1,2,1,2)= + mb*prop26**(-1)*prop58**(-1) * (  -1._dp/(za(p2,b2))
     &    /(za(p3,k5))*za(p3,p8)**2*za(b2,k6)*za(p7,p8)*zb(p1,p2)*zb(p4
     &    ,p8)*mbdk +1._dp/(za(p2,b2))/(za(p3,k5))*za(p3,p8)**2*za(b2,k6)*
     &    zb(p1,p2)*za345b(p7,p4)*mbdk -1._dp/(za(p2,b2))/(za(p3,k5))/(zb(
     &    k6,e6))*za(p3,p8)**2*za(p7,p8)*zb(p2,e6)*zb(p4,p8)*za26b(b2,
     &    p1)*mbdk +1._dp/(za(p2,b2))/(za(p3,k5))/(zb(k6,e6))*za(p3,p8)**2
     &    *zb(p2,e6)*za345b(p7,p4)*za26b(b2,p1)*mbdk )
      amp(1,2,1,2) = amp(1,2,1,2) + mb*prop167**(-1)*prop58**(-1) * (
     &     -1._dp/(za(p2,b2))/(za(p3,k5))/(zb(k6,e6))*za(p3,p8)**2*za(b2,
     &    p8)*zb(p1,e6)*zb(p4,p8)*za167b(p7,p2)*mbdk +1._dp/(za(p2,b2))/(
     &    za(p3,k5))/(zb(k6,e6))*za(p3,p8)**2*zb(p1,e6)*za167b(p7,p2)*
     &    za345b(b2,p4)*mbdk )
      amp(1,2,1,2) = amp(1,2,1,2) + mt**2*mb*prop167**(-1)*prop58**(-1)
     &  * (1._dp/(za(p2,b2))/(za(p3,k5))/(zb(k6,e6))*za(p3,p8)**2*za(b2,p7
     &    )*zb(p1,e6)*zb(p2,p4)*mbdk )

      amp(1,2,1,1)= + mb*prop26**(-1)*prop58**(-1) * (1._dp/(za(p2,b2))*
     &    za(p3,p8)*za(b2,k6)*za(p7,p8)*za(p8,k5)*zb(p1,p2)*zb(p4,p8)
     &     -1._dp/(za(p2,b2))*za(p3,p8)*za(b2,k6)*za(p8,k5)*zb(p1,p2)*
     &    za345b(p7,p4) +1._dp/(za(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p3,
     &    p8)*za(b2,k6)*za(p7,p8)*zb(p1,p2)*zb(p3,b8)*zb(p4,p8)*mbdk**2
     &     -1._dp/(za(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p3,p8)*za(b2,k6)*
     &    zb(p1,p2)*zb(p3,b8)*za345b(p7,p4)*mbdk**2 +1._dp/(za(p2,b2))/(
     &    zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*za(p3,p8)*za(p7,p8)*zb(p2,
     &    e6)*zb(p3,b8)*zb(p4,p8)*za26b(b2,p1)*mbdk**2 -1._dp/(za(p2,b2))
     &    /(zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*za(p3,p8)*zb(p2,e6)*zb(
     &    p3,b8)*za345b(p7,p4)*za26b(b2,p1)*mbdk**2 -1._dp/(za(p2,b2))/(
     &    zb(p8,b8))*za(p3,k5)*za(b2,k6)*za(p7,p8)*za(p8,k5)*zb(p1,p2)*
     &    zb(p4,p8)*zb(b8,k5) +1._dp/(za(p2,b2))/(zb(p8,b8))*za(p3,k5)*za(
     &    b2,k6)*za(p8,k5)*zb(p1,p2)*zb(b8,k5)*za345b(p7,p4) -1._dp/(za(p2
     &    ,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(p7,p8)*za(p8,k5)*
     &    zb(p2,e6)*zb(p4,p8)*zb(b8,k5)*za26b(b2,p1) )
      amp(1,2,1,1) = amp(1,2,1,1) + mb*prop26**(-1)*prop58**(-1) * (
     &   1._dp/(za(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(p8,k5)*zb(
     &    p2,e6)*zb(b8,k5)*za345b(p7,p4)*za26b(b2,p1) +1._dp/(za(p2,b2))/(
     &    zb(k6,e6))*za(p3,p8)*za(p7,p8)*za(p8,k5)*zb(p2,e6)*zb(p4,p8)*
     &    za26b(b2,p1) -1._dp/(za(p2,b2))/(zb(k6,e6))*za(p3,p8)*za(p8,k5)*
     &    zb(p2,e6)*za345b(p7,p4)*za26b(b2,p1) )
      amp(1,2,1,1) = amp(1,2,1,1) + mb*prop26**(-1)*invtwoptDp8 * (1._dp/(
     &    za(p2,b2))*za(p3,k5)*za(b2,k6)*za(p7,p8)*zb(p1,p2)*za345b(p8,
     &    p4) +1._dp/(za(p2,b2))/(zb(p8,b8))*za(p3,k5)*za(b2,k6)*zb(p1,p2)
     &    *za345b(p7,b8)*za345b(p8,p4) +1._dp/(za(p2,b2))/(zb(p8,b8))/(zb(
     &    k6,e6))*za(p3,k5)*zb(p2,e6)*za345b(p7,b8)*za345b(p8,p4)*
     &    za26b(b2,p1) +1._dp/(za(p2,b2))/(zb(k6,e6))*za(p3,k5)*za(p7,p8)*
     &    zb(p2,e6)*za345b(p8,p4)*za26b(b2,p1) )
      amp(1,2,1,1) = amp(1,2,1,1) + mb*prop167**(-1)*prop58**(-1) * (
     &   1._dp/(za(p2,b2))/(zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*za(p3,p8)*
     &    za(b2,p8)*zb(p1,e6)*zb(p3,b8)*zb(p4,p8)*za167b(p7,p2)*mbdk**2
     &     -1._dp/(za(p2,b2))/(zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*za(p3,p8
     &    )*zb(p1,e6)*zb(p3,b8)*za167b(p7,p2)*za345b(b2,p4)*mbdk**2 -
     &   1._dp/(za(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(b2,p8)*za(
     &    p8,k5)*zb(p1,e6)*zb(p4,p8)*zb(b8,k5)*za167b(p7,p2) +1._dp/(za(p2
     &    ,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(p8,k5)*zb(p1,e6)*
     &    zb(b8,k5)*za167b(p7,p2)*za345b(b2,p4) +1._dp/(za(p2,b2))/(zb(k6,
     &    e6))*za(p3,p8)*za(b2,p8)*za(p8,k5)*zb(p1,e6)*zb(p4,p8)*
     &    za167b(p7,p2) -1._dp/(za(p2,b2))/(zb(k6,e6))*za(p3,p8)*za(p8,k5)
     &    *zb(p1,e6)*za167b(p7,p2)*za345b(b2,p4) )
      amp(1,2,1,1) = amp(1,2,1,1) + mb*prop167**(-1)*invtwoptDp8 * (
     &   1._dp/(za(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*zb(p1,e6)*
     &    za167b(p7,p2)*za345b(b2,b8)*za345b(p8,p4) +1._dp/(za(p2,b2))/(
     &    zb(k6,e6))*za(p3,k5)*za(b2,p8)*zb(p1,e6)*za167b(p7,p2)*
     &    za345b(p8,p4) )
      amp(1,2,1,1) = amp(1,2,1,1) + mt**2*mb*prop26**(-1)*invtwoptDp8
     &  * (  -1._dp/(za(p2,b2))/(zb(p8,b8))*za(p3,k5)*za(b2,k6)*za(p7,p8)*
     &    zb(p1,p2)*zb(p4,b8) -1._dp/(za(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*
     &    za(p3,k5)*za(p7,p8)*zb(p2,e6)*zb(p4,b8)*za26b(b2,p1) )
      amp(1,2,1,1) = amp(1,2,1,1) + mt**2*mb*prop167**(-1)*prop58**(-1)
     &  * (  -1._dp/(za(p2,b2))/(zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*za(p3,
     &    p8)*za(b2,p7)*zb(p1,e6)*zb(p2,p4)*zb(p3,b8)*mbdk**2 +1._dp/(za(
     &    p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(b2,p7)*za(p8,k5)
     &    *zb(p1,e6)*zb(p2,p4)*zb(b8,k5) -1._dp/(za(p2,b2))/(zb(k6,e6))*
     &    za(p3,p8)*za(b2,p7)*za(p8,k5)*zb(p1,e6)*zb(p2,p4) )
      amp(1,2,1,1) = amp(1,2,1,1) + mt**2*mb*prop167**(-1)*invtwoptDp8
     &  * (1._dp/(za(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(b2,p7)*
     &    zb(p1,e6)*zb(p2,b8)*za345b(p8,p4) -1._dp/(za(p2,b2))/(zb(p8,b8))
     &    /(zb(k6,e6))*za(p3,k5)*za(b2,p7)*zb(p1,e6)*zb(p4,b8)*za345b(
     &    p8,p2) -1._dp/(za(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(
     &    b2,p8)*zb(p1,e6)*zb(p4,b8)*za167b(p7,p2) )

      amp(1,1,2,2)=czip

      amp(1,1,2,1)= + mb*prop26**(-1)*prop58**(-1) * (  -1._dp/(za(p8,b8))
     &    /(zb(p2,b2))*za(p2,k6)*za(p3,k5)*za(p7,p8)*za(b8,k5)*zb(p1,b2
     &    )*zb(p4,p8)*zb(p8,k5) +1._dp/(za(p8,b8))/(zb(p2,b2))*za(p2,k6)*
     &    za(p3,k5)*za(b8,k5)*zb(p1,b2)*zb(p8,k5)*za345b(p7,p4) +1._dp/(
     &    za(p8,b8))/(zb(p2,b2))/(zb(p3,k5))*za(p2,k6)*za(p3,b8)*za(p7,
     &    p8)*zb(p1,b2)*zb(p3,p8)*zb(p4,p8)*mbdk**2 -1._dp/(za(p8,b8))/(
     &    zb(p2,b2))/(zb(p3,k5))*za(p2,k6)*za(p3,b8)*zb(p1,b2)*zb(p3,p8
     &    )*za345b(p7,p4)*mbdk**2 +1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(p3,k5
     &    ))/(zb(k6,e6))*za(p3,b8)*za(p7,p8)*zb(p3,p8)*zb(p4,p8)*zb(b2,
     &    e6)*za26b(p2,p1)*mbdk**2 -1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(p3,
     &    k5))/(zb(k6,e6))*za(p3,b8)*zb(p3,p8)*zb(b2,e6)*za345b(p7,p4)*
     &    za26b(p2,p1)*mbdk**2 -1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(k6,e6))*
     &    za(p3,k5)*za(p7,p8)*za(b8,k5)*zb(p4,p8)*zb(b2,e6)*zb(p8,k5)*
     &    za26b(p2,p1) +1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(k6,e6))*za(p3,k5
     &    )*za(b8,k5)*zb(b2,e6)*zb(p8,k5)*za345b(p7,p4)*za26b(p2,p1) )
      amp(1,1,2,1) = amp(1,1,2,1) + mb*prop26**(-1)*invtwoptDp8 * (1._dp/(
     &    za(p8,b8))/(zb(p2,b2))*za(p2,k6)*za(p3,k5)*zb(p1,b2)*za345b(
     &    p7,p8)*za345b(b8,p4) +1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(k6,e6))*
     &    za(p3,k5)*zb(b2,e6)*za345b(p7,p8)*za345b(b8,p4)*za26b(p2,p1)
     &     )
      amp(1,1,2,1) = amp(1,1,2,1) + mb*prop167**(-1)*prop58**(-1) * (
     &   1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(p3,k5))/(zb(k6,e6))*za(p2,p8)*
     &    za(p3,b8)*zb(p1,e6)*zb(p3,p8)*zb(p4,p8)*za167b(p7,b2)*mbdk**2
     &     -1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(p3,k5))/(zb(k6,e6))*za(p3,b8
     &    )*zb(p1,e6)*zb(p3,p8)*za167b(p7,b2)*za345b(p2,p4)*mbdk**2 -
     &   1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(k6,e6))*za(p2,p8)*za(p3,k5)*za(
     &    b8,k5)*zb(p1,e6)*zb(p4,p8)*zb(p8,k5)*za167b(p7,b2) +1._dp/(za(p8
     &    ,b8))/(zb(p2,b2))/(zb(k6,e6))*za(p3,k5)*za(b8,k5)*zb(p1,e6)*
     &    zb(p8,k5)*za167b(p7,b2)*za345b(p2,p4) )
      amp(1,1,2,1) = amp(1,1,2,1) + mb*prop167**(-1)*invtwoptDp8 * (
     &   1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(k6,e6))*za(p3,k5)*zb(p1,e6)*
     &    za167b(p7,b2)*za345b(p2,p8)*za345b(b8,p4) )
      amp(1,1,2,1) = amp(1,1,2,1) + mt**2*mb*prop26**(-1)*invtwoptDp8
     &  * (  -1._dp/(za(p8,b8))/(zb(p2,b2))*za(p2,k6)*za(p3,k5)*za(p7,b8)*
     &    zb(p1,b2)*zb(p4,p8) -1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(k6,e6))*
     &    za(p3,k5)*za(p7,b8)*zb(p4,p8)*zb(b2,e6)*za26b(p2,p1) )
      amp(1,1,2,1) = amp(1,1,2,1) + mt**2*mb*prop167**(-1)*prop58**(-1)
     &  * (1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(p3,k5))/(zb(k6,e6))*za(p2,p7)
     &    *za(p3,b8)*zb(p1,e6)*zb(p3,p8)*zb(p4,b2)*mbdk**2 -1._dp/(za(p8,
     &    b8))/(zb(p2,b2))/(zb(k6,e6))*za(p2,p7)*za(p3,k5)*za(b8,k5)*
     &    zb(p1,e6)*zb(p4,b2)*zb(p8,k5) )
      amp(1,1,2,1) = amp(1,1,2,1) + mt**2*mb*prop167**(-1)*invtwoptDp8
     &  * (  -1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(k6,e6))*za(p2,p7)*za(p3,k5
     &    )*zb(p1,e6)*zb(p4,p8)*za345b(b8,b2) +1._dp/(za(p8,b8))/(zb(p2,b2
     &    ))/(zb(k6,e6))*za(p2,p7)*za(p3,k5)*zb(p1,e6)*zb(b2,p8)*
     &    za345b(b8,p4) -1._dp/(za(p8,b8))/(zb(p2,b2))/(zb(k6,e6))*za(p2,
     &    b8)*za(p3,k5)*zb(p1,e6)*zb(p4,p8)*za167b(p7,b2) -1._dp/(zb(p2,b2
     &    ))/(zb(k6,e6))*za(p2,p7)*za(p3,k5)*zb(p1,e6)*zb(p4,p8)*zb(b2,
     &    p8) )

      amp(1,1,1,2)= + mb*prop26**(-1)*prop58**(-1) * (1._dp/(za(p3,k5))/(
     &    zb(p2,b2))*za(p2,k6)*za(p3,p8)**2*za(p7,p8)*zb(p1,b2)*zb(p4,
     &    p8)*mbdk -1._dp/(za(p3,k5))/(zb(p2,b2))*za(p2,k6)*za(p3,p8)**2*
     &    zb(p1,b2)*za345b(p7,p4)*mbdk +1._dp/(za(p3,k5))/(zb(p2,b2))/(zb(
     &    k6,e6))*za(p3,p8)**2*za(p7,p8)*zb(p4,p8)*zb(b2,e6)*za26b(p2,
     &    p1)*mbdk -1._dp/(za(p3,k5))/(zb(p2,b2))/(zb(k6,e6))*za(p3,p8)**2
     &    *zb(b2,e6)*za345b(p7,p4)*za26b(p2,p1)*mbdk )
      amp(1,1,1,2) = amp(1,1,1,2) + mb*prop167**(-1)*prop58**(-1) * (
     &   1._dp/(za(p3,k5))/(zb(p2,b2))/(zb(k6,e6))*za(p2,p8)*za(p3,p8)**2*
     &    zb(p1,e6)*zb(p4,p8)*za167b(p7,b2)*mbdk -1._dp/(za(p3,k5))/(zb(p2
     &    ,b2))/(zb(k6,e6))*za(p3,p8)**2*zb(p1,e6)*za167b(p7,b2)*
     &    za345b(p2,p4)*mbdk )
      amp(1,1,1,2) = amp(1,1,1,2) + mt**2*mb*prop167**(-1)*prop58**(-1)
     &  * (1._dp/(za(p3,k5))/(zb(p2,b2))/(zb(k6,e6))*za(p2,p7)*za(p3,p8)**
     &    2*zb(p1,e6)*zb(p4,b2)*mbdk )

      amp(1,1,1,1)= + mb*prop26**(-1)*prop58**(-1) * (  -1._dp/(zb(p2,b2))
     &    *za(p2,k6)*za(p3,p8)*za(p7,p8)*za(p8,k5)*zb(p1,b2)*zb(p4,p8)
     &     +1._dp/(zb(p2,b2))*za(p2,k6)*za(p3,p8)*za(p8,k5)*zb(p1,b2)*
     &    za345b(p7,p4) -1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p2,
     &    k6)*za(p3,p8)*za(p7,p8)*zb(p1,b2)*zb(p3,b8)*zb(p4,p8)*mbdk**2
     &     +1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))*za(p2,k6)*za(p3,p8)*
     &    zb(p1,b2)*zb(p3,b8)*za345b(p7,p4)*mbdk**2 -1._dp/(zb(p2,b2))/(
     &    zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*za(p3,p8)*za(p7,p8)*zb(p3,
     &    b8)*zb(p4,p8)*zb(b2,e6)*za26b(p2,p1)*mbdk**2 +1._dp/(zb(p2,b2))
     &    /(zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*za(p3,p8)*zb(p3,b8)*zb(
     &    b2,e6)*za345b(p7,p4)*za26b(p2,p1)*mbdk**2 +1._dp/(zb(p2,b2))/(
     &    zb(p8,b8))*za(p2,k6)*za(p3,k5)*za(p7,p8)*za(p8,k5)*zb(p1,b2)*
     &    zb(p4,p8)*zb(b8,k5) -1._dp/(zb(p2,b2))/(zb(p8,b8))*za(p2,k6)*za(
     &    p3,k5)*za(p8,k5)*zb(p1,b2)*zb(b8,k5)*za345b(p7,p4) +1._dp/(zb(p2
     &    ,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(p7,p8)*za(p8,k5)*
     &    zb(p4,p8)*zb(b2,e6)*zb(b8,k5)*za26b(p2,p1) )
      amp(1,1,1,1) = amp(1,1,1,1) + mb*prop26**(-1)*prop58**(-1) * (
     &     -1._dp/(zb(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(p8,k5)*
     &    zb(b2,e6)*zb(b8,k5)*za345b(p7,p4)*za26b(p2,p1) -1._dp/(zb(p2,b2)
     &    )/(zb(k6,e6))*za(p3,p8)*za(p7,p8)*za(p8,k5)*zb(p4,p8)*zb(b2,
     &    e6)*za26b(p2,p1) +1._dp/(zb(p2,b2))/(zb(k6,e6))*za(p3,p8)*za(p8,
     &    k5)*zb(b2,e6)*za345b(p7,p4)*za26b(p2,p1) )
      amp(1,1,1,1) = amp(1,1,1,1) + mb*prop26**(-1)*invtwoptDp8 * (  -
     &   1._dp/(zb(p2,b2))*za(p2,k6)*za(p3,k5)*za(p7,p8)*zb(p1,b2)*za345b(
     &    p8,p4) -1._dp/(zb(p2,b2))/(zb(p8,b8))*za(p2,k6)*za(p3,k5)*zb(p1,
     &    b2)*za345b(p7,b8)*za345b(p8,p4) -1._dp/(zb(p2,b2))/(zb(p8,b8))/(
     &    zb(k6,e6))*za(p3,k5)*zb(b2,e6)*za345b(p7,b8)*za345b(p8,p4)*
     &    za26b(p2,p1) -1._dp/(zb(p2,b2))/(zb(k6,e6))*za(p3,k5)*za(p7,p8)*
     &    zb(b2,e6)*za345b(p8,p4)*za26b(p2,p1) )
      amp(1,1,1,1) = amp(1,1,1,1) + mb*prop167**(-1)*prop58**(-1) * (
     &     -1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*za(p2,p8
     &    )*za(p3,p8)*zb(p1,e6)*zb(p3,b8)*zb(p4,p8)*za167b(p7,b2)*
     &    mbdk**2 +1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*
     &    za(p3,p8)*zb(p1,e6)*zb(p3,b8)*za167b(p7,b2)*za345b(p2,p4)*
     &    mbdk**2 +1._dp/(zb(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p2,p8)*za(
     &    p3,k5)*za(p8,k5)*zb(p1,e6)*zb(p4,p8)*zb(b8,k5)*za167b(p7,b2)
     &     -1._dp/(zb(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*za(p8,k5)*
     &    zb(p1,e6)*zb(b8,k5)*za167b(p7,b2)*za345b(p2,p4) -1._dp/(zb(p2,b2
     &    ))/(zb(k6,e6))*za(p2,p8)*za(p3,p8)*za(p8,k5)*zb(p1,e6)*zb(p4,
     &    p8)*za167b(p7,b2) +1._dp/(zb(p2,b2))/(zb(k6,e6))*za(p3,p8)*za(p8
     &    ,k5)*zb(p1,e6)*za167b(p7,b2)*za345b(p2,p4) )
      amp(1,1,1,1) = amp(1,1,1,1) + mb*prop167**(-1)*invtwoptDp8 * (
     &     -1._dp/(zb(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p3,k5)*zb(p1,e6)*
     &    za167b(p7,b2)*za345b(p2,b8)*za345b(p8,p4) -1._dp/(zb(p2,b2))/(
     &    zb(k6,e6))*za(p2,p8)*za(p3,k5)*zb(p1,e6)*za167b(p7,b2)*
     &    za345b(p8,p4) )
      amp(1,1,1,1) = amp(1,1,1,1) + mt**2*mb*prop26**(-1)*invtwoptDp8
     &  * (1._dp/(zb(p2,b2))/(zb(p8,b8))*za(p2,k6)*za(p3,k5)*za(p7,p8)*zb(
     &    p1,b2)*zb(p4,b8) +1._dp/(zb(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(
     &    p3,k5)*za(p7,p8)*zb(p4,b8)*zb(b2,e6)*za26b(p2,p1) )
      amp(1,1,1,1) = amp(1,1,1,1) + mt**2*mb*prop167**(-1)*prop58**(-1)
     &  * (  -1._dp/(zb(p2,b2))/(zb(p3,k5))/(zb(p8,b8))/(zb(k6,e6))*za(p2,
     &    p7)*za(p3,p8)*zb(p1,e6)*zb(p3,b8)*zb(p4,b2)*mbdk**2 +1._dp/(zb(
     &    p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p2,p7)*za(p3,k5)*za(p8,k5)
     &    *zb(p1,e6)*zb(p4,b2)*zb(b8,k5) -1._dp/(zb(p2,b2))/(zb(k6,e6))*
     &    za(p2,p7)*za(p3,p8)*za(p8,k5)*zb(p1,e6)*zb(p4,b2) )
      amp(1,1,1,1) = amp(1,1,1,1) + mt**2*mb*prop167**(-1)*invtwoptDp8
     &  * (1._dp/(zb(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p2,p7)*za(p3,k5)*
     &    zb(p1,e6)*zb(p4,b8)*za345b(p8,b2) -1._dp/(zb(p2,b2))/(zb(p8,b8))
     &    /(zb(k6,e6))*za(p2,p7)*za(p3,k5)*zb(p1,e6)*zb(b2,b8)*za345b(
     &    p8,p4) +1._dp/(zb(p2,b2))/(zb(p8,b8))/(zb(k6,e6))*za(p2,p8)*za(
     &    p3,k5)*zb(p1,e6)*zb(p4,b8)*za167b(p7,b2) )


      amp_dkqg_tbqdk_g=0._dp
      do i=1,2
      do j=1,2
      do k=1,2
      do m=1,2
      amp_dkqg_tbqdk_g=amp_dkqg_tbqdk_g+abs(amp(i,j,k,m))**2
      enddo
      enddo
      enddo
      enddo

      s17=real(za(p1,p7)*zb(p7,p1))
      s34=real(za(p3,p4)*zb(p4,p3))
c--- add missing overall factors
      amp_dkqg_tbqdk_g=amp_dkqg_tbqdk_g
     &         /((s17-wmass**2)**2)
     &         /((s34-wmass**2)**2+(wmass*wwidth)**2)
     &         /(mt*twidth)**2


      return
      end
