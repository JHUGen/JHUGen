      function qg_tbqndk_ampanti(p1,p2,p3,p4,p7,
     & k6,e6,twonDpt,twonDp6,zanb)
      implicit none
      include 'types.f'
      real(dp):: qg_tbqndk_ampanti
      
c---- 4-flavour gvec routine
c---- averaged over initial colours and spins
c---- line is contracted with the vector n(mu)
c---- and routine allows for massive vectors p5,p6
c     q(-p1)+g(-p2)--> t(p3,p4,p5)+bb(p6)+q'(p7)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'masses.f'
      integer:: p1,p2,p3,p4,p7,k6,e6,i1,i2
      complex(dp):: resp,resm,za345b,za167b,zanb(mxpart,mxpart)
      real(dp):: c6,prop26,prop167,s167,twonDpt,twonDp6

c--- statement functions
      za167b(i1,i2)=za(i1,p1)*zb(p1,i2)+za(i1,p7)*zb(p7,i2)
     &             +za(i1,k6)*zb(k6,i2)+za(i1,e6)*zb(e6,i2)*c6
      za345b(i1,i2)=-za167b(i1,i2)-za(i1,p2)*zb(p2,i2)

      c6=mb**2/real(za(k6,e6)*zb(e6,k6))

      s167=real(za(p1,p7)*zb(p7,p1)
     &         +za(p1,k6)*zb(k6,p1)+c6*za(p1,e6)*zb(e6,p1)
     &         +za(p7,k6)*zb(k6,p7)+c6*za(p7,e6)*zb(e6,p7)
     &         +cplx1(mb**2))
      prop26=real(za(p2,k6)*zb(k6,p2)+c6*za(p2,e6)*zb(e6,p2))
      prop167=s167-mt**2

      resp= + prop26**(-1) * (  - 1/(za(k6,e6))*za(p2,p7)*za345b(p3,p1)
     &    *zanb(e6,p2)*mb + 1/(za(k6,e6))*za(p7,e6)*za345b(p3,p1)*mb*
     &    twonDp6 )
      resp = resp + prop167**(-1) * ( 1/(za(k6,e6))*za(p7,e6)*za167b(p2
     &    ,p1)*zanb(p3,p2)*mb + 1/(za(k6,e6))*za(p7,e6)*za167b(p3,p1)*
     &    mb*twonDpt + 1/(za(k6,e6))*za(p7,e6)*zanb(p3,p1)*mb*s167 - 
     &    1/(za(k6,e6))*za(p7,e6)*zanb(p3,p1)*mt**2*mb )

      resm= + prop26**(-1) * ( za(p2,p7)*za345b(p3,p1)*zanb(k6,p2) - 
     &    za(p7,k6)*za345b(p3,p1)*twonDp6 - za(k6,e6)*za345b(p3,p1)*
     &    zanb(p7,e6)*c6 - 1/(zb(k6,e6))*za345b(p3,p1)*zanb(p7,e6)*
     &    mb**2 )
      resm = resm + prop167**(-1) * (  - za(p7,k6)*za167b(p2,p1)*zanb(
     &    p3,p2) - za(p7,k6)*za167b(p3,p1)*twonDpt - za(p7,k6)*zanb(p3,
     &    p1)*s167 + za(p7,k6)*zanb(p3,p1)*mt**2 )


      qg_tbqndk_ampanti=abs(resp)**2+abs(resm)**2

      return
      end
