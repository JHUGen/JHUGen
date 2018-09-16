      subroutine toppaironshell(p,iswitch,m,mab,mba)
      implicit none
      include 'types.f'
c---                                                                      
c---     q(-p1) +qbar(-p2)=nu(p3)+e+(p4)+b(p5)+bbar(p6)+e-(p7)+nubar(p8) [+g(p9]]  
c---                                                                      
c--- iswitch= 0 for no gluon emission
c--- iswitch=+1 for gluon emission in top decay
c--- iswitch=-1 for gluon emission in anti-top decay
c--- the matrix m is the q-qbar initiated process
c--- and mab,mba are the color-ordered gluon initiated process
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4),
     & alt,ala,dot,s12,mt2,twoptDp1,twoptDp2
      complex(dp):: m(2,2,2),mab(2,2,2,2),mba(2,2,2,2),iza,izb
      integer:: p1,p2,t,a,e,eb,si,aa,bb,iswitch
      parameter(p1=1,p2=2,t=3,eb=4,a=5,e=6)
C-----matrix element for p1+p2 -> t+a where both t and a are on shell
C-----t rendered massless wrt eb, and a rendered massless wrt e
c--- statement functions
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
c--- end statement functions

C---zero all arrays
      m(:,:,:)=czip
      mab(:,:,:,:)=czip
      mba(:,:,:,:)=czip

      do si=1,4
      q(p1,si)=p(1,si)
      q(p2,si)=p(2,si)
      q(eb,si)=p(4,si)
      q(e,si)=p(7,si)
      if (iswitch == 1) then
      q(t,si)=p(3,si)+p(4,si)+p(5,si)+p(9,si)
      q(a,si)=p(6,si)+p(7,si)+p(8,si)
      elseif (iswitch == -1) then
      q(t,si)=p(3,si)+p(4,si)+p(5,si)
      q(a,si)=p(6,si)+p(7,si)+p(8,si)+p(9,si)
      elseif (iswitch == 0) then
      q(t,si)=p(3,si)+p(4,si)+p(5,si)
      q(a,si)=p(6,si)+p(7,si)+p(8,si)
      endif
      enddo
      mt2=mt**2
C---- now render "t" massless wrt to vector eb 
C---- now render "a" massless wrt to vector e 
      alt=mt**2/(2._dp*dot(q,t,eb))
      ala=mt**2/(2._dp*dot(q,a,e))
      s12=2._dp*dot(q,1,2)
      twoptDp1=2._dp*dot(q,t,p1)
      twoptDp2=2._dp*dot(q,t,p2)
      do si=1,4
      q(t,si)=q(t,si)-alt*q(eb,si)
      q(a,si)=q(a,si)-ala*q(e,si)
      enddo

      call spinoru(6,q,za,zb)

c      include 'qqb.f'
C----order of indices is polt,polq1,pola
      m(1,1,1)= - za(t,p1)*zb(e,p2)*izb(e,a)*mt*s12**(-1) + za(a,p1)*
     & zb(eb,p2)*izb(t,eb)*mt*s12**(-1)

      m(1,2,1)= - za(t,p2)*zb(e,p1)*izb(e,a)*mt*s12**(-1) + za(a,p2)*
     & zb(eb,p1)*izb(t,eb)*mt*s12**(-1)

      m(1,1,2)=za(e,p1)*zb(eb,p2)*iza(e,a)*izb(t,eb)*mt**2*s12**(-1) - 
     & za(t,p1)*zb(a,p2)*s12**(-1)

      m(1,2,2)=za(e,p2)*zb(eb,p1)*iza(e,a)*izb(t,eb)*mt**2*s12**(-1) - 
     & za(t,p2)*zb(a,p1)*s12**(-1)

      m(2,1,1)= - za(a,p1)*zb(t,p2)*s12**(-1) + za(eb,p1)*zb(e,p2)*iza(
     & t,eb)*izb(e,a)*mt**2*s12**(-1)

      m(2,2,1)= - za(a,p2)*zb(t,p1)*s12**(-1) + za(eb,p2)*zb(e,p1)*iza(
     & t,eb)*izb(e,a)*mt**2*s12**(-1)

      m(2,1,2)= - za(e,p1)*zb(t,p2)*iza(e,a)*mt*s12**(-1) + za(eb,p1)*
     & zb(a,p2)*iza(t,eb)*mt*s12**(-1)

      m(2,2,2)= - za(e,p2)*zb(t,p1)*iza(e,a)*mt*s12**(-1) + za(eb,p2)*
     & zb(a,p1)*iza(t,eb)*mt*s12**(-1)

c      include 'gg.f'
c--- Ordering of polarizations, poltop * polp * polq * polantitop
      mab(1,1,1,1)= + twoptDp1**(-1) * (  - za(t,p1)*za(t,p2)*zb(t,a)*
     &    izb(a,p1)*izb(a,p2)*mt - za(t,p1)*za(a,p2)*zb(t,a)*zb(a,eb)*
     &    izb(t,eb)*izb(a,p1)*izb(a,p2)*mt + za(t,p1)*za(eb,p2)*zb(a,eb
     &    )*izb(a,p1)*izb(a,p2)*alt*mt + za(t,p1)*za(p1,p2)*izb(a,p2)*
     &    mt + za(a,p2)*za(eb,p1)*zb(a,eb)**2*izb(t,eb)*izb(a,p1)*izb(a
     &    ,p2)*alt*mt + za(p1,p2)*zb(a,eb)*izb(t,eb)*izb(a,p1)*izb(a,p2
     &    )*mt*mt2 )
      mab(1,1,1,1) = mab(1,1,1,1) - za(t,p1)*izb(a,p2)*izb(p1,p2)*mt - 
     &    za(t,p2)*izb(a,p1)*izb(p1,p2)*mt - za(a,p1)*zb(a,eb)*izb(t,eb
     &    )*izb(a,p2)*izb(p1,p2)*mt - za(a,p2)*zb(a,eb)*izb(t,eb)*izb(a
     &    ,p1)*izb(p1,p2)*mt

      mab(1,2,1,1)= + twoptDp1**(-1) * (  - za(t,p2)**2*zb(e,p1)*zb(t,
     &    p1)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*mt + za(t,p2)*za(a,p2)*zb(
     &    t,p1)*zb(eb,p1)*iza(p1,p2)*izb(t,eb)*izb(p1,p2)*mt - za(t,p2)
     &    *za(eb,p2)*zb(e,p1)*zb(eb,p1)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*
     &    alt*mt + za(a,p2)*za(eb,p2)*zb(eb,p1)**2*iza(p1,p2)*izb(t,eb)
     &    *izb(p1,p2)*alt*mt )

      mab(1,1,2,1)= + twoptDp1**(-1) * (  - za(t,p1)**2*zb(e,p2)*zb(t,
     &    p2)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*mt + za(t,p1)*za(a,p1)*zb(
     &    t,p2)*zb(eb,p2)*iza(p1,p2)*izb(t,eb)*izb(p1,p2)*mt - za(t,p1)
     &    *za(eb,p1)*zb(e,p2)*zb(eb,p2)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*
     &    alt*mt + za(a,p1)*za(eb,p1)*zb(eb,p2)**2*iza(p1,p2)*izb(t,eb)
     &    *izb(p1,p2)*alt*mt )

      mab(1,2,2,1)= + twoptDp1**(-1) * ( za(t,a)*za(t,eb)*zb(eb,p1)*zb(
     &    eb,p2)*iza(t,p1)*iza(t,p2)*izb(t,eb)*alt*mt + za(t,a)*zb(eb,
     &    p1)*zb(p1,p2)*iza(t,p2)*izb(t,eb)*mt )
      mab(1,2,2,1) = mab(1,2,2,1) - za(t,a)*zb(eb,p1)*iza(t,p2)*iza(p1,
     & p2)*izb(t,eb)*mt - za(t,a)*zb(eb,p2)*iza(t,p1)*iza(p1,p2)*izb(t,
     &    eb)*mt

      mab(2,1,1,1)= + twoptDp1**(-1) * (  - za(eb,p1)*za(eb,p2)*zb(e,t)
     &    *zb(t,eb)*iza(t,eb)*izb(e,a)*izb(t,p1)*izb(t,p2)*alt*mt2 - 
     &    za(eb,p1)*za(p1,p2)*zb(e,t)*iza(t,eb)*izb(e,a)*izb(t,p2)*mt2
     &     )
      mab(2,1,1,1) = mab(2,1,1,1) + za(eb,p1)*zb(e,t)*iza(t,eb)*izb(e,a
     & )*izb(t,p2)*izb(p1,p2)*mt2 + za(eb,p2)*zb(e,t)*iza(t,eb)*izb(e,a
     &    )*izb(t,p1)*izb(p1,p2)*mt2

      mab(2,2,1,1)= + twoptDp1**(-1) * (  - za(t,p2)*za(a,p2)*zb(t,p1)
     &    **2*iza(p1,p2)*izb(p1,p2) + za(t,p2)*za(eb,p2)*zb(e,p1)*zb(t,
     &    p1)*iza(t,eb)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*mt2 - za(a,p2)*
     &    za(eb,p2)*zb(t,p1)*zb(eb,p1)*iza(p1,p2)*izb(p1,p2)*alt + za(
     &    eb,p2)**2*zb(e,p1)*zb(eb,p1)*iza(t,eb)*iza(p1,p2)*izb(e,a)*
     &    izb(p1,p2)*alt*mt2 )

      mab(2,1,2,1)= + twoptDp1**(-1) * (  - za(t,p1)*za(a,p1)*zb(t,p2)
     &    **2*iza(p1,p2)*izb(p1,p2) + za(t,p1)*za(eb,p1)*zb(e,p2)*zb(t,
     &    p2)*iza(t,eb)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*mt2 - za(a,p1)*
     &    za(eb,p1)*zb(t,p2)*zb(eb,p2)*iza(p1,p2)*izb(p1,p2)*alt + za(
     &    eb,p1)**2*zb(e,p2)*zb(eb,p2)*iza(t,eb)*iza(p1,p2)*izb(e,a)*
     &    izb(p1,p2)*alt*mt2 )

      mab(2,2,2,1)= + twoptDp1**(-1) * (  - za(t,a)*za(a,eb)*zb(e,p2)*
     &    zb(t,p1)*iza(t,eb)*iza(a,p1)*iza(a,p2)*izb(e,a)*mt2 + za(a,eb
     &    )**2*zb(e,p2)*zb(eb,p1)*iza(t,eb)*iza(a,p1)*iza(a,p2)*izb(e,a
     &    )*alt*mt2 )
      mab(2,2,2,1) = mab(2,2,2,1) - za(a,eb)*zb(e,p1)*iza(t,eb)*iza(a,
     & p2)*iza(p1,p2)*izb(e,a)*mt2 - za(a,eb)*zb(e,p2)*iza(t,eb)*iza(a,
     &    p1)*iza(p1,p2)*izb(e,a)*mt2

      mab(1,1,1,2)= + twoptDp1**(-1) * (  - za(e,p2)*za(t,p1)*zb(t,a)*
     &    zb(a,eb)*iza(e,a)*izb(t,eb)*izb(a,p1)*izb(a,p2)*mt2 + za(e,p2
     &    )*za(eb,p1)*zb(a,eb)**2*iza(e,a)*izb(t,eb)*izb(a,p1)*izb(a,p2
     &    )*alt*mt2 )
      mab(1,1,1,2) = mab(1,1,1,2) - za(e,p1)*zb(a,eb)*iza(e,a)*izb(t,eb
     & )*izb(a,p2)*izb(p1,p2)*mt2 - za(e,p2)*zb(a,eb)*iza(e,a)*izb(t,eb
     &    )*izb(a,p1)*izb(p1,p2)*mt2

      mab(1,2,1,2)= + twoptDp1**(-1) * ( za(e,p2)*za(t,p2)*zb(t,p1)*zb(
     &    eb,p1)*iza(e,a)*iza(p1,p2)*izb(t,eb)*izb(p1,p2)*mt2 + za(e,p2
     &    )*za(eb,p2)*zb(eb,p1)**2*iza(e,a)*iza(p1,p2)*izb(t,eb)*izb(p1
     &    ,p2)*alt*mt2 - za(t,p2)**2*zb(t,p1)*zb(a,p1)*iza(p1,p2)*izb(
     &    p1,p2) - za(t,p2)*za(eb,p2)*zb(a,p1)*zb(eb,p1)*iza(p1,p2)*
     &    izb(p1,p2)*alt )

      mab(1,1,2,2)= + twoptDp1**(-1) * ( za(e,p1)*za(t,p1)*zb(t,p2)*zb(
     &    eb,p2)*iza(e,a)*iza(p1,p2)*izb(t,eb)*izb(p1,p2)*mt2 + za(e,p1
     &    )*za(eb,p1)*zb(eb,p2)**2*iza(e,a)*iza(p1,p2)*izb(t,eb)*izb(p1
     &    ,p2)*alt*mt2 - za(t,p1)**2*zb(t,p2)*zb(a,p2)*iza(p1,p2)*izb(
     &    p1,p2) - za(t,p1)*za(eb,p1)*zb(a,p2)*zb(eb,p2)*iza(p1,p2)*
     &    izb(p1,p2)*alt )

      mab(1,2,2,2)= + twoptDp1**(-1) * (  - za(e,t)*za(t,eb)*zb(eb,p1)*
     &    zb(eb,p2)*iza(e,a)*iza(t,p1)*iza(t,p2)*izb(t,eb)*alt*mt2 - 
     &    za(e,t)*zb(eb,p1)*zb(p1,p2)*iza(e,a)*iza(t,p2)*izb(t,eb)*mt2
     &     )
      mab(1,2,2,2) = mab(1,2,2,2) + za(e,t)*zb(eb,p1)*iza(e,a)*iza(t,p2
     & )*iza(p1,p2)*izb(t,eb)*mt2 + za(e,t)*zb(eb,p2)*iza(e,a)*iza(t,p1
     &    )*iza(p1,p2)*izb(t,eb)*mt2

      mab(2,1,1,2)= + twoptDp1**(-1) * (  - za(e,p2)*za(t,p1)*zb(t,a)**
     &    2*iza(e,a)*izb(a,p1)*izb(a,p2)*mt + za(e,p2)*za(eb,p1)*zb(t,a
     &    )*zb(a,eb)*iza(e,a)*izb(a,p1)*izb(a,p2)*alt*mt )
      mab(2,1,1,2) = mab(2,1,1,2) - za(e,p1)*zb(t,a)*iza(e,a)*izb(a,p2)
     & *izb(p1,p2)*mt - za(e,p2)*zb(t,a)*iza(e,a)*izb(a,p1)*izb(p1,p2)*
     &    mt

      mab(2,2,1,2)= + twoptDp1**(-1) * (  - za(e,p2)*za(t,p2)*zb(t,p1)
     &    **2*iza(e,a)*iza(p1,p2)*izb(p1,p2)*mt - za(e,p2)*za(eb,p2)*
     &    zb(t,p1)*zb(eb,p1)*iza(e,a)*iza(p1,p2)*izb(p1,p2)*alt*mt + 
     &    za(t,p2)*za(eb,p2)*zb(t,p1)*zb(a,p1)*iza(t,eb)*iza(p1,p2)*
     &    izb(p1,p2)*mt + za(eb,p2)**2*zb(a,p1)*zb(eb,p1)*iza(t,eb)*
     &    iza(p1,p2)*izb(p1,p2)*alt*mt )

      mab(2,1,2,2)= + twoptDp1**(-1) * (  - za(e,p1)*za(t,p1)*zb(t,p2)
     &    **2*iza(e,a)*iza(p1,p2)*izb(p1,p2)*mt - za(e,p1)*za(eb,p1)*
     &    zb(t,p2)*zb(eb,p2)*iza(e,a)*iza(p1,p2)*izb(p1,p2)*alt*mt + 
     &    za(t,p1)*za(eb,p1)*zb(t,p2)*zb(a,p2)*iza(t,eb)*iza(p1,p2)*
     &    izb(p1,p2)*mt + za(eb,p1)**2*zb(a,p2)*zb(eb,p2)*iza(t,eb)*
     &    iza(p1,p2)*izb(p1,p2)*alt*mt )

      mab(2,2,2,2)= + twoptDp1**(-1) * ( za(e,t)*za(t,eb)*zb(t,p1)*zb(
     &    eb,p2)*iza(e,a)*iza(t,p1)*iza(t,p2)*alt*mt + za(e,t)*zb(t,p1)
     &    *zb(p1,p2)*iza(e,a)*iza(t,p2)*mt + za(e,t)*zb(p1,p2)*iza(e,a)
     &    *iza(t,p1)*iza(t,p2)*mt*mt2 + za(t,eb)*zb(a,p2)*zb(eb,p1)*
     &    iza(t,p1)*iza(t,p2)*alt*mt )
      mab(2,2,2,2) = mab(2,2,2,2) - za(e,t)*zb(t,p1)*iza(e,a)*iza(t,p2)
     & *iza(p1,p2)*mt - za(e,t)*zb(t,p2)*iza(e,a)*iza(t,p1)*iza(p1,p2)*
     &    mt - zb(a,p1)*iza(t,p2)*iza(p1,p2)*mt - zb(a,p2)*iza(t,p1)*
     &    iza(p1,p2)*mt

      mba(1,1,1,1)= + twoptDp2**(-1) * (  - za(t,p1)*za(t,p2)*zb(t,a)*
     &    izb(a,p1)*izb(a,p2)*mt - za(t,p2)*za(a,p1)*zb(t,a)*zb(a,eb)*
     &    izb(t,eb)*izb(a,p1)*izb(a,p2)*mt + za(t,p2)*za(eb,p1)*zb(a,eb
     &    )*izb(a,p1)*izb(a,p2)*alt*mt - za(t,p2)*za(p1,p2)*izb(a,p1)*
     &    mt + za(a,p1)*za(eb,p2)*zb(a,eb)**2*izb(t,eb)*izb(a,p1)*izb(a
     &    ,p2)*alt*mt - za(p1,p2)*zb(a,eb)*izb(t,eb)*izb(a,p1)*izb(a,p2
     &    )*mt*mt2 )
      mba(1,1,1,1) = mba(1,1,1,1) + za(t,p1)*izb(a,p2)*izb(p1,p2)*mt + 
     &    za(t,p2)*izb(a,p1)*izb(p1,p2)*mt + za(a,p1)*zb(a,eb)*izb(t,eb
     &    )*izb(a,p2)*izb(p1,p2)*mt + za(a,p2)*zb(a,eb)*izb(t,eb)*izb(a
     &    ,p1)*izb(p1,p2)*mt

      mba(1,2,1,1)= + twoptDp2**(-1) * (  - za(t,p2)**2*zb(e,p1)*zb(t,
     &    p1)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*mt + za(t,p2)*za(a,p2)*zb(
     &    t,p1)*zb(eb,p1)*iza(p1,p2)*izb(t,eb)*izb(p1,p2)*mt - za(t,p2)
     &    *za(eb,p2)*zb(e,p1)*zb(eb,p1)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*
     &    alt*mt + za(a,p2)*za(eb,p2)*zb(eb,p1)**2*iza(p1,p2)*izb(t,eb)
     &    *izb(p1,p2)*alt*mt )

      mba(1,1,2,1)= + twoptDp2**(-1) * (  - za(t,p1)**2*zb(e,p2)*zb(t,
     &    p2)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*mt + za(t,p1)*za(a,p1)*zb(
     &    t,p2)*zb(eb,p2)*iza(p1,p2)*izb(t,eb)*izb(p1,p2)*mt - za(t,p1)
     &    *za(eb,p1)*zb(e,p2)*zb(eb,p2)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*
     &    alt*mt + za(a,p1)*za(eb,p1)*zb(eb,p2)**2*iza(p1,p2)*izb(t,eb)
     &    *izb(p1,p2)*alt*mt )

      mba(1,2,2,1)= + twoptDp2**(-1) * ( za(t,a)*za(t,eb)*zb(eb,p1)*zb(
     &    eb,p2)*iza(t,p1)*iza(t,p2)*izb(t,eb)*alt*mt - za(t,a)*zb(eb,
     &    p2)*zb(p1,p2)*iza(t,p1)*izb(t,eb)*mt )
      mba(1,2,2,1) = mba(1,2,2,1) + za(t,a)*zb(eb,p1)*iza(t,p2)*iza(p1,
     & p2)*izb(t,eb)*mt + za(t,a)*zb(eb,p2)*iza(t,p1)*iza(p1,p2)*izb(t,
     &    eb)*mt

      mba(2,1,1,1)= + twoptDp2**(-1) * (  - za(eb,p1)*za(eb,p2)*zb(e,t)
     &    *zb(t,eb)*iza(t,eb)*izb(e,a)*izb(t,p1)*izb(t,p2)*alt*mt2 + 
     &    za(eb,p2)*za(p1,p2)*zb(e,t)*iza(t,eb)*izb(e,a)*izb(t,p1)*mt2
     &     )
      mba(2,1,1,1) = mba(2,1,1,1) - za(eb,p1)*zb(e,t)*iza(t,eb)*izb(e,a
     & )*izb(t,p2)*izb(p1,p2)*mt2 - za(eb,p2)*zb(e,t)*iza(t,eb)*izb(e,a
     &    )*izb(t,p1)*izb(p1,p2)*mt2

      mba(2,2,1,1)= + twoptDp2**(-1) * (  - za(t,p2)*za(a,p2)*zb(t,p1)
     &    **2*iza(p1,p2)*izb(p1,p2) + za(t,p2)*za(eb,p2)*zb(e,p1)*zb(t,
     &    p1)*iza(t,eb)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*mt2 - za(a,p2)*
     &    za(eb,p2)*zb(t,p1)*zb(eb,p1)*iza(p1,p2)*izb(p1,p2)*alt + za(
     &    eb,p2)**2*zb(e,p1)*zb(eb,p1)*iza(t,eb)*iza(p1,p2)*izb(e,a)*
     &    izb(p1,p2)*alt*mt2 )

      mba(2,1,2,1)= + twoptDp2**(-1) * (  - za(t,p1)*za(a,p1)*zb(t,p2)
     &    **2*iza(p1,p2)*izb(p1,p2) + za(t,p1)*za(eb,p1)*zb(e,p2)*zb(t,
     &    p2)*iza(t,eb)*iza(p1,p2)*izb(e,a)*izb(p1,p2)*mt2 - za(a,p1)*
     &    za(eb,p1)*zb(t,p2)*zb(eb,p2)*iza(p1,p2)*izb(p1,p2)*alt + za(
     &    eb,p1)**2*zb(e,p2)*zb(eb,p2)*iza(t,eb)*iza(p1,p2)*izb(e,a)*
     &    izb(p1,p2)*alt*mt2 )

      mba(2,2,2,1)= + twoptDp2**(-1) * (  - za(t,a)*za(a,eb)*zb(e,p1)*
     &    zb(t,p2)*iza(t,eb)*iza(a,p1)*iza(a,p2)*izb(e,a)*mt2 + za(a,eb
     &    )**2*zb(e,p1)*zb(eb,p2)*iza(t,eb)*iza(a,p1)*iza(a,p2)*izb(e,a
     &    )*alt*mt2 )
      mba(2,2,2,1) = mba(2,2,2,1) + za(a,eb)*zb(e,p1)*iza(t,eb)*iza(a,
     & p2)*iza(p1,p2)*izb(e,a)*mt2 + za(a,eb)*zb(e,p2)*iza(t,eb)*iza(a,
     &    p1)*iza(p1,p2)*izb(e,a)*mt2

      mba(1,1,1,2)= + twoptDp2**(-1) * (  - za(e,p1)*za(t,p2)*zb(t,a)*
     &    zb(a,eb)*iza(e,a)*izb(t,eb)*izb(a,p1)*izb(a,p2)*mt2 + za(e,p1
     &    )*za(eb,p2)*zb(a,eb)**2*iza(e,a)*izb(t,eb)*izb(a,p1)*izb(a,p2
     &    )*alt*mt2 )
      mba(1,1,1,2) = mba(1,1,1,2) + za(e,p1)*zb(a,eb)*iza(e,a)*izb(t,eb
     & )*izb(a,p2)*izb(p1,p2)*mt2 + za(e,p2)*zb(a,eb)*iza(e,a)*izb(t,eb
     &    )*izb(a,p1)*izb(p1,p2)*mt2

      mba(1,2,1,2)= + twoptDp2**(-1) * ( za(e,p2)*za(t,p2)*zb(t,p1)*zb(
     &    eb,p1)*iza(e,a)*iza(p1,p2)*izb(t,eb)*izb(p1,p2)*mt2 + za(e,p2
     &    )*za(eb,p2)*zb(eb,p1)**2*iza(e,a)*iza(p1,p2)*izb(t,eb)*izb(p1
     &    ,p2)*alt*mt2 - za(t,p2)**2*zb(t,p1)*zb(a,p1)*iza(p1,p2)*izb(
     &    p1,p2) - za(t,p2)*za(eb,p2)*zb(a,p1)*zb(eb,p1)*iza(p1,p2)*
     &    izb(p1,p2)*alt )

      mba(1,1,2,2)= + twoptDp2**(-1) * ( za(e,p1)*za(t,p1)*zb(t,p2)*zb(
     &    eb,p2)*iza(e,a)*iza(p1,p2)*izb(t,eb)*izb(p1,p2)*mt2 + za(e,p1
     &    )*za(eb,p1)*zb(eb,p2)**2*iza(e,a)*iza(p1,p2)*izb(t,eb)*izb(p1
     &    ,p2)*alt*mt2 - za(t,p1)**2*zb(t,p2)*zb(a,p2)*iza(p1,p2)*izb(
     &    p1,p2) - za(t,p1)*za(eb,p1)*zb(a,p2)*zb(eb,p2)*iza(p1,p2)*
     &    izb(p1,p2)*alt )

      mba(1,2,2,2)= + twoptDp2**(-1) * (  - za(e,t)*za(t,eb)*zb(eb,p1)*
     &    zb(eb,p2)*iza(e,a)*iza(t,p1)*iza(t,p2)*izb(t,eb)*alt*mt2 + 
     &    za(e,t)*zb(eb,p2)*zb(p1,p2)*iza(e,a)*iza(t,p1)*izb(t,eb)*mt2
     &     )
      mba(1,2,2,2) = mba(1,2,2,2) - za(e,t)*zb(eb,p1)*iza(e,a)*iza(t,p2
     & )*iza(p1,p2)*izb(t,eb)*mt2 - za(e,t)*zb(eb,p2)*iza(e,a)*iza(t,p1
     &    )*iza(p1,p2)*izb(t,eb)*mt2

      mba(2,1,1,2)= + twoptDp2**(-1) * (  - za(e,p1)*za(t,p2)*zb(t,a)**
     &    2*iza(e,a)*izb(a,p1)*izb(a,p2)*mt + za(e,p1)*za(eb,p2)*zb(t,a
     &    )*zb(a,eb)*iza(e,a)*izb(a,p1)*izb(a,p2)*alt*mt )
      mba(2,1,1,2) = mba(2,1,1,2) + za(e,p1)*zb(t,a)*iza(e,a)*izb(a,p2)
     & *izb(p1,p2)*mt + za(e,p2)*zb(t,a)*iza(e,a)*izb(a,p1)*izb(p1,p2)*
     &    mt

      mba(2,2,1,2)= + twoptDp2**(-1) * (  - za(e,p2)*za(t,p2)*zb(t,p1)
     &    **2*iza(e,a)*iza(p1,p2)*izb(p1,p2)*mt - za(e,p2)*za(eb,p2)*
     &    zb(t,p1)*zb(eb,p1)*iza(e,a)*iza(p1,p2)*izb(p1,p2)*alt*mt + 
     &    za(t,p2)*za(eb,p2)*zb(t,p1)*zb(a,p1)*iza(t,eb)*iza(p1,p2)*
     &    izb(p1,p2)*mt + za(eb,p2)**2*zb(a,p1)*zb(eb,p1)*iza(t,eb)*
     &    iza(p1,p2)*izb(p1,p2)*alt*mt )

      mba(2,1,2,2)= + twoptDp2**(-1) * (  - za(e,p1)*za(t,p1)*zb(t,p2)
     &    **2*iza(e,a)*iza(p1,p2)*izb(p1,p2)*mt - za(e,p1)*za(eb,p1)*
     &    zb(t,p2)*zb(eb,p2)*iza(e,a)*iza(p1,p2)*izb(p1,p2)*alt*mt + 
     &    za(t,p1)*za(eb,p1)*zb(t,p2)*zb(a,p2)*iza(t,eb)*iza(p1,p2)*
     &    izb(p1,p2)*mt + za(eb,p1)**2*zb(a,p2)*zb(eb,p2)*iza(t,eb)*
     &    iza(p1,p2)*izb(p1,p2)*alt*mt )

      mba(2,2,2,2)= + twoptDp2**(-1) * ( za(e,t)*za(t,eb)*zb(t,p2)*zb(
     &    eb,p1)*iza(e,a)*iza(t,p1)*iza(t,p2)*alt*mt - za(e,t)*zb(t,p2)
     &    *zb(p1,p2)*iza(e,a)*iza(t,p1)*mt - za(e,t)*zb(p1,p2)*iza(e,a)
     &    *iza(t,p1)*iza(t,p2)*mt*mt2 + za(t,eb)*zb(a,p1)*zb(eb,p2)*
     &    iza(t,p1)*iza(t,p2)*alt*mt )
      mba(2,2,2,2) = mba(2,2,2,2) + za(e,t)*zb(t,p1)*iza(e,a)*iza(t,p2)
     & *iza(p1,p2)*mt + za(e,t)*zb(t,p2)*iza(e,a)*iza(t,p1)*iza(p1,p2)*
     &    mt + zb(a,p1)*iza(t,p2)*iza(p1,p2)*mt + zb(a,p2)*iza(t,p1)*
     &    iza(p1,p2)*mt


      return
      end
