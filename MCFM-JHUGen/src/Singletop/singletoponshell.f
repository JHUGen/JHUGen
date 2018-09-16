      subroutine singletoponshell(q1,q2,q7,p,iswitch,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*                                                                      *
*     q(-p1)+g(-p2)--> t(p3,p4,p5)+pc(p6)+q'(p7)                       *
*                                                                      *
*     keeping polarization information for t,pc and gluon              *
*     iswitch= 0 for no gluon emission                                 *
*     iswitch=+1 for gluon emission in top decay                       *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4),
     & alt,alc,dot,s17,mt2,twoptDp2,twopcDp2
      complex(dp):: m(2,2,2),iza,izb,cprop
      integer:: p1,p2,t,p7,c,eb,si,aa,bb,iswitch,q1,q2,q7
      parameter(p1=1,p2=2,t=3,eb=4,c=5,p7=6)
C-----matrix element for p1+p2(g) -> t+c where both t and c are on shell
C-----t rendered massless wrt eb, and c rendered massless wrt p2
c--- statement functions
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
c--- end statement functions
c      write(6,*) 'iswitch',iswitch
c      call writeout(p)
c      pause
C---zero all arrays
      m(:,:,:)=czip
      do si=1,4
      q(p1,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      if (iswitch == 0) then
      q(t,si)=p(3,si)+p(4,si)+p(5,si)
      elseif (iswitch == 1) then
      q(t,si)=p(3,si)+p(4,si)+p(5,si)+p(8,si)
      endif
      q(eb,si)=p(4,si)
      q(c,si)=p(6,si)
      q(p7,si)=p(q7,si)
      enddo
      mt2=mt**2
C---- now render "t" massless wrt to vector eb
C---- now render "c" massless wrt to vector p2
      alt=mt2/(2._dp*dot(q,t,eb))
      alc=mb**2/(2._dp*dot(q,c,p2))
      s17=2._dp*dot(q,p1,p7)
      cprop=cplx2(s17-wmass**2,zip)
      twoptDp2=2._dp*dot(q,t,p2)
      twopcDp2=2._dp*dot(q,c,p2)
      do si=1,4
      q(t,si)=q(t,si)-alt*q(eb,si)
      q(c,si)=q(c,si)-alc*q(p2,si)
      enddo
      call spinoru(6,q,za,zb)
      
C----order of indices is polt,polg,polc
      m(1,1,2)= + cprop**(-1)*twoptDp2**(-1) * (  - za(t,p2)*za(t,p7)*
     &    zb(c,t)*zb(c,p1)*izb(c,p2) - za(t,p2)*za(eb,p7)*zb(c,eb)*zb(c
     &    ,p1)*izb(c,p2)*alt - za(t,p2)*za(p2,p7)*zb(c,p1) - za(p2,p7)*
     &    zb(c,eb)*zb(c,p1)*izb(c,p2)*izb(t,eb)*mt2 )

      m(1,2,2)= + cprop**(-1)*twoptDp2**(-1) * ( za(c,t)*za(t,p7)*zb(c,
     &    p1)*zb(t,p2)*iza(c,p2) + za(c,t)*za(eb,p7)*zb(c,p1)*zb(eb,p2)
     &    *iza(c,p2)*alt - za(c,p7)*zb(c,p1)*zb(eb,p2)*iza(c,p2)*izb(t,
     &    eb)*mt2 )
      m(1,2,2) = m(1,2,2) + cprop**(-1)*twopcDp2**(-1) * ( za(t,p7)*zb(
     &    c,p2)*zb(p1,p2) + za(t,p7)*zb(c,p2)*zb(p1,p2)*alc + za(t,p7)*
     &    zb(p1,p2)*iza(c,p2)*mb**2 )

      m(1,1,1)= + cprop**(-1)*twoptDp2**(-1) * (  - za(t,p2)*za(t,p7)*
     &    zb(c,t)*zb(p1,p2)*izb(c,p2)**2*mb - za(t,p2)*za(eb,p7)*zb(c,
     &    eb)*zb(p1,p2)*izb(c,p2)**2*alt*mb - za(t,p2)*za(p2,p7)*zb(p1,
     &    p2)*izb(c,p2)*mb - za(p2,p7)*zb(c,eb)*zb(p1,p2)*izb(c,p2)**2*
     &    izb(t,eb)*mt2*mb )

      m(1,2,1)= + cprop**(-1)*twoptDp2**(-1) * ( za(c,t)*za(t,p7)*zb(t,
     &    p2)*zb(p1,p2)*iza(c,p2)*izb(c,p2)*mb + za(c,t)*za(eb,p7)*zb(
     &    eb,p2)*zb(p1,p2)*iza(c,p2)*izb(c,p2)*alt*mb - za(c,p7)*zb(eb,
     &    p2)*zb(p1,p2)*iza(c,p2)*izb(c,p2)*izb(t,eb)*mt2*mb )

      
      m(2,1,2)= + cprop**(-1)*twoptDp2**(-1) * ( za(t,p7)*za(eb,p2)*zb(
     &    c,t)*zb(c,p1)*iza(t,eb)*izb(c,p2)*mt + za(eb,p2)*za(eb,p7)*
     &    zb(c,eb)*zb(c,p1)*iza(t,eb)*izb(c,p2)*alt*mt + za(eb,p2)*za(
     &    p2,p7)*zb(c,p1)*iza(t,eb)*mt + za(p2,p7)*zb(c,t)*zb(c,p1)*
     &    izb(c,p2)*mt )

      m(2,2,2)= + cprop**(-1)*twoptDp2**(-1) * (  - za(c,eb)*za(t,p7)*
     &    zb(c,p1)*zb(t,p2)*iza(c,p2)*iza(t,eb)*mt - za(c,eb)*za(eb,p7)
     &    *zb(c,p1)*zb(eb,p2)*iza(c,p2)*iza(t,eb)*alt*mt + za(c,p7)*zb(
     &    c,p1)*zb(t,p2)*iza(c,p2)*mt )
      m(2,2,2) = m(2,2,2) + cprop**(-1)*twopcDp2**(-1) * (  - za(eb,p7)
     &    *zb(c,p2)*zb(p1,p2)*iza(t,eb)*mt - za(eb,p7)*zb(c,p2)*zb(p1,
     &    p2)*iza(t,eb)*alc*mt - za(eb,p7)*zb(p1,p2)*iza(c,p2)*iza(t,eb
     &    )*mt*mb**2 )

      m(2,1,1)= + cprop**(-1)*twoptDp2**(-1) * ( za(t,p7)*za(eb,p2)*zb(
     &    c,t)*zb(p1,p2)*iza(t,eb)*izb(c,p2)**2*mt*mb + za(eb,p2)*za(eb
     &    ,p7)*zb(c,eb)*zb(p1,p2)*iza(t,eb)*izb(c,p2)**2*alt*mt*mb + 
     &    za(eb,p2)*za(p2,p7)*zb(p1,p2)*iza(t,eb)*izb(c,p2)*mt*mb + za(
     &    p2,p7)*zb(c,t)*zb(p1,p2)*izb(c,p2)**2*mt*mb )

      m(2,2,1)= + cprop**(-1)*twoptDp2**(-1) * (  - za(c,eb)*za(t,p7)*
     &    zb(t,p2)*zb(p1,p2)*iza(c,p2)*iza(t,eb)*izb(c,p2)*mt*mb - za(c
     &    ,eb)*za(eb,p7)*zb(eb,p2)*zb(p1,p2)*iza(c,p2)*iza(t,eb)*izb(c,
     &    p2)*alt*mt*mb + za(c,p7)*zb(t,p2)*zb(p1,p2)*iza(c,p2)*izb(c,
     &    p2)*mt*mb )

      return
      end
