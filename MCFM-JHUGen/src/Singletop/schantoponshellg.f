      subroutine schantoponshellg(q1,q2,q7,p,iswitch,mi,mf)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*                                                                      *
*     u(-p1)+d~(-p2)--> t(nu,eb,pb)+pc~(p6)+g(p7)                      *
*                                                                      *
*     amplitudes with radiation in initial state mi, final state mf    *
*      (returned separately since they do not interfere)               *
*                                                                      *
*     keeping polarization information for t                           *
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
     & alt,alc,dot,s12,s127,mt2,twopcDg,twoptDg
      complex(dp):: mi(2,2,2),mf(2,2,2),iza,izb,cprop,cprop7
      integer:: p1,p2,t,eb,c,g,si,i,j,k,iswitch,q1,q2,q7
      parameter(p1=1,p2=2,t=3,eb=4,c=5,g=6)
C-----matrix element for p1+p2(g) -> t+c where both t and b are on shell
C-----t rendered massless wrt eb, and p rendered massless wrt p2
c--- statement functions
      iza(i,j)=cone/za(i,j)
      izb(i,j)=cone/zb(i,j)
c--- end statement functions
c      write(6,*) 'iswitch',iswitch
c      call writeout(p)
c      pause
C---zero all arrays
      do i=1,2
      do j=1,2
      do k=1,2
      mi(i,j,k)=czip
      mf(i,j,k)=czip
      enddo
      enddo
      enddo
      do si=1,4
      q(p1,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      if (iswitch == 0) then
      q(t,si)=p(3,si)+p(4,si)+p(5,si)
      elseif (iswitch == 1) then
      q(t,si)=p(3,si)+p(4,si)+p(5,si)+p(7,si)
      endif
      q(eb,si)=p(4,si)
      q(c,si)=p(6,si)
      q(g,si)=p(q7,si)
      enddo
      mt2=mt**2
      s12=2._dp*dot(q,p1,p2)
      s127=s12+2._dp*dot(q,p1,g)+2._dp*dot(q,p2,g)
      twoptDg=2._dp*dot(q,t,g)
      twopcDg=2._dp*dot(q,c,g)
      cprop=cplx2(s12-wmass**2,wmass*wwidth)
      cprop7=cplx2(s127-wmass**2,wmass*wwidth)
C---- now render "t" massless wrt to vector eb
C---- now render "pb" massless wrt to vector p2
      alt=mt2/(2._dp*dot(q,t,eb))
      alc=mb**2/(2._dp*dot(q,c,p1))
      do si=1,4
      q(t,si)=q(t,si)-alt*q(eb,si)
      q(c,si)=q(c,si)-alc*q(p1,si)
      enddo
      call spinoru(6,q,za,zb)
c      call writeout(q)
c      pause
      
C----order of indices is polt,polg,polc
      mi(1,1,2)= + cprop7**(-1) * ( za(t,g)*zb(c,p1)*izb(g,p2) - za(t,
     &    p2)*zb(c,p1)*zb(p1,p2)*izb(g,p1)*izb(g,p2) )

      mi(1,1,1)= czip

      mi(1,2,2)= + cprop7**(-1) * ( za(t,p2)*za(p1,p2)*zb(c,p1)*iza(g,
     &    p1)*iza(g,p2) + za(t,p2)*zb(c,g)*iza(g,p1) )

      mi(1,2,1)= + cprop7**(-1) * ( za(t,p2)*zb(g,p1)*iza(g,p1)*izb(c,
     &    p1)*mb )

      
      mf(1,1,2)= + cprop**(-1) * (  - za(c,g)*za(t,p2)*zb(c,p1)*zb(c,p2
     &    )*izb(g,p2)*twopcDg**(-1) + za(t,g)*za(t,p2)*zb(c,p1)*zb(t,p2
     &    )*izb(g,p2)*twoptDg**(-1) + za(t,g)*za(g,p2)*zb(c,p1)*
     &    twoptDg**(-1) + za(t,g)*za(eb,p2)*zb(c,p1)*zb(eb,p2)*izb(g,p2
     &    )*alt*twoptDg**(-1) - za(t,p2)*za(g,p1)*zb(p1,p2)*iza(c,p1)*
     &    izb(g,p2)*twopcDg**(-1)*mb**2 + za(g,p2)*zb(c,p1)*zb(eb,p2)*
     &    izb(t,eb)*izb(g,p2)*mt2*twoptDg**(-1) )

      mf(1,1,1)= czip

      mf(1,2,2)= + cprop**(-1) * ( za(c,p1)*za(t,p2)*zb(c,g)*zb(c,p1)*
     &    iza(g,p1)*twopcDg**(-1) - za(t,p1)*za(t,p2)*zb(c,p1)*zb(t,g)*
     &    iza(g,p1)*twoptDg**(-1) + za(t,p1)*za(eb,p2)*zb(c,p1)*zb(g,eb
     &    )*iza(g,p1)*alt*twoptDg**(-1) + za(t,p2)*zb(c,g)*zb(g,p1)*
     &    twopcDg**(-1) + za(p1,p2)*zb(c,p1)*zb(g,eb)*iza(g,p1)*izb(t,
     &    eb)*mt2*twoptDg**(-1) )

      mf(1,2,1)= + cprop**(-1) * ( za(t,p2)*zb(g,p1)**2*izb(c,p1)*
     &    twopcDg**(-1)*mb )

      
      mi(2,1,2)= + cprop7**(-1) * ( za(g,eb)*zb(c,p1)*iza(t,eb)*izb(g,
     &    p2)*mt + za(eb,p2)*zb(c,p1)*zb(p1,p2)*iza(t,eb)*izb(g,p1)*
     &    izb(g,p2)*mt )

      mi(2,1,1)= czip

      mi(2,2,2)= + cprop7**(-1) * (  - za(eb,p2)*za(p1,p2)*zb(c,p1)*
     &    iza(t,eb)*iza(g,p1)*iza(g,p2)*mt - za(eb,p2)*zb(c,g)*iza(t,eb
     &    )*iza(g,p1)*mt )

      mi(2,2,1)= + cprop7**(-1) * (  - za(eb,p2)*zb(g,p1)*iza(t,eb)*
     &    iza(g,p1)*izb(c,p1)*mt*mb )

      
      mf(2,1,2)= + cprop**(-1) * ( za(c,g)*za(eb,p2)*zb(c,p1)*zb(c,p2)*
     &    iza(t,eb)*izb(g,p2)*mt*twopcDg**(-1) + za(t,p2)*za(g,eb)*zb(c
     &    ,p1)*zb(t,p2)*iza(t,eb)*izb(g,p2)*mt*twoptDg**(-1) + za(g,eb)
     &    *za(g,p2)*zb(c,p1)*iza(t,eb)*mt*twoptDg**(-1) + za(g,eb)*za(
     &    eb,p2)*zb(c,p1)*zb(eb,p2)*iza(t,eb)*izb(g,p2)*alt*mt*
     &    twoptDg**(-1) + za(g,p1)*za(eb,p2)*zb(p1,p2)*iza(c,p1)*iza(t,
     &    eb)*izb(g,p2)*mt*twopcDg**(-1)*mb**2 - za(g,p2)*zb(c,p1)*zb(t
     &    ,p2)*izb(g,p2)*mt*twoptDg**(-1) )

      mf(2,1,1)= czip

      mf(2,2,2)= + cprop**(-1) * (  - za(c,p1)*za(eb,p2)*zb(c,g)*zb(c,
     &    p1)*iza(t,eb)*iza(g,p1)*mt*twopcDg**(-1) + za(t,p2)*za(eb,p1)
     &    *zb(c,p1)*zb(t,g)*iza(t,eb)*iza(g,p1)*mt*twoptDg**(-1) - za(
     &    eb,p1)*za(eb,p2)*zb(c,p1)*zb(g,eb)*iza(t,eb)*iza(g,p1)*alt*mt
     &    *twoptDg**(-1) - za(eb,p2)*zb(c,g)*zb(g,p1)*iza(t,eb)*mt*
     &    twopcDg**(-1) + za(p1,p2)*zb(c,p1)*zb(t,g)*iza(g,p1)*mt*
     &    twoptDg**(-1) )

      mf(2,2,1)= + cprop**(-1) * (  - za(eb,p2)*zb(g,p1)**2*iza(t,eb)*
     &    izb(c,p1)*mt*twopcDg**(-1)*mb )

      
      return
      end
