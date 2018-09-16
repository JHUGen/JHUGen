      subroutine Wtoponshell(q1,q2,p,iswitch,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, May 2012                                     *
*                                                                      *
*     b(-c)+g(-p2)--> W^-(e-(e)+nu~(nb))+t(p5,p6,p7)                   *
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
     & alt,dot,mt2,twoptDp2,s34
      complex(dp):: m(2,2,2),iza,izb,cprop
      integer:: p2,t,c,e,nb,eb,si,aa,bb,iswitch,q1,q2
      parameter(c=1,p2=2,e=3,nb=4,t=5,eb=6)
C-----matrix element for b(p1)+p2(g) -> t+W where t is on shell
C-----t rendered massless wrt eb, and c massless
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
      q(c,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      q(e,si)=p(3,si)
      q(nb,si)=p(4,si)
      if (iswitch == 0) then
      q(t,si)=p(5,si)+p(6,si)+p(7,si)
      elseif (iswitch == 1) then
      q(t,si)=p(5,si)+p(6,si)+p(7,si)+p(8,si)
      endif
      q(eb,si)=p(6,si)
      enddo
      mt2=mt**2
C---- now render "t" massless wrt to vector eb
      alt=mt2/(two*dot(q,t,eb))
      s34=two*dot(q,e,nb)
      cprop=cplx2(s34-wmass**2,wmass*wwidth)
      twoptDp2=two*dot(q,t,p2)
      do si=1,4
      q(t,si)=q(t,si)-alt*q(eb,si)
      enddo
      call spinoru(6,q,za,zb)
      
C----order of indices is polt,polg,polc
      m(1,1,2)=czip

      m(1,2,2)=czip

      m(1,1,1)= + cprop**(-1)*twoptDp2**(-1) * ( za(e,t)*za(t,p2)*zb(c,
     &    t)*zb(c,nb)*izb(c,p2) + za(e,eb)*za(t,p2)*zb(c,nb)*zb(c,eb)*
     &    izb(c,p2)*alt + za(e,p2)*za(t,p2)*zb(c,nb) + za(e,p2)*zb(c,nb
     &    )*zb(c,eb)*izb(c,p2)*izb(t,eb)*mt2 )

      m(1,2,1)= + cprop**(-1)*twoptDp2**(-1) * ( za(e,c)*zb(c,nb)*zb(eb
     &    ,p2)*iza(c,p2)*izb(t,eb)*mt2 - za(e,t)*za(c,t)*zb(c,nb)*zb(t,
     &    p2)*iza(c,p2) - za(e,eb)*za(c,t)*zb(c,nb)*zb(eb,p2)*iza(c,p2)
     &    *alt )
      m(1,2,1) = m(1,2,1) + cprop**(-1) * ( za(e,t)*zb(nb,p2)*iza(c,p2)
     &     )

      
      m(2,1,2)=czip

      m(2,2,2)=czip

      m(2,1,1)= + cprop**(-1)*twoptDp2**(-1) * (  - za(e,t)*za(eb,p2)*
     &    zb(c,t)*zb(c,nb)*iza(t,eb)*izb(c,p2)*mt - za(e,eb)*za(eb,p2)*
     &    zb(c,nb)*zb(c,eb)*iza(t,eb)*izb(c,p2)*alt*mt - za(e,p2)*za(eb
     &    ,p2)*zb(c,nb)*iza(t,eb)*mt - za(e,p2)*zb(c,t)*zb(c,nb)*izb(c,
     &    p2)*mt )

      m(2,2,1)= + cprop**(-1)*twoptDp2**(-1) * (  - za(e,c)*zb(c,nb)*
     &    zb(t,p2)*iza(c,p2)*mt + za(e,t)*za(c,eb)*zb(c,nb)*zb(t,p2)*
     &    iza(c,p2)*iza(t,eb)*mt + za(e,eb)*za(c,eb)*zb(c,nb)*zb(eb,p2)
     &    *iza(c,p2)*iza(t,eb)*alt*mt )
      m(2,2,1) = m(2,2,1) + cprop**(-1) * (  - za(e,eb)*zb(nb,p2)*iza(c
     &    ,p2)*iza(t,eb)*mt )

      return
      end
