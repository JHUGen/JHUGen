      subroutine adecayg(p,pqq,pqb,pc,pg,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*     antitop decay  a --> e(pe)+nb(pnb)+bbar(pc)+g(pc)                *
*     with bottom and top masses (and radiation from t and b lines)    *
*     in massless spinor notation                                      *
*     pqq,pqb,pc are integer::s that point to                            *
*     the appropriate four-momenta in p                                *
*     pqq=quark                                                        *
*     pqb=antiquark                                                    *
*     pc=anti-bottom quark                                             *
*     pg points to the radiated gluon                                  *
*     q(c) is rendered massless wrt to pqb                             *
*     q(a) is rendered massless wrt to pqq                             *
*     returned m(apol,gpol,cpol)                                       *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,twopcDg,twopaDg,
     & alag,alc
      complex(dp):: m(2,2,2),cprop
      integer:: qq,qb,c,g,ag,nu,pc,pqq,pqb,pg
      parameter(ag=1,g=2,qq=3,qb=4,c=5)

      do nu=1,4
      q(ag,nu)=p(pqq,nu)+p(pqb,nu)+p(pc,nu)+p(pg,nu)
      q(g,nu)=p(pg,nu)
      q(qq,nu)=p(pqq,nu)
      q(qb,nu)=p(pqb,nu)
      q(c,nu)=p(pc,nu)
      enddo
      twopaDg=2._dp*dot(q,ag,g)
      twopcDg=2._dp*dot(q,c,g)
      alag=mt**2/(2._dp*dot(q,ag,qq))
      alc=mb**2/(2._dp*dot(q,c,qb))
      do nu=1,4
      q(ag,nu)=q(ag,nu)-alag*q(qq,nu)
      q(c,nu)=q(c,nu)-alc*q(qb,nu)
      enddo
      call spinoru(5,q,za,zb)
      sw=s(qq,qb)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is the m(apol,gpol,cpol)
      m(1,1,1)=czip

      m(1,2,1)=czip

      m(2,1,1)=czip

      m(2,2,1)= + cprop**(-1)*mb * ( za(qq,ag)*zb(qb,g)**2/zb(qb,c)*
     &    twopcDg**(-1) )

      m(1,1,2)= + cprop**(-1)*mt * ( za(qq,g)**2*zb(qb,c)/za(qq,ag)*
     &    twopaDg**(-1) )

      m(1,2,2)=czip

      m(2,1,2)= + cprop**(-1) * ( za(qq,g)*za(g,ag)*zb(qb,c)*
     &    twopaDg**(-1) - za(qq,ag)*za(g,ag)*zb(qb,c)*zb(c,ag)/zb(c,g)
     &    *twopaDg**(-1) )
      m(2,1,2) = m(2,1,2) + cprop**(-1)*mb**2 * (  - za(qq,ag)*za(qb,g)
     &    *zb(qb,c)/za(qb,c)/zb(c,g)*twopcDg**(-1) )
      m(2,1,2) = m(2,1,2) + cprop**(-1)*mt**2 * ( za(qq,g)*zb(qq,c)*zb(
     &    qb,c)/zb(qq,ag)/zb(c,g)*twopaDg**(-1) )

      m(2,2,2)= + cprop**(-1) * ( za(qq,ag)*za(c,ag)*zb(qb,c)*zb(g,ag)
     &    /za(c,g)*twopaDg**(-1) + za(qq,ag)*zb(qb,g)*zb(c,g)*
     &    twopcDg**(-1) )
      m(2,2,2) = m(2,2,2) + cprop**(-1)*mb**2 * ( za(qq,ag)*zb(qb,g)
     &    /za(c,g)*twopcDg**(-1) )
      m(2,2,2) = m(2,2,2) + cprop**(-1)*mt**2 * (  - za(qq,c)*zb(qq,g)*
     &    zb(qb,c)/za(c,g)/zb(qq,ag)*twopaDg**(-1) )

      return
      end
