      subroutine tdecayWg(p,pqq,pqb,pb,pg,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*     top decay  t --> q(pqq)+qb(pb)+b(pb)+g(pg)                       *
*     with bottom and top masses (radiation from decay products of W   * 
*     in massless spinor notation                                      *
*     pqq,pqb,pb,pg are integer::s that point to                         *
*     the appropriate four-momenta in p                                *
*     pqq=quark                                                        *
*     pqb=antiquark                                                    *
*     pb=bottom quark                                                  *
*     pg=gluon                                                         *
*     q(t) is rendered massless wrt to pqb                             *
*     q(b) is rendered massless wrt to pqq                             *
*     returned m(bpol,gpol,tpol)                                       *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,altg,alb
      complex(dp):: m(2,2,2),cprop
      integer:: qq,qb,b,g,tg,nu,pb,pqq,pqb,pg
      parameter(tg=1,g=2,qq=3,qb=4,b=5)
      do nu=1,4
      q(tg,nu)=p(pqq,nu)+p(pqb,nu)+p(pb,nu)+p(pg,nu)
      q(g,nu)=p(pg,nu)
      q(qq,nu)=p(pqq,nu)
      q(qb,nu)=p(pqb,nu)
      q(b,nu)=p(pb,nu)
      enddo
      altg=mt**2/(2._dp*dot(q,tg,qb))
      alb=mb**2/(2._dp*dot(q,b,qq))
      do nu=1,4
      q(tg,nu)=q(tg,nu)-altg*q(qb,nu)
      q(b,nu)=q(b,nu)-alb*q(qq,nu)
      enddo
      call spinoru(5,q,za,zb)
      sw=s(qq,qb)+s(qb,g)+s(qq,g)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is the m(bpol,gpol,tpol)
      m(1,1,1)= + za(qq,b)*zb(qq,qb)*zb(qb,tg)/zb(qq,g)/zb(qb,g)*
     &    cprop**(-1) + za(b,g)*zb(qb,tg)/zb(qq,g)*cprop**(-1)

      m(1,2,1)= - za(qq,qb)*za(qq,b)*zb(qb,tg)/za(qq,g)/za(qb,g)*
     &    cprop**(-1) - za(qq,b)*zb(g,tg)/za(qb,g)*cprop**(-1)

      m(2,1,1)= + mb * ( za(qq,g)*zb(qb,tg)/za(qq,b)/zb(qq,g)*
     &    cprop**(-1) )

      m(2,2,1)=czip

      m(1,1,2)=czip

      m(1,2,2)= + mt * (  - za(qq,b)*zb(qb,g)/za(qb,g)/zb(qb,tg)*
     &    cprop**(-1) )

      m(2,1,2)=czip

      m(2,2,2)=czip

      return
      end
