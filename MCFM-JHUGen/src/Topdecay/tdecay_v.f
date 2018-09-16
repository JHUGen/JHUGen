c--- File written by FORM program tdecay_v.frm on Thu Mar  1 14:02:55 CST 2012
      subroutine tdecay_v(p,pqq,pqb,pb,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*     Virtual corrections to                                           *
*     top decay  t --> q(pqq)+qb(pqb)+b(pb)                            *
*     with bottom and top masses (and no radiation)                    *
*     and with overall factor of (4 pi)^e/Gamma(1-e)                   *
*     in massless spinor notation                                      *
*     pe,pnb,pc are integer::s that point to                             *
*     the appropriate four-momenta in p                                *
*     pqq=quark                                                        *
*     pqb=antiquark                                                    *
*     pb=bottom quark                                                  *
*     q(t) is rendered massless wrt to pqb                             *
*     q(b) is rendered massless wrt to pqq                             *
*     returned m(bpol,tpol)                                            *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'qcdcouple.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,alt,alb,ctm,
     & nloratiotopdecay,corr
      complex(dp):: m(2,2),cprop,
     & c0L,c0R,c1L,c1R,c1Lon2,c1Ron2,iza,izb
      integer:: qq,qb,b,t,si,pb,pqq,pqb,aa,bb
      parameter(t=1,qq=3,qb=4,b=2)
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
      corr=nloratiotopdecay(mt,mb,wmass,wwidth)/(cf*ason2pi)
      do si=1,4
      q(t,si)=p(pqq,si)+p(pqb,si)+p(pb,si)
      q(qq,si)=p(pqq,si)
      q(qb,si)=p(pqb,si)
      q(b,si)=p(pb,si)
      enddo
      alt=mt**2/(2._dp*dot(q,t,qb))
      alb=mb**2/(2._dp*dot(q,b,qq))
      do si=1,4
      q(t,si)=q(t,si)-alt*q(qb,si)
      q(b,si)=q(b,si)-alb*q(qq,si)
      enddo
      call spinoru(4,q,za,zb)
      sw=s(qq,qb)
      call coefsdkmass(sw,mt,mb,ctm,c0L,c0R,c1L,c1R)
      c0L=c0L+cplx1(ctm)-corr
      c1Lon2=c1L/2._dp
      c1Ron2=c1R/2._dp
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is m(bpol,tpol)
      m(1,1)= + cprop**(-1)*mt**(-1)*mb * (  - za(qq,b)*zb(qq,t)*zb(qb,
     &    b)*izb(qq,b)*c1Lon2 )
      m(1,1) = m(1,1) + cprop**(-1) * (  - za(qq,b)*za(qb,b)*zb(qb,b)*
     &    iza(qb,t)*c1Ron2 - za(qq,b)*zb(qb,t)*c0L )
      m(1,1) = m(1,1) + cprop**(-1)*mt*mb * (  - za(qq,qb)*zb(qq,qb)*
     &    iza(qb,t)*izb(qq,b)*c0R )

      m(1,2)= + cprop**(-1)*mt**(-1) * (  - za(qq,b)*za(b,t)*zb(qb,b)*
     &    c1Ron2 )
      m(1,2) = m(1,2) + cprop**(-1)*mb * ( za(qq,b)*zb(qq,qb)*zb(qb,b)*
     &    izb(qq,b)*izb(qb,t)*c1Lon2 + za(qq,t)*zb(qq,qb)*izb(qq,b)*c0R
     &     )

      m(2,1)= + cprop**(-1)*mt**(-1) * (  - za(qq,b)*zb(qb,b)*zb(b,t)*
     &    c1Lon2 )
      m(2,1) = m(2,1) + cprop**(-1)*mb * ( za(qq,qb)*zb(qb,b)*iza(qb,t)
     &    *c1Ron2 )
      m(2,1) = m(2,1) + cprop**(-1)*mt * ( za(qq,qb)*zb(qb,b)*iza(qb,t)
     &    *c0R )

      m(2,2)= + cprop**(-1)*mt**(-1)*mb * (  - za(qq,t)*zb(qb,b)*c1Ron2
     &     )
      m(2,2) = m(2,2) + cprop**(-1) * (  - za(qq,b)*zb(qb,b)**2*izb(qb,
     &    t)*c1Lon2 - za(qq,t)*zb(qb,b)*c0R )

      return
      end
