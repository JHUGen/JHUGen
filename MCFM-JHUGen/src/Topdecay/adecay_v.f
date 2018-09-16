c--- File written by FORM program adecay_v.frm on Thu Mar  1 14:02:57 CST 2012
      subroutine adecay_v(p,pqq,pqb,pc,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*     Virtual corrections to                                           *
*     antitop decay  a --> qq(pqq)+qb(pqb)+bbar(pc)                    *
*     with bottom and top masses (and no radiation)                    *
*     in massless spinor notation                                      *
*     pqq,pqb,pc are integer::s that point to                            *
*     the appropriate four-momenta in p                                *
*     pqq=electron                                                     *
*     pqb=antineutrino                                                 *
*     pc=anti-bottom quark                                             *
*     q(c) is rendered massless wrt to pqb                             *
*     q(a) is rendered massless wrt to pqq                             *
*     returned m(apol,cpol)                                            *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'qcdcouple.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,ala,alc,ctm,
     & nloratiotopdecay,corr
      complex(dp):: m(2,2),cprop,
     & c0L,C0R,C1L,C1R,C1Lon2,C1Ron2,iza,izb
      integer:: qb,a,c,qq,si,pc,pqq,pqb,aa,bb
      parameter(a=1,qq=3,qb=4,c=2)
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)

      corr=nloratiotopdecay(mt,mb,wmass,wwidth)/(cf*ason2pi)
      do si=1,4
      q(a,si)=p(pqq,si)+p(pqb,si)+p(pc,si)
      q(qq,si)=p(pqq,si)
      q(qb,si)=p(pqb,si)
      q(c,si)=p(pc,si)
      enddo
      ala=mt**2/(2._dp*dot(q,a,qq))
      alc=mb**2/(2._dp*dot(q,c,qb))
      do si=1,4
      q(a,si)=q(a,si)-ala*q(qq,si)
      q(c,si)=q(c,si)-alc*q(qb,si)
      enddo
      call spinoru(4,q,za,zb)
      sw=s(qq,qb)
      call coefsdkmass(sw,mt,mb,ctm,c0L,c0R,c1L,C1R)
      c0L=c0L+cplx1(ctm)-corr
      c1Lon2=c1L/2._dp
      c1Ron2=c1R/2._dp
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is m(apol,cpol)
      m(1,1)= + cprop**(-1)*mt**(-1)*mb * ( za(qq,c)*zb(qb,a)*c1Ron2 )
      m(1,1) = m(1,1) + cprop**(-1) * ( za(qq,c)**2*zb(qb,c)*iza(qq,a)*
     &    c1Lon2 + za(qq,c)*zb(qb,a)*c0R )

      m(1,2)= + cprop**(-1)*mt**(-1) * ( za(qq,c)*zb(qb,c)*zb(c,a)*
     &    c1Ron2 )
      m(1,2) = m(1,2) + cprop**(-1)*mb * ( za(qq,qb)*za(qq,c)*zb(qb,c)*
     &    iza(qq,a)*iza(qb,c)*c1Lon2 + za(qq,qb)*zb(qb,a)*iza(qb,c)*c0R
     &     )

      m(2,1)= + cprop**(-1)*mt**(-1) * ( za(qq,c)*za(c,a)*zb(qb,c)*
     &    c1Lon2 )
      m(2,1) = m(2,1) + cprop**(-1)*mb * ( za(qq,c)*zb(qq,qb)*izb(qq,a)
     &    *c1Ron2 )
      m(2,1) = m(2,1) + cprop**(-1)*mt * ( za(qq,c)*zb(qq,qb)*izb(qq,a)
     &    *c0R )

      m(2,2)= + cprop**(-1)*mt**(-1)*mb * ( za(qq,c)*za(qb,a)*zb(qb,c)*
     &    iza(qb,c)*c1Lon2 )
      m(2,2) = m(2,2) + cprop**(-1) * ( za(qq,c)*zb(qq,c)*zb(qb,c)*izb(
     &    qq,a)*c1Ron2 + za(qq,a)*zb(qb,c)*c0L )
      m(2,2) = m(2,2) + cprop**(-1)*mt*mb * ( za(qq,qb)*zb(qq,qb)*iza(
     &    qb,c)*izb(qq,a)*c0R )

      return
      end
