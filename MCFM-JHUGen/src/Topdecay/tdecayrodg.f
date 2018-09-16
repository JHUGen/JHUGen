c--- File written by FORM program tdecayrodg.frm on Thu May 24 10:26:49 CDT 2012
      subroutine tdecayrodg(p,pnu,peb,pb,pbb,pem,pnb,pg,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, May 2012                                     *
*     top decay  t --> q(pnu)+qb(pb)+b(pb)+g(pg)                       *
*     with bottom and top masses (and radiation form t-b line)         *
*     in massless spinor notation                                      *
*     pnu,peb,pb,pg are integer::s that point to                         *
*     the appropriate four-momenta in p                                *
*     pnu=quark                                                        *
*     peb=antiquark                                                    *
*     pb=bottom quark                                                  *
*     pbb=anti-bottom quark                                            *
*     pem=e^-                                                          *
*     pnub=antineutrino                                                *
*     pg=gluon                                                         *
*     q(t) is rendered massless wrt to peb                             *
*     q(b) is rendered massless wrt to pg                              *
*     returned m(bpol,gpol,tpol)                                       *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,twoptDg,
     & alb,a(4),t(4),tpa(4),s34,betasq,be,bp,bm,rtbp
      complex(dp):: m(2,2,2),cprop,iza,izb
      integer:: qq,qb,b,g,pb,pnu,peb,pbb,pem,pnb,pg,aa,bb,si,k5,k6
      parameter(g=1,qq=2,qb=3,b=4,k5=5,k6=6)
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
C construct top and antitop momenta
      do si=1,4
      t(si)=p(pnu,si)+p(peb,si)+p(pb,si)+p(pg,si)
      a(si)=p(pem,si)+p(pnb,si)+p(pbb,si)
      tpa(si)=t(si)+a(si)
      enddo
      twoptDg=2._dp*(t(4)*p(pg,4)-t(1)*p(pg,1)-t(2)*p(pg,2)-t(3)*p(pg,3))
      alb=mb**2/(2._dp*dot(p,pb,pg))
C calculate betap
      s34=tpa(4)**2-tpa(1)**2-tpa(2)**2-tpa(3)**2
      betasq=1._dp-4._dp*mt**2/s34
      if (betasq >= 0._dp) then
        be=sqrt(betasq)
        bp=0.5_dp*(1._dp+be)
        bm=1._dp-bp
        rtbp=sqrt(bp)
      else
        write(6,*) 'betasq < 0 in tdecayrodg.f, betasq=',betasq
        call flush(6)
        stop
      endif
      do si=1,4
      q(g,si)=p(pg,si)
      q(qq,si)=p(pnu,si)
      q(qb,si)=p(peb,si)
      q(b,si)=p(pb,si)-alb*p(pg,si)
      q(k5,si)=(bp*t(si)-bm*a(si))/be
      q(k6,si)=(bp*a(si)-bm*t(si))/be
      enddo
      call spinoru(6,q,za,zb)
      sw=s(qq,qb)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is the m(bpol,gpol,tpol)
      m(1,1,1)= + cprop**(-1)*twoptDg**(-1) * ( za(qq,b)*za(k5,g)*zb(qb
     &    ,k5)*zb(k5,b)*izb(b,g)*bp*bp*rtbp**(-1) + za(qq,b)*za(k6,g)*
     &    zb(qb,k6)*zb(k5,b)*izb(b,g)*bp*bm*rtbp**(-1) )
      m(1,1,1) = m(1,1,1) + cprop**(-1) * ( za(qq,g)*za(b,g)*zb(qb,k5)*
     &    iza(b,g)*izb(g,b)*bp*rtbp**(-1) + za(qq,g)*za(b,g)*zb(qb,k5)*
     &    iza(b,g)*izb(g,b)*bp*alb*rtbp**(-1) )
      m(1,1,1) = m(1,1,1) + cprop**(-1)*mb**2 * ( za(qq,g)*zb(qb,k5)*
     &    iza(b,g)*izb(b,g)*izb(g,b)*bp*rtbp**(-1) )
      m(1,1,1) = m(1,1,1) + cprop**(-1)*mt**2*twoptDg**(-1) * ( za(qq,b
     &    )*za(k6,g)*zb(qb,b)*iza(k5,k6)*izb(b,g)*rtbp**(-1) )

      m(1,2,1)= + cprop**(-1)*twoptDg**(-1) * (  - za(qq,b)*za(k5,b)*
     &    zb(qb,k5)*zb(k5,g)*iza(b,g)*bp*bp*rtbp**(-1) - za(qq,b)*za(k6
     &    ,b)*zb(qb,k6)*zb(k5,g)*iza(b,g)*bp*bm*rtbp**(-1) - za(qq,b)*
     &    zb(qb,g)*zb(k5,g)*bp*rtbp**(-1) )
      m(1,2,1) = m(1,2,1) + cprop**(-1)*mt**2*twoptDg**(-1) * (  - za(
     &    qq,b)*za(k6,b)*zb(qb,g)*iza(k5,k6)*iza(b,g)*rtbp**(-1) )

      m(2,1,1)= + cprop**(-1)*twoptDg**(-1)*mb * (  - za(qq,g)*za(k5,g)
     &    *zb(qb,k5)*zb(k5,b)*iza(b,g)*izb(b,g)*bp*bp*rtbp**(-1) - za(
     &    qq,g)*za(k6,g)*zb(qb,k6)*zb(k5,b)*iza(b,g)*izb(b,g)*bp*bm*
     &    rtbp**(-1) )
      m(2,1,1) = m(2,1,1) + cprop**(-1)*mt**2*twoptDg**(-1)*mb * (  - 
     &    za(qq,g)*za(k6,g)*zb(qb,b)*iza(k5,k6)*iza(b,g)*izb(b,g)*
     &    rtbp**(-1) )

      m(2,2,1)= + cprop**(-1)*twoptDg**(-1)*mb * ( za(qq,g)*za(k5,b)*
     &    zb(qb,k5)*zb(k5,g)*iza(b,g)**2*bp*bp*rtbp**(-1) + za(qq,g)*
     &    za(k6,b)*zb(qb,k6)*zb(k5,g)*iza(b,g)**2*bp*bm*rtbp**(-1) + 
     &    za(qq,g)*zb(qb,g)*zb(k5,g)*iza(b,g)*bp*rtbp**(-1) )
      m(2,2,1) = m(2,2,1) + cprop**(-1)*mt**2*twoptDg**(-1)*mb * ( za(
     &    qq,g)*za(k6,b)*zb(qb,g)*iza(k5,k6)*iza(b,g)**2*rtbp**(-1) )

      m(1,1,2)= + cprop**(-1)*mt*twoptDg**(-1) * ( za(qq,b)*za(k5,g)*
     &    zb(qb,k5)*zb(k6,b)*izb(k5,k6)*izb(b,g)*bp*rtbp**(-1) + za(qq,
     &    b)*za(k5,g)*zb(qb,b)*izb(b,g)*bp*rtbp**(-1) + za(qq,b)*za(k6,
     &    g)*zb(qb,k6)*zb(k6,b)*izb(k5,k6)*izb(b,g)*bm*rtbp**(-1) )
      m(1,1,2) = m(1,1,2) + cprop**(-1)*mt * ( za(qq,g)*za(b,g)*zb(qb,
     &    k6)*iza(b,g)*izb(k5,k6)*izb(g,b)*rtbp**(-1) + za(qq,g)*za(b,g
     &    )*zb(qb,k6)*iza(b,g)*izb(k5,k6)*izb(g,b)*alb*rtbp**(-1) )
      m(1,1,2) = m(1,1,2) + cprop**(-1)*mt*mb**2 * ( za(qq,g)*zb(qb,k6)
     &    *iza(b,g)*izb(k5,k6)*izb(b,g)*izb(g,b)*rtbp**(-1) )

      m(1,2,2)= + cprop**(-1)*mt*twoptDg**(-1) * (  - za(qq,b)*za(k5,b)
     &    *zb(qb,k5)*zb(k6,g)*iza(b,g)*izb(k5,k6)*bp*rtbp**(-1) - za(qq
     &    ,b)*za(k5,b)*zb(qb,g)*iza(b,g)*bp*rtbp**(-1) - za(qq,b)*za(k6
     &    ,b)*zb(qb,k6)*zb(k6,g)*iza(b,g)*izb(k5,k6)*bm*rtbp**(-1) - 
     &    za(qq,b)*zb(qb,g)*zb(k6,g)*izb(k5,k6)*rtbp**(-1) )

      m(2,1,2)= + cprop**(-1)*mt*twoptDg**(-1)*mb * (  - za(qq,g)*za(k5
     &    ,g)*zb(qb,k5)*zb(k6,b)*iza(b,g)*izb(k5,k6)*izb(b,g)*bp*
     &    rtbp**(-1) - za(qq,g)*za(k5,g)*zb(qb,b)*iza(b,g)*izb(b,g)*bp*
     &    rtbp**(-1) - za(qq,g)*za(k6,g)*zb(qb,k6)*zb(k6,b)*iza(b,g)*
     &    izb(k5,k6)*izb(b,g)*bm*rtbp**(-1) )

      m(2,2,2)= + cprop**(-1)*mt*twoptDg**(-1)*mb * ( za(qq,g)*za(k5,b)
     &    *zb(qb,k5)*zb(k6,g)*iza(b,g)**2*izb(k5,k6)*bp*rtbp**(-1) + 
     &    za(qq,g)*za(k5,b)*zb(qb,g)*iza(b,g)**2*bp*rtbp**(-1) + za(qq,
     &    g)*za(k6,b)*zb(qb,k6)*zb(k6,g)*iza(b,g)**2*izb(k5,k6)*bm*
     &    rtbp**(-1) + za(qq,g)*zb(qb,g)*zb(k6,g)*iza(b,g)*izb(k5,k6)*
     &    rtbp**(-1) )

      return
      end
