c--- File written by FORM program tdecayrod_v.frm on Thu May 24 10:26:43 CDT 2012
      subroutine tdecayrod_v(p,nu,eb,b,bb,e,nb,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, May 2012                                     *
*     Virtual corrections to                                           *
*     top decay  t --> nu(nu)+eb(eb)+b(pb)                             *
*     anti-top decay  t --> em(e)+nb(nb)+bb(bb)                        *
*     with bottom and top masses (and no radiation)                    *
*     and with overall factor of (4 pi)^e/Gamma(1-e)                   *
*     in massless spinor notation                                      *
*     nu,eb,b,bb,e,nb are integer::s that point to                       *
*     the appropriate four-momenta in p                                *
*     xnu=neutrino or quark                                            *
*     xeb=positron or antiquark                                        *
*     pb=bottom quark                                                  *
*     q(t) is rendered massless in the rodrigo scheme                  *
*     q(b) is rendered massless wrt to nu                              *
*     returned m(bpol,tpol)                                            *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'qcdcouple.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),sw,alb,bp,bm,ctm,
     & t(4),a(4),tpa(4),s34,be,betasq,rtbp,dot,corr,nloratiotopdecay
      complex(dp):: m(2,2),cprop,c0L,c0R,c1L,C1R,c1Lon2,c1Ron2,
     & iza,izb
      integer:: nu,eb,b,e,nb,bb,aa
      integer:: si,xnu,xeb,xb,xk3,xk4
      parameter(xnu=1,xeb=2,xb=3,xk3=4,xk4=5)
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
      corr=nloratiotopdecay(mt,mb,wmass,wwidth)/(cf*ason2pi)
C construct top and antitop momenta
      do si=1,4
      t(si)=p(nu,si)+p(eb,si)+p(b,si)
      a(si)=p(e,si)+p(nb,si)+p(bb,si)
      tpa(si)=t(si)+a(si)
      enddo
C calculate betap
      s34=tpa(4)**2-tpa(1)**2-tpa(2)**2-tpa(3)**2
      betasq=1._dp-4._dp*mt**2/s34
      if (betasq >= 0._dp) then
        be=sqrt(betasq)
        bp=0.5_dp*(1._dp+be)
        bm=1._dp-bp
        rtbp=sqrt(bp)
      else
        write(6,*) 'betasq < 0 in tdecayrod_v.f, betasq=',betasq
        call flush(6)
        stop
      endif
      alb=mb**2/(2._dp*dot(p,b,nu))
      do si=1,4
      q(xnu,si)=p(nu,si)
      q(xeb,si)=p(eb,si)
      q(xb,si)=p(b,si)-alb*p(nu,si)
      q(xk3,si)=(bp*t(si)-bm*a(si))/be
      q(xk4,si)=(bp*a(si)-bm*t(si))/be
      enddo
      call spinoru(5,q,za,zb)
      sw=2._dp*dot(q,xeb,xnu)
      call coefsdkmass(sw,mt,mb,ctm,c0L,c0R,c1L,c1R)
      c0L=c0L+cplx1(ctm-corr)
      c1Lon2=c1L/2._dp
      c1Ron2=c1R/2._dp
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is m(bpol,tpol)
      m(1,1)= + cprop**(-1)*rtbp**(-1)*mt**(-1)*mb * ( za(xb,xnu)*zb(xb
     &    ,xeb)*zb(xnu,xk3)*izb(xb,xnu)*bp*c1Lon2 )
      m(1,1) = m(1,1) + cprop**(-1)*rtbp**(-1) * (  - za(xb,xnu)*za(xb,
     &    xk4)*zb(xb,xeb)*iza(xk3,xk4)*c1Ron2 + za(xb,xnu)*zb(xeb,xk3)*
     &    bp*c0L )
      m(1,1) = m(1,1) + cprop**(-1)*rtbp**(-1)*mt*mb * (  - za(xnu,xk4)
     &    *zb(xnu,xeb)*iza(xk3,xk4)*izb(xb,xnu)*c0R )

      m(1,2)= + cprop**(-1)*rtbp**(-1)*mt**(-1) * (  - za(xb,xnu)*za(xb
     &    ,xk3)*zb(xb,xeb)*bp*c1Ron2 )
      m(1,2) = m(1,2) + cprop**(-1)*rtbp**(-1)*mb * ( za(xb,xnu)*zb(xb,
     &    xeb)*zb(xnu,xk4)*izb(xb,xnu)*izb(xk3,xk4)*c1Lon2 - za(xnu,xk3
     &    )*zb(xnu,xeb)*izb(xb,xnu)*bp*c0R )
      m(1,2) = m(1,2) + cprop**(-1)*rtbp**(-1)*mt * ( za(xb,xnu)*zb(xeb
     &    ,xk4)*izb(xk3,xk4)*c0L )

      m(2,1)= + cprop**(-1)*rtbp**(-1)*mt**(-1) * (  - za(xb,xnu)*zb(xb
     &    ,xeb)*zb(xb,xk3)*bp*c1Lon2 )
      m(2,1) = m(2,1) + cprop**(-1)*rtbp**(-1)*mb * ( za(xnu,xk4)*zb(xb
     &    ,xeb)*iza(xk3,xk4)*c1Ron2 )
      m(2,1) = m(2,1) + cprop**(-1)*rtbp**(-1)*mt * ( za(xnu,xk4)*zb(xb
     &    ,xeb)*iza(xk3,xk4)*c0R )

      m(2,2)= + cprop**(-1)*rtbp**(-1)*mt**(-1)*mb * ( za(xnu,xk3)*zb(
     &    xb,xeb)*bp*c1Ron2 )
      m(2,2) = m(2,2) + cprop**(-1)*rtbp**(-1) * (  - za(xb,xnu)*zb(xb,
     &    xeb)*zb(xb,xk4)*izb(xk3,xk4)*c1Lon2 + za(xnu,xk3)*zb(xb,xeb)*
     &    bp*c0R )

      return
      end
