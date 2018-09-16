c--- File written by FORM program adecayrod_v.frm on Thu May 24 10:26:46 CDT 2012
      subroutine adecayrod_v(p,nu,eb,b,c,em,nb,m)
      implicit none
      include 'types.f'
      
C     anti-topdecay
C     c is rendered massless wrt to nb
C     a is rendered massless via the rodrigo scheme
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'qcdcouple.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),sw,alc,bp,bm,ctm,
     & t(4),a(4),tpa(4),s34,be,betasq,rtbp,dot,corr,nloratiotopdecay
      complex(dp):: m(2,2),cprop,iza,izb
      complex(dp):: c0L,c0R,c1L,C1R,c1Lon2,c1Ron2
      integer:: nu,eb,b,em,nb,c
      integer:: si,xnb,xem,xc,xk3,xk4,aa,bb
      parameter(xnb=1,xem=2,xc=3,xk3=4,xk4=5)
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
      corr=nloratiotopdecay(mt,mb,wmass,wwidth)/(cf*ason2pi)
C construct top and antitop momenta
      do si=1,4
      t(si)=p(nu,si)+p(eb,si)+p(b,si)
      a(si)=p(em,si)+p(nb,si)+p(c,si)
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
        write(6,*) 'betasq < 0 in adecayrod_v.f, betasq=',betasq
        call flush(6)
        stop
      endif
      alc=mb**2/(2._dp*dot(p,c,nb))
      do si=1,4
      q(xnb,si)=p(nb,si)
      q(xem,si)=p(em,si)
      q(xc,si)=p(c,si)-alc*p(nb,si)
C     k3=+bp/beta*t-bm/beta*a
C     k4=-bm/beta*t+bp/beta*a
C     t=bp*k3+bm*k4
C     a=bm*k3+bp*k4
      q(xk3,si)=(bp*t(si)-bm*a(si))/be
      q(xk4,si)=(bp*a(si)-bm*t(si))/be
      enddo
      call spinoru(5,q,za,zb)
      sw=2._dp*dot(q,xem,xnb)
      call coefsdkmass(sw,mt,mb,ctm,c0L,c0R,c1L,c1R)
      c0L=c0L+cplx1(ctm-corr)
      c1Lon2=c1L/2._dp
      c1Ron2=c1R/2._dp
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is m(apol,cpol)
      m(1,1)= + cprop**(-1)*mt**(-1)*mb * (  - za(xc,xem)*zb(xk4,xnb)*
     &    bp*c1Ron2*rtbp**(-1) )
      m(1,1) = m(1,1) + cprop**(-1) * ( za(xk3,xc)*za(xc,xem)*zb(xnb,xc
     &    )*iza(xk3,xk4)*c1Lon2*rtbp**(-1) - za(xc,xem)*zb(xk4,xnb)*bp*
     &    c0R*rtbp**(-1) )

      m(1,2)= + cprop**(-1)*mt**(-1) * (  - za(xc,xem)*zb(xk4,xc)*zb(
     &    xnb,xc)*bp*c1Ron2*rtbp**(-1) )
      m(1,2) = m(1,2) + cprop**(-1)*mb * ( za(xk3,xnb)*za(xc,xem)*zb(
     &    xnb,xc)*iza(xk3,xk4)*iza(xnb,xc)*c1Lon2*rtbp**(-1) - za(xnb,
     &    xem)*zb(xk4,xnb)*iza(xnb,xc)*bp*c0R*rtbp**(-1) )
      m(1,2) = m(1,2) + cprop**(-1)*mt * (  - za(xk3,xem)*zb(xnb,xc)*
     &    iza(xk3,xk4)*c0L*rtbp**(-1) )

      m(2,1)= + cprop**(-1)*mt**(-1) * (  - za(xk4,xc)*za(xc,xem)*zb(
     &    xnb,xc)*bp*c1Lon2*rtbp**(-1) )
      m(2,1) = m(2,1) + cprop**(-1)*mb * ( za(xc,xem)*zb(xk3,xnb)*izb(
     &    xk3,xk4)*c1Ron2*rtbp**(-1) )
      m(2,1) = m(2,1) + cprop**(-1)*mt * ( za(xc,xem)*zb(xk3,xnb)*izb(
     &    xk3,xk4)*c0R*rtbp**(-1) )

      m(2,2)= + cprop**(-1)*mt**(-1)*mb * (  - za(xk4,xnb)*za(xc,xem)*
     &    zb(xnb,xc)*iza(xnb,xc)*bp*c1Lon2*rtbp**(-1) )
      m(2,2) = m(2,2) + cprop**(-1) * ( za(xk4,xem)*zb(xnb,xc)*bp*c0L*
     &    rtbp**(-1) + za(xc,xem)*zb(xk3,xc)*zb(xnb,xc)*izb(xk3,xk4)*
     &    c1Ron2*rtbp**(-1) )
      m(2,2) = m(2,2) + cprop**(-1)*mt*mb * ( za(xnb,xem)*zb(xk3,xnb)*
     &    iza(xnb,xc)*izb(xk3,xk4)*c0R*rtbp**(-1) )

C---choice of auxiliary vector for a is xk3
C---choice of auxiliary vector for c is xnb
      return
      end
