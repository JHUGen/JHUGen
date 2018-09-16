c--- File written by FORM program adecayrodg.frm on Thu May 24 10:26:48 CDT 2012
      subroutine adecayrodg(p,pnu,peb,pb,pc,pem,pnb,pg,m)
      implicit none
      include 'types.f'
      
C     antitopdecay with radiation from t-b line
C     b is rendered massless wrt to pqq
C     t is rendered massless wrt to pqb
C     pqq,pqb,pc point to the anti-top decay products
C     pg points to the radiated gluon
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),dot,sw,twopaDg,
     & alc,a(4),t(4),tpa(4),s34,betasq,be,bp,bm,rtbp
      complex(dp):: m(2,2,2),cprop,iza,izb
      integer:: si,qq,qb,c,g,pnu,peb,pb,pc,pem,pnb,pg,aa,bb,k5,k6
      parameter(g=1,qq=2,qb=3,c=4,k5=5,k6=6)
C---statement functions
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
C---end statement functions
C construct top and antitop momenta
      do si=1,4
      t(si)=p(pnu,si)+p(peb,si)+p(pb,si)
      a(si)=p(pem,si)+p(pnb,si)+p(pc,si)+p(pg,si)
      tpa(si)=t(si)+a(si)
      enddo
      twopaDg=2._dp*(a(4)*p(pg,4)-a(1)*p(pg,1)-a(2)*p(pg,2)-a(3)*p(pg,3))
      alc=mb**2/(2._dp*dot(p,pc,pg))
C calculate betap
      s34=tpa(4)**2-tpa(1)**2-tpa(2)**2-tpa(3)**2
      betasq=1._dp-4._dp*mt**2/s34
      if (betasq >= 0._dp) then
        be=sqrt(betasq)
        bp=0.5_dp*(1._dp+be)
        bm=1._dp-bp
        rtbp=sqrt(bp)
      else
        write(6,*) 'betasq < 0 in adecayrodg.f, betasq=',betasq
        call flush(6)
        stop
      endif
      do si=1,4
      q(g,si)=p(pg,si)
      q(qq,si)=p(pem,si)
      q(qb,si)=p(pnb,si)
      q(c,si)=p(pc,si)-alc*p(pg,si)
      q(k5,si)=(bp*t(si)-bm*a(si))/be
      q(k6,si)=(bp*a(si)-bm*t(si))/be
      enddo
      call spinoru(6,q,za,zb)
      sw=s(qq,qb)
      cprop=cplx2(sw-wmass**2,wmass*wwidth)
C---order of polarizations is the m(apol,gpol,cpol)
      m(1,1,1)= + cprop**(-1)*mt*mb * (  - za(qq,k5)*za(k5,g)*zb(qb,g)*
     &    zb(k5,c)*iza(k5,k6)*izb(c,g)**2*bm*rtbp**(-1)*twopaDg**(-1)
     &     - za(qq,k6)*za(k5,g)*zb(qb,g)*zb(k6,c)*iza(k5,k6)*izb(c,g)**
     &    2*bp*rtbp**(-1)*twopaDg**(-1) - za(qq,g)*za(k5,g)*zb(qb,g)*
     &    iza(k5,k6)*izb(c,g)*rtbp**(-1)*twopaDg**(-1) + za(qq,g)*zb(qb
     &    ,g)*zb(k6,c)*izb(c,g)**2*bp*rtbp**(-1)*twopaDg**(-1) )

      m(1,2,1)= + cprop**(-1)*mt*mb * ( za(qq,k5)*za(k5,c)*zb(qb,g)*zb(
     &    k5,g)*iza(k5,k6)*iza(c,g)*izb(c,g)*bm*rtbp**(-1)*
     &    twopaDg**(-1) + za(qq,k6)*za(k5,c)*zb(qb,g)*zb(k6,g)*iza(k5,
     &    k6)*iza(c,g)*izb(c,g)*bp*rtbp**(-1)*twopaDg**(-1) - za(qq,c)*
     &    zb(qb,g)*zb(k6,g)*iza(c,g)*izb(c,g)*bp*rtbp**(-1)*
     &    twopaDg**(-1) )

      m(2,1,1)= + cprop**(-1)*mb * ( za(qq,k5)*za(k6,g)*zb(qb,g)*zb(k5,
     &    c)*izb(c,g)**2*bp*bm*rtbp**(-1)*twopaDg**(-1) + za(qq,k6)*za(
     &    k6,g)*zb(qb,g)*zb(k6,c)*izb(c,g)**2*bp*bp*rtbp**(-1)*
     &    twopaDg**(-1) + za(qq,g)*za(k6,g)*zb(qb,g)*izb(c,g)*bp*
     &    rtbp**(-1)*twopaDg**(-1) )
      m(2,1,1) = m(2,1,1) + cprop**(-1)*mt**2*mb * (  - za(qq,g)*zb(qb,
     &    g)*zb(k5,c)*izb(k5,k6)*izb(c,g)**2*rtbp**(-1)*twopaDg**(-1) )

      m(2,2,1)= + cprop**(-1)*mb * (  - za(qq,k5)*za(k6,c)*zb(qb,g)*zb(
     &    k5,g)*iza(c,g)*izb(c,g)*bp*bm*rtbp**(-1)*twopaDg**(-1) - za(
     &    qq,k6)*za(k6,c)*zb(qb,g)*zb(k6,g)*iza(c,g)*izb(c,g)*bp*bp*
     &    rtbp**(-1)*twopaDg**(-1) )
      m(2,2,1) = m(2,2,1) + cprop**(-1)*mt**2*mb * ( za(qq,c)*zb(qb,g)*
     &    zb(k5,g)*iza(c,g)*izb(k5,k6)*izb(c,g)*rtbp**(-1)*
     &    twopaDg**(-1) )

      m(1,1,2)= + cprop**(-1)*mt * ( za(qq,k5)*za(k5,g)*zb(qb,c)*zb(k5,
     &    c)*iza(k5,k6)*izb(c,g)*bm*rtbp**(-1)*twopaDg**(-1) + za(qq,k6
     &    )*za(k5,g)*zb(qb,c)*zb(k6,c)*iza(k5,k6)*izb(c,g)*bp*
     &    rtbp**(-1)*twopaDg**(-1) + za(qq,g)*za(k5,g)*zb(qb,c)*iza(k5,
     &    k6)*rtbp**(-1)*twopaDg**(-1) - za(qq,g)*zb(qb,c)*zb(k6,c)*
     &    izb(c,g)*bp*rtbp**(-1)*twopaDg**(-1) )

      m(1,2,2)= + cprop**(-1)*mt * (  - za(qq,k5)*za(k5,c)*zb(qb,c)*zb(
     &    k5,g)*iza(k5,k6)*iza(c,g)*bm*rtbp**(-1)*twopaDg**(-1) - za(qq
     &    ,k5)*zb(qb,g)*zb(c,g)*iza(k5,k6)*iza(c,g)*izb(g,c)*rtbp**(-1)
     &     - za(qq,k5)*zb(qb,g)*zb(c,g)*iza(k5,k6)*iza(c,g)*izb(g,c)*
     &    alc*rtbp**(-1) - za(qq,k6)*za(k5,c)*zb(qb,c)*zb(k6,g)*iza(k5,
     &    k6)*iza(c,g)*bp*rtbp**(-1)*twopaDg**(-1) + za(qq,c)*zb(qb,c)*
     &    zb(k6,g)*iza(c,g)*bp*rtbp**(-1)*twopaDg**(-1) )
      m(1,2,2) = m(1,2,2) + cprop**(-1)*mt*mb**2 * (  - za(qq,k5)*zb(qb
     &    ,g)*iza(k5,k6)*iza(c,g)**2*izb(g,c)*rtbp**(-1) )

      m(2,1,2)= + cprop**(-1) * (  - za(qq,k5)*za(k6,g)*zb(qb,c)*zb(k5,
     &    c)*izb(c,g)*bp*bm*rtbp**(-1)*twopaDg**(-1) - za(qq,k6)*za(k6,
     &    g)*zb(qb,c)*zb(k6,c)*izb(c,g)*bp*bp*rtbp**(-1)*twopaDg**(-1)
     &     - za(qq,g)*za(k6,g)*zb(qb,c)*bp*rtbp**(-1)*twopaDg**(-1) )
      m(2,1,2) = m(2,1,2) + cprop**(-1)*mt**2 * ( za(qq,g)*zb(qb,c)*zb(
     &    k5,c)*izb(k5,k6)*izb(c,g)*rtbp**(-1)*twopaDg**(-1) )

      m(2,2,2)= + cprop**(-1) * ( za(qq,k5)*za(k6,c)*zb(qb,c)*zb(k5,g)*
     &    iza(c,g)*bp*bm*rtbp**(-1)*twopaDg**(-1) + za(qq,k6)*za(k6,c)*
     &    zb(qb,c)*zb(k6,g)*iza(c,g)*bp*bp*rtbp**(-1)*twopaDg**(-1) + 
     &    za(qq,k6)*zb(qb,g)*zb(c,g)*iza(c,g)*izb(g,c)*bp*rtbp**(-1) + 
     &    za(qq,k6)*zb(qb,g)*zb(c,g)*iza(c,g)*izb(g,c)*bp*alc*
     &    rtbp**(-1) )
      m(2,2,2) = m(2,2,2) + cprop**(-1)*mb**2 * ( za(qq,k6)*zb(qb,g)*
     &    iza(c,g)**2*izb(g,c)*bp*rtbp**(-1) )
      m(2,2,2) = m(2,2,2) + cprop**(-1)*mt**2 * (  - za(qq,c)*zb(qb,c)*
     &    zb(k5,g)*iza(c,g)*izb(k5,k6)*rtbp**(-1)*twopaDg**(-1) )

      return
      end
