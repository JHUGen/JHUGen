      subroutine tdecayrod(p,nu,eb,b,bb,e,nb,g,m)
      implicit none
      include 'types.f'
      
C     topdecay
C     b is rendered massless wrt to nu
C     t is rendered massless wrt the rodrigo scheme
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),q(mxpart,4),sw,alb,bp,bm,
     & t(4),a(4),tpa(4),s34,be,betasq,rtbp,dot
      complex(dp):: m(2,2),cprop
      integer:: nu,eb,b,e,nb,bb,g
      integer:: si,xnu,xeb,xb,xk3,xk4
      parameter(xnu=1,xeb=2,xb=3,xk3=4,xk4=5)

C construct top and antitop momenta
      do si=1,4
      t(si)=p(nu,si)+p(eb,si)+p(b,si)
      if (g == 0) then
      a(si)=p(e,si)+p(nb,si)+p(bb,si)
      else
      a(si)=p(e,si)+p(nb,si)+p(bb,si)+p(g,si)
      endif
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
        write(6,*) 'betasq < 0 in tdecayrod.f, betasq=',betasq
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
      cprop=cplx2(sw-wmass**2,wmass*wwidth)

C---order of polarizations is m(bpol,tpol)
      m(1,1)= za(xb,xnu)*zb(xeb,xk3)*rtbp/cprop
      m(1,2)= mt*za(xb,xnu)*zb(xeb,xk4)/zb(xk3,xk4)/(rtbp*cprop)
      m(2,1)= czip
      m(2,2)= czip


c--debug for comparison with stable top
c      m(1,1)= mt*twidth/gwsq
c      m(1,2)= czip
c      m(2,1)= czip
c      m(2,2)= mt*twidth/gwsq
c--debug
      return
      end
