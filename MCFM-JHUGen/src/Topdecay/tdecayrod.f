      subroutine tdecayrod(p,nu,eb,b,bb,e,nb,g,m)
      implicit none
C     topdecay
C     b is rendered massless wrt to nu
C     t is rendered massless wrt the rodrigo scheme
      include 'constants.f'
      include 'zprods_decl.f'
      include 'masses.f'
      double precision p(mxpart,4),q(mxpart,4),sw,alb,bp,bm,
     & t(4),a(4),tpa(4),s34,be,betasq,rtbp,dot
      double complex m(2,2),cprop
      integer nu,eb,b,e,nb,bb,g
      integer si,xnu,xeb,xb,xk3,xk4
      parameter(xnu=1,xeb=2,xb=3,xk3=4,xk4=5)

C construct top and antitop momenta
      do si=1,4
      t(si)=p(nu,si)+p(eb,si)+p(b,si)
      if (g .eq. 0) then
      a(si)=p(e,si)+p(nb,si)+p(bb,si)
      else
      a(si)=p(e,si)+p(nb,si)+p(bb,si)+p(g,si)
      endif
      tpa(si)=t(si)+a(si)
      enddo

C calculate betap
      s34=tpa(4)**2-tpa(1)**2-tpa(2)**2-tpa(3)**2
      betasq=1d0-4d0*mt**2/s34
      if (betasq .ge. 0d0) then
        be=dsqrt(betasq)
        bp=0.5d0*(1d0+be)
        bm=1d0-bp
        rtbp=sqrt(bp)
      else
        write(6,*) 'betasq < 0 in tdecayrod.f, betasq=',betasq
        call flush(6)
        stop
      endif

      alb=mb**2/(2d0*dot(p,b,nu))

      do si=1,4
      q(xnu,si)=p(nu,si)
      q(xeb,si)=p(eb,si)
      q(xb,si)=p(b,si)-alb*p(nu,si)
      q(xk3,si)=(bp*t(si)-bm*a(si))/be
      q(xk4,si)=(bp*a(si)-bm*t(si))/be
      enddo
      call spinoru(5,q,za,zb)
      sw=2d0*dot(q,xeb,xnu)
      cprop=dcmplx(sw-wmass**2,wmass*wwidth)

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
