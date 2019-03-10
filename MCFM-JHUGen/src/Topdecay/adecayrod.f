      subroutine adecayrod(p,nu,eb,b,c,em,nb,g,m)
      implicit none
C     anti-topdecay
C     c is rendered massless wrt to nb
C     a is rendered massless via the rodrigo scheme
      include 'constants.f'
      include 'zprods_decl.f'
      include 'masses.f'
      double precision p(mxpart,4),q(mxpart,4),sw,alb,bp,bm,
     & t(4),a(4),tpa(4),s34,be,betasq,rtbp,dot
      double complex m(2,2),cprop
      integer nu,eb,b,em,nb,c,g
      integer si,xnb,xem,xc,xk3,xk4
      parameter(xnb=1,xem=2,xc=3,xk3=4,xk4=5)
C construct top and antitop momenta
      do si=1,4
      if (g .eq. 0) then
      t(si)=p(nu,si)+p(eb,si)+p(b,si)
      else
      t(si)=p(nu,si)+p(eb,si)+p(b,si)+p(g,si)
      endif
      a(si)=p(em,si)+p(nb,si)+p(c,si)
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
        write(6,*) 'betasq < 0 in adecayrod.f, betasq=',betasq
        call flush(6)
        stop
      endif
      alb=mb**2/(2d0*dot(p,c,nb))

      do si=1,4
      q(xnb,si)=p(nb,si)
      q(xem,si)=p(em,si)
      q(xc,si)=p(c,si)-alb*p(nb,si)
C     k3=+bp/beta*t-bm/beta*a
C     k4=-bm/beta*t+bp/beta*a
C     t=bp*k3+bm*k4
C     a=bm*k3+bp*k4
      q(xk3,si)=(bp*t(si)-bm*a(si))/be
      q(xk4,si)=(bp*a(si)-bm*t(si))/be
      enddo

      call spinoru(5,q,za,zb)
      
      sw=2d0*dot(q,xnb,xem)
      cprop=dcmplx(sw-wmass**2,wmass*wwidth)

C---order of polarizations is m(apol,cpol)
C---choice of auxiliary vector for a is xk3
C---choice of auxiliary vector for c is xnb
      m(1,1)= czip
      m(2,1)= czip
      m(1,2)= -mt*za(xk3,xem)/za(xk3,xk4)*zb(xnb,xc)/(rtbp*cprop)
      m(2,2)= +za(xk4,xem)*zb(xnb,xc)*rtbp/cprop

c--debug for comparison with stable top
c      m(1,1)= mt*twidth/gwsq
c      m(1,2)= czip
c      m(2,1)= czip
c      m(2,2)= mt*twidth/gwsq
c--debug
      return
      end


