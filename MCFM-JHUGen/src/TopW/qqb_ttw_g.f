      subroutine qqb_ttw_g(pin,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     March, 2012.                                                     *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  g*  + W +g(p11)
c                           |     |
c                           |      --> nu(p9)+e^+(p10)
c                           |
c                           ---> t(p3+p4+p5)+t(p6+p7+p8)
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'plabel.f'
      integer j,k,nu
      double precision pin(mxpart,4),p(mxpart,4),q(mxpart,4),dot,
     & msq(-nf:nf,-nf:nf),trodmsqm,fac,bp,bm,beta,betasq,s56,tmass
      double precision qqbWbbg,qbqWbbg,qgWbbq,gqWbbq,gqbWbbqb,qbgWbbqb
      double complex mtop(2,2),manti(2,2)
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---setup the decay matrix elements common for all crossings
      call tdecayrod(pin,3,4,5,6,7,8,0,mtop)
      call adecayrod(pin,3,4,5,6,7,8,0,manti)

C--define concatenated momenta
      do nu=1,4
      p(1,nu)=pin(1,nu)
      p(2,nu)=pin(2,nu)
      p(3,nu)=pin(9,nu)
      p(4,nu)=pin(10,nu)
      p(5,nu)=pin(3,nu)+pin(4,nu)+pin(5,nu)
      p(6,nu)=pin(6,nu)+pin(7,nu)+pin(8,nu)
      p(7,nu)=pin(11,nu)
      enddo

      tmass=mt
      fac=gsq**3*gwsq**6/4d0*32d0/(mt*twidth)**4

c--- include factor for hadronic decays of W
      if (plabel(3) .eq. 'pp') fac=2d0*xn*fac
      if (plabel(7) .eq. 'pp') fac=2d0*xn*fac

C calculate betap
      s56=2d0*mt**2+2d0*dot(p,5,6)
      betasq=1d0-4d0*mt**2/s56
      if (betasq .ge. 0d0) then
        beta=dsqrt(betasq)
        bp=0.5d0*(1d0+beta)
        bm=1d0-bp
      else
        write(6,*) 'betasq < 0 in qqb_ttw_g.f, betasq=',betasq
        call flush(6)
        stop
      endif


C---Create rodrigo momenta
      do j=1,7
      do nu=1,4
      if (j.eq.5) then
      q(j,nu)=(bp/beta)*p(5,nu)-(bm/beta)*p(6,nu)
      elseif (j.eq.6) then
      q(j,nu)=(bp/beta)*p(6,nu)-(bm/beta)*p(5,nu)
      else
      q(j,nu)=p(j,nu)
      endif
      enddo
      enddo
C--- and calculate their spinor products     
      call spinoru(7,q,za,zb)

c--- q-qb and qb-q
      qqbWbbg =+trodmsqm(2,1,5,6,7,4,3,pin,tmass,mtop,manti)*fac*aveqq
      qbqWbbg =+trodmsqm(1,2,5,6,7,4,3,pin,tmass,mtop,manti)*fac*aveqq

c--- q-g and g-q
      qgWbbq  =+trodmsqm(7,1,5,6,2,4,3,pin,tmass,mtop,manti)*fac*aveqg
      gqWbbq  =+trodmsqm(7,2,5,6,1,4,3,pin,tmass,mtop,manti)*fac*aveqg

c--- g-qb and qb-g
      gqbWbbqb=+trodmsqm(2,7,5,6,1,4,3,pin,tmass,mtop,manti)*fac*aveqg
      qbgWbbqb=+trodmsqm(1,7,5,6,2,4,3,pin,tmass,mtop,manti)*fac*aveqg

      do j=-nf,nf
      do k=-nf,nf

      if     ((j .gt. 0) .and. (k .lt. 0)) then
      msq(j,k)=Vsq(j,k)*qqbWbbg

      elseif ((j .lt. 0) .and. (k .gt. 0)) then
      msq(j,k)=Vsq(j,k)*qbqWbbg

      elseif ((j .gt. 0) .and. (k .eq. 0)) then
      msq(j,k)=Vsum(j)*qgWbbq 
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
      msq(j,k)=Vsum(j)*qbgWbbqb

      elseif ((j .eq. 0) .and. (k .gt. 0)) then
      msq(j,k)=Vsum(k)*gqWbbq

      elseif ((j .eq. 0) .and. (k .lt. 0)) then
      msq(j,k)=Vsum(k)*gqbWbbqb
      endif

      enddo
      enddo

      return
      end



