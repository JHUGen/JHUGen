      subroutine qqb_ttw_g(pin,msq)
      implicit none
      include 'types.f'
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

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'plabel.f'
      integer:: j,k,nu
      real(dp):: pin(mxpart,4),p(mxpart,4),q(mxpart,4),dot,
     & msq(-nf:nf,-nf:nf),trodmsqm,fac,bp,bm,beta,betasq,s56,tmass
      real(dp):: qqbWbbg,qbqWbbg,qgWbbq,gqWbbq,gqbWbbqb,qbgWbbqb
      complex(dp):: mtop(2,2),manti(2,2)
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
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
      fac=gsq**3*gwsq**6/4._dp*32._dp/(mt*twidth)**4

c--- include factor for hadronic decays of W
      if (plabel(3) == 'pp') fac=2._dp*xn*fac
      if (plabel(7) == 'pp') fac=2._dp*xn*fac

C calculate betap
      s56=2._dp*mt**2+2._dp*dot(p,5,6)
      betasq=1._dp-4._dp*mt**2/s56
      if (betasq >= 0._dp) then
        beta=sqrt(betasq)
        bp=0.5_dp*(1._dp+beta)
        bm=1._dp-bp
      else
        write(6,*) 'betasq < 0 in qqb_ttw_g.f, betasq=',betasq
        call flush(6)
        stop
      endif


C---Create rodrigo momenta
      do j=1,7
      do nu=1,4
      if (j==5) then
      q(j,nu)=(bp/beta)*p(5,nu)-(bm/beta)*p(6,nu)
      elseif (j==6) then
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

      if     ((j > 0) .and. (k < 0)) then
      msq(j,k)=Vsq(j,k)*qqbWbbg

      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=Vsq(j,k)*qbqWbbg

      elseif ((j > 0) .and. (k == 0)) then
      msq(j,k)=Vsum(j)*qgWbbq
      elseif ((j < 0) .and. (k == 0)) then
      msq(j,k)=Vsum(j)*qbgWbbqb

      elseif ((j == 0) .and. (k > 0)) then
      msq(j,k)=Vsum(k)*gqWbbq

      elseif ((j == 0) .and. (k < 0)) then
      msq(j,k)=Vsum(k)*gqbWbbqb
      endif

      enddo
      enddo

      return
      end



