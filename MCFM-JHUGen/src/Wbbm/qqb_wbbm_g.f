      subroutine qqb_wbbm_g(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     October, 2010.                                                      *
************************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  g*  + W +g(p7)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+bb(p6)
c   with mass for the b and the bbar

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'heavyflav.f'
      integer:: j,k,nu
      real(dp):: p(mxpart,4),q(mxpart,4),
     & msq(-nf:nf,-nf:nf),rodmsqm,fac,bp,bm,beta,s56,mQsq,mq
      real(dp):: qqbWbbg,qbqWbbg,qgWbbq,gqWbbq,gqbWbbqb,qbgWbbqb

C--- set up the correct mass, according to 'flav'
      if     (flav == 6) then
        mQsq=mt**2
      elseif (flav == 5) then
        mQsq=mb**2
      elseif (flav == 4) then
        mQsq=mc**2
      else
        write(6,*) 'Wrong flavour in qqb_wbbm_v.f: flav=',flav
        call flush(6)
        stop
      endif
      mq=sqrt(mQsq)

C----Initialize whole array to zero
      msq(:,:)=zero

      fac=gsq**3*gw**4/four*32._dp

      s56=two*mQsq
     & +two*(+p(5,4)*p(6,4)-p(5,1)*p(6,1)-p(5,2)*p(6,2)-p(5,3)*p(6,3))
      beta=sqrt(one-four*mQsq/s56)
      bp=half+half*beta
      bm=half-half*beta

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

      call spinoru(7,q,za,zb)

c--- q-qb and qb-q
      qqbWbbg =+rodmsqm(2,1,5,6,7,4,3,mq)*fac*aveqq
      qbqWbbg =+rodmsqm(1,2,5,6,7,4,3,mq)*fac*aveqq

c--- q-g and g-q
      qgWbbq  =+rodmsqm(7,1,5,6,2,4,3,mq)*fac*aveqg
      gqWbbq  =+rodmsqm(7,2,5,6,1,4,3,mq)*fac*aveqg

c--- g-qb and qb-g
      gqbWbbqb=+rodmsqm(2,7,5,6,1,4,3,mq)*fac*aveqg
      qbgWbbqb=+rodmsqm(1,7,5,6,2,4,3,mq)*fac*aveqg

      do j=-(flav-1),(flav-1)
      do k=-(flav-1),(flav-1)

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


      function rodmsqm(j1,j2,j3,j4,j5,j6,j7,mQ)
      implicit none
      include 'types.f'
      real(dp):: rodmsqm

c matrix element squared summed over colors and spins
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6,j7,h3,h4,h5
      complex(dp):: qcda(2,2,2),qcdb(2,2,2),qedi(2,2,2),qedf(2,2,2)
      real(dp):: prop,mQ
c---calculate the W propagator
      prop=((s(j6,j7)-wmass**2)**2+(wmass*wwidth)**2)

C---These two calls exploit the symmetry under 1<->2,3<->4,6<->7,za<->zb
C---and overall sign change
      call Wbb(j1,j2,j3,j4,j5,j6,j7,mQ,za,zb,1,qedi,qedf,qcda,qcdb)
      call Wbb(j2,j1,j4,j3,j5,j7,j6,mQ,zb,za,2,qedi,qedf,qcda,qcdb)

      rodmsqm=zip
      do h3=1,2
      do h4=1,2
      do h5=1,2

      rodmsqm=rodmsqm+
     & V*xn/eight*(abs(qcda(h3,h4,h5))**2+abs(qcdb(h3,h4,h5))**2)
     &+V/(eight*xn)*(abs(qedi(h3,h4,h5))**2+abs(qedf(h3,h4,h5))**2
     &-two*(abs(qedi(h3,h4,h5)+qedf(h3,h4,h5)))**2)
      enddo
      enddo
      enddo

      rodmsqm=rodmsqm/prop
      return
      end
