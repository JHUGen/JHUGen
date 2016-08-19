      subroutine qqb_wbbm_g(p,msq)
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
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'heavyflav.f'
      integer j,k,nu
      double precision p(mxpart,4),q(mxpart,4),
     . msq(-nf:nf,-nf:nf),rodmsqm,fac,bp,bm,beta,s56,mQsq,mq
      double precision qqbWbbg,qbqWbbg,qgWbbq,gqWbbq,gqbWbbqb,qbgWbbqb

C--- set up the correct mass, according to 'flav'
      if     (flav .eq. 6) then
        mQsq=mt**2
      elseif (flav .eq. 5) then
        mQsq=mb**2
      elseif (flav .eq. 4) then
        mQsq=mc**2
      else
        write(6,*) 'Wrong flavour in qqb_wbbm_v.f: flav=',flav
        call flush(6)
        stop
      endif
      mq=dsqrt(mQsq)

C----Initialize whole array to zero
      msq(:,:)=0d0

      fac=gsq**3*gw**4/4d0*32d0

      s56=2d0*mQsq
     & +2d0*(+p(5,4)*p(6,4)-p(5,1)*p(6,1)-p(5,2)*p(6,2)-p(5,3)*p(6,3))
      beta=sqrt(1d0-4d0*mQsq/s56)
      bp=0.5d0+0.5d0*beta
      bm=0.5d0-0.5d0*beta

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


      double precision function rodmsqm(j1,j2,j3,j4,j5,j6,j7,mQ)
      implicit none
c matrix element squared summed over colors and spins
      include 'constants.f'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,h3,h4,h5
      double complex qcda(2,2,2),qcdb(2,2,2),qedi(2,2,2),qedf(2,2,2)
      double precision prop,mQ
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
     & V*xn/eight*(cdabs(qcda(h3,h4,h5))**2+cdabs(qcdb(h3,h4,h5))**2)
     &+V/(eight*xn)*(cdabs(qedi(h3,h4,h5))**2+cdabs(qedf(h3,h4,h5))**2
     &-two*(cdabs(qedi(h3,h4,h5)+qedf(h3,h4,h5)))**2)
      enddo
      enddo
      enddo

      rodmsqm=rodmsqm/prop
      return
      end
