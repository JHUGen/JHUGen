      subroutine gen7m(r,p,m3,m4,m5,wt3,*)
      implicit none
      include 'types.f'
c----generate 7 dimensional phase space weight and vectors p(mxpart,4)
c----           p1+p2+(p3+p4+p5)+(p6+p7+p8)+p9=0
c----  where (p3+p4+p5) has mass m3, (p6+p7+p8) has mass m4 and p9, m5
c----and x1 and x2 given nineteen random numbers

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'phasemin.f'
      include 'x1x2.f'
      integer:: nu,j
      real(dp):: r(mxdim),wt3
      real(dp):: p(mxpart,4),
     & p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),p9(4),m3,m4,m5
      real(dp):: pswt,xjac,tau,y

      include 'energy.f'

      do nu=1,4
      do j=1,mxpart
      p(j,nu)=0._dp
      enddo
      enddo

      wt3=0._dp
      tau=exp(log(taumin)*r(18))
      y=0.5_dp*log(tau)*(1._dp-2._dp*r(19))
      xjac=log(taumin)*tau*log(tau)

      xx(1)=sqrt(tau)*exp(+y)
      xx(2)=sqrt(tau)*exp(-y)

c---if x's out of normal range alternative return
      if   ((xx(1) > 1._dp)
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) return 1

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half

      if (kcase==ktottth) then
      m3=mt
      m4=mt
      m5=hmass
      endif

      call phase7m_alt(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,m3,m4,m5,pswt)

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=p7(nu)
      p(8,nu)=p8(nu)
      p(9,nu)=p9(nu)
      enddo
      wt3=xjac*pswt

      if (wt3 == 0._dp) return 1

      return
      end
