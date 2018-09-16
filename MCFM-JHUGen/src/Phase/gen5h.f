c--- Generates 2->4 phase space with (12) a Breit-Wigner around
c--- the Higgs mass
      subroutine gen5h(r,p,wt5,*)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'debug.f'
      include 'masses.f'
      include 'phasemin.f'
      include 'kprocess.f'
      include 'breit.f'
      include 'x1x2.f'
      integer:: nu
      real(dp):: r(mxdim)
      real(dp):: wt5,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: p(mxpart,4),rtshat
      real(dp):: pswt,xjac
      real(dp):: s12,wt12,ymax,yave
      include 'energy.f'

      wt5=0._dp

      if (kcase.ne.kHZZqgI) then
c--- this is the usual case
        call breitw(r(9),0._dp,sqrts**2,hmass,hwidth,s12,wt12)
      else
c--- this is the HZZqgI case (cf. gen4handc.f)
        if (hmass < mass2+mass3-hwidth*5._dp) then
          call breitw(r(9),0._dp,sqrts**2,hmass,10._dp,s12,wt12)
        else
          call breitw(r(9),0._dp,sqrts**2,hmass,hwidth,s12,wt12)
        endif
      endif

      rtshat=sqrt(s12)
      ymax=log(sqrts/rtshat)
      yave=ymax*(two*r(10)-1._dp)
      xjac=two*ymax*wt12

      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)

      if   ((xx(1) > 1._dp)
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) then
c      write(6,*) 'problems with xx(1),xx(2) in gen4h',xx(1),xx(2)
      return 1
      endif

      p1(4)=-0.5_dp*xx(1)*sqrts
      p1(1)=0._dp
      p1(2)=0._dp
      p1(3)=-0.5_dp*xx(1)*sqrts

      p2(4)=-0.5_dp*xx(2)*sqrts
      p2(1)=0._dp
      p2(2)=0._dp
      p2(3)=+0.5_dp*xx(2)*sqrts

      call phase5h(r,p1,p2,p3,p4,p5,p6,p7,pswt)

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=p7(nu)
      enddo

      wt5=xjac*pswt/sqrts**2

      if (debug) write(6,*) 'wt5 in gen5h',wt5
      return

c 999  return 1
      end

