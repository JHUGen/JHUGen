      subroutine gen4hvv(r,p,wt4,*)
      implicit none
      include 'types.f'
c--- Generates 2->4 phase space with (3456) a Breit-Wigner around
c--- the Higgs mass, and subsequent decays to either
c--- Z(34) Z(56) or W(54) W(36)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'debug.f'
      include 'masses.f'
      include 'phasemin.f'
      include 'breit.f'
      include 'interference.f'
      include 'ipsgen.f'
      include 'x1x2.f'
      include 'energy.f'
      integer:: nu
      real(dp):: r(mxdim)
      real(dp):: wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p(mxpart,4),rtshat
      real(dp):: pswt,xjac,s34,s56,s45,s36,wt_ww,wt_zz
      real(dp):: s3456,wt3456,ymax,yave

      wt4=0._dp

      call breitw(r(9),0._dp,sqrts**2,hmass,hwidth,s3456,wt3456)

      rtshat=sqrt(s3456)
      ymax=log(sqrts/rtshat)
      yave=ymax*(two*r(10)-1._dp)
      xjac=two*ymax*wt3456

      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)

      if   ((xx(1) > 1._dp)
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) then
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

      if (ipsgen == 1) then
c--- generating Z(34) Z(56)
        bw34_56=.true.
        mass2=zmass
        mass3=zmass
        width2=zwidth
        width3=zwidth
        call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*999)
      else
c--- generating W(54) W(36)
        bw34_56=.false.
        mass2=wmass
        mass3=wmass
        width2=wwidth
        width3=wwidth
        call phase4(r,p1,p2,p5,p4,p3,p6,pswt,*999)
      endif

      s34=2._dp*(p3(4)*p4(4)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3))
      s56=2._dp*(p5(4)*p6(4)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3))
      s36=2._dp*(p3(4)*p6(4)-p3(1)*p6(1)-p3(2)*p6(2)-p3(3)*p6(3))
      s45=2._dp*(p5(4)*p4(4)-p5(1)*p4(1)-p5(2)*p4(2)-p5(3)*p4(3))

c--- weighting to suppress poorly-sampled regions;
c--- wt_zz must also suppress possible photon pole for Z(34), in addition to BW
      wt_zz=(((s34-zmass**2)**2+(zmass*zwidth)**2)
     &      *((s56-zmass**2)**2+(zmass*zwidth)**2))
     &      *(s34/zmass**2)**2
      wt_ww=(((s36-wmass**2)**2+(wmass*wwidth)**2)
     &      *((s45-wmass**2)**2+(wmass*wwidth)**2))
      if (bw34_56) then
        xjac=xjac*(wt_ww/(wt_zz+wt_ww))
      else
        xjac=xjac*(wt_zz/(wt_zz+wt_ww))
      endif

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=0._dp
      enddo

      wt4=xjac*pswt/sqrts**2

      if (debug) write(6,*) 'wt4 in gen4hvv',wt4
      return

 999  return 1
      end

