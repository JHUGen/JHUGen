      subroutine gen4vv(r,p,wt4,*)
      implicit none
      include 'types.f'
c--- Generates 2->4 phase space with either
c--- Z(34) Z(56) or W(54) W(36)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'debug.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'phasemin.f'
      include 'breit.f'
      include 'interference.f'
      include 'ipsgen.f'
      include 'x1x2.f'
      include 'energy.f'
      include 'nproc.f'
      integer:: nu,icount
      real(dp):: r(mxdim)
      real(dp):: wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p(mxpart,4)
      real(dp):: pswt,xjac,s34,s56,s45,s36,wt_ww,wt_zz
      real(dp):: tau,x1mx2,surd
      real(dp):: lntaum

      wt4=0._dp

      lntaum=log(taumin)
      tau=exp(lntaum*(one-r(9)))
      xjac=-lntaum*tau

      x1mx2=two*r(10)-one
      surd=sqrt(x1mx2**2+four*tau)

      xx(1)=half*(+x1mx2+surd)
      xx(2)=half*(-x1mx2+surd)

      xjac=xjac*two/surd

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

      wt4=xjac*pswt

      if (debug) write(6,*) 'wt4 in gen4vv',wt4
      return

 999  return 1
      end

