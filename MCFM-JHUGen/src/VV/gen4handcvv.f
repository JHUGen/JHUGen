      subroutine gen4handcvv(r,p,wt4,*)
c--- Generates 2->4 phase space with (3456) a Breit-Wigner around
c--- the Higgs mass, and subsequent decays to either
c--- Z(34) Z(56) or W(54) W(36);  width of the BW is artificially enhanced to
c--- better map the continuum contribution
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'masses.f'
      include 'phasemin.f'
      include 'breit.f'
      include 'interference.f'
      include 'ipsgen.f'
      include 'x1x2.f'
      include 'energy.f'
      include 'first.f'
      integer nu
      double precision r(mxdim)
      double precision wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p(mxpart,4),rtshat
      double precision pswt,xjac,pswidth,s34,s56,s45,s36,wt_ww,wt_zz
      double precision s3456,wt3456,ymax,yave
      
      wt4=0d0

      if (first) then
c--- width to use in generation of PS: if too far from threshold, just
c--- use a width of 10 GeV in the B.W. to sample PS adequately
        if (hmass .lt. mass2+mass3-hwidth*5d0) then
          pswidth=10d0
        else
          pswidth=hwidth
        endif
      endif
      
      call breitw(r(9),0d0,sqrts**2,hmass,pswidth,s3456,wt3456)
            
      rtshat=dsqrt(s3456)
      ymax=dlog(sqrts/rtshat)
      yave=ymax*(two*r(10)-1d0)
      xjac=two*ymax*wt3456
           
      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
        return 1 
      endif

      p1(4)=-0.5d0*xx(1)*sqrts
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=-0.5d0*xx(1)*sqrts
      
      p2(4)=-0.5d0*xx(2)*sqrts
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=+0.5d0*xx(2)*sqrts

      if (ipsgen .eq. 1) then
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
      
      s34=2d0*(p3(4)*p4(4)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3))
      s56=2d0*(p5(4)*p6(4)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3))
      s36=2d0*(p3(4)*p6(4)-p3(1)*p6(1)-p3(2)*p6(2)-p3(3)*p6(3))
      s45=2d0*(p5(4)*p4(4)-p5(1)*p4(1)-p5(2)*p4(2)-p5(3)*p4(3))

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
      p(7,nu)=0d0
      enddo 
      
      wt4=xjac*pswt/sqrts**2
            
      if (debug) write(6,*) 'wt4 in gen4handcvv',wt4
      return

 999  return 1
      end

