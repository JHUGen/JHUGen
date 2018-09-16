      subroutine genVjtaucut(r,p,wt,*)
      implicit none
      include 'types.f'
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4+p5+p6
c----
c---- with p6 generated using tau as a variable of integration,
c---- with minimum value taucut
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'limits.f'
      include 'vegas_common.f'
      include 'phasemin.f'
      include 'nodecay.f'
      include 'breit.f'
      include 'x1x2.f'
      include 'taucut.f'
      include 'energy.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'debug.f'
      logical pass
      real(dp)::r(mxdim),p(mxpart,4),ph(mxpart,4),p3(4),p4(4),wt,t0
      real(dp), parameter:: wt0=one

      p(:,:)=zip

c--- generate p1+p2 -> p3+p4+p5
      call gen3(r,p,wt,*999)
!      call gen_njets(r,1,p,wt,*999)
c--- flip momentum convention
      p(1:2,:)=-p(1:2,:)
      
!      t0=taucut
      t0=1.e-12_dp
!      t0=1.e-8_dp
      
c--- generate extra parton
c--- note that r(ndim+1) is a uniform random variable not adapted by VEGAS
      call genpartonVj(p,r(8),r(9),r(10),r(ndim+1),t0,ph,wt,pass)
      if (pass .eqv. .false.) goto 999

c--- need to Lorentz boost p(3,:) and p(4,:) from
c--- p(3,:)+p(4,:) to ph(3,:)
      call boostx(p(3,:),p(3,:)+p(4,:),ph(3,:),p3)
      call boostx(p(4,:),p(3,:)+p(4,:),ph(3,:),p4)
      
c--- translate to momenta to be returned
      p(1,:)=-ph(1,:)
      p(2,:)=-ph(2,:)
      p(3,:)=p3(:)
      p(4,:)=p4(:)
      p(5,:)=ph(4,:)
      p(6,:)=ph(5,:)
      p(7,:)=zip
      
c--- symmetrize phase space
      if (rand() > 0.5_dp) then
        p(5,:)=ph(5,:)
        p(6,:)=ph(4,:)
      endif
      
c      call writeout(p)
c      pause
      
      wt=wt*wt0
      
      xx(1)=-two*p(1,4)/sqrts
      xx(2)=-two*p(2,4)/sqrts
      
c      call writeout(p)
c      pause
      
      if   ((xx(1) > 1._dp) 
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) then
      if (debug) write(6,*) 'problems with xx(1),xx(2) in genVjtaucut',
     & xx(1),xx(2)  
      return 1 
      endif
          
      return

  999 wt=zip
      return 1

      end
