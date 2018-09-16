      subroutine gen3taucut(r,p,wt,*)
      implicit none
      include 'types.f'
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4+p5
c----
c---- with p5 generated using tau as a variable of integration,
c---- with minimum value taucut
      
      include 'constants.f'
      include 'mxpart.f'
      include 'limits.f'
      include 'vegas_common.f'
      include 'phasemin.f'
      include 'breit.f'
      include 'x1x2.f'
      include 'taucut.f'
      include 'energy.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'debug.f'
      real(dp):: r(mxdim),p(mxpart,4),ph(mxpart,4),Q(4),p3(4),p4(4)
      real(dp):: wbw,wt,wtdk,Qsq,rtshat,mass,Qsqmin,Qsqmax
      real(dp), parameter:: wt0=one/twopi

      p(:,:)=zip

      wt=zip

c--- generate invariant mass of Q=p3+p4 
      if (n3==0) then
         wbw=one
         Qsqmin=max(wsqmin,one) ! ensure minimum value of m(34)>1 GeV
         Qsqmax=min(wsqmax,sqrts**2*0.999d0) ! ensure maximum value a bit below s
         call pick(2,Qsq,Qsqmin,Qsqmax,r(1),wbw)
      elseif (n3==1) then 
         call breitw(r(1),wsqmin,wsqmax,mass3,width3,Qsq,wbw)
      endif
      
      rtshat=sqrt(Qsq)

c--- generate pa+pb -> Q
      call genQ(rtshat,r(2),taucut,p,wt)
      
c      write(6,*) 'rtshat=',rtshat
c      write(6,*) 'p1',p(1,4),p(1,1),p(1,2),p(1,3)
c      write(6,*) 'p2',p(2,4),p(2,1),p(2,2),p(2,3)
c      write(6,*) 'p3',p(3,4),p(3,1),p(3,2),p(3,3)
c      write(6,*) 'p4',p(4,4),p(4,1),p(4,2),p(4,3)
      
c--- generate extra parton
c--- note that r(ndim+1) is a uniform random variable not adapted by VEGAS
      call genparton(2,p,r(3),r(4),r(5),r(ndim+1),taucut,ph,wt)

c--- decay Q
      if     (hdecaymode == 'bqba') then
        mass=mb
      elseif (hdecaymode == 'tlta') then
        mass=mtau
      else
        mass=zip
      endif
      Q(:)=ph(3,:)
      call phi3m(r(6),r(7),Q,p3,p4,mass,mass,wtdk,*999)

c      write(6,*) 'ph1',ph(1,4),ph(1,1),ph(1,2),ph(1,3)
c      write(6,*) 'ph2',ph(2,4),ph(2,1),ph(2,2),ph(2,3)
c      write(6,*) 'ph3',ph(3,4),ph(3,1),ph(3,2),ph(3,3)
c      write(6,*) 'ph4',ph(4,4),ph(4,1),ph(4,2),ph(4,3)
      
c--- translate to momenta to be returned
      p(1,:)=-ph(1,:)
      p(2,:)=-ph(2,:)
      p(3,:)=p3(:)
      p(4,:)=p4(:)
      p(5,:)=ph(4,:)
      p(6,:)=zip
      
      xx(1)=-two*p(1,4)/sqrts
      xx(2)=-two*p(2,4)/sqrts
      
      wt=wt0*wtdk*wbw*wt
c--- the factor below will be added later, in lowint
      wt=wt*xx(1)*xx(2)*sqrts**2
      
c      call writeout(p)
c      pause
      
      if   ((xx(1) > 1._dp) 
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) then
        if (debug) write(6,*) 'problems with xx(1),xx(2) in gen3taucut',
     & xx(1),xx(2)  
        return 1 
      endif
          
      return

  999 wt=zip
      return 1

      end
