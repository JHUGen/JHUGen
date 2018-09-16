      subroutine genVHjtaucut(r,p,wt,*)
      implicit none
      include 'types.f'
C---generate two particle phase space and x1,x2 integration
C---p1+p2 --> p3+p4+p5+p6+p7
c----
c---- with p7 generated using tau as a variable of integration,
c---- with minimum value taucut
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
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
      include 'zerowidth.f'
      include 'npart.f'
      include 'kprocess.f'
      real(dp):: r(mxdim),p(mxpart,4),ph(mxpart,4),p3(4),p4(4),
     & p5(4),p6(4),p34(4),p56(4),p3456(4),wt34,wt3456,Qsqmin,Qsqmax,
     & p5678(4),p78(4),p7(4),p8(4),wt5678,wt56,wt78
      real(dp):: wbw,wt,wtdk,Qsq,rtshat,mass
      real(dp), parameter:: wt0=one/twopi**3
      real(dp), parameter:: Qsqmincut = 0.1_dp

      p(:,:)=zip

      wt=zip

      Qsqmin=max(m3456min**2,Qsqmincut)
      Qsqmax=sqrts**2
      if (zerowidth) then
        Qsqmin=max(Qsqmin,(n2*mass2+n3*mass3)**2)
      endif

c--- generate invariant mass of Q=p3+p4 
c--- linearly
c      Qsq=(Qsqmax-Qsqmin)*r(1)+Qsqmin
c      wbw=(Qsqmax-Qsqmin)
c--- logarithmically
      Qsq=Qsqmin*(Qsqmax/Qsqmin)**r(1)
      wbw=log(Qsqmax/Qsqmin)*Qsq
     
      rtshat=sqrt(Qsq)

c--- generate pa+pb -> Q
      call genQ(rtshat,r(2),taucut,p,wt)
      
c--- generate extra parton
c--- note that r(ndim+1) is a uniform random variable not adapted by VEGAS
      call genparton(2,p,r(3),r(4),r(5),r(ndim+1),taucut,ph,wt)

c--- decay Q into W+H
      p3456(:)=ph(3,:)
      mass2=hmass
      width2=hwidth
      if     (kcase==kWH1jet) then 
        mass3=wmass
        width3=wwidth
      elseif (kcase==kZH1jet) then
        mass3=zmass
        width3=zwidth
      else
        write(6,*) 'Unexpected kcase in genVHjtaucut.f: ',kcase
        stop
      endif
      call phi1_2(r(6),r(7),r(8),r(9),p3456,p56,p34,wt3456,*999)
      call phi3m0(r(10),r(11),p34,p3,p4,wt34,*999)
      
c--- translate to momenta to be returned
      p(1,:)=-ph(1,:)
      p(2,:)=-ph(2,:)
      p(3,:)=p3(:)
      p(4,:)=p4(:)
      
      if     (hdecaymode == 'bqba') then
        mass=mb
        call phi3m(r(12),r(13),p56,p5,p6,mass,mass,wtdk,*999)
        p(5,:)=p5(:)
        p(6,:)=p6(:)
        p(7,:)=ph(4,:)
        p(8,:)=zip
      elseif (hdecaymode == 'tlta') then
        mass=mtau
        call phi3m(r(12),r(13),p56,p5,p6,mass,mass,wtdk,*999)
        p(5,:)=p5(:)
        p(6,:)=p6(:)
        p(7,:)=ph(4,:)
        p(8,:)=zip
      elseif (hdecaymode == 'gaga') then
        mass=zip
        call phi3m0(r(12),r(13),p56,p5,p6,wtdk,*999)
        p(5,:)=p5(:)
        p(6,:)=p6(:)
        p(7,:)=ph(4,:)
        p(8,:)=zip
      elseif (hdecaymode == 'wpwm') then
         npart=7
         
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth
        p5678(:)=p56(:)
        call phi1_2(r(12),r(13),r(14),r(15),p5678,p56,p78,wt5678,*999)
        call phi3m0(r(16),r(17),p56,p5,p6,wt56,*999)
        call phi3m0(r(18),r(19),p78,p7,p8,wt78,*999)
        wtdk=wt5678*wt56*wt78/twopi**2
        p(5,:)=p5(:)
        p(6,:)=p6(:)
        p(7,:)=p7(:)
        p(8,:)=p8(:)
        p(9,:)=ph(4,:)
        p(10,:)=zip
      else
        write(6,*) 'genVHjtaucut: unexpected hdecaymode=',hdecaymode
        stop
      endif
      
      xx(1)=-two*p(1,4)/sqrts
      xx(2)=-two*p(2,4)/sqrts
      
      wt=wt0*wt3456*wt34*wtdk*wbw*wt
c--- the factor below will be added later, in lowint
      wt=wt*xx(1)*xx(2)*sqrts**2
      
c      call writeout(p)
c      pause
      
      if (p(3,4) .ne. p(3,4)) then
!        write(6,*) 'p(3,4) NaN'
        return 1
      endif
      
      if   ((xx(1) > 1._dp) 
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) then
!        write(6,*) 'problems with xx(1),xx(2) in genVHjtaucut',xx(1),xx(2)  
      return 1 
      endif
      
      return

  999 wt=zip
      return 1

      end
