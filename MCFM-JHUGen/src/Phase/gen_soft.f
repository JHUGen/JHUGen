      subroutine gen_soft(r,njets,p,wt,*)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'limits.f'
      include 'xmin.f'
      include 'reset.f'
      include 'breit.f'
      include 'kpart.f'
      include 'x1x2.f'
      include 'first.f'
      include 'energy.f'
c---- generate phase space for 2-->2+n process
c---- with (34) being a vector boson and 5,..,4+n the jets
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
      real(dp):: r(mxdim)
      real(dp):: p(mxpart,4),p3(4),p34(4),psumjet(4),pcm(4),Q(4)
      real(dp):: wt,tmp(4)
      real(dp):: hmin,hmax,delh,h,pt,etamax,etamin
      real(dp):: y,sinhy,coshy,phi,mv2,wtbw,mjets
      real(dp):: ybar,ptsumjet2,ycm,sumpst,q0st,rshat
      real(dp):: costh,sinth,dely
      real(dp):: ptjetmin,etajetmin,etajetmax,pbreak
      real(dp):: plstar,estar,plstarsq,y5starmax,y5starmin
      integer:: j,nu,njets,ijet
      logical:: gen15,gen16,gen17,gen25,gen26,gen27,gen56,gen57,gen67
      logical:: xxerror
      data xxerror/.false./
      save ptjetmin,etajetmin,etajetmax,pbreak,xxerror
!$omp threadprivate(ptjetmin,etajetmin,etajetmax,pbreak,xxerror)

      if (first .or. reset) then
        first=.false.
        reset=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
        if (kpart==kreal) then
c--- if we're generating phase space for real emissions, then we need
c--- to produce partons spanning the whole phase space pt>0,eta<10;
c--- in this case, pbreak=ptjetmin simply means that we
c--- generate pt approx. 1/x for pt > pbreak and
c--- pt approx. uniformly for pt < pbreak
          pbreak=ptjetmin
          ptjetmin=zip
          etajetmax=20._dp
        else
c--- for lord and virt, the partons produced here can be generated
c--- right up to the jet cut boundaries and there is no need for pbreak
          pbreak=zip
        endif
c--- in case this routine is used for very small values of ptjetmin
        if (ptjetmin < 5._dp) pbreak=5._dp
      endif        

      do nu=1,4
        do j=1,4+njets
          p(j,nu)=zip
        enddo
        psumjet(nu)=zip
        pcm(nu)=zip
      enddo 

      wt=2._dp*pi
        
      gen15=.false.  
      gen16=.false.  
      gen17=.false.  
      gen25=.true.  
      gen26=.false.  
      gen27=.false.  
      gen56=.false.  
      gen57=.false.  
      gen67=.false.  
            
      do ijet=1,njets
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
        wt=wt/16._dp/pi**3
c        xmin=2._dp/sqrts
c        xmax=1._dp/ptjetmin
        hmin=1._dp/sqrt((sqrts/2._dp)**2+pbreak**2)
        hmax=1._dp/sqrt(ptjetmin**2+pbreak**2)
        delh=hmax-hmin

        if (gen57 .and. (ijet == njets)) then
c--- generate last jet collinear to jet 1
           r(ijet)=r(1)
           r(njets+ijet)=max(r(njets+1)-1.e-4_dp,zip)
           r(2*njets+ijet)=max(r(2*njets+1)-1.e-4_dp,zip)
        endif
        if (gen67 .and. (ijet == njets)) then
c--- generate last jet collinear to jet 2
           r(ijet)=r(2)
           r(njets+ijet)=max(r(njets+2)-1.e-4_dp,zip)
           r(2*njets+ijet)=max(r(2*njets+2)-1.e-4_dp,zip)
        endif        
        if (gen56 .and. (ijet == 2)) then
c--- generate second jet collinear to jet 1
           r(ijet)=r(1)
           r(njets+ijet)=max(r(njets+1)-1.e-4_dp,zip)
           r(2*njets+ijet)=max(r(2*njets+1)-1.e-4_dp,zip)
        endif

        h=hmin+r(ijet)*delh
c        if (ijet == njets) then
c--- generate last jet with a very low pt
c          h=hmax-r(ijet)*hmin/1d6          
c        else
c--- generate other jets with a much higher pt
c          h=hmin+r(ijet)*(hmax-hmin*101._dp)
c        endif
c        if (ijet == 1) then
c--- generate first jet with a very low pt
c          h=hmax-r(ijet)*hmin/1d6          
c        else
c--- generate other jets with a much higher pt
c          h=hmin+r(ijet)*(hmax-hmin*101._dp)
c        endif
c        if (ijet == 2) then
c--- generate second jet with a very low pt
c          h=hmax-r(ijet)*hmin/1d6          
c        else
c--- generate other jets with a much higher pt
c          h=hmin+r(ijet)*(hmax-hmin*101._dp)
c        endif
        pt=sqrt(1._dp/h**2-pbreak**2)
        etamax=sqrts/2._dp/pt
        if (etamax**2 <= 1._dp) then
            write(6,*) 'etamax**2 <= 1._dp in gen_njets.f',etamax**2 
            wt=zip
            return 1
        endif
        etamax=log(etamax+sqrt(etamax**2-1._dp))
        
        etamax=min(etamax,etajetmax)
        y=etamax*(2._dp*r(njets+ijet)-1._dp)

        if     ((gen17 .or. gen15 .or. gen16)
     &    .and. (ijet == njets)) then
c--- generate last jet collinear to beam 1
          y=+etamax*(1._dp-r(njets+ijet)/1d3)       
        elseif ((gen27 .or. gen25 .or. gen26)
     &    .and. (ijet == njets)) then
c--- generate last jet collinear to beam 2
          y=-etamax*(1._dp-r(njets+ijet)/1d3)       
c        elseif (gen16 .and. (ijet == 2)) then
c--- generate second jet collinear to beam 1
c          y=+etamax*(1._dp-r(njets+ijet)/1d3)       
c        elseif (gen26 .and. (ijet == 2)) then
c--- generate second jet collinear to beam 2
c          y=-etamax*(1._dp-r(njets+ijet)/1d3) 
c        elseif (gen15 .and. (ijet == 1)) then
c--- generate first jet collinear to beam 1
c          y=+etamax*(1._dp-r(njets+ijet)/1d3)       
c        elseif (gen25 .and. (ijet == 1)) then
c--- generate first jet collinear to beam 2
c          y=-etamax*(1._dp-r(njets+ijet)/1d3) 
        else
          y=etamax*(2._dp*r(njets+ijet)-1._dp)/3._dp     
        endif        
        wt=wt*2._dp*etamax
        
        sinhy=sinh(y)
        coshy=sqrt(1._dp+sinhy**2)
        
        p(4+ijet,4)=pt*coshy
        wt=wt*delh/h**3
        
        phi=2._dp*pi*r(2*njets+ijet)
        wt=wt*2._dp*pi
        
        p(4+ijet,1)=pt*cos(phi)
        p(4+ijet,2)=pt*sin(phi)
        p(4+ijet,3)=pt*sinhy
        
        do nu=1,4
          psumjet(nu)=psumjet(nu)+p(4+ijet,nu)
        enddo
      enddo
      
c--- swap jets 5 and 7
      if (gen15 .or. gen25) then
        do nu=1,4
          tmp(nu)=p(4+njets,nu)
          p(4+njets,nu)=p(5,nu)
          p(5,nu)=tmp(nu)
        enddo
      endif
      
c--- swap jets 6 and 7      
      if (gen16 .or. gen26) then
        do nu=1,4
          tmp(nu)=p(4+njets,nu)
          p(4+njets,nu)=p(6,nu)
          p(6,nu)=tmp(nu)
        enddo
      endif
      
c--- now generate Breit-Wigner        
      call breitw(r(3*njets+1),wsqmin,wsqmax,mass3,width3,mv2,wtbw)
      wt=wt*wtbw/2._dp/pi
c--- invariant mass of jets
      mjets=psumjet(4)**2-psumjet(1)**2-psumjet(2)**2-psumjet(3)**2
      mjets=sqrt(abs(mjets))
      
      ybar=0.5_dp*log((psumjet(4)+psumjet(3))/(psumjet(4)-psumjet(3)))
      ptsumjet2=psumjet(1)**2+psumjet(2)**2
      plstarsq=((sqrts**2-mv2-mjets**2)**2
     & -4._dp*(mjets**2*mv2+ptsumjet2*sqrts**2))/(4._dp*sqrts**2)
      if (plstarsq <= zip) then
        wt=zip
        return 1
      endif
      plstar=sqrt(plstarsq)
      Estar=sqrt(plstarsq+ptsumjet2+mjets**2)
      y5starmax=0.5_dp*log((Estar+plstar)/(Estar-plstar))
      y5starmin=-y5starmax

      etamax=ybar-y5starmin
      etamin=ybar-y5starmax
      dely=etamax-etamin
      ycm=etamin+r(3*njets+2)*dely     
      sinhy=sinh(ycm)
      coshy=sqrt(1._dp+sinhy**2)
      
c--- now make the initial state momenta
      sumpst=ptsumjet2+(psumjet(3)*coshy-psumjet(4)*sinhy)**2
      q0st=sqrt(mv2+sumpst)
      rshat=q0st+sqrt(mjets**2+sumpst)
      pcm(4)=rshat*coshy
      pcm(3)=rshat*sinhy
            
      xx(1)=(pcm(4)+pcm(3))/sqrts
      xx(2)=(pcm(4)-pcm(3))/sqrts
      
      if   ((xx(1)*xx(2) > 1._dp) .and. (xxerror .eqv. .false.)) then
c        xxerror=.true.
c        write(6,*) 'gen_njets: xx(1)*xx(2),xx(1),xx(2)',
c     &   xx(1)*xx(2),xx(1),xx(2)  
        wt=zip
        return 1
      endif

      if   ((xx(1) > 1._dp) .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin).or. (xx(2) < xmin)) then
c         xxerror=.true.
c         write(6,*) 'gen_njets: xx(1),xx(2),xmin',xx(1),xx(2),xmin  
         wt=zip
         return 1
      endif 
      
      wt=wt*dely
      do j=1,4
        Q(j)=pcm(j)-psumjet(j)
      enddo
      
      p(1,4)=-xx(1)*sqrts/2._dp
      p(1,3)=p(1,4)
      p(2,4)=-xx(2)*sqrts/2._dp
      p(2,3)=-p(2,4)
      
      wt=wt*rshat/(sqrts**2*q0st)
c--- decay boson into leptons, in boson rest frame
      costh=2._dp*r(3*njets+3)-1._dp
      sinth=sqrt(1._dp-costh**2)
      phi=2._dp*pi*r(3*njets+4)
      p34(4)=sqrt(mv2)/2._dp
      p34(1)=p34(4)*sinth*cos(phi)
      p34(2)=p34(4)*sinth*sin(phi)
      p34(3)=p34(4)*costh
      
c--- boost into lab frame    
      call boost(sqrt(mv2),Q,p34,p3)
      do j=1,4
      p(3,j)=p3(j)
      p(4,j)=Q(j)-p(3,j)
      enddo

      wt=wt/8._dp/pi
            
      return
      end
      
      
      
      
      
      
      
      
      
      
      
