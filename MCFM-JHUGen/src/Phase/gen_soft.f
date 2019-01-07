      subroutine gen_soft(r,njets,p,wt,*)
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'limits.f'
      include 'xmin.f'
      include 'reset.f'
      include 'breit.f'
      include 'part.f'
      include 'x1x2.f'
      include 'first.f'
      include 'energy.f'
c---- generate phase space for 2-->2+n process
c---- with (34) being a vector boson and 5,..,4+n the jets
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
      double precision r(mxdim)
      double precision p(mxpart,4),p3(4),p34(4),psumjet(4),pcm(4),Q(4)
      double precision wt,tmp(4)
      double precision hmin,hmax,delh,h,pt,etamax,etamin
      double precision y,sinhy,coshy,phi,mv2,wtbw,mjets
      double precision ybar,ptsumjet2,ycm,sumpst,q0st,rshat
      double precision costh,sinth,dely
      double precision ptjetmin,etajetmin,etajetmax,pbreak
      double precision plstar,estar,plstarsq,y5starmax,y5starmin
      integer j,nu,njets,ijet
      logical gen15,gen16,gen17,gen25,gen26,gen27,gen56,gen57,gen67
      logical xxerror
      data xxerror/.false./
      save ptjetmin,etajetmin,etajetmax,pbreak,xxerror
!$omp threadprivate(ptjetmin,etajetmin,etajetmax,pbreak,xxerror)

      if (first .or. reset) then
        first=.false.
        reset=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
        if (part .eq. 'real') then
c--- if we're generating phase space for real emissions, then we need
c--- to produce partons spanning the whole phase space pt>0,eta<10;
c--- in this case, pbreak=ptjetmin simply means that we
c--- generate pt approx. 1/x for pt > pbreak and
c--- pt approx. uniformly for pt < pbreak
          pbreak=ptjetmin
          ptjetmin=0d0
          etajetmax=20d0
        else
c--- for lord and virt, the partons produced here can be generated
c--- right up to the jet cut boundaries and there is no need for pbreak
          pbreak=0d0
        endif
c--- in case this routine is used for very small values of ptjetmin
        if (ptjetmin .lt. 5d0) pbreak=5d0
      endif        

      do nu=1,4
        do j=1,4+njets
          p(j,nu)=0d0
        enddo
        psumjet(nu)=0d0
        pcm(nu)=0d0
      enddo 

      wt=2d0*pi
        
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
        wt=wt/16d0/pi**3
c        xmin=2d0/sqrts
c        xmax=1d0/ptjetmin
        hmin=1d0/dsqrt((sqrts/2d0)**2+pbreak**2)
        hmax=1d0/dsqrt(ptjetmin**2+pbreak**2)
        delh=hmax-hmin

        if (gen57 .and. (ijet .eq. njets)) then
c--- generate last jet collinear to jet 1
           r(ijet)=r(1)
           r(njets+ijet)=max(r(njets+1)-1d-4,0d0)
           r(2*njets+ijet)=max(r(2*njets+1)-1d-4,0d0)
        endif
        if (gen67 .and. (ijet .eq. njets)) then
c--- generate last jet collinear to jet 2
           r(ijet)=r(2)
           r(njets+ijet)=max(r(njets+2)-1d-4,0d0)
           r(2*njets+ijet)=max(r(2*njets+2)-1d-4,0d0)
        endif        
        if (gen56 .and. (ijet .eq. 2)) then
c--- generate second jet collinear to jet 1
           r(ijet)=r(1)
           r(njets+ijet)=max(r(njets+1)-1d-4,0d0)
           r(2*njets+ijet)=max(r(2*njets+1)-1d-4,0d0)
        endif

        h=hmin+r(ijet)*delh
c        if (ijet .eq. njets) then
c--- generate last jet with a very low pt
c          h=hmax-r(ijet)*hmin/1d6          
c        else
c--- generate other jets with a much higher pt
c          h=hmin+r(ijet)*(hmax-hmin*101d0)
c        endif
c        if (ijet .eq. 1) then
c--- generate first jet with a very low pt
c          h=hmax-r(ijet)*hmin/1d6          
c        else
c--- generate other jets with a much higher pt
c          h=hmin+r(ijet)*(hmax-hmin*101d0)
c        endif
c        if (ijet .eq. 2) then
c--- generate second jet with a very low pt
c          h=hmax-r(ijet)*hmin/1d6          
c        else
c--- generate other jets with a much higher pt
c          h=hmin+r(ijet)*(hmax-hmin*101d0)
c        endif
        pt=dsqrt(1d0/h**2-pbreak**2)
        etamax=sqrts/2d0/pt
        if (etamax**2 .le. 1d0) then
            write(6,*) 'etamax**2 .le. 1d0 in gen_njets.f',etamax**2 
            wt=0d0
            return 1
        endif
        etamax=dlog(etamax+dsqrt(etamax**2-1d0))
        
        etamax=min(etamax,etajetmax)
        y=etamax*(2d0*r(njets+ijet)-1d0)

        if     ((gen17 .or. gen15 .or. gen16)
     .    .and. (ijet .eq. njets)) then
c--- generate last jet collinear to beam 1
          y=+etamax*(1d0-r(njets+ijet)/1d3)       
        elseif ((gen27 .or. gen25 .or. gen26)
     .    .and. (ijet .eq. njets)) then
c--- generate last jet collinear to beam 2
          y=-etamax*(1d0-r(njets+ijet)/1d3)       
c        elseif (gen16 .and. (ijet .eq. 2)) then
c--- generate second jet collinear to beam 1
c          y=+etamax*(1d0-r(njets+ijet)/1d3)       
c        elseif (gen26 .and. (ijet .eq. 2)) then
c--- generate second jet collinear to beam 2
c          y=-etamax*(1d0-r(njets+ijet)/1d3) 
c        elseif (gen15 .and. (ijet .eq. 1)) then
c--- generate first jet collinear to beam 1
c          y=+etamax*(1d0-r(njets+ijet)/1d3)       
c        elseif (gen25 .and. (ijet .eq. 1)) then
c--- generate first jet collinear to beam 2
c          y=-etamax*(1d0-r(njets+ijet)/1d3) 
        else
          y=etamax*(2d0*r(njets+ijet)-1d0)/3d0     
        endif        
        wt=wt*2d0*etamax
        
        sinhy=dsinh(y)
        coshy=dsqrt(1d0+sinhy**2)
        
        p(4+ijet,4)=pt*coshy
        wt=wt*delh/h**3
        
        phi=2d0*pi*r(2*njets+ijet)
        wt=wt*2d0*pi
        
        p(4+ijet,1)=pt*dcos(phi)
        p(4+ijet,2)=pt*dsin(phi)
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
      wt=wt*wtbw/2d0/pi
c--- invariant mass of jets
      mjets=psumjet(4)**2-psumjet(1)**2-psumjet(2)**2-psumjet(3)**2
      mjets=dsqrt(dabs(mjets))
      
      ybar=0.5d0*dlog((psumjet(4)+psumjet(3))/(psumjet(4)-psumjet(3)))
      ptsumjet2=psumjet(1)**2+psumjet(2)**2
      plstarsq=((sqrts**2-mv2-mjets**2)**2
     . -4d0*(mjets**2*mv2+ptsumjet2*sqrts**2))/(4d0*sqrts**2)
      if (plstarsq .le. 0d0) then
        wt=0d0
        return 1
      endif
      plstar=dsqrt(plstarsq)
      Estar=dsqrt(plstarsq+ptsumjet2+mjets**2)
      y5starmax=0.5d0*dlog((Estar+plstar)/(Estar-plstar))
      y5starmin=-y5starmax

      etamax=ybar-y5starmin
      etamin=ybar-y5starmax
      dely=etamax-etamin
      ycm=etamin+r(3*njets+2)*dely     
      sinhy=dsinh(ycm)
      coshy=dsqrt(1d0+sinhy**2)
      
c--- now make the initial state momenta
      sumpst=ptsumjet2+(psumjet(3)*coshy-psumjet(4)*sinhy)**2
      q0st=dsqrt(mv2+sumpst)
      rshat=q0st+dsqrt(mjets**2+sumpst)
      pcm(4)=rshat*coshy
      pcm(3)=rshat*sinhy
            
      xx(1)=(pcm(4)+pcm(3))/sqrts
      xx(2)=(pcm(4)-pcm(3))/sqrts
      
      if   ((xx(1)*xx(2) .gt. 1d0) .and. (xxerror .eqv. .false.)) then
c        xxerror=.true.
c        write(6,*) 'gen_njets: xx(1)*xx(2),xx(1),xx(2)',
c     .   xx(1)*xx(2),xx(1),xx(2)  
        wt=0d0
        return 1
      endif

      if   ((xx(1) .gt. 1d0) .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin).or. (xx(2) .lt. xmin)) then
c         xxerror=.true.
c         write(6,*) 'gen_njets: xx(1),xx(2),xmin',xx(1),xx(2),xmin  
         wt=0d0
         return 1
      endif 
      
      wt=wt*dely
      do j=1,4
        Q(j)=pcm(j)-psumjet(j)
      enddo
      
      p(1,4)=-xx(1)*sqrts/2d0
      p(1,3)=p(1,4)
      p(2,4)=-xx(2)*sqrts/2d0
      p(2,3)=-p(2,4)
      
      wt=wt*rshat/(sqrts**2*q0st)
c--- decay boson into leptons, in boson rest frame
      costh=2d0*r(3*njets+3)-1d0
      sinth=dsqrt(1d0-costh**2)
      phi=2d0*pi*r(3*njets+4)
      p34(4)=dsqrt(mv2)/2d0
      p34(1)=p34(4)*sinth*dcos(phi)
      p34(2)=p34(4)*sinth*dsin(phi)
      p34(3)=p34(4)*costh
      
c--- boost into lab frame    
      call boost(dsqrt(mv2),Q,p34,p3)
      do j=1,4
      p(3,j)=p3(j)
      p(4,j)=Q(j)-p(3,j)
      enddo

      wt=wt/8d0/pi
            
      return
      end
      
      
      
      
      
      
      
      
      
      
      
