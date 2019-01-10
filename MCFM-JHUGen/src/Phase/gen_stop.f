      subroutine gen_stop(r,njets,p,wt,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'limits.f'
      include 'xmin.f'
      include 'zerowidth.f'
      include 'process.f'
      include 'reset.f'
      include 'part.f'
      include 'x1x2.f'
      include 'energy.f'
      include 'notag.f'
      include 'first.f'
c---- Generate phase space for 2-->2+n process
c---- with (345) being a top and 6,..,5+n the jets
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(3n-4), where n is the number
c----  of final state particles)
c---- This routine has a minimum of 4 final state particles, hence
c---- the twopi**2 correction factor is given by the ratio of
c---- (1/twopi)**(3n-4) present in the phase space and the factor
c---- of [(1/twopi)**2]**(n-1) from the number of branchings
c---- For the specific case 'ttdkay' where one of the jets is
c---- associated with the top quark decay, we must add an extra
c---- factor of (1/twopi) since the number of jets generated is
c---- larger than the value of 'njets' passed 
      double precision r(mxdim)
      double precision p(mxpart,4),psumjet(4),pcm(4),Q(4)
      double precision wt,wt0,wtbg
      double precision hmin,hmax,delh,h,pt,etamax,etamin
      double precision y,sinhy,coshy,phi,mv2,wtbw,mjets
      double precision ybar,ptsumjet2,ycm,sumpst,q0st,rshat,dely
      double precision ptjetmin,etajetmin,etajetmax,pbreak
      double precision plstar,estar,plstarsq,y5starmax,y5starmin,mtrans
      double precision bm(4),wp(4),nn(4),ep(4),pbg(4),g(4),wtwp,wtepnn
      integer j,nu,njets,ijet,in
      logical oldzerowidth,xxerror
      parameter(wt0=1d0/twopi**2)
      data xxerror/.false./
      save ptjetmin,etajetmin,etajetmax,pbreak,xxerror
!$omp threadprivate(ptjetmin,etajetmin,etajetmax,pbreak,xxerror)

      if (first .or. reset) then
        first=.false.
        reset=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
        if (notag .eq. 1) then
c--- for the t-channel calculation, the default behaviour 
c--- in chooser.f is to set notag=1 for this process; in that case,
c--- the calculation is inclusive of all additional jets, so that
c--- partons must be generated without any explicit cuts on pt and eta
c--- necessary for t-channel process with notag=1
          pbreak=ptjetmin
          ptjetmin=0d0
          etajetmax=10d0
          if (pbreak .lt. 20d0) pbreak=20d0 ! in case ptjetmin not set
        else
c--- proceed as before (the old default)      
          if (part .eq. 'real') then
c---   if we're generating phase space for real emissions, then we need
c---   to produce partons spanning the whole phase space pt>0,eta<10;
c---   in this case, pbreak=ptjetmin simply means that we
c---   generate pt approx. 1/x for pt > pbreak and
c---   pt approx. uniformly for pt < pbreak
            pbreak=ptjetmin
            ptjetmin=0d0
            etajetmax=10d0
          else
c---   for lord and virt, the partons produced here can be generated
c---   right up to the jet cut boundaries and there is no need for pbreak
            pbreak=0d0
          endif
c--- in case this routine is used for very small values of ptjetmin
          if (ptjetmin .lt. 5d0) pbreak=5d0
      endif
      endif        

      do nu=1,4
        do j=1,5+njets
          p(j,nu)=0d0
        enddo
        psumjet(nu)=0d0
        pcm(nu)=0d0
      enddo 

      wt=2d0*pi
            
      do ijet=1,njets
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
        wt=wt/16d0/pi**3
c        xmin=2d0/sqrts
c        xmax=1d0/ptjetmin
        hmin=1d0/dsqrt((sqrts/2d0)**2+pbreak**2)
        hmax=1d0/dsqrt(ptjetmin**2+pbreak**2)
        delh=hmax-hmin
        h=hmin+r(ijet)*delh
        pt=dsqrt(1d0/h**2-pbreak**2)
        etamax=sqrts/2d0/pt
        if (etamax**2 .le. 1d0) then
            write(6,*) 'etamax**2 .le. 1d0 in gen_stop.f',etamax**2 
            wt=0d0
            return 1
        endif
        etamax=dlog(etamax+dsqrt(etamax**2-1d0))
        
        etamax=min(etamax,etajetmax)
        y=etamax*(2d0*r(njets+ijet)-1d0)
        wt=wt*2d0*etamax
        
        sinhy=dsinh(y)
        coshy=dsqrt(1d0+sinhy**2)
        
        p(5+ijet,4)=pt*coshy
        wt=wt*delh/h**3
        
        phi=2d0*pi*r(2*njets+ijet)
        wt=wt*2d0*pi
        
        p(5+ijet,1)=pt*dcos(phi)
        p(5+ijet,2)=pt*dsin(phi)
        p(5+ijet,3)=pt*sinhy

c--- generate a b-quark as particle 6 for s-channel processes
        if (((case .eq. 't_bbar') .or. (case .eq. 'tdecay'))
     &   .and. (ijet .eq. 1)) then
          mtrans=dsqrt(pt**2+mb**2)
          p(5+ijet,4)=mtrans*coshy
          p(5+ijet,3)=mtrans*sinhy      
      endif
        
        do nu=1,4
          psumjet(nu)=psumjet(nu)+p(5+ijet,nu)
        enddo
      enddo
      
c--- now generate Breit-Wigner, but always with zero width
      oldzerowidth=zerowidth
      zerowidth=.true.    
      call breitw(one,wsqmin,wsqmax,mt,twidth,mv2,wtbw)
      zerowidth=oldzerowidth
      wt=wt*wtbw
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
      ycm=etamin+r(3*njets+1)*dely     
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
        xxerror=.true.
        write(6,*) 'gen_stop: xx(1)*xx(2),xx(1),xx(2)',
     .   xx(1)*xx(2),xx(1),xx(2)  
      endif

      if   ((xx(1) .gt. 1d0) .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin).or. (xx(2) .lt. xmin)) then
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

c--- If we're calculating top decay then generate the additional jet
c--- for the real contribution here, after the decay      
      if ( ((case .eq. 'ttdkay') .or. (case .eq. 'tdecay'))
     .     .and. (part .eq. 'real') ) then
        in=3*njets+2
        call phi1_2(r(in),r(in+1),r(in+2),r(in+3),Q,pbg,wp,wtwp,*999)
        in=in+4
          call phi3m(r(in),r(in+1),pbg,bm,g,mb,zip,wtbg,*999)
          call phi3m0(r(in+2),r(in+3),wp,nn,ep,wtepnn,*999)
            wt=wt0*wt*wtwp*wtbg*wtepnn/twopi
        do nu=1,4
          p(7,nu)=g(nu)
        enddo
      else
        call phi1_2m(mb,r(3*njets+2),r(3*njets+3),r(3*njets+4),zip,
     .  Q,bm,wp,wtwp,*999)
        call phi3m0(r(3*njets+5),r(3*njets+6),wp,nn,ep,wtepnn,*999)
        wt=wt0*wt*wtwp*wtepnn
      endif

      do nu=1,4
        p(3,nu)=nn(nu)
        p(4,nu)=ep(nu)
        p(5,nu)=bm(nu)
      enddo
      
                             
      return
      
  999 wt=0d0
      return 1
      
      end
      
      
      
      
      
      
      
      
      
      
      
