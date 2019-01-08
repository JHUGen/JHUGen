      subroutine gen_Vphotons_jets_dkrad2(r,nphots,njets,p,wt,*)
c---- generate phase space for 2-->2+nphots+njets process
c----   with (3+4+photon+photon) being a vector boson
c----   and 5,..,4+nphots-2 the photons (if nphots>2)
c----   and 5+nphots,..,4+nphots+njets the jets
c----
c---- note: as denotd above, the final two photons (3+nphots) and
c----       (4+nphots) are generated in the decay of the vector boson
c-----
c---- This routine is adapted from gen_njets.f
c----
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'limits.f'
      include 'xmin.f'
      include 'leptcuts.f'
      include 'reset.f'
      include 'breit.f'
      include 'part.f'
      include 'process.f'
      include 'ipsgen.f'
      include 'x1x2.f'
      include 'first.f'
      double precision r(mxdim)
      double precision p(mxpart,4),p3(4),psumjet(4),pcm(4),Q(4)
      double precision wt
      double precision hmin,hmax,delh,h,pt,etamax,etamin
      double precision y,sinhy,coshy,phi,mv2,wtbw,mjets
      double precision ybar,ptsumjet2,ycm,sumpst,q0st,rshat,dely
      double precision ptjetmin,etajetmin,etajetmax,pbreak
      double precision ptmin_part,etamax_part,pbreak_part
      double precision plstar,estar,plstarsq,y5starmax,y5starmin
      double precision p4(4),p5(4),p6(4),wt35,wt46,wt3456,p35(4),p46(4)
      double precision ymax,yave
      integer j,nu,nphots,njets,nphotsjets,ijet
      logical xxerror,flatreal
      include 'energy.f'
      parameter(flatreal=.false.)
      data xxerror/.false./
      save ptjetmin,etajetmin,etajetmax,pbreak,xxerror
!$omp threadprivate(ptjetmin,etajetmin,etajetmax,pbreak)
!$omp threadprivate(xxerror)
      
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
        if ((ptjetmin .lt. 5d0) .and. (part .ne. 'real')) pbreak=5d0
c--- for processes in which it is safe to jet ptmin to zero at NLO
      if ((part .eq. 'real') .and. (pbreak .lt. 1d-8)) pbreak=5d0
      endif        


c--- total number of photons and jets generated outside of boson decay
      nphotsjets=nphots-2+njets

      do nu=1,4
        do j=1,4+nphots+njets
          p(j,nu)=0d0
        enddo
        Q(nu)=0d0 
        psumjet(nu)=0d0
        pcm(nu)=0d0
      enddo 

      wt=2d0*pi
      
      if     (nphotsjets .lt. 0) then
        write(6,*) 'nphotsjets < 0 in gen_phots_jets_dkrad2.f'
      write(6,*) 'nphots = ',nphots
      write(6,*) 'njets  = ',njets
        stop
      elseif (nphotsjets .eq. 0) then
        goto 77 ! no loop to perform (nphots=2, njets=0)
      endif               
                   
      do ijet=1,nphotsjets
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
        wt=wt/16d0/pi**3
c        xmin=2d0/sqrts
c        xmax=1d0/ptjetmin

        if (ijet .le. nphots-2) then
          if     (nphots .eq. 3) then
            ptmin_part=min(gammpt,gammpt2,gammpt3)
          elseif (nphots .eq. 2) then
            ptmin_part=min(gammpt,gammpt2)
          elseif (nphots .eq. 1) then
            ptmin_part=gammpt
          else
            write(6,*) 'Unexpected # photons gen_Vphotons_jets_dkrad2:',
     &                  nphots
            stop
          endif
          etamax_part=gammrap
          pbreak_part=0d0
          if (part .eq. 'real') then
c--- cannot generate exactly to match, since dipoles transform photon
          ptmin_part=0d0
          etamax_part=20d0
          pbreak_part=gammpt
        endif
        else
          ptmin_part=ptjetmin
          etamax_part=etajetmax
          pbreak_part=pbreak
        endif

        if ((flatreal) .and. (part .eq. 'real')) then
c--- generate flat pt for the real contribution
          pt=r(ijet)*((sqrts/2d0)-ptmin_part)+ptmin_part
          wt=wt*((sqrts/2d0)-ptmin_part)*pt
        else
c--- favour small pt region 
          hmin=1d0/dsqrt((sqrts/2d0)**2+pbreak_part**2)
          hmax=1d0/dsqrt(ptmin_part**2+pbreak_part**2)
          delh=hmax-hmin
          h=hmin+r(ijet)*delh        
          pt=dsqrt(1d0/h**2-pbreak_part**2)
          wt=wt*delh/h**3
        endif

        etamax=sqrts/2d0/pt
        if (etamax**2 .le. 1d0) then
          write(6,*) 'etamax**2 .le. 1d0 in gen_phots_jets.f',etamax**2 
          wt=0d0
          return 1
        endif
        etamax=dlog(etamax+dsqrt(etamax**2-1d0))
        
        etamax=min(etamax,etamax_part)
        y=etamax*(2d0*r(nphotsjets+ijet)-1d0)
        wt=wt*2d0*etamax
        
        sinhy=dsinh(y)
        coshy=dsqrt(1d0+sinhy**2)
        
        p(6+ijet,4)=pt*coshy
        
        phi=2d0*pi*r(2*nphotsjets+ijet)
        wt=wt*2d0*pi
        
        p(6+ijet,1)=pt*dcos(phi)
        p(6+ijet,2)=pt*dsin(phi)
        p(6+ijet,3)=pt*sinhy
        
        do nu=1,4
          psumjet(nu)=psumjet(nu)+p(6+ijet,nu)
        enddo
      enddo
   
   77 continue
      
      call breitw(r(3*nphotsjets+1),wsqmin,sqrts**2,
     &  mass3,width3,mv2,wtbw)
      wt=wt*wtbw/2d0/pi

c--- catch special case
      if (nphotsjets .eq. 0) then
        rshat=dsqrt(mv2)
        ymax=dlog(sqrts/rshat)
        yave=ymax*(two*r(2)-1d0)     
        wt=wt*two*ymax/sqrts**2         
        xx(1)=rshat/sqrts*dexp(+yave)
        xx(2)=rshat/sqrts*dexp(-yave)
        Q(4)=(xx(1)+xx(2))/2d0*sqrts
        Q(3)=(xx(1)-xx(2))/2d0*sqrts
c        etamax=8d0
c        etamin=-8d0
        goto 89
      endif

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
      if (Estar/plstar-1d0 .lt. 1d-8) then
        y5starmax=20d0
      else
        y5starmax=0.5d0*dlog((Estar+plstar)/(Estar-plstar))
      endif
      y5starmin=-y5starmax

      etamax=ybar-y5starmin
      etamin=ybar-y5starmax

      dely=etamax-etamin
      ycm=etamin+r(3*nphotsjets+2)*dely     
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
        write(6,*) 'gen_njets: xx(1)*xx(2),xx(1),xx(2)',
     .   xx(1)*xx(2),xx(1),xx(2)  
      endif

      wt=wt*dely
      do j=1,4
        Q(j)=pcm(j)-psumjet(j)
      enddo
  
      wt=wt*rshat/(sqrts**2*q0st)
      
   89 continue 
      
      if   ((xx(1) .gt. 1d0) .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin).or. (xx(2) .lt. xmin)) then
         wt=0d0
         return 1
      endif 
      
      p(1,4)=-xx(1)*sqrts/2d0
      p(1,3)=p(1,4)
      p(2,4)=-xx(2)*sqrts/2d0
      p(2,3)=-p(2,4)
      
c--- now decay Q into 3,4,5,6
      call phi1_2(r(3*nphotsjets+3),r(3*nphotsjets+4),
     & r(3*nphotsjets+5),r(3*nphotsjets+6),Q,p35,p46,wt3456,*999)
      call phi3m0(r(3*nphotsjets+7),r(3*nphotsjets+8),
     & p35,p3,p5,wt35,*999)
      call phi3m0(r(3*nphotsjets+9),r(3*nphotsjets+10),
     & p46,p4,p6,wt46,*999)
      wt=wt*wt3456*wt35*wt46/twopi**2

      do nu=1,4
      if (r(3*nphotsjets+11) .lt. 0.5d0) then
        p(3,nu)=p3(nu)
        p(4,nu)=p4(nu)
      else
        p(3,nu)=p4(nu)
        p(4,nu)=p3(nu)
      endif
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      enddo

      return
      
  999 wt=0d0
      return 1     
      
      end
      
      
      
      
      
      
      
      
      
      
      
