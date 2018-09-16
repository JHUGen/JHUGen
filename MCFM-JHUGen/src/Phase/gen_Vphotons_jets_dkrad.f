      subroutine gen_Vphotons_jets_dkrad(r,nphots,njets,p,wt,*)
      implicit none
      include 'types.f'
c---- generate phase space for 2-->2+nphots+njets process
c----   with (3+4+photon) being a vector boson
c----   and 5,..,4+nphots-1 the other photons (if nphots>1)
c----   and 5+nphots,..,4+nphots+njets the jets
c----
c---- note: as denotd above, the final photon (4+nphots) is generated 
c----       in the decay of the vector boson
c-----
c---- This routine is adapted from gen_njets.f
c----
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
c----
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mxdim.f'
      include 'limits.f'
      include 'xmin.f'
      include 'leptcuts.f'
      include 'reset.f'
      include 'breit.f'
      include 'kpart.f'
      include 'kprocess.f'
      include 'ipsgen.f'
      include 'x1x2.f'
      include 'first.f'
      real(dp):: r(mxdim)
      real(dp):: p(mxpart,4),p3(4),psumjet(4),pcm(4),Q(4)
      real(dp):: wt,dot,s346
      real(dp):: hmin,hmax,delh,h,pt,etamax,etamin
      real(dp):: y,sinhy,coshy,phi,mv2,wtbw,mjets
      real(dp):: ybar,ptsumjet2,ycm,sumpst,q0st,rshat,dely
      real(dp):: ptjetmin,etajetmin,etajetmax,pbreak
      real(dp):: ptmin_part,etamax_part,pbreak_part
      real(dp):: plstar,estar,plstarsq,y5starmax,y5starmin
      real(dp):: p45(4),p4(4),p5(4),wt45,wt345
      integer:: j,nu,nphots,njets,nphotsjets,ijet
      logical:: xxerror,flatreal
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
        if (kpart==kreal) then
c--- if we're generating phase space for real emissions, then we need
c--- to produce partons spanning the whole phase space pt>0,eta<10;
c--- in this case, pbreak=ptjetmin simply means that we
c--- generate pt approx. 1/x for pt > pbreak and
c--- pt approx. uniformly for pt < pbreak
          pbreak=ptjetmin
          ptjetmin=0._dp
          etajetmax=20._dp
        else
c--- for lord and virt, the partons produced here can be generated
c--- right up to the jet cut boundaries and there is no need for pbreak
          pbreak=0._dp
        endif
c--- in case this routine is used for very small values of ptjetmin
        if ((ptjetmin < 5._dp) .and. (kpart.ne.kreal)) pbreak=5._dp
c--- for processes in which it is safe to jet ptmin to zero at NLO
      if ((kpart==kreal) .and. (pbreak < 1.e-8_dp)) pbreak=5._dp
      endif        


c--- total number of photons and jets generated outside of boson decay
      nphotsjets=nphots-1+njets

      do nu=1,4
        do j=1,4+nphotsjets
          p(j,nu)=0._dp
        enddo
        Q(nu)=0._dp 
        psumjet(nu)=0._dp
        pcm(nu)=0._dp
      enddo 

      wt=2._dp*pi
                       
      do ijet=1,nphotsjets
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
        wt=wt/16._dp/pi**3
c        xmin=2._dp/sqrts
c        xmax=1._dp/ptjetmin

        if (ijet <= nphots-1) then
          if     (nphots == 3) then
            ptmin_part=min(gammpt,gammpt2,gammpt3)
          elseif (nphots == 2) then
            ptmin_part=min(gammpt,gammpt2)
          elseif (nphots == 1) then
            ptmin_part=gammpt
          else
            write(6,*) 'Unexpected # photons gen_Vphotons_jets_dkrad: ',
     &                  nphots
            stop
          endif
          etamax_part=gammrap
          pbreak_part=0._dp
          if (kpart==kreal) then
c--- cannot generate exactly to match, since dipoles transform photon
          ptmin_part=0._dp
          etamax_part=20._dp
          pbreak_part=gammpt
          endif
        else
          ptmin_part=ptjetmin
          etamax_part=etajetmax
          pbreak_part=pbreak
        endif

        if ((flatreal) .and. (kpart==kreal)) then
c--- generate flat pt for the real contribution
          pt=r(ijet)*((sqrts/2._dp)-ptmin_part)+ptmin_part
          wt=wt*((sqrts/2._dp)-ptmin_part)*pt
        else
c--- favour small pt region 
          hmin=1._dp/sqrt((sqrts/2._dp)**2+pbreak_part**2)
          hmax=1._dp/sqrt(ptmin_part**2+pbreak_part**2)
          delh=hmax-hmin
          h=hmin+r(ijet)*delh        
          pt=sqrt(1._dp/h**2-pbreak_part**2)
          wt=wt*delh/h**3
        endif

        etamax=sqrts/2._dp/pt
        if (etamax**2 <= 1._dp) then
          write(6,*) 'etamax**2 <= 1._dp in gen_phots_jets.f',etamax**2 
          wt=0._dp
          return 1
        endif
        etamax=log(etamax+sqrt(etamax**2-1._dp))
        
        etamax=min(etamax,etamax_part)
        y=etamax*(2._dp*r(nphotsjets+ijet)-1._dp)
        wt=wt*2._dp*etamax
        
        sinhy=sinh(y)
        coshy=sqrt(1._dp+sinhy**2)
        
        p(5+ijet,4)=pt*coshy
        
        phi=2._dp*pi*r(2*nphotsjets+ijet)
        wt=wt*2._dp*pi
        
        p(5+ijet,1)=pt*cos(phi)
        p(5+ijet,2)=pt*sin(phi)
        p(5+ijet,3)=pt*sinhy
        
        do nu=1,4
          psumjet(nu)=psumjet(nu)+p(5+ijet,nu)
        enddo
      enddo
      
c--- now generate Breit-Wigner        
      call breitw(r(3*nphotsjets+1),wsqmin,sqrts**2,
     & mass3,width3,mv2,wtbw)
      wt=wt*wtbw/2._dp/pi
c--- invariant mass of jets
      mjets=psumjet(4)**2-psumjet(1)**2-psumjet(2)**2-psumjet(3)**2
      mjets=sqrt(abs(mjets))
      
      ybar=0.5_dp*log((psumjet(4)+psumjet(3))/(psumjet(4)-psumjet(3)))
      ptsumjet2=psumjet(1)**2+psumjet(2)**2
      plstarsq=((sqrts**2-mv2-mjets**2)**2
     & -4._dp*(mjets**2*mv2+ptsumjet2*sqrts**2))/(4._dp*sqrts**2)
      if (plstarsq <= 0._dp) then
        wt=0._dp
        return 1
      endif
      plstar=sqrt(plstarsq)
      Estar=sqrt(plstarsq+ptsumjet2+mjets**2)
      y5starmax=0.5_dp*log((Estar+plstar)/(Estar-plstar))
      y5starmin=-y5starmax

      etamax=ybar-y5starmin
      etamin=ybar-y5starmax
      dely=etamax-etamin
      ycm=etamin+r(3*nphotsjets+2)*dely     
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
        xxerror=.true.
        write(6,*) 'gen_njets: xx(1)*xx(2),xx(1),xx(2)',
     &   xx(1)*xx(2),xx(1),xx(2)  
      endif

      if   ((xx(1) > 1._dp) .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin).or. (xx(2) < xmin)) then
         wt=0._dp
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
      
c--- now decay Q into 3,4,5
      call phi1_2m(zip,r(3*nphotsjets+3),r(3*nphotsjets+4),
     & r(3*nphotsjets+5),zip,Q,p3,p45,wt345,*999)     
      call phi3m0(r(3*nphotsjets+6),r(3*nphotsjets+7),
     & p45,p4,p5,wt45,*999)
      wt=wt*wt345*wt45/twopi     

      do nu=1,4
      if (r(3*nphotsjets+8) < 0.5_dp) then
        p(3,nu)=p3(nu)
        p(4,nu)=p4(nu)
      else
        p(3,nu)=p4(nu)
        p(4,nu)=p3(nu)
      endif
      p(5,nu)=p5(nu)
      enddo
            
c--- veto PS regions for W+2 photons and ipsgen=3
c      if (kcase==kW_2gam) then
c        if (ipsgen == 3) then
c          s346=2._dp*(dot(p,3,4)+dot(p,3,6)+dot(p,4,6))
c          if (abs(sqrt(s346)-mass3) < 5._dp*width3) return 1
c        endif
c      endif

      return
      
  999 wt=0._dp
      return 1     
      
      end
      
      
      
      
      
      
      
      
      
      
      
