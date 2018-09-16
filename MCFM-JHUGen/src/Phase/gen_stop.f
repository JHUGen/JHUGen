      subroutine gen_stop(r,njets,p,wt,*)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'mxdim.f'
      include 'limits.f'
      include 'xmin.f'
      include 'zerowidth.f'
      include 'kprocess.f'
      include 'reset.f'
      include 'kpart.f'
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
      real(dp):: r(mxdim)
      real(dp):: p(mxpart,4),psumjet(4),pcm(4),Q(4)
      real(dp):: wt,wt0,wtbg
      real(dp):: hmin,hmax,delh,h,pt,etamax,etamin
      real(dp):: y,sinhy,coshy,phi,mv2,wtbw,mjets
      real(dp):: ybar,ptsumjet2,ycm,sumpst,q0st,rshat,dely
      real(dp):: ptjetmin,etajetmin,etajetmax,pbreak
      real(dp):: plstar,estar,plstarsq,y5starmax,y5starmin,mtrans
      real(dp):: bm(4),wp(4),nn(4),ep(4),pbg(4),g(4),wtwp,wtepnn
      integer:: j,nu,njets,ijet,in
      logical:: oldzerowidth,xxerror
      parameter(wt0=1._dp/twopi**2)
      data xxerror/.false./
      save ptjetmin,etajetmin,etajetmax,pbreak,xxerror
!$omp threadprivate(ptjetmin,etajetmin,etajetmax,pbreak,xxerror)

      if (first .or. reset) then
        first=.false.
        reset=.false.
        call read_jetcuts(ptjetmin,etajetmin,etajetmax)
        if (notag == 1) then
c--- for the t-channel calculation, the default behaviour 
c--- in chooser.f is to set notag=1 for this process; in that case,
c--- the calculation is inclusive of all additional jets, so that
c--- partons must be generated without any explicit cuts on pt and eta
c--- necessary for t-channel process with notag=1
          pbreak=ptjetmin
          ptjetmin=0._dp
          etajetmax=10._dp
          if (pbreak < 20._dp) pbreak=20._dp ! in case ptjetmin not set
        else
c--- proceed as before (the old default)      
          if (kpart==kreal) then
c---   if we're generating phase space for real emissions, then we need
c---   to produce partons spanning the whole phase space pt>0,eta<10;
c---   in this case, pbreak=ptjetmin simply means that we
c---   generate pt approx. 1/x for pt > pbreak and
c---   pt approx. uniformly for pt < pbreak
            pbreak=ptjetmin
            ptjetmin=0._dp
            etajetmax=10._dp
          else
c---   for lord and virt, the partons produced here can be generated
c---   right up to the jet cut boundaries and there is no need for pbreak
            pbreak=0._dp
          endif
c--- in case this routine is used for very small values of ptjetmin
          if (ptjetmin < 5._dp) pbreak=5._dp
      endif
      endif        

      do nu=1,4
        do j=1,5+njets
          p(j,nu)=0._dp
        enddo
        psumjet(nu)=0._dp
        pcm(nu)=0._dp
      enddo 

      wt=2._dp*pi
            
      do ijet=1,njets
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
        wt=wt/16._dp/pi**3
c        xmin=2._dp/sqrts
c        xmax=1._dp/ptjetmin
        hmin=1._dp/sqrt((sqrts/2._dp)**2+pbreak**2)
        hmax=1._dp/sqrt(ptjetmin**2+pbreak**2)
        delh=hmax-hmin
        h=hmin+r(ijet)*delh
        pt=sqrt(1._dp/h**2-pbreak**2)
        etamax=sqrts/2._dp/pt
        if (etamax**2 <= 1._dp) then
            write(6,*) 'etamax**2 <= 1._dp in gen_stop.f',etamax**2 
            wt=0._dp
            return 1
        endif
        etamax=log(etamax+sqrt(etamax**2-1._dp))
        
        etamax=min(etamax,etajetmax)
        y=etamax*(2._dp*r(njets+ijet)-1._dp)
        wt=wt*2._dp*etamax
        
        sinhy=sinh(y)
        coshy=sqrt(1._dp+sinhy**2)
        
        p(5+ijet,4)=pt*coshy
        wt=wt*delh/h**3
        
        phi=2._dp*pi*r(2*njets+ijet)
        wt=wt*2._dp*pi
        
        p(5+ijet,1)=pt*cos(phi)
        p(5+ijet,2)=pt*sin(phi)
        p(5+ijet,3)=pt*sinhy

c--- generate a b-quark as particle 6 for s-channel processes
        if (((kcase==kt_bbar) .or. (kcase==ktdecay))
     &   .and. (ijet == 1)) then
          mtrans=sqrt(pt**2+mb**2)
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
      ycm=etamin+r(3*njets+1)*dely     
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
        write(6,*) 'gen_stop: xx(1)*xx(2),xx(1),xx(2)',
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

c--- If we're calculating top decay then generate the additional jet
c--- for the real contribution here, after the decay      
      if ( ((kcase==kttdkay) .or. (kcase==ktdecay))
     &     .and. (kpart==kreal) ) then
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
     &  Q,bm,wp,wtwp,*999)
        call phi3m0(r(3*njets+5),r(3*njets+6),wp,nn,ep,wtepnn,*999)
        wt=wt0*wt*wtwp*wtepnn
      endif

      do nu=1,4
        p(3,nu)=nn(nu)
        p(4,nu)=ep(nu)
        p(5,nu)=bm(nu)
      enddo
      
                             
      return
      
  999 wt=0._dp
      return 1
      
      end
      
      
      
      
      
      
      
      
      
      
      
