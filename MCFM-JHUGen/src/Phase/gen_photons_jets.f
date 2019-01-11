      subroutine gen_photons_jets(r,nphots,njets,p,wt,*)
c----  WARNING: although mostly written in a generic manner, this
c----           routine has been tailored for direct photon (dirgam)
c----           and photon+two jets (gamjet) production

c---- generate phase space for 2-->nphots+njets process
c----   and 3,..,2+nphots the photons
c----   and 2+nphots,..,2+nphots+njets the jets
c----
c---- This routine is adapted from gen_Vphotons_jets.f
c----
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
c----
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'leptcuts.f'
      include 'reset.f'
      include 'xmin.f'
      include 'part.f'
      include 'process.f'
      include 'x1x2.f'
      include 'first.f'
      include 'energy.f'
      double precision r(mxdim),p(mxpart,4),psumjet(4),wt
      double precision hmin,hmax,delh,h,xpt,etamax
      double precision y,sinhy,coshy,phi,xjac,tmp(mxpart,4)
      double precision ptjetmin,etajetmin,etajetmax,pbreak
      double precision ptmin_part,etamax_part,pbreak_part,swap,lastmax
      integer j,nu,nphots,njets,nphotsjets,ijet,icount
      logical flatreal
      integer,save:: perm=0
      integer,parameter:: perm2_3(0:1)=(/3,4/)
      integer,parameter:: perm2_4(0:1)=(/4,3/)
      integer,parameter:: perm3_3(0:5)=(/3,3,4,4,5,5/)
      integer,parameter:: perm3_4(0:5)=(/4,5,3,5,3,4/)
      integer,parameter:: perm3_5(0:5)=(/5,4,5,3,4,3/)
      integer,parameter:: perm4_3(0:23)=(/3,3,3,3,3,3,4,4,4,4,4,4,
     & 5,5,5,5,5,5,6,6,6,6,6,6/)
      integer,parameter:: perm4_4(0:23)=(/4,4,5,5,6,6,3,3,5,5,6,6,
     & 3,3,4,4,6,6,3,3,4,4,5,5/)
      integer,parameter:: perm4_5(0:23)=(/5,6,4,6,4,5,5,6,3,6,3,5,
     & 4,6,3,6,3,4,4,5,3,5,3,4/)
      integer,parameter:: perm4_6(0:23)=(/6,5,6,4,5,4,6,5,6,3,5,3,
     & 6,4,6,3,4,3,5,4,5,3,4,3/)
      save ptjetmin,etajetmin,etajetmax,pbreak
!$omp threadprivate(ptjetmin,etajetmin,etajetmax,pbreak,perm)
      
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


c--- total number of photons and jets
      nphotsjets=nphots+njets

      do nu=1,4
        do j=1,4+nphotsjets
          p(j,nu)=0d0
        enddo
        psumjet(nu)=0d0
      enddo 

      wt=2d0*pi
                       
     
      icount=0
c--- loop over all the jets and photons except the last one
      do ijet=1,nphotsjets-1
      
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
      wt=wt/16d0/pi**3
c      xmin=2d0/sqrts
c      xmax=1d0/ptjetmin

      if (ijet .le. nphots) then
        if     ((nphots .eq. 3) .or. (nphots .eq. 4)) then
          ptmin_part=min(gammpt,gammpt2,gammpt3)
        elseif (nphots .eq. 2) then
          ptmin_part=min(gammpt,gammpt2)
        elseif (nphots .eq. 1) then
          ptmin_part=gammpt
        else
          write(6,*) 'Unexpected # photons in gen_photons_jets: ',nphots
          stop
        endif
        etamax_part=gammrap
c        pbreak_part=0d0
        if (part .eq. 'real') then
c--- cannot generate exactly to match, since dipoles transform photon
          etamax_part=20d0
c          ptmin_part=0d0
c          pbreak_part=gammpt
        endif
      else
        if (part .eq. 'lord') then
c--- generate jets according to phase-space boundaries
          ptmin_part=ptjetmin
          etamax_part=etajetmax
c          pbreak_part=pbreak
        else
c--- cannot generate exactly to match
          etamax_part=20d0
c          ptmin_part=0d0
c--- attempt to generate jets to balance photon for dirgam,
          if (case .eq. 'dirgam') then
            ptmin_part=gammpt
c          else
c            pbreak_part=pbreak
          endif
        endif
      endif

c      if ((flatreal) .and. (part .eq. 'real')) then
cc--- generate flat pt for the real contribution
c         pt=r(icount+1)*((sqrts/2d0)-ptmin_part)+ptmin_part
c         wt=wt*((sqrts/2d0)-ptmin_part)*pt
c      else
cc--- favour small pt region 
c         hmin=1d0/dsqrt((sqrts/2d0)**2+pbreak_part**2)
c         hmax=1d0/dsqrt(ptmin_part**2+pbreak_part**2)
c         delh=hmax-hmin
c         h=hmin+r(icount+1)*delh        
c         pt=dsqrt(1d0/h**2-pbreak_part**2)
c         wt=wt*delh/h**3
c      endif
      
      if (part .eq. 'real') then
        call genpt(r(icount+1),ptmin_part,.false.,xpt,xjac)
      else
        call genpt(r(icount+1),ptmin_part,.true.,xpt,xjac)
      endif
      wt=wt*xjac
 
      etamax=sqrts/2d0/xpt
      if (etamax**2 .le. 1d0) then
        write(6,*) 'etamax**2 .le. 1d0 in gen_photons_jets.f',etamax**2
        wt=0d0
        return 1
      endif
      etamax=dlog(etamax+dsqrt(etamax**2-1d0))
      
      etamax=min(etamax,etamax_part)
      y=etamax*(2d0*r(icount+2)-1d0)
      wt=wt*2d0*etamax
      
      sinhy=dsinh(y)
      coshy=dsqrt(1d0+sinhy**2)
      
      p(2+ijet,4)=xpt*coshy
      
      phi=2d0*pi*r(icount+3)
      wt=wt*2d0*pi
      
      p(2+ijet,1)=xpt*dcos(phi)
      p(2+ijet,2)=xpt*dsin(phi)
      p(2+ijet,3)=xpt*sinhy
      
      do nu=1,4
        psumjet(nu)=psumjet(nu)+p(2+ijet,nu)
      enddo
     
      icount=icount+3
      
c--- end of loop over photons and jets      
      enddo

c---- have generated nphotsjets-1 massless momenta, last pt is fixed but rapidity is unconstrained
      xpt=dsqrt(psumjet(1)**2+psumjet(2)**2)
      etamax=sqrts/2d0/xpt
      if (etamax**2 .le. 1d0) then
c         write(6,*) 'etamax**2 .le. 1d0 in gen_photons_jets.f',etamax**2 
         wt=0d0
         return 1
      endif
      etamax=dlog(etamax+dsqrt(etamax**2-1d0))
        
      if ((part .eq. 'lord') .or. (part .eq. 'virt')) then
c--- rapidity constrained: need to check identity of last parton
        if (njets .eq. 0) then
          lastmax=gammrap
        else
          lastmax=etajetmax
        endif
      else
c--- rapidity not constrained      
        lastmax=20d0
      endif
           
      etamax=min(etamax,lastmax)
      
      y=etamax*(2d0*r(3*nphotsjets-2)-1d0)
      wt=wt*2d0*etamax
        
      sinhy=dsinh(y)
      coshy=dsqrt(1d0+sinhy**2)

      p(2+nphotsjets,4)=xpt*coshy
               
      p(2+nphotsjets,1)=-psumjet(1)
      p(2+nphotsjets,2)=-psumjet(2)
      p(2+nphotsjets,3)=xpt*sinhy
       
      do nu=1,4
      psumjet(nu)=psumjet(nu)+p(2+nphotsjets,nu)
      enddo

c--- now make the initial state momenta
      xx(1)=(psumjet(4)+psumjet(3))/sqrts
      xx(2)=(psumjet(4)-psumjet(3))/sqrts
      if (abs(psumjet(4)/psumjet(3)-1d0) .lt. 1d-12) then
        wt=0d0
        return 1
      endif
      
      if   ((xx(1) .gt. 1d0) .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin).or. (xx(2) .lt. xmin)) then
         wt=0d0
         return 1
      endif 
      
      p(1,4)=-0.5d0*xx(1)*sqrts
      p(1,1)=0d0
      p(1,2)=0d0
      p(1,3)=-0.5d0*xx(1)*sqrts

      p(2,4)=-0.5d0*xx(2)*sqrts
      p(2,1)=0d0
      p(2,2)=0d0
      p(2,3)=+0.5d0*xx(2)*sqrts

      wt=wt/sqrts**2
      
      do j=3+nphotsjets,mxpart
      do nu=1,4
      p(j,nu)=0d0
      enddo
      enddo
      
c--- symmetrize phase-space by permuting photon momenta appropriately
      if     (nphots .eq. 2) then
        perm=int(r(3*nphotsjets-1)*2d0)
        if (perm .lt. 0) perm=0
        if (perm .gt. 1) perm=1
        tmp(3,:)=p(3,:)  
        tmp(4,:)=p(4,:)  
        p(3,:)=tmp(perm2_3(perm),:)
        p(4,:)=tmp(perm2_4(perm),:)
c        perm=mod(perm+1,2)
      elseif (nphots .eq. 3) then
        perm=int(r(3*nphotsjets-1)*6d0)
        if (perm .lt. 0) perm=0
        if (perm .gt. 5) perm=5
        tmp(3,:)=p(3,:)  
        tmp(4,:)=p(4,:)  
        tmp(5,:)=p(5,:)  
        p(3,:)=tmp(perm3_3(perm),:)
        p(4,:)=tmp(perm3_4(perm),:)
        p(5,:)=tmp(perm3_5(perm),:)
c        perm=mod(perm+1,6)
      elseif (nphots.eq.4) then 
        perm=int(r(3*nphotsjets-1)*24d0)
        if (perm .lt. 0) perm=0
        if (perm .gt. 23) perm=23
        tmp(3,:)=p(3,:)  
        tmp(4,:)=p(4,:)  
        tmp(5,:)=p(5,:)  
        tmp(6,:)=p(6,:)  
        p(3,:)=tmp(perm4_3(perm),:)
        p(4,:)=tmp(perm4_4(perm),:)
        p(5,:)=tmp(perm4_5(perm),:)
        p(6,:)=tmp(perm4_6(perm),:)
c        perm=mod(perm+1,24)
      elseif (nphots .ge. 5) then
        write(6,*) 'Symmetrize gen_photons_jets.f before proceeding'
        write(6,*) 'nphots = ',nphots
        stop      
      endif
      
      return
      end
      
      
      
      
      
      
      
      
      
      
      
