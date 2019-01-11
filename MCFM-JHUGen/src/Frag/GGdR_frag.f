      subroutine GGdR_frag(x,i,D,str) 
c--- LO fragmentation functions from:
c---   A.~Gehrmann-De Ridder, E.~W.~N.~Glover,
c---   %``Final state photon production at LEP,''
c---   Eur.\ Phys.\ J.\  {\bf C7}, 29-48 (1999).  [hep-ph/9806316].
      implicit none 
      include 'frag.f'
      include 'qcdcouple.f'
      include 'ewcouple.f' 
      include 'constants.f'
      include 'ewcharge.f'
      double precision x,mu_zero,omx
      double precision pref_ew,pref_ewas
      double precision Pqga,Pqq,P1loop
      double precision aewon2pi,log_mu,ddilog
      double precision lx,dix,lmx,D,Dnp,PxP,PxD
      double precision GGdR_NonP,GGdR_PxP,GGdR_DxP
      integer pow,str,i
      double precision GdRG_NLO_frag
      external GdRG_NLO_frag 
      
c--- Safety cuts::
c---    x < 0.9999999 to avoid x->1 singularity
c---    x > 0.0001
      if((x .gt. 0.9999999d0) .or. (x .lt. 0.0001d0)) then 
         D=0d0 
         return 
      endif

c--- set up logarithms 
      omx=one-x
      lx=dlog(x)
      dix=ddilog(1-x) 
      lmx=dlog(omx) 
      pow=2

c--- no gluon fragmentation
      if (i.eq.0) then 
        D=0d0 
        return 
      endif

      aewon2pi=esq/(fourpi*twopi) 
  
c--- prefactors 
      pref_ew=aewon2pi*Q(i)**2
      if     (str.eq.0) then 
c----- LO only, set alpha_s terms= 0 
         pref_ewas = 0d0 
         mu_zero=0.14d0 
      elseif (str.eq.1) then 
         D=GdRG_NLO_frag(x,i) 
         return 
!         write(6,*) 'NLO not currently implemented.'
!         stop 
      else
         write(6,*) 'Unrecognised order in Frag'  
      endif

      log_mu=dlog(frag_scale**2/mu_zero**2)

      
c-- splitting functions      
      Pqga=(one+omx**2)/x
c      P1loop=Cf*(-1d0/2d0+9d0/2d0*x+(-8d0+x/2d0)*lx+2d0*x*lmx
c     &     +(1d0-x/2d0)*lx**2+(lmx**2+4d0*lx*lmx+8d0*dix
c     &     -4d0*pi**2/3d0)*Pqga)

    
!----- Non perturbative pieces 
      Dnp=GGdR_NonP(x,str)

!      Pxp=GGdR_PxP(x) 
!      PxD=GGdR_DxP(x)      

      D=pref_ew*Dnp+log_mu*pref_ew*Pqga

      return 
      end
       
      
      
      double precision function GGdR_NonP(x,str) 
!-----returns the non-pertubrbative pieces str=0,LO 1,NLO 
      implicit none 
      include 'constants.f' 
      double precision x,Pqga
      double precision omx,lmx
      integer str

      omx=one-x
      lmx=dlog(omx**2)
      Pqga=(one+omx**2)/x

      if     (str.eq.0) then 
c-- LO Non Pert eq. 2.13 with prefactor removed 
         GGdR_NonP=-Pqga*lmx-13.26d0 
      elseif (str.eq.1) then
         GGdR_NonP=-Pqga*lmx+20.8d0*omx-11.07d0 
      else
         write(6,*) 'Unrecognised Non P input in GGdR' 
         stop
      endif

      return 
      end


