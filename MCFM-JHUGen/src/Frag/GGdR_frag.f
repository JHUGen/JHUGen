      subroutine GGdR_frag(x,i,D,str) 
      implicit none
      include 'types.f'
c--- LO fragmentation functions from:
c---   A.~Gehrmann-De Ridder, E.~W.~N.~Glover,
c---   %``Final state photon production at LEP,''
c---   Eur.\ Phys.\ J.\  {\bf C7}, 29-48 (1999).  [hep-ph/9806316].
       
      include 'frag.f'
      include 'qcdcouple.f'
      include 'ewcouple.f' 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcharge.f'
      real(dp):: x,mu_zero,omx
      real(dp):: pref_ew,pref_ewas
      real(dp):: Pqga,Pqq,P1loop
      real(dp):: aewon2pi,log_mu,ddilog
      real(dp):: lx,dix,lmx,D,Dnp,PxP,PxD
      real(dp):: GGdR_NonP,GGdR_PxP,GGdR_DxP
      integer:: pow,str,i
      real(dp):: GdRG_NLO_frag
      external GdRG_NLO_frag 
      
c--- Safety cuts::
c---    x < 0.9999999 to avoid x->1 singularity
c---    x > 0.0001
      if((x > 0.9999999_dp) .or. (x < 0.0001_dp)) then 
         D=0._dp 
         return 
      endif

c--- set up logarithms 
      omx=one-x
      lx=log(x)
      dix=ddilog(1-x) 
      lmx=log(omx) 
      pow=2

c--- no gluon fragmentation
      if (i==0) then 
        D=0._dp 
        return 
      endif

      aewon2pi=esq/(fourpi*twopi) 
  
c--- prefactors 
      pref_ew=aewon2pi*Q(i)**2
      if     (str==0) then 
c----- LO only, set alpha_s terms= 0 
         pref_ewas = 0._dp 
         mu_zero=0.14_dp 
      elseif (str==1) then 
         D=GdRG_NLO_frag(x,i) 
         return 
!         write(6,*) 'NLO not currently implemented.'
!         stop 
      else
         write(6,*) 'Unrecognised order in Frag'  
      endif

      log_mu=log(frag_scale**2/mu_zero**2)

      
c-- splitting functions      
      Pqga=(one+omx**2)/x
c      P1loop=Cf*(-1._dp/2._dp+9._dp/2._dp*x+(-8._dp+x/2._dp)*lx+2._dp*x*lmx
c     &     +(1._dp-x/2._dp)*lx**2+(lmx**2+4._dp*lx*lmx+8._dp*dix
c     &     -4._dp*pi**2/3._dp)*Pqga)

    
!----- Non perturbative pieces 
      Dnp=GGdR_NonP(x,str)

!      Pxp=GGdR_PxP(x) 
!      PxD=GGdR_DxP(x)      

      D=pref_ew*Dnp+log_mu*pref_ew*Pqga

      return 
      end
       
      
      
      function GGdR_NonP(x,str) 
      implicit none
      include 'types.f'
      real(dp):: GGdR_NonP
!-----returns the non-pertubrbative pieces str=0,LO 1,NLO 
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      real(dp):: x,Pqga
      real(dp):: omx,lmx
      integer:: str

      omx=one-x
      lmx=log(omx**2)
      Pqga=(one+omx**2)/x

      if     (str==0) then 
c-- LO Non Pert eq. 2.13 with prefactor removed 
         GGdR_NonP=-Pqga*lmx-13.26_dp 
      elseif (str==1) then
         GGdR_NonP=-Pqga*lmx+20.8_dp*omx-11.07_dp 
      else
         write(6,*) 'Unrecognised Non P input in GGdR' 
         stop
      endif

      return 
      end


