      subroutine getggHWWamps(p,Mloop_bquark,Mloop_tquark)
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process gg->Higgs->ZZ; there are:
c---        Mloop_bquark(h1,h2,h34,h56)   top quark mass=mt
c---        Mloop_tquark(h1,h2,h34,h56)   bottom quark mass=mb
c---
c--- The overall factor on the amplitude is:
c---
c---      (2d0*gwsq*gsq/(16d0*pisq)*gwsq/2d0)**2 * delta(a,b)
c---
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f' 
      include 'spinzerohiggs_anomcoupl.f'
      double precision p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac
      double precision mfsq,tau,tauinv,rt,rescale
      double complex Mloop_tquark(2,2),Mloop_bquark(2,2),
     & fachiggs,amphiggs,f,e3De4,higgsprop,props
      double complex anomhwwamp,qlI3,C0,a1,a3,a1_4gen,a3_4gen


      Mloop_tquark(:,:)=czip
      Mloop_bquark(:,:)=czip
      if(hmass.lt.zip) then
         return
      endif

      call spinoru(6,p,za,zb)
      
      fachiggs=higgsprop(s(1,2))

      e3De4=-2d0*anomhwwamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)
     & /(s(3,4)*s(5,6))      
      
      

c    Couplings for point-like interactions
      a1 = ghg2+ghg3*s(1,2)/4d0/LambdaBSM**2
      a3 = -2d0*ghg4/s(1,2)
      a1_4gen= ghg2_4gen+ghg3_4gen*s(1,2)/4d0/LambdaBSM**2
      a3_4gen= -2d0*ghg4_4gen/s(1,2)

      
      
      
! c--- fill amplitudes with contributions of Higgs: SM top loop
!       mfsq=mt**2
!       tau=s(1,2)/(4d0*mfsq)
!       tauinv=1d0/tau
! 
!       if (tau .le. 1d0) then
!          f=dcmplx(dasin(sqrt(tau))**2)
!       elseif (tau .gt. 1d0) then
!          rt=sqrt(1d0-tauinv)
!          f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
!       else
!          f=czip
!       endif
!       amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4     
!       Mloop_tquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
!       Mloop_tquark(1,2)=czip
!       Mloop_tquark(2,1)=czip
!       Mloop_tquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)

      
c        SM top
c        kappa couplings (quark loop interaction)      
      mfsq=mt**2
      C0=qlI3(zip,zip,s(1,2),mfsq,mfsq,mfsq,musq,0)
      amphiggs=mfsq/(-2d0)*im*e3De4*
     &   (-s(1,2)*C0*kappa_top*(1d0-4d0*mfsq/s(1,2))
     &   +2d0*kappa_top
     &   +s(1,2)*C0*kappa_tilde_top )*za(1,2)/zb(1,2)
      Mloop_tquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_tquark(1,2)=czip
      Mloop_tquark(2,1)=czip
      Mloop_tquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)     

c        ghg couplings (point-like interaction)
      Mloop_tquark(2,2)=Mloop_tquark(2,2) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1*zb(1,2)/za(1,2)+a3*0.5d0*zb(1,2)**2)     

      Mloop_tquark(1,1)=Mloop_tquark(1,1) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1*za(1,2)/zb(1,2)-a3*0.5d0*za(1,2)**2)

     
      
c        4th generation top
c        kappa couplings (quark loop interaction)
      mfsq=mt_4gen**2
      C0=qlI3(zip,zip,s(1,2),mfsq,mfsq,mfsq,musq,0)
      amphiggs=mfsq/(-2d0)*im*e3De4*
     &   (-s(1,2)*C0*kappa_4gen_top*(1d0-4d0*mfsq/s(1,2))
     &   +2d0*kappa_4gen_top
     &   +s(1,2)*C0*kappa_tilde_4gen_top )*za(1,2)/zb(1,2)
      Mloop_tquark(1,1)=Mloop_tquark(1,1) + 
     &                  fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_tquark(2,2)=Mloop_tquark(2,2) + 
     &                  fachiggs*amphiggs*zb(1,2)/za(2,1) 

c        ghg couplings (point-like interaction)
      Mloop_tquark(2,2)=Mloop_tquark(2,2) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)     

      Mloop_tquark(1,1)=Mloop_tquark(1,1) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1_4gen*za(1,2)/zb(1,2)-a3_4gen*0.5d0*za(1,2)**2)    
      
      
! c--- fill amplitudes with contributions of Higgs: bottom loop
!       mfsq=mb**2
!       tau=s(1,2)/(4d0*mfsq)
!       tauinv=1d0/tau
! 
!       if (tau .le. 1d0) then
!          f=dcmplx(dasin(sqrt(tau))**2)
!       elseif (tau .gt. 1d0) then
!          rt=sqrt(1d0-tauinv)
!          f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
!       else
!          f=czip
!       endif
!       amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4
! 
!       Mloop_bquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
!       Mloop_bquark(1,2)=czip
!       Mloop_bquark(2,1)=czip
!       Mloop_bquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)


      
c        SM bot
c        kappa couplings (quark loop interaction)      
      mfsq=mb**2
      C0=qlI3(zip,zip,s(1,2),mfsq,mfsq,mfsq,musq,0)
      amphiggs=mfsq/(-2d0)*im*e3De4*
     &   (-s(1,2)*C0*kappa_bot*(1d0-4d0*mfsq/s(1,2))
     &   +2d0*kappa_bot
     &   +s(1,2)*C0*kappa_tilde_bot )*za(1,2)/zb(1,2)
      Mloop_bquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_bquark(1,2)=czip
      Mloop_bquark(2,1)=czip
      Mloop_bquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)     

c        ghg couplings (point-like interaction)
      Mloop_bquark(2,2)=Mloop_tquark(2,2) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1*zb(1,2)/za(1,2)+a3*0.5d0*zb(1,2)**2)     

      Mloop_bquark(1,1)=Mloop_tquark(1,1) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1*za(1,2)/zb(1,2)-a3*0.5d0*za(1,2)**2)      


c        4th generation bot
c        kappa couplings (quark loop interaction)
      mfsq=mb_4gen**2
      C0=qlI3(zip,zip,s(1,2),mfsq,mfsq,mfsq,musq,0)
      amphiggs=mfsq/(-2d0)*im*e3De4*
     &   (-s(1,2)*C0*kappa_4gen_bot*(1d0-4d0*mfsq/s(1,2))
     &   +2d0*kappa_4gen_bot
     &   +s(1,2)*C0*kappa_tilde_4gen_bot )*za(1,2)/zb(1,2)
      Mloop_bquark(1,1)=Mloop_bquark(1,1) + 
     &                  fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_bquark(2,2)=Mloop_bquark(2,2) + 
     &                  fachiggs*amphiggs*zb(1,2)/za(2,1) 

c        ghg couplings (point-like interaction)
      Mloop_bquark(2,2)=Mloop_bquark(2,2) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)     

      Mloop_bquark(1,1)=Mloop_bquark(1,1) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1_4gen*za(1,2)/zb(1,2)-a3_4gen*0.5d0*za(1,2)**2)          
      
      
      
      
      
c--- Rescale for width study
      if((keep_smhiggs_norm).and.(anom_higgs)) then 
         rescale=chi_higgs**2 
         Mloop_tquark(:,:)=Mloop_tquark(:,:)*rescale
         Mloop_bquark(:,:)=Mloop_bquark(:,:)*rescale
      endif
      
      props=
     &  s(3,4)/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
     & *s(5,6)/dcmplx(s(5,6)-wmass**2,wmass*wwidth)
c      props=cone ! debug phase
      Mloop_bquark(:,:)=Mloop_bquark(:,:)*props
      Mloop_tquark(:,:)=Mloop_tquark(:,:)*props
     
      return
      end
      

      
      
      
      
      
      
      
      
      
      
      
      subroutine getggH2WWamps(p,Mloop_bquark,Mloop_tquark)
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process gg->Higgs->ZZ; there are:
c---        Mloop_bquark(h1,h2,h34,h56)   top quark mass=mt
c---        Mloop_tquark(h1,h2,h34,h56)   bottom quark mass=mb
c---
c--- The overall factor on the amplitude is:
c---
c---      (2d0*gwsq*gsq/(16d0*pisq)*gwsq/2d0)**2 * delta(a,b)
c---
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f' 
      include 'spinzerohiggs_anomcoupl.f'
      double precision p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac
      double precision mfsq,tau,tauinv,rt,rescale
      double complex Mloop_tquark(2,2),Mloop_bquark(2,2),
     & fachiggs,amphiggs,f,e3De4,higgs2prop,props
      double complex anomhwwamp,qlI3,C0,a1,a3,a1_4gen,a3_4gen
      

      Mloop_tquark(:,:)=czip
      Mloop_bquark(:,:)=czip
      if(h2mass.lt.zip) then
         return
      endif

      call spinoru(6,p,za,zb)
      
      fachiggs=higgs2prop(s(1,2))

      e3De4=-2d0*anomhwwamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)
     & /(s(3,4)*s(5,6))  

c    Couplings for point-like interactions
      a1 = gh2g2+gh2g3*s(1,2)/4d0/Lambda2BSM**2
      a3 = -2d0*gh2g4/s(1,2)
      a1_4gen= gh2g2_4gen+gh2g3_4gen*s(1,2)/4d0/Lambda2BSM**2
      a3_4gen= -2d0*gh2g4_4gen/s(1,2)


      
      
c--- fill amplitudes with contributions of Higgs: top loop
!       mfsq=mt**2
!       tau=s(1,2)/(4d0*mfsq)
!       tauinv=1d0/tau
!       if (tau .le. 1d0) then
!          f=dcmplx(dasin(sqrt(tau))**2)
!       elseif (tau .gt. 1d0) then
!          rt=sqrt(1d0-tauinv)
!          f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
!       else
!          f=czip
!       endif
!       amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4
!       Mloop_tquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
!       Mloop_tquark(1,2)=czip
!       Mloop_tquark(2,1)=czip
!       Mloop_tquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)

      
c        SM top
c        kappa2 couplings (quark loop interaction)      
      mfsq=mt**2
      C0=qlI3(zip,zip,s(1,2),mfsq,mfsq,mfsq,musq,0)
      amphiggs=mfsq/(-2d0)*im*e3De4*
     &   (-s(1,2)*C0*kappa2_top*(1d0-4d0*mfsq/s(1,2))
     &   +2d0*kappa2_top
     &   +s(1,2)*C0*kappa2_tilde_top )*za(1,2)/zb(1,2)
      Mloop_tquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_tquark(1,2)=czip
      Mloop_tquark(2,1)=czip
      Mloop_tquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)     

c        ghg couplings (point-like interaction)
      Mloop_tquark(2,2)=Mloop_tquark(2,2) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1*zb(1,2)/za(1,2)+a3*0.5d0*zb(1,2)**2)     

      Mloop_tquark(1,1)=Mloop_tquark(1,1) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1*za(1,2)/zb(1,2)-a3*0.5d0*za(1,2)**2)

     
      
c        4th generation top
c        kappa2 couplings (quark loop interaction)
      mfsq=mt_4gen**2
      C0=qlI3(zip,zip,s(1,2),mfsq,mfsq,mfsq,musq,0)
      amphiggs=mfsq/(-2d0)*im*e3De4*
     &   (-s(1,2)*C0*kappa2_4gen_top*(1d0-4d0*mfsq/s(1,2))
     &   +2d0*kappa2_4gen_top
     &   +s(1,2)*C0*kappa2_tilde_4gen_top )*za(1,2)/zb(1,2)
      Mloop_tquark(1,1)=Mloop_tquark(1,1) + 
     &                  fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_tquark(2,2)=Mloop_tquark(2,2) + 
     &                  fachiggs*amphiggs*zb(1,2)/za(2,1) 

c        ghg couplings (point-like interaction)
      Mloop_tquark(2,2)=Mloop_tquark(2,2) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)     

      Mloop_tquark(1,1)=Mloop_tquark(1,1) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1_4gen*za(1,2)/zb(1,2)-a3_4gen*0.5d0*za(1,2)**2)    
            
      
      
      
c--- fill amplitudes with contributions of Higgs: bottom loop
!       mfsq=mb**2
!       tau=s(1,2)/(4d0*mfsq)
!       tauinv=1d0/tau
!       if (tau .le. 1d0) then
!          f=dcmplx(dasin(sqrt(tau))**2)
!       elseif (tau .gt. 1d0) then
!          rt=sqrt(1d0-tauinv)
!          f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
!       else
!          f=czip
!       endif
!       amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4
!       Mloop_bquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
!       Mloop_bquark(1,2)=czip
!       Mloop_bquark(2,1)=czip
!       Mloop_bquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)


c        SM bot
c        kappa2 couplings (quark loop interaction)      
      mfsq=mb**2
      C0=qlI3(zip,zip,s(1,2),mfsq,mfsq,mfsq,musq,0)
      amphiggs=mfsq/(-2d0)*im*e3De4*
     &   (-s(1,2)*C0*kappa2_bot*(1d0-4d0*mfsq/s(1,2))
     &   +2d0*kappa2_bot
     &   +s(1,2)*C0*kappa2_tilde_bot )*za(1,2)/zb(1,2)
      Mloop_bquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_bquark(1,2)=czip
      Mloop_bquark(2,1)=czip
      Mloop_bquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)     

c        ghg couplings (point-like interaction)
      Mloop_bquark(2,2)=Mloop_tquark(2,2) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1*zb(1,2)/za(1,2)+a3*0.5d0*zb(1,2)**2)     

      Mloop_bquark(1,1)=Mloop_tquark(1,1) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1*za(1,2)/zb(1,2)-a3*0.5d0*za(1,2)**2)      


c        4th generation bot
c        kappa2 couplings (quark loop interaction)
      mfsq=mb_4gen**2
      C0=qlI3(zip,zip,s(1,2),mfsq,mfsq,mfsq,musq,0)
      amphiggs=mfsq/(-2d0)*im*e3De4*
     &   (-s(1,2)*C0*kappa2_4gen_bot*(1d0-4d0*mfsq/s(1,2))
     &   +2d0*kappa2_4gen_bot
     &   +s(1,2)*C0*kappa2_tilde_4gen_bot )*za(1,2)/zb(1,2)
      Mloop_bquark(1,1)=Mloop_bquark(1,1) + 
     &                  fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_bquark(2,2)=Mloop_bquark(2,2) + 
     &                  fachiggs*amphiggs*zb(1,2)/za(2,1) 

c        ghg couplings (point-like interaction)
      Mloop_bquark(2,2)=Mloop_bquark(2,2) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)     

      Mloop_bquark(1,1)=Mloop_bquark(1,1) + 
     &   fachiggs*s(1,2)/3d0/(-2d0)*im*e3De4*
     &   (a1_4gen*za(1,2)/zb(1,2)-a3_4gen*0.5d0*za(1,2)**2)          
      
      
            
      
      
      
      
      
c--- Rescale for width study
      if((keep_smhiggs_norm).and.(anom_higgs)) then 
         rescale=chi_higgs**2 
         Mloop_tquark(:,:)=Mloop_tquark(:,:)*rescale
         Mloop_bquark(:,:)=Mloop_bquark(:,:)*rescale
      endif
      
      props=
     &  s(3,4)/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
     & *s(5,6)/dcmplx(s(5,6)-wmass**2,wmass*wwidth)
c      props=cone ! debug phase
      Mloop_bquark(:,:)=Mloop_bquark(:,:)*props
      Mloop_tquark(:,:)=Mloop_tquark(:,:)*props
     
      return
      end
