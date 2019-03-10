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
      double complex FFa1, FFa2, FFa3
      double complex aa1,aa2,aa3
      double precision shat,q3_q3,q4_q4
      double complex ghW1_dyn,ghw2_dyn,ghw3_dyn,ghw4_dyn

      Mloop_tquark(:,:)=czip
      Mloop_bquark(:,:)=czip
      if(hmass.lt.zip) then
         return
      endif

      call spinoru(6,p,za,zb)
      
      fachiggs=higgsprop(s(1,2))

c--- fill amplitudes with contributions of Higgs: top loop
      mfsq=mt**2
      tau=s(1,2)/(4d0*mfsq)
      tauinv=1d0/tau

      if (tau .le. 1d0) then
         f=dcmplx(dasin(sqrt(tau))**2)
      elseif (tau .gt. 1d0) then
         rt=sqrt(1d0-tauinv)
         f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
      else
         f=czip
      endif

      IF( AllowAnomalousCouplings ) THEN

      shat  = s(1,2)
      q3_q3 = s(3,4)
      q4_q4 = s(5,6)

c--- q^2-dependent couplings
      call HVVSpinZeroDynCoupl(ghw1_dyn,1,1,q3_q3,q4_q4,shat,.true.)
      call HVVSpinZeroDynCoupl(ghw2_dyn,2,1,q3_q3,q4_q4,shat,.true.)
      call HVVSpinZeroDynCoupl(ghw3_dyn,3,1,q3_q3,q4_q4,shat,.true.)
      call HVVSpinZeroDynCoupl(ghw4_dyn,4,1,q3_q3,q4_q4,shat,.true.)

      aa1 =ghw1_dyn*wmass**2/shat
     &     + (s(1,2)-s(3,4)-s(5,6))/shat*
     &       (ghw2_dyn
     &       +ghw3_dyn*(s(1,2)-s(3,4)-s(5,6))/4d0/LambdaBSM**2)
      aa2 =-2d0*ghw2_dyn
     &     -ghw3_dyn*(s(1,2)-s(3,4)-s(5,6))/2d0/LambdaBSM**2
      aa3 =-2d0*ghw4_dyn

      aa1 = aa1 / wmass**2 
      aa2 = aa2 / wmass**2 
      aa3 = aa3 / wmass**2

c------ FORM FACTORS FOR ANOMALOUS COUPLINGS
c L1L2
      FFa1 = za(3,5)*zb(4,6)*shat
      FFa2 = -0.5d0*za(3,5)**2*zb(3,6)*zb(4,5) 
     &       -0.5d0*za(3,5)*za(3,6)*zb(3,6)*zb(4,6) 
     &       -0.5d0*za(3,5)*za(4,5)*zb(4,5)*zb(4,6) 
     &       -0.5d0*za(3,6)*za(4,5)*zb(4,6)**2
      FFa3 =  0.5d0*za(3,4)*za(5,6)*zb(4,6)**2 
     &       -0.5d0*za(3,5)**2*zb(3,4)*zb(5,6)

      FFa3 = FFa3 * (0d0,-1d0)!  phase convention to match JHUGen

      e3De4=-2d0*( aa1*FFa1 + aa2*FFa2 + aa3*FFa3 )/(s(3,4)*s(5,6))
      
      ELSE

      e3De4=2d0*za(3,5)*zb(6,4)/(s(3,4)*s(5,6))

      ENDIF

      amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4
      Mloop_tquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_tquark(1,2)=czip
      Mloop_tquark(2,1)=czip
      Mloop_tquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)

c--- fill amplitudes with contributions of Higgs: bottom loop
      mfsq=mb**2
      tau=s(1,2)/(4d0*mfsq)
      tauinv=1d0/tau

      if (tau .le. 1d0) then
         f=dcmplx(dasin(sqrt(tau))**2)
      elseif (tau .gt. 1d0) then
         rt=sqrt(1d0-tauinv)
         f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
      else
         f=czip
      endif
      amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4

      Mloop_bquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_bquark(1,2)=czip
      Mloop_bquark(2,1)=czip
      Mloop_bquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)

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
      double complex FFa1, FFa2, FFa3
      double complex aa1,aa2,aa3
      double precision shat,q3_q3,q4_q4
      double complex ghW1_dyn,ghw2_dyn,ghw3_dyn,ghw4_dyn

      Mloop_tquark(:,:)=czip
      Mloop_bquark(:,:)=czip
      if(h2mass.lt.zip) then
         return
      endif

      call spinoru(6,p,za,zb)
      
      fachiggs=higgs2prop(s(1,2))

c--- fill amplitudes with contributions of Higgs: top loop
      mfsq=mt**2
      tau=s(1,2)/(4d0*mfsq)
      tauinv=1d0/tau

      if (tau .le. 1d0) then
         f=dcmplx(dasin(sqrt(tau))**2)
      elseif (tau .gt. 1d0) then
         rt=sqrt(1d0-tauinv)
         f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
      else
         f=czip
      endif

      IF( AllowAnomalousCouplings ) THEN

      shat  = s(1,2)
      q3_q3 = s(3,4)
      q4_q4 = s(5,6)

c--- q^2-dependent couplings
      call HVVSpinZeroDynCoupl(ghw1_dyn,1,2,q3_q3,q4_q4,shat,.true.)
      call HVVSpinZeroDynCoupl(ghw2_dyn,2,2,q3_q3,q4_q4,shat,.true.)
      call HVVSpinZeroDynCoupl(ghw3_dyn,3,2,q3_q3,q4_q4,shat,.true.)
      call HVVSpinZeroDynCoupl(ghw4_dyn,4,2,q3_q3,q4_q4,shat,.true.)

      aa1 =ghw1_dyn*wmass**2/shat
     &     + (s(1,2)-s(3,4)-s(5,6))/shat*
     &       (ghw2_dyn
     &       +ghw3_dyn*(s(1,2)-s(3,4)-s(5,6))/4d0/Lambda2BSM**2)
      aa2 =-2d0*ghw2_dyn
     &     -ghw3_dyn*(s(1,2)-s(3,4)-s(5,6))/2d0/Lambda2BSM**2
      aa3 =-2d0*ghw4_dyn

      aa1 = aa1 / wmass**2 
      aa2 = aa2 / wmass**2 
      aa3 = aa3 / wmass**2

c------ FORM FACTORS FOR ANOMALOUS COUPLINGS
c L1L2
      FFa1 = za(3,5)*zb(4,6)*shat
      FFa2 = -0.5d0*za(3,5)**2*zb(3,6)*zb(4,5) 
     &       -0.5d0*za(3,5)*za(3,6)*zb(3,6)*zb(4,6) 
     &       -0.5d0*za(3,5)*za(4,5)*zb(4,5)*zb(4,6) 
     &       -0.5d0*za(3,6)*za(4,5)*zb(4,6)**2
      FFa3 =  0.5d0*za(3,4)*za(5,6)*zb(4,6)**2 
     &       -0.5d0*za(3,5)**2*zb(3,4)*zb(5,6)

      FFa3 = FFa3 * (0d0,-1d0)!  phase convention to match JHUGen

      e3De4=-2d0*( aa1*FFa1 + aa2*FFa2 + aa3*FFa3 )/(s(3,4)*s(5,6))
      
      ELSE

      e3De4=2d0*za(3,5)*zb(6,4)/(s(3,4)*s(5,6))

      ENDIF

      amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4
      Mloop_tquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_tquark(1,2)=czip
      Mloop_tquark(2,1)=czip
      Mloop_tquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)

c--- fill amplitudes with contributions of Higgs: bottom loop
      mfsq=mb**2
      tau=s(1,2)/(4d0*mfsq)
      tauinv=1d0/tau

      if (tau .le. 1d0) then
         f=dcmplx(dasin(sqrt(tau))**2)
      elseif (tau .gt. 1d0) then
         rt=sqrt(1d0-tauinv)
         f=-0.25d0*(dcmplx(log((1d0+rt)/(1d0-rt)))-im*pi)**2
      else
         f=czip
      endif
      amphiggs=mfsq*(cone+(cone-dcmplx(tauinv))*f)*im*e3De4

      Mloop_bquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_bquark(1,2)=czip
      Mloop_bquark(2,1)=czip
      Mloop_bquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)

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
