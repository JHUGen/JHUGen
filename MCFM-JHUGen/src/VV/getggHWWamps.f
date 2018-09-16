      subroutine getggHWWamps(p,Mloop_bquark,Mloop_tquark)
      implicit none
      include 'types.f'
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process gg->Higgs->ZZ; there are:
c---        Mloop_bquark(h1,h2,h34,h56)   top quark mass=mt
c---        Mloop_tquark(h1,h2,h34,h56)   bottom quark mass=mb
c---
c--- The overall factor on the amplitude is:
c---
c---      (2._dp*gwsq*gsq/(16._dp*pisq)*gwsq/2._dp)**2 * delta(a,b)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f' 
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),msqgg,fac
      real(dp):: mfsq,tau,tauinv,rt,rescale
      complex(dp):: Mloop_tquark(2,2),Mloop_bquark(2,2),
     & fachiggs,amphiggs,f,e3De4,higgsprop,props

      call spinoru(6,p,za,zb)
      
      fachiggs=higgsprop(s(1,2))

c--- fill amplitudes with contributions of Higgs: top loop
      mfsq=mt**2
      tau=s(1,2)/(4._dp*mfsq)
      tauinv=1._dp/tau

      if (tau <= 1._dp) then
         f=cplx1(asin(sqrt(tau))**2)
      elseif (tau > 1._dp) then
         rt=sqrt(1._dp-tauinv)
         f=-0.25_dp*(cplx1(log((1._dp+rt)/(1._dp-rt)))-im*pi)**2
      else
         f=czip
      endif
      e3De4=2._dp*za(3,5)*zb(6,4)/(s(3,4)*s(5,6))
      amphiggs=mfsq*(cone+(cone-cplx1(tauinv))*f)*im*e3De4
      Mloop_tquark(1,1)=fachiggs*amphiggs*za(1,2)/zb(2,1)
      Mloop_tquark(1,2)=czip
      Mloop_tquark(2,1)=czip
      Mloop_tquark(2,2)=fachiggs*amphiggs*zb(1,2)/za(2,1)

c--- fill amplitudes with contributions of Higgs: bottom loop
      mfsq=mb**2
      tau=s(1,2)/(4._dp*mfsq)
      tauinv=1._dp/tau

      if (tau <= 1._dp) then
         f=cplx1(asin(sqrt(tau))**2)
      elseif (tau > 1._dp) then
         rt=sqrt(1._dp-tauinv)
         f=-0.25_dp*(cplx1(log((1._dp+rt)/(1._dp-rt)))-im*pi)**2
      else
         f=czip
      endif
      amphiggs=mfsq*(cone+(cone-cplx1(tauinv))*f)*im*e3De4

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
     &  s(3,4)/cplx2(s(3,4)-wmass**2,wmass*wwidth)
     & *s(5,6)/cplx2(s(5,6)-wmass**2,wmass*wwidth)
c      props=cone ! debug phase
      Mloop_bquark(:,:)=Mloop_bquark(:,:)*props
      Mloop_tquark(:,:)=Mloop_tquark(:,:)*props
     
      return
      end
      
