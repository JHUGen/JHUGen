      subroutine coefsdk(s12,mtsq,ct,cv,c1)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      include 'alfacut.f'
      include 'includect.f'
      real(dp):: cv,ct,c1,s12,mtsq,ddilog,rsq,omrsq,eta,Kfun,
     & wlog,rlog,mulog

      if (scheme =='dred') then
C------        eta=0 4.e-_dphel
         eta=0._dp
      elseif (scheme == 'tH-V') then
C------       eta=1 t'Hooft Veltman
         eta=1._dp
      endif

      rsq=s12/mtsq
      omrsq=1._dp-rsq
      wlog=log(omrsq)
      rlog=log(rsq)

c--- epinv stands for (4*pi)^ep/Gamma(1-ep)/ep  (as usual)
c----  ct is the integrated counter-term, including alpha-dependence,
c----  see Eq. (10) of arXiv:1102.1967
      mulog=log(musq/mtsq)
      if (includect) then
      ct=(epinv2*epinv+epinv*mulog+0.5_dp*mulog**2)
     & +(epinv+mulog)*(2.5_dp-2._dp*wlog)
     & +25._dp/4._dp+0.5_dp*(1._dp/omrsq**2-8._dp/omrsq+7._dp)*rlog
     & +0.5_dp/omrsq+2._dp*ddilog(omrsq)-5._dp*pisqo6
     & -5._dp*wlog+2._dp*wlog**2+eta/2._dp
     & -2._dp*log(aff)**2-(3.5_dp-4._dp*aff+aff**2/2._dp)*log(aff)
     & +2._dp*(1._dp-aff)*rsq/omrsq*log(rsq/(1._dp-aff+rsq*aff))
      else
      ct=0._dp
      endif

C---- this routine has been constructed from 
C---- %\cite{Gottschalk:1980rv}
C---- \bibitem{Gottschalk:1980rv}
C---- T.~Gottschalk,
C---- %``Chromodynamic Corrections To Neutrino Production Of Heavy Quarks,''
C---- Phys.\ Rev.\ D {\bf 23}, 56 (1981).
C---- %%CITATION = PHRVA,D23,56;%%
C----- Adapted from Eqs.(A8,A9)

      Kfun=1._dp/rsq*wlog
      c1=2._dp*Kfun      
      cv=-(epinv2*epinv+epinv*mulog+0.5_dp*mulog**2)
     & -(epinv+mulog)*(2.5_dp-2._dp*wlog)
     & -0.5_dp*(11._dp+eta)-pisqo6-2._dp*ddilog(rsq)
     &  +3._dp*wlog-2._dp*wlog**2-Kfun
      
      return
      end

