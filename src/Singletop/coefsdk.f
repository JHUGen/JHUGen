      subroutine coefsdk(s12,mtsq,ct,cv,c1)
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      include 'alfacut.f'
      include 'includect.f'
      double precision cv,ct,c1,s12,mtsq,ddilog,rsq,omrsq,eta,Kfun,
     & wlog,rlog,mulog

      if (scheme .eq.'dred') then
C------        eta=0 4d-hel
         eta=0d0
      elseif (scheme .eq. 'tH-V') then
C------       eta=1 t'Hooft Veltman
         eta=1d0
      endif

      rsq=s12/mtsq
      omrsq=1d0-rsq
      wlog=dlog(omrsq)
      rlog=dlog(rsq)

c--- epinv stands for (4*pi)^ep/Gamma(1-ep)/ep  (as usual)
c----  ct is the integrated counter-term, including alpha-dependence,
c----  see Eq. (10) of arXiv:1102.1967
      mulog=dlog(musq/mtsq)
      if (includect) then
      ct=(epinv2*epinv+epinv*mulog+0.5d0*mulog**2)
     & +(epinv+mulog)*(2.5d0-2d0*wlog)
     & +25d0/4d0+0.5d0*(1d0/omrsq**2-8d0/omrsq+7d0)*rlog
     & +0.5d0/omrsq+2d0*ddilog(omrsq)-5d0*pisqo6
     & -5d0*wlog+2d0*wlog**2+eta/2d0
     & -2d0*dlog(aff)**2-(3.5d0-4d0*aff+aff**2/2d0)*dlog(aff)
     & +2d0*(1d0-aff)*rsq/omrsq*dlog(rsq/(1d0-aff+rsq*aff))
      else
      ct=0d0
      endif

C---- this routine has been constructed from 
C---- %\cite{Gottschalk:1980rv}
C---- \bibitem{Gottschalk:1980rv}
C---- T.~Gottschalk,
C---- %``Chromodynamic Corrections To Neutrino Production Of Heavy Quarks,''
C---- Phys.\ Rev.\ D {\bf 23}, 56 (1981).
C---- %%CITATION = PHRVA,D23,56;%%
C----- Adapted from Eqs.(A8,A9)

      Kfun=1d0/rsq*wlog
      c1=2d0*Kfun      
      cv=-(epinv2*epinv+epinv*mulog+0.5d0*mulog**2)
     & -(epinv+mulog)*(2.5d0-2d0*wlog)
     & -0.5d0*(11d0+eta)-pisqo6-2d0*ddilog(rsq)
     &  +3d0*wlog-2d0*wlog**2-Kfun
      
      return
      end

