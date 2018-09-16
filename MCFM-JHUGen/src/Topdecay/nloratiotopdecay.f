      function nloratiotopdecay(mt,mb,mw,gamw)
      implicit none
      include 'types.f'
      real(dp):: nloratiotopdecay
************************************************************************
*     Authors: R.K. Ellis and J. Campbell, February 2012               *
*                                                                      *
*     ratio NLO/LO for the width of the top quark, including the       *
*     effect of the bottom quark mass.                                 *
*                                                                      *
*     Formula is taken from Eq. (27) of                                *
*                                                                      *
*     \bibitem{Czarnecki:1990kv}                                       *
*     A.~Czarnecki,                                                    *
*     ``QCD corrections to the decay t ---> W b                        *
*       in dimensional regularization,''                               *
*     Phys.\ Lett.\ B {\bf 252}, 467 (1990).                           *
*     %%CITATION = PHLTA,B252,467;%%                                   *
*                                                                      *
*     Formula has been improved in order to have a                     *
*     smooth mb -> 0 limit (in agreement with topwidth.f)              *
*                                                                      *
************************************************************************
      
C formula taken from Eq.(27) of 
c---  %\cite{Czarnecki:1990kv}
c---  \bibitem{Czarnecki:1990kv} 
c---  A.~Czarnecki,
c---  %``QCD corrections to the decay t ---> W b in dimensional regularization,''
c---  Phys.\ Lett.\ B {\bf 252}, 467 (1990).
c---- %%CITATION = PHLTA,B252,467;%%
      
      include 'zerowidth.f'
      real(dp):: mb,mt,mt1,mw,om,omsq,lo,ho,be,besq,
     & Gamma0,Gamma0int,asGamma1,asGamma1int,dgauss,
     & xlo,xhi,xi,gamw,ga
      real(dp):: cachemass,cacheratio,tiny
      data cachemass,cacheratio,tiny/0._dp,0._dp,1.e-8_dp/
      common/transfer/mt1,besq,xi,ga
      save cachemass,cacheratio,tiny
!$omp threadprivate(cachemass,cacheratio,tiny,/transfer/)
      external Gamma0int,asGamma1int

c--- check to see if result has already been computed
      if (abs(mt*mw-mb-cachemass) < tiny) then
        nloratiotopdecay=cacheratio
        return
      endif

      mt1=mt
      om=mw/mt
      be=mb/mt
      ga=gamw/mw
      xi=(mt/mw)**2
      besq=be**2
      omsq=om**2

      if (zerowidth) then
        lo=Gamma0(mt,besq,omsq)
        ho=asGamma1(mt,besq,omsq)
      else
        xlo=0._dp
        xhi=(1._dp-be)**2
        lo=dgauss(Gamma0int,xlo,xhi,tiny)
        ho=dgauss(asGamma1int,xlo,xhi,tiny)     
      endif

      nloratiotopdecay=ho/lo


c--- set-up caching variables
      cachemass=mt*mw-mb
      cacheratio=nloratiotopdecay

      return
      end

