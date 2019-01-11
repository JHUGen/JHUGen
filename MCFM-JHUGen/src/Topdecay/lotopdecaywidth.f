      double precision function lotopdecaywidth(mt,mb,mw,gamw)
************************************************************************
*     Authors: R.K. Ellis and J. Campbell, February 2012               *
*                                                                      *
*     LO width of the top quark, including the effect of the           *
*     bottom quark mass and off shellness                              *
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
*                                                                      *
************************************************************************
      implicit none
      include 'zerowidth.f'
      double precision mb,mt,mt1,mw,om,omsq,be,besq,dgauss,Gamma0,
     & xlo,xhi,xi,gamw,ga,Gamma0int
      double precision cachemass,cachewidth,tiny
      data cachemass,cachewidth,tiny/0d0,0d0,1d-8/
      common/transfer/mt1,besq,xi,ga
      save cachemass,cachewidth,tiny
!$omp threadprivate(cachemass,cachewidth,tiny,/transfer/)
      external Gamma0int


c--- check to see if result has already been computed
      if (abs(mt*mw-mb-cachemass) .lt. tiny) then
        lotopdecaywidth=cachewidth
        return
      endif
       
      mt1=mt
      om=mw/mt
      be=mb/mt
      besq=be**2

      if (zerowidth) then
      omsq=om**2
      lotopdecaywidth=Gamma0(mt,besq,omsq)
      else
      xlo=0d0
      xhi=(1d0-be)**2
      ga=gamw/mw
      xi=(mt/mw)**2
      lotopdecaywidth=dgauss(Gamma0int,xlo,xhi,tiny)
      endif

c--- set-up caching variables
      cachemass=mt*mw-mb
      cachewidth=lotopdecaywidth

      return
      end



