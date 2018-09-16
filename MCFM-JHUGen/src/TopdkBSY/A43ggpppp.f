      function A43ggpppp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: A43ggpppp

C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken form arXiv:1101.5947 [hep-ph], Eq. (5,52)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      real(dp):: mt2,beta,xlog
      complex(dp):: BSYALggpppphp,BSYARggpppphp,BSYALslggpppp
      complex(dp):: BSYA0ggpppp,lnrat,Vslc
      integer:: e1,p2,p3,e4

      mt2=mt**2
      beta=sqrt(1d0-4d0*mt2/s(p2,p3))
      xlog=log((1d0-beta)/(1d0+beta))/beta
C---- arXiv:1101.5947 [hep-ph], Eq. (52)
c----   with an extra minus sign to account for A4;3 -> -Vslc, c.f. Eq. (53))
      Vslc=-cplx1(epinv)*
     & (+s(1,p2)/s(p2,p3)*(lnrat(musq,-s(1,p2))+lnrat(mt2,-s(1,p2)))
     &  +s(1,p3)/s(p2,p3)*(lnrat(musq,-s(1,p3))+lnrat(mt2,-s(1,p3)))
     & +lnrat(musq,-s(p2,p3))+(1d0-2d0*mt2/s(p2,p3))*xlog)

      A43ggpppp=(
     &  Vslc*(BSYA0ggpppp(e1,p2,p3,e4,za,zb,zab,zba)
     &       +BSYA0ggpppp(e1,p3,p2,e4,za,zb,zab,zba))
     & +BSYALggpppphp(e1,p2,p3,e4,za,zb,zab,zba)
     & +BSYALggpppphp(e1,p3,p2,e4,za,zb,zab,zba)
     & +BSYALslggpppp(e1,p2,e4,p3,za,zb,zab,zba)
     & +BSYALslggpppp(e1,p3,e4,p2,za,zb,zab,zba)
     & +BSYARggpppphp(e1,p2,p3,e4,za,zb,zab,zba)
     & +BSYARggpppphp(e1,p3,p2,e4,za,zb,zab,zba))

      return
      end

