      double complex function A1slcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (56)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      double precision mt2,beta,xlog
      double complex BSYA0qqppmp,BSYA1slcqqppmp,Vslc,lnrat
      integer e1,p2,p3,e4
      mt2=mt**2 
      beta=sqrt(1d0-4d0*mt2/s(p2,p3))
      xlog=log((1d0-beta)/(1d0+beta))/beta
c--- NB: added a minus sign in front of the double pole compared to Eq. (56)
      Vslc=-dcmplx(epinv**2)-dcmplx(epinv)*(1d0
     & +lnrat(musq,-s(p2,p3))+(1d0-2d0*mt2/s(p2,p3))*xlog)
      A1slcqqppmp=Vslc*BSYA0qqppmp(e1,p2,p3,e4,za,zb,zab,zba)
     &                +BSYA1slcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)  
      
      return
      end
