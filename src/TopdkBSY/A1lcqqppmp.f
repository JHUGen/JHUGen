      double complex function A1lcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- taken from arXiv:1101.5947 [hep-ph], Eq. (55)
      include 'constants.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      double precision mt2
      double complex BSYA0qqppmp,BSYA1lcqqppmp,Vlc,lnrat
      integer e1,p2,p3,e4
      mt2=mt**2 
c--- NB: added a minus sign in front of the double pole compared to Eq. (55)
      Vlc=-dcmplx(epinv**2)+dcmplx(epinv)*(8d0/3d0
     & -(lnrat(musq,-s(1,p2))+lnrat(mt2,-s(1,p2))))
      A1lcqqppmp=Vlc*BSYA0qqppmp(e1,p2,p3,e4,za,zb,zab,zba)
     &              +BSYA1lcqqppmp(e1,p2,p3,e4,za,zb,zab,zba)  
      
      return
      end
