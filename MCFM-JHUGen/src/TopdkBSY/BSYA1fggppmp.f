      double complex function BSYA1fggppmp(e1,p2,p3,e4,za,zb,zab,zba)
      implicit none
C-----Authors: John Campbell and Keith Ellis, November 2011
C---- arXiv:1101.5947 [hep-ph], (limit of) Eq. (101)
C---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'constants.f'
      include 'zprods_decl.f'
      include 'zabprods_decl.f'
      integer e1,p2,p3,e4
 
      BSYA1fggppmp=czip
      return
      end

