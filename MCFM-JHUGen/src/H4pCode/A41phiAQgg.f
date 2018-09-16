      double complex function A41phiAQggmpmm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer j1,j2,j3,j4
      double complex A1phiAQggmpmmL,A1phiAQggmpmmR,A1phiAQggmpmmF

c--- implementation of arXiv:0906.0008, Eq.(2.21)
      A41phiAQggmpmm=
     .     +A1phiAQggmpmmL(j1,j2,j3,j4,za,zb)
     .     -A1phiAQggmpmmR(j1,j2,j3,j4,za,zb)/xn**2
     .     +A1phiAQggmpmmF(j1,j2,j3,j4,za,zb)*dfloat(nflav)/xn

      return
      end
      
      double complex function A41phiAQggmppp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer j1,j2,j3,j4
      double complex A1phiAQggmpppL,A1phiAQggmpppR,A1phiAQggmpppF

c--- implementation of arXiv:0906.0008, Eq.(2.21)

      A41phiAQggmppp=
     .     +A1phiAQggmpppL(j1,j2,j3,j4,za,zb)
     .     -A1phiAQggmpppR(j1,j2,j3,j4,za,zb)/xn**2
     .     +A1phiAQggmpppF(j1,j2,j3,j4,za,zb)*dfloat(nflav)/xn

      return
      end
      
      double complex function A41phiAQggmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer j1,j2,j3,j4
      double complex A1phiAQggmpmpL,A1phiAQggmpmpR,A1phiAQggmpmpF

c--- implementation of arXiv:0906.0008, Eq.(2.21)
      A41phiAQggmpmp=
     .     +A1phiAQggmpmpL(j1,j2,j3,j4,za,zb)
     .     -A1phiAQggmpmpR(j1,j2,j3,j4,za,zb)/xn**2
     .     +A1phiAQggmpmpF(j1,j2,j3,j4,za,zb)*dfloat(nflav)/xn

      return
      end
      
      double complex function A41phiAQggmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      integer j1,j2,j3,j4
      double complex A1phiAQggmppmL,A1phiAQggmppmR,A1phiAQggmppmF

c--- implementation of arXiv:0906.0008, Eq.(2.21)
      A41phiAQggmppm=
     .     +A1phiAQggmppmL(j1,j2,j3,j4,za,zb)
     .     -A1phiAQggmppmR(j1,j2,j3,j4,za,zb)/xn**2
     .     +A1phiAQggmppmF(j1,j2,j3,j4,za,zb)*dfloat(nflav)/xn

      return
      end
      
