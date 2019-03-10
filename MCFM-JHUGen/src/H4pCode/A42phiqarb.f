      double complex function A42phiqarbmppm(j1,j2,j3,j4,za,zb)
C     arXiv:09060008v1, Eq.2.16
      implicit none 
      include 'constants.f'
      include 'nflav.f'
      include 'zprods_decl.f'
      double complex Alcphiqarbmppm,Alcphiqarbmpmp,Aslcphiqarbmppm,
     . Afphiqarbmppm
      integer j1,j2,j3,j4
      A42phiqarbmppm=
     &           Alcphiqarbmpmp(j1,j2,j4,j3,za,zb)
     & +1d0/xnsq*(Alcphiqarbmpmp(j1,j2,j4,j3,za,zb)
     &           +Alcphiqarbmppm(j1,j2,j3,j4,za,zb))
     & +1d0/xnsq*Aslcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &    -dfloat(nflav)/xn*Afphiqarbmppm(j1,j2,j3,j4,za,zb)
      return
      end

      double complex function A42phiqarbmpmp(j1,j2,j3,j4,za,zb)
C     based on arXiv:09060008v1, Eq.2.17
c       but with slight modifications traced back to equation 2.8 of
c        Z. Bern, L. Dixon, D. Kosower & S. Weinzierl, hep-ph/9610370
c       
c        a) interchanged 3- and 4+ in first term
c        b) changed sign of second line [ i.e. +1/N^2*(Alc+Alc)]

      implicit none 
      include 'constants.f'
      include 'nflav.f'
      include 'zprods_decl.f'
      double complex Alcphiqarbmpmp,Alcphiqarbmppm,Aslcphiqarbmppm,
     . Afphiqarbmpmp
      integer j1,j2,j3,j4
      A42phiqarbmpmp=
     &           Alcphiqarbmppm(j1,j2,j4,j3,za,zb)
     & +1d0/xnsq*(Alcphiqarbmpmp(j1,j2,j3,j4,za,zb)
     &           +Alcphiqarbmppm(j1,j2,j4,j3,za,zb))
     & -1d0/xnsq*Aslcphiqarbmppm(j1,j2,j4,j3,za,zb)
     &    -dfloat(nflav)/xn*Afphiqarbmpmp(j1,j2,j3,j4,za,zb)
      return
      end

