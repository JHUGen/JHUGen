      double complex function A41phiqarbmppm(j1,j2,j3,j4,za,zb)
C     arXiv:09060008v1, Eq.2.14
      implicit none 
      include 'constants.f'
      include 'nflav.f'
      include 'zprods_decl.f'
      double complex Alcphiqarbmppm,Alcphiqarbmpmp,Aslcphiqarbmppm,
     . Afphiqarbmppm
      integer j1,j2,j3,j4
      A41phiqarbmppm=
     &                      Alcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &            -2d0/xnsq*(Alcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &                      +Alcphiqarbmpmp(j1,j2,j4,j3,za,zb))
     &            -1d0/xnsq*Aslcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &    +dfloat(nflav)/xn*Afphiqarbmppm(j1,j2,j3,j4,za,zb)
      return
      end

      double complex function A41phiqarbmpmp(j1,j2,j3,j4,za,zb)
C     arXiv:09060008v1, Eq.2.15
      implicit none 
      include 'constants.f'
      include 'nflav.f'
      include 'zprods_decl.f'
      double complex Alcphiqarbmpmp,Alcphiqarbmppm,Aslcphiqarbmppm,
     . Afphiqarbmpmp
      integer j1,j2,j3,j4
      A41phiqarbmpmp=
     &                      Alcphiqarbmpmp(j1,j2,j3,j4,za,zb)
     &            -2d0/xnsq*(Alcphiqarbmpmp(j1,j2,j3,j4,za,zb)
     &                      +Alcphiqarbmppm(j1,j2,j4,j3,za,zb))
     &            +1d0/xnsq*Aslcphiqarbmppm(j1,j2,j4,j3,za,zb)
     &    +dfloat(nflav)/xn*Afphiqarbmpmp(j1,j2,j3,j4,za,zb)
      return
      end
