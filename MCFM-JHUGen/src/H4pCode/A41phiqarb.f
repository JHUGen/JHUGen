      function A41phiqarbmppm(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41phiqarbmppm
C     arXiv:09060008v1, Eq.2.14
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nflav.f'
      include 'zprods_decl.f'
      complex(dp):: Alcphiqarbmppm,Alcphiqarbmpmp,Aslcphiqarbmppm,
     & Afphiqarbmppm
      integer:: j1,j2,j3,j4
      A41phiqarbmppm=
     &                      Alcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &            -2._dp/xnsq*(Alcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &                      +Alcphiqarbmpmp(j1,j2,j4,j3,za,zb))
     &            -1._dp/xnsq*Aslcphiqarbmppm(j1,j2,j3,j4,za,zb)
     &    +real(nflav,dp)/xn*Afphiqarbmppm(j1,j2,j3,j4,za,zb)
      return
      end

      function A41phiqarbmpmp(j1,j2,j3,j4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A41phiqarbmpmp
C     arXiv:09060008v1, Eq.2.15
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nflav.f'
      include 'zprods_decl.f'
      complex(dp):: Alcphiqarbmpmp,Alcphiqarbmppm,Aslcphiqarbmppm,
     & Afphiqarbmpmp
      integer:: j1,j2,j3,j4
      A41phiqarbmpmp=
     &                      Alcphiqarbmpmp(j1,j2,j3,j4,za,zb)
     &            -2._dp/xnsq*(Alcphiqarbmpmp(j1,j2,j3,j4,za,zb)
     &                      +Alcphiqarbmppm(j1,j2,j4,j3,za,zb))
     &            +1._dp/xnsq*Aslcphiqarbmppm(j1,j2,j4,j3,za,zb)
     &    +real(nflav,dp)/xn*Afphiqarbmpmp(j1,j2,j3,j4,za,zb)
      return
      end
