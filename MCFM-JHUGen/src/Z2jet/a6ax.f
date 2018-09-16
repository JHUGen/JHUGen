      function a6ax(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'types.f'
      complex(dp):: a6ax
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: fmt,fzip
C     Implementation of BDKW Nucl Phys. 513, 3 (1998)
C     Eqn. (12.13). Their function is multipled by i*4*pi^2
      a6ax=
     &    +2._dp*(fmt(s(j1,j2),s(j3,j4),s(j5,j6))
     &     -fzip(s(j1,j2),s(j3,j4),s(j5,j6)))/s(j5,j6)
     & *(zb(j6,j3)*za(j4,j2)*za(j2,j5)/za(j1,j2)
     & -zb(j6,j1)*zb(j1,j3)*za(j4,j5)/zb(j1,j2))
     &    +2._dp*(fmt(s(j3,j4),s(j1,j2),s(j5,j6))
     &     -fzip(s(j3,j4),s(j1,j2),s(j5,j6)))/s(j5,j6)
     & *(zb(j6,j1)*za(j2,j4)*za(j4,j5)/za(j3,j4)
     & -zb(j6,j3)*zb(j3,j1)*za(j2,j5)/zb(j3,j4))
      return
      end

