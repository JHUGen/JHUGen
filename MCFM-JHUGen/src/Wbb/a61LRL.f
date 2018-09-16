      function a61LRL(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a61LRL

C---  Corresponds to all outgoing
*     q(j1,-)+Q(j3,+)+e(j6,-)+q~(j4)+Q~(j2)+e~(j5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a61
      real(dp):: prop
      prop=s(5,6)/sqrt((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
c---Note interchange of za,zb to effect complex conjugation
      a61LRL=a61('pp',j1,j2,j3,j4,j5,j6,zb,za)*prop
      return
      end

