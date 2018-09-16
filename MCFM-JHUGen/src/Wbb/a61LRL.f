      double complex function a61LRL(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
C---  Corresponds to all outgoing
*     q(j1,-)+Q(j3,+)+e(j6,-)+q~(j4)+Q~(j2)+e~(j5)
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer j1,j2,j3,j4,j5,j6
      double complex a61
      double precision prop
      prop=s(5,6)/dsqrt((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
c---Note interchange of za,zb to effect complex conjugation
      a61LRL=a61('pp',j1,j2,j3,j4,j5,j6,zb,za)*prop
      return
      end

