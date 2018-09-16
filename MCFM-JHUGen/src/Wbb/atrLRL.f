      function atrLRL(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: atrLRL

c---atrLRL is the amplitude for
c---q+(-p4)+Q-(-p2)+l+(-p5) ---> q-(p1)+Q+(p3)+l-(p6)
c---All outgoing particles are except p3 left-handed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: atree
      real(dp):: prop
      prop=s(5,6)/sqrt((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
c---Note interchange of za,zb to effect complex conjugation
      atrLRL=atree('pp',j1,j2,j3,j4,j5,j6,zb,za)*prop
c      write(6,*) 'atrLRL',atrLRL
      return
      end
