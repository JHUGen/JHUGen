      function a6f(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6f
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************

c---Atreepm is the amplitude for
c---q-(-p4)+Q-(-p2)+l-(-p5) ---> q+(p1)+Q+(p3)+l+(p6)
c---All outgoing particles are right-handed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: atree,virt,Lnrat
      character*2 st
      virt=epinv+Lnrat(musq,-s(j2,j3))+2._dp
c---???continuation
      a6f=atree(st,j1,j2,j3,j4,j5,j6,za,zb)*virt
      return
      end


