      function a6g(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6g

************************************************************************
*     Author: R.K. Ellis                                               *
*     August, 1999.                                                    *
*     implementation of Eq. (6.1) of BDKW hep-ph/9708239               *
*     character string st can take various values                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      character*9 st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6treeg,vvg,fcc,fsc

      a6g=a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)*vvg(st,j1,j2,j3,j4,j5,j6)
      a6g=a6g+fcc(st,j1,j2,j3,j4,j5,j6,za,zb)
     &       +fsc(st,j1,j2,j3,j4,j5,j6,za,zb)

      end

