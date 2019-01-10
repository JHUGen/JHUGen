      double complex function a6virtLL(st,j1,j2,j3,j4,j5,j6,za,zb)
****************************************************
* virtual amplitude for
* 0->q(p1)+qb(p2)+gam(p3)+gam(p4)+lb*p5)+l(p6)
* where both photon are coming from lepton line
****************************************************
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*9 st
      integer j1,j2,j3,j4,j5,j6
      double complex a6treeg,vll
c-----
      a6virtLL = a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)
     .              *vll(st,j5,j6,j3,j4,j1,j2)
c-----
      return
      end


      double complex function vll(st,j1,j2,j3,j4,j5,j6)
c-----divergent part
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'epinv2.f'
      integer j1,j2,j3,j4,j5,j6
      double complex Lnrat,xl12
      character*9 st
c-----
      xl12=Lnrat(musq,-s(j1,j2))
c-----
      vll = -4D0
     .      -(epinv*epinv2+epinv*xl12+half*xl12**2)
     .      -(3D0/2D0)*(epinv+xl12)
c-----
      return
      end
