      function A51mpmpp(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      complex(dp):: A51mpmpp
      complex(dp):: Vf,Vs,Ff,Fs,L0,L1,L2,Ls1,miA5tree,lnrat
      integer:: j1,j2,j3,j4,j5
C----Eq.(9,10) of hep-ph/9302280v1 of BDK multiplied by 16*pi^2*(-i)
C--- to give (16*pi^2)*(-i)*A^{[1/2]}_{5;1}
      miA5tree=za(j1,j3)**4
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*za(j5,j1))
      Vf=cplx1(-five/two*epinv-two)
     & -half*(lnrat(musq,-s(j3,j4))+lnrat(musq,-s(j5,j1)))
      Vs=-Vf/three+cplx1(two/nine)
      Ff=
     & -za(j1,j3)**2*za(j4,j1)*zb(j2,j4)**2
     & /(za(j4,j5)*za(j5,j1)*s(j5,j1)**2)
     & *Ls1(-s(j2,j3),-s(j5,j1),-s(j3,j4),-s(j5,j1))

     & +za(j1,j3)**2*za(j5,j3)*zb(j2,j5)**2
     & /(za(j3,j4)*za(j4,j5)*s(j3,j4)**2)
     & *Ls1(-s(j1,j2),-s(j3,j4),-s(j5,j1),-s(j3,j4))

     & -za(j1,j3)**3*L0(-s(j3,j4),-s(j5,j1))
     & *(za(j1,j5)*zb(j5,j2)*za(j2,j3)-za(j3,j4)*zb(j4,j2)*za(j2,j1))
     & /(two*za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*za(j5,j1)*s(j5,j1))
      Fs=
     & -za(j1,j2)*za(j2,j3)*za(j3,j4)*(za(j4,j1)*zb(j2,j4))**2
     & /(za(j4,j5)*za(j5,j1)*za(j2,j4)**2*s(j5,j1)**2)
     & *(two*Ls1(-s(j2,j3),-s(j5,j1),-s(j3,j4),-s(j5,j1))
     & +L1(-s(j2,j3),-s(j5,j1))+L1(-s(j3,j4),-s(j5,j1)))

     & +za(j3,j2)*za(j2,j1)*za(j1,j5)*(za(j5,j3)*zb(j2,j5))**2
     & /(za(j5,j4)*za(j4,j3)*za(j2,j5)**2*s(j3,j4)**2)
     & *(two*Ls1(-s(j1,j2),-s(j3,j4),-s(j5,j1),-s(j3,j4))
     & +L1(-s(j1,j2),-s(j3,j4))+L1(-s(j5,j1),-s(j3,j4)))

     &+two*za(j2,j3)**2*(za(j4,j1)*zb(j2,j4))**3*L2(-s(j2,j3),-s(j5,j1))
     & /(three*za(j4,j5)*za(j5,j1)*za(j2,j4)*s(j5,j1)**3)
     &-two*za(j2,j1)**2*(za(j5,j3)*zb(j2,j5))**3*L2(-s(j1,j2),-s(j3,j4))
     & /(three*za(j5,j4)*za(j4,j3)*za(j2,j5)*s(j3,j4)**3)

     & +L2(-s(j3,j4),-s(j5,j1))/s(j5,j1)**3
     * *(za(j1,j3)*zb(j2,j4)*zb(j2,j5)/(three*za(j4,j5))
     & *(za(j1,j5)*zb(j5,j2)*za(j2,j3)-za(j3,j4)*zb(j4,j2)*za(j2,j1))
     & +two*(za(j1,j2)*za(j3,j4))**2*za(j4,j1)*zb(j2,j4)**3
     & /(three*za(j4,j5)*za(j5,j1)*za(j2,j4))
     & -two*(za(j3,j2)*za(j1,j5))**2*za(j5,j3)*zb(j2,j5)**3
     & /(three*za(j5,j4)*za(j4,j3)*za(j2,j5)))

     & +za(j1,j3)**3*L0(-s(j3,j4),-s(j5,j1))
     & /(six*za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*za(j5,j1)*s(j5,j1))
     & *(za(j1,j5)*zb(j5,j2)*za(j2,j3)-za(j3,j4)*zb(j4,j2)*za(j2,j1))

     & +(zb(j2,j4)*zb(j2,j5))**2
     & /(three*zb(j1,j2)*zb(j2,j3)*zb(j3,j4)*za(j4,j5)*zb(j5,j1))

     & -za(j1,j2)*za(j4,j1)**2*zb(j2,j4)**3
     & /(three*za(j4,j5)*za(j5,j1)*za(j2,j4)
     & *zb(j2,j3)*zb(j3,j4)*s(j5,j1))

     & +za(j3,j2)*za(j5,j3)**2*zb(j2,j5)**3
     & /(three*za(j5,j4)*za(j4,j3)*za(j2,j5)
     & *zb(j2,j1)*zb(j1,j5)*s(j3,j4))

     & +za(j1,j3)**2*zb(j2,j4)*zb(j2,j5)
     & /(six*s(j3,j4)*za(j4,j5)*s(j5,j1))

      A51mpmpp=-(Vf+Vs)*miA5tree-(Ff+Fs)
      return
      end
