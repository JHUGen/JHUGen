      double complex function A51mmppp(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'epinv.f'
      include 'scale.f'
      double complex Vf,Vs,Ff,Fs,L0,L2,miA5tree,lnrat
      integer j1,j2,j3,j4,j5
C----Eq.(8,10) of hep-ph/9302280v1 of BDK multiplied by 16*pi^2*(-i)
C--- to give (16*pi^2)*(-i)*A^{[1/2]}_{5;1}
      miA5tree=za(j1,j2)**4
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*za(j5,j1))
      Vf=dcmplx(-2.5d0*epinv-2d0)
     & -0.5d0*(lnrat(musq,-s(j2,j3))+lnrat(musq,-s(j5,j1)))
      Vs=-Vf/3d0+dcmplx(2d0/9d0)
      Ff=-za(j1,j2)**2*L0(-s(j2,j3),-s(j5,j1))
     & *(za(j2,j3)*zb(j3,j4)*za(j4,j1)+za(j2,j4)*zb(j4,j5)*za(j5,j1))
     & /(2d0*za(j2,j3)*za(j3,j4)*za(j4,j5)*za(j5,j1)*s(j5,j1))
      Fs=-Ff/3d0
     & -zb(j3,j4)*za(j4,j1)*za(j2,j4)*zb(j4,j5)*L2(-s(j2,j3),-s(j5,j1))
     & *(za(j2,j3)*zb(j3,j4)*za(j4,j1)+za(j2,j4)*zb(j4,j5)*za(j5,j1))
     & /(3d0*za(j3,j4)*za(j4,j5)*s(j5,j1)**3)
     & -za(j3,j5)*zb(j3,j5)**3
     & /(3d0*zb(j1,j2)*zb(j2,j3)*za(j3,j4)*za(j4,j5)*zb(j5,j1))
     & +za(j1,j2)*zb(j3,j5)**2
     & /(3d0*zb(j2,j3)*za(j3,j4)*za(j4,j5)*zb(j5,j1))
     & +za(j1,j2)*zb(j3,j4)*za(j4,j1)*za(j2,j4)*zb(j4,j5)
     & /(6d0*s(j2,j3)*za(j3,j4)*za(j4,j5)*s(j5,j1))
      A51mmppp=-(Vf+Vs)*miA5tree-(Ff+Fs)
      
c      write(6,*) 'j1,j2,j3,j4,j5,miA5tree',j1,j2,j3,j4,j5,miA5tree
      return
      end 
