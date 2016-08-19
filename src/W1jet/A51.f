      double complex function A51(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer j1,j2,j3,j4,j5
      double complex Vcc,Fcc,Vsc,Fsc,l12,l23,L0,L1,Lsm1,A5lom
      double complex lnrat

C    -i * A5tree  
      A5lom =-za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5))
      l12=lnrat(musq,-s(j1,j2))
      l23=lnrat(musq,-s(j2,j3))

C--leading N
      Vcc=
     . -(epinv**2+epinv*l12+0.5d0*l12**2)
     . -(epinv**2+epinv*l23+0.5d0*l23**2)
     . -2d0*(epinv+l23)-4d0      

      Fcc=za(j3,j4)**2/(za(j1,j2)*za(j2,j3)*za(j4,j5))
     . *(Lsm1(-s(j1,j2),-s(j4,j5),-s(j2,j3),-s(j4,j5))
     . -2d0*za(j3,j1)*zb(j1,j5)*za(j5,j4)/za(j3,j4)
     .   *L0(-s(j2,j3),-s(j4,j5))/s(j4,j5))

      Vsc =0.5d0*(epinv+l23)+1d0
      Fsc =za(j3,j4)*za(j3,j1)*zb(j1,j5)*za(j5,j4)
     . /(za(j1,j2)*za(j2,j3)*za(j4,j5))*L0(-s(j2,j3),-s(j4,j5))/s(j4,j5)
     . +0.5d0*(za(j3,j1)*zb(j1,j5))**2*za(j4,j5)
     . /(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s(j4,j5))/s(j4,j5)**2

      A51=(Vcc+Vsc)*A5Lom+Fcc+Fsc

      return  
      end
