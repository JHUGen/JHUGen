      double complex function A52(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      integer j1,j2,j3,j4,j5
      double complex Vcc,Fcc,Vsc,Fsc,l12,l45,L0,L1,Lsm1,A5lom
      double complex lnrat

      l12=lnrat(musq,-s(j1,j2))  
      l45=lnrat(musq,-s(j4,j5))  
C    -i * A5tree  
      A5lom=za(j2,j4)**2/(za(j2,j3)*za(j3,j1)*za(j4,j5))
      Vcc=-(epinv**2+epinv*l12+0.5d0*l12**2)
     . -2d0*(epinv+l45)-4d0

C--subleading N
      Fcc=-za(j2,j4)**2/(za(j2,j3)*za(j3,j1)*za(j4,j5))
     . *Lsm1(-s(j1,j2),-s(j4,j5),-s(j1,j3), -s(j4,j5))
     . +za(j2,j4)*(za(j1,j2)*za(j3,j4)-za(j1,j4)*za(j2,j3))
     .  /(za(j2,j3)*za(j1,j3)**2*za(j4,j5))
     . *Lsm1(-s(j1,j2),-s(j4,j5),-s(j2,j3),-s(j4,j5))
     . +2d0*zb(j1,j3)*za(j1,j4)*za(j2,j4)/(za(j1,j3)*za(j4,j5))
     . *L0(-s(j2,j3),-s(j4,j5))/s(j4,j5)

      Vsc=0.5d0*(epinv+l45)+0.5d0
      Fsc=za(j1,j4)**2*za(j2,j3)/(za(j1,j3)**3*za(j4,j5))
     . *Lsm1(-s(j1,j2),-s(j4,j5),-s(j2,j3),-s(j4,j5))
     . -0.5d0*(za(j4,j1)*zb(j1,j3))**2*za(j2,j3)/(za(j1,j3)*za(j4,j5))
     . *L1(-s(j4,j5),-s(j2,j3))/s(j2,j3)**2
     . +za(j1,j4)**2*za(j2,j3)*zb(j3,j1)/(za(j1,j3)**2*za(j4,j5))
     . *L0(-s(j4,j5),-s(j2,j3))/s(j2,j3)
     . -za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j5)/za(j1,j3)
     . *L1(-s(j4,j5),-s(j1,j2))/s(j1,j2)**2
     . -za(j2,j1)*zb(j1,j3)*za(j3,j4)*za(j1,j4)/(za(j1,j3)**2*za(j4,j5))
     . *L0(-s(j4,j5),-s(j1,j2))/s(j1,j2)
     . -0.5d0*zb(j3,j5)*(zb(j1,j3)*zb(j2,j5)+zb(j2,j3)*zb(j1,j5))
     . /(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*zb(j4,j5))

      A52=(Vcc+Vsc)*A5lom+Fcc+Fsc

      return  
      end


