      double complex function Vpole(sij)
      implicit none
c---  DKS Eq. 2.12
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      double precision sij
      double complex Lnrat,xl12
        
      xl12=Lnrat(-sij,musq)

      Vpole=-epinv*epinv2+epinv*(-1.5d0+xl12)
     .   -0.5d0*xl12**2+1.5d0*xl12-3.5d0

      return
      end
