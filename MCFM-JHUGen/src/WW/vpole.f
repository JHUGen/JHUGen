      function Vpole(sij)
      implicit none
      include 'types.f'
      complex(dp):: Vpole
      
c---  DKS Eq. 2.12
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      real(dp):: sij
      complex(dp):: Lnrat,xl12
        
      xl12=Lnrat(-sij,musq)

      Vpole=-epinv*epinv2+epinv*(-1.5_dp+xl12)
     &   -0.5_dp*xl12**2+1.5_dp*xl12-3.5_dp

      return
      end
