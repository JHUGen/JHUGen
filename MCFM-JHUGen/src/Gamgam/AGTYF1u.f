      function AGTYF1u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.18
      include 'types.f'
      real(dp):: AGTYF1u
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYF1u=(5._dp/36._dp*Lx**2+ (1._dp/3*Lu
     & -10._dp/27-1._dp/3._dp*Ly )*Lx
     & +1._dp/3._dp*Ly**2+ (20._dp/27-2._dp/3._dp*Lu )*Ly
     & -20._dp/27._dp*Lu-13._dp/108._dp*pisq+1._dp/3._dp*Lu**2)
     & *(t**2+s**2)/(s*t)
      return
      end
