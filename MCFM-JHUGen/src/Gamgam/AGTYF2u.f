      function AGTYF2u(s,t,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.22
      include 'types.f'
      include 'zeta.f'
      real(dp):: AGTYF2u,s,t,Lu

      AGTYF2u=(4._dp*pisq/27._dp+32._dp*Lu**2/9._dp-160._dp*Lu/27._dp)
     & *(t**2 + s**2)/(s*t)
      return
      end
