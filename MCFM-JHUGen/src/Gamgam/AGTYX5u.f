      function AGTYX5u(s,t,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.16
      include 'types.f'
      real(dp):: AGTYX5u
      real(dp)::s,t,Lu
      AGTYX5u=32._dp/9._dp*Lu**2*(t**2+s**2)/(s*t)
      return
      end
