      function AGTYX3u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.14
      include 'types.f'
      real(dp):: AGTYX3u
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYX3u= (1._dp/ 18*Lx**2+ (-2._dp/ 9*Ly+2._dp/ 9*Lu )*Lx
     & -4._dp/9*Ly*Lu+1._dp/18*pisq+2._dp/9*Lu**2+2._dp/9*Ly**2)
     & *(t**2+s**2)/(s*t)
      return
      end
