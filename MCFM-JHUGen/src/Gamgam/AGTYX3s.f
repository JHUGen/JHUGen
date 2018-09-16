      function AGTYX3s(t,u,Lx,Ly,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.6
      include 'types.f'
      real(dp):: AGTYX3s,X3sx,t,u,Lx,Ly,Ls
      AGTYX3s=+X3sx(t,u,Lx,Ly,Ls)+X3sx(u,t,Ly,Lx,Ls)
      return
      end

      function X3sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      real(dp):: X3sx
      real(dp)::t,u,Lx,Ly,Ls

      X3sx=1._dp/18._dp*(2*Ls+Ly+Lx)**2*(t/u)
! + \t \leftrightarrow u \
      return
      end
