      function AGTYX5s(t,u,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.8
      include 'types.f'
      real(dp):: AGTYX5s,X5sx,t,u,Ls
      AGTYX5s=X5sx(t,u,Ls)+X5sx(u,t,Ls)
      return
      end

      function X5sx(t,u,Ls)
      implicit none
      include 'types.f'
      include 'zeta.f'
      real(dp):: X5sx
      real(dp)::t,u,Ls

      X5sx= (32._dp/9._dp*Ls**2+32._dp/9._dp*pisq)*t/u
!       + t leftrightarrow u
      return
      end
