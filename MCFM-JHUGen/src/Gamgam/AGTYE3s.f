      function AGTYE3s(t,u,Lx,Ly,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.9
      include 'types.f'
      real(dp):: AGTYE3s,E3sx
      real(dp):: t,u,Lx,Ly,Ls
      AGTYE3s=E3sx(t,u,Lx,Ly,Ls)
     &       +E3sx(u,t,Ly,Lx,Ls)
      return
      end


      function E3sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      real(dp):: E3sx
      include 'zeta.f'
      real(dp)::t,u,Lx,Ly,Ls

      E3sx=(16._dp/9*Lx**3+ (-76._dp/9+8._dp/3._dp*Ls)*Lx**2
     & +16._dp/9*pisq*Lx+8._dp/9*Ly**3
     & +(4._dp/3._dp*Ls-2._dp/9._dp)*Ly**2
     & +(8._dp/9*pisq+4*Ls-10._dp)*Ly
     & -1._dp/3._dp*pisq*Ls-202._dp/27*Ls+19._dp/9._dp*pisq
     & -2._dp/9*zeta3+3401._dp/162._dp)*(t/u)
     & + (16._dp/9*pisq*Lx+16._dp/9*Ly**3+
     & (8._dp/3._dp*Ls-52._dp/9._dp)*Ly**2
     & + (-76._dp/9._dp+8._dp/3._dp*Ls )*Ly+8._dp/9._dp*pisq )
!    + \t \leftrightarrow u \
      return
      end

