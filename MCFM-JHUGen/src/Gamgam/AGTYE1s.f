      function AGTYE1s(t,u,Lx,Ly,Ls)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.5
      include 'types.f'
      real(dp):: AGTYE1s,E1sx
      real(dp):: t,u,Lx,Ly,Ls
      AGTYE1s=E1sx(t,u,Lx,Ly,Ls)
     &       +E1sx(u,t,Ly,Lx,Ls)
      return
      end


      function E1sx(t,u,Lx,Ly,Ls)
      implicit none
      include 'types.f'
      real(dp):: E1sx
      include 'zeta.f'
      real(dp)::t,u,Lx,Ly,Ls

      E1sx= (22._dp/ 9*Lx**3+ (-76._dp/
     & 9+4*Ls+2._dp/3._dp*Ly )*Lx**2+ (1._dp/
     & 3*Ly**2+Ly+16._dp/ 9*pisq-7._dp/3._dp )*Lx
     & +11._dp/ 9*Ly**3+ (7._dp/ 9+2*Ls
     & )*Ly**2+ (6*Ls+8._dp/ 9*pisq-37._dp/3._dp
     & )*Ly
     & +19._dp/ 9*pisq-328._dp/ 27*Ls+3401._dp/ 162-2._dp/
     & 9*zeta3-1._dp/3._dp*pisq*Ls )*(t/u)
     & + (-46._dp/ 9*Lx**2+ (2._dp/3._dp*Ly+2._dp/
     & 3*Ly**2-76._dp/ 9 )*Lx
     & +22._dp/ 9*Ly**3+4*Ly**2*Ls+ (16._dp/
     & 9*pisq+4*Ls )*Ly+8._dp/ 9*pisq )
!+ \t \leftrightarrow u \
      return
      end
