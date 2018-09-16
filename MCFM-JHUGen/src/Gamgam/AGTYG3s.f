      function AGTYG3s(t,u,Lx,Ly,Ls)
      include 'types.f'
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.3
      real(dp):: AGTYG3s,G3sx,t,u,Lx,Ly,Ls
      AGTYG3s=+G3sx(t,u,Lx,Ly,Ls)+G3sx(u,t,Ly,Lx,Ls)
      return
      end

      function G3sx(t,u,Lx,Ly,Ls)
      include 'types.f'
      real(dp):: G3sx
      include 'zeta.f'
      real(dp)::t,u,Lx,Ly,Ls

      G3sx= (2*Lx**4+2*Lx**3*Ly+22._dp/
     & 3*Lx**3+11._dp/
     & 3*Lx**2*Ls+2*Lx**2*Ly**2+10*Lx**2*Ly+68._dp/
     & 9*Lx**2+13*Lx**2*pisq
     & +20._dp/3._dp*Ly**2*Lx+6*Lx*Ly*pisq+100._dp/
     & 9*Lx*Ly+22._dp/3._dp*Lx*Ly*Ls+110._dp/
     & 9*Lx*Ls+50._dp/3._dp*pisq*Lx+50._dp/
     & 9*Ly**2
     & +2*Ly**2*pisq+8._dp/3._dp*pisq*Ly+110._dp/
     & 9*Ly*Ls+2._dp+121._dp/18*Ls**2-11._dp/3*pisq*Ls
     & +1._dp/ 2*pi**4+13._dp/ 2*pisq)*(t/u)
     & +2*Lx*(Lx**3+Lx**2+4*pisq*Lx+2*pisq)*(t/u)**2
     & +1._dp/2._dp*Lx**2*(Lx**2+4*pisq )*(t/u)**3
     & + (4*Lx**3*Ly+8*Lx**3-2*Lx**2*pisq+26._dp/
     & 3*Lx**2*Ly+2*Ly**2*Lx+8*Lx*Ly*pisq+20._dp/
     & 3*Lx*Ly+22._dp/3._dp*Lx*Ls
     & +18*pisq*Lx-4._dp/3._dp*Ly**3+20._dp/3._dp*Ly**2+22._dp/
     & 3*Ly**2*Ls+8*Ly**2*pisq+6*pisq )

!  + \t \leftrightarrow u \
      return
      end
