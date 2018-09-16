      function AGTYG2s(t,u,Lx,Ly,Ls)
      include 'types.f'
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.2
      real(dp):: AGTYG2s,G2sx,t,u,Lx,Ly,Ls
      AGTYG2s=+G2sx(t,u,Lx,Ly,Ls)+G2sx(u,t,Ly,Lx,Ls)
      return
      end

      function G2sx(t,u,Lx,Ly,Ls)
      include 'types.f'
      real(dp):: G2sx
      include 'zeta.f'
      real(dp)::t,u,Lx,Ly,Ls

      G2sx= (-10*Lx**4-106._dp/
     & 3*Lx**3-8*Lx**3*Ly-2*Lx**2*Ly**2-52*Lx**2*pisq
     & -44._dp/3._dp*Lx**2*Ls+6*Lx**2-40._dp/3._dp*Lx**2*Ly
     & -4*Ly**3*Lx-56._dp/
     & 3*Ly**2*Lx-32*Lx*Ly*pisq+8*Lx*Ly-80*pisq*Lx+140._dp/3._dp*Lx
     & -20._dp/3._dp*Ly**3-22._dp/3._dp*Ly**2*Ls
     & -20*Ly**2-6*Ly**2*pisq-18*pisq*Ly-22*Ly*Ls+140._dp/3._dp*Ly-4._dp
     & +154._dp/3._dp*Ls-40*pisq )*(t/u)
     & -8*Lx
     & *(Lx**3+Lx**2+4*pisq*Lx+2*pisq)*(t/u)**2
     & -2*Lx**2*(Lx**2+4*pisq)*(t/u)**3
     & +(-8*Lx**3*Ly-32*Lx**3+4*Lx**2*pisq-4*Lx**2*Ly**2-52._dp/
     & 3*Lx**2*Ly-8*Ly**2*Lx-40._dp/3._dp*Lx*Ly
     & -32*Lx*Ly*pisq-76*pisq*Lx-44._dp/
     & 3*Lx*Ls-4*Ly**4+8._dp/3._dp*Ly**3+8._dp/
     & 3*Ly**2-44._dp/3._dp*Ly**2*Ls-32*Ly**2*pisq
     & +28*Ly+4._dp-24*pisq)
!+ \t \leftrightarrow u \
      return
      end



