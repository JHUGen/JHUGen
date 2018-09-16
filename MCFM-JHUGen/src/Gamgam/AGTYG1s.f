      function AGTYG1s(t,u,Lx,Ly)
      include 'types.f'
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.1
      real(dp):: AGTYG1s,G1sx,t,u,Lx,Ly
      AGTYG1s=+G1sx(t,u,Lx,Ly)+G1sx(u,t,Ly,Lx)
      return
      end

      function G1sx(t,u,Lx,Ly)
      implicit none
      include 'types.f'
      real(dp):: G1sx
      include 'zeta.f'
      real(dp)::t,u,Lx,Ly

      G1sx=
     & (14*Lx**4+28*Lx**3+8*Lx**2*Ly**2+56*Lx**2*pisq-48*Lx**2
     & +12*Lx**2*Ly+32*Lx*Ly*pisq+80*pisq*Lx
     & +2*Ly**4+12*Ly**3-10*Ly**2+8*Ly**2*pisq+26*pisq+24*pisq*Ly
     * -84*Ly+102._dp)*(t/u)
     & +8*Lx*(Lx**3+Lx**2+4*pisq*Lx+2*pisq)*(t/u)**2
     & +2*Lx**2*(Lx**2+4*pisq)*(t/u)**3
     & +(32*Lx**3+8*Lx**2*Ly**2+80*pisq*Lx+32*Lx*Ly*pisq+8*Ly**2*Lx
     & +8*Ly**4+32*Ly**2*pisq-32*Ly**2-4._dp-56*Ly+24*pisq)
!     + \t \leftrightarrow u \
      return
      end
