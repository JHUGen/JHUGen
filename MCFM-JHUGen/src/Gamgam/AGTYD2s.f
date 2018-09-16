      function AGTYD2s(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.8
      implicit none
      include 'types.f'
      real(dp):: AGTYD2s,D2sx
      real(dp):: s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv

      AGTYD2s=
     &  D2sx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li4z)
     & +D2sx(s,u,t,Ly,Lx,Ls,Li2y,Li3y,Li4y,Li2x,Li3x,Li4x,Li4zinv)
      return
      end

      function D2sx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z)
      implicit none
      include 'types.f'
      real(dp):: D2sx
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Ls,
     & Li4x,Li4y,Li4z,Li3x,Li3y,Li2x,Li2y
      D2sx=
     & (48*Li4z-16*Li4x+24*Li4y+
     & (56*Ly-8*Lx+20._dp)*Li3x
     & +(8*Lx-12._dp+16*Ly )*Li3y
     & +(16._dp/3._dp*pisq-20*Lx-12*Ly-8*Lx**2+4*Ly**2)*Li2x
     & +1._dp/3._dp*Lx**4+ (-8*Ly-70._dp/9
     & )*Lx**3+ (-4*pisq+286._dp/
     & 9-16*Ly+14*Ly**2-44._dp/3._dp*Ls )*Lx**2
     & + (-22._dp/
     & 9*pisq+4*Ly**3-8*pisq*Ly-6*Ly**2
     & )*Lx-44._dp/9*Ly**3+ (-4._dp/3._dp*pisq+35._dp/
     & 9-22._dp/3._dp*Ls )*Ly**2
     & +(57._dp-26._dp/9*pisq-72*zeta3-22*Ls)*Ly
     & +479._dp/9*zeta3+19._dp/60*pi**4
     & -52*zeta3*Ls+1141._dp/27*Ls-215._dp/18._dp*pisq
     & -43417._dp/324._dp+23._dp/6._dp*pisq*Ls)*(t/u)
     & +(6*Lx**2+(-12._dp-12*Ly)*Lx+6*Ly**2+12*Ly+6*pisq)*(t/s)**2
     & -6*Lx**2*(t/u)**2
     & +(16*Li4y+48*Li3x*Ly+64*Li3y-8*Ly**2*Li2y-64*Li2x*Lx
     & -4._dp/3._dp*Lx**4+ (-20._dp/3._dp*pisq+6*Ly**2
     & )*Lx**2+ (-24*Ly**2+ (-16._dp/
     & 3*pisq-14._dp)*Ly-148._dp/9._dp*pisq)*Lx
     & -112._dp/9*Ly**3+ (-44._dp/3._dp*Ls+298._dp/9._dp)*Ly**2
     & +(538._dp/9._dp-48*zeta3-44._dp/3._dp*Ls)*Ly
     & -8*zeta3-1._dp/3._dp*pi**4+61._dp/9*pisq )

!     & & & + \t \leftrightarrow u \
      return
      end

