      function AGTYD1s(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.4
      implicit none
      include 'types.f'
      real(dp):: AGTYD1s,D1sx
      real(dp):: s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv

      AGTYD1s=
     &  D1sx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li4z)
     & +D1sx(s,u,t,Ly,Lx,Ls,Li2y,Li3y,Li4y,Li2x,Li3x,Li4x,Li4zinv)
      return
      end

      function D1sx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z)
      implicit none
      include 'types.f'
      real(dp):: D1sx
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Ls,
     & Li4x,Li4y,Li4z,Li3x,Li3y,Li2x,Li2y
      D1sx=
     & (96*Li4z-48*Li4x+52*Li4y+
     & (124*Ly-8*Lx+46._dp)*Li3x
     & +(-16*Lx**2+ (-46._dp+8*Ly)*Lx
     & +6*Ly**2-30*Ly+4._dp/3._dp*pisq)*Li2x
     & + (-30._dp+28*Ly+36*Lx )*Li3y+1._dp/
     & 2*Lx**4+ (-56._dp/3._dp*Ly-100._dp/9 )*Lx**3
     & + (-125._dp/3._dp*Ly+39*Ly**2+214._dp/
     & 9+3*pisq-22*Ls )*Lx**2
     & + (14._dp/3._dp*Ly**3-73._dp/3._dp*Ly**2+
     & (-24*pisq+4._dp)*Ly+155._dp/9*pisq+148._dp/3)*Lx
     & -74._dp/9*Ly**3+ (10._dp/3._dp*pisq-55._dp/9-11*Ls )*Ly**2
     & +(-33*Ls+136._dp/9*pisq-140*zeta3+241._dp/3._dp)*Ly
     & -43417._dp/324._dp+23._dp/
     & 6*pisq*Ls-52*zeta3*Ls-173._dp/18*pisq+227._dp/
     & 180*pi**4+1834._dp/27*Ls+515._dp/9*zeta3
     & )*(t/u)
     & + (14*Lx**2+ (-28._dp-28*Ly
     & )*Lx+14*Ly**2+28*Ly+14*pisq
     & )*(t/s)**2-5*Lx**2*(t/u)**2
     & + (-8*Li4y+24*Lx*Li3x+
     & (144*Lx+120._dp)*Li3y-120*Li2x*Lx-20*Ly**2*Li2y
     & -8._dp/3._dp*Lx**3*Ly+ (20*Ly**2+472._dp/9
     & )*Lx**2+ (-182._dp/3._dp*Ly**2+ (-40._dp/
     & 3*pisq-104._dp/3._dp )*Ly+898._dp/9 )*Lx
     & -3*Ly**4-184._dp/9*Ly**3+ (-22*Ls-8._dp/
     & 3*pisq )*Ly**2+ (-22*Ls+56._dp/
     & 9*pisq-136*zeta3 )*Ly
     & +4._dp/3._dp*pi**4+148._dp/9*pisq+2._dp-12*zeta3 )
!     & + \t \leftrightarrow u \
      return
      end
