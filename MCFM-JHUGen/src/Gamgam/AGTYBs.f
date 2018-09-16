      function AGTYBs(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,
     & Li4z,Li4zinv)
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.2
      implicit none
      include 'types.f'
      real(dp):: Bsx,AGTYBs
      real(dp)::s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,
     & Li4z,Li4zinv
      AGTYBs=Bsx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li2y,Li3y,Li4y,Li4z)
     &      +Bsx(s,u,t,Ly,Lx,Ls,Li2y,Li3y,Li2x,Li3x,Li4x,Li4zinv)
      return
      end


      function Bsx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li2y,Li3y,Li4y,Li4z)
      implicit none
      include 'types.f'
      real(dp):: Bsx
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Ls,
     & Li4y,Li4z,Li3x,Li3y,Li2x,Li2y
      Bsx= (-112*Li4z-88*Li4y+
     & (-128*Ly+48*Lx-64._dp)*Li3x
     & + (-16*Ly-16*Lx+12._dp)*Li3y+
     & (12*Ly-4*Ly**2+8*Lx**2-8*pisq+64*Lx
     & )*Li2x
     & +2._dp/3._dp*Lx**4+56._dp/3._dp*Lx**3*Ly+
     & (44*Ly-4*pisq+2._dp-32*Ly**2 )*Lx**2
     & +(-4*Ly**3-8._dp-32*zeta3-80._dp/3._dp*pisq
     & +6*Ly**2+56._dp/3._dp*pisq*Ly )*Lx
     & +Ly**4+6*Ly**3+ (-10._dp/3._dp*pisq-5._dp)*Ly**2
     & +(-39._dp-18*pisq+144*zeta3)*Ly
     & +3*Ls+187._dp/4-4*pisq*Ls+4._dp/45._dp*pi**4
     & -5*pisq-20*zeta3+48*zeta3*Ls)*(t/u)
     & +(-12*Lx**2+(24*Ly+24._dp)*Lx
     & -12*Ly**2-24*Ly-12*pisq)*(t/s)**2
     & +8*Lx**2*(t/u)**2
     & +(-80*Li4y+32*Lx*Li3x+
     & (-128*Lx-152._dp)*Li3y
     & +152*Li2x*Lx
     & +8*Ly**2*Li2y+ (-16*Ly**2-24._dp)*Lx**2+
     & (60*Ly**2+(28._dp+32._dp/3._dp*pisq)*Ly-58._dp)*Lx
     & +14._dp/3._dp*Ly**4+44._dp/3._dp*Ly**3+8._dp/3._dp*Ly**2*pisq
     & +(96*zeta3-32._dp/3._dp*pisq)*Ly
     & +32._dp/45._dp*pi**4+16*zeta3-86._dp/3._dp*pisq-2._dp)
!     & & & + \t \leftrightarrow u \
      return
      end

