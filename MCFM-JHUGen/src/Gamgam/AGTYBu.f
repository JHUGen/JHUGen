      function AGTYBu(s,t,u,Lx,Ly,Lu,Li2x,Li3x,Li4x,Li3y,Li4y,Li4z)
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.13
      implicit none
      include 'types.f'
      real(dp):: AGTYBu
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Lu,
     & Li4x,Li4y,Li4z,Li3x,Li3y,Li2x

      AGTYBu=
     & -12*Lx**2*(t**2+s**2)/u**2+24*Lx*(t**2-s**2)/u**2
     & +8*Ly**2*(s/t)**2+(-16*Lx*Ly+8*Ly**2+8*pisq+8*Lx**2 )*(t/s)**2
     & +(44*Li4z+44*Li4y+(-16*Lx-56*Ly+26._dp)*Li3x-88*Li3y*Lx
     & +(-6*Lx**2+(12*Ly-26._dp)*Lx-6*pisq)*Li2x+5*Lx**4
     & +(-20*Ly+3._dp)*Lx**3
     & +(-3._dp/2._dp+10._dp/3._dp*pisq-28*Ly+Ly**2)*Lx**2
     & +(-28._dp/3._dp*Ly**3+40*Ly**2+(20._dp/3._dp*pisq+3._dp)*Ly
     & -1._dp/3._dp*pisq-47._dp/2+72*zeta3)*Lx
     & +14._dp/3._dp*Ly**4-80._dp/3._dp*Ly**3
     & +(-3._dp-10._dp/3._dp*pisq)*Ly**2
     & +(47._dp-56*zeta3-28._dp/3._dp*pisq)*Ly
     & -13._dp/2._dp*pisq-46*zeta3-4*pisq*Lu-28._dp/9._dp*pi**4
     & +3*Lu+48*zeta3*Lu+187._dp/4._dp)*(t**2+s**2)/(s*t)

     & +(-44*Li4z+44*Li4y+112*Li4x
     & +(-32*Lx+38._dp-24*Ly)*Li3x+(24*Lx-48*Ly+76._dp)*Li3y
     & +(-2*Lx**2+(-38._dp+4*Ly)*Lx+76*Ly-4*Ly**2+6*pisq)*Li2x
     & +1._dp/3._dp*Lx**4+(-4._dp/3._dp*Ly-3._dp)*Lx**3
     & +(-6*pisq-Ly**2+7._dp/2._dp-16*Ly)*Lx**2
     & +(-8*Ly**3+54*Ly**2+(12*pisq-7._dp)*Ly
     & -7._dp/3._dp*pisq-56*zeta3+31._dp/2._dp)*Lx
     & -4._dp/3._dp*Ly**2*pisq+(24*zeta3-86._dp/3._dp*pisq)*Ly
     & -206._dp/45._dp*pi**4-38*zeta3+7._dp/2._dp*pisq)
     & *(t**2-s**2)/(s*t)

     & +(80*Li4z+80*Li4y+(-96*Ly-32*Lx+152._dp)*Li3x-160*Li3y*Lx
     & +(-8*Lx**2+(-152._dp+16*Ly)*Lx-8*pisq)*Li2x
     & +8*Lx**4+(44._dp/3._dp-32*Ly )*Lx**3
     & +(16._dp/3._dp*pisq-4*Ly**2-104*Ly-24._dp)*Lx**2
     & +(-16*Ly**3+72*Ly**2+(16*pisq-8._dp)*Ly
     & +4._dp/3._dp*pisq+128*zeta3-58._dp)*Lx
     & +8*Ly**4-48*Ly**3+(8._dp-8._dp/3._dp*pisq)*Ly**2
     & +(116._dp-56._dp/3._dp*pisq-96*zeta3)*Ly
     & -4._dp-76._dp/3._dp*pisq-120*zeta3-24._dp/5._dp*pi**4)
      return
      end
