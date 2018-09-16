      function AGTYAu(s,t,u,Lx,Ly,Li2x,Li3x,Li4x,Li4y,Li3y,Li4z)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.12
      include 'types.f'
      real(dp):: AGTYAu
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,
     & Li4x,Li4y,Li4z,Li3x,Li3y,Li2x

      AGTYAu=(24*pisq-48*Lx*Ly+24*Ly**2+24*Lx**2)*(t/s)**2
     & +24*Ly**2*(s/t)**2
     & +((64*Ly+32._dp/3._dp)*Li3x+64*Li3y*Lx
     & -32._dp/3*Li2x*Lx
     & +(-8._dp+16*Ly**2)*Lx**2
     & +((-64._dp/3._dp*pisq+16._dp)*Ly
     & -16._dp/9._dp*pisq+24._dp-64*zeta3+32._dp/3._dp*Ly**3
     & -32._dp/3._dp*Ly**2 )*Lx+64._dp/9*Ly**3
     & +(-48._dp+64*zeta3+32._dp/9*pisq)*Ly
     & +(-16._dp/3._dp*pisq-16._dp)*Ly**2-16._dp/3._dp*Ly**4
     & -8*pisq+32._dp/3._dp*zeta3
     & +88._dp/45*pi**4 )*(t**2+s**2)/(s*t)
     & +(-128*Li4x-128*Li4y+128*Li4z+
     & (64*Ly+32._dp/3._dp )*Li3x
     & +(128*Ly+64._dp/3._dp-64*Lx)*Li3y
     & +(-32._dp/3._dp*Lx+64._dp/3._dp*Ly+64._dp/3._dp*pisq)*Li2x
     & +16._dp/3._dp*Lx**4-64._dp/3._dp*Lx**3*Ly+
     & (16*Ly**2-8._dp+32._dp/3._dp*pisq )*Lx**2
     & + (32._dp/3._dp*Ly**2+32._dp/
     & 3*Ly**3+64*zeta3-16._dp/9*pisq+24._dp+16*Ly
     & )*Lx+88._dp/15*pi**4-32._dp/3._dp*zeta3
     & + (-32._dp/9*pisq-64*zeta3
     & )*Ly+16*Ly**2*pisq-8*pisq )*(t**2-s**2)/(s*t)
     & + (-32._dp/3._dp*Li3x-64._dp/3._dp*Li3y+
     & (32._dp/3._dp*Lx-64._dp/3._dp*Ly )*Li2x
     & + (-32._dp/3._dp*Ly**2-64._dp/3._dp+16._dp/9*pisq
     & )*Lx+32._dp/9*pisq*Ly+32._dp/3._dp*zeta3)*(t**2-s**2)/u**2
     & +((64*Ly-416._dp/3._dp)*Li3x+64*Li3y*Lx+416._dp/3*Li2x*Lx
     & +(64*Ly+16._dp+16*Ly**2)*Lx**2
     & + (-160._dp/3._dp*Ly**2+32._dp/3._dp*Ly**3+
     & (-80._dp/3._dp-64._dp/3._dp*pisq )*Ly+208._dp/
     & 9*pisq-64*zeta3+80._dp/3._dp )*Lx+320._dp/
     & 9*Ly**3
     & + (-160._dp/3._dp+160._dp/9*pisq+64*zeta3
     & )*Ly+ (-16._dp/3._dp*pisq+80._dp/3._dp
     & )*Ly**2-16._dp/3._dp*Ly**4+88._dp/45*pi**4-200._dp/9*pisq
     & -416._dp/3._dp*zeta3 )
      return
      end

