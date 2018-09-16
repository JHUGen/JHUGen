      function AGTYD2u(s,t,u,Lx,Ly,Lu,Li2x,Li3x,Li4x,Li3y,Li4y,Li4z)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.19
      include 'types.f'
      real(dp):: AGTYD2u
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Lu,
     & Li4x,Li4y,Li4z,Li3x,Li3y,Li2x
      AGTYD2u=
     & (-6*pisq-6*Ly**2+12*Lx*Ly-6*Lx**2
     & )*(t/s)**2 -6*Ly**2*(s/t)**2
     & +6*Lx**2*(t**2+s**2)/u**2-12*Lx*(t**2-s**2)/u**2
     & + (-4*Li4y-4*Li4z+
     & (-4._dp-4*Lx+36*Ly
     & )*Li3x+28*Li3y*Lx
     & + (6*Lx**2-12*Lx*Ly+6*pisq+4*Lx
     & )*Li2x+20*Ly**3+ (107._dp/
     & 3+3*Lx**2-22*LU-30*Lx )*Ly**2
     & + (4*Lx**3+24*Lx**2+
     & (-4*pisq+22*LU-107._dp/3._dp
     & )*Lx-57._dp+36*zeta3+22*LU-2._dp/3._dp*pisq
     & )*Ly
     & -Lx**4-19._dp/3._dp*Lx**3+ (-11*LU+107._dp/6+1._dp/
     & 3*pisq )*Lx**2+
     & (-11*LU-32*zeta3+57._dp/2-20._dp/3._dp*pisq
     & )*Lx
     & -43417._dp/324._dp-43._dp/6*pisq*LU
     & -52*zeta3*LU+515._dp/9*zeta3+251._dp/
     & 9*pisq+1141._dp/27*LU+65._dp/36._dp*pi**4
     & )*(t**2+s**2)/(s*t)
     & +(-20*Li4y-48*Li4x+20*Li4z+
     & (12*Ly-16._dp+12*Lx )*Li3x
     & +(4*Ly**2+ (-4*Lx-32._dp)*Ly
     & +2*Lx**2+16*Lx-10._dp/3._dp*pisq)*Li2x
     & + (-12*Lx+24*Ly-32._dp
     & )*Li3y+4*Ly**3*Lx+ (-94._dp/
     & 3*Lx+Lx**2+4._dp/3._dp*pisq )*Ly**2
     & + (46._dp/3._dp*Lx**2+ (-20._dp/3._dp*pisq+22._dp/
     & 3*LU-251._dp/9 )*Lx-12*zeta3+62._dp/
     & 3*pisq )*Ly
     & -13._dp/9*Lx**3+ (251._dp/18+3*pisq-11._dp/
     & 3*LU )*Lx**2+ (-16._dp/
     & 9*pisq+24*zeta3-57._dp/2+11*LU )*Lx
     & +16*zeta3+251._dp/18*pisq-11._dp/
     & 3*pisq*LU+94._dp/45*pi**4 )*(t**2-s**2)/(s*t)
     & + (-16*Li4y-16*Li4z+
     & (-64._dp+48*Ly )*Li3x+48*Li3y*Lx
     & +
     & (-16*Lx*Ly+64*Lx+8*Lx**2+8*pisq
     & )*Li2x+272._dp/9*Ly**3
     & + (4*Lx**2+344._dp/9-136._dp/3._dp*Lx-88._dp/
     & 3*LU )*Ly**2+ (8*Lx**3+184._dp/
     & 3*Lx**2+ (88._dp/3._dp*LU-344._dp/9-16._dp/
     & 3*pisq )*Lx
     & +88._dp/3._dp*LU+32._dp/9*pisq+48*zeta3-1076._dp/
     & 9 )*Ly-2*Lx**4-112._dp/9*Lx**3+
     & (-44._dp/3._dp*LU+298._dp/9 )*Lx**2
     & + (-76._dp/9*pisq-44._dp/3._dp*LU+538._dp/
     & 9-48*zeta3 )*Lx+48*pisq+48*zeta3-44._dp/
     & 3*pisq*LU+92._dp/45*pi**4 )
      return
      end
