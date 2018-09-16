      function AGTYCu(s,t,u,Lx,Ly,Lu,Li2x,Li3x,Li4x,Li4y,Li3y,Li4z)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.14
      include 'types.f'
      real(dp):: AGTYCu
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Lu,
     & Li4x,Li4y,Li4z,Li3x,Li3y,Li2x

      AGTYCu= (-Lx*Ly+1._dp/2._dp*pisq+1._dp/
     & 2*Ly**2+1._dp/ 2*Lx**2 )*(t/s)**2 +1._dp/
     & 2*Ly**2*(s/t)**2
     & -5._dp/ 4*Lx**2*(t**2+s**2)/u**2+5._dp/ 2*Lx*(t**2-s**2)/u**2
     & + (-14*Li4y-14*Li4z+ (-59._dp/
     & 6+11*Lx-31*Ly )*Li3x-9*Li3y*Lx
     & + (-4*Lx**2+ (59._dp/ 6+8*Ly
     & )*Lx-4*pisq )*Li2x-1._dp/ 4*Lx**4+
     & (4._dp/3._dp*Ly+13._dp/ 18 )*Lx**3
     & + (-7._dp/ 2*Ly**2-22._dp/3._dp*Ly-5._dp/
     & 2*pisq+11._dp/ 4*Lu+14._dp/3._dp )*Lx**2
     & + (-4._dp/3._dp*Ly**3+55._dp/ 2*Ly**2+ (35._dp/
     & 6*pisq-33._dp/ 2*Lu )*Ly+2*zeta3+55._dp/
     & 3*Lu-847._dp/ 54+1._dp/ 18*pisq )*Lx
     & +2._dp/3._dp*Ly**4-55._dp/3._dp*Ly**3+ (-1._dp/
     & 2*pisq+33._dp/ 2*Lu )*Ly**2+ (847._dp/
     & 27-28._dp/ 9*pisq+5*zeta3-110._dp/3._dp*Lu
     & )*Ly
     & -1142._dp/ 81+61._dp/ 9*zeta3-13*Lu-166._dp/
     & 9*pisq+47._dp/360._dp*pi**4-2*zeta3*Lu+143._dp/
     & 18*pisq*Lu+121._dp/ 12*Lu**2 )*(t**2+s**2)/(s*t)
     & +
     & (20*Li4x+14*Li4y-14*Li4z+
     & (-7*Ly+19._dp/ 2-Lx )*Li3x
     & + (-14*Ly+7*Lx+19._dp)*Li3y+
     & (-2*Lx**2-19._dp/ 2*Lx+2._dp/3._dp*pisq+19*Ly
     & )*Li2x
     & -1._dp/ 4*Lx**4+ (2._dp/3._dp*Ly+13._dp/ 18
     & )*Lx**3+ (-3._dp/ 2*Ly**2-38._dp/ 9-3._dp/
     & 2*pisq-10*Ly+11._dp/ 4*Lu )*Lx**2
     & + (-4._dp/3._dp*Ly**3+39._dp/ 2*Ly**2+ (76._dp/
     & 9+13._dp/ 6*pisq-11._dp/ 2*Lu )*Ly+77._dp/
     & 36*pisq-31._dp/ 6-12*zeta3 )*Lx
     & -3._dp/ 2*Ly**2*pisq+ (7*zeta3-79._dp/
     & 6*pisq )*Ly-19._dp/ 2*zeta3-38._dp/
     & 9*pisq-5._dp/ 6*pi**4+11._dp/ 4*pisq*Lu
     & )*(t**2-s**2)/(s*t)
     & + (-24*Li4y-24*Li4z+
     & (22._dp-60*Ly+20*Lx
     & )*Li3x-20*Li3y*Lx
     & + (-8*Lx**2+ (-22._dp+16*Ly
     & )*Lx-8*pisq )*Li2x-2._dp/3._dp*Lx**4+
     & (4._dp/3._dp*Ly+59._dp/ 9 )*Lx**3
     & + (-35*Ly-4._dp/3._dp*pisq-575._dp/
     & 36-6*Ly**2+11*Lu )*Lx**2
     & + (-8._dp/3._dp*Ly**3+131._dp/3._dp*Ly**2+ (26._dp/
     & 3*pisq-22*Lu+178._dp/ 9 )*Ly+11*Lu-637._dp/
     & 18+24*zeta3+53._dp/ 9*pisq )*Lx
     & +4._dp/3._dp*Ly**4-262._dp/ 9*Ly**3+ (-178._dp/
     & 9+2._dp/3._dp*pisq+22*Lu )*Ly**2+ (637._dp/
     & 9-28*zeta3-67._dp/ 9*pisq-22*Lu )*Ly
     & -18*zeta3+2._dp/ 45*pi**4+11*pisq*Lu-71._dp/
     & 3*pisq )
      return
      end
