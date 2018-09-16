      function AGTYD1u(s,t,u,Lx,Ly,Lu,Li2x,Li3x,Li4x,Li3y,Li4y,Li4z)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.15
      include 'types.f'
      real(dp):: AGTYD1u
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Lu,
     & Li4x,Li4y,Li4z,Li3x,Li3y,Li2x

      AGTYD1u=
     & (10*Lx*Ly-5*pisq-5*Ly**2-5*Lx**2
     & )*(t/s)**2 -5*Ly**2*(s/t)**2
     & +14*Lx**2*(t**2+s**2)/u**2-28*Lx*(t**2-s**2)/u**2
     & + (-2*Li4y-2*Li4z+
     & (-8._dp-10*Lx+90*Ly
     & )*Li3x+70*Li3y*Lx
     & +(11*Lx**2+ (8._dp-22*Ly
     & )*Lx+11*pisq )*Li2x-11._dp/
     & 6*Lx**4+ (28._dp/3._dp*Ly-29._dp/3._dp )*Lx**3
     & + (-33._dp/ 2*Lu+53._dp/ 6+47*Ly+13._dp/
     & 2*Ly**2 )*Lx**2+ (22._dp/
     & 3*Ly**3-75*Ly**2+ (-65._dp/
     & 3-15*pisq+33*Lu )*Ly
     & +389._dp/ 6-60*zeta3-31._dp/3._dp*pisq-33._dp/ 2*Lu
     & )*Lx-11._dp/3._dp*Ly**4+50*Ly**3+ (8._dp/
     & 3*pisq+65._dp/3._dp-33*Lu )*Ly**2
     & + (50*zeta3+29._dp/3._dp*pisq-389._dp/
     & 3+33*Lu )*Ly+587._dp/ 9*zeta3+326._dp/
     & 9*pisq-52*zeta3*Lu+61._dp/36._dp*pi**4
     & +1834._dp/27*Lu-43417._dp/324._dp-38._dp/3._dp*pisq*Lu
     & )*(t**2+s**2)/(s*t)
     & +(-96*Li4x-50*Li4y+50*Li4z+
     & (18*Lx+26*Ly-38._dp)*Li3x
     & + (52*Ly-26*Lx-76._dp)*Li3y+
     & (5*Lx**2+ (-2*Ly+38._dp
     & )*Lx-76*Ly+2*Ly**2-13._dp/3._dp*pisq
     & )*Li2x
     & +1._dp/3._dp*Lx**4+ (-2._dp/3._dp*Ly-13._dp/ 9
     & )*Lx**3+ (269._dp/ 18+6*pisq-11._dp/
     & 2*Lu+28*Ly+7._dp/ 2*Ly**2 )*Lx**2
     & + (20._dp/3._dp*Ly**3-66*Ly**2+ (-269._dp/
     & 9-31._dp/3._dp*pisq+11*Lu )*Ly+8._dp/
     & 9*pisq+52*zeta3+33._dp/ 2*Lu-31._dp/ 2 )*Lx
     & +11._dp/3._dp*Ly**2*pisq+ (-26*zeta3+122._dp/
     & 3*pisq )*Ly+38*zeta3+178._dp/
     & 45*pi**4+269._dp/ 18*pisq-11._dp/ 2*pisq*Lu
     & )*(t**2-s**2)/(s*t)
     & + (8*Li4y+8*Li4z+
     & (-120._dp+168*Ly-24*Lx
     & )*Li3x+120*Li3y*Lx
     & + (20*Lx**2+ (120._dp-40*Ly
     & )*Lx+20*pisq )*Li2x-8._dp/3._dp*Lx**4+
     & (40._dp/3._dp*Ly-184._dp/ 9 )*Lx**3
     & + (472._dp/ 9+14*Ly**2+122*Ly-22*Lu
     & )*Lx**2+ (40._dp/3._dp*Ly**3-370._dp/3._dp*Ly**2+
     & (-76._dp/3._dp*pisq+44*Lu-320._dp/ 9 )*Ly
     & -112._dp/ 9*pisq+898._dp/ 9-112*zeta3-22*Lu
     & )*Lx-20._dp/3._dp*Ly**4+740._dp/ 9*Ly**3+ (320._dp/
     & 9-44*Lu )*Ly**2
     & + (44*Lu+218._dp/ 9*pisq+104*zeta3-1796._dp/
     & 9 )*Ly+96*zeta3+104._dp/
     & 45*pi**4+4._dp-22*pisq*Lu+60*pisq )
      return
      end
