      function AGTYG3u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.11
      include 'types.f'
      real(dp):: AGTYG3u
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYG3u= (2*Lx**4+ (-8*Ly+2._dp)*Lx**3
     & +(12*Ly**2-6*Ly+4*pisq)*Lx**2
     & +(-8*Ly**3+6*Ly**2-8*pisq*Ly+2*pisq)*Lx
     & +2*Ly**4-2*Ly**3+4*Ly**2*pisq-2*pisq*Ly+2*pi**4)*(t/s)**2
     & +(2*Ly**4-2*Ly**3
     & +8*pisq*Ly**2-4*pisq*Ly )*(s/t)**2
     & + (1._dp/ 2*Lx**4-2*Lx**3*Ly+
     & (pisq+3*Ly**2 )*Lx**2+
     & (-2*Ly**3-2*pisq*Ly )*Lx
     & +1._dp/ 2*Ly**4+Ly**2*pisq+1._dp/ 2*pi**4
     & )*(t/s)**3 +(1._dp/ 2*Ly**4 + 2*pisq*Ly**2)
     & *(s/t)**3
     & + (Lx**4+ (11._dp/3._dp-5*Ly
     & )*Lx**3+ (11._dp/ 6*Lu+59._dp/ 9+11*Ly**2+9._dp/
     & 2*pisq-58._dp/3._dp*Ly )*Lx**2
     & + (-12*Ly**3+36*Ly**2+
     & (-14*pisq-218._dp/ 9-11*Lu )*Ly+110._dp/
     & 9*Lu+41._dp/3._dp*pisq )*Lx
     & +6*Ly**4-24*Ly**3+ (11*Lu+218._dp/
     & 9+14*pisq )*Ly**2+ (-26*pisq-220._dp/
     & 9*Lu )*Ly
     & +121._dp/ 18*Lu**2+2*pi**4+11._dp/ 2*pisq*Lu+59._dp/
     & 9*pisq+2._dp)*(t**2+s**2)/(s*t)
     & + (Lx**4+ (11._dp/3._dp-5*Ly
     & )*Lx**3+ (1._dp-38._dp/3._dp*Ly+5._dp/
     & 2*pisq+9*Ly**2+11._dp/ 6*Lu )*Lx**2
     & + (-6*Ly**3+38._dp/3._dp*Ly**2+
     & (-4*pisq-11._dp/3._dp*Lu-2._dp)*Ly
     & +11._dp/3._dp*pisq)*Lx
     & -3*Ly**2*pisq+2*pisq*Ly+1._dp/
     & 6*pisq*(9*pisq+11*Lu-6._dp)
     & )*(t**2-s**2)/(s*t)
     & + ( (20._dp/3._dp-4*Ly )*Lx**3+
     & (12*Ly**2-92._dp/3._dp*Ly+20._dp/3._dp+22._dp/
     & 3*Lu+2*pisq )*Lx**2
     & + (-16*Ly**3+52*Ly**2+ (-44._dp/
     & 3*Lu-80._dp/3._dp-16*pisq )*Ly+22._dp/
     & 3*Lu+38._dp/3._dp*pisq )*Lx
     & +8*Ly**4-104._dp/3._dp*Ly**3+ (16*pisq+80._dp/
     & 3+44._dp/3._dp*Lu )*Ly**2+ (-104._dp/
     & 3*pisq-44._dp/3._dp*Lu )*Ly
     & +2._dp/3._dp*pisq*(3*pisq+10._dp+11*Lu))
      return
      end


