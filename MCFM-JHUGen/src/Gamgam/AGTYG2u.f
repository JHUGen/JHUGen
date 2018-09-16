      function AGTYG2u(s,t,Lx,Ly,Lu)
      include 'types.f'
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.10
      real(dp):: AGTYG2u
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYG2u= (-8*Lx**4+ (32*Ly-8._dp)*Lx**3
     & +(-48*Ly**2+24*Ly-16*pisq)*Lx**2
     & +(32*Ly**3-24*Ly**2+32*pisq*Ly-8*pisq)*Lx
     & -8*Ly**4+8*Ly**3-16*Ly**2*pisq+8*pisq*Ly-8*pi**4
     & )*(t/s)**2
     & -(8*Ly**4-8*Ly**3
     & +32*pisq*Ly**2-16*pisq*Ly )*(s/t)**2
     & + (-2*Lx**4+8*Lx**3*Ly+
     & (-4*pisq-12*Ly**2 )*Lx**2+
     & (8*Ly**3+8*pisq*Ly )*Lx
     & -2*Ly**4-4*Ly**2*pisq-2*pi**4 )*(t/s)**3
     & -(2*Ly**4 + 8*pisq*Ly**2) *(s/t)**3
     & + (-5*Lx**4+ (-21._dp+26*Ly)*Lx**3+
     & (-11*Lu-7._dp-50*Ly**2-13*pisq+79*Ly)*Lx**2
     & + (48*Ly**3-111*Ly**2+
     & (6._dp+44*pisq+22*Lu )*Ly+140._dp/
     & 3-30*pisq-11*Lu )*Lx
     & -24*Ly**4+74*Ly**3+
     & (-6._dp-56*pisq-22*Lu )*Ly**2+(-280._dp/3._dp+85*pisq+22*Lu)*Ly
     & +154._dp/3._dp*Lu-4._dp+7*pisq-11*pisq*Lu-8*pi**4)
     & *(t**2+s**2)/(s*t)
     & + (-5*Lx**4+ (22*Ly-43._dp/3._dp)*Lx**3
     & +(121._dp/3*Ly-11*pisq-36*Ly**2+13._dp-11._dp/3._dp*Lu)*Lx**2
     & + (24*Ly**3-121._dp/3._dp*Ly**2+ (-26._dp+22._dp/
     & 3*Lu+20*pisq )*Ly-52._dp/3._dp*pisq+11*Lu)*Lx
     & +12*Ly**2*pisq-5*pisq*Ly-1._dp/3._dp*pisq*(18*pisq+11*Lu-3._dp))
     & *(t**2-s**2)/(s*t)
     & +(-4*Lx**4+ (-88._dp/3._dp+24*Ly)*Lx**3
     & +(-44._dp/3._dp*Lu+8._dp/3._dp-56*Ly**2
     & -12*pisq+340._dp/3._dp*Ly )*Lx**2
     & + (64*Ly**3-164*Ly**2+ (64._dp/
     & 3+48*pisq+88._dp/3._dp*Lu )*Ly+28._dp-124._dp/
     & 3*pisq-44._dp/3._dp*Lu )*Lx
     & -32*Ly**4+328._dp/3._dp*Ly**3
     & +(-64._dp/3-64*pisq-88._dp/3._dp*Lu)*Ly**2
     & +(-56._dp+364._dp/3._dp*pisq+88._dp/3._dp*Lu )*Ly
     & +8._dp+8._dp/3._dp*pisq-44._dp/3._dp*pisq*Lu-8*pi**4
     & )
      return
      end
