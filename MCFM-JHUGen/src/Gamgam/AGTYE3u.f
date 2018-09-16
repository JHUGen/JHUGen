      function AGTYE3u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.20
      include 'types.f'
      real(dp):: AGTYE3u
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYE3u= (-8._dp/3._dp*Ly**3+(4*Lu-26._dp/3._dp+4*Lx)*Ly**2
     & +(-4*Lx**2+(26._dp/3._dp-4*Lu)*Lx
     & +4._dp/3._dp*pisq+10._dp-4*Lu)*Ly
     & +4._dp/3._dp*Lx**3+(2*Lu-13._dp/3._dp)*Lx**2
     & +(4._dp/3._dp*pisq+2*Lu-5._dp)*Lx
     & -2._dp/9._dp*zeta3+5._dp/3._dp*pisq*Lu+3401._dp/162._dp
     & -56._dp/9._dp*pisq-202._dp/27._dp*Lu)*(t**2+s**2)/(s*t)
     & +(4._dp/3._dp*Ly**2*Lx+ (-4._dp/3._dp*Lx**2
     & +(-4._dp/3._dp*Lu+74._dp/9._dp)*Lx-4._dp/3._dp*pisq)*Ly
     & +4._dp/9*Lx**3+(-37._dp/9._dp+2._dp/3._dp*Lu)*Lx**2
     & +(-2*Lu+5._dp+4._dp/9*pisq)*Lx
     & +1._dp/9._dp*pisq*(-37._dp+6*Lu))*(t**2-s**2)/(s*t)
     & +(-32._dp/9*Ly**3
     & +(16._dp/3._dp*Lu+16._dp/3._dp*Lx-104._dp/9)*Ly**2
     & +(-16._dp/3._dp*Lx**2+(-16._dp/3._dp*Lu+104._dp/9._dp)*Lx
     & +152._dp/9-16._dp/3._dp*Lu+16._dp/9*pisq)*Ly
     & +16._dp/9*Lx**3+(8._dp/3._dp*Lu-52._dp/9._dp)*Lx**2
     & +(8._dp/3._dp*Lu-76._dp/9+16._dp/9._dp*pisq)*Lx
     & +4._dp/3._dp*pisq*(2*Lu-7._dp))
      return
      end

