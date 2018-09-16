      function AGTYE4u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.21
      include 'types.f'
      include 'zeta.f'
      real(dp):: AGTYE4u
      real(dp):: s,t,Lx,Ly,Lu
      AGTYE4u=(8._dp/3._dp*Lx**3+(8._dp*Lu-26._dp/3._dp-8._dp*Ly)*Lx**2
     & +(8._dp*Ly**2+(52._dp/3._dp-16._dp*Lu)*Ly+8._dp*Lu
     & +8._dp/3._dp*pisq-10._dp)*Lx
     & -16._dp/3._dp*Ly**3+(-52._dp/3._dp+16._dp*Lu)*Ly**2
     & +(8._dp/3._dp*pisq-16._dp*Lu+20._dp)*Ly
     & -112._dp/9._dp*pisq-908._dp/27._dp*Lu-4._dp/9._dp*zeta3
     & +22._dp/3._dp*pisq*Lu+3401._dp/81._dp)*(t**2+s**2)/(s*t)
     & +(8._dp/9._dp*Lx**3+(8._dp/3._dp*Lu-74._dp/9._dp
     & -8._dp/3._dp*Ly)*Lx**2
     & +(8._dp/3._dp*Ly**2+(-16._dp/3._dp*Lu+148._dp/9._dp)*Ly
     & +8._dp/9._dp*pisq+10._dp-8._dp*Lu)*Lx
     & -8._dp/3._dp*pisq*Ly+2._dp/9._dp*pisq
     & *(12._dp*Lu-37._dp))*(t**2-s**2)/(s*t)
     & +(32._dp/9._dp*Lx**3+(-104._dp/9._dp-32._dp/3._dp*Ly
     & +32._dp/3._dp*Lu)*Lx**2
     & +(32._dp/3._dp*Ly**2+(-64._dp/3._dp*Lu+208._dp/9._dp)*Ly
     & -152._dp/9._dp+32._dp/3._dp*Lu+32._dp/9._dp*pisq)*Lx
     & -64._dp/9._dp*Ly**3+(-208._dp/9._dp+64._dp/3._dp*Lu)*Ly**2
     & +(-64._dp/3._dp*Lu+304._dp/9._dp+32._dp/9._dp*pisq)*Ly
     & +8._dp/3._dp*pisq*(-7._dp+4._dp*Lu))

      return
      end

