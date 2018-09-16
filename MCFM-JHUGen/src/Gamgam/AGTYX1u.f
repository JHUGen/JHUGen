      function AGTYX1u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.12
      include 'types.f'
      real(dp):: AGTYX1u
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYX1u= (Lx**3+ (2*Lu-4*Ly+1._dp
     & )*Lx**2+ (6*Ly**2+ (-4._dp-4*Lu
     & )*Ly+pisq+2*Lu-14._dp/3._dp )*Lx
     & -4*Ly**3+ (4._dp+4*Lu )*Ly**2+
     & (-4*Lu+28._dp/3._dp-4*pisq
     & )*Ly+2*pisq*Lu-28._dp/3._dp*Lu+pisq)*(t**2+s**2)/(s*t)
     & + (1._dp/3._dp*Lx**3
     & +(-4._dp/3._dp*Ly-1._dp+2._dp/3._dp*Lu)*Lx**2
     & + (4._dp/3._dp*Ly**2+
     & (-4._dp/3._dp*Lu+2._dp)*Ly+1._dp/
     & 3*pisq-2*Lu )*Lx
     & +1._dp/3._dp*pisq*(2*Lu+3._dp)
     & )*(t**2-s**2)/(s*t)
     & + (4._dp/3._dp*Lx**3+ (-16._dp/3._dp*Ly+8._dp/
     & 3*Lu+4._dp/3._dp )*Lx**2+ (8*Ly**2+
     & (-16._dp/3._dp*Lu-16._dp/3._dp )*Ly+4._dp/
     & 3*pisq+8._dp/3._dp*Lu )*Lx
     & -16._dp/3._dp*Ly**3+ (16._dp/3._dp+16._dp/3._dp*Lu)*Ly**2
     & +(-16._dp/3._dp*Lu-16._dp/3._dp*pisq)*Ly
     & +4._dp/3._dp*pisq*(2*Lu+1._dp))
      return
      end

