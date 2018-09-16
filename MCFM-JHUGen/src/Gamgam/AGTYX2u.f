      function AGTYX2u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.13
      include 'types.f'
      real(dp):: AGTYX2u
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYX2u= (-1._dp/ 6*Lx**3+ (-1._dp/
     & 3*Lu+4._dp/3._dp*Ly-10._dp/ 9 )*Lx**2+
     & (-3*Ly**2+ (2*Lu+40._dp/ 9 )*Ly-7._dp/
     & 6*pisq-31._dp/ 9*Lu )*Lx
     & +2*Ly**3+ (-2*Lu-40._dp/ 9
     & )*Ly**2+ (2*pisq+62._dp/ 9*Lu
     & )*Ly-22._dp/ 9*Lu**2-pisq*Lu-10._dp/ 9*pisq
     & )*(t**2+s**2)/(s*t)
     & + (-1._dp/ 6*Lx**3+ (2._dp/3._dp*Ly-1._dp/
     & 3*Lu )*Lx**2+ (-2._dp/3._dp*Ly**2-1._dp/
     & 6*pisq+2._dp/3._dp*Ly*Lu )*Lx-1._dp/
     & 3*pisq*Lu )*(t**2-s**2)/(s*t)
     & + (-2._dp/3._dp*Lx**3+ (8._dp/3._dp*Ly-4._dp/
     & 3*Lu-2._dp/3._dp )*Lx**2+ (-4*Ly**2+
     & (8._dp/3._dp*Lu+8._dp/3._dp )*Ly-2._dp/3._dp*pisq-4._dp/
     & 3*Lu )*Lx
     & +8._dp/3._dp*Ly**3+ (-8._dp/3._dp-8._dp/3._dp*Lu
     & )*Ly**2+ (8._dp/3._dp*Lu+8._dp/3._dp*pisq
     & )*Ly-2._dp/3._dp*pisq*(2*Lu+1._dp))
      return
      end
