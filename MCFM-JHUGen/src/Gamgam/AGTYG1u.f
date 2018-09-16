      function AGTYG1u(s,t,Lx,Ly)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.9
      include 'types.f'
      real(dp):: AGTYG1u
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly

      AGTYG1u=(8*Lx**4+(-32*Ly+8._dp)*Lx**3
     & +(48*Ly**2-24*Ly+16*pisq)*Lx**2
     & +(-32*Ly**3+24*Ly**2-32*pisq*Ly+8*pisq)*Lx
     & +8*Ly**4-8*Ly**3+16*Ly**2*pisq-8*pisq*Ly+8*pi**4
     & )*(t/s)**2
     & +(8*Ly**4-8*Ly**3
     & +32*pisq*Ly**2-16*pisq*Ly )*(s/t)**2
     & + (2*Lx**4-8*Lx**3*Ly+
     & (4*pisq+12*Ly**2 )*Lx**2+
     & (-8*Ly**3-8*pisq*Ly)*Lx
     & +2*Ly**4+4*Ly**2*pisq+2*pi**4)*(t/s)**3
     & +(2*Ly**4+8*pisq*Ly**2)*(s/t)**3
     & +(8*Lx**4+(-32*Ly+20._dp)*Lx**3
     & +(16*pisq+56*Ly**2-66*Ly-29._dp)*Lx**2
     & + (-48*Ly**3+78*Ly**2+
     & (58._dp-32*pisq )*Ly-42._dp+20*pisq)*Lx
     & +24*Ly**4-52*Ly**3+ (56*pisq-58._dp)*Ly**2
     & +(-66*pisq+84._dp)*Ly
     & +102._dp+8*pi**4-29*pisq )*(t**2+s**2)/(s*t)
     & + (6*Lx**4+ (-24*Ly+8._dp)*Lx**3
     & +(36*Ly**2-19._dp-30*Ly+12*pisq)*Lx**2
     & +(-24*Ly**3+30*Ly**2+(38._dp-24*pisq)*Ly+8*pisq+42._dp)*Lx
     & -12*Ly**2*pisq+2*pisq*Ly+3*pisq*(2*pisq-3._dp))*(t**2-s**2)/(s*t)
     & +(8*Lx**4+(32._dp-32*Ly)*Lx**3
     & +(-104*Ly+64*Ly**2-32._dp+16*pisq)*Lx**2
     & +(-64*Ly**3+120*Ly**2+(-32*pisq+64._dp)*Ly-56._dp+32*pisq)*Lx
     & +32*Ly**4-80*Ly**3+64*(pi-1._dp)*(pi+1._dp)*Ly**2
     & +(-104*pisq+112._dp)*Ly-8._dp+8*pi**4-32*pisq)
      return
      end


