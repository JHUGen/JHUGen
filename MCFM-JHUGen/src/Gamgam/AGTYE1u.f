      function AGTYE1u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.16
      include 'types.f'
      real(dp):: AGTYE1u
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYE1u= (11._dp/6*Lx**3+(3*Lu-23._dp/6-6*Ly )*Lx**2
     & +(7*Ly**2+(-6*Lu+20._dp/3._dp)*Ly+11._dp/6*pisq-22._dp/3+3*Lu)*Lx
     & -14._dp/3._dp*Ly**3+(-20._dp/3._dp+6*Lu)*Ly**2
     & +(-6*Lu+44._dp/3._dp+4._dp/3._dp*pisq)*Ly
     & -2._dp/9*zeta3+8._dp/3._dp*pisq*Lu-121._dp/18*pisq
     &  +3401._dp/162._dp-328._dp/27._dp*Lu)*(t**2+s**2)/(s*t)
     & +(11._dp/18._dp*Lx**3+(-83._dp/18._dp+Lu-2*Ly)*Lx**2
     & +(2*Ly**2+(83._dp/9-2*Lu)*Ly-3*Lu+5._dp+11._dp/18._dp*pisq)*Lx
     & -2*pisq*Ly+1._dp/18*pisq
     & *(-83._dp+18*Lu))*(t**2-s**2)/(s*t)
     & + (22._dp/9*Lx**3+ (-46._dp/
     & 9+4*Lu-8*Ly )*Lx**2+ (28._dp/
     & 3*Ly**2+ (80._dp/9-8*Lu )*Ly+22._dp/
     & 9*pisq-76._dp/9+4*Lu )*Lx
     & -56._dp/9*Ly**3+ (-80._dp/9+8*Lu
     & )*Ly**2+ (152._dp/9-8*Lu+16._dp/9*pisq
     & )*Ly+2*pisq*(-5._dp+2*Lu))
      return
      end
