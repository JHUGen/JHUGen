      function AGTYCs(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,
     & Li4z,Li4zinv)
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.3
      implicit none
      include 'types.f'
      real(dp):: Csx,AGTYCs
      real(dp)::s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,
     & Li4z,Li4zinv
      AGTYCs=Csx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li4z)
     &      +Csx(s,u,t,Ly,Lx,Ls,Li2y,Li3y,Li4y,Li2x,Li3x,Li4x,Li4zinv)
      return
      end


      function Csx(s,t,u,Lx,Ly,Ls,Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li4z)
      implicit none
      include 'types.f'
      real(dp):: Csx
      include 'zeta.f'
      real(dp)::s,t,u,Lx,Ly,Ls,
     & Li4x,Li4y,Li4z,Li3x,Li3y,Li2x,Li2y
      Csx= (-20*Li4z+28*Li4x+
     & (-28*Ly-10*Lx+1._dp/3._dp )*Li3x
     & + (6*Lx**2+ (-1._dp/
     & 3-4*Ly )*Lx-2*Ly**2+58._dp/3._dp*Ly+4._dp/
     & 3*pisq )*Li2x
     & +(58._dp/
     & 3-12*Lx-12*Ly)*Li3y-1._dp/
     & 6*Lx**4+ (10._dp/3._dp*Ly+13._dp/ 9 )*Lx**3
     & + (-9*Ly**2-1._dp/3._dp*pisq+4._dp/ 9+11._dp/
     & 2*Ls+13*Ly )*Lx**2
     & + (-8._dp/3._dp*Ly**3+50._dp/3._dp*Ly**2+ (17._dp/
     & 3*pisq-28._dp/3._dp+11*Ls )*Ly-563._dp/
     & 27-233._dp/36._dp*pisq+55._dp/3._dp*Ls )*Lx
     & + (-2._dp/3._dp*pisq+80._dp/ 9 )*Ly**2+
     & (-284._dp/ 27-299._dp/36._dp*pisq+26*zeta3+55._dp/
     & 3*Ls )*Ly-209._dp/36._dp*pisq*Ls
     & -2*zeta3*Ls+121._dp/ 12*Ls**2-13*Ls-1142._dp/
     & 81-197._dp/360._dp*pi**4+461._dp/36._dp*pisq-55._dp/
     & 18*zeta3 )*(t/u)
     & + (-5._dp/ 4*Lx**2+ (5._dp/ 2+5._dp/
     & 2*Ly )*Lx-5._dp/ 4*Ly**2-5._dp/ 2*Ly-5._dp/
     & 4*pisq )*(t/s)**2+1._dp/ 2*Lx**2*(t/u)**2
     & + (24*Li4y-20*Lx*Li3x+
     & (-40*Lx-22._dp)*Li3y+22*Li2x*Lx+8*Ly**2*Li2y
     & +4._dp/3._dp*Lx**3*Ly+ (-6*Ly**2-575._dp/36._dp
     & )*Lx**2+ (46._dp/3._dp*Ly**2+ (73._dp/
     & 12+4*pisq )*Ly-637._dp/ 18 )*Lx
     & +1._dp/3._dp*Ly**4+59._dp/ 9*Ly**3+ (2._dp/
     & 3*pisq+11*Ls )*Ly**2+
     & (44*zeta3-4._dp/ 9*pisq+11*Ls )*Ly
     & -38._dp/ 45*pi**4+77._dp/ 72*pisq+2*zeta3 )
!     & + \t \leftrightarrow u \
      return
      end
